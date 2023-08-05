#####################################################################
# Functions written to implement 3D CGR for DNA sequences 
# By: Stephanie Young <syoung2@sdsu.edu> <syoung49@its.jnj.com>
#####################################################################

library(stringr)
library(hypervolume)
library(Rcpp)

#=====================================================================
# Convert a DNA sequence to 3D CGR 
#=====================================================================

switch.V <- 
  function(x, ...){
    sapply(x, function(EXPR) switch(EXPR, ...))
  }

read_dict <- 
  function(x){
    # The order of this function must stay fixed as this is the order used in 
    # the C++ implementation of the chaos game. 
    switch(
      x, 
      R = c("A", "G"),
      Y = c("C", "T"),
      K = c("G", "T"),
      M = c("A", "C"),
      S = c("C", "G"),
      W = c("A", "T"),
      B = c("C", "G", "T"),
      D = c("A", "G", "T"),
      H = c("A", "C", "T"),
      V = c("A", "C", "G"), 
      N = c("A", "C", "G", "T")
    )
  }

# Coordinate table
## Function to make a table containing 3DCGR coordinates
CGR_table <- 
  function(A_, C_, G_, T_){
    unknown_reads <- c("R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "N")
    output <- matrix(nrow = 4 + length(unknown_reads), ncol = length(A_)) 
    
    output[1:4, ] <- rbind(A_, C_, G_, T_) 

    idx <- purrr::partial(switch.V, A = 1, C = 2, G = 3, T = 4)
    
    row_index <- 5
    for(u_r in unknown_reads){
      output[row_index,] <- output[idx(read_dict(u_r)),] |> colMeans() |> t()
      row_index <- row_index + 1
    }
    
    output <- output |> as.data.frame()
    output[["Nucleotide"]] <- c("A", "C", "G", "T", unknown_reads)
    
    return(output)
  }

coord1 <- CGR_table(
  (c(0, 0, 0) - (1/2)) * 2*sqrt(1/3),
  (c(1, 0, 1) - (1/2)) * 2*sqrt(1/3),
  (c(0, 1, 1) - (1/2)) * 2*sqrt(1/3),
  (c(1, 1, 0) - (1/2)) * 2*sqrt(1/3)
)

coord2 <- CGR_table(
  c(0, 0, 1),
  c(-sqrt(2)/3, sqrt(6)/3, -1/3),
  c(-sqrt(2)/3, -sqrt(6)/3, -1/3),
  c(2*sqrt(2) / 3, 0, -1/3)
)

sourceCpp("./R/chaosgame.cpp") 

# (recursive implementation)
seq_to_hypercomplex_cg <- 
  function(dna_seq, CGR_coord = coord1, df = F, axes = c("i", "j", "k")){ 
    # The input to this function is the DNA sequence and 
    # a table of the CGR coordinates to be used
    if(ncol(CGR_coord) != (length(axes) + 1))
      stop("Number of axes needs to match that of the CGR_coord") 
    
    if(length(dna_seq) == 1) dna_seq <- str_split(toupper(dna_seq), "")[[1]]
    
    # Remove gaps
    dna_seq <- dna_seq[which(dna_seq != "-")]
    
    # Chaos game 
    cg <- chaosGame(dna_seq, as.matrix(CGR_coord[,1:length(axes)]))
    colnames(cg) <- axes
    return(cg)
  }

#=====================================================================
# Functions for shape signature methods
#=====================================================================

sourceCpp("./R/features.cpp")

## Functions for building histograms for multiple features.
## Not really used, but may be helpful for future implementations using 
## joint probability distribution information 

### One histogram for each feature
hist_df_unpaired <- function(df, breaks){
  if(length(breaks) != ncol(df)) "hist_df_unpaired being used incorrectly"
  sapply(1:ncol(df), function(x) hist(df[,x], breaks[[x]], plot = F)$counts)
}

### Functions for building joint histograms for multiple features
cutting_df <- function(df, breaks){
  # df is a data frame with n columns 
  # breaks is a data frame with n columns
  output <- data.frame(matrix(ncol = 0, nrow = nrow(df)))
  for(i in 1:ncol(df)){
    output[[paste0("interval_", i)]] <- cut(df[,i], breaks[,i], include.lowest = T)
  }
  
  output |> group_by(across(starts_with("interval_"))) |> 
    summarise(n = n(), .groups = "drop")
}

hist_paired <- function(list_df, bin_count, return_bins = F){
  # Inputs are a list of features that have been column bound;
  # Each element in a list is an organism
  
  # cut points  
  cut_lb <- apply(sapply(list_df, function(x) apply(x, 2, min)), 1, min)
  cut_ub <- apply(sapply(list_df, function(x) apply(x, 2, max)), 1, max)
  cut_points <- 
    sapply(
      1:ncol(list_df[[1]]), 
      function(i) seq(cut_lb[i], cut_ub[i], length = bin_count + 1)
    ) 
  
  output <- 
    lapply(
      seq_along(list_df), 
      function(x) cutting_df(list_df[[x]], breaks = cut_points)
    )
  
  cut_point_centers <- apply(cut_points, 2, zoo::rollmean, k = 2) |> 
    as_tibble() |> 
    rename_with(.f = partial(gsub, pattern = "V", replacement = "interval_")) 
  
  output <- 
    Reduce(function(...) full_join(by = grep("interval_", colnames(output[[1]]), value = T), ...), output) |> 
    mutate(across(where(is.numeric), replace_na, 0)) |> 
    mutate(across(contains("interval_"), .fns = function(x) cut_point_centers[[cur_column()]][x]))
  
  counts <- output |> select(!starts_with("interval_")) |> as.matrix() |> t()
  rownames(counts) <- names(list_df)
  if(!return_bins) return(counts)
  else return(list(counts, output |> select(starts_with("interval"))))
}

#=====================================================================
# Function to implement feature signature methods 
# by building histograms for features (such as angles and edge lengths) 
# and/or coordinates
#=====================================================================


feature_histograms <- 
  function(features, bin_count, return_bins = F, return_bin_centers = F, drop_empty = F){
    if(max(unlist(features)) == pi & min(unlist(features)) == -pi){
      # Histograms for angles may not need to be uniformly distributed. 
      # At least on one instance, which this is accounting for, the angles 
      # immediately adjacent to pi and -pi cannot be attained in the CGR. 
      # This removes the unattainable angles next to pi and -pi out of the set of 
      # histogram breaks. 
      remove_pi <- function(x){
        new_x <- x[which(abs(pi - x) > 3e-08)]
        new_x <- new_x[which(abs(-pi - new_x) > 3e-08)]
        return(new_x)
      }
      not_pi <- remove_pi(unlist(features))
      features_breaks <- seq(min(not_pi), max(not_pi), length.out = bin_count - 1)
      features_breaks <- c(-pi, features_breaks, pi)
    } else {
      features_breaks <- seq(min(unlist(features)), max(unlist(features)), length.out = bin_count + 1)
    }
    # Get histogram bin counts
    features_tabulation <-
      sapply(
        features,
        function(x) hist(x, breaks = features_breaks, plot = F)$counts
      )
    bin_centers <- zoo::rollmean(features_breaks, k = 2)
    
    # Remove bins where it's zero across all organisms
    if(drop_empty){
      tabulations.rowSums <- rowSums(features_tabulation)
      features_tabulation <- features_tabulation[which(tabulations.rowSums != 0),,drop = F]
      bin_centers <- bin_centers[which(tabulations.rowSums != 0)]
    }
    
    output <- t(features_tabulation)
    if(return_bins) output <- list(output, features_breaks)
    if(return_bin_centers) output <- list(output, bin_centers)
    return(output)
  }

# feature_histograms <- function(dna_features, bin_count, return_bins = F, drop_empty = F){
#   features_breaks <- seq(min(unlist(dna_features)), max(unlist(dna_features)),
#                          length.out = bin_count + 1)
# 
#   # Get histogram bin counts
#   features_tabulation <-
#     sapply(dna_features,
#            function(x) hist(x, breaks = features_breaks, plot = F)$counts
#     )
#   bin_centers <- zoo::rollmean(features_breaks, k = 2)
# 
#   # Remove bins where it's zero across all organisms
#   if(drop_empty){
#     tabulations.rowSums <- rowSums(features_tabulation)
#     features_tabulation <- features_tabulation[which(tabulations.rowSums != 0),,drop = F]
#     bin_centers <- bin_centers[which(tabulations.rowSums != 0)]
#   }
# 
#   if(!return_bins)  return(t(features_tabulation))
#   else return(list(t(features_tabulation), bin_centers))
# }

coordinate_histograms <- 
  function(cgr_coords, bin_count, return_bins = F, drop_empty = F){
    # Remove columns with constant data
    cgr_coords <- lapply(cgr_coords, 
                         function(x) preProcess(x, method = "zv") |>  predict(x))
    
    # Remove the origin point if it is included
    cgr_coords <- lapply(cgr_coords,
                         function(x){
                           if(all.equal(x[1,,drop = T], c(0, 0, 0), check.attributes = F)) return(x[-1,,drop = F])
                           else return(x)
                         }
    )
    
    # Remove row of zeros if there is one 
    if(all(sapply(cgr_coords, function(x) all(x[1,]==0)))){
      cgr_coords <- lapply(cgr_coords, function(x) x[-1,])
    }
    
    # Get histogram bounds
    hist_lb <- apply(do.call(rbind, cgr_coords), 2, min)
    hist_ub <- apply(do.call(rbind, cgr_coords), 2, max)
    
    # Get histogram intervals for each axis
    hist_breaks <- list()
    for(i in seq_along(hist_lb)){
      insert_me <- unique(seq(hist_lb[i], hist_ub[i], length.out = bin_count + 1))
      if(length(insert_me) < 2) insert_me <- c(insert_me, insert_me + 1) 
      hist_breaks[[i]] <- insert_me
    }
    bin_centers <- sapply(hist_breaks, partial(.f = zoo::rollmean, k = 2))
    
    # Get histogram bin counts
    tabulations <- lapply(cgr_coords, 
                          function(x) hist_df_unpaired(x, breaks = hist_breaks))
    tabulations <- sapply(tabulations, unlist)
    colnames(tabulations) <- names(cgr_coords)
    
    # Remove bins where it's zero across all organisms 
    if(drop_empty){
      tabulations.rowSums <- rowSums(tabulations)
      tabulations <- tabulations[which(tabulations.rowSums != 0),, drop = F]
      bin_centers <- bin_centers[which(tabulations.rowSums != 0),, drop = F]
    }
    
    if(!return_bins)  return(t(tabulations)) 
    else return(list(t(tabulations), bin_centers))
  }

feature_signature <- function(dna_features, bin_count){
  features_breaks <- seq(min(unlist(dna_features)), max(unlist(dna_features)),
                         length.out = bin_count + 1)

  # Get histogram bin counts
  features_tabulation <- sapply(dna_features,
                                function(x) hist(x, breaks = features_breaks, plot = F)$counts)

  # Calculate z-scores within organisms
  features_tabulation <-
    preProcess(features_tabulation, method = c("center", "scale")) |>
    predict(features_tabulation)

  return(t(features_tabulation))
}

cgr_distance <- function(seq_cg, frac = 1/15, cs = T, ...){
  combine_distances <- function(dist_list, p = 1){
    n <- length(dist_list)
    Reduce(`+`, lapply(dist_list, function(x){(1/n) * x^p}))^(1/p)
  }
  center_scale <- function(hist_mat,force_pos  = F){
    output <- 
      if(force_pos) {
        output <- apply(hist_mat, 2, function(x) {(x- mean(x)) / sd(x)})
        output <- output + abs(min(output))
      } else {
        output <- apply(hist_mat, 2, function(x) {(x- mean(x)) / sd(x)})
      }
    return(output)
  }
  if(cs) cen_scal <- center_scale
  else cen_scal <- identity
  
  bins <- ceiling(frac * (mean(sapply(seq_cg, nrow)) - 1))
  
  combine_distances(list(
    feature_histograms(lapply(seq_cg, function(x) orientedangle1(x[-1,], v =c(1, 0, 0))), bin_count = bins) |> cen_scal() |> dist() |> as.matrix(),
    feature_histograms(lapply(seq_cg, function(x) orientedangle1(x[-1,], v =c(0, 1, 0))), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, function(x) orientedangle1(x[-1,], v =c(0, 0, 1))), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, by3roworientedangle4, v =c(1, 0, 0)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, by3roworientedangle4, v =c(0, 1, 0)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, by3roworientedangle4, v =c(0, 0, 1)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, orienteddistance1, v =c(1, 0, 0)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, orienteddistance1, v =c(0, 1, 0)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, orienteddistance1, v =c(0, 0, 1)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, by3roworienteddistance4, v =c(1, 0, 0)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, by3roworienteddistance4, v =c(0, 1, 0)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    feature_histograms(lapply(seq_cg, by3roworienteddistance4, v =c(0, 0, 1)), bin_count = bins) |> cen_scal() |> dist()|> as.matrix(),
    (coordinate_histograms(seq_cg,  bin_count = bins) |> split_coord_histograms()) [[1]]|> cen_scal() |>  dist()|> as.matrix(),
    (coordinate_histograms(seq_cg,  bin_count = bins) |> split_coord_histograms()) [[2]]|> cen_scal() |>  dist()|> as.matrix(),
    (coordinate_histograms(seq_cg,  bin_count = bins) |> split_coord_histograms()) [[3]]|> cen_scal()|>  dist()|> as.matrix())
  ) 
}

#=====================================================================
# Function to calculate distances between histograms 
#=====================================================================

kdist <- function(input_hist) {
  # Input is a histogram with each row associated with an organism 
  # and and each column corresponds to a histogram bin
  n_organisms <- nrow(input_hist)
  seq_lengths <- rowSums(input_hist)
  Fij <- matrix(0, nrow = n_organisms, ncol = n_organisms)
  rownames(Fij) <- rownames(input_hist)
  colnames(Fij) <- rownames(input_hist)
  res <- matrix(0, nrow = n_organisms, ncol = n_organisms)
  a <- log(1.1)
  b <- log(0.1) - a
  for (i in 1:(n_organisms-1)) {
    for (j in (i + 1):n_organisms) {
      denom <- min(seq_lengths[i], seq_lengths[j])
      Fij[i, j] <- sum(apply(input_hist[c(i, j),], 2, min) / denom)
    }
  }
  Fij <- Fij + t(Fij)
  res <- (log(0.1 + Fij) - a) / b
  res[which(res < 0)] <- 0
  diag(Fij) <- 1
  diag(res) <- 0
  return(res)
}


kdist2 <- function(input_hist) {
  # Input is a histogram with each row associated with an organism 
  # and and each column corresponds to a histogram bin
  n_organisms <- nrow(input_hist)
  seq_lengths <- rowSums(input_hist)
  Fij <- matrix(0, nrow = n_organisms, ncol = n_organisms)
  rownames(Fij) <- rownames(input_hist)
  colnames(Fij) <- rownames(input_hist)
  res <- matrix(0, nrow = n_organisms, ncol = n_organisms)
  for (i in 1:(n_organisms-1)) {
    for (j in (i + 1):n_organisms) {
      denom <- abs(seq_lengths[i]/ncol(input_hist)-seq_lengths[j]/ncol(input_hist))
      Fij[i, j] <- sum(apply(input_hist[c(i, j),], 2, min) / max(denom,   min(seq_lengths)/ (2*ncol(input_hist))))
    }
  }
  Fij <- Fij + t(Fij)
  a <- max(Fij + .1) 
  b <- log(0.1) - a
  res <- (log(0.1 + Fij) - a) / b
  res[which(res < 0)] <- 0
  diag(res) <- 0
  return(res)
}

kdist5 <- function(input_hist) {
  # Input is a histogram with each row associated with an organism 
  # and and each column corresponds to a histogram bin
  n_organisms <- nrow(input_hist)
  seq_lengths <- rowSums(input_hist)
  Fij <- matrix(0, nrow = n_organisms, ncol = n_organisms)
  rownames(Fij) <- rownames(input_hist)
  colnames(Fij) <- rownames(input_hist)
  res <- matrix(0, nrow = n_organisms, ncol = n_organisms)
  a <- log(1.1)
  b <- log(0.1) - a
  for (i in 1:(n_organisms-1)) {
    for (j in (i + 1):n_organisms) {
      denom <- max(seq_lengths[i], seq_lengths[j])
      og_numerator <- apply(input_hist[c(i, j),], 2, min)
      og_numerator[which(og_numerator == 0)] <- abs(seq_lengths[i] - seq_lengths[j])
      Fij[i, j] <- sum(apply(input_hist[c(i, j),], 2, min) / denom)
    }
  }
  Fij <- Fij + t(Fij)
  res <- (log(0.1 + Fij) - a) / b
  res[which(res < 0)] <- 0
  diag(res) <- 0
  return(res)
}

#=====================================================================
# Function to implement coordinate signature method
#=====================================================================

coordinate_signature <- function(cgr_coords, bin_count){
  # Remove columns with constant data
  cgr_coords <- lapply(cgr_coords, 
                       function(x) preProcess(x, method = "zv") |>  predict(x))
  
  # Remove row of zeros if there is one 
  if(all(sapply(cgr_coords, function(x) all(x[1,]==0)))){
    cgr_coords <- lapply(cgr_coords, function(x) x[-1,])
  }
  
  # Get histogram bounds
  hist_lb <- apply(do.call(rbind, cgr_coords), 2, min)
  hist_ub <- apply(do.call(rbind, cgr_coords), 2, max)
  
  # Get histogram intervals for each axis
  hist_breaks <- list()
  for(i in seq_along(hist_lb)){
    insert_me <- unique(seq(hist_lb[i], hist_ub[i], length.out = bin_count + 1))
    if(length(insert_me) < 2) insert_me <- c(insert_me, insert_me + 1) 
    hist_breaks[[i]] <- insert_me
  }
  # Get histogram bin counts
  tabulations <- lapply(cgr_coords, 
                        function(x) hist_df_unpaired(x, breaks = hist_breaks))
  tabulations <- sapply(tabulations, unlist)
  colnames(tabulations) <- names(cgr_coords)
  
  # Calculate z-scores within organisms 
  tabulations <- preProcess(tabulations, method = c("center","scale")) |> 
    predict(tabulations)
  
  return(t(tabulations))
}

#=====================================================================
# Function to implement volume intersection method
#=====================================================================

volume_intersection_tanimoto <- function(sequence_list, bandwidth = 0.003, hv_args = list(), vi_args = list()){
  output <- matrix(1, nrow = length(sequence_list), ncol = length(sequence_list))
  
  cg4_list <- lapply(sequence_list, function(dna_seq) seq_to_hypercomplex_cg4(dna_seq) |> (\(x)(x[-1,-1]))())  
  
  hypervolumes1 <- 
    lapply(seq_along(sequence_list), 
           function(i) 
             do.call(hypervolume_gaussian, 
                     c(list(cg4_list[[i]], sd.count = 4,  
                            kde.bandwidth = estimate_bandwidth(data=cg4_list[[i]], method = "fixed", value = bandwidth),
                            name = names(sequence_list)[[i]]), hv_args)))
  
  pairwise_combn <- combn(1:length(sequence_list), 2, simplify = F)
  
  output <- 
    lapply(pairwise_combn, 
           function(i){
             test1 <- hypervolumes1[[i[1]]]
             test2 <- hypervolumes1[[i[2]]]
             check <- do.call(hypervolume_set, c(list(hv1 = test1, hv2 = test2, check.memory = F), vi_args))
             tanimoto <- check@HVList$Intersection@Volume / 
               (test1@Volume + test2@Volume - check@HVList$Intersection@Volume)
             data.frame(Var1 = c(i), 
                        Var2 = c(rev(i)), 
                        tanimoto = rep(tanimoto, 2))})
  
  output <- do.call(rbind, output)
  output <- 
    rbind(output, data.frame(Var1 = 1:length(sequence_list),
                             Var2 = 1:length(sequence_list), 
                             tanimoto = rep(1, length(sequence_list)))) |>
    arrange(Var1, Var2)
  output <- output |> 
    pivot_wider(id_cols = Var1, names_from = Var2, values_from = tanimoto) |>
    column_to_rownames("Var1") |> 
    as.matrix()
  
  colnames(output) <- names(sequence_list)
  rownames(output) <- names(sequence_list)
  return(output)
}


#=====================================================================
# Plotting functions
#=====================================================================

