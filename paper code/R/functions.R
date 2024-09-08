#####################################################################
# Functions written to implement 3D CGR for DNA sequences 
# By: Stephanie Young <syoung2@sdsu.edu> <syoung49@its.jnj.com>
#####################################################################

library(stringr)
library(zoo)
library(hypervolume)
library(Rcpp)
library(parallel)

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
seq_to_cgr <- 
  function(dna_seq, CGR_coord = coord1, df = F, axes = c("i", "j", "k")){ 
    # The input to this function is the DNA sequence and 
    # a table of the CGR coordinates to be used
    if(ncol(CGR_coord) != (length(axes) + 1))
      stop("Number of axes needs to match that of the CGR_coord") 
    
    dna_seq <- toupper(dna_seq)
    if(length(dna_seq) == 1) dna_seq <- str_split(dna_seq, "")[[1]]
    
    # Remove gaps
    dna_seq <- dna_seq[which(dna_seq != "-")]
    
    # Chaos game 
    cg <- chaosGame(dna_seq, as.matrix(CGR_coord[,1:length(axes)]))
    colnames(cg) <- axes
    return(cg)
  }

#=====================================================================
# Functions for assist in building histograms
#=====================================================================

sourceCpp("./R/features.cpp")

## Functions for building histograms for multiple features.
## Paired histogram currently not used, but may be helpful for future 
## implementations using joint probability distribution information 

### One histogram for each feature
hist_df_unpaired <- 
  function(df, breaks){
    if(length(breaks) != ncol(df)) "hist_df_unpaired being used incorrectly"
    sapply(1:ncol(df), function(x) hist(df[,x], breaks[[x]], plot = F)$counts)
  }

### Functions for building joint histograms for multiple features
cutting_df <- 
  function(df, breaks){
    # df is a data frame with n columns 
    # breaks is a data frame with n columns
    output <- data.frame(matrix(ncol = 0, nrow = nrow(df)))
    for(i in 1:ncol(df)){
      output[[paste0("interval_", i)]] <- cut(df[,i], breaks[,i], include.lowest = T)
    }
    
    output |> group_by(across(starts_with("interval_"))) |> 
      summarise(n = n(), .groups = "drop")
  }

hist_paired <- 
  function(list_df, bin_count, return_bins = F){
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

zv <- function(input_matrix, dim = 2){
  # identify zero-variance columns
  zv <- apply(input_matrix, dim, function(x){ length(unique(x)) == 1})
  return(zv)
}

common_function <- function(list_of_vectors, FUNC, ...){
  # A faster way to find a common min/max across all lists of vectors
  output <- sapply(list_of_vectors, FUNC, ...)
  output <- FUNC(output)
  output
}

feature_histograms <- 
  function(features, bin_count, return_breaks = F, return_bin_centers = F){
    if(common_function(features, max) == pi | common_function(features, min) == -pi){
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
      features_breaks <- seq(common_function(not_pi, min), common_function(not_pi, max), length.out = bin_count - 1)
      features_breaks <- c(-pi, features_breaks, pi)
    } else {
      features_breaks <- seq(common_function(features, min), common_function(features, max), length.out = bin_count + 1)
    }
    # Get histogram bin counts (parameters set to replicate hist() function)
    features_tabulation <-
      sapply(
        features,
        # function(x) hist(x, breaks = features_breaks, plot = F)$counts # This version is slower
        function(x) tabulate(findInterval(x, vec = features_breaks, rightmost.closed = T, all.inside = T, left.open = T), bin_count)
      )
    bin_centers <- zoo::rollmean(features_breaks, k = 2)
    
    output <- t(features_tabulation)
    if(return_breaks) output <- list(output, features_breaks)
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
  function(cgr_coords, bin_count, return_breaks = F, return_bin_centers = F){
    
    # Remove columns with constant data
    const <- do.call("rbind", cgr_coords) |> zv()
    cgr_coords <- lapply(cgr_coords, function(x) x[, !const, drop = F])
    
    # Remove the origin if it is included
    cgr_coords <- 
      lapply(cgr_coords,
             function(x){
               if(all(x[1,] == 0)) return(x[-1, , drop = F])
               else return(x)
             }
      )
    
    # Get histogram bounds
    ncols =  dim(cgr_coords[[1]])[2]
    hist_lb <- vector(length = ncols, mode = "numeric")
    hist_ub <- vector(length = ncols, mode = "numeric")
    for(col in 1:ncols){
      hist_lb[col] <- common_function(lapply(cgr_coords, function(x) x[,col]), min)
      hist_ub[col] <- common_function(lapply(cgr_coords, function(x) x[,col]), max)
    }
    
    # Get histogram intervals for each axis
    hist_breaks <- vector(mode = "list", length = length(hist_lb))
    for(i in seq_along(hist_lb)){
      breaks_i <- unique(seq(hist_lb[i], hist_ub[i], length.out = bin_count + 1))
      if(length(breaks_i) < 2) breaks_i <- c(breaks_i, breaks_i + 1) 
      hist_breaks[[i]] <- breaks_i
    }
    bin_centers <- sapply(hist_breaks, partial(.f = zoo::rollmean, k = 2))
    
    # Get histogram bin counts
    tabulations <- lapply(cgr_coords, 
                          function(x) hist_df_unpaired(x, breaks = hist_breaks))
    tabulations <- sapply(tabulations, unlist)
    colnames(tabulations) <- names(cgr_coords)

    output <- t(tabulations)
    if(return_breaks) output <- list(output, hist_breaks)
    if(return_bin_centers) output <- list(output, bin_centers)
    return(output)
  }

# feature_signature <- function(dna_features, bin_count){
#   features_breaks <- seq(min(unlist(dna_features)), max(unlist(dna_features)),
#                          length.out = bin_count + 1)
# 
#   # Get histogram bin counts
#   features_tabulation <- sapply(dna_features,
#                                 function(x) hist(x, breaks = features_breaks, plot = F)$counts)
# 
#   # Calculate z-scores within organisms
#   features_tabulation <-
#     preProcess(features_tabulation, method = c("center", "scale")) |>
#     predict(features_tabulation)
# 
#   return(t(features_tabulation))
# }

# #=====================================================================
# # Function to implement coordinate signature method
# #=====================================================================
# 
# coordinate_signature <- function(cgr_coords, bin_count){
#   # Remove columns with constant data
#   cgr_coords <- lapply(cgr_coords, 
#                        function(x) preProcess(x, method = "zv") |>  predict(x))
#   
#   # Remove row of zeros if there is one 
#   if(all(sapply(cgr_coords, function(x) all(x[1,]==0)))){
#     cgr_coords <- lapply(cgr_coords, function(x) x[-1,])
#   }
#   
#   # Get histogram bounds
#   hist_lb <- apply(do.call(rbind, cgr_coords), 2, min)
#   hist_ub <- apply(do.call(rbind, cgr_coords), 2, max)
#   
#   # Get histogram intervals for each axis
#   hist_breaks <- list()
#   for(i in seq_along(hist_lb)){
#     insert_me <- unique(seq(hist_lb[i], hist_ub[i], length.out = bin_count + 1))
#     if(length(insert_me) < 2) insert_me <- c(insert_me, insert_me + 1) 
#     hist_breaks[[i]] <- insert_me
#   }
#   # Get histogram bin counts
#   tabulations <- lapply(cgr_coords, 
#                         function(x) hist_df_unpaired(x, breaks = hist_breaks))
#   tabulations <- sapply(tabulations, unlist)
#   colnames(tabulations) <- names(cgr_coords)
#   
#   # Calculate z-scores within organisms 
#   tabulations <- preProcess(tabulations, method = c("center","scale")) |> 
#     predict(tabulations)
#   
#   return(t(tabulations))
# }

#=====================================================================
# Functions to calculate distance based on shape signature
#=====================================================================

split_coord_histograms <- 
  function(histograms, num_axes){
    histograms.ncol <- ncol(histograms)
    starts <- (0:(num_axes - 1)) * (histograms.ncol / num_axes) + 1
    ends <-  (1:num_axes) * (histograms.ncol / num_axes)
    output <- lapply(1:num_axes, function(i) histograms[,starts[i]:ends[i]])
    return(output) 
  }

cgr_distance <- 
  function(seq_cg, frac = 1/15, cs = T, power = 1, normalize = F){
    
    combine_distances <- function(dist_list, p){
      # Power mean: arithmetic p = 1, geometric p = 0, harmonic p = -1, quadratic p = 2
      # As p -> infinity, power mean -> max; 
      # as p -> -infinity, power mean -> min. 
      n <- length(dist_list)
      Reduce(`+`, lapply(dist_list, function(x){(1/n) * x^p}))^(1/p)
    }
    
    center_scale <- function(hist_mat){
      apply(hist_mat, 2, function(x) {(x- mean(x)) / sd(x)}) 
    }
    
    cen_scal <- identity
    norm <- identity
    if(cs) cen_scal <- center_scale
    if(normalize) norm <- function(x){x / max(x)}
    
    bins <- ceiling(frac * (mean(sapply(seq_cg, nrow)) - 1))
    
    shape_features <- list(
      lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(1, 0, 0))),
      lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(0, 1, 0))),
      lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(0, 0, 1))),
      lapply(seq_cg, by3roworientedangle4, v = c(1, 0, 0)),
      lapply(seq_cg, by3roworientedangle4, v = c(0, 1, 0)),
      lapply(seq_cg, by3roworientedangle4, v = c(0, 0, 1)),
      lapply(seq_cg, orienteddistance1, v = c(1, 0, 0)),
      lapply(seq_cg, orienteddistance1, v = c(0, 1, 0)),
      lapply(seq_cg, orienteddistance1, v = c(0, 0, 1)),
      lapply(seq_cg, by3roworienteddistance4, v = c(1, 0, 0)),
      lapply(seq_cg, by3roworienteddistance4, v = c(0, 1, 0)),
      lapply(seq_cg, by3roworienteddistance4, v = c(0, 0, 1))
    ) |> 
      lapply(feature_histograms, bin_count = bins) |> 
      lapply(function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
    
    coord_features <- 
      coordinate_histograms(seq_cg,  bin_count = bins) |> 
      split_coord_histograms(ncol(seq_cg[[1]])) |>
      lapply(function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
    
    combine_distances(
      c(shape_features, coord_features), 
      p = power
    )
  }

cgr_distance2 <- 
  function(seq_cg, frac = 1/15, cs = T, power = 1, normalize = F){
    
    combine_distances <- function(dist_list, p){
      # Power mean: arithmetic p = 1, geometric p = 0, harmonic p = -1, quadratic p = 2
      # As p -> infinity, power mean -> max; 
      # as p -> -infinity, power mean -> min. 
      n <- length(dist_list)
      Reduce(`+`, lapply(dist_list, function(x){(1/n) * x^p}))^(1/p)
    }
    
    center_scale <- function(hist_mat){
      apply(hist_mat, 1, function(x) {(x- mean(x)) / sd(x)}) |> t()
    }
    
    cen_scal <- identity
    norm <- identity
    if(cs) cen_scal <- center_scale
    if(normalize) norm <- function(x){x / max(x)}
    
    bins <- ceiling(frac * (mean(sapply(seq_cg, nrow)) - 1))
    
    reference_coordinates <- list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
    feature_functions <- list(function(x, ...) orientedangle1(x[-1,], ...),
                              by3roworientedangle4,
                              orienteddistance1,
                              by3roworienteddistance4)
    shape_histograms <- vector(mode = "list", length = 12)
    i <- 1
    for(f in seq_along(feature_functions)){
      for(v in seq_along(reference_coordinates)){
        shape_histograms[[i]] <- lapply(seq_cg, feature_functions[[f]], v = reference_coordinates[[v]])
        shape_histograms[[i]] <- feature_histograms(shape_histograms[[i]] , bin_count = bins)
        shape_histograms[[i]] <- shape_histograms[[i]]  |> cen_scal() |> dist() |> as.matrix() |> norm()
        i <- i + 1
      }
    }
    
    coord_features <- 
      coordinate_histograms(seq_cg,  bin_count = bins) |> 
      split_coord_histograms(ncol(seq_cg[[1]])) 
    coord_features <- 
      lapply(coord_features, function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
    combine_distances(
      c(shape_histograms, coord_features), 
      p = power
    )
  }

cgr_distance3<- 
  function(seq_cg, frac = 1/15, cs = T, power = 1, normalize = F){
    
    combine_distances <- function(dist_list, p){
      # Power mean: arithmetic p = 1, geometric p = 0, harmonic p = -1, quadratic p = 2
      # As p -> infinity, power mean -> max; 
      # as p -> -infinity, power mean -> min. 
      n <- length(dist_list)
      Reduce(`+`, lapply(dist_list, function(x){(1/n) * x^p}))^(1/p)
    }
    
    center_scale <- function(hist_mat){
      sc1 <- apply(hist_mat, 1, function(x) {(x- mean(x)) / sd(x)}) |> t()
      apply(sc1, 2, function(x) {(x- mean(x)) / sd(x)})
    }
    
    cen_scal <- identity
    norm <- identity
    if(cs) cen_scal <- center_scale
    if(normalize) norm <- function(x){x / max(x)}
    
    bins <- ceiling(frac * (mean(sapply(seq_cg, nrow)) - 1))
    
    shape_features <- list(
      lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(1, 0, 0))),
      lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(0, 1, 0))),
      lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(0, 0, 1))),
      lapply(seq_cg, by3roworientedangle4, v = c(1, 0, 0)),
      lapply(seq_cg, by3roworientedangle4, v = c(0, 1, 0)),
      lapply(seq_cg, by3roworientedangle4, v = c(0, 0, 1)),
      lapply(seq_cg, orienteddistance1, v = c(1, 0, 0)),
      lapply(seq_cg, orienteddistance1, v = c(0, 1, 0)),
      lapply(seq_cg, orienteddistance1, v = c(0, 0, 1)),
      lapply(seq_cg, by3roworienteddistance4, v = c(1, 0, 0)),
      lapply(seq_cg, by3roworienteddistance4, v = c(0, 1, 0)),
      lapply(seq_cg, by3roworienteddistance4, v = c(0, 0, 1))
    ) |> 
      lapply(feature_histograms, bin_count = bins) |> 
      lapply(function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
    
    coord_features <- 
      coordinate_histograms(seq_cg,  bin_count = bins) |> 
      split_coord_histograms(ncol(seq_cg[[1]])) |>
      lapply(function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
    
    combine_distances(
      c(shape_features, coord_features), 
      p = power
    )
  }

cgr_distance4 <- function(seq_cg, frac = 1/15, cs = T, power = 1, normalize = F){
  
  combine_distances <- function(dist_list, p){
    # Power mean: arithmetic p = 1, geometric p = 0, harmonic p = -1, quadratic p = 2
    # As p -> infinity, power mean -> max; 
    # as p -> -infinity, power mean -> min. 
    n <- length(dist_list)
    Reduce(`+`, lapply(dist_list, function(x){(1/n) * x^p}))^(1/p)
  }
  
  center_scale <- function(hist_mat){
    apply(hist_mat, 2, function(x) {(x- mean(x)) / sd(x)}) 
  }
  
  cen_scal <- identity
  norm <- identity
  if(cs) cen_scal <- center_scale
  if(normalize) norm <- function(x){x / max(x)}
  
  bins <- ceiling(frac * (mean(sapply(seq_cg, nrow)) - 1))
  
  shape_features <- list(
    lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(1, 0, 0))),
    lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(0, 1, 0))),
    lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(0, 0, 1))),
    lapply(seq_cg, by3roworientedangle4, v = c(1, 0, 0)),
    lapply(seq_cg, by3roworientedangle4, v = c(0, 1, 0)),
    lapply(seq_cg, by3roworientedangle4, v = c(0, 0, 1)),
    lapply(seq_cg, orienteddistance1, v = c(1, 0, 0)),
    lapply(seq_cg, orienteddistance1, v = c(0, 1, 0)),
    lapply(seq_cg, orienteddistance1, v = c(0, 0, 1)),
    lapply(seq_cg, by3roworienteddistance4, v = c(1, 0, 0)),
    lapply(seq_cg, by3roworienteddistance4, v = c(0, 1, 0)),
    lapply(seq_cg, by3roworienteddistance4, v = c(0, 0, 1))
  ) |> 
    lapply(feature_histograms, bin_count = bins) |> 
    lapply(function(x) t(apply(x, 1, function(x) x / sum(x)))) |>
    lapply(function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
  
  coord_features <- 
    coordinate_histograms(seq_cg,  bin_count = bins) |> 
    split_coord_histograms(ncol(seq_cg[[1]])) |>
    lapply(function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
  
  combine_distances(
    c(shape_features, coord_features), 
    p = power
  )
}

cgr_distance5 <- function(seq_cg, frac = 1/15, cs = T, power = 1, normalize = F){
  
  combine_distances <- function(dist_list, p){
    # Power mean: arithmetic p = 1, geometric p = 0, harmonic p = -1, quadratic p = 2
    # As p -> infinity, power mean -> max; 
    # as p -> -infinity, power mean -> min. 
    n <- length(dist_list)
    Reduce(`+`, lapply(dist_list, function(x){(1/n) * x^p}))^(1/p)
  }
  
  center_scale <- function(hist_mat){
    sc1 <- apply(hist_mat, 1, function(x) {(x- mean(x)) / sd(x)}) |> t()
    apply(sc1, 2, function(x) {(x- mean(x)) / sd(x)})
  }
  
  cen_scal <- identity
  norm <- identity
  if(cs) cen_scal <- center_scale
  if(normalize) norm <- function(x){x / max(x)}
  
  bins <- ceiling(frac * (mean(sapply(seq_cg, nrow)) - 1))
  
  shape_features <- list(
    lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(1, 0, 0))),
    lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(0, 1, 0))),
    lapply(seq_cg, function(x) orientedangle1(x[-1,], v = c(0, 0, 1))),
    lapply(seq_cg, by3roworientedangle4, v = c(1, 0, 0)),
    lapply(seq_cg, by3roworientedangle4, v = c(0, 1, 0)),
    lapply(seq_cg, by3roworientedangle4, v = c(0, 0, 1)),
    lapply(seq_cg, orienteddistance1, v = c(1, 0, 0)),
    lapply(seq_cg, orienteddistance1, v = c(0, 1, 0)),
    lapply(seq_cg, orienteddistance1, v = c(0, 0, 1)),
    lapply(seq_cg, by3roworienteddistance4, v = c(1, 0, 0)),
    lapply(seq_cg, by3roworienteddistance4, v = c(0, 1, 0)),
    lapply(seq_cg, by3roworienteddistance4, v = c(0, 0, 1))
  ) |> 
    lapply(feature_histograms, bin_count = bins) |> 
    lapply(function(x) t(apply(x, 1, function(x) x / sum(x)))) |>
    lapply(function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
  
  coord_features <- 
    coordinate_histograms(seq_cg,  bin_count = bins) |> 
    split_coord_histograms(ncol(seq_cg[[1]])) |>
    lapply(function(x) x |> cen_scal() |> dist() |> as.matrix() |> norm())
  
  combine_distances(
    c(shape_features, coord_features), 
    p = power
  )
}




#=====================================================================
# Function to implement volume intersection method
#=====================================================================

vol_int_tan <- 
  function(cg_list, bandwidth = 0.003, hv_args = list(), vi_args = list()){
    # This implementation does not use parallelization and is not recommended 
    output <- matrix(1, nrow = length(cg_list), ncol = length(cg_list))

    hypervolumes1 <- 
      lapply(seq_along(cg_list), 
             function(i) 
               do.call(hypervolume_gaussian, 
                       c(list(cg_list[[i]], sd.count = 4,  
                              kde.bandwidth = estimate_bandwidth(data=cg_list[[i]], method = "fixed", value = bandwidth),
                              name = names(cg_list)[[i]]), hv_args)))
    
    pairwise_combn <- combn(1:length(cg_list), 2, simplify = F)
    
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
      rbind(output, data.frame(Var1 = 1:length(cg_list),
                               Var2 = 1:length(cg_list), 
                               tanimoto = rep(1, length(cg_list)))) |>
      arrange(Var1, Var2)
    output <- output |> 
      pivot_wider(id_cols = Var1, names_from = Var2, values_from = tanimoto) |>
      column_to_rownames("Var1") |> 
      as.matrix()
    
    colnames(output) <- names(cg_list)
    rownames(output) <- names(cg_list)
    return(output)
  }


volume_intersection_tanimoto <- 
  function(cg_list, bandwidth = 0.003, hv_args = list(), vi_args = list()){
    output <- matrix(1, nrow = length(cg_list), ncol = length(cg_list))
    
    cl <- makeCluster(min(detectCores(), 4) - 1)

    ## set up each worker.
    clusterEvalQ(cl, {library(hypervolume)})
    
    clusterExport(
      cl, 
      c("cg_list", "bandwidth"),
      envir = environment()
    )
    
    hypervolumes1 <- 
      parLapply(
        cl, 
        seq_along(cg_list), 
        function(i) 
          do.call(hypervolume_gaussian, 
                  c(list(cg_list[[i]], sd.count = 4,  
                         kde.bandwidth = estimate_bandwidth(data=cg_list[[i]], method = "fixed", value = bandwidth),
                         name = names(cg_list)[[i]]), hv_args)
          )
      )
    
    gc()
    pairwise_combn <- combn(1:length(cg_list), 2, simplify = F)
    
    output <-
      parLapply(
        cl, 
        pairwise_combn, 
        function(i){
          test1 <- hypervolumes1[[i[1]]]
          test2 <- hypervolumes1[[i[2]]]
          check <- do.call(hypervolume_set, c(list(hv1 = test1, hv2 = test2, check.memory = F), vi_args))
          tanimoto <- check@HVList$Intersection@Volume / 
            (test1@Volume + test2@Volume - check@HVList$Intersection@Volume)
          data.frame(Var1 = c(i), Var2 = c(rev(i)), tanimoto = rep(tanimoto, 2))
        }
      )
    
    stopCluster(cl)
    output <- do.call(rbind, output)
    output <- rbind(output, data.frame(Var1 = 1:length(cg_list),
                                       Var2 = 1:length(cg_list), 
                                       tanimoto = rep(1, length(cg_list)))) |>
      arrange(Var1, Var2)
    output <- output |> 
      pivot_wider(id_cols = Var1, names_from = Var2, values_from = tanimoto) |>
      column_to_rownames("Var1") |> 
      as.matrix()
    
    colnames(output) <- names(cg_list)
    rownames(output) <- names(cg_list)
    return(output)
  }

#=====================================================================
# Plotting functions
#=====================================================================

