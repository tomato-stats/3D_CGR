#####################################################################
# Functions written to implement 3D CGR for DNA sequences 
# By: Stephanie Young <syoung2@sdsu.edu> <syoung49@its.jnj.com>
#####################################################################

library(stringr)
library(hypervolume)
library(parallel)

#=====================================================================
# Convert a DNA sequence to 3D CGR 
#=====================================================================

# Coordinate table
## Function to make a table containing 3DCGR coordinates
CGR_table <- function(A_, T_, C_, G_){
  output <- rbind(A_, T_, C_, G_) |> as.data.frame() 
  colnames(output) <- c("i", "j", "k")
  rownames(output)
  output[["Nucleotide"]] <- c("A", "T", "C", "G")
  output
}

coord1 <- CGR_table(
  (c(0, 0, 0) - (1/2)) * 2*sqrt(1/3),
  (c(1, 1, 0) - (1/2)) * 2*sqrt(1/3),
  (c(1, 0, 1) - (1/2)) * 2*sqrt(1/3),
  (c(0, 1, 1) - (1/2)) * 2*sqrt(1/3)
)

coord2 <- CGR_table(
  c(0, 0, 1),
  c(2*sqrt(2) / 3, 0, -1/3),
  c(-sqrt(2)/3, sqrt(6)/3, -1/3),
  c(-sqrt(2)/3, -sqrt(6)/3, -1/3)
)

# (recursive implementation)
seq_to_hypercomplex_cg4 <- function(dna_seq){ 
  # The input to this function is the DNA sequence and 
  # a table of the CGR coordinates to be used
  if(length(dna_seq) == 1) dna_seq <- str_split(dna_seq, "")[[1]]
  A <- (c(0, 0, 0) - (1/2)) * 2*sqrt(1/3)
  T <- (c(1, 1, 0) - (1/2)) * 2*sqrt(1/3) 
  C <- (c(1, 0, 1) - (1/2)) * 2*sqrt(1/3)
  G <- (c(0, 1, 1) - (1/2)) * 2*sqrt(1/3)
  cg <- data.frame(i = 0,
                   j = 0, 
                   k = 0)
  for(i in seq_along(dna_seq)){
    last_obs <- cg[nrow(cg), ]
    letter_values <- get(dna_seq[i])
    cg <- rbind(cg, 
                data.frame(i = mean(c(last_obs[[1]], letter_values[1])),
                           j = mean(c(last_obs[[2]], letter_values[2])),
                           k = mean(c(last_obs[[3]], letter_values[3])))
    )
  }
  cg[["r"]] <- 0
  cg <- cg[, c("r", "i", "j", "k")]
  return(cg)
}

# (non-recursive implementation)
seq_to_hypercomplex_cg4_nr <- function(dna_seq, CGR_coord = coord1){
  # The input to this function is the DNA sequence and 
  # a table of the CGR coordinates to be used
  if(length(dna_seq) == 1) dna_seq <- str_split(dna_seq, "")[[1]]
  temp <- CGR_coord[match(dna_seq, CGR_coord$Nucleotide), ]
  r <- 1/2
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, c("temp"), envir = environment())

  cg <- 
    parLapply(
      cl,
      1:nrow(temp), 
      function(j){
        k <- 1:j
        rowSums(r * rep((1-r)^(k-1), each = 3) * t(temp[j-k+1,1:3]))
      }
    )
  cg <- do.call(rbind, cg) |> as.data.frame()
  cg[["r"]] <- 0
  cg <- cg[, c("r", "i", "j", "k")]
  cg <- rbind(data.frame(r = 0, i = 0, j = 0, k = 0), cg)
  return(cg)
}


#=====================================================================
# Functions necessary for shape signature methods
#=====================================================================

## Functions related to angular signature

pt_dist <- function(p1, p2){
  sqrt(sum((p1 - p2)^2))
}

angle_between_3pts <- function(p1, pv, p3){
  pv1 <- pt_dist(pv, p1)
  pv3 <- pt_dist(pv, p3)
  p13 <- pt_dist(p1, p3)
  acos_input <- round((pv1^2 + pv3^2 - p13^2) / (2 * pv1 * pv3), digits = 10)
  acos(acos_input)
}

by3rowangle <- function(df){
  i <- 3:nrow(df)
  angles <- sapply(i, function(x) angle_between_3pts(df[x,], df[x-1,], df[x-2,]))
}

## Functions related to edge signature

by3rowdistance <- function(df){
  i <- 3:nrow(df)
  angles <- sapply(i, function(x) pt_dist(df[x,], df[x-2,]))
}

## Functions related to coordinate signature 

hist_df_unpaired <- function(df, breaks){
  if(length(breaks) != ncol(df)) "hist_df_unpaired being used incorrectly"
  sapply(1:ncol(df), function(x) hist(df[,x], breaks[[x]], plot = F)$counts)
}

#=====================================================================
# Function to implement feature signature methods 
# (applicable to angles and edges)
#=====================================================================
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

#=====================================================================
# Function to implement coordinate signature method
#=====================================================================

coordinate_signature <- function(cgr_coords, bin_count){
  cgr_coords <- lapply(cgr_coords, 
                       function(x) preProcess(x, method = "zv") |>  predict(x))
  
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
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, 
                c("str_split", 
                  "seq_to_hypercomplex_cg4", "seq_to_hypercomplex_cg4_nr"), 
                envir = environment())
  
  cg4_list <- parLapply(cl, sequence_list, function(dna_seq) seq_to_hypercomplex_cg4_nr(dna_seq) |> (\(x)(x[-1,-1]))())  
  
  clusterEvalQ(cl, {
    ## set up each worker.
    library(hypervolume)
  })
  
  clusterExport(cl, 
                c("sequence_list", "cg4_list", "bandwidth"),
                envir = environment())
  
  hypervolumes1 <- 
    parLapply(cl, seq_along(sequence_list), 
              function(i) 
                do.call(hypervolume_gaussian, c(list(cg4_list[[i]], sd.count = 4,  
                                                     kde.bandwidth = estimate_bandwidth(data=cg4_list[[i]], method = "fixed", value = bandwidth),
                                                     name = names(sequence_list)[[i]]), hv_args)))
  
  gc()
  pairwise_combn <- combn(1:length(sequence_list), 2, simplify = F)
  
  
  output <- parLapply(cl, pairwise_combn, 
                      function(i){
                        test1 <- hypervolumes1[[i[1]]]
                        test2 <- hypervolumes1[[i[2]]]
                        check <- do.call(hypervolume_set, c(list(hv1 = test1, hv2 = test2, check.memory = F), vi_args))
                        tanimoto <- check@HVList$Intersection@Volume / 
                          (test1@Volume + test2@Volume - check@HVList$Intersection@Volume)
                        data.frame(Var1 = c(i), 
                                   Var2 = c(rev(i)), 
                                   tanimoto = rep(tanimoto, 2))})
  
  stopCluster(cl)
  output <- do.call(rbind, output)
  output <- rbind(output, data.frame(Var1 = 1:length(sequence_list),
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
