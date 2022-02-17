#####################################################################
# Functions written to implement 3D CGR for DNA sequences 
# By: Stephanie Young <syoung2@sdsu.edu> <syoung49@its.jnj.com>
#####################################################################

library(stringr)
library(hypervolume)

#=====================================================================
# Convert a DNA sequence to 3D CGR 
#=====================================================================

# (recursive implementation)
seq_to_hypercomplex_cg4 <- function(dna_seq){ 
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
seq_to_hypercomplex_cg4_nr <- function(dna_seq){ 
  if(length(dna_seq) == 1) dna_seq <- str_split(dna_seq, "")[[1]]
  A <- (c(0, 0, 0) - (1/2)) * 2*sqrt(1/3) 
  T <- (c(1, 1, 0) - (1/2)) * 2*sqrt(1/3) 
  C <- (c(1, 0, 1) - (1/2)) * 2*sqrt(1/3)
  G <- (c(0, 1, 1) - (1/2)) * 2*sqrt(1/3)
  r <- 1/2
  cg <- data.frame(i = 0,
                   j = 0, 
                   k = 0)
  y <- function(i){
    S_i = dna_seq[i]
    if(S_i == "A") return(A)
    if(S_i == "T") return(T)
    if(S_i == "C") return(C)
    if(S_i == "G") return(G)
    else return("error")
  }
  y.Vectorize <- Vectorize(y, "i")
  i <- seq_along(dna_seq)
  cg <- lapply(i, 
               function(j){ 
                 k <- 1:j
                 rowSums(r * rep((1-r)^(k-1), each = 3) * y.Vectorize(j-k+1))
               }
  )
  cg <- do.call(rbind, cg) |> as.data.frame()
  colnames(cg) <- c("i", "j", "k")
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

volume_intersection_method <- function(sequence_list, bandwidth = 0.003, ...){
  output <- matrix(1, nrow = length(sequence_list), ncol = length(sequence_list))
  
  hypervolumes1 <- seq_to_hypercomplex_cg4(sequence_list[[1]]) |> (\(x)(x[-1,-1]))() |> 
    (\(x)(hypervolume_gaussian(data= x, sd = 4, 
                               kde.bandwidth = estimate_bandwidth(data=x, method = "fixed", value = bandwidth, ...),
                               name = names(sequence_list)[[1]])))()
  
  for(i in 2:length(sequence_list)){
    append_me <- seq_to_hypercomplex_cg4(sequence_list[[i]]) |> (\(x)(x[-1,-1]))() |> 
      (\(x)(hypervolume_gaussian(data= x, sd = 4, 
                                 kde.bandwidth = estimate_bandwidth(data=x, method = "fixed", value = bandwidth, ...),
                                 name = names(sequence_list)[[i]])))()
    hypervolumes1 <- hypervolume_join(hypervolumes1, append_me)
    gc()
  }
  
  for(i in 1:(length(sequence_list)-1)){
    for(j in (i+1):length(sequence_list)){
      test1 <- hypervolumes1[[i]]
      test2 <- hypervolumes1[[j]]
      check <- hypervolume_set(test1, test2, check.memory = F)
      tanimoto <- check@HVList$Intersection@Volume / 
        (test1@Volume + test2@Volume - check@HVList$Intersection@Volume)
      output[i, j] <- tanimoto
      output[j, i] <- tanimoto
      gc()
    }
  }
  colnames(output) <- names(sequence_list)
  rownames(output) <- names(sequence_list)
  par(mai = c(0, 0, 0, 0), lwd = 1)
  print(plot(hclust(as.dist(1-output), method = "complete" ), 
             axes = F, xlab = NULL, ylab = NULL, main = NULL, 
             sub = NULl, ann = F, cex = .75))
  return(1-output)
}



