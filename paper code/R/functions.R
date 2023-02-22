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

# Coordinate table
## Function to make a table containing 3DCGR coordinates
CGR_table <- function(A_, C_, G_, T_){
  output <- rbind(A_, C_, G_, T_) |> as.data.frame() 
  colnames(output) <- c("i", "j", "k")
  rownames(output)
  output[["Nucleotide"]] <- c("A", "C", "G", "T")
  output <- 
    rbind(
      output, 
      data.frame(output[output$Nucleotide %in% c("A", "G"),1:3] |> colMeans() |> t(), Nucleotide ="R"),
      data.frame(output[output$Nucleotide %in% c("C", "T"),1:3] |> colMeans() |> t(), Nucleotide ="Y"),
      data.frame(output[output$Nucleotide %in% c("G", "T"),1:3] |> colMeans() |> t(), Nucleotide ="K"),
      data.frame(output[output$Nucleotide %in% c("A", "C"),1:3] |> colMeans() |> t(), Nucleotide ="M"),
      data.frame(output[output$Nucleotide %in% c("C", "G"),1:3] |> colMeans() |> t(), Nucleotide ="S"),
      data.frame(output[output$Nucleotide %in% c("A", "T"),1:3] |> colMeans() |> t(), Nucleotide ="W"),
      data.frame(output[output$Nucleotide %in% c("C", "G", "T"),1:3] |> colMeans() |> t(), Nucleotide ="B"),
      data.frame(output[output$Nucleotide %in% c("A", "G", "T"),1:3] |> colMeans() |> t(), Nucleotide ="D"),
      data.frame(output[output$Nucleotide %in% c("A", "C", "T"),1:3] |> colMeans() |> t(), Nucleotide ="H"),
      data.frame(output[output$Nucleotide %in% c("A", "C", "G"),1:3] |> colMeans() |> t(), Nucleotide ="V"),
      data.frame(i = 0, j = 0, k = 0, Nucleotide = "N")
    )
  output
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
seq_to_hypercomplex_cg4 <- function(dna_seq, CGR_coord = coord1, df = F){ 
  # The input to this function is the DNA sequence and 
  # a table of the CGR coordinates to be used
  if(length(dna_seq) == 1) dna_seq <- str_split(toupper(dna_seq), "")[[1]]
  
  dna_seq <- dna_seq[which(dna_seq != "-")]
  
  cg <- chaosGame(dna_seq, as.matrix(CGR_coord[,1:3]))
  if(df){
    cg <- data.frame(cg)
    cg[["r"]] <- 0
    
    cg <- cg[, c("r", "i", "j", "k")]
  }
  else colnames(cg) <- c("i", "j", "k")
  return(cg)
}

#=====================================================================
# Functions necessary for shape signature methods
#=====================================================================

sourceCpp("./R/features.cpp")

distance_from_vertices <- function(points, CGR_coord = coord1){
  if(!is.matrix(points)) points <- as.matrix(points)
  if(all(points[1,] == 0)) points <- points[-1,]
  if(all(points[,1] == 0)) points <- points[,-1]
  distance_matrix(points, as.matrix(coord1[1:4,1:3]))
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

feature_histograms <- function(dna_features, bin_count){
  features_breaks <- seq(min(unlist(dna_features)), max(unlist(dna_features)), 
                         length.out = bin_count + 1)
  
  # Get histogram bin counts
  features_tabulation <- sapply(dna_features, 
                                function(x) hist(x, breaks = features_breaks, plot = F)$counts)
  
  return(t(features_tabulation))
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

