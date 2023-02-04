#####################################################################
# Example DNA sequences analyzed by 3D CGR / paper reproduction
# This code does not run in parallel and can take a while to run.
# By: Stephanie Young <syoung2@sdsu.edu> <syoung49@its.jnj.com>
#####################################################################

#=====================================================================
# Loading data and packages -- 
#   assumes "3D_CGR/paper code" is the current working directory
#=====================================================================

source("./data/sequences.R")
source("../R/functions.R")

library(reshape2)
library(msa)
library(seqinr)
library(plotly)
library(stringr)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(caret)
library(parallel)

set.seed(0)

#=====================================================================
# Clustal Omega results for comparison
#=====================================================================

clustal_omega_tree <- function(dna_seq){
  msa(dna_seq, type = "dna", method = "ClustalOmega") |> 
    msaConvert(type = "seqinr::alignment") |> 
    dist.alignment(gap = T) |> 
    as.matrix() |> as.dist() |> 
    hclust(method = "complete") |>
    plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)
}

clustal_omega_tree(beta_seq)
clustal_omega_tree(nadh)
clustal_omega_tree(sim_fas)

#=====================================================================
# Shape signature results 
#=====================================================================

  
beta_cg <- lapply(beta_seq, seq_to_hypercomplex_cg4)  
nadh_cg <- lapply(nadh, seq_to_hypercomplex_cg4)  
sim_cg <- lapply(sim_fas, seq_to_hypercomplex_cg4)  
  
# Angular signature

beta_features <- lapply(beta_cg, by3rowangle)
nadh_features <- lapply(nadh_cg, by3rowangle)
sim_features <- lapply(sim_cg, by3rowangle)

feature_signature(beta_features, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)
feature_signature(nadh_features, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)
feature_signature(sim_features, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)

# Edge signature

beta_features <- lapply(beta_cg, by3rowdistance)
nadh_features <- lapply(nadh_cg, by3rowdistance)
sim_features <- lapply(sim_cg, by3rowdistance)

feature_signature(beta_features, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)
feature_signature(nadh_features, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)
feature_signature(sim_features, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)

# Coordinate signature

coordinate_signature(beta_cg, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)
coordinate_signature(nadh_cg, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)
coordinate_signature(sim_cg, bin_count = 6110) |>
  dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)

#=====================================================================
# Volume intersection method results 
#=====================================================================

#This part of the code takes the longest to run. 
# Reduce samples.per.point and num.points.max to improve runtime or eliminate
# the parameters altogether to use default values of the hypervolume package
volume_intersection_tanimoto(beta_seq, 
                             hv_args = list(samples.per.point = 5000), 
                             vi_args = list(num.points.max = 150000)) %>% 
  print(plot(hclust(as.dist(1-.), method = "complete" ), 
           axes = F, xlab = NULL, ylab = NULL, main = NULL, 
           sub = NULl, ann = F, cex = .75))

volume_intersection_tanimoto(nadh, 
                             hv_args = list(samples.per.point = 5000), 
                             vi_args = list(num.points.max = 150000)) %>% 
  print(plot(hclust(as.dist(1-.), method = "complete" ), 
             axes = F, xlab = NULL, ylab = NULL, main = NULL, 
             sub = NULl, ann = F, cex = .75))

volume_intersection_tanimoto(sim_fas, 
                             hv_args = list(samples.per.point = 5000), 
                             vi_args = list(num.points.max = 150000)) %>% 
  print(plot(hclust(as.dist(1-.), method = "complete" ), 
             axes = F, xlab = NULL, ylab = NULL, main = NULL, 
             sub = NULl, ann = F, cex = .75))


#=====================================================================
# Demonstrative plots in the paper
#=====================================================================

# Code for 2D chaos game example 
seq_to_2Dcg <- function(dna_seq, plot=F, ...){
  A <- c(0, 0)
  C <- c(0, 1)
  G <- c(1, 1)
  T <- c(1, 0)
  if(length(dna_seq) == 1) dna_seq <- str_split(dna_seq, "")[[1]]
  cg_x <- c(1/2)
  cg_y <- c(1/2)
  for(i in seq_along(dna_seq)){
    cg_x <- c(cg_x, mean(c(cg_x[length(cg_x)], get(dna_seq[i])[1])))
    cg_y <- c(cg_y, mean(c(cg_y[length(cg_y)], get(dna_seq[i])[2])))
  }
  xy <- data.frame(cg_x, cg_y, nucleotide = c("-", dna_seq)) |> 
    rownames_to_column() |> 
    mutate(position = as.numeric(rowname)-1,
           subsequence = "") 
  
  for(i in 2:nrow(xy)){
    xy[["subsequence"]][i] <- paste(xy[["subsequence"]][i-1], xy[["nucleotide"]][i] )
  }
  
  return(xy)
}

# Figure 1.2 An example of 2D chaos game representation for ATGG
rows <- 5
seq_to_2Dcg(beta_seq[[1]]) |> head(rows) |> 
  ggplot(aes(label = subsequence, x = cg_x, y = cg_y)) + 
  geom_text(position = position_nudge(x = -.05)) +
  geom_path(aes(col = as.numeric(rowname)), 
            arrow = arrow(ends = "last", type = "closed", 
                          length = unit(c(rep(.15, rows), 0.25), "inches")), 
            alpha = .75) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() + xlab(NULL) + ylab(NULL) +
  guides(col = "none") + 
  annotate(geom = "text", 
           label = c("A", "C", "G", "T"), 
           x = c(0, 0, 1, 1), y = c(0, 1, 1, 0), size = 7) 


# Figure 2.1. 3D CGR of the MT-NADH for the Macaca fascicular
for_plotly <- seq_to_hypercomplex_cg4(nadh[[1]]) %>% rownames_to_column() %>% 
  dplyr::rename(position = rowname) %>% mutate(position = as.numeric(position))
fig <- plot_ly(for_plotly, x = ~i, y = ~j, z = ~k, 
               size = I(300), color = ~position, opacity = .85)
fig <- fig %>% add_markers()
fig


# Figure 3.2 Sequence alignment of the CDS of theβ-globin gene for the chimp, human, and rabbit
alignment_to_vector <- function(alignment){
  output <- alignment$seq
  names(output) <- alignment$nam
  return(output)
}
msa(beta_seq[c("chimp", "human", "rabbit")], type = "dna", method = "ClustalOmega") %>% 
  msaConvert(., type = "seqinr::alignment") |> 
  alignment_to_vector() |> DNAStringSet() |> 
  ggmsa::ggmsa( 1, 52, 
                color = "Shapely_NT", font = "DroidSansMono", char_width = 0.5, 
                seq_name = T, consensus_views = T  ) 
msa(beta_seq[c("chimp", "human", "rabbit")], type = "dna", method = "ClustalOmega") %>% 
  msaConvert(., type = "seqinr::alignment") |> 
  alignment_to_vector() |> DNAStringSet() |> 
  ggmsa::ggmsa( 53, 105, 
                color = "Shapely_NT", font = "DroidSansMono", char_width = 0.5,
                seq_name = T, consensus_views = T  ) 
