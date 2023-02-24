#####################################################################
# Comparing the effect of mutation accumulation between
# substitution and deletion on the distance from the parent sequence
# This code is takes a while to run.  You can shorten the run time a bit
# by decreasing the number of bin counts/bandwidths checked.
# By: Stephanie Young <syoung2@sdsu.edu> <syoung49@its.jnj.com>
#####################################################################

set.seed(0)

parent_dna <- beta_seq[[2]]

# Code to mutate a single nucleotide
mtte <- function(nucleotide){
  sample(setdiff(c("A", "T", "C", "G"), nucleotide), 1)
}

mtte.Vectorize <- Vectorize(mtte, "nucleotide")

# Select mutation locations

intervals <- 2
num_mutations <- seq(0, 76, by = intervals)

mutation_locations <- c()
for(i in seq_along(num_mutations)[-1]){
  new_locations <- sample(setdiff(1:nchar(parent_dna), mutation_locations), intervals)
  mutation_locations <- c(mutation_locations, c(new_locations))
}

agg_sub_fas <- c(`Substitutions = 0` = parent_dna)
agg_del_fas <- c(`Deletions = 0` = parent_dna)

for(i in seq_along(num_mutations)[-1]){
  # Deletion mutations
  previous_del <- str_split(agg_del_fas[i - 1], "")[[1]]
  del_seq <- previous_del
  del_seq[mutation_locations[(1:intervals) + (i - 2) * intervals]] <- "-"
  agg_del_fas <- c(agg_del_fas, paste(del_seq, collapse = ""))
  names(agg_del_fas)[i] <- paste("Deletions =", num_mutations[i])
  
  # Substitution mutations
  previous_sub <- str_split(agg_sub_fas[i - 1], "")[[1]]
  sub_seq <- previous_sub
  sub_seq[mutation_locations[(1:intervals) + (i - 2) * intervals]] <- 
    mtte.Vectorize(previous_sub[mutation_locations[(1:intervals) + (i - 2) * intervals]])
  agg_sub_fas <- c(agg_sub_fas, paste(sub_seq, collapse = ""))
  names(agg_sub_fas)[i] <- paste("Substitutions =", num_mutations[i])
}

agg_fas <- c(agg_sub_fas[-1], gsub("-","",agg_del_fas)[-1])

all_mutations <- c(`Parent` = parent_dna, agg_fas)

#=====================================================================
# Mutation aggregation effects under Clustal Omega
#=====================================================================

clustal_output <- msa(all_mutations, type = "dna", method = "ClustalOmega") |> 
  msaConvert(type = "seqinr::alignment") |> 
  dist.alignment(gap = T) |> 
  as.matrix() %>% reshape2::melt()
clustal_output <- 
  clustal_output |> 
  filter(Var1 == "Parent" & Var2 != "Parent") |> 
  mutate(`Mutation count` = as.numeric(gsub("(Substitutions|Deletions) = ", "", Var2)),
         `Mutation percent` = `Mutation count`/nchar(parent_dna),
         `Mutation type` = gsub(" = [0-9]*", "", Var2)) 

clustal_substitutions <- clustal_output |> filter(`Mutation type` == "Substitutions") 
clustal_deletions <- clustal_output |> filter(`Mutation type` == "Deletions") 

# Figure 3.5(a)
# Plots of the distance from the parent sequence as a function of the 
# percent of DNA that is mutated. 
clustal_output %>% 
  ggplot(aes(x=`Mutation percent`, y = value,col = `Mutation type`, pch = `Mutation type`)) + 
  geom_point() + 
  geom_line(stat = "smooth", method="loess",
            method.args = list(degree = 1), se = F, alpha = .25, lwd = 2, span = .75) +
  theme_minimal() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  ylab("Distance") + theme(legend.position = "bottom")

# Figure 3.6(a) 
# Distance from parent sequence due to deletion mutation v 
# distance from parent sequence due to substitution mutation, 
# colored by proportion of sequence that has been mutated.
full_join(clustal_substitutions, clustal_deletions, 
          by = c("Var1", "Mutation count", "Mutation percent"), 
          suffix=c(".S", ".D")) |>
  ggplot() + 
  geom_point(aes(x=value.S, y=value.D, col = `Mutation count`), alpha = .75, size = 3) + 
  theme_minimal() +  
  geom_abline(slope = 1, intercept = 0, lty = "dashed") +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  xlab("Distance from parent sequence\ndue to substitution mutation") + 
  ylab("Distance from parent sequence\ndue to deletion mutation") +
  scale_color_continuous(low = "blue", high = "orange")

#=====================================================================
# Mutation aggregation effects under volume intersection method
#=====================================================================

bandwidth_options <- 
  c(seq(0.0001, 0.0009, by = 0.0001),
    seq(0.001, 0.009, by = 0.001),
    seq(0.01, 0.1, by = 0.01))

volume_results <- 
  lapply(bandwidth_options, 
         function(x){ 
           volume_intersection_tanimoto(all_mutations, x)
         })

volume_results_df <- 
  lapply(seq_along(bandwidth_options), 
         function(x){
           melt(volume_results[x], value.name = "tanimoto") |> 
             filter(Var2 == "Parent" & Var1 != "Parent") |> 
             mutate(bandwidth_value = bandwidth_options[x],
                    `Mutation type` =  gsub(" = [0-9]*", "", Var1),
                    `Mutation count` = as.numeric(gsub("(Substitutions|Deletions) = ","", Var1)),
                    `Mutation percent` = `Mutation count`/nchar(parent_dna))
         }) 

volume_results_df <- do.call(rbind, volume_results_df)

sub_vol <- volume_results_df |> filter(`Mutation type` == "Substitutions")
del_vol <- volume_results_df |> filter(`Mutation type` == "Deletions")

# Figure 3.5(b)
# Distance from parent sequence as a function of percent of DNA that is mutated 
rbind(sub_vol, del_vol) |>   
  ggplot() + 
  geom_point(aes(x=`Mutation percent`, y = 1-tanimoto, col = bandwidth_value)) + 
  facet_wrap(~`Mutation type`) +
  scale_color_continuous("Bandwidth", low = "blue", high = "orange", trans = "log") +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  theme_minimal() + scale_x_log10() + 
  ylab("Distance") 


# Figure 3.6 (b)
# Distance from parent sequence due to deletion mutation v 
# distance from parent sequence due to substitution mutation, 
# colored by proportion of sequence that has been mutated.
full_join(sub_vol, del_vol, 
          by = c("Mutation count", "Mutation percent", "bandwidth_value"), suffix = c(".S", ".D")) |>
  ggplot() + 
  geom_point(aes(x=1-tanimoto.S, y=1-tanimoto.D, col = `Mutation percent`), alpha = .5) + 
  scale_color_continuous("Mutation\npercent", low = "blue", high = "orange") + 
  theme_minimal() +  
  geom_abline(slope = 1, intercept = 0, lty =2) +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  xlab("Distance from parent sequence\ndue to substitution mutation") + 
  ylab("Distance from parent sequence\ndue to deletion mutation")

# Figure 3.7 (a)
# Distance from parent sequence due to deletion mutation v 
# distance from parent sequence due to substitution mutation, 
# colored by bandwidth
full_join(sub_vol, del_vol, 
          by = c("Mutation count", "Mutation percent", "bandwidth_value"), suffix = c(".S", ".D")) |>
  ggplot() + 
  geom_point(aes(x=1-tanimoto.S, y=1-tanimoto.D, col = bandwidth_value), alpha = .5) + 
  scale_color_continuous("Bandwidth", low = "blue", high = "orange", trans = "log") + 
  theme_minimal() +  
  geom_abline(slope = 1, intercept = 0, lty =2) +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  xlab("Distance from parent sequence\ndue to substitution mutation") + 
  ylab("Distance from parent sequence\ndue to deletion mutation")



#=====================================================================
# Mutation aggregation effects under shape signature methods
#=====================================================================

# Function to combine deletion and substitution results 
combine_sub_del <- function(bin_counts, sub_results, del_results, sig = "", DISTFUNC = dist){
  sub_del_comb <- data.frame()
  for(i in 1:length(sub_results)){
    i_count <- bin_counts[i]
    i_processed_substitutions <- 
      sub_results[[i]] |> DISTFUNC() |> as.matrix() |> as.data.frame() |> 
      rownames_to_column() |> 
      pivot_longer(cols = !contains("rowname"), values_to = "distance") |> 
      filter(rowname != name & grepl("= 0$", rowname)) |> 
      mutate(`Mutation count` = as.numeric(gsub("[a-zA-Z]* = *([0-9]*)$", "\\1", name)), 
             info = "processed", 
             bin_count = i_count) 
    
    i_processed_deletions <- 
      del_results[[i]] |> DISTFUNC() |> as.matrix() |> as.data.frame() |> 
      rownames_to_column() |> 
      pivot_longer(cols = !contains("rowname"), values_to = "distance") |> 
      filter(rowname != name & grepl("= 0$", rowname)) |> 
      mutate(`Mutation count` = as.numeric(gsub("[a-zA-Z]* = *([0-9]*)$", "\\1", name)), 
             info = "processed", 
             bin_count = i_count) 
    
    sub_del_comb <- rbind(sub_del_comb, 
                              full_join(i_processed_substitutions, 
                                        i_processed_deletions, 
                                        by = c("Mutation count", "info", "bin_count")) |> 
                                mutate(info = "processed", signature = sig))
  }
  return(sub_del_comb)
}

bin_counts <- 
  c(seq(5, 95, by = 10), 
    seq(100, 900, by = 100), 
    seq(1000, 9000, by = 1000), 
    seq(10000, 90000, by = 10000))

# Angular signature --------------------------------------------------

check_substitutions <- list()
for(i in bin_counts){
  print(i)
  temp <- 
    feature_signature(lapply(agg_sub_fas, 
                             function(x) by3rowangle(seq_to_hypercomplex_cg4(x))),
                      bin_count = i)
  check_substitutions <- append(check_substitutions, list(temp))
}

check_deletions<- list()
for(i in bin_counts){
  print(i)
  temp <- 
    feature_signature(lapply(gsub("-", "", agg_del_fas), 
                             function(x) by3rowangle(seq_to_hypercomplex_cg4(x))),
                      bin_count = i)
  check_deletions <- append(check_deletions, list(temp))
}

sub_del_comb_ang <- combine_sub_del(bin_counts, 
                                    check_substitutions, check_deletions, 
                                    sig = "angular")


# Edge signature --------------------------------------------------

check_substitutions <- list()
for(i in bin_counts){
  print(i)
  temp <- 
    feature_signature(lapply(agg_sub_fas, 
                             function(x) by3rowdistance(seq_to_hypercomplex_cg4(x))),
                      bin_count = i)
  check_substitutions <- append(check_substitutions, list(temp))
}

check_deletions<- list()
for(i in bin_counts){
  print(i)
  temp <- 
    feature_signature(lapply(gsub("-", "", agg_del_fas), 
                             function(x) by3rowdistance(seq_to_hypercomplex_cg4(x))),
                      bin_count = i)
  check_deletions <- append(check_deletions, list(temp))
}

sub_del_comb_edg <- combine_sub_del(bin_counts, 
                                    check_substitutions, check_deletions,
                                    sig = "edge")

# Coordinate signature -----------------------------------------------

check_substitutions <- list()
for(i in bin_counts){
  print(i)
  temp <- 
    coordinate_signature(lapply(agg_sub_fas, seq_to_hypercomplex_cg4),
                         bin_count = i)
  check_substitutions <- append(check_substitutions, list(temp))
}

check_deletions<- list()
for(i in bin_counts){
  print(i)
  temp <- 
    coordinate_signature(lapply(gsub("-", "", agg_del_fas), seq_to_hypercomplex_cg4),
                         bin_count = i)
  check_deletions <- append(check_deletions, list(temp))
}

sub_del_comb_coord <- combine_sub_del(bin_counts, 
                                      check_substitutions, check_deletions, 
                                      sig = "coordinate")

# Figure 3.5 (c)
# Distance from parent sequence as a function of percent of DNA that is mutated 
# by substitution
rbind(sub_del_comb_ang, 
      sub_del_comb_edg, 
      sub_del_comb_coord) |> 
  filter(bin_count %in% c(15, 45, 100, 500, 1000, 5000, 10000)) |>
  ggplot(aes(x=`Mutation count`/nchar(parent_dna), y = distance.x)) + 
  geom_point(aes(col = bin_count)) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  ylab("Distance induced by substitution") + xlab("Mutation percent") +
  theme_minimal() + facet_grid(~signature, scales = "free_y") + scale_x_log10() + 
  scale_color_continuous("bin\ncount", low = "blue", high = "orange", trans = "log", labels = function(x) sprintf("%.f", x)) +
  geom_line(stat = "smooth", method = "loess", aes(col = bin_count, group = bin_count), method.args = list(degree = 1), se = F, alpha = .25, lwd = 2, span = .75) +
  theme(legend.position = "bottom")


# Figure 3.5 (d)
# Distance from parent sequence as a function of percent of DNA that is mutated 
# by deletion
rbind(sub_del_comb_ang, 
      sub_del_comb_edg, 
      sub_del_comb_coord) |> 
  filter(bin_count %in% c(15, 45, 100, 500, 1000, 5000, 10000)) |>
  ggplot(aes(x=`Mutation count`/nchar(parent_dna), y = distance.y)) + 
  geom_point(aes(col = bin_count)) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  ylab("Distance induced by deletion") + xlab("Mutation percent") +
  theme_minimal() + facet_grid(~signature, scales = "free_y") + scale_x_log10() + 
  scale_color_continuous("bin\ncount", low = "blue", high = "orange", trans = "log", labels = function(x) sprintf("%.f", x)) +
  geom_line(stat = "smooth", method = "loess", aes(col = bin_count, group = bin_count), method.args = list(degree = 1), se = F, alpha = .25, lwd = 2, span = .75) +
  theme(legend.position = "bottom")




















