#####################################################################
# Comparing the effect of mutation accumulation between
# substitution and deletion on the distance from the parent sequence
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

#=====================================================================
# Mutation aggregation effects under Clustal Omega
#=====================================================================

agg_fas <- c(agg_sub_fas[-1], gsub("-","",agg_del_fas)[-1])

temp <- c(`Parent` = parent_dna, agg_fas)
clustal_output <- msa(temp, type = "dna", method = "ClustalOmega") |> 
  msaConvert(type = "seqinr::alignment") |> 
  dist.alignment(gap = T) |> 
  as.matrix() %>% reshape2::melt()
clustal_output <- 
  clustal_output |> 
  filter(Var1 == "Parent" & Var2 != "Parent") |> 
  mutate(`Mutation count` = as.numeric(gsub("(Substitutions|Deletions) = ", "", Var2)),
         `Mutation percent` = `Mutation count`/nchar(parent_seq),
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

