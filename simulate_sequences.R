# Simulating my own damn sequences
set.seed(0)

# Substitution mutation
mtte <- function(nucleotide){
  # Give the function vector of nucleotides and it will replace each one
  sapply(nucleotide, function(x) sample(setdiff(c("A", "C", "G", "T"), x), 1))
}

# Insertion mutation
insert_at <- function(a, pos, insert){
  # Give the function a vector of the sequence of nucleotides and it will insert 
  # new nucleotides (insert) at every position (pos)
  if(min(pos) == 1){
    result <- vector("list",2*length(pos))
    result[c(FALSE,TRUE)] <- split(a, cumsum(seq_along(a) %in% (pos)))
    result[c(TRUE,FALSE)] <- insert[order(pos)]
  } else{
    result <- vector("list",2*length(pos)+1)
    result[c(TRUE, FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos)))
    result[c(FALSE, TRUE)] <- insert[order(pos)]
  }
  unlist(result)
}

# Transposition mutation
transposition <- function(orig_seq, bp_length, orig_loc, new_loc){
  # Give the function a vector of the sequence of nucleotides and it transpose it,
  # moving it from one location to another
  if(orig_loc > length(orig_seq) - bp_length + 1){
    print("There is an error with the input")
    return(NULL)
  } else {
    cut_index <- orig_loc : (orig_loc + bp_length -1)
    cut_seq <- orig_seq[cut_index]
    output <- insert_at(orig_seq, new_loc, list(cut_seq))
    if(new_loc > orig_loc) output <- output[-cut_index]
    else output <- output[-(cut_index + bp_length)]
    return(output)
  }
}

# Substitution mutation
substitution <- function(orig_seq, mut_loc){
  # Give the function a vector of the sequence of nucleotides it will mutate all the
  # nucleotides at mut_loc
  orig_seq[mut_loc] <- mtte(orig_seq[mut_loc])
  return(orig_seq)
}

# deletion mutation
deletion <- function(orig_seq, mut_loc){
  orig_seq[mut_loc] <- "-"
  return(orig_seq)
}

# insertion mutation
insertion <- function(orig_seq, mut_loc){
  elements <- sample(c("A", "C", "G", "T"), length(mut_loc), replace = T) |> as.list()
  insert_at(
    orig_seq, 
    mut_loc, 
    elements
  )
}

# duplication mutation
duplication <- function(orig_seq, mut_loc, times = 1){
  duplicate_this <- orig_seq[mut_loc]
  insert_this <- list(rep(duplicate_this, times))
  insert_at(orig_seq, mut_loc[1], insert_this)
}

sample_interval <- function(x, size, interval_length){
  # Sample interval_length contiguous/adjacent elements from x size times.
  # Does not allow overlapping of sampled intervals
  # Does not allow repeat sample of the same interval
  indices <- 1:(length(x) - interval_length + 1)
  sampled_indices <- sample(indices, 1)
  sampled_indices <- sampled_indices : (sampled_indices + interval_length - 1)
  times <- 1
  while(times < size){
    remaining_viable_indices <- setdiff(1:length(x), sampled_indices)
    check <- lead(remaining_viable_indices, interval_length - 1) - remaining_viable_indices
    remaining_viable_indices <- remaining_viable_indices[!(is.na(check) | check >= interval_length)]
    if(length(remaining_viable_indices) > 0){
      next_sample <- sample(remaining_viable_indices, 1)
      next_sample <- next_sample : (next_sample + interval_length - 1)
      sampled_indices <- c(sampled_indices, next_sample)
      times <- times + 1
    } else{
      print("Ran out of nucleotides to sample")
      times <- size
    }
  }
  return(sampled_indices)
}


dna_mutation <- function(orig_seq, num_mutations, interval_size = 1, lineage = T, mutation_locations = NULL){
  if(length(orig_seq == 1)) orig_seq <- str_split(orig_seq, "")[[1]]
  del_seq <- sub_seq <- ins_seq <- dup_seq <- orig_seq 
  
  # The first generation sequences
  del_line <- sub_line <- ins_line <- dup_line <- 
    vector(length = length(num_mutations) + 1)
  del_line[1] <- sub_line[1] <- ins_line[1] <- dup_line[1] <- paste0(orig_seq, collapse = "")
  names(del_line)[1] <- "Deletions = 0"
  names(sub_line)[1] <- "Substitutions = 0"
  names(ins_line)[1] <- "Insertions = 0"
  names(dup_line)[1] <- "Duplications = 0"

  if(lineage){
    # Mutations in a lineage are aggregative
    # Under this case, (num_)mutations must increase 
    if(any(diff(num_mutations) < 0)) print("Incorrect use of function. If lineage = T, then the num_mutations input must always be increasing.")
    # intervals <- num_mutations - c(0, lag(num_mutations)[-1]) Looks like I never needed this
    if(is.null(mutation_locations)) 
      mutation_locations <- sample_interval(1:length(orig_seq), last(num_mutations), interval_size)
    
    # Some last generation sequences
    sub_line.last <- sub_seq
    sub_line.last[mutation_locations] <- mtte(sub_seq[mutation_locations])
    
    ins_line.last <- insert_at(
      orig_seq, 
      mutation_locations[seq(1, length(mutation_locations), by = interval_size)], 
      lapply(1:last(num_mutations), function(x) sample(c("A", "C", "G", "T"), interval_size, replace = T))
    )
    full_ins_helper <- insert_at(
      orig_seq, 
      mutation_locations[seq(1, length(mutation_locations), by = interval_size)], 
      lapply(1:last(num_mutations), function(x) rep(x, each = interval_size))
    ) |> (\(x) gsub("[a-zA-Z]", "0", x))()
    

    dup_line.last <- insert_at(
      orig_seq, 
      mutation_locations[seq(1, length(mutation_locations), by = interval_size)], 
      lapply(mutation_locations[seq(1, length(mutation_locations), by = interval_size)], function(x) orig_seq[x : (x+interval_size-1)])
    )
    
    for(i in seq_along(num_mutations)){
      # Deletion mutations
      del_seq[mutation_locations[1:(interval_size * num_mutations[i])]] <- "-"
      del_line[i+1] <- paste(del_seq, collapse = "")
      names(del_line)[i+1] <- paste("Deletions =", num_mutations[i], "; Length =", interval_size)

      # Insertion mutations
      ins_seq <- ins_line.last[which(full_ins_helper <= num_mutations[i])]
      ins_line[i+1] <- paste(ins_seq, collapse = "")
      names(ins_line)[i+1] <- paste("Insertions =", num_mutations[i], "; Length =", interval_size)
      
      # Duplication mutations
      dup_seq <- dup_line.last[which(full_ins_helper <= num_mutations[i])]
      dup_line[i+1] <- paste(dup_seq, collapse = "")
      names(dup_line)[i+1] <- paste("Duplications =", num_mutations[i], "; Length =", interval_size)
      
      # Substitution mutations
      sub_seq[mutation_locations[1:(interval_size * num_mutations[i])]] <- sub_line.last[mutation_locations[1:(interval_size * num_mutations[i])]]
      sub_line[i+1] <- paste(sub_seq, collapse = "")
      names(sub_line)[i+1] <- paste("Substitutions =", num_mutations[i], "; Length =", interval_size)
    }
    print(mutation_locations)
  }
  else{
    for(i in seq_along(num_mutations)){
      if(is.null(mutation_locations) )
         mutation_locations <- sample_interval(1:length(orig_seq), num_mutations[i], interval_size)
      del_seq <- sub_seq <- ins_seq <- dup_seq <- orig_seq
      
      # Deletion mutations
      del_seq[mutation_locations[1:(interval_size * num_mutations[i])]] <- "-"
      del_line[i + 1] <- paste(del_seq, collapse = "")
      names(del_line)[i + 1] <- paste("Deletions =", num_mutations[i], "; Length =", interval_size)
      
      # Insertion mutations
      ins_seq <- insert_at(
        orig_seq, 
        mutation_locations[seq(1, length(mutation_locations), by = interval_size)], 
        lapply(1:num_mutations[i], function(x) sample(c("A", "C", "G", "T"), interval_size, replace = T))
      )
      ins_line[i + 1] <- paste(ins_seq, collapse = "")
      names(ins_line)[i + 1] <- paste("Insertions =", num_mutations[i], "; Length =", interval_size)
      
      # Duplication mutations
      dup_seq <- insert_at(
        orig_seq, 
        mutation_locations[seq(1, length(mutation_locations), by = interval_size)], 
        lapply( mutation_locations[seq(1, length(mutation_locations), by = interval_size)] , function(x) orig_seq[x:(x + interval_size - 1)])
      )
      dup_line[i + 1] <- paste(dup_seq, collapse = "")
      names(dup_line)[i+1] <- paste("Duplications =", num_mutations[i], "; Length =", interval_size)

      # Substitution mutations
      sub_seq[mutation_locations] <- mtte(sub_seq[mutation_locations])
      sub_line[i + 1] <- paste(sub_seq, collapse = "")
      names(sub_line)[i + 1] <- paste("Substitutions =", num_mutations[i], "; Length = ", interval_size)
      print(mutation_locations)
      mutation_locations <- NULL
    }
  }
  
  return(list(del_line, ins_line, sub_line, dup_line))
}


#yinish_fas

set.seed(100)

yin_fas <- 
  c(
    sim_fas["A"],
    dna_mutation(
      sim_fas["A"], 
      num_mutations = c(2, 2, 5, 5, 10, 10), 
      lineage = F
    )[[3]] |> (\(x) x[-1])(),
    sim_fas["B"],
    dna_mutation(
      sim_fas["B"], 
      num_mutations = c(2, 2, 5, 5, 10, 10), 
      lineage = F
    )[[3]] |> (\(x) x[-1])(),
    dna_mutation(
      sim_fas["B"],
      num_mutations = c(1),
      interval_size = 5,
      mutation_locations = c(51:55),
      lineage = F) |> sapply(`[`, 2) |> (\(x) x[grepl("Insertion|Deletion", names(x))])(),
    dna_mutation(
      sim_fas["B"],
      num_mutations = c(1),
      interval_size = 5,
      mutation_locations = c(101:105),
      lineage = F) |> sapply(`[`, 2)  |> (\(x) x[grepl("Insertion|Deletion", names(x))])(),
    transposition(
      str_split(sim_fas["B"], "")[[1]], bp_length = 5, orig_loc = 50, new_loc = 150) |> 
      paste(collapse = ""),
    transposition(
      str_split(sim_fas["B"], "")[[1]], bp_length = 5, orig_loc = 50, new_loc = 250) |> 
      paste(collapse = "")
  ) |> toupper() |> (\(x) gsub("-", "", x))()
names(yin_fas) <-
  c(
    "A", 
    "A/\u2206 2NT/1", "A/\u2206 2NT/2", 
    "A/\u2206 5NT/1", "A/\u2206 5NT/2", 
    "A/\u2206 10NT/1", "A/\u2206 10NT/2", 
    "B", 
    "B/\u2206 2NT/1", "B/\u2206 2NT/2", 
    "B/\u2206 5NT/1", "B/\u2206 5NT/2", 
    "B/\u2206 10NT/1", "B/\u2206 10NT/2", 
    "B/-5NT/51:55", "B/+5NT/51:55",
    "B/-5NT/101:105", "B/+5NT/101:105",
    "B/\U00B1 5NT/50:54 -> 150:154", "B/\U00B1 5NT/50:54 -> 250:254"
  )
clustal_omega(yin_fas[-(1:7)]) |> upgma() |> plot()
clustal_omega(yin_fas[-(1:7)]) |> bionj() |> plot()
cgr_distance(lapply(yin_fas[-(1:7)], seq_to_cgr), frac = .15) |> upgma() |> plot() 
cgr_distance(lapply(yin_fas[-(1:7)], seq_to_cgr), frac = 1) |> upgma() |> plot() 
cgr_distance(lapply(yin_fas[-(1:7)], seq_to_cgr), frac = .15) |> bionj() |> plot() 
cgr_distance(lapply(yin_fas[-(1:7)], seq_to_cgr), frac = 1) |> bionj() |>  plot() 




set.seed(0)
A.1 <- dna_mutation(my_fas[1], num_mutations = c(2, 5, 10), lineage = T)[[3]]
A.2to4 <- dna_mutation(my_fas[1], num_mutations = c(2, 5, 10), lineage = F)[[3]]
B.1to6 <- dna_mutation(my_fas[2], num_mutations = c(2, 2, 5, 5, 10, 10), lineage = F)[[3]]
B.7 <- dna_mutation(my_fas[2], num_mutations = 1, interval_size = 5, lineage = F, mutation_locations = 51:55)[[2]]
B.8 <- dna_mutation(my_fas[2], num_mutations = 1, interval_size = 5, lineage = F, mutation_locations = 101:105)[[2]]
B.9 <- dna_mutation(my_fas[2], num_mutations = 1, interval_size = 5, lineage = F, mutation_locations = 51:55)[[1]]
B.10 <- dna_mutation(my_fas[2], num_mutations = 1, interval_size = 5, lineage = F, mutation_locations = 101:105)[[1]]
B.11 <- transposition(str_split(my_fas[2], "")[[1]], bp_length = 5, orig_loc = 50, new_loc = 150) |> paste(collapse = "")
B.11.1 <- transposition(str_split(B.12, "")[[1]], bp_length = 5, orig_loc = 300, new_loc = 500) |> paste(collapse = "")
B.12 <- transposition(str_split(my_fas[2], "")[[1]], bp_length = 5, orig_loc = 50, new_loc = 250) |> paste(collapse = "")

new_fas <-
  c(
    my_fas[1], 
    `A.01/2NT Substitutions` = A.1[[which(grepl("Substitutions = 2", names(A.1)))]],
    `A.01.01/2NT+3NT=5NT Substitutions` = A.1[[which(grepl("Substitutions = 5", names(A.1)))]],
    `A.01.01.01/2NT+3NT+5NT=10NT Substitutions` = A.1[[which(grepl("Substitutions = 10", names(A.1)))]],
    `A.02/2NT Substitutions` = A.2to4[[which(grepl("Substitutions = 2", names(A.1)))]],
    `A.03/5NT Substitutions` = A.2to4[[which(grepl("Substitutions = 5", names(A.1)))]],
    `A.04/10NT Substitutions` = A.2to4[[which(grepl("Substitutions = 10", names(A.1)))]],
    my_fas[2],
    `B.01/2NT Substitutions` = B.1to6[[2]],
    `B.02/2NT Substitutions` = B.1to6[[3]],
    `B.03/5NT Substitutions` = B.1to6[[4]],
    `B.04/5NT Substitutions` = B.1to6[[5]],
    `B.05/10NT Substitutions` = B.1to6[[6]],
    `B.06/10NT Substitutions` = B.1to6[[7]],
    `B.07/+5NT Insertion/51:55` = B.7[[2]],
    `B.08/+5NT Insertion/101:105` = B.8[[2]],
    `B.09/-5NT Deletion/51:55` = B.9[[2]],
    `B.10/-5NT Deletion/101:105` = B.10[[2]],
    `B.11/5NT Transposition/50->150` = B.11,
    `B.11.1/5NT+5NT Transposition/50->150,300->500` = B.11.1,
    `B.12/5NT Transpotition/50->250` = B.12
  )

























set.seed(0)
A.1 <- dna_mutation(my_fas[1], num_mutations = c(2, 4, 6), lineage = T)[[1]]
A.2 <- dna_mutation(my_fas[1], num_mutations = c(2, 4, 6), lineage = T)[[2]]
A.3 <- dna_mutation(my_fas[1], num_mutations = c(2, 4, 6), lineage = T)[[3]]
B.1to6 <- dna_mutation(my_fas[2], num_mutations = c(2, 2, 5, 5, 10, 10), lineage = F)[[3]]
B.7 <- dna_mutation(my_fas[2], num_mutations = 1, interval_size = 5, lineage = F, mutation_locations = 51:55)[[2]]
B.8 <- dna_mutation(my_fas[2], num_mutations = 1, interval_size = 5, lineage = F, mutation_locations = 101:105)[[2]]
B.9 <- dna_mutation(my_fas[2], num_mutations = 1, interval_size = 5, lineage = F, mutation_locations = 51:55)[[1]]
B.10 <- dna_mutation(my_fas[2], num_mutations = 1, interval_size = 5, lineage = F, mutation_locations = 101:105)[[1]]
B.11 <- transposition(str_split(my_fas[2], "")[[1]], bp_length = 5, orig_loc = 50, new_loc = 150) |> paste(collapse = "")
B.11.1 <- transposition(str_split(B.11, "")[[1]], bp_length = 5, orig_loc = 300, new_loc = 500) |> paste(collapse = "")
B.12 <- transposition(str_split(my_fas[2], "")[[1]], bp_length = 5, orig_loc = 60, new_loc = 250) |> paste(collapse = "")

new_fas <-
  c(
    my_fas[1], 
    `A-01/-2NT Deletions` = A.1[[which(grepl("Deletions = 2", names(A.1)))]],
    `A-01-01/-2NT-2NT=4NT Deletions` = A.1[[which(grepl("Deletions = 4", names(A.1)))]],
    `A-01-01-01/-2NT-2NT-2NT=6NT Deletions` = A.1[[which(grepl("Deletions = 6", names(A.1)))]],
    `A-02/+2NT Insertions` = A.2[[which(grepl("Insertions = 2", names(A.2)))]],
    `A-02-01/+2NT+2NT=4NT Insertions` = A.2[[which(grepl("Insertions = 4", names(A.2)))]],
    `A-02-01-01/+2NT+2NT+2=6NT Insertions` = A.2[[which(grepl("Insertions = 6", names(A.2)))]],
    `A-03/2NT Substitutions` = A.3[[which(grepl("Substitutions = 2", names(A.3)))]],
    `A-03-01/2NT+2NT=4NT Substitutions` = A.3[[which(grepl("Substitutions = 4", names(A.3)))]],
    `A-03-01-01/2NT+2NT+2NT=6 Substitutions` = A.3[[which(grepl("Substitutions = 6", names(A.3)))]],
    my_fas[2],
    `B-01/2NT Substitutions` = B.1to6[[2]],
    `B-02/2NT Substitutions` = B.1to6[[3]],
    `B-03/5NT Substitutions` = B.1to6[[4]],
    `B-04/5NT Substitutions` = B.1to6[[5]],
    `B-05/10NT Substitutions` = B.1to6[[6]],
    `B-06/10NT Substitutions` = B.1to6[[7]],
    `B-07/+5NT Insertion/51:55` = B.7[[2]],
    `B-08/+5NT Insertion/101:105` = B.8[[2]],
    `B-09/-5NT Deletion/51:55` = B.9[[2]],
    `B-10/-5NT Deletion/101:105` = B.10[[2]],
    `B-11/5NT Transposition/50->150` = B.11,
    `B-11-1/5NT+5NT Transposition/50->150,300->500` = B.11.1,
    `B-12/5NT Transposition/60->250` = B.12
  )


new_fas.align <- 
  c(
    `B-07/+5NT Insertion/51:55` = insert_at(str_split(B.7[[2]], "")[[1]], c(101+5, 150+5, 250+5, 500+5), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5))) |> (\(x)paste(x, collapse = ""))(),
    `B-08/+5NT Insertion/101:105` = insert_at(str_split(B.8[[2]], "")[[1]], c(51, 150+5, 250+5, 500+5), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5))) |> (\(x)paste(x, collapse = ""))(),

    `B/Original` = insert_at(str_split(my_fas[2], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-01/2NT Substitutions` = insert_at(str_split(B.1to6[[2]], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-02/2NT Substitutions` = insert_at(str_split(B.1to6[[3]], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-03/5NT Substitutions` = insert_at(str_split(B.1to6[[4]], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-04/5NT Substitutions` = insert_at(str_split(B.1to6[[5]], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-05/10NT Substitutions` = insert_at(str_split(B.1to6[[6]], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-06/10NT Substitutions` = insert_at(str_split(B.1to6[[7]], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),

    `B-09/-5NT Deletion/51:55` = insert_at(str_split(B.9[[2]], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-10/-5NT Deletion/101:105` = insert_at(str_split(B.10[[2]], "")[[1]], c(51, 101, 150, 250, 500), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),

    `B-11/5NT Transposition/50->150` =  insert_at(insert_at(str_split( B.11, "")[[1]], c(50), list(rep("-", 5))), c(51, 101, 255, 505), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-11-1/5NT+5NT Transposition/50->150,300->500` = insert_at(insert_at(str_split( B.11.1, "")[[1]], c(50), list(rep("-", 5))), c(51, 101, 250, 305), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))(),
    `B-12/5NT Transposition/60->250` = insert_at(insert_at(str_split( B.12, "")[[1]], c(60), list(rep("-", 5))), c(51, 101, 150, 505), list(rep("-", 5), rep("-", 5), rep("-", 5), rep("-", 5)))|> (\(x)paste(x, collapse = ""))()
  )
do.call("rbind", lapply(gsub("-","N",new_fas.align), function(x) str_split(tolower(x), "")[[1]]))  |> as.DNAbin() |> ape::dist.gene() |> sqrt() |> ape::bionj() |> plot()
image(do.call("rbind", lapply(gsub("-","N",new_fas.align), function(x) str_split(tolower(x), "")[[1]]))  |> as.DNAbin() )

set.seed(0)
test <- c(
  dna_mutation(my_fas[1], num_mutations = 1:10, lineage = T)[[3]],
  dna_mutation(my_fas[1], num_mutations = 1:10, lineage = T)[[1]][-1],
  dna_mutation(my_fas[2], num_mutations = 1:10, lineage = T)[[3]]
)
msa(test, type = "dna", method = "ClustalOmega") |> 
  msaConvert(type = "ape::DNAbin") |> ape::dist.gene() |> sqrt() |> sqrt() |> 
  ape::bionj() |> 
  plot()


check_alignment <- function(dna_seq){
  msa(dna_seq, type = "dna", method = "ClustalOmega") |> 
    msaConvert(type = "seqinr::alignment") |>
    (\(x){
      output <- do.call("cbind", x$seq |> str_split(""))
      colnames(output) <- x$nam
      return(output)
    })()
}

# 
# dna_mutation <- function(orig_seq, num_mutations, lineage = T){
#   
#   # The first generation sequences
#   del_line <- c(`Deletions = 0` = orig_seq)
#   sub_line <- c(`Substitutions = 0` = orig_seq)
#   ins_line <- c(`Insertions = 0` = orig_seq)
#   del_seq <- sub_seq <- ins_seq <- orig_seq <- str_split(orig_seq, "")[[1]]
#   
#   if(lineage){
#     # Mutations in a lineage are aggregative
#     # Under this case, (num_)mutations must increase 
#     if(any(diff(num_mutations) < 0)) print("Incorrect use of function.")
#     intervals <- num_mutations - c(0, lag(num_mutations)[-1])
#     mutation_locations <- sample(1:length(orig_seq), last(num_mutations), replace = F)
#     
#     # Some last generation sequences
#     sub_line.last <- sub_seq
#     sub_line.last[mutation_locations] <- mtte(sub_seq[mutation_locations])
#     ins_line.last <- insert_at(orig_seq, mutation_locations, sample(c("A", "C", "G", "T"), last(num_mutations), replace = T))
#     full_ins_helper <- insert_at(orig_seq, mutation_locations, 1:last(num_mutations)) 
#     
#     for(i in seq_along(num_mutations)){
#       # Deletion mutations
#       del_seq[mutation_locations[1:num_mutations[i]]] <- "-"
#       del_line <- c(del_line, paste(del_seq, collapse = ""))
#       names(del_line)[i+1] <- paste("Deletions =", num_mutations[i])
#       
#       # Insertion mutations
#       ins_seq <- ins_line.last[full_ins_helper %in% c("A", "C", "G", "T", 1:num_mutations[i])]
#       ins_line <- c(ins_line, paste(ins_seq, collapse = ""))
#       names(ins_line)[i+1] <- paste("Insertions =", num_mutations[i])
#       
#       # Substitution mutations
#       sub_seq[mutation_locations[1:num_mutations[i]]] <- sub_line.last[mutation_locations[1:num_mutations[i]]]
#       sub_line <- c(sub_line, paste(sub_seq, collapse = ""))
#       names(sub_line)[i+1] <- paste("Substitutions =", num_mutations[i])
#     }
#   }
#   else{
#     for(i in seq_along(num_mutations)){
#       mutation_locations <- sample(1:length(orig_seq), num_mutations[i], replace = F)
#       del_seq <- sub_seq <- ins_seq <- orig_seq
#       
#       # Deletion mutations
#       del_seq[mutation_locations] <- "-"
#       del_line <- c(del_line, paste(del_seq, collapse = ""))
#       names(del_line)[i + 1] <- paste("Deletions =", num_mutations[i])
#       
#       # Insertion mutations
#       ins_seq <- insert_at(ins_seq, mutation_locations, sample(c("A", "C", "G", "T"), length(mutation_locations), replace = T))
#       ins_line <- c(ins_line, paste(ins_seq, collapse = ""))
#       names(ins_line)[i + 1] <- paste("Insertions =", num_mutations[i])
#       
#       # Substitution mutations
#       sub_seq[mutation_locations] <- mtte(sub_seq[mutation_locations])
#       sub_line <- c(sub_line, paste(sub_seq, collapse = ""))
#       names(sub_line)[i + 1] <- paste("Substitutions =", num_mutations[i])
#     }
#   }
#   return(list(del_line, sub_line, ins_line))
# }
