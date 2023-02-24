freq_to_prob <- function(input_hist){
  # Input is a matrix where each row is an organism
  max_sequence_length <- rowSums(input_hist) |> max()
  sequence_length_deficit <- max_sequence_length - rowSums(input_hist) 
  add_these_to_all <- matrix(sequence_length_deficit / max_sequence_length, ncol= ncol(input_hist), nrow = nrow(input_hist)) 
  input_hist <- input_hist + add_these_to_all
  output <- input_hist / rowSums(input_hist)
  rownames(output) <- rownames(input_hist)
  return(output)
}

modify_dist <- function(input_hist, DISTFUNC){
  output <- matrix(0, ncol = nrow(input_hist), nrow = nrow(input_hist))
  for(i in 1:(nrow(input_hist)-1)){
    for(j in (i + 1):nrow(input_hist)){
      output[i, j] <- DISTFUNC(input_hist[i,,drop = F], input_hist[j,,drop = F])
    }
  }
  rownames(output) <- rownames(input_hist)
  colnames(output) <- rownames(input_hist)
  return(t(output) + output)
}



emdist:emd2d(freq_to_prob)
philentrophy::KL
test

kl <- function(x){
  output <- philentropy::KL(x)
  output <- (log(output + .1) - log(1.1)) / (log(0.1) - log(1.1))
  rownames(output) <- rownames(x)
  colnames(output) <- rownames(x)
  return(output)
}
