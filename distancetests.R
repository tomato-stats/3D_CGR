freq_to_prob <- function(input_hist){
  # Input is a matrix where each row is an organism
  max_sequence_length <- rowSums(input_hist) |> max()
  sequence_length_deficit <- rowSums(input_hist) - max_sequence_length
  add_these_to_all <- matrix(sequence_length_deficit / max_sequence_length, ncol= ncol(input_hist), nrow = nrow(input_hist)) 
  input_hist <- input_hist + add_these_to_all
  input_hist / rowSums(input_hist)
}

modify_dist <- function(){
  
  
}



emdist:emd2d(freq_to_prob)
philentrophy::KL
test

test