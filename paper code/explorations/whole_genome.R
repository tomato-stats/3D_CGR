library(Rcpp)
sourceCpp("./R/chaosgame.cpp")

library(seqinr)

sars <- read.fasta("./paper code/data/sarscov2.fasta")
sars <- lapply(sars, function(x) x[x!="-"] |> toupper())

sars_seq <- lapply(sars |> chaosGame(baseCoords = as.matrix(coord1[,1:3])))

temp_sig <- coordinate_signature(sars_seq, bin_count = 6110) 
rownames(temp_sig) <- apply(sars_info, 1, paste0, collapse = "")
  temp_sig |> dist() |> hclust(method = "complete") |>
  plot(axes = F, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, ann = F)

library(kmer)

# Function to assess speeds for different number of lengths 
  
time_kmer <- function(seq, ...){
  times <- system.time(kcount(seq, ...))
  data.frame(user.self = times[1], sys.self = times[2], elapsed = times[3], 
    user.child = times[4], sys.child = times[5])
}

test_kmer <- function(test_sequence){
  lengths <- 10^(2:5)
  output <- data.frame()
  for(i in lengths){
    shortened_sequence <- lapply(sars, function(x) x[1:i])
    append_me <- data.frame(length = i, time_kmer(shortened_sequence, ...))
    output <- rbind(output, append_me)
  }
  return(output)
}
kmer_times <- test_kmer(sars)
kmer_times |> pivot_longer(cols = -1) |> filter(!is.na(value)) |> ggplot(aes(x = length, y = value )) + geom_point() + geom_line() + facet_wrap(~name)

build_features <- function(seq, ...){
  seq_cg <- lapply(seq |> chaosGame(baseCoords = as.matrix(coord1[,1:3])))
  seq_vdist <- distance_from_vertices(seq_cg)
  coordinate_signature(seq_vdist, ...)
}

time_cg <- function(){
  times <- system.time(build_features(seq, ...))
  data.frame(user.self = times[1], sys.self = times[2], elapsed = times[3], 
             user.child = times[4], sys.child = times[5])
}

test_cg <- function(test_sequence){
  lengths <- 10^(2:5)
  output <- data.frame()
  for(i in lengths){
    shortened_sequence <- lapply(sars, function(x) x[1:i])
    append_me <- data.frame(length = i, time_cg(shortened_sequence, ...))
    output <- rbind(output, append_me)
  }
  return(output)
}