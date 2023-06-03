library(transport)

first <- wpp(coordinates=cbind(1:3, 0), c(1, .5, 0))
second <- wpp(coordinates=cbind(1:3, 0), c(1, 0, .5))


temp <- beta_cg |> hist_paired(bin_count = 100, return_bins = T) 
first <- wpp(coordinates = temp[[2]], mass = temp[[1]][1,])
second <- wpp(coordinates = temp[[2]], mass = temp[[1]][4,])
my_wasserstein(first, second)

my_wasserstein <- function (a, b, p = 1, tplan = NULL, costm = NULL, prob = TRUE, all_eqdist = F, normalize = F, ...) 
{
  default <- FALSE
  if (is.null(tplan)) {
    argus <- list(...)
    argus$fullreturn <- FALSE
    if (default) {
      allargus <- c(list(a = a, b = b, costm = costm^p), 
                    argus)
      tplan <- do.call(transport, allargus)
    }
    else {
      allargus <- c(list(a = a, b = b, p = p), argus)
      tplan <- do.call(transport, allargus)
    }
  }
  K <- dim(tplan)[1]
  if (K == 0) {
    return(0)
  }
  wpsum <- function(x, m = rep(1, length(x)), pp) {
    mmax <- max(m)
    xmax <- max(abs(x))
    if (mmax == 0 || xmax == 0) {
      return(0)
    }
    if (pp == 1) {
      return(mmax * xmax * sum((m/mmax) * (x/xmax)))
    }
    else if (length(unique(m)) == 1 && unique(m) == 1) {
      return(xmax * sum((x/xmax)^pp)^(1/pp))
    }
    else {
      return((mmax)^(1/pp) * xmax * sum((m/mmax) * (x/xmax)^pp)^(1/pp))
    }
  }
  if (default) {
    dd <- costm[cbind(tplan$from, tplan$to)]
  }
  else {
    if (is(a, "pgrid") && is(b, "pgrid")) {
      gg <- expand.grid(a$generator)
      orig <- gg[tplan$from, , drop = FALSE]
      dest <- gg[tplan$to, , drop = FALSE]
    }
    else if (is(a, "pp") && is(b, "pp")) {
      orig <- a$coordinates[tplan$from, , drop = FALSE]
      dest <- b$coordinates[tplan$to, , drop = FALSE]
    }
    else if (is(a, "wpp") && is(b, "wpp")) {
      orig <- a$coordinates[tplan$from, , drop = FALSE]
      dest <- b$coordinates[tplan$to, , drop = FALSE]
    }
    else {
      stop("a and b must be both of the same class among 'pgrid', 'pp', 'wpp'")
    }
    dd <- apply(orig - dest, 1, wpsum, pp = 2)
  }
  if(all_eqdist) dd[which(dd !=0)] <- .1#dd <- rep(1, length(dd))
  # Make each deletion the average substitution mutation cost 
  else{ 
    dd[which(tplan$from == 1)] <- .1#coalesce(mean(dd[intersect(which(tplan$from != 1), which(dd!=0))]), .2)
    if(length(intersect(which(tplan$from != 1), which(dd!=0))) > 0)
      dd[intersect(which(tplan$from != 1), which(dd!=0))] <- median(dd[intersect(which(tplan$from != 1), which(dd!=0))])
  }
  #print(mean(dd[intersect(which(tplan$from != 1), which(dd!=0))]))
  #print(median(dd[intersect(which(tplan$from != 1), which(dd!=0))]))
  res <- wpsum(dd, tplan$mass, p)
  if (prob) {
    if (default) {
      res <- res/(sum(a)^(1/p))
    }
    else if (is(a, "pgrid") || is(a, "wpp")) {
      res <- res/(a$totmass^(1/p))
    }
    else {
      res <- res/(a$N^(1/p))
    }
  }
  if(normalize){
    a <- log(1.1)
    b <- log(0.1) - a
    res <- (log(0.1 + (1-res)) - a) / (b);
  }
  return(res)
}

all_wasserstein <- function(){
  temp <- beta_cg |> hist_paired(bin_count = 1000, return_bins = T) 
  all_lengths <- temp[[1]] |> rowSums() 
  max_length <- max(all_lengths)
  
  distances <- matrix(0, ncol = length(all_lengths), nrow = length(all_lengths))
  for(i in 1:(length(all_lengths)-1)){
    for(j in (i + 1):length(all_lengths)){
      first  <- wpp(coordinates = rbind(temp[[2]], 0), mass = c(temp[[1]][i,], max_length - all_lengths[i]))
      second <- wpp(coordinates = rbind(temp[[2]], 0), mass = c(temp[[1]][j,], max_length - all_lengths[j]))
      distances[i, j] <- my_wasserstein(first, second, all_eqdist = F, prob = F)
    }
  }
  distances <- t(distances) + distances
  rownames(distances) <- names(beta_seq)
}

all_wasserstein_cg <- function(list_coords){
  # Each coordinate gets count = 1
  all_lengths <- sapply(list_coords, nrow)
  max_length <- max(all_lengths)
  
  distances <- matrix(0, ncol = length(all_lengths), nrow = length(all_lengths))
  for(i in 1:(length(all_lengths)-1)){
    for(j in (i + 1):length(all_lengths)){
      first  <- wpp(coordinates = rbind(list_coords[[i]], 0), mass = c(rep(1, nrow(list_coords[[i]])), max_length - all_lengths[i]))
      second <- wpp(coordinates = rbind(list_coords[[j]], 0), mass = c(rep(1, nrow(list_coords[[j]])), max_length - all_lengths[j]))
      distances[i, j] <- my_wasserstein(first, second, all_eqdist = F, prob = F)
    }
  }
  distances <- t(distances) + distances
  rownames(distances) <- names(list_coords)
  return(distances)
}






temp <- nadh_cg
all_lengths <- sapply(temp, nrow) - 1
max_length <- max(all_lengths)

distances <- matrix(0, ncol = length(all_lengths), nrow = length(all_lengths))
for(i in 1:(length(all_lengths)-1)){
  for(j in (i + 1):length(all_lengths)){
    first  <- wpp(coordinates = temp[[i]], mass = c(max_length - all_lengths[i], rep(1, all_lengths[i])) )
    second <- wpp(coordinates = temp[[j]], mass = c(max_length - all_lengths[j], rep(1, all_lengths[j])) )
    distances[i, j] <- my_wasserstein(first, second, all_eqdist = F, normalize = T)
  }
}
distances <- t(distances) + distances
rownames(distances) <- names(nadh)
