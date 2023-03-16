library(transport)

first <- wpp(coordinates=cbind(1:3, 0), c(1, .5, 0))
second <- wpp(coordinates=cbind(1:3, 0), c(1, 0, .5))


temp <- beta_cg |> hist_paired(bin_count = 100, return_bins = T) 
first <- wpp(coordinates = temp[[2]], mass = temp[[1]][1,])
second <- wpp(coordinates = temp[[2]], mass = temp[[1]][4,])
my_wasserstein(first, second)

my_wasserstein <- function (a, b, p = 1, tplan = NULL, costm = NULL, prob = TRUE, all_eqdist = F, ...) 
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
  if(all_eqdist) dd <- rep(1, length(dd))
  
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
  return(res)
}

