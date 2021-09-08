proxReorder <- function(d) {

  prox <- as.dist(d)
  rprox <- range(prox)
  if (rprox[2] != rprox[1]) {
    prox <- (prox - rprox[1]) / (rprox[2] - rprox[1])
  }
  score <- apply(as.matrix(prox), 1, max)
  o <- DendSer::dser(-prox, -score, cost = DendSer::costLS)
  res <- d[o, o]
  class(res)<- class(d)
  res
}

