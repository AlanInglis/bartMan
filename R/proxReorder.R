proxReorder <- function(d) {

  prox <- as.dist(d)
  score <- apply(as.matrix(prox), 1, max)
  o <- DendSer::dser(-prox, -score, cost = DendSer::costLS)
  res <- d[o, o]
  class(res)<- class(d)
  res

}

