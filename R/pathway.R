pathway.var.test <- function(exp1, exp2, dag, alpha) {
  l <- .procParams(exp1, exp2, dag)
  exp1 <- l$exp1
  exp2 <- l$exp2

  res           <- .runPathwayVarTest(exp1, exp2, l$graph, alpha)
  ns            <- nodes(res$graph)
  res$cli.moral <- lapply(res$cli.moral, function(ixs) ns[ixs])

  res
}

pathway.mean.test <- function(exp1, exp2, dag, alpha, perm.num=10000) {
  l <- .procParams(exp1, exp2, dag)
  exp1 <- l$exp1
  exp2 <- l$exp2

  exp1.num  <- nrow(exp1)
  exps      <- rbind(exp1, exp2)
  exps.num  <- nrow(exp1) + nrow(exp1)

  path.test <- .runPathwayVarTest(exp1, exp2, l$graph, alpha)
  cli.moral <- path.test$cli.moral
  check     <- path.test$check

  stat.obs  <- .hote(exp1, exp2, check, cli.moral)

  stat.perm <- vector("numeric", perm.num)
  for (i in seq_len(perm.num)) {
    ind          <- sample(exps.num)
    exp1.perm    <- exps[ind[1:exp1.num],]
    exp2.perm    <- exps[ind[(exp1.num+1):exps.num],]
    stat.perm[i] <- .hote(exp1.perm, exp2.perm, check, cli.moral)
  }

  list(alpha.obs=sum(stat.perm >= stat.obs) / perm.num, graph=path.test$graph)
}

.runPathwayVarTest <- function(exp1, exp2, graph, alpha) {
  cliques <- graph$cli.moral
  maxCliqueSize <- max(sapply(cliques, length))
  if (nrow(exp1) <= maxCliqueSize)
    stop("exp1 should have more than ", maxCliqueSize, " rows (samples)")
  else if (nrow(exp2) <= maxCliqueSize)
    stop("exp2 should have more than ", maxCliqueSize, " rows (samples)")

  cov <- .estimateCov(exp1, exp2)

  s1.hat <- qpIPF(cov$s1, cliques)
  s2.hat <- qpIPF(cov$s2, cliques)
  s.hat  <- qpIPF(cov$s,  cliques)

  k1.hat <- solve(s1.hat)
  k2.hat <- solve(s2.hat)
  k.hat  <- solve(s.hat)

  k1.det <- det(k1.hat)
  k2.det <- det(k2.hat)
  k.det  <- det(k.hat)

  lambda.obs  <- nrow(exp1)*log(k1.det/k.det) + nrow(exp2)*log(k2.det/k.det)
  arcs        <- (sum(graph$adj.moral)/2) + ncol(exp1)
  lambda.theo <- qchisq(0.95, arcs)
  alpha.obs   <- 2 * min(pchisq(lambda.obs, arcs), 1-pchisq(lambda.obs, arcs))
  check       <- alpha.obs <= alpha

  list(alpha.obs=alpha.obs, cli.moral=cliques, check=check, graph=graph$moral, lambda.obs=lambda.obs, lambda.theo=lambda.theo)
}
