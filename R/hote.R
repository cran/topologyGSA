.hote <- function(exp1, exp2, exact, cliques=NULL) {
  exp1.num <- nrow(exp1)
  exp2.num <- nrow(exp2)
  gene.num <- ncol(exp1)

  exp1.bar <- colMeans(exp1)
  exp2.bar <- colMeans(exp2)
  exp1.s   <- cov(exp1)
  exp2.s   <- cov(exp2)

  exp.diff <- exp1.bar - exp2.bar
  if (!is.null(cliques)) {
    exp1.s <- qpIPF(exp1.s, cliques)
    exp2.s <- qpIPF(exp2.s, cliques)
  }

  if (exact) {
    s  <- ((exp1.num-1)*exp1.s + (exp2.num-1)*exp2.s) / (exp1.num + exp2.num - 2)
    t2 <- ((exp1.num*exp2.num) / (exp1.num+exp2.num)) * (exp.diff %*% solve(s) %*% exp.diff)

    if (is.null(cliques)) {
      c <- exp1.num + exp2.num - gene.num - 1
      t.obs <- t2 * c / (gene.num * (exp1.num + exp2.num - 2))
      alpha.obs <- 1 - pf(t.obs, gene.num, c)
      list(alpha.obs=alpha.obs, t.obs=t.obs)
    } else {
      as.numeric(t2)
    }
  } else {
    s <- exp1.s/exp1.num + exp2.s/exp2.num
    as.numeric(exp.diff %*% solve(s) %*% exp.diff)
  }
}
