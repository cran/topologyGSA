clique.var.test <- function(exp1, exp2, dag, alpha) {
  l <- .checkParams(exp1, exp2, dag)
  exp1 <- l$exp1
  exp2 <- l$exp2
  dag  <- l$dag
  
  .runCliqueVarTest(exp1, exp2, dag, alpha)
}

clique.mean.test <- function(exp1, exp2, dag, alpha, perm.num=1000) {
  l <- .checkParams(exp1, exp2, dag)
  exp1 <- l$exp1
  exp2 <- l$exp2
  dag  <- l$dag

  exp.all    <- rbind(exp1, exp2)
  cli.test   <- .runCliqueVarTest(exp1, exp2, dag, alpha)
  check      <- cli.test$check
  cliques    <- cli.test$cliques
  clique.num <- length(cliques)

  alpha.obs  <- vector("numeric", clique.num)
  t.obs      <- vector("numeric", clique.num)
  for (i in 1:clique.num) {
    cli      <- unlist(cliques[i])
    exp1.cli <- exp1[,cli]
    exp2.cli <- exp2[,cli]

    
    if (length(cli) != 1) {
      if (check[i])
        r <- .mult.test(exp1.cli, exp2.cli, perm.num)
      else
        r <- .hote(exp1.cli, exp2.cli, TRUE)

      alpha.obs[i] <- r$alpha.obs
      t.obs[i]     <- r$t.obs
    } else {
      r            <- t.test(exp1.cli, exp2.cli)
      alpha.obs[i] <- r$p.value
      t.obs[i]     <- r$statistic
    }
  }

  list(alpha.obs=as.numeric(alpha.obs), cliques=cliques, check=check, graph=cli.test$graph, t.obs=t.obs)
}

.runCliqueVarTest <- function(exp1, exp2, dag, alpha) {
  cov <- .estimateCov(exp1, exp2)

  pg         <- .processGraph(dag)
  cliques    <- pg$cli.tg$maxCliques
  clique.num <- length(cliques)

  alpha.obs  <- rep(0,     clique.num)
  lambda.obs <- rep(0,     clique.num)
  check      <- rep(FALSE, clique.num)

  for (i in 1:clique.num) {
    cli <- unlist(cliques[i])
    p   <- length(cli)

    s1.hat <- cov$s1[cli, cli]
    s2.hat <- cov$s2[cli, cli]
    s.hat  <- cov$s[cli, cli]

    if (p == 1) {
      s1.det <- s1.hat
      s2.det <- s2.hat
      s.det  <- s.hat
    } else {
      s1.det <- det(s1.hat)
      s2.det <- det(s2.hat)
      s.det  <- det(s.hat)
    }

    lambda.obs[i] <- nrow(exp1)*log(s.det/s1.det) + nrow(exp2)*log(s.det/s2.det)
    alpha.obs[i]  <- 2 * min(pchisq(lambda.obs[i], p*(p+1)/2), 1-pchisq(lambda.obs[i], p*(p+1)/2))

    if (alpha.obs[i] <= alpha)
      check[i] <- TRUE
  }

  list(alpha.obs=alpha.obs, cliques=cliques, check=check, graph=pg$tg, lambda.obs=lambda.obs)
}
