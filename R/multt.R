.mult.test <- function(exp1, exp2, perm.num) {
  exp1.num <- nrow(exp1)
  exps     <- rbind(exp1, exp2)
  exps.num <- nrow(exps)

  t.obs <- .hote(exp1, exp2, FALSE)

  stat.perm <- vector("numeric", perm.num)
  for (i in 1:perm.num) {
    ind          <- sample(exps.num)
    exp1.perm    <- exps[ind[1:exp1.num],]
    exp2.perm    <- exps[ind[(exp1.num+1):exps.num],]
    stat.perm[i] <- .hote(exp1.perm, exp2.perm, FALSE)
  }

  alpha.obs <- sum(stat.perm >= t.obs) / perm.num
  list(alpha.obs=alpha.obs, t.obs=t.obs)
}
