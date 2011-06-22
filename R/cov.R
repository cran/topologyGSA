.estimateCov <- function(exp1, exp2) {
  exp1.num <- nrow(exp1)
  exp2.num <- nrow(exp2)
  gene.num <- ncol(exp1)

  exp1.mean        <- apply(exp1, 2, mean)
  names(exp1.mean) <- NULL
  exp1.scal        <- exp1 - matrix(rep(exp1.mean, exp1.num*gene.num), exp1.num, gene.num, byrow=T)
  s1               <- cov(exp1.scal)

  exp2.mean        <- apply(exp2, 2, mean)
  names(exp2.mean) <- NULL
  exp2.scal        <- exp2 - matrix(rep(exp2.mean, exp2.num*gene.num), exp2.num, gene.num, byrow=T)
  s2               <- cov(exp2.scal)

  s <- (s1*(exp1.num-1)+ s2*(exp2.num-1)) / (exp1.num+exp2.num-2)
  list(s1=s1, s2=s2, s=s)
}
