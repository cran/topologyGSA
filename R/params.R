.procParams <- function(exp1, exp2, dag) {
  if (!is.matrix(exp1) && !is.data.frame(exp1))
    stop("exp1 is not a matrix or a data frame.")
  else if (!is.matrix(exp2) && !is.data.frame(exp2))
    stop("exp2 is not a matrix or a data frame.")
  else if (ncol(exp1) != ncol(exp2))
    stop("exp1 and exp2 differ in the number of columns (genes)")
  else if (any(colnames(exp1) != colnames(exp2)))
    stop("exp1 and exp2 differ in the column names (gene names)")
  else if (nrow(exp1) < 3)
    stop("exp1 should have at least 3 rows (samples)")
  else if (nrow(exp2) < 3)
    stop("exp2 should have at least 3 rows (samples)")

  common <- intersect(colnames(exp1), nodes(dag))
  if (length(common) < 3)
    stop("need at least 3 genes in common between expression and dag")

  exp1  <- exp1[,common,drop=FALSE]
  exp2  <- exp2[,common,drop=FALSE]
  dag   <- subGraph(common, dag)
  graph <- .processGraph(dag)
  
  list(exp1=exp1, exp2=exp2, graph=graph)
}
