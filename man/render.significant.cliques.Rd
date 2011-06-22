\name{render.significant.cliques}
\alias{render.significant.cliques}
\title{Graphical display of the significant cliques}

\description{This function highlights the significant cliques on the moralized and triangulated graph. Different colors correspond to different values of -log(pvalue)
}
\usage{
render.significant.cliques(info)
}

\arguments{
\item{info}{the result of \code{\link{clique.var.test}} or
  \code{\link{clique.mean.test}}.
  }
}

\references{Massa, M.S., Chiogna, M., Romualdi, C. (2010). Gene set
  analysis exploiting the topology of a pathway. BMC Systems Biology,
  4:121
  \url{http://www.biomedcentral.com/1752-0509/4/121}
}
\seealso{\code{\link{clique.var.test}}, \code{\link{clique.mean.test}}.
}
\examples{
data(examples)

info <- clique.var.test(exp1, exp2, dag_bcell, 0.1)
\dontrun{render.significant.cliques(info)}
}
