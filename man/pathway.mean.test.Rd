\name{pathway.mean.test}
\alias{pathway.mean.test}
\title{Testing the means of a pathway.}

\description{Test the equality of means in two experimental conditions for a pathway,
conditioning on the result of the test on the concentration matrices
}
\usage{pathway.mean.test(exp1,exp2,dag,alpha,perm.num=10000)}

\arguments{
\item{exp1}{a matrix with n1 individuals (rows) in the first experimental condition and p genes (columns).}
\item{exp2}{a matrix with n2 individuals (rows) in the second experimental condition and p genes (columns).
The genes in the two experimental conditions must be the same.}
\item{dag}{graphNEL object, a directed acyclic graph (DAG) corresponding to the pathway of interest. See package \code{gRbase} for more details.}
\item{perm.num}{number of permutations of the test.}
\item{alpha}{significance level of the test.}
}
\details{The function tests  the equality of the means of  two
  experimental conditions, conditioning on the result of
  \code{\link{pathway.var.test}}.

  The expression data may contain some genes differing
  from those in the pathway: in such case the function automatically takes
  the intersection between the two gene sets.

  This function requires \code{gRBase} and \code{qpgraph} packages.
}

\value{

\item{alpha.obs}{P-value of the test}
\item{graph}{A graphNEL object reprenting the moral graph.}

}


\author{M. Sofia Massa,  Gabriele Sales}
\references{Massa, M.S., Chiogna, M., Romualdi, C. (2010). Gene set analysis exploiting the topology of a pathway. BMC Systems Biology, 4:121
\url{http://www.biomedcentral.com/1752-0509/4/121}
}
\seealso{\code{\link{pathway.var.test}}, \code{\link{clique.var.test}},
\code{\link{clique.mean.test}},
}
\examples{
data(examples)

pathway.mean.test(exp1, exp2, dag_bcell, 0.05, 100)
}
