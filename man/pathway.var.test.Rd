\name{pathway.var.test}
\alias{pathway.var.test}
\title{Testing the concentration matrices of a pathway.}

\description{Test the equality of concentration matrices in two experimental conditions for a pathway
}
\usage{
pathway.var.test(exp1,exp2,dag,alpha)
}

\arguments{
\item{exp1}{a matrix with n1 individuals (rows) in the first experimental condition and p genes (columns).}
\item{exp2}{a matrix with n2 individuals (rows) in the second experimental condition and p genes (columns).
The genes in the two experimental conditions must be the same.}
\item{dag}{graphNEL object, directed acyclic graph (DAG) corresponding to the pathway of interest. See package \code{gRbase} for more details.}
\item{alpha}{significance level of the test.}
}
\details{The graph of a pathway is first converted into a DAG and then
  into a moral graph. The data is modelled with two Gaussian graphical
  models with zero mean and graph provided by the moral graph. The
  function tests the equality of the two concentration matrices
  (inverse of the covariance matrices).

  The expression data may contain some genes differing
  from those in the pathway: in such case the function automatically takes
  the intersection between the two gene sets.

  This function requires \code{gRBase} and \code{qpgraph} packages.
}
\value{
\item{lambda.obs}{the observed value of the test statistic.}
\item{lambda.theo}{the theoretical value of the test statistic.}
\item{alpha.obs}{the observed value of the significance level of the test.}
\item{cli.moral}{the cliques of the moral graph.}
\item{check}{logical, a vector with a TRUE for a significant pathway and
  a FALSE otherwise.}
\item{graph}{a graphNEL object representing the moral graph.}
}

\author{M. Sofia Massa, Gabriele Sales}

\references{Massa, M.S., Chiogna, M., Romualdi, C. (2010). Gene set
  analysis exploiting the topology of a pathway. BMC Systems Biology,
  4:121
  \url{http://www.biomedcentral.com/1752-0509/4/121}

Lauritzen, S.L. (1996). Graphical models. Clarendon Press, Oxford.
}
\seealso{\code{\link{pathway.mean.test}}, \code{\link{clique.var.test}},
\code{\link{clique.mean.test}}.
}
\examples{
data(examples)

pathway.var.test(exp1, exp2, dag_bcell, 0.05)
}
