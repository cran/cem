\name{L1.profile}
\alias{L1.profile}
\title{Calculates L1 distance for different coarsenings}
\description{
Calculates L1 distance for different coarsenings
}
\usage{
L1.profile(group, data, drop = NULL, min.cut=2, max.cut=6, ref.cut=3, 
 weights, plot=TRUE, add=FALSE, col="red", lty=1)
}
\arguments{
  \item{group}{the group variable}
  \item{data}{the data}
  \item{drop}{a vector of variable names in the data frame to ignore}
  \item{min.cut}{minimum number of cut points per variable}
  \item{max.cut}{maximum number of cut points per variable}
  \item{ref.cut}{reference number of cut points per variable}
  \item{weights}{weights}
  \item{plot}{plot a graph?}
  \item{add}{add graph to an existing plot? Makes sense only if \code{plot} is \code{TRUE}}
  \item{col}{draw in specified color}
  \item{lty}{draw using specified lty}
}
\details{
  The L1 measure depends on the coarsening chosen to calculate it, and
  as such the comparison of different matching solutions may differ
  depending on this somewhat arbitrary choice.  This function computes
  L1 for a very wide range of possible
  coarsenings.  The point of this function is that if
  one matching solution has a lower L1 than another, then it dominates
  without regard to the choice of coarsening.  A graphic display conveys
  the results succinctly.  (The logic is similar to
  that for ROC curves used for classification algorithms.)
  (This degree of coarsening should remain fixed for different CEM
  runs.)

This function first tries to move the cutpoints of the support of all
variables together from \code{min.cut} to \code{max.cut}. Then, fixing a
coarsening into \code{ref.cut} number of cut points on all variables, it
changes the coarsening one variable at time, again from \code{min.cut}
to \code{max.cut}.
 
Non numeric variables are not grouped, i.e. no coarsening occurs.

A \code{plot} method exists for the returned object.
}
\value{
  An invisible object of class \code{L1profile} which contains a named list of coarsening.
}
\references{Stefano Iacus, Gary King, Giuseppe Porro, ``Matching for
  Casual Inference Without Balance Checking: Coarsened Exact Matching,''
  http://gking.harvard.edu/files/abs/cem-abs.shtml}
\author{Stefano Iacus, Gary King, and Giuseppe Porro}
\examples{
data(LL)
for(i in c(4:6,10:12))
 LL[[i]] <- factor(LL[[i]])

L1.profile(LL$treated,LL, drop=c("treated","re78"), max.cut=10)

if(require(MatchIt)){
 m2 <- matchit(treated ~ black + hispanic + married + nodegree + u74 + u75 + education +
  age + re74 + re75, data=LL, distance="logit")
 
 L1.profile(LL$treated,LL, drop=c("treated","re78"), max.cut=10, 
  weights=m2$w, add=TRUE, col="green", lty=2)
}

m1 <- cem("treated", LL, drop="re78")

L1.profile(LL$treated,LL, drop=c("treated","re78"), max.cut=10, 
 weights=m1$w, add=TRUE, col="blue", lty=3)

legend(5, 0.9, legend=c("raw data", "pscore", "cem"), lty=1:3, col=c("red", "green", "blue"))
}
\keyword{ datagen }
