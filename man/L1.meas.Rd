\name{L1.meas}
\alias{L1.meas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates L1 distance between multidimensional histograms}
\description{
Evaluates L1 distance between multidimensional histograms
}
\usage{
L1.meas(group, data, drop=NULL, breaks = NULL, weights)
}
\arguments{
  \item{group}{the group variable}
  \item{data}{the data}
  \item{drop}{a vector of variable names in the data frame to ignore}
  \item{breaks}{a list of vectors of cutpoints; if not specified,
  automatic choice will be made}
  \item{weights}{weights}
}
\details{
This function calculates the L1 distance on the k-dimensional histogram.

If \code{breaks} is not specified, the Scott automated bin calculation
is used (which coarsens less than Sturges, which used in
\code{\link{cem}}).  Please refer to \code{\link{cem}} help page. In
this case, breaks are used to calculate the L1 measure.

If \code{breaks} is missing, the default rule to calculate cutpoints
is the Scott's rule. 

This code also calculate the Local Common Support (LCS) measure, which is
the proportion of non empty k-dimensional cells of the histogram which
contain at least one observation per group.
}
\value{
  An object of class \code{L1.meas} which is a list with the following fields
  \item{L1}{The numerical value of the L1 measure}
  \item{breaks}{A list of cutpoints used to calculate the L1 measure}
  \item{LCS}{The numerical value of the Local Common Support proportion}
}
\references{Stefano Iacus, Gary King, Giuseppe Porro, ``Matching for
  Casual Inference Without Balance Checking,''
  http://gking.harvard.edu/files/abs/cem-abs.shtml}
\author{Stefano Iacus, Gary King, and Giuseppe Porro}
\examples{
data(LL)
L1.meas(LL$treated,LL, drop=c("treated","re78"))
}
\keyword{ datagen }
