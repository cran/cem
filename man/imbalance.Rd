\name{imbalance}
\alias{imbalance}
\title{Calculates several imbalance measures}
\description{
Calculates several imbalance measures for the original
and matched data sets}
\usage{
imbalance(group, data, drop=NULL, breaks = NULL, weights)
}
\arguments{
 \item{group}{the group variable}
  \item{data}{the data}
  \item{drop}{a vector of variable names in the data frame to ignore}
  \item{breaks}{a list of vectors of cutpoints used to calculate the L1
  measure. See Details.}
  \item{weights}{weights}
}
\details{
This function calculates several imbalance measures.
For numeric variables, the difference in means (under the column
 \code{statistic}), the difference in quantiles and the L1 measure is
 calculated. For categorical variables the L1 measure and the
 Chi-squared distance (under column \code{statistic}) is calculated.
 Column \code{type} reports either \code{(diff)} or \code{(Chi2)} to
 indicate the type of statistic being calculated.

If \code{breaks} is not specified, the Scott automated bin calculation
is used (which coarsens less than Sturges, which used in
\code{\link{cem}}).  Please refer to \code{\link{cem}} help page. In
this case, breaks are used to calculate the L1 measure.

This function also calculate the global L1 imbalance measure. 
If \code{breaks} is missing, the default rule to calculate cutpoints
is the Scott's rule. 
See \code{\link{L1.meas}} help page for details.  
}
\value{
  An object of class \code{imbalance} which is a list with the following
  two elements
  \item{tab}{Table of imbalance measures}
  \item{L1}{The global L1 measure of imbalance}
}
\references{Stefano Iacus, Gary King, Giuseppe Porro, ``Matching for
  Casual Inference Without Balance Checking: Coarsened Exact Matching,''
  http://gking.harvard.edu/files/abs/cem-abs.shtml}
\author{Stefano Iacus, Gary King, and Giuseppe Porro}
\examples{ 

data(LL)

todrop <- c("treated","re78")
   
imbalance(LL$treated, LL, drop=todrop)

# cem match: automatic bin choice
mat <- cem(treatment="treated", data=LL, drop="re78")
}
\keyword{ datagen }
