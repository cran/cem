\name{att}
\alias{att}
\title{Example of ATT estimation from CEM output}
\description{
An example of ATT estimation from CEM output
}
\usage{
att(obj, formula, data, model="lm", family="binomial")
}
\arguments{
  \item{obj}{a \code{cem.atch} or \code{multicem} object}
  \item{data}{a single data.frame or a list of data.frame's in case of \code{multicem}}
  \item{formula}{formula type specification of model. See Details.}
  \item{model}{either \code{lm} or \code{glm}. See Details.}
  \item{family}{used if model is \code{glm}, otherwise ignored.}
}
\details{
Argument \code{data} must be a single data frame or a list of (mulitply
imputed) data frames. 

Argument \code{model} can be \code{lm} or \code{glm} if the outcome
variable in the ATT estimation is, e.g., a binary outcome. If the
outcome is \code{y} and the treatment variable is \code{T}, then a
\code{formula} like \code{y ~ T} is enough to estimate the ATT: it is
just the coefficient of \code{T}. User can add covariates to span any
remaining imbalance after the match, such as \code{y ~ T + age + sex},
to adjust for variables \code{age} and \code{sex}.

In the case of multiply imputed datasets, the model is applied to each
single matched data and the ATT and is the standard error estimated
using the standard formulas for combining results of multiply imputed
data.  } 
\value{
  A matrix of estimates with their standard error, or a list in
  the case of \code{multicem}.  } 
\references{Stefano Iacus, Gary King, Giuseppe Porro, ``Matching for
  Casual Inference Without Balance Checking,''
  http://gking.harvard.edu/files/abs/cem-abs.shtml}
\author{Stefano Iacus, Gary King, and Giuseppe Porro}

\examples{ 
data(LL)

# cem match: automatic bin choice
mat <- cem(treatment="treated",data=LL, drop="re78")
mat
mat$k2k

# ATT estimate
att(mat, re78~treated,  data=LL)


# reduce the match into k2k using euclidean distance within cem strata
mat2 <- k2k(mat, LL, "euclidean", 1)
mat2
mat2$k2k

# ATT estimate after k2k
att(mat2, re78~treated, data=LL)

# example with missing data
# using multiply imputated data
# we use Amelia for multiple imputation
if(require(Amelia)){
 data(LL)
 n <- dim(LL)[1]
 k <- dim(LL)[2]

# we generate missing values in 30% of the rows of LL data
# randomly in one colum per row
 LL1 <- LL
 idx <- sample(1:n, .3*n)
 invisible(sapply(idx, function(x) LL1[x,sample(2:k,1)] <<- NA))


 imputed <- amelia(LL1)[1:5] 

 mat <- cem("treated", datalist=imputed, data=LL1, drop="re78")

 print(mat)
  
 att(mat, re78 ~ treated, data=imputed)
}
}
\keyword{multivariate}

