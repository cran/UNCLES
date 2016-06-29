\name{fuzzystretch}
\alias{fuzzystretch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fuzzystretch}
\description{
Stretches the fuzzy values in each row of the matrix X such the one
stays one, zeros stays zero, x0 stays x0, values between zero and x0
decreased, and values between x0 and one increased.
   
If a vector was given for x0 then it should have elements with the same
number of rows in X. If a single value of x0 was given then it is used
for all of the rows of X.

The default value of x0 is the mean of the non-zero elements of the
corresponding row, we recommend using this value because it preserve
the criterion of (sum of fuzzy values for a single row is unity) to a
good level.
}
\usage{
fuzzystretch(X, x0 = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, x0 = -1) 
{
    return(NULL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
