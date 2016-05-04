\name{mnplots}
\alias{mnplots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{M-N Plots}
\description{
Generate M-N plots and select clusters based on them
}
\usage{
mnplots(X, B, type, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{\dots}{
}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, B, MCs = numeric(), type = "A", Vmse = numeric(), 
    fig = 0, corner = c(0, 1), setsP = seq(from = 1, to = floor(length(X)/2)), 
    setsN = seq(from = floor(length(X)/2) + 1, to = length(X)), 
    removedtype = "perc", removedval = .Machine$double.eps, Xtype = "data", 
    mseCache = numeric(), subplotdim = numeric(), wsets = rep(1, 
        length(X)), GDM = numeric()) 
{
    return(NULL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
