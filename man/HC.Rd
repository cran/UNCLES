\name{HC}
\alias{HC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{HC}
\description{
Perform hierarchical clustering
}
\usage{
HC(X, K, distancemetric = "euclidean", method = "ward.D2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
\item{K}{
%%     ~~Describe \code{K} here~~
}
\item{distancemetric}{
Any distance method which can be used by the function "dist"

Default: "euclidean"
}
\item{method}{
The hierarchical clustering method / algorithm. Acceptable
values are: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
"mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

Default: "ward.D2"
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
