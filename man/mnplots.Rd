\name{mnplots}
\alias{mnplots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{M-N Plots}
\description{
Generate M-N plots and select clusters based on them
}
\usage{
mnplots(unclesResult, MCs = 10, corner = c(0, 1),
    removedtype = 'abs', removedval = 1,
    Vmse = numeric(), mseCache = numeric(), doplot = FALSE, subplotdim = numeric(),
    subplotind = 1:MCs, minimiseDistance = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{unclesResult}{
The result of the "uncles" function as it is, or a
similarly structred list as explained here.

This argument is a list which must include two named
elements at least:

- unclesResult$X: A list of the datasets as matrices
or data frames. If a single dataset is provided, it
may be provided as it is rather than as a list containing
a single matrix or data frame.

- unclesResult$B: This can simply be a list of partitions, or
can be a multi-dimensional list array of partitions. Any of the
partitions, for example unclesResult$B[[i]], is a binary
partition of M rows representing M genes, and K columns
representing K clusters. The element unclesResult$B[[i]][j,k]
is the binary membership (either 1 or 0) of the (j)th gene in
the (k)th cluster as per the (i)th partition.

The partitions can have different numbers of clusters (K values)
represented by their columns, but they all have to have the same
numbers of genes, represented by rows. Moreover, the gene
represented by the (j)th row in one of the partitions has to be
the same gene represented by the (j)th row in all of the other
partitions. In other words, the rows of the partitions must be
aligned.

If this argument is passed as the output of the "uncles"
function, it will be a 4D list array of binary partitions with
the dimensions of:

(T)x(NBP1)x(NBP2)x(NKs).

where (T) is the number of the CoPaM final trials; (NBP1) is
the number of the different values of the parameter of the
binarisation technique if the UNCLES type is "A", and is the
number of the different values of the parameter of the positive
binarisation technique if the UNCLES type is "B"; (NBP2) is 1
if the UNCLES type is "A", and is the number of the different
values of the parameter of the negative binarisation technique
if UNCLES type is "B"; (NKs) is the number of the different
numbers of clusters (K values). For example: if 5 trials of the
final CoPaM were considered, UNCLES type "A" was used with a
DTB binarisation technique whose parameter delta ranges from
0.0 to 1.0 with steps of 0.1, and 4 different K values were
considered (e.g. K = 4, 8, 12, and 16), then the dimensions of
unclesResult$B will be (5x11x1x4).

Other optional named elements of the argument unclesResult
include:

- unclesResult$GDM: Gene-dataset logical matrix of M rows representing
M genes and L columns representing L datasets. A value of 1 in an
element of this matrix indicates that the corresponding gene is
found in the corresponding dataset, i.e. it is represented by some
probe(s) in that dataset.

Default: All ones (all considered genes are found in all of the datasets).

- unclesResult$params$type: The type of UNCLES, 'A' or 'B'. Default: 'A'.

- unclesResult$params$setsP: For UNCLES type 'B', these are the datasets
considered in the positive set of datasets. See the description of the
argument "setsP" of the function "uncles".

- unclesResult$params$setsN: For UNCLES type 'B', these are the datasets
considered in the negative set of datasets. See the description of the
argument "setsN" of the function "uncles".

- unclesResult$params$wsets: For L datasets, this is a vector of L
numeric values representing the relative weights of the datasets. The
vector does not have to be normalised as it will be normalised within
the "uncles" function. Valid examples for 5 datasets include:

wsets = c(0.2, 0.2, 0.2, 0.2, 0.2)

wsets = rep(1, 5)

wsets = c(4, 4, 4, 4, 4)

wsets = c(1, 2, 2, 0, 1)

wsets = c(0.2, 0.3, 0, 0.4, 0.4)

Note that the first three examples result in the same weighting, 
which is to treat all datasets equally. If the weight of a dataset
was set to zero, this implies excluding it of the analysis.

Default: numeric() # which will be read as equal weights for all
datasets.
}
  \item{MCs}{
The number of clusters to be selected by the M-N scatter plots
technique. This is also the number of iterations, as in each
iteration one cluster is selected.

Default: 10.
}
  \item{corner}{
The coordinates of point at the unity-normalised M-N plots
which is considered as the reference point from which the 
distance is measured for the points of all of the clusters.
Better clusters are those which are closer to this corner.
Its default is the top-left corner of the plot with the
coordinates of (0.0, 1.0). As the horizontal axis represents
the dispersion within the cluster and the vertical axis
represents the size of the cluster, clusters closer to that
top-left corner minimise dispersion while maximize their size.

If the reference was moved a bit towards the right on the
horizontal axis (e.g. to become at (0.2, 1.0)), wider clusters
will be selected. While if it was moved towards the left (e.g.
to become at (-0.2, 1.0)), tighter clusters will be selected.

Default: c(0.0, 1.0).
}
\item{removedtype}{
This is either 'perc' or 'abs'. When a cluster is selected as
the best cluster (closest to the "corner" argument), how do we
identify the other clusters which overlap with it?

Read the description of the argument "removedval" below for details.

Default: 'abs'
}
\item{removedval}{
A numeric value indicating the minimum amount of overlap between the
cluster selected as the best cluster in the current iteration and
the other clusters for these other clusters to be removed before the
following iteration.

If "removedtype" is 'perc', "removedval" represents the percentage
of the overlap out of the smaller cluster between the two clusters
being compared. "removedval" in this case should be in the range
0.0 to 1.0. For example, if "removedval" is 0.25, the overlap between
the two clusters has to be at least 25\% of the smaller cluster of
the two to consider it a significant overlap.

If "removedtype" is 'abs', "removedval" represents the minimum number
of genes in the overlap to consider it as a significant overlap.
"removedval" in this case should be an integer greater than zero.

Default: 1. As the default of "removedtype" is 'abs', this means that
if a single gene was found in the overlap, the overlap is considered
significant.
}
\item{Vmse}{
%%     ~~Describe \code{Vmse} here~~
}
\item{mseCache}{
%%     ~~Describe \code{mseCache} here~~
}
\item{doplot}{
If TRUE, the function plots the M-N plots in addition to providing
the calculated results in the output. If FALSE, it just calculates
the results and provides them without plotting.
}
\item{subplotdim}{
%%     ~~Describe \code{subplotdim} here~~
}
\item{subplotind}{
%%     ~~Describe \code{subplotind} here~~
}
\item{minimiseDistance}{
%%     ~~Describe \code{minimiseDistance} here~~
}
}
\examples{
# This is the simplist way to apply UNCLES and MN plots.
# Just pass the datasets to the "uncles" function and then pass
# the UNCLES result to the "mnplots" function.
# Both functions will use default values for all other arguments.
#
# Define three random gene expression datasets for 1000 genes.
# The number of samples in the datasets are 6, 4, and 9, respectively.
#
# X = list()
# X[[1]] = matrix(rnorm(6000), 1000, 6)
# X[[2]] = matrix(rnorm(4000), 1000, 4)
# X[[3]] = matrix(rnorm(9000), 1000, 9)
#
# unclesResult <- uncles(X)
# mnResult <- mnplots(unclesResult)
#
# The clusters will be available in the form of a partition matrix in the variable:
# mnResult$B;
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~M-N scatter plots }% use one of  RShowDoc("KEYWORDS")
