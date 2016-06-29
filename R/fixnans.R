fixnans <- function(X, type = 'spline') {
    fixrow <- function(rowin, type) {
        #Length of the row
        M = length(rowin)

        #If one not-nan value exists, make it a constant vector with it, all nans, return zeros
        if (sum(is.na(rowin)) == M - 1) {
            return(rep(rowin[!is.na(rowin)], M));
        } else if (sum(is.na(rowin)) == M) {
            return(rep(0, M));
        }

        #For now, this is the case
        rowout = rowin

        #Generating binary vector locating the known points
        known = !is.na(rowin)

        #Extracting 2 time series from (t) for the known & the unknown points
        tknown = which(known)
        tunknown = which(!known)

        #If no unknowns are in the row, return it as it is
        if (length(tunknown) == 0) {
            return(rowin)
        }

        #Extracting the known x values
        xknown = rowin[known]

        #Using the interpolation methods to calculate the unknown points
        if (type == "spline") {
            xunknown = spline(tknown, xknown, xout = tunknown)$y
        } else {
            stop("Unsupported type of interpolation")
        }

        rowout[tunknown] = xunknown

        return(rowout)
    }

    if (is.vector(X)) {
        return(fixrow(X, type))
    }

    N = nrow(X)
    M = ncol(X)
    Xout = matrix(0, N, M)

    for (i in 1:N) {
        Xout[i,] = fixrow(X[i,], type)
    }

    return(Xout)
}
