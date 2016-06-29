fuzzystretch <- function(X, x0 = -1) {
    M = nrow(X)
    N = ncol(X)
    if (x0 == -1) {
        x0 = rep(0, M)
        for (i in 1:M) {
            xrow = X[i,]
            x0[i] = mean(xrow[xrow > 0])
        }
        x0[x0 == 1] = 0.5
    } else {
        if (length(x0) == 1) {
            x0 = rep(x0, M)
        } else if (length(x0) != M) {
            stop("x0 must be a single value or a vector with elements equal in number to the number of rows of X")
        }
    }

    y = matrix(, M, N)
    for (i in 1:M) {
        xrow = X[i,]
        xt = xrow
        xt[xrow < x0[i]] = (pi * xrow[xrow < x0[i]]) / (2 * x0[i]) - pi / 2
        xt[xrow >= x0[i]] = (xrow[xrow >= x0[i]] - x0[i]) * pi / (2 * (1 - x0[i]))

        yt = rep(0, N)
        yt[xrow < x0[i]] = x0[i] + x0[i] * sin(xt[xrow < x0[i]])
        yt[xrow >= x0[i]] = x0[i] + (1 - x0[i]) * sin(xt[xrow >= x0[i]])

        y[i,] = yt
    }

    return(y)
}