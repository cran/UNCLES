KA <- function(X, K, distancemetric = "euclidean") {
    M = nrow(X)
    Dists = as.matrix(dist(X, method = distancemetric))

    ResultInd = rep(0, K)
    Xmean = colMeans(X)
    Dmean = as.matrix(pdist::pdist(X, Xmean))

    ResultInd[1] = which(Dmean == min(Dmean))[1]

    for (k in 1:(K - 1)) {
        D = apply(as.matrix(Dists[, ResultInd[1:k]]), 1, min)
        C = rep(0, M)
        for (i in 1:M) {
            if (all(ResultInd == i)) {
                next
            }
            C[i] = sum(pmax(D - Dists[, i], 0))
        }
        ResultInd[k + 1] = which(C == max(C))[1]
    }
    return(X[ResultInd,])
}