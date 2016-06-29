clustDist <- function(U1, U2, X = matrix(, 0, 0), criterion = "direct_euc") {
    if (criterion == "direct_euc") {
        return(as.matrix(pdist::pdist(U1, U2)))
    } else if (criterion == "centers_euc") {
        # find the centers
        N = ncol(X)
        centres1 = (U1 %*% X) / (rowSums(U1) %*% rep(1, N))
        centres2 = (U2 %*% X) / (rowSums(U2) %*% rep(1, N))

        # find the distances between the centres
        D = as.matrix(pdist::pdist(centres1, centres2))
        m = max(D[!is.na(D)])
        if (!isempty(m)) {
            D[is.na(D)] = m + 1
        }
        return(D)
    } else if (criterion == "union_std") {
        N = ncol(X)
        K1 = nrow(U1)
        K2 = nrow(U2)
        D = matrix(, K1, K2)

        for (i in 1:K1) {
            for (j in 1:K2) {
                uU = pmax(U1[i,], U2[j,])
                uCentre = (uU %*% X) / (rowSums(uU) %*% rep(1, N))
                dists = as.matrix(pdist::pdist(X, uCentre))
                D[i, j] = (t(dists) %*% uU) / sum(uU)
            }
        }
        return(D)
    } else {
        stop("Unsupported distances")
    }
}