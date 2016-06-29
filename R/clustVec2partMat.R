clustVec2partMat <- function(C) {
    labels = unique(C[C != 0])
    if (length(labels) == 0) {
        return(rep(0, length(C)))
    }

    K = max(labels)
    N = length(C)
    U = matrix(FALSE, K, N)
    for (i in 1:N) {
        U[labels[labels == C[i]], i] = TRUE
    }

    if (sum(rowSums(U) == 0) > 0) {
        U = U[ - which(rowSums(U) == 0),]
    }
    return(U)
}