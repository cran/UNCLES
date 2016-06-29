isValidBPM <- function(U) {
    s1 = colSums(U)
    s2 = rowSums(U)

    if (any(U != 1 && U != 0)) {
        return(FALSE)
    }

    if (max(s1) != 1 || max(s2) != 1) {
        return(FALSE)
    }

    if (min(s2) <= 0) {
        return(FALSE)
    }

    return(TRUE)
}