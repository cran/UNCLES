binarise <- function(U, K, technique = "DTB", parameter = 0.0) {
    technique = tolower(technique)
    if (technique == "union" || technique == "ub") {
        return(U > 0)
    } else if (technique == "intersection" || technique == "ib") {
        return(U == 1)
    } else if (technique == "max" || technique == "mvb") {
        return(U == matrix(1, K, 1) %*% apply(U, 2, max))
    } else if (technique == "valuethresh" || technique == "value" || technique == "vtb") {
        return(U >= parameter)
    } else if (technique == "stdthresh" || technique == "std") {
        return((matrix(1, K, 1) %*% apply(U, 2, sd) >= parameter) & (U == matrix(1, K, 1) %*% apply(U, 2, max)))
    } else if (technique == "difference" || technique == "dtb") {
        if (is.vector(U) || nrow(U) == 1) {
            diff = rep(1, length(U))
        } else {
            Us = apply(U, 2, sort)
            diff = Us[nrow(Us),] - Us[nrow(Us) - 1,]
        }
        return(((matrix(1, K, 1) %*% diff) >= parameter) & (U == matrix(1, K, 1) %*% apply(U, 2, max)))
    } else if (technique == "top" || technique == "tb") {
        return(U >= (matrix(1, K, 1) %*% apply(U, 2, max) - parameter))
    } else {
        stop("Unknown binarisation technique")
    }
}
