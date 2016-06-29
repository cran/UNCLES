closestToSquareFactors <- function(n) {
    f = factors(n)
    d = abs(f - sqrt(n))
    j = which(d == min(d))
    sf = c(0, 0);

    if (f[j] < sqrt(n)) {
        sf[1] = f[j]
        sf[2] = n / sf[1]
    } else {
        sf[2] = f[j]
        sf[1] = n / sf[2]
    }

    return(sf)
}