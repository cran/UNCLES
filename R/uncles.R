uncles <- function(X, type = 'A', Ks = c(4, 8, 12, 16),
methods = list(kmeansKA, list(HC, method = "ward.D2"), SOMs),
methodsDetailed = list(), inparams = list(), normalise = 0,
samplesIDs = numeric(), flipSamples = list(), U = list(),
UType = 'PM', Xn = list(), relabel_technique = "minmin",
binarisation_technique = "DTB", binarisation_param = seq(0, 1, 0.1),
setsP = numeric(), setsN = numeric(),
dofuzzystretch = FALSE, wsets = numeric(), wmethods = numeric(),
GDM = numeric(), CoPaMforDatasetTrials = 1, CoPaMfinaltrials = 1) {
    ####################### Helper functions (sortClusters) & (calcMc) #######################
    sortClusters <- function(CoPaM, Mc, minGenesInClust = 7) {
        K = ncol(Mc)
        Cf = rep(0, K)

        for (i in seq(nrow(Mc), 1, -1)) {
            M = sort(Mc[i,], decreasing = TRUE)
            C = match(M, Mc[i,])
            Cf[M >= minGenesInClust & Cf == 0] = C[M >= minGenesInClust & Cf == 0]
            if (i > 1) {
                Mc[i - 1, Cf[Cf != 0]] = Inf
            }
        }

        Cf[Cf == 0] = setdiff(1:K, Cf)
        return(CoPaM[Cf,])
    }
    calcMc <- function(B) {
        if (!is.list(B)) {
            B = list(B)
        }

        Nb = length(B)
        K = ncol(B[[1]])
        Mc = matrix(, Nb, K)

        for (i in 1:Nb) {
            Mc[i,] = colSums(B[[i]])
        }

        return(Mc)
    }

    ####################### Some parameters' fixing #######################
    if (TRUE) {
        if (!is.list(X) || is.data.frame(X)) {
            X = list(X)
        }

        params = inparams
        params$type = type
        params$Ks = Ks
        params$methods = methods
        params$normalise = normalise
        params$samplesIDs = samplesIDs
        params$flipSamples = flipSamples
        params$relabel_technique = relabel_technique
        params$setsP = setsP
        params$setsN = setsN
        params$dofuzzystretch = dofuzzystretch
        params$wmethods = wmethods
        params$wsets = wsets
        params$binarisation_technique = binarisation_technique
        params$binarisation_param = binarisation_param
        params$CoPaMforDatasetTrials = CoPaMforDatasetTrials
        params$CoPaMfinaltrials = CoPaMfinaltrials

        if (isempty(params$samplesIDs)) {
            params$samplesIDs = list()
            for (l in 1:length(X)) {
                params$samplesIDs[[l]] = 1:ncol(X[[l]])
            }
        }

        if (!is.list(params$normalise)) {
            params$normalise = rep(list(params$normalise), length(X))
        } else if (length(params$normalise) == 1) {
            params$normalise = rep(params$normalise, length(X))
        }

        if (isempty(methodsDetailed)) {
            params$methods = rep(list(params$methods), length(X))
        } else {
            params$methods = methodsDetailed
        }

        if (isempty(params$setsP)) {
            if (toupper(params$type) == "A" || length(X) == 1) {
                params$setsP = 1:length(X)
            } else if (toupper(params$type) == "B") {
                params$setsP = 1:(ceiling(length(X) / 2))
            }
        }

        if (isempty(params$setsN)) {
            if (toupper(params$type) == "A" || length(X) == 1) {
                params$setsN = numeric()
            } else if (toupper(params$type) == "B") {
                params$setsN = seq(from = ceiling(length(X) / 2) + 1, to = length(X))
            }
        }

        X = X[c(params$setsP, params$setsN)]
        if (!isempty(Xn)) {
            Xn = Xn[c(params$setsP, params$setsN)]
            if (length(X) != length(Xn)) {
                stop("X and Xn do not match in number of elements")
            }
        }

        if (!isempty(U)) {
            U = U[c(params$setsP, params$setsN),]
        }

        params$samplesIDs = params$samplesIDs[c(params$setsP, params$setsN)]
        params$normalise = params$normalise[c(params$setsP, params$setsN)]
        params$methods = params$methods[c(params$setsP, params$setsN)]

        if (isempty(params$wsets)) {
            params$wsets = rep(1, length(X))
        } else {
            params$wsets = params$wsets[c(params$setsP, params$setsN)]
        }

        if (isempty(params$wmethods)) {
            params$wmethods = rep(1, length(params$methods[[1]]))
        }

        if (!isempty(GDM)) {
            GDM = GDM[, c(params$setsP, params$setsN)]
        }

        params$setsP = 1:length(params$setsP)
        if (length(X) > length(params$setsP)) {
            params$setsN = (length(params$setsP) + 1):(length(X))
        }

        params$Ds = matrix(0, length(params$Ks), 2)
        for (d in 1:length(params$Ks)) {
            params$Ds[d,] = closestToSquareFactors(params$Ks[d])
        }

        if (toupper(params$type) == "B") {
            if (!is.list(params$binarisation_param)) {
                params$binarisation_param = list(params$binarisation_param, params$binarisation_param)
            }
            if (length(params$binarisation_param) == 1) {
                params$binarisation_param = list(binarisation_param[[1]], binarisation_param[[1]])
            } else if (length(params$binarisation_param) != 2) {
                stop("Invalid format of binarisation_param")
            }
            if (is.list(params$binarisation_technique)) {
                params$binarisation_technique = as.vector(binarisation_technique)
            }
            if (length(params$binarisation_technique) == 1) {
                params$binarisation_technique = c(params$binarisation_technique, params$binarisation_technique)
            } else if (length(params$binarisation_technique) != 2) {
                stop("Invalid format of binarisation_technique")
            }
        }

        if (isempty(GDM)) {
            GDM = matrix(TRUE, nrow(X[[1]]), length(c(setsP, setsN)))
        } else {
            GDM = GDM[, c(setsP, setsN)]
        }

        params$L = length(X)
        params$NKs = length(params$Ks)
    }

    ####################### Normalise, flip samples, and combine replicates #######################
    if (isempty(Xn)) {
        Xn = list()
        for (l in 1:params$L) {
            Xtmp = X[[l]]

            # If the normalisation (6) is needed, do it before flipping samples
            if (any(params$normalise[[l]] == 6)) {
                Xtmp = normaliseMatrix(Xtmp, 6)
            }

            # Flip samples
            if (!isempty(params$flipSamples) && !isempty(params$flipSamples[[l]])) {
                Xtmp[, params$flipSamples[[l]] == 1] = 1 / Xtmp[, params$flipSamples[[l]] == 1]
                Xtmp[, params$flipSamples[[l]] == 2] = -Xtmp[, params$flipSamples[[l]] == 2]
            }

            # Combine replicates
            uniqueSamples = sort(unique(params$samplesIDs[[l]]))
            if (is.data.frame(X[[l]]))
                X[[l]] = data.frame(matrix(, nrow(X[[l]]), sum(uniqueSamples > 0)), row.names = row.names(X[[l]]))
            else
                X[[l]] = matrix(, nrow(X[[l]]), sum(uniqueSamples > 0))
            ss = 1
            for (s in 1:length(uniqueSamples)) {
                if (uniqueSamples[s] > 0) {
                    X[[l]][, ss] = apply(matrix(Xtmp[, params$samplesIDs[[l]] == uniqueSamples[s]]), 1, median)
                    ss = ss + 1
                }
            }

            # Normalise
            Xn[[l]] = normaliseMatrix(X[[l]], params$normalise[[l]])
        }
    }

    ####################### Clustering #######################
    if (isempty(U)) {
        UType = "PM"
        U = matrix(list(), params$L, params$NKs)
        for (l in 1:params$L) {
            for (ki in 1:params$NKs) {
                sprintf("Dataset %i, K = %i", l, params$Ks[ki])
                U[[l, ki]] = clusterDataset(Xn[[l]], params$Ks[ki], params$Ds[ki,], params$methods[[l]])
            }
        }
    } else {
        if (is.vector(U)) {
            U = as.matrix(U);
        }
        if (nrow(U) != params$L) {
            stop("U and X do not match")
        }
    }

    ####################### Calculate a CoPaM for each dataset at each K #######################
    CoPaMsFine = matrix(list(), nrow(U), ncol(U))
    for (l in 1:params$L) {
        for (ki in 1:params$NKs) {
            CoPaMsFineTmp = rep(list(), CoPaMforDatasetTrials)
            for (t in 1:CoPaMforDatasetTrials) {
                if (tolower(UType) == "pm") {
                    CoPaMsFineTmp[[t]] = generateCoPaM(U[[l, ki]], relabel_technique = params$relabel_technique, w = params$wmethods)
                } else if (tolower(UType) == "idx") {
                    stop("Unsupported in this version of the package. Email the maintainer Basel Abu-Jamous for details.")
                }
            }

            CoPaMsFine[[l, ki]] = generateCoPaM(CoPaMsFineTmp, relabel_technique = params$relabel_technique)

            if (params$dofuzzystretch) {
                CoPaMsFine[[l, ki]] = fuzzystretch(CoPaMsFine[[l, ki]])
            }
        }
        sprintf("CoPaM fine dataset %i", l)
    }

    ####################### Calculate the final CoPaM for each K #######################
    CoPaMs = matrix(list(), CoPaMfinaltrials, params$NKs)
    CoPaMsP = matrix(list(), CoPaMfinaltrials, params$NKs)
    CoPaMsN = matrix(list(), CoPaMfinaltrials, params$NKs)
    for (t in 1:CoPaMfinaltrials) {
        for (ki in 1:params$NKs) {
            if (toupper(params$type) == "A") {
                CoPaMs[[t, ki]] = generateCoPaM(CoPaMsFine[, ki], relabel_technique = params$relabel_technique, w = params$wsets, GDM = GDM)
            } else if (toupper(params$type) == "B") {
                CoPaMsP[[t, ki]] = generateCoPaM(CoPaMsFine[params$setsP, ki], relabel_technique = params$relabel_technique, w = params$wsets[params$setsP], GDM = GDM[, params$setsP])
                CoPaMsN[[t, ki]] = generateCoPaM(CoPaMsFine[params$setsN, ki], relabel_technique = params$relabel_technique, w = params$wsets[params$setsN], GDM = GDM[, params$setsN])
            } else {
                stop("Invalid UNCLES type. It has to be either A or B")
            }
        }
        sprintf("CoPaM final trial %i", t)
    }
    if (toupper(params$type) == "A") {
        params$CoPaMs = CoPaMs
    } else if (toupper(params$type) == "B") {
        params$CoPaMsP = CoPaMsP
        params$CoPaMsN = CoPaMsN
    }

    ####################### Binarisation, sorting clusters, and performing UNCLES A or B #######################
    if (toupper(params$type) == "A") {
        Bs = array(list(), c(CoPaMfinaltrials, length(params$binarisation_param), 1, params$NKs))
        Mc = matrix(list(), CoPaMfinaltrials, params$NKs)
    } else if (toupper(params$type) == "B") {
        Bs = array(list(), c(CoPaMfinaltrials, length(params$binarisation_param[[1]]), length(params$binarisation_param[[2]]), params$NKs))
        Mc = matrix(list(), CoPaMfinaltrials, params$NKs)
    }

    for (t in 1:CoPaMfinaltrials) {
        for (ki in 1:params$NKs) {
            if (toupper(params$type) == "A") {
                ##################### Perform UNCLES A
                # Pre-sorting binarisation
                for (i in 1:length(params$binarisation_param)) {
                    Bs[[t, i, 1, ki]] = t(binarise(CoPaMs[[t, ki]], params$Ks[ki], technique = params$binarisation_technique, params$binarisation_param[i]))
                }
                Mc[[t, ki]] = calcMc(Bs[t,,, ki])

                # Sorting
                CoPaMs[[t, ki]] = sortClusters(CoPaMs[[t, ki]], Mc[[t, ki]], 7)

                # Post-sorting binarisation
                for (i in 1:length(params$binarisation_param)) {
                    Bs[[t, i, 1, ki]] = t(binarise(CoPaMs[[t, ki]], params$Ks[ki], technique = params$binarisation_technique, params$binarisation_param[i]))
                }
                Mc[[t, ki]] = calcMc(Bs[t,,, ki])
            } else if (toupper(params$type) == "B") {
                ##################### Perform UNCLES B
                # Pre-sorting binarisation
                Bs_P = rep(list(), length(params$binarisation_param[[1]]))
                for (i in 1:length(params$binarisation_param[[1]])) {
                    Bs_P[[i]] = t(binarise(CoPaMsP[[t, ki]], params$Ks[ki], technique = params$binarisation_technique[1], params$binarisation_param[[1]][i]))
                }
                Mc_P = calcMc(Bs_P)

                Bs_N = rep(list(), length(params$binarisation_param[[2]]))
                for (i in 1:length(params$binarisation_param[[2]])) {
                    if (tolower(params$binarisation_technique[2]) == "dtb" && params$binarisation_param[[2]][i] == 0) {
                        params$binarisation_param[[2]][i] = .Machine$double.eps
                    }

                    Bs_N[[i]] = t(binarise(CoPaMsN[[t, ki]], params$Ks[ki], technique = params$binarisation_technique[2], params$binarisation_param[[2]][i]))
                }
                Mc_N = calcMc(Bs_N)

                #Sorting
                CoPaMsP[[t, ki]] = sortClusters(CoPaMsP[[t, ki]], Mc_P, 7)
                CoPaMsN[[t, ki]] = sortClusters(CoPaMsN[[t, ki]], Mc_P, 7)

                # Post-sorting binarisation
                Bs_P = rep(list(), length(params$binarisation_param[[1]]))
                for (i in 1:length(params$binarisation_param[[1]])) {
                    Bs_P[[i]] = t(binarise(CoPaMsP[[t, ki]], params$Ks[ki], technique = params$binarisation_technique[1], params$binarisation_param[[1]][i]))
                }
                Mc_P = calcMc(Bs_P)

                Bs_N = rep(list(), length(params$binarisation_param[[2]]))
                for (i in 1:length(params$binarisation_param[[2]])) {
                    if (tolower(params$binarisation_technique[2]) == "dtb" && params$binarisation_param[[2]][i] == 0) {
                        params$binarisation_param[[2]][i] = .Machine$double.eps
                    }

                    Bs_N[[i]] = t(binarise(CoPaMsN[[t, ki]], params$Ks[ki], technique = params$binarisation_technique[2], params$binarisation_param[[2]][i]))
                }
                Mc_N = calcMc(Bs_N)

                # UNCLES type 'B'
                for (i in 1:length(params$binarisation_param[[1]])) {
                    for (j in 1:length(params$binarisation_param[[2]])) {
                        Bs[[t, i, j, ki]] = Bs_P
                        Bs[[t, i, j, ki]][apply(Bs_N, 1, any),] = matrix(FALSE, sum(apply(Bs_N, 1, any)), params$Ks[ki])
                    }
                }

                # Fill Mc
                Mc[[t, ki]] = rep(list(), params$Ks[ki])
                for (k in 1:params$Ks[ki]) {
                    Mc[[t, ki]][[k]] = matrix(0, length(params$binarisation_param[[1]]), length(params$binarisation_param[[2]]))
                    for (i in 1:length(params$binarisation_param[[1]])) {
                        for (j in 1:length(params$binarisation_param[[2]])) {
                            Mc[[t, ki]][[k]][i, j] = sum(Bs[[t, i, j, ki]][, k])
                        }
                    }
                }
            }
        }
        sprintf("UNCLES trial %i", t)
    }

    ####################### Return the result #######################
    return(list(B = Bs, Mc = Mc, params = params, X = X, Xn = Xn, U = U, GDM = GDM))
}
