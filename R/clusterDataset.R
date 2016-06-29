clusterDataset <- function(X, K, D = 0, methods = list(kmeansKA)) {
    if (!is.list(methods)) {
        methods = list(methods)
    }

    U = list()

    for (ms in 1:length(methods)) {
        if (!is.list(methods[[ms]])) {
            methods[[ms]] = list(methods[[ms]])
        }
        # Get the method's function
        meth = methods[[ms]][[1]]

        # Get the name of the output variable if provided
        if (is.null(methods[[ms]]$outputVariable)) {
            outputVariable = NULL
        } else {
            outputVariable = methods[[ms]]$outputVariable
            methods[[ms]]$outputVariable <- NULL
        }

        # Get the parameters and run the method
        if (length(meth) > 1) {
            params = methods[[ms]][2:length(methods[[ms]])]
            res = meth(X, K, params)
        } else {
            params = NULL
            res = meth(X, K)
        }

        # Get the output variable from the result if relevant
        if (!is.null(outputVariable)) {
            res = res[[outputVariable]]
        }

        # Convert the result to a partition matrix if not already so
        if (is.vector(res) && max(res) > 1) {
            res = clustVec2partMat(res)
        }

        # Store the result
        U[[ms]] = res
    }

    return(U)
}