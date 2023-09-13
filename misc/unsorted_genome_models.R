# for now I will just try to create a model that will fit 1n and 2n peaks indipendently

nls_2peak <- function(x, y, estKmercov, max_iterations){

    model2 <- NULL

  cat("trying nls_2peak standard algorithm\n")
    try(model2 <- nls(y ~ N * (dnbinom(x, size = kmercov_1n / bias1, mu = kmercov_1n) +
                               dnbinom(x, size = kmercov_2n / bias2, mu = kmercov_2n)),
                    start = list(kmercov_1n=estKmercov, kmercov_2n=(2 * estKmercov), bias1 = 1, bias2 = 2, N = sum(y)),
                    control = list(minFactor=1e-12, maxiter=max_iterations)))

    if(is.null(model2)){
        cat("retrying nls_2peak with port algorithm\n")
        try(model2 <- nls(y ~ N * (dnbinom(x, size = kmercov_1n / bias1, mu = kmercov_1n) +
                                         dnbinom(x, size = kmercov_2n / bias2, mu = kmercov_2n)),
                                start = list(kmercov_1n=estKmercov, kmercov_2n=(2 * estKmercov), bias1 = 1, bias2 = 2, N = sum(y)),
                          algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)))
    }

    return(model2)
}


# for now I will just try to create a model that will fit 1n and 2n peaks indipendently
nls_2peak_conditional <- function(x, y, k, estKmercov, estLength, rEst, biasEst, max_iterations){
    model2 = NULL

  cat("trying nls_2peak_conditional standard algorithm\n")

    try(model2 <- nls(y ~ ((2*(1-(1-rEst)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                              ((1-rEst)^k)        * dnbinom(x, size = kmercov*2 / bias, mu = kmercov2)) * estLength,
                      start = list(kmercov=estKmercov, kmercov2 = (estKmercov * 2), bias = biasEst),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}


# for now I will just try to create a model that will fit 1n and 2n peaks indipendently

nls_2peak_kmer_explicit <- function(x, y, k, estKmercov, estLength, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_kmer_explicit standard algorithm\n")

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-r)^k)        * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2)) * length,
                      start = list(r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}

