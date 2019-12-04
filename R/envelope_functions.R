F_norm <- function(V, Y, resid, n, p, s, r, v1, v0, U1, U0, q,
                   prior_diff, Lambda0, alpha, nu, L) {

    G <- V[, 1:s, drop=FALSE]
    if(r > 0) {
        G0 <- V[, (s+1):(s+r)]
        YG0 <- Y %*%  G0
        G0part <- as.numeric((n + r + v0 - 1)/2 *
                             determinant(crossprod(YG0) + U0,
                                         logarithm = TRUE)$modulus)

        if(r + s < p){
            YV <- Y %*% V
            sig2part <- (n*(p-s-r)/2 + alpha) *
                log(sum(Y^2)/2 - sum(YV^2)/2 + nu)
        } else {
            sig2part <- 0
        }
        
    } else {
        G0part <- 0

        if(r + s < p){ 
            YV <- Y %*% V
            sig2part <- (n*(p-s-r)/2 + alpha) *
                log(sum(Y^2) / 2 - sum(YV^2)/2 + nu)
        } else {
            sig2part <- 0
        }
        
    }


    rG <- resid %*% G
    Gpart <- as.numeric((n + s + v1 + 1 - q)/2 *
                        determinant(crossprod(rG) +
                                    t(prior_diff %*% G) %*% Lambda0 %*%
                                    (prior_diff %*% G) + U1, logarithm = TRUE)$modulus)
    
    ## Minimize the negative log likelihood
    1/n * (Gpart + G0part + sig2part) + L * sum(sqrt(rowSums(V[, 1:s]^2)))

}

## A <- t(Y - X %*% beta_hat) %*% (Y - X %*% beta_hat) +
## t(Beta0 - beta_hat) %*% Lambda0 %*% (Beta0 - beta_hat)

## compute gradient of negative log-likelihood
dF_norm <- function(V, Y, resid, n, p, s, r, q, v1, v0, U1, U0,
                    prior_diff, Lambda0, alpha, nu, L) {

    G <- V[, 1:s, drop=FALSE]

    rG <- resid %*% G

    Gpart <- (n + s + v1 + 1 - q) * (t(resid) %*% rG + t(prior_diff) %*% Lambda0 %*%
                                     (prior_diff  %*% G)) %*%
        solve(crossprod(rG) + t(prior_diff %*% G) %*% Lambda0 %*%
              (prior_diff %*% G) + U1)

    if(r > 0) {
        G0 <- V[, (s+1):(s+r)]
        YG0 <- Y %*%  G0
        G0part <- (n + r + v0 - 1) * t(Y) %*% (Y %*% G0) %*%
            solve(crossprod(YG0) + U0)
    } else {
        G0part <- matrix(0, nrow=nrow(V), ncol=0)
    }

    if(r + s < p){
        YV <- Y %*% V
        sig2part <- (n*(p-s-r) + 2*alpha) /
            (sum(Y^2)/2 - sum(YV^2)/2 + nu) * t(Y) %*% YV
    } else {
        sig2part <- matrix(0, nrow=p, ncol=s+r)
    }
    
    dV <- cbind(Gpart, matrix(0, nrow=p, ncol=r)) +
        cbind(matrix(0, nrow=p, ncol=s), G0part) +
        sig2part

    if(L > 0) {
        spart <- L * V[, 1:s]  / matrix(sqrt(rowSums(V[, 1:s]^2)),
                                        nrow=p, ncol=s, byrow=FALSE)
        Lpenalty <- cbind(spart, matrix(0, nrow=p, ncol=r))
    } else {
        Lpenalty <- 0
    }

    
    ## negative ll grad
    -1/n * dV + Lpenalty

    
}

## compute gradient of negative log-likelihood?
dFi_norm <- function(Vi, Y, YN, indices, s_indices, r_indices,
                     RN, n, p, s, r, q, v1, v0, U1, U0,
                     prior_diff, Lambda0, alpha, nu, L,
                     AVs_part, AVr_part,
                     LPN, rVr2, pVr2, VMVinv_r,
                     rVs2, pVs2, VMVinv_s) {

    
    Vs1 <- Vi[, 1:length(s_indices), drop=FALSE]
    if(length(r_indices) > 0)
        Vr1 <- Vi[, (length(s_indices)+1):length(indices), drop=FALSE]
    else
        Vr1 <- matrix(nrow=nrow(Vs1), ncol=0)

    if(length(r_indices) > 0) {

        ## More efficient than commented part

        rVr1 <- RN %*% Vr1
        pVr1 <- LPN %*% Vr1

        W <- (crossprod(rVr1, rVr2) + crossprod(pVr1, pVr2))
        Vr1AVr1inv <- solve(crossprod(rVr1) + crossprod(rVr1) -
                            W %*% VMVinv_r %*% t(W))

        AVr1 <- crossprod(RN, rVr1) + crossprod(LPN, pVr1) -
            tcrossprod(AVr_part, W)

        G0part <- (n + r + v0 - 1) * AVr1 %*% Vr1AVr1inv

        
    } else {
        G0part <- matrix(0, nrow=nrow(Vi), ncol=length(r_indices))
    }

    if(r + s < p){
        suppressMessages(browser())
        YV <- YN %*% Vi
        sig2part <- (n*(p-s-r) + 2*alpha) /
            (sum(Y^2)/2 - sum(YV^2)/2 - sum((Y %*% V[, -indices])^2)/2 + nu) * t(YN) %*% YV

        
    } else {
        sig2part <- matrix(0, nrow=nrow(Vi), ncol=ncol(Vi))
    }


    if(length(s_indices) > 0) {
        
        ## More efficient than commented part

        rVs1 <- RN %*% Vs1
        pVs1 <- LPN %*% Vs1

        W <- (crossprod(rVs1, rVs2) + crossprod(pVs1, pVs2))

        Vs1AVs1inv <- solve(crossprod(rVs1) + crossprod(pVs1) -
                            W %*% VMVinv_s %*% t(W))

        AVs1 <- crossprod(RN, rVs1) + crossprod(LPN, pVs1) -
            tcrossprod(AVs_part, W)
        
        Gpart <- (n + s + v1 + 1 - q) * AVs1 %*% Vs1AVs1inv

    } else {
        Gpart <- matrix(0, nrow=ncol(NullVfixed),
                        ncol=length(s_indices))
    }

    if( L > 0) {
        warning("L not fully supported in coord mode")
        Vstar <- V
        Vstar[, indices] <- NullVfixed %*% Vi

        tmp <- Vstar  / matrix(sqrt(rowSums(Vstar^2)),
                               nrow=p, ncol=(s+r), byrow=FALSE)

        Lpenalty <- L * (t(NullVfixed) %*% tmp)[, indices]
        
    } else {
        Lpenalty <- 0
    }

    
    dV <- 1/n * cbind(Gpart, G0part) + sig2part -
        Lpenalty
    
    
    ## negative ll grad
    -dV

}




######################### Covariance Regression ##################


#' Title
#'
#' @param V current basis for the subspace of relevant variation
#' @param Y n x p data matrix
#' @param resid n x p residual matrix
#' @param SigInvXList list of length n with entries E[Sigma(X_i)^{-1}]
#' @param etaSigXList list of length n with entries E[eta^T * Sigma(X_i)^{-1}]
#' @param prior_diff 
#'
#' 
#' @return
#' @export
#'
#' @examples
F_cov_reg <- function(V, Y, resid, n, p, s, r, q, SigInvXList, etaSigXList, 
                      U1, U0, v1, v0, prior_diff, Lambda0, alpha, nu, L) {

    G <- V[, 1:s, drop=FALSE]
    if(r > 0) {
        G0 <- V[, (s+1):(s+r)]
        YG0 <- Y %*%  G0
        G0part <- as.numeric((n + r + v0 - 1)/2 *
                             determinant(crossprod(YG0) + U0,
                                         logarithm = TRUE)$modulus)

        if(r + s < p){
            YV <- Y %*% V
            sig2part <- (n*(p-s-r)/2 + alpha) *
                log(sum(Y^2)/2 - sum(YV^2)/2 + nu)
        } else {
            sig2part <- 0
        }
        
    } else {
        G0part <- 0

        if(r + s < p){ 
            YV <- Y %*% V
            sig2part <- (n*(p-s-r)/2 + alpha) *
                log(sum(Y^2) / 2 - sum(YV^2)/2 + nu)
        } else {
            sig2part <- 0
        }
        
    }

    
    ##

    Gpart <- 0
    for(i in 1:n) {
        YG <- Y[i, ] %*% G
        YGSigInvYG <- tcrossprod(YG, (YG %*% SigInvXList[[i]]))
        Gpart <- Gpart + 1/2*(YGSigInvYG -  2 * tcrossprod((X[i, ] %*% etaSigXList[[i]]), YG))
    }
    
    ## Minimize the negative log likelihood
    1/n * (Gpart + G0part + sig2part) + L * sum(sqrt(rowSums(V[, 1:s]^2)))

}

dF_cov_reg <- function(V, Y, resid, n, p, s, r, q, SigInvXList, etaSigXList, 
                      U1, U0, v1, v0, prior_diff, Lambda0, alpha, nu, L) {


    ## For covariance regression
    G <- V[, 1:s, drop=FALSE]
    Gpart <- matrix(0, nrow=p, ncol=s)
    for(i in 1:n) {
        YG <- Y[i, ] %*% G
        Gpart <- Gpart +
            crossprod(Y[i, , drop=FALSE], YG) %*% SigInvXList[[i]] -
            crossprod(Y[i, , drop=FALSE], X[i, ] %*% etaSigXList[[i]])
    }

    if(r > 0) {
        G0 <- V[, (s+1):(s+r)]
        YG0 <- Y %*%  G0
        G0part <- (n + r + v0 - 1) * t(Y) %*% (Y %*% G0) %*%
            solve(crossprod(YG0) + U0)
    } else {
        G0part <- matrix(0, nrow=nrow(V), ncol=0)
    }

    if(r + s < p){
        YV <- Y %*% V
        sig2part <- (n*(p-s-r) + 2*alpha) /
            (sum(Y^2)/2 - sum(YV^2)/2 + nu) * t(Y) %*% YV
    } else {
        sig2part <- matrix(0, nrow=p, ncol=s+r)
    }
    
    dV <- cbind(Gpart, matrix(0, nrow=p, ncol=r)) +
        cbind(matrix(0, nrow=p, ncol=s), G0part) +
        sig2part

    if(L > 0) {
        spart <- L * V[, 1:s]  / matrix(sqrt(rowSums(V[, 1:s]^2)),
                                        nrow=p, ncol=s, byrow=FALSE)
        Lpenalty <- cbind(spart, matrix(0, nrow=p, ncol=r))
    } else {
        Lpenalty <- 0
    }

    
    ## negative ll grad
    -1/n * dV + Lpenalty
    
}

covariance_regression_estep <- function(YV, X,  sm, method="vb") {

    print("Starting e-step sampling...")
    s  <- ncol(YV)
    q  <- ncol(X)
    n  <- nrow(YV)
    
    data_list <- list(s=s, q=q, n=n, X=X, Y=YV)
    if(method == "vb") {
        stan_fit <- rstan::vb(sm, data=data_list)
    } else {
        stan_fit  <- rstan::sampling(sm, data=data_list)
    }
    print("Finished e-step sampling...")
    
    nsamples  <- dim(stan_fit)[[1]]
    
    SigInvList  <- list(nrow(YV))
    etaSigInvList  <- list(nrow(YV))
    samples <- rstan::extract(stan_fit)
    for(i in 1:nrow(X)) {
        SigInvSamples  <- lapply(1:nsamples, function(s) {
            A  <- samples$A[s, ,]
            L  <- samples$gamma[s, ,] %*%  X[i, ]
            SigInv  <- solve(tcrossprod(L) + A)
        })
        SigInvList[[i]]  <- Reduce(`+`, SigInvSamples)/nsamples
        eta  <- samples$beta
        etaSigInvSamples  <- lapply(1:nsamples, function(s) {
            t(eta[s, , ]) %*%  SigInvSamples[[s]]
        })
        
        etaSigInvList[[i]]  <- Reduce(`+`, etaSigInvSamples)/nsamples
    }
    
    list(SigInvList = SigInvList, etaSigInvList = etaSigInvList,
         samples=samples)
    
    
}


covariance_regression_mstep <- function(Vcur, searchParams, maxIters, pars) {


    F <- function(Vcur) do.call(F_cov_reg, c(list(V=Vcur), pars))
    dF <- function(Vcur) do.call(dF_cov_reg, c(list(V=Vcur), pars))

    V <- Vcur
    
    max_count <- Inf
    count <- 1

    Vprev <- matrix(Inf, nrow=nrow(V), ncol=ncol(V))
    Fcur <- F(V)
    Fprev <- Fcur + abs(Fcur)
    Finit <- Fprev
    tol_f <- 1e-8
    tol_v <- 1e-8
    
    Fprev <- F(V)
    Vprev <- V

    print(sprintf("------ F(V) = %f --------", F(V)))

    Vfit <- rstiefel::optStiefel(
                          F,
                          dF,
                          method = "bb",
                          Vinit = V,
                          verbose = TRUE,
                          maxIters = maxIters,
                          maxLineSearchIters = 25,
                          searchParams = searchParams
                      )

    list(V = Vfit, Fcur = F(Vfit))
        
    
    ## Fcur <- F(V)        
    ## count <- count + 1    

}






