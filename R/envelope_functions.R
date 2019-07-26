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
dFi_norm <- function(Vi, Y, s_indices, r_indices,
                resid, n, p, s, r, q, v1, v0, U1, U0,
                prior_diff, Lambda0, alpha, nu, L) {

    Vs1 <- Vi[, 1:length(s_indices), drop=FALSE]
    if(length(r_indices) > 0)
        Vr1 <- Vi[, (length(s_indices)+1):length(indices), drop=FALSE]
    else
        Vr1 <- matrix(nrow=nrow(Vs1), ncol=0)

    if(length(r_indices) > 0) {
        ## For speed
        if(n >= p) {
            M <- t(Y) %*% Y
            MV2 <- M %*% Vr2
            if(ncol(Vr2) > 0)
                VMVinv <- solve(t(Vr2) %*% M %*% Vr2)
            else
                VMVinv <- matrix(nrow=0, ncol=0)

            A <- t(NullVfixed) %*%
                (M - MV2 %*% VMVinv %*% t(MV2)) %*%
                NullVfixed
            
            G0part <- (n + r + v0 - 1) * A %*% Vr1 %*%
                solve(t(Vr1) %*% A %*% Vr1)
        } else {

            ## More efficient than commented part

            rVr1 <- RN %*% Vr1
            pVr1 <- LPN %*% Vr1

            W <- (crossprod(rVr1, rVr2) + crossprod(pVr1, pVr2))
            Vr1AVr1inv <- solve(crossprod(rVr1) + crossprod(rVr1) -
                                W %*% VMVinv_r %*% t(W))

            AVr1 <- crossprod(RN, rVr1) + crossprod(LPN, pVr1) -
                tcrossprod(AVr_part, W)

            G0part <- (n + r + v0 - 1) * AVr1 %*% Vr1AVr1inv

        }
        
    } else {
        G0part <- matrix(0, nrow=nrow(Vi), ncol=length(r_indices))
    }

    if(r + s < p){

        YV <- YN %*% Vi
        sig2part <- (n*(p-s-r) + 2*alpha) /
            (sum(Y^2)/2 - sum(YV^2)/2 - sum((Y %*% V[, -indices])^2)/2 + nu) * t(YN) %*% YV

        
    } else {
        sig2part <- matrix(0, nrow=nrow(Vi), ncol=ncol(Vi))
    }


    if(length(s_indices) > 0) {
        
        if(n >= p) {
            M <- crossprod(resid, resid)  +
                crossprod(prior_diff, Lambda0 %*% prior_diff) 
            MV2 <- M %*% Vs2
            if(ncol(Vs2) > 0)
                VMVinv <- solve(t(Vs2) %*% M %*% Vs2)
            else
                VMVinv <- matrix(ncol=0, nrow=0)
            
            A <- t(NullVfixed) %*%
                (M - MV2 %*% VMVinv %*% t(MV2)) %*%
                NullVfixed
            
            Gpart <- (n + s + v1 + 1 - q) * (A %*% Vs1 %*% solve(t(Vs1) %*% A %*% Vs1))
        } else {

            ## More efficient than commented part

            rVs1 <- RN %*% Vs1
            pVs1 <- LPN %*% Vs1

            W <- (crossprod(rVs1, rVs2) + crossprod(pVs1, pVs2))

            Vs1AVs1inv <- solve(crossprod(rVs1) + crossprod(pVs1) -
                                W %*% VMVinv_s %*% t(W))

            AVs1 <- crossprod(RN, rVs1) + crossprod(LPN, pVs1) -
                tcrossprod(AVs_part, W)
            
            Gpart <- (n + s + v1 + 1 - q) * AVs1 %*% Vs1AVs1inv
        }

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


F_cov_reg <- function(V, Y, resid, prior_diff, SigInvXList, etaSigXList) {

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
        YGSiginvYG <- crossprod(YG, (YG %*% SigInvXList[[i]]))
        Gpart <- Gpart + 1/2*(YGSigInvYG -  2 * tcrossprod((x[i, ] %*% etaSigXList[[i]]), YG))
    }
    
    ## Minimize the negative log likelihood
    1/n * (Gpart + G0part + sig2part) + L * sum(sqrt(rowSums(V[, 1:s]^2)))

}

dF_cov_reg <- function(V, Y_tilde, resid_tilde, prior_diff_tilde,
                       SigInvXList, etaSigXList) {

    G <- V[, 1:s, drop=FALSE]

    Gpart <- matrix(0, nrow=p, ncol=s)
    for(i in 1:n) {
        YG <- Y[i, ] %*% G
        Gpart <- crossprod(Y[i,], YG) %*% SigInvXList[[i]] -
            crossprod(Y[i, ], X[i, ] %*% etaSigXList[[i]])
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













