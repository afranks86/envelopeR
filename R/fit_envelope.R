#' Empirical Bayes inference for the envelope model, as described in Franks et al (2018).  
#'
#' `fit_envelope` fits an envelope model to data.
#' N = sample size
#' P = number of features
#' S = subspace of mean variation
#' R = orthogonal supace of shared covariation
#' S + R <= P.  P - S - R dimensions have isotropic covariance sigma^2
#' D scales the variances
#' 
#' @param Y an n x p matrix.  Rows are independent observations of a p-variate normal.  
#' @param X an n x q matrix of q-covariates.  
#' @param D a diagonal matrix
#' @param s dimension of the subspace of relevant variation
#' @param r subspace of irrelevant non-isotropic variation?
#' @param Vinit todo
#' @return A list including: a semi-orthogonal matrix, V, the basis for the subspace of relevant variation.
#' @author Alexander Franks
#' @examples
#' 
#' @export fit_envelope
fit_envelope <- function(Y, X, distn = "normal", ...){

    if(distn == "normal") {
        normal_fit <- optimize_envelope(Y, X, ...)
        normal_fit 
        
    } else if (distn == "normal_kd") {
        cook_fit <- optimize_envelope_kd(Y, X, ...)
        cook_fit
        
    } else if (distn == "cook") {
        cook_fit <- optimize_envelope_cook(Y, X, ...)
        cook_fit
        
    } else if (distn == "cook_kd") {
        cook_fit <- optimize_envelope_cook_kd(Y, X, ...)
        cook_fit
        
    } else if (distn == "laplace") {
        

    } else if (distn == "") {

        
    } else {
        stop("distn but must be normal, laplace or t")
    }




}

#' N = sample size
#' P = number of features
#' S = subspace of mean variation
#' R = orthogonal supace of shared covariation
#' S + R <= P.  P - S - R dimensions have isotropic covariance sigma^2
#' D scales the variances
#' 
#' @author Alexander Franks
#' @examples
#' 
#' @export optimize_envelope
optimize_envelope <- function(Y, X, D = diag(ncol(Y)),
                              s=2, r=0,
                              Vinit = "OLS",
                              Lambda0 = t(X) %*% X,
                              prior_counts=0,
                              Beta0 = matrix(0, nrow=ncol(X), ncol=ncol(Y)),
                              v1=0, v0 = 0,
                              U1=matrix(0, nrow=s, ncol=s),
                              U0=matrix(0, nrow=r, ncol=r),
                              alpha=0, nu = 0, L=0,
                              center=TRUE, maxIters=1000, use_py=FALSE,
                              searchParams=NULL, ...){


    Y <- Y %*% D
    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(X)

    intercept <- rep(0, ncol(Y))
    if(center) {
        intercept <- colMeans(Y)
        Y <- scale(Y, scale=FALSE)
    } else {
        intercept <- rep(0, ncol(Y))
    }
    
    Lambda0 <- Lambda0 * prior_counts / n

    beta_hat <- solve(t(X) %*% X + Lambda0) %*% (t(X) %*% Y + Lambda0 %*% Beta0)
    resid <- Y - X %*% beta_hat

    if(is.character(Vinit)) {
        if(Vinit == "OLS") {

            Ginit <- svd(beta_hat)$v[, 1:min(s, q), drop=FALSE]
            if(s > q) {
                Ginit <- cbind(Ginit, NullC(Ginit)[, 1:(s-q), drop=FALSE])
            }

            res_proj <- (diag(p) - Ginit %*% t(Ginit)) %*% t(resid)
            
            if(r > 0) {
                Uinit <- svd(res_proj)$u[, 1:r, drop=FALSE]
            } else
                Uinit <- matrix(nrow=p, ncol=0)
            
            Vinit <- cbind(Ginit, Uinit)
            
        } else if(Vinit == "COV") {
            Vinit <- svd(Y)$v[, 1:(s+r), drop=FALSE]

        } else {
            warn("Randomly initializing")
            Vinit <- rustiefel(ncol(Y), s+r)
        }
    } else if(is.null(Vinit)) {
        Vinit <- rustiefel(ncol(Y), s+r)
    } else if( !is.matrix(Vinit) ) {
        stop("Vinit must be a semi-orthogonal matrix or string")
    } 
    
    if (s <= 0){
        stop("Require s > 0")
    }
    if (r < 0){
        stop("Require r >= 0")
    }
    if(s + r > p) {
        stop("Require s + r <= p")
    }

    if(norm(t(Vinit) %*% Vinit - diag(s+r)) > 1e-6) {
        stop("Vinit not semi-orthogonal")
    }
    
    evals <- eigen(t(X) %*% X)$values
    if (evals[length(evals)] < 1e-6) {
        warning("Perfect Colinearity in X")
    }



    
    ## shrinkage_coefs <- lw_shrinkage(resid)

    ## U0 <- diag(shrinkage_coefs$m, nrow=r)
    ## U1 <- diag(shrinkage_coefs$m, nrow=s)

    ## shrinkage_factor <- shrinkage_coefs$b2 / shrinkage_coefs$d2
    ## v1 <- shrinkage_factor / (1 - shrinkage_factor) * n
    ## v0 <- v1
    
    
    prior_diff <- Beta0 - beta_hat 
    ## compute negative log-likelihood
    F <- function(V) {

        G <- V[, 1:s, drop=FALSE]
        if(r > 0) {
            G0 <- V[, (s+1):(s+r)]
            YG0 <- Y %*%  G0
            rG0 <- resid %*%  G0 ## to remove
            G0part <- (n + r + v0 - 1)/2 * determinant((t(YG0)) %*% (YG0) + U0,
                                                   logarithm = TRUE)$modulus

            if(r + s < p){
                YV <- Y %*% V
                sig2part <- (n*(p-s-r)/2 + alpha) *
                    log(sum(Y^2)/2 - sum(YV^2)/2 + nu)
            } else {
                sig2part <- 0
            }
            
        } else {
            G0part <- 0
            YV <- Y %*% V
            sig2part <- (n*(p-s-r)/2 + alpha) *
                log(sum(Y^2) / 2 - sum(YV^2)/2 + nu)
        }

        rG <- resid %*% G
        Gpart <- (n + s + v1 + 1 - q)/2 *
            determinant(t(rG) %*% rG +  t(prior_diff %*% G) %*% Lambda0 %*%
                        (prior_diff %*% G) + U1, logarithm = TRUE)$modulus
        
        ## Minimize the negative log likelihood
        Gpart + G0part + sig2part 
        
    }

    ## A <- t(Y - X %*% beta_hat) %*% (Y - X %*% beta_hat) +
    ## t(Beta0 - beta_hat) %*% Lambda0 %*% (Beta0 - beta_hat)
    
    ## compute gradient of negative log-likelihood
    dF <- function(V) {

        G <- V[, 1:s, drop=FALSE]

        rG <- resid %*% G

        Gpart <- (n + v1 - q) * (t(resid) %*% rG + t(prior_diff) %*% Lambda0 %*%
                                 (prior_diff  %*% G)) %*%
            solve(t(rG) %*% rG + t(prior_diff %*% G) %*% Lambda0 %*%
                  (prior_diff %*% G) + U1)

        if(r > 0) {
            G0 <- V[, (s+1):(s+r)]
            YG0 <- Y %*%  G0
            G0part <- (n + v0) * t(Y) %*% (Y %*% G0) %*%
                solve(t(YG0) %*% YG0 + U0)
        } else {
            G0part <- matrix(0, nrow=nrow(V), ncol=0)
        }

        if(r + s < p){
            YV <- Y %*% V
            sig2part <- (n*(p-s-r)/2 + alpha) /
                (sum(Y^2)/2 - sum(YV^2)/2 + nu) * t(Y) %*% YV
        } else {
            sig2part <- matrix(0, nrow=p, ncol=s+r)
        }

        dV <- cbind(Gpart, matrix(0, nrow=p, ncol=r)) +
            cbind(matrix(0, nrow=p, ncol=s), G0part) +
            sig2part

        ## negative ll grad
        -dV
        
    }
    browser()
    print(sprintf("------ F(V) = %f --------", F(Vinit)))
    print("Fitting Stiefel manifold")
    if (use_py) {

        Vfit <- optim_envelope_py(Vinit,
                                  as.integer(n),
                                  as.integer(p),
                                  as.integer(s),
                                  as.integer(r),
                                  as.integer(q),
                                  U0, U1, alpha, v0, v1, S, nu, A, L,
                                  maxiters = maxIters)
    }
    else {
        ##source("test_opt.R")
        
        Vfit <- optStiefel(
            F,
            dF,
            method = "bb",
            Vinit = Vinit,
            verbose = TRUE,
            maxIters = maxIters,
            maxLineSearchIters = 20,
            searchParams = searchParams
            )
    }

    Yproj <- Y %*% Vfit[, 1:s]
    eta_hat_env <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(Vfit[, 1:s])

    list(V=Vfit, intercept=intercept, beta_ols=beta_hat,
         beta_env=beta_env, eta_hat=eta_hat_env, F=F, dF=dF)
    
}


optimize_envelope_kd <- function(Y, X, D = diag(ncol(Y)),
                                 s=2, r=0,
                                 Vinit = "OLS",
                                 Lambda0 = t(X) %*% X,
                                 prior_counts=0,
                                 Beta0 = matrix(0, nrow=ncol(X), ncol=ncol(Y)),
                                 v1=0, v0 = 0,
                                 U1=matrix(0, nrow=s, ncol=s),
                                 U0=matrix(0, nrow=r, ncol=r),
                                 alpha=0, nu=0, nchunks = 1, L=0,
                                 center=TRUE, maxIters=1000, use_py=FALSE,
                                 searchParams=NULL, ...){
    
    Y <- Y %*% D
    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(X)

    intercept <- rep(0, ncol(Y))
    if(center) {
        intercept <- colMeans(Y)
        Y <- scale(Y, scale=FALSE)
    } else {
        intercept <- rep(0, ncol(Y))
    }

    L0 <- chol(Lambda0) * sqrt(prior_counts / n)
    Lambda0 <- Lambda0 * prior_counts / n

    

    beta_hat <- solve(t(X) %*% X + Lambda0) %*% (t(X) %*% Y + Lambda0 %*% Beta0)
    resid <- Y - X %*% beta_hat

    if(is.character(Vinit)) {
        if(Vinit == "OLS") {

            Ginit <- svd(beta_hat)$v[, 1:min(s, q), drop=FALSE]
            if(s > q) {
                Ginit <- cbind(Ginit, NullC(Ginit)[, 1:(s-q), drop=FALSE])
            }

            res_proj <- (diag(p) - Ginit %*% t(Ginit)) %*% t(resid)
            
            if(r > 0) {
                Uinit <- svd(res_proj)$u[, 1:r, drop=FALSE]
            } else
                Uinit <- matrix(nrow=p, ncol=0)
            
            Vinit <- cbind(Ginit, Uinit)
            
        } else if(Vinit == "COV") {
            Vinit <- svd(Y)$v[, 1:(s+r), drop=FALSE]

        } else {
            warn("Randomly initializing")
            Vinit <- rustiefel(ncol(Y), s+r)
        }
    } else if(is.null(Vinit)) {
        Vinit <- rustiefel(ncol(Y), s+r)
    } else if( !is.matrix(Vinit) ) {
        stop("Vinit must be a semi-orthogonal matrix or string")
    } 
    
    if (s <= 0){
        stop("Require s > 0")
    }
    if (r < 0){
        stop("Require r >= 0")
    }
    if(s + r > p) {
        stop("Require s + r <= p")
    }

    if(norm(t(Vinit) %*% Vinit - diag(s+r)) > 1e-6) {
        stop("Vinit not semi-orthogonal")
    }
    
    evals <- eigen(t(X) %*% X)$values
    if (evals[length(evals)] < 1e-6) {
        warning("Perfect Colinearity in X")
    }

    ## m1 <- median(svd(resid)$d)^2
    ## m0 <- median(svd(Y)$d)^2
    
    ## shrinkage_coefs <- lw_shrinkage(resid)

    ## U0 <- diag(shrinkage_coefs$m, nrow=r)
    ## U1 <- diag(shrinkage_coefs$m, nrow=s)

    ## shrinkage_factor <- shrinkage_coefs$b2 / shrinkage_coefs$d2
    ## v1 <- shrinkage_factor / (1 - shrinkage_factor) * n
    ## v0 <- v1
    
    
    prior_diff <- Beta0 - beta_hat 
    ## compute negative log-likelihood



    ## if(n >= p) {
    ##     M <- t(Y) %*% Y
    ##     MV2 <- M %*% Vr2
    ##     if(ncol(Vr2) > 0)
    ##         VMVinv <- solve(t(Vr2) %*% M %*% Vr2)
    ##     else
    ##         VMVinv <- matrix(nrow=0, ncol=0)

    ##     A <- t(NullVfixed) %*%
    ##         (M - MV2 %*% VMVinv %*% t(MV2)) %*%
    ##         NullVfixed
        
    ## } else {

    ##     YNR1 <- Y %*% NullVfixed %*% Vr1
    ##     YR2 <- (Y %*% Vr2)
    ##     M <- t(YNR1) %*% YNR1
    ##     if(ncol(Vr2) > 0)
    ##         VMVinv <- solve(t(YR2) %*% YR2)
    ##     else
    ##         VMVinv <- matrix(nrow=0, ncol=0)

    ##     VAVinv <- solve((M - t(YNR1) %*% YR2 %*%
    ##                      VMVinv %*% t(YR2) %*% YNR1))

    ## }
    

    F_factory <- function(V, Y_tilde, resid_tilde, prior_diff_tilde) {
        G <- V[, 1:s, drop=FALSE]
        if(r > 0) {
            G0 <- V[, (s+1):(s+r)]
            YG0 <- Y_tilde %*%  G0
            rG0 <- resid_tilde %*%  G0 ## to remove
            G0part <- (n + r + v0 - 1)/2 * determinant((t(YG0)) %*% (YG0) + U0,
                                                       logarithm = TRUE)$modulus

            if(r + s < p){
                YV <- Y_tilde %*% V
                sig2part <- (n*(p-s-r)/2 + alpha) *
                    log(sum(Y_tilde^2)/2 - sum(YV^2)/2 + nu)
            } else {
                sig2part <- 0
            }
            
        } else {
            G0part <- 0
            YV <- Y_tilde %*% V
            sig2part <- (n*(p-s-r)/2 + alpha) *
                log(sum(Y_tilde^2) / 2 - sum(YV^2)/2 + nu)
        }


        rG <- resid_tilde %*% G
        Gpart <- (n + s + v1 + 1 - q)/2 *
            determinant(t(rG) %*% rG +  t(prior_diff_tilde %*% G) %*% Lambda0 %*%
                        (prior_diff_tilde %*% G) + U1, logarithm = TRUE)$modulus
        
        ## Minimize the negative log likelihood
        Gpart + G0part + sig2part
    }
    
    F <- function(V) { F_factory(V, Y, resid, prior_diff) }

    ## A <- t(Y - X %*% beta_hat) %*% (Y - X %*% beta_hat) +
    ## t(Beta0 - beta_hat) %*% Lambda0 %*% (Beta0 - beta_hat)
    
    ## compute gradient of negative log-likelihood
    dF <- function(V) {

        G <- V[, 1:s, drop=FALSE]

        rG <- resid %*% G

        Gpart <- (n + s + v1 + 1 - q) * (t(resid) %*% rG + t(prior_diff) %*% Lambda0 %*%
                                 (prior_diff  %*% G)) %*%
            solve(t(rG) %*% rG + t(prior_diff %*% G) %*% Lambda0 %*%
                  (prior_diff %*% G) + U1)

        if(r > 0) {
            G0 <- V[, (s+1):(s+r)]
            YG0 <- Y %*%  G0
            G0part <- (n + r + v0 - 1) * t(Y) %*% (Y %*% G0) %*%
                solve(t(YG0) %*% YG0 + U0)
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

        ## negative ll grad
        -dV
        
    }
    
    V <- Vinit
    ## (Until convergence do)
    
    max_count <- 10
    count <- 1

    ## if(nchunks == ncol(V)) {
    ##     stop("nchunks = ncol(V).  Don't use _kd opt.")
    ## }

    while(TRUE & count < max_count) {
        

        
        print(sprintf("------ F(V) = %f --------", F(V)))

        ## V <- V %*% as.matrix(Matrix::bdiag(list(rustiefel(s, s), rustiefel(r,r))))

        df <- dF(V)
        rotS <- svd(df[, 1:s])$v
        rotR <- svd(df[, (s+1):(s+r)])$v
        V <- V %*% as.matrix(Matrix::bdiag(list(rotS, rotR)))
        
        ## indices_mat <- suppressWarnings(matrix(1:ncol(V),
        ##                                        ncol=min(nchunks, ncol(V)),
        ##                                        byrow=TRUE))

        
        
        indices_mat <- suppressWarnings(
            matrix(sample(ncol(V)),
                   ncol=min(round(ncol(V) / nchunks), ncol(V)),
                   byrow=TRUE))

        ## indices_mat <- cbind(sample(s, round(ncol(V), nchunks), replace=TRUE),
        ##                      sample(r, round(ncol(V), nchunks), replace=TRUE))
        
        for(i in 1:nrow(indices_mat)) {

            indices <- sort(indices_mat[i, ])
            
            s_indices <- indices[indices <= s]
            r_indices <- indices[indices > s]
            
            Vfixed <- V[, -indices, drop=FALSE]
            NullVfixed <- NullC(Vfixed)        

            ## Fixed part
            Vr2 <- V[, setdiff((s+1):(s+r), r_indices), drop=FALSE] 
            Vs2 <- V[, setdiff(1:s, s_indices), drop=FALSE] 

            YN <- Y %*% NullVfixed
            RN <- resid %*% NullVfixed
            LPN <- L0 %*% prior_diff %*% NullVfixed

            rVr2 <- resid %*% Vr2
            pVr2 <- L0 %*% prior_diff %*% Vr2

            if(ncol(Vr2) > 0)                            
                VMVinv_r <- solve(t(rVr2) %*% rVr2 + t(pVr2) %*% pVr2)
            else
                VMVinv_r <- matrix(ncol=0, nrow=0)
            
            rVs2 <- resid %*% Vs2
            pVs2 <- L0 %*% prior_diff %*% Vs2

            if(ncol(Vs2) > 0)                            
                VMVinv_s <- solve(t(rVs2) %*% rVs2 + t(pVs2) %*% pVs2)
            else
                VMVinv_s <- matrix(ncol=0, nrow=0)

            AVs_part <- (t(RN) %*% (rVs2) + t(LPN) %*% pVs2) %*% VMVinv_s
            AVr_part <- (t(RN) %*% (rVr2) + t(LPN) %*% pVr2) %*% VMVinv_r
                
            ## Vi is a (p-nchunks) x (s+r) dim matrix
                
            Fi <- function(Vi) {
                Vcur <-  V
                Vcur[, s_indices] <- NullVfixed %*% Vi[, 1:length(s_indices), drop=FALSE]
                if(length(r_indices) > 0)
                    Vcur[, r_indices] <- NullVfixed %*% Vi[, (length(s_indices)+1):length(indices), drop=FALSE]
                F(Vcur)

            }

            ## TO DO : speed up matrix computations by ordering
            ## compute gradient of negative log-likelihood
            dFi <- function(Vi) {

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

                        W <- (t(rVr1) %*% (rVr2) + t(pVr1) %*% pVr2)
                        Vr1AVr1inv <- solve(t(rVr1) %*% rVr1 + t(pVr1) %*% pVr1 -
                                            W %*% VMVinv_r %*% t(W))

                        AVr1 <- t(RN) %*% rVr1 + t(LPN) %*% pVr1 -
                            (AVr_part %*% t(W))

                        G0part <- (n + r + v0 - 1) * AVr1 %*% Vr1AVr1inv

                        
                        ## YNR1 <- Y %*% NullVfixed %*% Vr1
                        ## YR2 <- (Y %*% Vr2)
                        ## M <- t(YNR1) %*% YNR1
                        ## if(ncol(Vr2) > 0)
                        ##     VMVinv <- solve(t(YR2) %*% YR2)
                        ## else
                        ##     VMVinv <- matrix(nrow=0, ncol=0)

                        ## VAVinv <- solve((M - t(YNR1) %*% YR2 %*%
                        ##                  VMVinv %*% t(YR2) %*% YNR1))

                        ## G0part2 <- (n + r + v0 - 1) *
                        ##     (t(Y %*% NullVfixed) %*% YNR1 - t(Y %*% NullVfixed)
                        ##         %*% YR2 %*% VMVinv %*% t(YR2) %*% YNR1) %*% VAVinv
                    }
                    
                } else {
                    G0part <- matrix(0, nrow=nrow(Vi), ncol=0)
                }

                if(r + s < p){
                    ## Double check
                    YV <- YN %*% Vi
                    sig2part <- (n*(p-s-r) + 2*alpha) /
                        (sum(Y^2)/2 - sum(YV^2)/2 - sum((Y %*% V[, -indices])^2)/2 + nu) * t(YN) %*% YV

                    
                } else {
                    sig2part <- matrix(0, nrow=ncol(NullVfixed), ncol=s+r)
                }


                if(length(s_indices) > 0) {
                    
                    if(n >= p) {
                        M <- t(resid) %*% resid  +
                            t(prior_diff) %*% Lambda0 %*% prior_diff 
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

                        W <- (t(rVs1) %*% (rVs2) + t(pVs1) %*% pVs2)

                        Vs1AVs1inv <- solve(t(rVs1) %*% rVs1 + t(pVs1) %*% pVs1 -
                                     W %*% VMVinv_s %*% t(W))

                        AVs1 <- t(RN) %*% rVs1 + t(LPN) %*% pVs1 -
                            (AVs_part %*% t(W))
                        
                        Gpart <- (n + s + v1 + 1 - q) * AVs1 %*% Vs1AVs1inv

                        ## M <- t(resid) %*% resid  +
                        ##     t(prior_diff) %*% Lambda0 %*% prior_diff 
                        ## MV2 <- M %*% Vs2

                        ## if(ncol(Vs2) > 0)                            
                        ##     VMVinv <- solve(t(Vs2) %*% M %*% Vs2)
                        ## else
                        ##     VMVinv <- matrix(ncol=0, nrow=0)

                        ## A <- t(NullVfixed) %*%
                        ##     (M - MV2 %*% VMVinv %*% t(MV2)) %*%
                        ##     NullVfixed
                        
                        ## Gpart_tmp <- (n + s + v1 + 1 - q) * (A %*% Vs1 %*% solve(t(Vs1) %*% A %*% Vs1))

                        
                    }

                } else {
                    Gpart <- matrix(0, nrow=ncol(NullVfixed),
                                    ncol=length(s_indices))
                }

                dV <- cbind(Gpart, G0part) + sig2part
                
                ## negative ll grad
                -dV

            }
            browser()
            Vi_fit <- rstiefel::optStiefel(
                Fi,
                dFi,
                method = "bb",
                Vinit = t(NullVfixed) %*% V[, indices],
                verbose = TRUE,
                maxIters = maxIters,
                maxLineSearchIters = 25,
                searchParams = searchParams
            )

            V[, indices] <- NullVfixed %*% Vi_fit

        }
        count <- count + 1
    }
    
    Yproj <- Y %*% V[, 1:s]
    eta_hat_env <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(V[, 1:s])

    
    list(V=V, intercept=intercept, beta_ols=beta_hat,
         beta_env=beta_env, eta_hat=eta_hat_env, F=F, dF=dF)
    
}

#' N = sample size
#' P = number of features
#' S = subspace of mean variation
#' R = orthogonal supace of shared covariation
#' S + R <= P.  P - S - R dimensions have isotropic covariance sigma^2
#' D scales the variances
#' 
#' @author Alexander Franks
#' @examples
#' 
#' @export optimize_envelope_cook
optimize_envelope_cook <- function(Y, X, D = diag(ncol(Y)),
                              s=2, r=0,
                              Vinit = "OLS",
                              Lambda0 = t(X) %*% X,
                              prior_counts=0,
                              Beta0 = matrix(0, nrow=ncol(X), ncol=ncol(Y)),
                              v1=0, v0 = 0,
                              U1=matrix(0, nrow=s, ncol=s),
                              U0=matrix(0, nrow=r, ncol=r),
                              alpha=0, nchunks = 0, L=0,
                              center=TRUE, maxIters=1000, use_py=FALSE,
                              searchParams=NULL, ...) {

    Y <- Y %*% D
    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(X)

    intercept <- rep(0, ncol(Y))
    if(center) {
        intercept <- colMeans(Y)
        Y <- scale(Y, scale=FALSE)
    } else {
        intercept <- rep(0, ncol(Y))
    }
    
    Lambda0 <- Lambda0 * prior_counts / n
    if(is.character(Vinit)) {
        if(Vinit == "OLS") {
            beta_hat <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Y
            residual <- Y - X %*% beta_hat
            Vinit <- svd(beta_hat)$v[, 1:min(s, q), drop=FALSE]
            res_proj <- (diag(p) - Vinit %*% t(Vinit)) %*% t(residual)
            if( r > 0 ) {
                Uinit <- svd(res_proj)$u[, 1:r]
            } else
                Uinit <- matrix(nrow=p, ncol=0)
            
            if(q < (s+r)) {
                Vinit <- cbind(Vinit, Uinit)
                if(s > q) {
                    Vinit <- cbind(Vinit, NullC(Vinit)[, 1:(s-q), drop=FALSE])
                }
            }
        } else if(Vinit == "COV") {
            Vinit <- svd(Y)$v[, 1:min((s+r), q), drop=FALSE]
            if(q < (s+r)) {
                Vinit <- cbind(Vinit, NullC(Vinit)[, 1:(s+r-q), drop=FALSE])
            } 
        } else {
            warn("Randomly initializing")
            Vinit <- rustiefel(ncol(Y), s+r)
        }
    } else if(is.null(Vinit)) {
        Vinit <- rustiefel(ncol(Y), s+r)
    } else if( !is.matrix(Vinit) ) {
        stop("Vinit must be a semi-orthogonal matrix or string")
    } 
    
    if (s <= 0){
        stop("Require s > 0")
    }
    if (r < 0){
        stop("Require r >= 0")
    }
    if(s + r > p) {
        stop("Require s + r <= p")
    }
    
    if(norm(t(Vinit) %*% Vinit - diag(s+r)) > 1e-6) {
        stop("Vinit not semi-orthogonal")
    }
    
    evals <- eigen(t(X) %*% X)$values
    if (evals[length(evals)] < 1e-6) {
        warning("Perfect Colinearity in X")
    }
    
    beta_hat <- solve(t(X) %*% X + Lambda0) %*% (t(X) %*% Y + Lambda0 %*% Beta0)
    
    resid <- Y - X %*% beta_hat

    prior_diff <- Beta0 - beta_hat 

    evals <- eigen(t(X) %*% X)$values
    if (evals[length(evals)] < 1e-6) {
        browser()
    }

    betaHat <- t(Y) %*% X %*% solve(t(X) %*% X)

    res <- (Y - X %*% t(betaHat))

    MU <- t(Y) %*% Y / n  ##diag(p)

    ## shrinakge estimation
    MUvecs <- eigen(MU)$vectors[, (getRank(MU) + 1):p]

    M <- t(resid) %*% resid / n
    SigmaHat <- M

    Mvecs <- eigen(M)$vectors[, (getRank(M) + 1):p]

    ## M <- M + 0.01*diag(nrow(M))
    ## MUinv <- M + 0.01*diag(nrow(MUinv))
    print("Inverting...")
    MUinv <- solve(MU)
    Vinit <- Vinit[, 1:s, drop=FALSE]
    F <- function(V) {
        Vs <- V[, 1:s, drop=FALSE]
        RV <- resid %*% Vs
        determinant(t(RV) %*% RV / n, logarithm = TRUE)$modulus +
        determinant(t(Vs) %*% MUinv %*% Vs, logarithm = TRUE)$modulus
    }

    dF <- function(V) {
        Vs <- V[, 1:s, drop=FALSE]
        RV <- resid %*% Vs
        2 * t(resid) %*% RV %*% solve(t(RV) %*% RV) +
            2 * MUinv %*% Vs %*% solve(t(Vs) %*% MUinv %*% Vs)
    }

    print("Fitting Stiefel manifold")
    Vfit <- optStiefel(
        F,
        dF,
        method = "bb",
        Vinit = Vinit,
        verbose = TRUE,
        maxIters = maxIters,
        maxLineSearchIters = 20
    )
    
    Yproj <- Y %*% Vfit[, 1:s]
    eta_hat_env <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(Vfit[, 1:s])

    list(V=Vfit, intercept=intercept, beta_ols=beta_hat,
         beta_env=beta_env, eta_hat=eta_hat_env, F=F, dF=dF)
    

}

optimize_envelope_cook_kd <- function(Y, X, D = diag(ncol(Y)),
                              s=2, r=0,
                              Vinit = "OLS",
                              Lambda0 = t(X) %*% X,
                              prior_counts=0,
                              Beta0 = matrix(0, nrow=ncol(X), ncol=ncol(Y)),
                              v1=0, v0 = 0,
                              U1=matrix(0, nrow=s, ncol=s),
                              U0=matrix(0, nrow=r, ncol=r),
                              alpha=0, nchunks = 4, L=0,
                              center=TRUE, maxIters=1000, use_py=FALSE,
                              searchParams=NULL, ...) {

    Y <- Y %*% D
    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(X)

    intercept <- rep(0, ncol(Y))
    if(center) {
        intercept <- colMeans(Y)
        Y <- scale(Y, scale=FALSE)
    } else {
        intercept <- rep(0, ncol(Y))
    }
    
    Lambda0 <- Lambda0 * prior_counts / n

    if(is.character(Vinit)) {
        if(Vinit == "OLS") {
            beta_hat <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Y
            residual <- Y - X %*% beta_hat
            Vinit <- svd(beta_hat)$v[, 1:min(s, q), drop=FALSE]
            res_proj <- (diag(p) - Vinit %*% t(Vinit)) %*% t(residual)
            if( r > 0 ) {
                Uinit <- svd(res_proj)$u[, 1:r]
            } else
                Uinit <- matrix(nrow=p, ncol=0)
            
            if(q < (s+r)) {
                Vinit <- cbind(Vinit, Uinit)
                if(s > q) {
                    Vinit <- cbind(Vinit, NullC(Vinit)[, 1:(s-q), drop=FALSE])
                }
            }
        } else if(Vinit == "COV") {
            Vinit <- svd(Y)$v[, 1:min((s+r), q), drop=FALSE]
            if(q < (s+r)) {
                Vinit <- cbind(Vinit, NullC(Vinit)[, 1:(s+r-q), drop=FALSE])
            } 
        } else {
            warn("Randomly initializing")
            Vinit <- rustiefel(ncol(Y), s+r)
        }
    } else if(is.null(Vinit)) {
        Vinit <- rustiefel(ncol(Y), s+r)
    } else if( !is.matrix(Vinit) ) {
        stop("Vinit must be a semi-orthogonal matrix or string")
    } 
    
    if (s <= 0){
        stop("Require s > 0")
    }
    if (r < 0){
        stop("Require r >= 0")
    }
    if(s + r > p) {
        stop("Require s + r <= p")
    }

    if(norm(t(Vinit) %*% Vinit - diag(s+r)) > 1e-6) {
        stop("Vinit not semi-orthogonal")
    }
    
    evals <- eigen(t(X) %*% X)$values
    if (evals[length(evals)] < 1e-6) {
        warning("Perfect Colinearity in X")
    }
    
    beta_hat <- solve(t(X) %*% X + Lambda0) %*% (t(X) %*% Y + Lambda0 %*% Beta0)
    
    resid <- Y - X %*% beta_hat

    prior_diff <- Beta0 - beta_hat 

    evals <- eigen(t(X) %*% X)$values
    if (evals[length(evals)] < 1e-6) {
        browser()
    }

    betaHat <- t(Y) %*% X %*% solve(t(X) %*% X)

    res <- (Y - X %*% t(betaHat))

    MU <- t(Y) %*% Y / n  ##diag(p)

    ## shrinakge estimation
    MUvecs <- eigen(MU)$vectors[, (getRank(MU) + 1):p]

    M <- t(resid) %*% resid / n
    SigmaHat <- M

    Mvecs <- eigen(M)$vectors[, (getRank(M) + 1):p]

    ## M <- M + 0.01*diag(nrow(M))
    ## MUinv <- M + 0.01*diag(nrow(MUinv))
    print("Inverting...")
    MUinv <- solve(MU)
    Vinit <- Vinit[, 1:s]

    V <- Vinit

    for(iter in 1:3) {

        indices_mat <- matrix(sample(ncol(Vinit)), ncol=nchunks, byrow=TRUE)
        F <- function(V) {
            RV <- resid %*% V
            determinant(t(RV) %*% RV / n, logarithm = TRUE)$modulus +
                                                              determinant(t(V) %*% MUinv %*% V, logarithm = TRUE)$modulus
        }
        
        dF <- function(V) {
            RV <- resid %*% V
            2 * t(resid) %*% RV %*% solve(t(RV) %*% RV) +
                2 * MUinv %*% V %*% solve(t(Vs) %*% MUinv %*% V)
        }

        print(sprintf("------ F(V) = %f --------", F(V)))


    browser()
    for(i in 1:nrow(indices_mat)) {

        indices <- indices_mat[i, ]
        print(sprintf("COUNT: %i, Indices: %s", i, paste0(indices, collapse=", ")))

        V1 <- V[, indices]
        V2 <- V[, -indices]
        NullV2 <- NullC(V2)        

        ## Part 1
        MV2 <- M %*% V2
        A <-  t(NullV2) %*% (M - MV2 %*% solve(t(V2) %*% MV2) %*% t(MV2)) %*% NullV2

        ## Part 2
        MUinvV2 <- MUinv %*% V2
        MV2 <- M %*% V2
        B <-  t(NullV2) %*% (MUinv - MUinvV2 %*% solve(t(V2) %*% MUinvV2) %*% t(MUinvV2)) %*% NullV2
        
        ## log det: t(V1) %*% (M - MV2 %*% VMVing %*% t(MV2)) %*% V1
        ## V1 = NulLV2 %*% V1star

        Fi <- function(Vi) {
            determinant(t(Vi) %*% A %*% Vi,
                        logarithm = TRUE)$modulus +
            determinant(t(Vi) %*% B %*% Vi, logarithm = TRUE)$modulus
        }

        dFi <- function(Vi) {
            2 * A %*% Vi %*% solve(t(Vi) %*% A %*% Vi) +
                2 * B %*% Vi %*% solve(t(Vi) %*% B %*% Vi)
        }

        browser()
        print("Fitting Stiefel manifold")
        Vfit <- optStiefel(
            Fi,
            dFi,
            method = "bb",
            Vinit = t(NullV2) %*% V1,
            verbose = TRUE,
            maxIters = maxIters,
            maxLineSearchIters = 20
        )

        V[, indices] <- NullV2 %*% Vfit

    }



    }
    print(sprintf("------ F(V) = %f --------", F(V)))
    Yproj <- Y %*% V[, 1:s]
    eta_hat_env <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(V[, 1:s])

    list(V=V, intercept=intercept, beta_ols=beta_hat,
         beta_env=beta_env, eta_hat=eta_hat_env, F=F, dF=dF)
    

}




#' Ledoit-Wolfe Shrinkage Estimatiaon of Sigma
#'
#' 
lw_shrinkage <- function(Y) {

    ## matrix dimensions
    p <- ncol(Y)
    n <- nrow(Y)

    ## sample covariance matrix
    
    m <- sum(Y^2)/ (n*p)

    trace_S2 <- sum(svd(Y)$d^4)/n^2
    
    d2 <- (trace_S2 + m^2 *p - 2*m*sum(svd(Y)$d^2)/n) / p

    dispersion_vec <- sapply(1:nrow(Y), function(i) {

        disp <- sum(svd(Y[i, ])$d^4) + 
            trace_S2 -
            2/n * tr((Y[i, , drop=FALSE] %*% t(Y)) %*%
                     (Y %*% t(Y[i, , drop=FALSE])))
        
        disp / p

    })

    
    b2 <- min(sum(dispersion_vec) / n^2, d2)

    list(b2=b2, d2=d2, m=m)
    
}




if(FALSE) {

    library(tidyverse)
    library(mvtnorm)
    library(rstiefel)
    source("utility_functions.R")


    n <- 100
    p <- 1000
    s <- 10
    r <- 5
    q <- 5
    
    ## n <- 1000
    ## p <- 500
    ## s <- 50
    ## r <- 10
    ## q <- 10

    X <- matrix(rnorm(n*q), nrow=n, ncol=q)
    
    V <- rustiefel(p, s)
    eta <- matrix(rnorm(q*s, sd=2), nrow=q)
    beta <- eta %*% t(V)

    U <- NullC(V) %*% rustiefel(p-s, r)
    YV <- X %*% eta + rmvnorm(n, mean= rep(0, s),
                              sigma = diag(sort(rexp(s, 1/4), decreasing=TRUE)))

    YU <- rmvnorm(n, mean= rep(0, r),
                  sigma = diag(sort(rexp(r, 1/64), decreasing=TRUE)))
    
    Y <- YV %*% t(V) + YU %*% t(U) +
        matrix(rnorm(n * p, sd=0.1), nrow=n, ncol=p)


    ## ## Variance fit tests
    ## res <- fit_envelope(Y, X, s=s, r=r, Vinit=rustiefel(p, s+r),
    ##                     use_py=FALSE,
    ##                     maxIters=1000,
    ##                     U1=0.1*diag(s),
    ##                     U0=0.1*diag(r),
    ##                     searchParams=list(rho=0.1, eta=0.9))
    
    ## res <- fit_envelope(Y, X, s=s, r=r, Vinit=rustiefel(p, s+r),
    ##                     use_py=FALSE, maxIters=10000,
    ##                     U1=0.1*diag(s), U0=0.1*diag(r),
    ##                     searchParams=list(rho=0.1, eta=0.2))

    ## res <- fit_envelope(Y, X, s=2*s, r=0,
    ##                     use_py=FALSE, maxIters=3000,
    ##                     U1=0.1*diag(2*s), U0=0.1*diag(r),
    ##                     searchParams=list(rho=0.1, eta=0.2))

    ## res <- fit_envelope(Y, X, s=s, r=r,
    ##                     Vinit=cbind(res$V, NullC(res$V)[, 1:r]), use_py=FALSE,
    ##                     maxIters=20000, U1=1*diag(s), U0=1*diag(r),
    ##                     searchParams=list(rho=1e-4, eta=0.9))

    res <- fit_envelope(Y, X, s=s, r=r,
                        distn="normal",
                        Vinit="OLS", use_py=FALSE,
                        maxIters=10000,
                        searchParams=list(rho=1e-4, eta=0.9))
    res$F(cbind(V, U))
    res$F(rustiefel(p, s+r))
    res$F(res$V)

    tr(res$V[, 1:s] %*% t(res$V[, 1:s]) %*% V %*% t(V)) / s
    tr(res$V[, (s+1):(s+r)] %*% t(res$V[, (s+1):(s+r)]) %*%
       U %*% t(U)) / r
    tr(res$V %*% t(res$V) %*% cbind(V, U) %*% t(cbind(V, U))) / (s+r)

    res$V

    sum((res$beta_ols - beta)^2)
    sum((res$beta_env - beta)^2)


    
    res_kd <- fit_envelope(Y, X, s=s, r=r,
                        distn="normal_kd",
                        nchunks=1,
                        Vinit="OLS", use_py=FALSE,
                        maxIters=10, 
                        searchParams=list(rho=1e-4, eta=0.9))
    res_kd$F(cbind(V, U))
    res_kd$F(rustiefel(p, s+r))
    res_kd$F(res_kd$V)

    tr(res_kd$V[, 1:s] %*% t(res_kd$V[, 1:s]) %*% V %*% t(V)) / s
    tr(res_kd$V[, (s+1):(s+r)] %*% t(res_kd$V[, (s+1):(s+r)]) %*%
       U %*% t(U)) / r
    tr(res_kd$V %*% t(res_kd$V) %*% cbind(V, U) %*% t(cbind(V, U))) / (s+r)

    res_kd$V

    sum((res_kd$beta_ols - beta)^2)
    sum((res_kd$beta_env - beta)^2)

    res$F(res$V)
    res$F(res$V %*% ))

    res <- fit_envelope(Y, X, s=s, r=r,
                        distn="cook",
                        Vinit="OLS", use_py=FALSE,
                        maxIters=10000, 
                        searchParams=list(rho=1e-4, eta=0.9))

    res_kd <- fit_envelope(Y, X, s=s, r=r,
                        distn="cook_kd",
                        Vinit="OLS", use_py=FALSE,
                        maxIters=1000, 
                        searchParams=list(rho=1e-4, eta=0.9))
    Vhat1 <- res_kd$V[, 1:s, drop=FALSE]
    res_kd$F(res_kd$V)
    res_kd$F(V)

    
    tr(t(V) %*% Vhat1 %*% (t(Vhat1) %*% V)) / s
    tr(t(U) %*% Vhat1 %*% (t(Vhat1) %*% U))     
    tr(t(V) %*% Vhat2 %*% (t(Vhat2) %*% V)) / 2
    tr(t(V) %*% Vhat %*% (t(Vhat) %*% V)) / 4
    
    Vinit <- cbind(res_kd$V, NullC(res_kd$V)[, 1:r])
    res_kd <- fit_envelope(Y, X, s=50, r=10,
                        distn="cook_kd",
                        Vinit=Vinit, use_py=FALSE,
                        maxIters=1000, 
                        searchParams=list(rho=1e-4, eta=0.9))

    res_kd$F(res_kd$V)
    res_kd$F(V)
    



    
    res <- fit_envelope(Y, X, s=1, r=1,
                        Vinit="OLS", use_py=FALSE,
                        maxIters=1000, 
                        searchParams=list(rho=1e-4, eta=0.9))
    

    
    Vhat1 <- res$V[, 1:s, drop=FALSE]
    NV <- NullC(Vhat1)

    res2 <- fit_envelope(Y %*% (diag(p) - Vhat1 %*% t(Vhat1)), X, s=1, r=19,
                         distn="cook",
                         Vinit="OLS", use_py=FALSE,
                         maxIters=1000, 
                         searchParams=list(rho=1e-4, eta=0.9))
    Vhat2 <- res2$V[, 1:2]

    t(Vhat1) %*% Vhat2
    Vhat <- cbind(Vhat1, Vhat2)
    tr(t(V) %*% Vhat1 %*% (t(Vhat1) %*% V)) / s
    tr(t(U) %*% Vhat1 %*% (t(Vhat1) %*% U))     
    tr(t(V) %*% Vhat2 %*% (t(Vhat2) %*% V)) / 2
    tr(t(V) %*% Vhat %*% (t(Vhat) %*% V)) / 4


    beta_hat <- solve(t(X) %*% X) %*% (t(X) %*% Y)
    resid <- Y - X %*% beta_hat
    res_proj <- (diag(p) - Vhat %*% t(Vhat)) %*% t(resid)
    Uhat <- svd(res_proj)$u[, 1:r]
        
    res3 <- fit_envelope(Y, X, s=s, r=r,
                        Vinit=cbind(Vhat, Uhat), use_py=FALSE,
                        maxIters=1000, 
                        searchParams=list(rho=1e-4, eta=0.9))
    res3$F(res3$V)
    res3$F(cbind(V, U))



    
    res2 <- fit_envelope(Y, X, s=s, r=r,
                         Vinit=cbind(V, U), use_py=FALSE,
                         maxIters=1000, U1=0.1*diag(s), U0=0.1*diag(r),
                         searchParams=list(rho=0.1, eta=0.9))

    names(res)
    Vhat <- res$V[, 1:s]
    Uhat <- res$V[, (s+1):(s+r)]

    res$F(res$V)
    res$F(cbind(V, U))
    res$F(cbind(V))

    tr((t(Vhat) %*% V) %*% (t(V) %*% Vhat)) / s
    tr((t(Uhat) %*% U) %*% (t(U) %*% Uhat)) / r

    #V contained in res$V but not first s columns thereof
    
    tr((t(V) %*% res$V) %*%  (t(res$V) %*%  V) ) / s
    ## but not first s columns thereof
    tr((t(Vhat) %*% V) %*% (t(V) %*% Vhat)) / s
    
    norm(res$dF(res$V))
    norm(res$dF(cbind(V, U)))

    resid <- Y - X %*% res$beta_ols
    
    tr((t(U) %*% res$V) %*%  (t(res$V) %*%  U) ) / r
    tr((t(cbind(V, U)) %*% res$V) %*%  (t(res$V) %*%  cbind(V, U)) ) / (s+r)


    sum((beta - res$beta_ols)^2)
    sum((beta - res$beta_env)^2)

    Z <- Y %*% res$V
    beta_new <- solve(t(X) %*% X) %*% t(X) %*% Z %*% t(res$V)
    sum((beta - beta_new)^2)

    Z <- Y %*% Vhat
    beta_new <- solve(t(X) %*% X) %*% t(X) %*% Z %*% t(Vhat)
    
    sum((beta - beta_new)^2)

    Z <- Y %*% V
    beta_new <- solve(t(X) %*% X) %*% t(X) %*% Z %*% t(V)
    sum((beta - beta_new)^2)


    
    p <- 4000
    n <- 20
    Sigma = diag(sqrt(p:1))
    
    library(mvtnorm)
    Y2 <- rmvnorm(n, sigma = Sigma)

    shrinkage_coefs <- lw_shrinkage(Y2)
    shrinkage_coefs$b2 / shrinkage_coefs$d2



    A1 <- S - Sigma
    tr(A1 %*% t(A1))/p

    A2 <- Sstar - Sigma
    tr(A2 %*% t(A2))/p


    
}
