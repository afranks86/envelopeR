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
        
    } else if (distn == "covreg") {

        cook_fit <- optimize_envelope_covreg(Y, X, ...)
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
                              alpha=0, nu=0, nchunks = 1, L=0,
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

    L0 <- chol(Lambda0) * sqrt(prior_counts / n)
    Lambda0 <- Lambda0 * prior_counts / n

    beta_hat <- solve(t(X) %*% X + Lambda0) %*% (t(X) %*% Y + Lambda0 %*% Beta0)
    resid <- Y - X %*% beta_hat
    
    if(is.character(Vinit)) {
        if(Vinit == "OLS") {

            Ginit <- svd(beta_hat)$v[, 1:min(s, q), drop=FALSE]
            ## if(s > q) {
            ##     Ginit <- cbind(Ginit, NullC(Ginit)[, 1:(s-q), drop=FALSE])
            ## }

            res_proj <- (diag(p) - Ginit %*% t(Ginit)) %*% t(resid)
            Uinit <- svd(res_proj)$u[, 1:(s + r - q), drop=FALSE] %*%
                                 rustiefel(s+r-q, s+r-q)
            
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

    ## chunks can only be as large as s+r
    if(nchunks > s + r)
        nchunks <- s + r
    
    prior_diff <- Beta0 - beta_hat 
    ## compute negative log-likelihood

    
        
    pars <- list(Y=Y, resid=resid, n=n, p=p, s=s, r=r,
                 U1=U1, U0=U0, v1=v1, v0=v0, q=q,
                 prior_diff=prior_diff,
                 Lambda0=Lambda0, L=L, alpha=alpha, nu=nu)
        
    F <- function(Vcur) do.call(F_norm, c(list(V=Vcur), pars))
    dF <- function(Vcur) do.call(dF_norm, c(list(V=Vcur), pars))

    
    V <- Vinit

    
    max_count <- Inf
    count <- 1


    Vprev <- matrix(Inf, nrow=nrow(V), ncol=ncol(V))
    Fcur <- F(V)
    Fprev <- Fcur + abs(Fcur)
    Finit <- Fprev
    tol_f <- 1e-8
    tol_v <- 1e-8
    
    ## tol_x <- sqrt(sum((Vprev - V)^2)/n)
    ## tol_f <- (F(Vprev) - F(V)) / (abs(F(Vprev)) + 1)

    while((Fprev-Fcur) / (abs(Finit - Fcur) + 1) > tol_f & sqrt(sum((Vprev - V)^2)/n) > tol_v & count < max_count) {

        Fprev <- F(V)
        Vprev <- V
        
        ## df <- dF(V)
        ## rotS <- svd(df[, 1:s])$v
        ## if(r > 0) {
        ##     rotR <- svd(df[, (s+1):(s+r)])$v
        ##     V <- V %*% as.matrix(Matrix::bdiag(list(rotS, rotR)))
        ## } else {
        ##     V <- V %*% rotS
        ## }
        
        ## indices_mat <- suppressWarnings(matrix(1:ncol(V),
        ##                                        ncol=min(nchunks, ncol(V)),
        ##                                        byrow=TRUE))


        if(nchunks == 1) {
            print(sprintf("------ F(V) = %f --------", F(V)))
            start <- Sys.time()
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

            V <- Vfit
        }
        else {
            indices_mat <- suppressWarnings(
                matrix(sample(ncol(V)),
                       ncol=min(round(ncol(V) / nchunks), ncol(V)),
                       byrow=TRUE))
            
            ## indices_mat <- cbind(sample(s, round(ncol(V), nchunks), replace=TRUE),
            ##                      sample(r, round(ncol(V), nchunks), replace=TRUE))

            start <- Sys.time()
            for(i in 1:nrow(indices_mat)) {

                indices <- sort(indices_mat[i, ])
                
                s_indices <- indices[indices <= s]
                r_indices <- indices[indices > s]
                
                Vfixed <- V[, -indices, drop=FALSE]
                NullVfixed <- NullC(Vfixed)        

                ## Fixed part
                if(r > 0) {
                    Vr2 <- V[, setdiff((s+1):(s+r), r_indices), drop=FALSE]
                } else {
                    Vr2 <- matrix(0, nrow=nrow(V), ncol=0)
                }
                Vs2 <- V[, setdiff(1:s, s_indices), drop=FALSE] 

                YN <- Y %*% NullVfixed
                RN <- resid %*% NullVfixed
                PN <- prior_diff %*% NullVfixed
                LPN <- L0 %*% PN
                NV <- t(NullVfixed) %*% V

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
                pars_i <- pars
                pars_i$Y <- YN
                pars_i$resid <- RN
                pars_i$prior_diff <- PN
                
                Fi <- function(Vi)  {
                    if(length(s_indices) > 0)
                        NV[, s_indices] <- Vi[, 1:length(s_indices), drop=FALSE]
                    if(length(r_indices) > 0)
                        NV[, r_indices] <- Vi[, (length(s_indices)+1):length(indices), drop=FALSE]
                    do.call(F_norm, c(list(V=NV), pars_i))
                }

                pars_dfi <- c(pars_i,
                              list(AVs_part=AVs_part, AVr_part=AVr_part,
                                   LPN = LPN,
                                   rVr2=rVr2, pVr2=pVr2, VMVinv_r=VMVinv_r,
                                   rVs2=rVs2, pVs2=pVs2,
                                   VMVinv_s = VMVinv_s, indices,
                                   s_indices=s_indices, r_indices=r_indices,
                                   YN=YN, RN=RN))
                pars_dfi$Y <- Y
                pars_dfi$resid <- NULL
                                 

                ## dFi_norm function in envelope_functions.R
                dFi <- function(Vi) { do.call(dFi_norm,
                                            c(list(Vi=Vi),
                                              pars_dfi))
                }


                print(sprintf("------ F(V) = %f --------", F(V)))
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
        }

        Fcur <- F(V)        
        print(sprintf("F(V) = %f, time = %s", Fcur, Sys.time() - start))
        count <- count + 1


    }

    Yproj <- Y %*% V[, 1:s]
    eta_hat_env <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(V[, 1:s])

    
    list(V=V, intercept=intercept, beta_ols=beta_hat,
         beta_env=beta_env, eta_hat=eta_hat_env, F=F, dF=dF)
    
}


optimize_envelope_covreg <- function(Y, X,
                                     D = diag(ncol(Y)),
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

    L0 <- chol(Lambda0) * sqrt(prior_counts / n)
    Lambda0 <- Lambda0 * prior_counts / n

    beta_hat <- solve(t(X) %*% X + Lambda0) %*% (t(X) %*% Y + Lambda0 %*% Beta0)
    resid <- Y - X %*% beta_hat

    if(is.character(Vinit)) {
        if(Vinit == "OLS") {

            Ginit <- svd(beta_hat)$v[, 1:min(s, q), drop=FALSE]
            ## if(s > q) {
            ##     Ginit <- cbind(Ginit, NullC(Ginit)[, 1:(s-q), drop=FALSE])
            ## }

            if(s + r - q > 0) {
                res_proj <- (diag(p) - Ginit %*% t(Ginit)) %*% t(resid)
                Uinit <- svd(res_proj)$u[, 1:(s + r - q), drop=FALSE] %*%
                                     rustiefel(s+r-q, s+r-q)
            } else {
                Uinit  <- matrix(0, nrow=p, ncol=0)
            }
            
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

    ## chunks can only be as large as s+r
    if(nchunks > s + r)
        nchunks <- s + r
    
    prior_diff <- Beta0 - beta_hat 
    ## compute negative log-likelihood

    ## print("Compiling Stan Model...")
    ## sm  <- rstan::stan_model(file="src/stan_files/cov_regression.stan")
    sm  <- NULL

    ## Need to initialize these
    eta  <- beta_hat %*% Vinit
    z <- (Y - X %*% beta_hat) %*% Vinit
    sig_inv  <- solve(t(z) %*% z)
    
    SigInvXList  <-  lapply(1:n, function(i) sig_inv)
    etaSigXList  <-  lapply(1:n, function(i) eta %*% sig_inv)

    pars <- list(Y=Y, resid=resid, n=n, p=p, s=s, r=r,
                 SigInvXList=SigInvXList, etaSigXList = etaSigXList,
                 U1=U1, U0=U0, v1=v1, v0=v0, q=q,
                 prior_diff=prior_diff,
                 Lambda0=Lambda0, L=L, alpha=alpha, nu=nu)


    F <- function(Vcur) do.call(F_cov_reg, c(list(V=Vcur), pars))
    dF <- function(Vcur) do.call(dF_cov_reg, c(list(V=Vcur), pars))

    V  <- Vinit
    max_count <- Inf
    count <- 1

    Vprev <- matrix(Inf, nrow=nrow(V), ncol=ncol(V))
    Fcur <- F(V)
    Fprev <- Fcur + abs(Fcur)
    Finit <- Fprev
    tol_f <- 1e-8
    tol_v <- 1e-8

  ## (Fprev-Fcur) / (abs(Finit - Fcur) + 1) > tol_f &
    while(sqrt(sum((Vprev - V)^2)/n) > tol_v & count < max_count) {

              start  <- Sys.time()
              
              Fprev <- F(V)
              Vprev <- V

              ## E - Step
              YV  <- Y %*% V
              estep  <- covariance_regression_estep(YV=YV, X=X,
                                                    method="covreg",
                                                    sm=sm)
              

              SigInvXList  <-  estep$SigInvList
              etaSigXList  <-  estep$etaSigInvList
              
              pars <- list(Y=Y, resid=resid, n=n, p=p, s=s, r=r, q=q,
                           SigInvXList=SigInvXList, etaSigXList = etaSigXList,
                           U1=U1, U0=U0, v1=v1, v0=v0,
                           prior_diff=prior_diff,
                           Lambda0=Lambda0, alpha=alpha, nu=nu, L=L)

              ## M - Step
                            
              mstep  <- covariance_regression_mstep(V,
                                                    searchParams=searchParams,
                                                    maxIters=maxIters,
                                                    pars)

              V  <- mstep$V
              Fcur  <- mstep$Fcur
              print(sprintf("F(V) = %f, time = %s", Fcur, Sys.time() - start))
              count <- count + 1

    }

    Yproj <- Y %*% V[, 1:s]
    eta_hat_env <- solve(t(X) %*% X + Lambda0) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(V[, 1:s])

    
    list(V=V, intercept=intercept, beta_ols=beta_hat,
         beta_env=beta_env, eta_hat=eta_hat_env, F=F, dF=dF,
         covariance_list = estep)
    
}


optimize_envelope_cook <- function(Y, X, D = diag(ncol(Y)),
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
    
    L0 <- chol(Lambda0) * sqrt(prior_counts / n)
    Lambda0 <- Lambda0 * prior_counts / n

    beta_hat <- solve(t(X) %*% X + Lambda0) %*% (t(X) %*% Y + Lambda0 %*% Beta0)
    residual <- Y - X %*% beta_hat

    if(is.character(Vinit)) {
        if(Vinit == "OLS") {

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

    ## res <- cv.glmnet(X, Y, family="mgaussian", alpha=1)
    ## lambda_min <- res$lambda.min

    ## res <- glmnet(X, Y, family="mgaussian", alpha=1, lambda=lambda_min)
    
    prior_diff <- Beta0 - beta_hat 

    evals <- eigen(t(X) %*% X)$values
    if (evals[length(evals)] < 1e-6) {
        browser()
    }

    MU <- t(Y) %*% Y / n  ##diag(p)

    ## shrinakge estimation
    MUvecs <- eigen(MU)$vectors[, (getRank(MU) + 1):p]

    Mvecs <- eigen(M)$vectors[, (getRank(M) + 1):p]

    U0ing <- solve(U0)
    if(n < p)
        MUinv <- U0inv - U0inv %*% t(Y) %*% (diag() + Y %*% U0inv %*% t(Y))
    Vinit <- Vinit[, 1:s]

    V <- Vinit

    for(iter in 1:3) {



        G <- V[, 1:s, drop=FALSE]
        if(r > 0) {
            G0 <- V[, (s+1):(s+r)]
            YG0 <- Y_tilde %*%  G0
            G0part <- as.numeric((n + r + v0 - 1)/2 *
                                 determinant(crossprod(YG0) + U0,
                                             logarithm = TRUE)$modulus)

            if(r + s < p){
                YV <- Y_tilde %*% V
                sig2part <- (n*(p-s-r)/2 + alpha) *
                    log(sum(Y_tilde^2)/2 - sum(YV^2)/2 + nu)
            } else {
                sig2part <- 0
            }
            
        } else {
            G0part <- 0

            if(r + s < p){ 
                YV <- Y_tilde %*% V
                sig2part <- (n*(p-s-r)/2 + alpha) *
                    log(sum(Y_tilde^2) / 2 - sum(YV^2)/2 + nu)
            } else {
                sig2part <- 0
            }

            
            indices_mat <- matrix(sample(ncol(Vinit)), ncol=nchunks, byrow=TRUE)
            
            F <- function(V) {
                RV <- residual %*% V
                det1 <- (n + s + v1 + 1 - q)/2 *determinant(crossprod(RV), logarithm = TRUE)$modulus
                det2 <- (n + r + v0 - 1)/2 %*% determinant(t(V) %*% MUinv %*% V, logarithm = TRUE)$modulus
                det1 + det2
            }
            
            dF <- function(V) {
                RV <- resid %*% V
                2 * t(resid) %*% RV %*% solve(t(RV) %*% RV) +
                    2 * MUinv %*% V %*% solve(t(Vs) %*% MUinv %*% V)
            }

            print(sprintf("------ F(V) = %f --------", F(V)))


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

                print("Fitting Stiefel manifold")
                Vfit <- optStiefel(
                    Fi,
                    dFi,
                    method = "bb",
                    Vinit = t(NullV2) %*% V1,
                    verbose = FALSE,
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

}


if(FALSE) {

    library(tidyverse)
    library(mvtnorm)
    library(rstiefel)
    library(glmnet)
    source("utility_functions.R")
    source("envelope_functions.R")

    n <- 200
    p <- 800
    s <- 40
    r <- 0
    q <- 5
    
    ## n <- 1000
    ## p <- 500
    ## s <- 50
    ## r <- 10
    ## q <- 10

    X <- matrix(rnorm(n*q), nrow=n, ncol=q)
    
    V <- rustiefel(p, s)
    V <- rbind(diag(s), matrix(0, nrow=p-s, ncol=s))
    
    eta <- matrix(rnorm(q*s, sd=0.1), nrow=q)
    beta <- eta %*% t(V)

    YV <- X %*% eta + rmvnorm(n, mean= rep(0, s),
                              sigma = diag(sort(rexp(s, 1/8), decreasing=TRUE)))

    if(r > 0) {
        U <- NullC(V) %*% rustiefel(p-s, r)
        YU <- rmvnorm(n, mean= rep(0, r),
                      sigma = diag(sort(rexp(r, 1/8), decreasing=TRUE)))
    } else {
        U <- matrix(nrow=p, ncol=0)
        YU <- matrix(0, nrow=n, ncol=0)
    }
    
    Y <- YV %*% t(V) + YU %*% t(U) +
        matrix(rnorm(n * p, sd=3), nrow=n, ncol=p)


    res <- fit_envelope(Y, X, s=s, r=0,
                        distn="normal",
                        nchunks=1,
                        Vinit="COV", use_py=FALSE,
                        maxIters=10000,
                        searchParams=list(rho=1e-5, eta=0.9),
                        L = 50,
                        prior_counts = 50)
    res$F(cbind(V, U))
    res$F(rustiefel(p, s+r))
    res$F(res$V)
    tr(res$V %*% t(res$V) %*% cbind(V, U) %*% t(cbind(V, U)))  / (s+r)

    sum((res$beta_ols - beta)^2)
    sum((res$beta_env - beta)^2)
    mean(rowSums(res$V^2) < 1e-10)
    mean(rowSums(cbind(V, U)) < 1e-10)
    
    #######################################################################

    ## Covariance Regression model

    ## Generate some fake data
    library(mvtnorm)
    P <- 50
    N <- 50
    S <- 20
    R <- 20
    Q <- 2
    
    ## Generate envelope
    V <- rustiefel(P, S)
    Vperp <- NullC(V)

    ## regression parameters

    ## eta small relative to evals gives based improvements in MSE

    ## eta <- matrix(c(-1, 1), ncol=1)

    eta <- matrix(c(-10, 10, 10, 2), nrow=2, ncol=2)
    eta <- matrix(rnorm(S^2, 0, 1), nrow=S)
    Beta <- t(eta) %*% t(V)




    X <- matrix(rnorm(N*Q), nrow=N, ncol=Q)
    
    eta <- matrix(rnorm(Q*S, sd=2), nrow=Q)

    beta <- matrix(rnorm(Q*S, sd=2), nrow=Q)


    get_cov <- function(X, beta) {
        beta_X <- X %*% beta
        sig_x <- t(beta_X) %*%  beta_X
        sig_x
    }


    SigInvXList <- lapply(1:N, function(i) solve(get_cov(X[i, ], beta) + diag(S)))
    etaSigXList <- lapply(1:N, function(i) eta %*% SigInvXList[[i]])
    
    Y <- sapply(1:N, function(i)
        rmvnorm(1,
                mean = X[i, ] %*% eta %*% t(V),
                sigma = V %*% solve(SigInvXList[[i]]) %*% t(V) + diag(P))) %>% t
    svd(Y)$d

    
    ## covariates

    X <- matrix(rmvnorm(N, rep(0, S), sigma=diag(rep(5^2, S))), ncol=S)
                                        #X <- rmvnorm(N, mean=c(0,0), sigma=diag(c(5,5)))

    ## covariance matrices
    psi <- diag(sort(rexp(S, 1 / 100), decreasing = TRUE))
    ## psi <- diag(c(10, 10))
    ## psi.0 <- diag(rep(100, P-S))
    psi.0 <- diag(sort(rexp(P - S, 1 / 100), decreasing = TRUE))

    ecdf(diag(psi.0))(diag(psi))

    Sigma <- V %*% psi %*% t(V) + Vperp %*% psi.0 %*% t(Vperp)

    Y <- sapply(1:N, function(i) rmvn(1, X[i, ] %*% Beta, Sigma)) %>% t

    betaHat <- solve(t(X) %*% X) %*% t(X) %*% Y

    R <- 2
    betav <- svd(betaHat)$v
    Vinit <- cbind(betav, NullC(betav)[, 1:R])
    prior_counts <- 10



    
    X <- matrix(rnorm(n*q), nrow=n, ncol=q)
    
    V <- rustiefel(p, s)
    V <- rbind(diag(s), matrix(0, nrow=p-s, ncol=s))
    
    eta <- matrix(rnorm(q*s, sd=0.1), nrow=q)
    beta <- eta %*% t(V)

    YV <- X %*% eta + rmvnorm(n, mean= rep(0, s),
                              sigma = diag(sort(rexp(s, 1/8), decreasing=TRUE)))

    if(r > 0) {
        U <- NullC(V) %*% rustiefel(p-s, r)
        YU <- rmvnorm(n, mean= rep(0, r),
                      sigma = diag(sort(rexp(r, 1/8), decreasing=TRUE)))
    } else {
        U <- matrix(nrow=p, ncol=0)
        YU <- matrix(0, nrow=n, ncol=0)
    }
    
    Y <- YV %*% t(V) + YU %*% t(U) +
        matrix(rnorm(n * p, sd=3), nrow=n, ncol=p)


    res <- fit_envelope(Y, X, s=s, r=0,
                        distn="normal",
                        nchunks=1,
                        Vinit="COV", use_py=FALSE,
                        maxIters=10000,
                        searchParams=list(rho=1e-5, eta=0.9),
                        L = 50,
                        prior_counts = 50)
    res$F(cbind(V, U))
    res$F(rustiefel(p, s+r))
    res$F(res$V)
    tr(res$V %*% t(res$V) %*% cbind(V, U) %*% t(cbind(V, U)))  / (s+r)

    sum((res$beta_ols - beta)^2)
    sum((res$beta_env - beta)^2)
    mean(rowSums(res$V^2) < 1e-10)
    mean(rowSums(cbind(V, U)) < 1e-10)


    #################################################################
    
    
    res_kd <- fit_envelope(Y, X, s=s, r=r,
                           distn="normal",
                           nchunks=3,
                           prior_counts = 0,
                           Vinit="OLS", use_py=FALSE,
                           maxIters=1000, 
                           searchParams=list(rho=1e-5, eta=0.9))

    sum((res_kd$beta_ols - beta)^2)
    sum((res_kd$beta_env - beta)^2)

    
    res_cv <- cv.glmnet(X, Y, family="mgaussian", alpha=0, intercept=FALSE)
    res <- glmnet(X, Y, family="mgaussian", lambda=res_cv$lambda.min, alpha=0, intercept=FALSE)
    sum((do.call(cbind, coef(res))[-1, ] - beta)^2)
    
    
    res_kd$F(cbind(V, U))
    res_kd$F(rustiefel(p, s+r))
    res_kd$F(res_kd$V)
    res_kd$F(res$V)
    res$F(cbind(V, U))
    

    vbeta <- svd(beta)$v
    
    tr(res_kd$V[, 1:s] %*% t(res_kd$V[, 1:s]) %*% V %*% t(V)) / s
    tr(res_kd$V[, 1:s] %*% t(res_kd$V[, 1:s]) %*% vbeta %*% t(vbeta)) / ncol(vbeta)
    tr(res$V[, 1:s] %*% t(res$V[, 1:s]) %*% V %*% t(V)) / s
    tr(res$V[, 1:s] %*% t(res$V[, 1:s]) %*% vbeta %*% t(vbeta)) / ncol(vbeta)
    
    tr(res_kd$V[, (s+1):(s+r)] %*% t(res_kd$V[, (s+1):(s+r)]) %*%
       U %*% t(U)) / r
    tr(res_kd$V %*% t(res_kd$V) %*% cbind(V, U) %*% t(cbind(V, U))) / (s+r)

    res_kd$F(res_kd$V)
    res_kd$F(res$V)
    


    res <- fit_envelope(Y, X, s=s, r=r,
                        distn="normal",
                        Vinit="OLS", use_py=FALSE,
                        maxIters=10000,
                        searchParams=list(rho=1e-4, eta=0.9))
    

    
    Yproj <- Y %*% res$V
    eta_hat_env <- solve(t(X) %*% X) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(res_kd$V)
    sum((beta - beta_env)^2)


    
    res_kd$V

    sum((res_kd$beta_ols - beta)^2)
    sum((res_kd$beta_env - beta)^2)

    res_kd$F(res_kd$V)
    res_kd$F(res$V)

    res$F(res_kd$V)

    res <- fit_envelope(Y, X, s=s, r=r,
                        distn="cook",
                        Vinit="COV", use_py=FALSE,
                        maxIters=10000, 
                        searchParams=list(rho=1e-4, eta=0.9))

    sum((res$beta_ols - beta)^2)
    sum((res$beta_env - beta)^2)

    res_kd <- fit_envelope(Y, X, s=s, r=r,
                           distn="cook_kd",
                           Vinit="OLS", use_py=FALSE,
                           maxIters=1000, 
                           searchParams=list(rho=1e-4, eta=0.9))
    Vhat1 <- res_kd$V[, 1:s, drop=FALSE]
    res_kd$F(res_kd$V)
    res_kd$F(V)

    Yproj <- Y %*% res$V
    eta_hat_env <- solve(t(X) %*% X) %*% t(X) %*% Yproj
    beta_env <- eta_hat_env %*% t(res$V)
    sum((beta - beta_env)^2)

    
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

    n <- 200
    p <- 10
    s <- 2
    r <- 0
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

    YV <- X %*% eta + rmvnorm(n, mean= rep(0, s),
                              sigma = diag(sort(rexp(s, 1/4), decreasing=TRUE)))

    if(r > 0) {
        U <- NullC(V) %*% rustiefel(p-s, r)
        YU <- rmvnorm(n, mean= rep(0, r),
                      sigma = diag(sort(rexp(r, 1/64), decreasing=TRUE)))
    } else {
        U <- matrix(nrow=p, ncol=0)
        YU <- matrix(0, nrow=n, ncol=0)
    }
    
    Y <- YV %*% t(V) + YU %*% t(U) +
        matrix(rnorm(n * p, sd=2), nrow=n, ncol=p)




    envelope_fit <- fit_envelope(Y, X,
                                 distn = "normal",
                                 s=s, r=0,
                                 nchunks=1,
                                 Vinit=rustiefel(p, s+r), use_py=FALSE,
                                 prior_counts = 0, maxIters=100,
                                 search_params=list(rho=1e-5, eta=0.9))


    






    x <- matrix(rnorm(100), ncol=10, nrow=10)
    Sig <- t(x) %*% x
    g <- rustiefel(10, 2)
    g2 <- NullC(g)
    
    a <- det(t(g) %*% Sig %*% g)
    b <- det(t(g) %*% solve(Sig) %*% g)
    a * det(t(g2) %*% solve(Sig) %*% g2)
    
    
    log(det(Sig))
    log(det(t(g) %*% Sig %*% g)) + log(det(t(g2) %*% Sig %*% g2))


}




