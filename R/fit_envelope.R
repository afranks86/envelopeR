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


