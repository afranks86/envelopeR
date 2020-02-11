library(ggforce)
library(ggrepel)

## Optimal threshold from Gavish, Donoho 2014
getRank <- function(Y) {

  svals <- svd(Y)$d

  m <- max(nrow(Y), ncol(Y))
  n <- min(nrow(Y), ncol(Y))
  
  if(m==n) {
    rank <- sum(svals > 2.858*median(svals))
  } else {
    beta <- n/m
    omeg <- 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
    rank <- sum(svals > omeg*median(svals))
  }

  rank
}


## samples: an n x s x s x nsamples array
#' Title
#'
#' @param V 
#' @param samples 
#' @param n1 
#' @param n2 
#' @param view 
#' @param nlabeled 
#' @param to_plot 
#' @param obs_names 
#' @param R 
#' @param label_size 
#' @param plot_type 
#' @param col_values 
#' @param ... 
#'
#' @return
#' @export create_plots
#'
#' @examples
create_plots <- function(V, samples, n1, n2=NULL, view=c(1,2), nlabeled=20,
                         to_plot=NULL, obs_names=NULL, R=ncol(V),
                         labels=rownames(V),
                         label_size=2, plot_type="both", col_values=NULL, ...) {

  rownames(V)  <- labels
  rotation_result <- rotate_basis(V, samples, n1, n2)

  if(is.null(to_plot)) {
    nobs  <- 2
    to_plot <- c(n1, n2)
  }
  else
    nobs  <- length(to_plot)

  Vstar <- rotation_result$rotV[, view]
  O <- rotation_result$rotMat[, view]

  nsamps  <- dim(samples)[4]

  Osamps_proj <- array(dim = c(2, 2, nobs, nsamps))
  omegaSamps_proj <- array(dim = c(2, nobs, nsamps))
  cov_proj <- array(dim = c(nobs, 2, 2, nsamps))

  for(i in 1:nsamps) {
    for(k in 1:length(to_plot)) {

      cov_proj_ik <- t(O) %*% samples[to_plot[k], , , i]  %*% O
      cov_proj[k, , , i]  <- cov_proj_ik
      eig <- eigen(cov_proj_ik)
      Osamps_proj[, , k, i] <- eig$vectors
      lambda <- eig$values
      omegaSamps_proj[, k, i] <- lambda/(lambda+1)
    }
  }

  obs_to_plot  <- 1:length(to_plot)
  names(obs_to_plot)  <-  nms[to_plot]

  posterior_plot <- posteriorPlot(cov_proj,
                                  Osamps_proj, omegaSamps_proj,
                                  nsamps=nsamps,
                                  obs_to_plot=obs_to_plot,
                                  probRegion=0.95, legend=TRUE, ...)



  biplot <- covarianceBiplot(Vstar, cov_proj, obs_to_plot=obs_to_plot,
                             nlabeled=40, label_size=label_size,
                             col_values=col_values)


  if(plot_type == "both") {
    if(is.null(dev.list()))
      dev.new(width=14, height=7)

    posterior_plot + biplot
  }
  else if(plot_type == "posterior") {
    if(is.null(dev.list()))
      dev.new()
    posterior_plot
  } else if(plot_type == "biplot" ) {
    if(is.null(dev.list()))
      dev.new()
    biplot
  } else {
    stop("Invalid plot_type")
  }

}

#' Title
#'
#' @param V 
#' @param samples 
#' @param n1 
#' @param n2 
#'
#' @return
#' @export rotate_basis
#'
#' @examples
rotate_basis <- function(V, samples, n1=1, n2=NULL) {

  S <- ncol(V)

  pmPsi1 <- apply(samples[n1, , , ], 1:2, mean)

  if(is.null(n2)) {

    ## Compare subspace that explains largest source of variability
    O <- svd(pmPsi1[1:S, 1:S])$u

  } else {
    ## Compare largest difference between 2 obs
    pmPsi2 <- apply(samples[n2, , , ], 1:2, mean)

    O <- svd(pmPsi1[1:S, 1:S] - pmPsi2[1:S, 1:S])$u

  }

  Vstar <- V[, 1:S] %*% O
  list(rotV=Vstar, rotMat=O)

}

#' Title
#'
#' @param covSamps 
#' @param Osamps 
#' @param OmegaSamps 
#' @param nsamps 
#' @param obs_to_plot 
#' @param probRegion 
#' @param hline 
#' @param ymax 
#' @param type 
#' @param plotPoints 
#' @param polar 
#' @param legend 
#' @param col_values 
#'
#' @return
#' @export posteriorPlot
#'
#' @examples
posteriorPlot <- function(covSamps, Osamps, OmegaSamps, nsamps, obs_to_plot,
                          probRegion=0.95, hline=NULL,  ymax=NULL, type = "mag",
                          plotPoints=TRUE, polar=FALSE, legend=TRUE,
                          col_values=NULL, alpha=1) {

  ngroups <- length(obs_to_plot)

  group_names <- names(obs_to_plot)
  if(is.null(group_names))
    group_names = factor(1:ngroups)

  if(type=="mag") {
    ylab <- expression(lambda[1])
  }
  else if(type=="logmag") {
    ylab <- expression("log"[2]~"(" ~ lambda[1] ~ ")")
  } else if (type =="logratio") {
    ylab <- expression("log"[2]~"(" ~ lambda[1]/lambda[2] ~ ")")
  } else {
    ylab <- expression("(" ~ lambda[1]/lambda[2] ~ ")")
  }

  if(is.null(ymax))
    ymax <- 1.1*max(OmegaSamps/(1-OmegaSamps))

  ##plot(0, 0, xlim=c(-pi/2, pi/2),
  ##     ylim=c(0, ymax), cex=0, xlab=, ylab=ylab, xaxt="n", cex.axis=cex.axis, cex.lab=1.5)
  ## axis(1, at=seq(-pi/2, pi/2, by=pi/4), labels=expression(-pi/2, -pi/4, 0, pi/4, pi/2), cex.axis=cex.axis, cex.lab=1.5)

  group_pts_angle <- group_pts_eval <- group_type <- c()
  count <- 1
  for(g in obs_to_plot) {

    pmPsi <- apply(covSamps[g, , , ], 1:2, mean)

    eigPsi <- eigen(pmPsi)
    pmValues <- eigPsi$values
    pmVectors <- eigPsi$vectors
    maxIndex <- which.max(pmValues)

    hp <- getHullPoints(nsamps, pmPsi, OmegaSamps[, g, ], Osamps[, , g, ],
                        type=type, probRegion=probRegion)
    pts <- hp$pts
    hullPoints <- hp$hullPoints

    group_pts_angle <- c(group_pts_angle, pts[1, ])
    group_pts_eval <- c(group_pts_eval, pts[2, ])
    group_type <- c(group_type, rep(group_names[count], length(pts[2, ])))
    count <- count + 1

  }

  posterior_summaries <- tibble(angle = group_pts_angle, eval = group_pts_eval, Group = factor(group_type))

  p <- ggplot(posterior_summaries) +
    geom_point(aes(x=angle, y=eval, col=Group), alpha=alpha) +
    theme_bw(base_size=20) + xlim(c(-pi/2, pi/2)) + ylim(c(0, ymax))

  if(!is.null(hline))
    p <- p + geom_hline(yintercept=hline, lty=2)
  if(!legend) {
    p <- p + theme(legend.position = "none")
  } else{
    p <- p + theme(legend.position = "top", legend.title=element_blank())
  }
  if(!is.null(col_values))
    p <- p + scale_color_manual(values=col_values, alpha=alpha)

  p + ylab(ylab) + xlab(expression("angle, acos("~U[1]^T*V[1]~")")) +
      theme(legend.title=element_blank())

}


#' Title
#'
#' @param Vsub 
#' @param projected_samples 
#' @param obs_to_plot 
#' @param nlabeled 
#' @param legend 
#' @param label_size 
#' @param col_values 
#'
#' @return
#' @export covarianceBiplot
#'
#' @examples
covarianceBiplot <- function(Vsub, projected_samples, obs_to_plot=1:dim(projected_samples)[1],
                             nlabeled=20, legend=TRUE, label_size=2, col_values=NULL) {

  if(ncol(Vsub) != 2) {
    stop("Please provide 2-dimensional subspace")
  }

  if(is.null(rownames(Vsub)))
    rownames(Vsub) <- 1:nrow(Vsub)

  npos <- round(nlabeled/4)
  nneg <- round(nlabeled/4)

  xlimits <- c(-1.1, 1.1)*max(abs(Vsub))
  ylimits <- c(-1.1, 1.1)*max(abs(Vsub))

  p <- ggplot(as_tibble(Vsub), colnames=c("V1", "V2")) +
    geom_point(aes(x=V1, y=V2), size=0.5, col="light grey") +
    theme_bw(base_size=20) + xlim(xlimits) + ylim(ylimits) +
    theme(legend.text=element_text(size=15)) + xlab("") + ylab("")


  pos_x_indices <- order(Vsub[, 1], decreasing=TRUE)[1:npos]
  neg_x_indices <- order(Vsub[, 1], decreasing=FALSE)[1:nneg]
  pos_y_indices <- setdiff(order(Vsub[, 2], decreasing=TRUE), c(pos_x_indices, neg_x_indices))[1:npos]
  neg_y_indices <- setdiff(order(Vsub[, 2], decreasing=FALSE), c(pos_x_indices, pos_y_indices))[1:nneg]

  lambda_max <- lambda_min <- angle <- c()

  nobs <- length(obs_to_plot)
  for(k in obs_to_plot) {

    pmPsi <- apply(projected_samples[k, , , ], 1:2, mean)

    eigK <- eigen(pmPsi)
    lambda <- eigK$values
    evecs <- eigK$vectors

    maxIndex <- which.max(lambda)
    lamRatio <- lambda[maxIndex]/lambda[-maxIndex]

    lambda_max <- c(lambda_max, lambda[maxIndex])
    lambda_min <- c(lambda_min, lambda[-maxIndex])
    angle = c(angle, atan(evecs[2, maxIndex]/evecs[1, maxIndex]))

  }

  obs_names <- names(obs_to_plot)
  if(is.null(obs_names))
    obs_names = factor(1:nobs)

  ellipse_tibble <- tibble(lambda_max=lambda_max,
                           lambda_min=lambda_min,
                           angle=angle,
                           Group=obs_names)

  ellipse_tibble <- ellipse_tibble %>%
    mutate(sd_max = sqrt(lambda_max / max(lambda_max))) %>%
    mutate(sd_min = sqrt(lambda_min / max(lambda_max)))


  ## Plot Ellipses
  ## 95% contour
  p <- p + geom_ellipse(data=ellipse_tibble,
                        aes(x0 = 0, y0 = 0, a = sd_max*2*max(xlimits)*0.25,
                            b = sd_min*2*max(xlimits)*0.25,
                            angle = angle, color=Group), size=1.1)

  ## 50% contour
  p <- p + geom_ellipse(data=ellipse_tibble,
                        aes(x0 = 0, y0 = 0,
                            a = sd_max*0.67*max(xlimits)*0.25,
                            b = sd_min*0.67*max(xlimits)*.25,
                            angle = angle, color=Group), size=1.1)

  if(!is.null(col_values))
    p <- p + scale_color_manual(values=col_values)

  ## Add labels
  all_indices <- c(pos_x_indices, neg_x_indices, pos_y_indices, neg_y_indices)
  label_data <- tibble(x=Vsub[all_indices, 1], y=Vsub[all_indices, 2],
                       label=rownames(Vsub)[all_indices],
                       type=rep(c("pos_x", "neg_x", "pos_y", "neg_y"), each=nlabeled/4))

  p <- p + geom_point(data=label_data, aes(x=x, y=y), col="red", size=1.5)

  p <- p + geom_label_repel(
    data  = subset(label_data, type=="pos_x"),
    aes(x=x, y=y, label=label),
    nudge_x = 0.5,
    force=1,
    segment.size  = 0.5,
    segment.color = "grey50",
    direction     = "y",
    hjust         = 0,
    size = label_size
  ) +
    geom_label_repel(
      data = subset(label_data, type=="neg_x"),
      aes(x=x, y=y, label=label),
      force=1,
      nudge_x = -0.5,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction     = "y",
      hjust         = 1,
      size= label_size
    ) +
    geom_label_repel(
      data = subset(label_data, type=="pos_y"),
      aes(x=x, y=y, label=label),
      force=1,
      nudge_y = 0.05,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction = "both",
      size = label_size
    ) +
    geom_label_repel(
      data = subset(label_data, type=="neg_y"),
      aes(x=x, y=y, label=label),
      force=1,
      nudge_y = -0.05,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction     = "both",
      size=label_size
    )



  p + theme(legend.position = "top", legend.title=element_blank())

}

getHullPoints <- function(nsamps, pmPsi, OmegaSamps, Osamps, type="mag",
                          probRegion=0.95) {

  rot  <- eigen(pmPsi)$vectors
  ## rot  <- diag(2)

  PointsList <- lapply(1:nsamps, function(i) {
        LambdaSamp <- OmegaSamps[, i]/(1-OmegaSamps[, i])
        maxIndex <- which.max(LambdaSamp)

        if(type == "mag")
            yval <- LambdaSamp[maxIndex]
        else{
            yval <- LambdaSamp[maxIndex]/LambdaSamp[-maxIndex]
        }

        O1 <- Osamps[, maxIndex, i]
        angle <- atan(O1[2]/O1[1])

        O1_rot  <- rot  %*% O1
        angle_rot  <- atan(O1_rot[2]/O1_rot[1])

        c(angle, yval, angle_rot)

    })

    pts <- simplify2array(PointsList)
    allPts <- pts
    if(type == "logratio") {
        allPts[2, ] <- log2(allPts[2, ])
        pts[2, ] <- log2(pts[2, ])
    }

    numPtsToRemove <- round(nsamps*(1-probRegion))
    while(numPtsToRemove > 0) {
        hullPoints <- chull(pts[3, ], pts[2, ])
        if(length(hullPoints) > numPtsToRemove) {
            hullPoints <- sample(hullPoints, numPtsToRemove)
            pts <- pts[, -hullPoints]
            numPtsToRemove <- 0
        } else{
            pts <- pts[, -hullPoints]
            numPtsToRemove <- numPtsToRemove - length(hullPoints)
        }
    }



    hullPoints <- chull(pts[1, ], pts[2, ])

    list(allPts=allPts[1:2, ], pts=pts[1:2, ], hullPoints=hullPoints)

}
