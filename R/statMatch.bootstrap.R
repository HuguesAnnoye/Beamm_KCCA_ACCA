# Peut-être à modifier à réfléchir. Surtout en ce qui concerne Pred_ZC
#' @title Statistical Matching using bootstrap
#' @description Make a boostrap statistical matching using the function statMatch
#' @param df.rec a data.frame, the recipient data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param opts a list, supplied using a function \code{statMatch.***.options()}, containing all the options.
#' @param boot_fold a positive integer, with number of folds
#' @param length_boot a positive integer, with length of each bootstrap fold
#' @param n_cores a positive integer, with number of cores
#' @param print.boot a logical value, print tuned parameter or not
#' @param method a character vector with the name of statmatch method that have to be used
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom dplyr intersect setdiff
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#'
#' @export

statMatch.bootstrap <- function(df.rec,
                                df.don,
                                don.weights = NULL,
                                boot_fold = 10,
                                length_boot = 1000,
                                opts = NULL,
                                n_cores = 3,
                                print.boot = TRUE,
                                names.CV = NULL,
                                names.NCV = NULL,
                                zero.constraints = NULL,
                                method = method) {
  if (is.null(names.CV)) {
    names.CV <- intersect(colnames(df.rec), colnames(df.don))
  }
  if (is.null(names.NCV)) {
    names.NCV <- setdiff(colnames(df.don), names.CV)
  }
  cl <- makeCluster(n_cores, outfile = "")
  registerDoParallel(cl)
  # out_boot <- foreach(fold = 1:boot_fold,.combine=list, .packages = c("kernlab", "Hmisc", "MFAg","BEAMM.statMatch")) %dopar% {
  out_boot <- foreach(fold = 1:boot_fold, .packages = c("kernlab", "Hmisc", "BEAMM.statMatch")) %dopar% {
    n_b <- sample(nrow(df.don), length_boot, replace = TRUE, prob = don.weights)
    df.don_b <- df.don[n_b, ]
    if (is.null(don.weights)) {
      weights_b <- NULL
    } else {
      weights_b <- NULL
    } # don.weights[n_b]}
    # out_res <-BEAMM.statMatch.CCAorKCCA(df.rec,df.don_b,
    out_res <- statMatch(df.rec = df.rec,
                         df.don = df.don_b,
      don.weights = weights_b,
      opts = opts,
      names.CV = names.CV,
      names.NCV = names.NCV,
      zero.constraints = zero.constraints,
      method = method
    )
    if (print.boot) {
      if (method == "KCCA") {
        message("Tune cat : h: ", out_res$phase1$tune$results$h[1], ", hx : ", out_res$phase1$tune$results$hx[1],", hy : ", out_res$phase1$tune$results$hy[1], " and g : ", out_res$phase1$tune$results$g[1])
        message("Tune : h : ", out_res$phase2$tune$results$h[1], ", hx : ", out_res$phase2$tune$results$hx[1], ", hy : ", out_res$phase2$tune$results$hy[1], " and g : ", out_res$phase2$tune$results$g[1])
      } else  if (method == "CCA"){
        message("Tune cat h :", out_res$phase1$tune$results$h[1])
        message("Tune h :", out_res$phase2$tune$results$h[1])
      }
    }
    return(out_res)
  }
  stopCluster(cl)
  if (print.boot) {
    message("Taille 1 :", length(out_boot))
    message("Taille 2 :", boot_fold)
  }
  Pred_boot <- NULL
  is.fact <- sapply(df.don[, names.NCV], is.factor)
  for (i in 1:length(names.NCV)) {
    Pred <- NULL
    Pred_ZC <- NULL
    vi <- NULL
    tune <- NULL
    tune_cat <- NULL
    for (fold in 1:boot_fold) {
      if (opts$print.details) message("Fold (b) : ", fold)
      out1 <- out_boot[[fold]]
      tune <- rbind(tune, out1$phase1$tune$results[1,])
      tune_cat <- rbind(tune_cat, out1$phase2$tune$results[1,])
      if (is.null(Pred)) {
        Pred <- out1$phase2$df.match[, names.NCV[i]]
        # level0 <- levels(Pred)
      } else {
        Pred <- data.frame(Pred, out1$phase2$df.match[, names.NCV[i]])
      }
      if (is.element(names.NCV[i], zero.constraints) == TRUE) {
        if (is.null(Pred_ZC)) {
          Pred_ZC <- out1$phase2$df.match[, paste0("ZC_", names.NCV[i])]
          Pred_ZC <- as.numeric(levels(Pred_ZC))[Pred_ZC]
        } else {
          Pred_ZC_tmp <- out1$phase2$df.match[, paste0("ZC_", names.NCV[i])]
          Pred_ZC_tmp <- as.numeric(levels(Pred_ZC_tmp))[Pred_ZC_tmp]
          Pred_ZC <- cbind(Pred_ZC, Pred_ZC_tmp)
        }
      }
      out1 <- NULL
    }
    if (is.fact[names.NCV[i]] == TRUE) {
      if (is.null(vi)) {
        ni <- sample(1:boot_fold, nrow(Pred), replace = TRUE)
        vi <- cbind(1:nrow(Pred), ni)
        Pred0 <- Pred[vi]
      } else {
        Pred0 <- Pred[vi]
      }
      Pred0 <- as.factor(Pred0)
      if (is.null(Pred_boot)) {
        Pred_boot <- Pred0
      } else {
        Pred_boot <- data.frame(Pred_boot, Pred0)
      }
    } else {
      Pred0 <- rowMeans(Pred)
      if (is.element(names.NCV[i], zero.constraints) == TRUE) {
        Pred_ZC0 <- rowMeans(Pred_ZC)
        Pred0[Pred_ZC0 >= 0.5] <- 0
      }
      if (is.null(Pred_boot)) {
        Pred_boot <- Pred0
      } else {
        Pred_boot <- data.frame(Pred_boot, Pred0)
      }
    }
  }
  colnames(Pred_boot) <- names.NCV
  out <- list(
    boot = out_boot,
    pred = Pred_boot,
    tune = tune,
    tune_cat = tune_cat,
    final = data.frame(df.rec, Pred_boot)
  )
  return(out)
}
