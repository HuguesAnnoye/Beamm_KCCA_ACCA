#' @title kernel_1var_PDP
#' @description Plot partial dependence plot for one variable with kernel optimized by operating only with matrices
#' @param df_don a data.frame
#' @param num_of_h number of h for tuning part
#' @param h_min minimum value of h for tuning part
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param names.ZC  character vector with the names of the zero constraint variables used for statistical matching.
#' @param names_weights a numeric vector, the sample weights of the donor data set.
#' @param method a character with the name of statistical matching method.
#' @param options a list, supplied using the function \code{statMatch.MLP.options()}, containing all the options.
#' @param input a character with the name of one common variables
#' @param prediction a character with the name of one non-common variable
#' @param res a logical to check whether the statistical matching is conducted
#' @return A plot, ...
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate summarize group_by
#' @importFrom tidyr crossing
#' @importFrom ggplot2 ggplot geom_line theme_bw aes labs
#' @importFrom plotly plot_ly layout
#' @importFrom methods show


kernel_1var_PDP <- function(df_don, num_of_h = 10, h_min=0.5, names.CV, names.NCV, names.ZC, names_weights, method, options, input, prediction, res = NULL){

  #Check the dimension of input
  if(length(input)!=1)
    stop("BEAMM.KCCAACCA :: input should include one variable")

  weights2 <- df_don[[names_weights]]

  #Prepare data. Three possibilities: 1. The prediction variable already belongs to the common variables.
  # 2. The prediction variable does not belong to the common variables and we need to run the matching
  # 3. The prediction variable does not belong to the common variables but the statistical matching has already been performed on it.
  if(is.null(res)){
    if(input %in% names.CV && prediction %in% names.CV){
      x_s <- select(df_don, input)
      au <- mutate(x_s, .fitted = df_don[[prediction]],names_weights=df_don[[names_weights]])
    } else {
      # Initialization
      folds <- cut(seq(1, nrow(df_don)), breaks = 5, labels = F)
      res_matching <- NULL
      best_models <- NULL
      best_models_cat <- NULL

      # Creation of a donnor and a receiver data set
      fold <- 1
      idx.TEST <- folds == fold
      idx.TRAIN <- !idx.TEST
      df.don <- df_don[idx.TRAIN, c(names.CV, names.NCV)]
      df.rec <- df_don[idx.TEST, names.CV]
      df.rec[, names.NCV] <- NULL

      res <- statMatch(df.rec,
                       df.don,
                       xdon.weights = weights2[idx.TRAIN],
                       names.CV = names.CV,
                       names.NCV = names.NCV,
                       zero.constraints = names.ZC,
                       opts = options,
                       method = method
      )
      grid <- select(res$match, input)
      au <-  mutate(grid, .fitted = res$match[[prediction]])
    }
  }else {
    grid <- select(res$match, input)
    au <- mutate(grid, .fitted = res$match[[prediction]])
  }

  y_hat <- select(df_don, prediction)


  #Tuning bandwidth
  if (is.numeric(df_don[[input]])==TRUE){
    pred <- NULL
    bandwidth <- NULL
    h_max <- (max(x_s)-min(x_s))/2
    X <- as.numeric(unlist(x_s))
    Y <- as.numeric(unlist(y_hat))
    h <- seq(from=h_min, to = h_max, length.out = num_of_h)
    for (i in 1:length(h)) {
      message("h_", i, " of ", length(h), ", h =", h[i])
      kern<-exp(-(((outer(X,X,FUN = "-"))/h[i])^2))
      diag(kern)<-0
      num<-kern%*%(weights2*Y)
      den<-kern%*%weights2
      obj_func<-sum(weights2*(Y-num/den)^2)
      pred <- rbind(pred, obj_func)
      bandwidth <- rbind(bandwidth, h[i])
    }
    result <- cbind(bandwidth, pred)
    # Order results in increasing order
    idx <- order(result[,2])
    h_final <- result[idx[1],1]
    message("h_final is ", h_final)
  }


  #Calculate weighted mean of prediction
  if (is.numeric(df_don[[input]])==TRUE){
    X <- au[[input]]
    Y <- au$.fitted

    kern <- exp(-((outer(X,X,FUN = "-"))/h_final)^2)
    #diag(kern) <- 0 Not performed at this step.
    num <- kern%*%(weights2*Y)
    den <- kern%*%weights2

    au$.fitted2 <- num/den

    pd <- au %>%
      group_by(au[[input]]) %>%
      summarize(yhat = mean(.fitted2))
  } else {
    #Group by x_s and calculate the mean of prediction
    sum_weight <- sum(au$names_weights)
    au$.fitted <- au$.fitted*au$names_weights
    pd <- au %>%
      group_by(au[[input]]) %>%
      summarize(yhat = sum(.fitted)/sum_weight)
  }

  #Partial dependence plot
  p<- pd %>%
    ggplot(aes(x=`au[[input]]`, y=yhat, group = 1)) +
    geom_line(size = 1) +
    labs(title = "Kernel partial dependence plot",
         y = c(prediction),
         x = c(input)) +
    theme_bw()
}


