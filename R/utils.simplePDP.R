#' @title simple_PDP
#' @description Plot partial dependence plot for one or two variables without kernel
#' @param df_don a data.frame
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param names.ZC  character vector with the names of the zero constraint variables used for statistical matching.
#' @param names_weights a numeric vector, the sample weights of the donor data set.
#' @param grid_point number of grid points for estimation
#' @param method a character with the name of statistical matching method.
#' @param options a list, supplied using the function \code{statMatch.MLP.options()}, containing all the options.
#' @param input a character with the name of one or two common variables
#' @param prediction a character with the name of one non-common variable
#' @return A plot, ...
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate summarize group_by
#' @importFrom tidyr crossing
#' @importFrom ggplot2 ggplot scale_colour_gradient geom_point geom_line geom_tile theme_bw labs aes
#' @importFrom plotly plot_ly layout
#' @importFrom methods show
#'
#' @export


simple_PDP<- function(df_don, names.CV, names.NCV, names.ZC, names_weights, grid_point, method, options, input, prediction){

  #Check the dimension of input
  if(length(input)!=1 && length(input)!=2)
    stop("BEAMM.KCCAACCA :: input should include one or two variables")

  weights2 <- as.vector(df_don[[names_weights]])

  #Prepare data
  if(length(input)==1){ #PDP for 1 variable

    #Prepare data
    x_s <- select(df_don, input)   # grid where we want partial dependencies
    x_c <- select(df_don, -input, -names.NCV)  # other predictors

    # if the training dataset is large, use a subsample of x_c instead
    grid <- crossing(x_s, x_c)
    sum_weight <- sum(x_c[[names_weights]])

    #Statistical matching
    res <- statMatch(grid,
                     df_don,
                     don.weights = weights2,
                     names.CV = names.CV,
                     names.NCV = names.NCV,
                     zero.constraints = names.ZC,
                     opts = options,
                     method = method
    )

    #Add the prediction to grid
    au <- mutate(grid, .fitted = res$match[[prediction]])

    #Partial dependence plot for one variable
      #Group by x_s and calculate the mean of prediction
      au$.fitted <- au$.fitted*au[[names_weights]]
      pd <- au %>%
        group_by(au[[input]]) %>%
        summarize(yhat = sum(.fitted)/sum_weight)
      #Partial dependence plot
      p <- pd %>%
        ggplot(aes(`au[[input]]`, yhat, group = 1)) +
        geom_line(size = 1) +
        labs(title = "",
             y = c(prediction),
             x = c(input)) +
        theme_bw()
    } else if(length(input)==2){ #PDP for 2 variables

      if (is.numeric(df_don[[input[1]]])==TRUE){
        if (is.numeric(df_don[[input[2]]])==TRUE){
          grid_1 <- seq(from=min(df_don[[input[1]]])+0.001, to =max(df_don[[input[1]]]-0.001), length.out = grid_point)
          grid_2 <- seq(from=min(df_don[[input[2]]])+0.001, to =max(df_don[[input[2]]]-0.001), length.out = grid_point)

          x_s <- crossing(grid_1, grid_2)
          idx_cont <- c(1,2)
          } else {
            grid_1 <- seq(from=min(df_don[[input[1]]])+0.001, to =max(df_don[[input[1]]]-0.001), length.out = grid_point)
            grid_2 <- levels(df_don[[input[2]]])

            x_s <- crossing(grid_1, grid_2)
            idx_cont <- 1
            }
        } else if (is.numeric(df_don[[input[2]]])==TRUE){
          grid_1 <- levels(df_don[[input[1]]])
          grid_2 <- seq(from=min(df_don[[input[2]]])+0.001, to =max(df_don[[input[2]]]-0.001), length.out = grid_point)

          x_s <- crossing(grid_1, grid_2)
          idx_cont <- 2
          } else {
            x_s <- select(df_don, input)
            idx_cont <-0
            }
      names(x_s)[names(x_s) == "grid_1"] <- input[1]
      names(x_s)[names(x_s) == "grid_2"] <- input[2]

      x_c <- select(df_don, -input, -names.NCV)  # other predictors
      grid <- crossing(x_s, x_c)
      sum_weight <- sum(x_c[[names_weights]])

      #Statistical matching
      res <- statMatch(grid,
                       df_don,
                       don.weights = weights2,
                       names.CV = names.CV,
                       names.NCV = names.NCV,
                       zero.constraints = names.ZC,
                       opts = options,
                       method = method
      )

      #Add the prediction to grid
      au <- mutate(grid, .fitted = res$match[[prediction]])

      #Partial dependence plot for two variables

        #Group by x_s and calculate the mean of prediction

        au$.fitted2 <- au$.fitted*au[[names_weights]]
        pd <- au %>%
          group_by(au[[input[1]]], au[[input[2]]]) %>%
          summarize(yhat = sum(.fitted2)/sum_weight, .groups = 'drop')

        #Simple partial dependence plot
        if (length(idx_cont)==1){
          if (idx_cont!=0){
            #1 categorical variable and 1 continuous variable
            pd %>%
              ggplot(aes(`au[[input[1]]]`, `au[[input[2]]]`, yhat, color=yhat)) +
              geom_line(size=4) +
              scale_colour_gradient(low = "blue", high = "yellow") +
              labs(title = "Simple partial dependence plot",
                   y = c(input[2]),
                   x = c(input[1])) +
              theme_bw()
          }else{
            #2 categorical variables
            pd %>%
              ggplot(aes(`au[[input[1]]]`, `au[[input[2]]]`, yhat, color=yhat)) +
              geom_point(size=4) +
              scale_colour_gradient(low = "blue", high = "yellow") +
              labs(title = "Simple partial dependence plot",
                   y = c(input[2]),
                   x = c(input[1])) +
              theme_bw()
          }
        }else{
          #2 continuous variables
          pd %>%
            ggplot(aes(`au[[input[1]]]`, `au[[input[2]]]`, yhat, color=yhat)) +
            geom_tile(size=4) +
            scale_colour_gradient(low = "blue", high = "yellow") +
            labs(title = "Simple partial dependence plot",
                 y = c(input[2]),
                 x = c(input[1])) +
            theme_bw()
          }

    }else stop("BEAMM.KCCAACCA :: input should include one or two variables")

}
