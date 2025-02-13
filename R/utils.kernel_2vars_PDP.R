#' @title kernel_2vars_PDP
#' @description Plot partial dependence plot for one or two variables with kernel
#' @param df_don a data.frame
#' @param h_final preset value for h instead of tuning
#' @param num_of_h number of h for tuning part
#' @param h_min minimum value of h for tuning part
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param names.ZC  character vector with the names of the zero constraint variables used for statistical matching.
#' @param names_weights a numeric vector, the sample weights of the donor data set.
#' @param grid_point number of grid points for estimation
#' @param method a character with the name of statistical matching method.
#' @param options a list, supplied using the function \code{statMatch.MLP.options()}, containing all the options.
#' @param input a character with the name of one or two common variables
#' @param prediction a character with the name of one non-common variable
#' @param res a logical to check whether the statistical matching is conducted
#' @return A plot, ...
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate summarize group_by
#' @importFrom tidyr crossing
#' @importFrom ggplot2 ggplot scale_colour_gradient geom_point geom_tile geom_line theme_bw aes labs
#' @importFrom plotly plot_ly layout
#' @importFrom methods show
#' @importFrom readr write_csv


kernel_2vars_PDP <-
  function(df_don,
           h_final = NULL,
           num_of_h = 10,
           h_min = 0.5,
           names.CV,
           names.NCV = NULL,
           names.ZC,
           names_weights,
           grid_point,
           method,
           options,
           input,
           prediction,
           res = NULL) {
    #Check the dimension of input
    if (length(input) != 2)
      stop("BEAMM.KCCAACCA :: input should include two variables")
    if (is.null(names_weights)) {
      weights2 <- matrix(1, ncol = 1, nrow = nrow(df_don))
    } else{
      weights2 <- df_don[[names_weights]]
    }

    #Statistical matching
    if (is.null(res)) {
      if (input[1] %in% names.CV && input[2] %in% names.CV) {
        if (prediction %in% names.CV) {
          grid <- select(df_don, input)
          au <- mutate(grid, .fitted = df_don[[prediction]])
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

          res <- statMatch(
            df.rec,
            df.don,
            don.weights = weights2,
            names.CV = names.CV,
            names.NCV = names.NCV,
            zero.constraints = names.ZC,
            opts = options,
            method = method
          )
          grid <- select(res$match, input)
          au <-  mutate(grid, .fitted = res$match[[prediction]])
        }
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

        res <- statMatch(
          df.rec,
          df.don,
          don.weights = weights2,
          names.CV = names.CV,
          names.NCV = names.NCV,
          zero.constraints = names.ZC,
          opts = options,
          method = method
        )
        grid <- select(res$match, input)
        au <-  mutate(grid, .fitted = res$match[[prediction]])
      }
    } else {
      grid <- select(res$match, input)
      au <- mutate(grid, .fitted = res$match[[prediction]])
    }

    #Tuning phase
    if (is.numeric(df_don[[input[1]]]) == TRUE) {
      if (is.numeric(df_don[[input[2]]]) == TRUE) {
        idx_cont <- c(1, 2)
        if (length(h_final) != 2 | !is.numeric(h_final)) {
          message("No valid h given. Tuning...")
          h_final <-
            PDP_tuning_2vars(
              df_don,
              names.NCV = NULL,
              input,
              prediction,
              num_of_h = num_of_h,
              h_min = h_min,
              names_weights
            )
        }
      } else {
        idx_cont <- 1
        if (!is.numeric(h_final) | length(h_final) != 1) {
          message("No valid h given. Tuning...")
          h_final <-
            PDP_tuning_1var(
              df_don,
              names.NCV = NULL,
              input[1],
              prediction,
              num_of_h = num_of_h,
              h_min = h_min,
              names_weights
            )
        }
      }
    } else if (is.numeric(df_don[[input[2]]]) == TRUE) {
      idx_cont <- 2
      if (!is.numeric(h_final) | length(h_final) != 1) {
        message("No valid h given. Tuning...")
        h_final <-
          PDP_tuning_1var(
            df_don,
            names.NCV = NULL,
            input[2],
            prediction,
            num_of_h = num_of_h,
            h_min = h_min,
            names_weights
          )
      }
    } else
      idx_cont <- 0

    #Construct partial dependence plot & frequency plot
    if (length(idx_cont) == 1) {
      if (idx_cont != 0) {
        #1 categorical variable and 1 continuous variable
        #Create grid
        grid_cont <-
          seq(
            from = min(au[[input[idx_cont]]]) + 0.001,
            to = max(au[[input[idx_cont]]] - 0.001),
            length.out = grid_point
          )
        grid_cat <- levels(au[[input[-idx_cont]]])

        pd_graph <- crossing(grid_cont, grid_cat)

        #Calculate kernelized weighted mean
        pd_graph$.fitted1 <- 1
        density <- 1
        for (i in 1:nrow(pd_graph)) {
          X_cat <- pd_graph$grid_cat[i]
          X_cont <- pd_graph$grid_cont[i]
          #Y <- au$.fitted[i]

          au_new <- cbind(au,weights2)
          au_new <- au_new %>%
            subset(au_new[[input[-idx_cont]]] == X_cat)

          kern <- exp(-((X_cont - au_new[[input[idx_cont]]]) / h_final) ^
                        2)
          num <- sum(kern * au_new$.fitted * au_new$weights2)
          den <- sum(kern * au_new$weights2)
          density[i] <- den
          pd_graph$.fitted1[i] <- num / den
        }
        #Group by x_s
        pd <- pd_graph %>%
          group_by(pd_graph$grid_cont, pd_graph$grid_cat) %>%
          dplyr::summarize(yhat = mean(.fitted1), .groups = 'drop')

        fig_pdp_2d <- pd %>%
          ggplot(aes(`pd_graph$grid_cont`, `pd_graph$grid_cat`, yhat, color =
                       yhat)) +
          geom_line(size = 3) +
          scale_colour_gradient(low = "blue", high = "yellow") +
          labs(title = "Kernel partial dependence plot",
               y = c(input[2]),
               x = c(input[1])) +
          theme_bw()

        fig_pdp_3d <- pd %>%
          plot_ly(
            x = pd_graph$grid_cont,
            y = pd_graph$grid_cat,
            z = pd_graph$.fitted1,
            color = density,
            type = "scatter3d",
            mode = "markers"
          )
        fig_pdp_3d <- fig_pdp_3d %>% add_markers()

        fig_pdp_3d <-
          fig_pdp_3d %>% layout(scene = list(
            xaxis = list(title = input[idx_cont]),
            yaxis = list(title = input[-idx_cont]),
            zaxis = list(title = prediction)
          ))
        # fig3d <- fig3d %>% add_markers

        #Frequency table
        au$input_cont <- 1
        d <-
          c(seq(
            from = min(au[[input[idx_cont]]]),
            to = max(au[[input[idx_cont]]]),
            by = (max(au[[input[idx_cont]]]) - min(au[[input[idx_cont]]])) / nrow(df_don) ^
              (2 / 5)
          ), max(au[[input[idx_cont]]]))
        au$input_cont <-
          cut(au[[input[idx_cont]]], breaks = d, include.lowest = TRUE)

        freq <- au %>%
          group_by(au$input_cont, au[[input[-idx_cont]]]) %>%
          dplyr::summarize(n = n(), .groups = 'drop')
        write_csv(freq, "frequency table.csv")
        #Plot of frequency table
        freq$input_cont <- as.numeric(freq$`au$input_cont`)
        freq$input_cont <- min(au[[input[idx_cont]]]) +
          (freq$input_cont / max(freq$input_cont)) * (max(au[[input[idx_cont]]]) - min(au[[input[idx_cont]]]))


        fig_freq <- freq %>%
          ggplot(aes(`au[[input[-idx_cont]]]`, `input_cont`, n, color = n)) +
          geom_line(size = 3) +
          scale_colour_gradient(low = "blue", high = "yellow") +
          labs(title = "Frequency table kernel PDP",
               y = c(input[[idx_cont]]),
               x = c(input[[-idx_cont]])) +
          theme_bw()
        #View(freq)


      }
      else {
        #2 categorical variables
        au$.fitted2 <- au$.fitted * weigths2
        pd <- au %>%
          group_by(au[[input[1]]], au[[input[2]]]) %>%
          dplyr::summarize(yhat = sum(.fitted2) / sum(weights2),
                           .groups = 'drop')
        #Plot partial dependence plot
        fig_pdp_2d <- pd %>%
          ggplot(aes(`au[[input[1]]]`, `au[[input[2]]]`, yhat, color = yhat)) +
          geom_point(alpha = 0.5, size = 3) +
          scale_colour_gradient(low = "blue", high = "yellow") +
          labs(title = "Kernel partial dependence plot",
               y = c(input[2]),
               x = c(input[1])) +
          theme_bw()

        fig_pdp_3d <- pd %>%
          plot_ly(
            x = pd_graph$grid_cont,
            y = pd_graph$grid_cat,
            z = pd_graph$.fitted1,
            color = `au[[input[1]]]`,
            type = "scatter3d",
            mode = "markers"
          )
        fig_pdp_3d <- fig_pdp_3d %>% add_markers()

        fig_pdp_3d <- fig_pdp_3d %>%
          layout(scene = list(
            xaxis = list(title = input[idx_cont]),
            yaxis = list(title = input[-idx_cont]),
            zaxis = list(title = prediction)
          ))
        #Frequency table
        freq <- au %>%
          group_by(au[[input[1]]], au[[input[2]]]) %>%
          dplyr::summarize(n = n(), .groups = 'drop')
        #write_csv(freq, "frequency table.csv")
        #Plot of frequency table
        fig_freq <- freq %>%
          ggplot(aes(`au[[input[1]]]`, `au[[input[2]]]`, n, color = n)) +
          geom_point(size = 3) +
          scale_colour_gradient(low = "blue", high = "yellow") +
          labs(title = "Frequency table kernel PDP",
               y = c(input[2]),
               x = c(input[1])) +
          theme_bw()
        #View(freq)
      }
    } else {
      #2 continuous variables
      #Create grid
      grid_1 <-
        seq(
          from = min(au[[input[1]]]) + 0.001,
          to = max(au[[input[1]]] - 0.001),
          length.out = grid_point
        )
      grid_2 <-
        seq(
          from = min(au[[input[2]]]) + 0.001,
          to = max(au[[input[2]]] - 0.001),
          length.out = grid_point
        )

      pd_graph <- crossing(grid_1, grid_2)
      #Calculate kernelized weighted mean
      pd_graph$.fitted1 <- 1
      density <- 1
      for (i in 1:nrow(pd_graph)) {
        X1 <- pd_graph$grid_1[i]
        X2 <- pd_graph$grid_2[i]
        #Y <- au$.fitted[i]

        kern <-
          exp((-((X1 - au[[input[1]]]) / h_final[1]) ^ 2) + (-((X2 - au[[input[2]]]) /
                                                                 h_final[2]) ^ 2))
        num <- sum(kern * au$.fitted * weights2)
        den <- sum(kern * weights2)
        density[i] <- den
        pd_graph$.fitted1[i] <- num / den
      }
      #Group by x_s
      pd <- pd_graph %>%
        group_by(pd_graph$grid_1, pd_graph$grid_2) %>%
        dplyr::summarize(yhat = mean(.fitted1), .groups = 'drop')

      fig_pdp_2d <- pd %>%
        ggplot(aes(`pd_graph$grid_1`, `pd_graph$grid_2`, yhat, color = yhat)) +
        geom_tile(size = 4) +
        scale_colour_gradient(low = "blue", high = "yellow") +
        labs(title = "Kernel partial dependence plot",
             y = c(input[2]),
             x = c(input[1])) +
        theme_bw()

      fig_pdp_3d <- pd %>%
        plot_ly(
          x = pd_graph$grid_1,
          y = pd_graph$grid_2,
          z = pd_graph$.fitted1,
          color = density,
          type = "scatter3d",
          mode = "markers"
        )
      fig_pdp_3d <- fig_pdp_3d %>% add_markers()

      fig_pdp_3d <-
        fig_pdp_3d %>% layout(scene = list(
          xaxis = list(title = input[1]),
          yaxis = list(title = input[2]),
          zaxis = list(title = prediction)
        ))
      #Frequency table
      au$input_1 <- 1
      au$input_2 <- 1
      d1 <-
        c(seq(
          from = min(au[[input[1]]]),
          to = max(au[[input[1]]]),
          by = (max(au[[input[1]]]) - min(au[[input[1]]])) / nrow(df_don) ^ (2 /
                                                                               5)
        ), max(au[[input[1]]]))
      d2 <-
        c(seq(
          from = min(au[[input[2]]]),
          to = max(au[[input[2]]]),
          by = (max(au[[input[2]]]) - min(au[[input[2]]])) / nrow(df_don) ^ (2 /
                                                                               5)
        ), max(au[[input[2]]]))

      au$input_1 <-
        cut(au[[input[1]]], breaks = d1, include.lowest = TRUE)
      au$input_2 <-
        cut(au[[input[2]]], breaks = d2, include.lowest = TRUE)

      freq <- au %>%
        group_by(au$input_1, au$input_2) %>%
        dplyr::summarize(n = n(), .groups = 'drop')
      write_csv(freq, "frequency table.csv")
      #Plot of frequency table
      freq$input_1 <- as.numeric(freq$`au$input_1`)
      freq$input_1 <- min(au[[input[1]]]) +
        (freq$input_1 / max(freq$input_1)) * (max(au[[input[1]]]) - min(au[[input[1]]]))

      freq$input_2 <- as.numeric(freq$`au$input_2`)
      freq$input_2 <- min(au[[input[2]]]) +
        (freq$input_2 / max(freq$input_2)) * (max(au[[input[2]]]) - min(au[[input[2]]]))


      fig_freq <- freq %>%
        ggplot(aes(`input_1`, `input_2`, n, color = n)) +
        geom_tile(size = 3) +
        scale_colour_gradient(low = "blue", high = "yellow") +
        labs(title = "Frequency table kernel PDP",
             y = c(input[2]),
             x = c(input[1])) +
        theme_bw()

    }
    return(list(fig_pdp_2d, fig_pdp_3d, fig_freq))
  }

PDP_tuning_1var <-
  function(df,
           names.NCV,
           input,
           prediction,
           num_of_h = 10,
           h_min = 0.5,
           names_weights) {
    # Prepare data
    mt <- select(df,-names.NCV)
    mt <- cbind(mt, prediction = df[[prediction]])

    x_s <- select(mt, input)
    y_hat <- mt[, "prediction"]

    if (is.null(names_weights)) {
      weight <- matrix(1, ncol = 1, nrow = nrow(df))
    } else {
      weight <- mt[, names_weights]
    }

    message("\n Tuning h ...")

    # Perform cross-validation
    pred <- NULL
    bandwidth <- NULL
    h_max <- (max(x_s) - min(x_s)) / 2
    h <- seq(from = h_min,
             to = h_max,
             length.out = num_of_h)
    Y <- as.numeric(unlist(y_hat))
    X <- as.numeric(unlist(x_s))
    for (i in 1:length(h)) {
      message("h_", i, " of ", length(h), "=", h[i])
      obj_func <- 0
      kern <- exp(-(((
        outer(X, X, FUN = "-")
      ) / h[i]) ^ 2))
      diag(kern) <- 0
      num <- kern %*% (weight * Y)
      den <- kern %*% weight
      obj_func <- sum(weight * (Y - num / den) ^ 2)
      pred <- rbind(pred, obj_func)
      bandwidth <- rbind(bandwidth, h[i])
    }
    result <- cbind(bandwidth, pred)
    # Order results in increasing order
    idx <- order(result[, 2])
    h_final <- result[idx[1], 1]
    return(h_final)
  }

PDP_tuning_2vars <-
  function(df,
           names.NCV,
           input,
           prediction,
           num_of_h = 10,
           h_min = 0.5,
           names_weights) {
    # Prepare data
    mt <- select(df,-names.NCV)
    mt <- cbind(mt, prediction = df[[prediction]])

    #Prepare data
    x_s <-
      select(mt, input)   # grid where we want partial dependencies
    y_hat <- mt$prediction
    if (is.null(names_weights)) {
      weight <- matrix(1, ncol = 1, nrow = nrow(df))
    } else {
      weight <- mt[, names_weights]
    }

    message("\n Tuning h ...")

    # Perform cross-validation
    pred <- NULL
    bandwidth <- NULL

    #Interval of bandwidth
    h1_max <- (max(x_s[, 1]) - min(x_s[, 1])) / 2
    h1 <- seq(from = h_min,
              to = h1_max,
              length.out = num_of_h)
    h2_max <- (max(x_s[, 2]) - min(x_s[, 2])) / 2
    h2 <- seq(from = h_min,
              to = h2_max,
              length.out = num_of_h)
    for (i in 1:length(h1)) {
      message("h1[", i, "]=", h1[i])
      for (j in 1:length(h2)) {
        message("h2[", j, "]=", h2[j])
        kern <-
          exp(-((outer(
            x_s[, 1], x_s[, 1], FUN = "-"
          ) / h1[i]) ^ 2) - ((outer(
            x_s[, 2], x_s[, 2], FUN = "-"
          ) / h2[i]) ^ 2))
        diag(kern) <- 0
        num <- kern %*% (weight * y_hat)
        den <- kern %*% weight
        obj_func <- sum(weight * (y_hat - num / den) ^ 2)
        pred <- rbind(pred, obj_func)
        bandwidth <- rbind(bandwidth, list(h1[i], h2[j]))
      }
    }
    result <- cbind(bandwidth, pred)
    # Order results in increasing order
    idx <- order(as.numeric(result[, 3]))
    h1_final <- as.numeric(result[idx[1], 1])
    h2_final <- as.numeric(result[idx[1], 2])
    h_final <- c(h1_final, h2_final)
    return(h_final)
  }
