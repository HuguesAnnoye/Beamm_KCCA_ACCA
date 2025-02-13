#' @title statMatch.MLP.phase2.model
#' @description Initialize the MMLP using keras (phase2)
#' @param n_X an integer, number of explanatory variables.
#' @param n_Y an integer, number of dependent variables.
#' @param n_hidden_layers an integer, number of hidden layers.
#' @param n_units a numeric vector, number of units in each hidden layer.
#' @param penL1 a double, value of the LASSO penalty.
#' @param penL2 a double, value of the RIDGE penalty.
#' @param learning_rate a double, value of the learning rate.
#'
#' @return A keras model.
#'
#' @importFrom magrittr %>%
#' @importFrom keras layer_input layer_dense regularizer_l1_l2 keras_model compile optimizer_adam
#' @noRd

statMatch.MLP.phase2.model <- function(n_X, n_Y, n_hidden_layers, n_units, penL1, penL2, learning_rate) {
  # Define architecture of the neural network
  inputs <- layer_input(shape = n_X)
  outputs <- inputs
  for (j in seq_len(n_hidden_layers)) {
    outputs <- outputs %>%
      layer_dense(
        units = n_units[j],
        activation = 'relu',
        activity_regularizer = regularizer_l1_l2(l1 = penL1, l2 = penL2)
      )
  }
  outputs <- outputs %>% layer_dense(units = n_Y)

  # Initial model
  model <- keras_model(inputs = inputs, outputs = outputs)
  model <- model %>% compile(loss = "mean_squared_error",
                             weighted_metrics = list("mse"),
                             optimizer = optimizer_adam(learning_rate = learning_rate))

  return(model)
}
