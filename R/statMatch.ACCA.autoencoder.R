#' statMatch.ACCA.autoencoder
#' @description Initialize the Autoencoder using keras
#' @param n_var an integer, number of variables in the original space.
#' @param n_lat an integer, number of variables in the latent space.
#' @param n_hidden_layers an integer, number of hidden layers.
#' @param n_units a numeric vector, number of units in each hidden layer.
#' @param penL1 a double, value of the LASSO penalty.
#' @param penL2 a double, value of the RIDGE penalty.
#' @param lr a double, value of the learning rate.
#'
#' @return A keras model.
#'
#' @importFrom magrittr %>%
#' @importFrom keras layer_input layer_dense regularizer_l1_l2 keras_model compile optimizer_adam
#' @noRd
#'
statMatch.ACCA.autoencoder <- function(n_var, n_lat, n_hidden_layers, n_units, penL1, penL2, lr) {

  ### Define architecture of the autoencoder of X
  inputs <- layer_input(shape = n_var)
  # ENCODER X
  encoder <- inputs
  for (j in seq_len(n_hidden_layers)) {
    encoder <- encoder %>%
      layer_dense(
        units = n_units[j],
        activation = 'relu',
        activity_regularizer = regularizer_l1_l2(l1 = penL1, l2 = penL2)
      )
  }
  encoder <- encoder %>% layer_dense(units = n_lat)
  # DECODER X
  decoder <- encoder
  for (j in order(-seq_len(n_hidden_layers))) {
    decoder <- decoder %>%
      layer_dense(
        units = n_units[j],
        activation = 'relu',
        activity_regularizer = regularizer_l1_l2(l1 = penL1, l2 = penL2)
      )
  }
  decoder <- decoder %>% layer_dense(units = n_var)
  # AUTOENCODER X
  model <- keras_model(inputs = inputs, outputs = decoder)
  model <- model %>% compile(loss = "mean_squared_error",
                             weighted_metrics = list("mse"),
                             optimizer = optimizer_adam(learning_rate = lr))
  out <- list(
    model = model,
    inputs = inputs,
    encoder = encoder,
    decoder = decoder
  )

  return(out)
}
