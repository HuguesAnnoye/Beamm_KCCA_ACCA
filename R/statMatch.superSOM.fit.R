#' @title statMatch.superOM.fit
#' @description Fit superOM for statistical matching
#' @param data a list containing three data.frames: \code{mtx.don.CV.raw}, \code{mtx.don.NCV.raw} and \code{mtx.rec.CV.raw}.
#' @param grid_xdim an integer, the number of columns of the grid of neurons.
#' @param grid_ydim an integer, the number of rows of the grid of neurons.
#' @param weight.CV a double, the weight of the superOM layer containing the common variables.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param opts a list, supplied using the function \code{statMatch.superOM.options()}, containing all the options.
#' @return An object of class "kohonen".
#' @noRd

statMatch.superOM.fit <- function(data, grid_xdim, grid_ydim, weight.CV, don.weights, opts){
  list.raw <- list(layer1 = data$mtx.don.CV.raw,
                   layer2 = data$mtx.don.NCV.raw)
  def.grid <- kohonenBEAMM::somgrid(xdim = grid_xdim,
                                    ydim = grid_ydim,
                                    topo = opts$topo,
                                    neighbourhood.fct = opts$neighbourhood.fct,
                                    toroidal = opts$toroidal)
  fit <- kohonenBEAMM::supersom(list.raw,
                                grid = def.grid,
                                user.weights = c(weight.CV, 1 - weight.CV),
                                obs.weights = don.weights,
                                rlen = opts$rlen,
                                maxNA.fraction = opts$maxNA.fraction,
                                dist.fcts = opts$dist.fcts,
                                normalizeDataLayers = FALSE,
                                keep.data = TRUE,
                                mode = "batch")
  return(fit)
}
