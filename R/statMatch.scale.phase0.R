#' @title statMatch.scale
#' @description Scale a statMatch_data object
#' @param data a numeric matrix.
#' @param wt a numeric vector, sample weights.
#' @param method a character string, the method for scaling.
#' @param scale_NCV a logical, whether or not the non-common variables are scaled.
#' @return A numeric matrix.
#' @noRd

statMatch.scale.phase0 <- function(data, wt, method = c("z-score", "min-max"), scale_NCV = T) {
  method <- match.arg(method)
  if (method == "z-score") {
    data$mtx.don.CV.raw <- scale_mtx(data$mtx.don.CV.raw, wt)
    data$mtx.don.CV.raw_cont <- scale_mtx(data$mtx.don.CV.raw_cont, wt)
    if (scale_NCV) {
      data$mtx.don.NCV.raw <- scale_mtx(data$mtx.don.NCV.raw, wt)
      data$mtx.don.NCV.raw_cat <- scale_mtx(data$mtx.don.NCV.raw_cat, wt)
      data$mtx.don.NCV.raw_cont <- scale_mtx(data$mtx.don.NCV.raw_cont, wt)
    }
    data$mtx.rec.CV.raw <- scale_mtx(data$mtx.rec.CV.raw,
                                     center = attr(data$mtx.don.CV.raw, "scaled:center"),
                                     scale = attr(data$mtx.don.CV.raw, "scaled:scale"))
  } else if (method == "min-max") {
    data$mtx.don.CV.raw <- normalize_mtx(data$mtx.don.CV.raw)
    data$mtx.don.CV.raw_cont <- scale_mtx(data$mtx.don.CV.raw_cont, wt)
    if (scale_NCV) {
      data$mtx.don.NCV.raw <- normalize_mtx(data$mtx.don.NCV.raw)
      data$mtx.don.NCV.raw_cat <- normalize_mtx(data$mtx.don.NCV.raw_cat)
      data$mtx.don.NCV.raw_cont <- normalize_mtx(data$mtx.don.NCV.raw_cont)
      }
    data$mtx.rec.CV.raw <- normalize_mtx(data$mtx.rec.CV.raw,
                                 min = attr(data$mtx.don.CV.raw, "min"),
                                 max = attr(data$mtx.don.CV.raw, "max"))
  }
  return(data)
}
