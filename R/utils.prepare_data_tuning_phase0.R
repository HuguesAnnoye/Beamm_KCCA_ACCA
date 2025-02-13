prepare_data_tuning_phase0 <- function(df.don, names.CV, names.NCV, names.NCV_cat=NULL, names.NCV_cont=NULL, don.weights, opts, method = c("all", "MLP_phase1", "CCA"),
                                comp.mtx = F, names.CV.cat = NULL, comp.mtx.index = F) {
  method <- match.arg(method)
  n.don <- nrow(df.don)
  n.NCV <- length(names.NCV)
  if (is.null(names.NCV_cat))   names.NCV_cat <- names.NCV[unlist(lapply(df.don %>% select(all_of(names.NCV)), is.factor))]
  if (is.null(names.NCV_cont))   names.NCV_cont <- names.NCV[unlist(lapply(df.don %>% select(all_of(names.NCV)), !is.factor))]

  # Shuffle data before cross-validation (reorder matrices randomly)
  idx <- sample(n.don)
  df.don <- df.don[idx, ]
  don.weights <- don.weights[idx]
  # Tranform data.frame into matrix (categorical variables into dummies)
  if (method=="CCA"){
    mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)), n_1_levels = TRUE)
    mtx.don.NCV.raw <- df2mtx(df.don %>% select(all_of(names.NCV)), n_1_levels = TRUE)
    mtx.don.CV.raw_cont <- df2mtx(df.don %>% select(all_of(c(names.CV,names.NCV_cat))))
    mtx.don.NCV.raw_cat <- df2mtx(df.don %>% select(all_of(names.NCV_cat)), n_1_levels = FALSE)
    mtx.don.NCV.raw_cont <- df2mtx(df.don %>% select(all_of(names.NCV_cont)), n_1_levels = FALSE)
    mtx.don.CV.n <- df2mtx(df.don %>% select(all_of(names.CV)), n_1_levels = FALSE)
    mtx.don.NCV.n <- df2mtx(df.don %>% select(all_of(names.NCV)), n_1_levels = FALSE)
    names.NCV.raw <- colnames(mtx.don.NCV.raw)
  }else {
  mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)))
  mtx.don.CV.raw_cont <- df2mtx(df.don %>% select(all_of(c(names.CV,names.NCV_cat))))
  mtx.don.NCV.raw <- df2mtx(df.don %>% select(all_of(names.NCV)))
  mtx.don.NCV.raw_cat <- df2mtx(df.don %>% select(all_of(names.NCV_cat)))
  mtx.don.NCV.raw_cont <- df2mtx(df.don %>% select(all_of(names.NCV_cont)))
  names.NCV.raw <- colnames(mtx.don.NCV.raw)
  }
  # Prepare data for cross-validation and find compatible individuals (if required)
  folds <- cut(seq(1, n.don), breaks = opts$n_fold, labels = F)
  out <- list(raw_data = list(),
              weights = list())
  if (comp.mtx) out$comp.idx <- list()
  for (fold in seq_len(opts$n_fold)) {
    idx.TEST <- folds == fold
    idx.TRAIN <- !idx.TEST
    # Prepare donor's weights
    out$weights[[fold]] <- don.weights[idx.TRAIN]
    # Prepare CV and NCV
    if (method == "all") {
      out$raw_data[[fold]] <- list(mtx.don.CV.raw = mtx.don.CV.raw[idx.TRAIN, , drop = F],
                                   mtx.don.CV.raw_cont = mtx.don.CV.raw_cont[idx.TRAIN, , drop = F],
                                   mtx.don.NCV.raw = mtx.don.NCV.raw[idx.TRAIN, , drop = F],
                                   mtx.don.NCV.raw_cat = mtx.don.NCV.raw_cat[idx.TRAIN, , drop = F],
                                   mtx.don.NCV.raw_cont = mtx.don.NCV.raw_cont[idx.TRAIN, , drop = F],
                                   mtx.rec.CV.raw = mtx.don.CV.raw[idx.TEST, , drop = F],
                                   mtx.rec.CV.raw2 = mtx.don.CV.raw[idx.TEST, , drop = F]
      )
    } else if (method == "MLP_phase1") {
      out$raw_data[[fold]] <- list(mtx.don.CV.raw = mtx.don.CV.raw[idx.TRAIN, , drop = F],
                                   mtx.rec.CV.raw = mtx.don.CV.raw[idx.TEST, , drop = F])
      out$raw_data[[fold]]$mtx.don.NCV.raw <- list()
      for (i in seq_len(n.NCV)) {
        sel <- startsWith(names.NCV.raw, names.NCV[i])
        out$raw_data[[fold]]$mtx.don.NCV.raw[[i]] <- mtx.don.NCV.raw[idx.TRAIN, sel, drop = F]
      }
    } else if (method == "CCA") {
      out$raw_data[[fold]] <- list(mtx.don.CV.raw = mtx.don.CV.raw[idx.TRAIN, , drop = F],
                                   mtx.don.CV.raw_cont = mtx.don.CV.raw_cont[idx.TRAIN, , drop = F],
                                   mtx.don.NCV.raw = mtx.don.NCV.raw[idx.TRAIN, , drop = F],
                                   mtx.don.NCV.raw_cat = mtx.don.NCV.raw_cat[idx.TRAIN, , drop = F],
                                   mtx.don.NCV.raw_cont = mtx.don.NCV.raw_cont[idx.TRAIN, , drop = F],
                                   mtx.rec.CV.raw = mtx.don.CV.raw[idx.TEST, , drop = F],
                                   mtx.rec.CV.raw2 = mtx.don.CV.raw[idx.TEST, , drop = F],
                                   mtx.don.NCV.n = mtx.don.NCV.n[idx.TRAIN, , drop = F],
                                   mtx.don.CV.n = mtx.don.CV.n[idx.TRAIN, , drop = F],
                                   mtx.rec.CV.n = mtx.don.CV.n[idx.TEST, , drop = F])
    }
    # Scale CV/NCV (if necessary)
    if (opts$scaling != "no")
      out$raw_data[[fold]] <- statMatch.scale.phase0(data = out$raw_data[[fold]],
                                              wt = out$weights[[fold]],
                                              method = opts$scaling,
                                              scale_NCV = (method == "all"))
    # Compute compatibility matrix between donor and receiver
    if (comp.mtx)
      out$comp.idx[[fold]] <- compute_mtx_compatibilities(df.don[idx.TEST, ],
                                                          df.don[idx.TRAIN, ],
                                                          names.CV.cat,
                                                          opts$print.details,
                                                          index = comp.mtx.index)
  }
  out$mtx.don.NCV.raw <- mtx.don.NCV.raw
  if (method == "CCA") { out$mtx.don.NCV.n <- mtx.don.NCV.n }
  return(out)
}
