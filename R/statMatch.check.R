statMatch.check <- function(df.rec, df.don, names.CV, names.NCV){

  if (!is.data.frame(df.rec))
    stop("BEAMM.statMatch :: input 'df.rec' should be a data.frame")

  if (!is.data.frame(df.don))
    stop("BEAMM.statMatch :: input 'df.don' should be a data.frame")

  if (length(intersect(colnames(df.rec), colnames(df.don))) == 0)
    stop("BEAMM.statMatch :: 'df.rec' and 'df.don' do not have any variable in common.
         Statistical matching cannot be performed without common variables.")

  if (length(intersect(names.CV, names.NCV)) > 0)
    stop("BEAMM.statMatch :: the following variables cannot be specified in both
         names.CV and names.NCV : ", paste(intersect(names.CV, names.NCV),collapse = ", "),".")

  if (!is.null(names.CV)){
    if (!is.character(names.CV))
      stop("BEAMM.statMatch :: input 'names.CV' should be a character vector")
    wrong.CV.rec <- !names.CV %in% colnames(df.rec)
    if (any(wrong.CV.rec))
      stop("BEAMM.statMatch :: the following common variable(s) are not available in 'df.rec' : ",
           paste(names.CV[wrong.CV.rec],collapse = ", "),".")
    wrong.CV.don <- !names.CV %in% colnames(df.don)
    if (any(wrong.CV.don))
      stop("BEAMM.statMatch :: the following common variable(s) are not available in 'df.don' : ",
           paste(names.CV[wrong.CV.don],collapse = ", "),".")
  }

  if (!is.null(names.NCV)){
    if (!is.character(names.NCV))
      stop("BEAMM.statMatch :: input 'names.NCV' should be a character vector")
    wrong.NCV.don <- names.NCV[!names.NCV %in% colnames(df.don)]
    if (length(wrong.NCV.don)>0)
      stop("BEAMM.statMatch :: the following non-common variable(s) are not available in 'df.don' : ",
           paste(wrong.NCV.don,collapse = ", "),".")
  }

  for (varname in names.CV){
    if (class(df.don[, varname, drop = T]) != class(df.rec[, varname, drop = T]))
      stop("BEAMM.statMatch :: the class of variable ", varname, " is different in df.rec and df.don.")
    if (!setequal(levels(df.don[, varname, drop = T]), levels(df.rec[, varname, drop = T])))
      stop("BEAMM.statMatch :: the levels of the factor variable ", varname, " are different in df.rec and df.don.")
  }

}

check.KCCA.options <- function(opts) {
  if (is.null(opts))
    opts <- statMatch.KCCA.options()
  else if (!inherits(opts, "KCCA.options"))
    stop("Argument 'opts' should be an object of class 'KCCA.options'")
  return(opts)
}

check.superOM.options <- function(opts) {
  if (is.null(opts))
    opts <- statMatch.superOM.options()
  else if (!inherits(opts, "superOM.options"))
    stop("Argument 'opts' should be an object of class 'superOM.options'")
  return(opts)
}

check.options <- function(opts, method) {
  cl  <- paste0(method, ".options")
  if (is.null(opts))
    opts <- eval(parse(text = paste0("statMatch.", cl, "()")))
  else if (!inherits(opts, cl))
    stop("Argument 'opts' should be an object of class '", cl, "'")
  return(opts)
}

check.options.one.phase <- function(opts, method) {
  cl  <- paste0(method, ".options.one.phase")
  if (is.null(opts))
    opts <- eval(parse(text = paste0("statMatch.", cl, "()")))
  else if (!inherits(opts, cl))
    stop("Argument 'opts' should be an object of class '", cl, "'")
  return(opts)
}

check.don.weights <- function(df.don, don.weights) {
  if (!is.null(don.weights)) {
    if (!is.vector(don.weights, mode = "numeric") | nrow(df.don) != length(don.weights)) {
      stop("Argument 'don.weights' should be a numeric vector of length ", nrow(df.don), ".")
    } else {
      don.weights <- don.weights
    }
  } else {
    don.weights <- rep(1, nrow(df.don))
  }
  return(don.weights)
}

check.names_data <- function(df.rec, df.don, names.CV, names.NCV, zero.constraints) {
  # Check whether df.rec and df.don are correctly specified as data.frame
  if (!is.data.frame(df.rec)) stop("Input 'df.rec' should be a data.frame")
  if (!is.data.frame(df.don)) stop("Input 'df.don' should be a data.frame")
  # Check whether names.CV exists and it is correctly specified
  if (is.null(names.CV)) names.CV <- intersect(names(df.rec), names(df.don))
  else {
    if (!is.character(names.CV)) stop("Input 'names.CV' should be a character vector.")
    wrong.CV.rec <- !names.CV %in% colnames(df.rec)
    if (any(wrong.CV.rec))
      stop("Variable(s) in 'names.CV' not available in 'df.rec' : ",
           paste(names.CV[wrong.CV.rec], collapse = ", "), ".")
    wrong.CV.don <- !names.CV %in% colnames(df.don)
    if (any(wrong.CV.don))
      stop("Variable(s) in 'names.CV' not available in 'df.don' : ",
           paste(names.CV[wrong.CV.don], collapse = ", "), ".")
  }
  # Check whether names.NCV exists and it is correctly specified
  if (is.null(names.NCV)) names.NCV <- setdiff(names(df.don), names.CV)
  else {
    if (!is.character(names.NCV)) stop("Input 'names.NCV' should be a character vector.")
    wrong.NCV.don <- !names.NCV %in% colnames(df.don)
    if (any(wrong.NCV.don))
      stop("Variable(s) in 'names.NCV' not available in 'df.don' : ",
           paste(names.NCV[wrong.NCV.don], collapse = ", "), ".")
  }
  # Check if df.rec and df.don have common variables
  if (length(intersect(colnames(df.rec), colnames(df.don))) == 0)
    stop("'df.rec' and 'df.don' do not have any common variable.")
  # Check if any variable is specified in both names.CV and names.NCV
  if (length(intersect(names.CV, names.NCV)) > 0)
    stop("The same variable cannot be specified in both names.CV and names.NCV : ",
         paste(intersect(names.CV, names.NCV),collapse = ", "),".")
  # Check whether the common variables have the same class and levels
  for (varname in names.CV) {
    if (class(df.don[, varname, drop = T]) != class(df.rec[, varname, drop = T]))
      stop("Variable ", varname, " has different classes in df.rec and df.don.")
    if (!setequal(levels(df.don[, varname, drop = T]), levels(df.rec[, varname, drop = T])))
      stop("Variable ", varname, " has different levels in df.rec and df.don.")
  }

  out <- list(
    names.CV = names.CV,
    names.NCV = names.NCV
  )

  return(out)
}

check.zero_constraints <- function(df.don, zero.constraints, names.NCV) {
  if (!is.null(zero.constraints)) {
    check <- setdiff(zero.constraints, names.NCV)
    if (length(check) > 0)
      stop("Variable(s) in 'zero.constraints' not in 'names.NCV': ", paste0(check, collapse = ", "))
    check <- colSums(df.don[, zero.constraints, drop = F] == 0) == 0
    if (any(check))
      stop("Problem in 'zero.constraints', no zeros in variable(s): ", paste0(names(check)[check], collapse = ", "))
    df.tmp <- df.don %>%
      select(all_of(zero.constraints)) %>%
      mutate_all(list(~ factor(ifelse(. == 0, 1, 0)))) %>%
      select_all(list(~ paste0("ZC_", .)))
    names.NCV <- c(names.NCV, names(df.tmp))
    df.don <- cbind(df.don, df.tmp)
  }
  return(df.don)
}
