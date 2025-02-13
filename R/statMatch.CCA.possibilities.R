#' @title Matrix of possibilities for CCA and KCCA
#' @description Do prediction
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @param data_cat_rec Receiver data frame with only categorical variable
#' @param data_cat_don Donor data frame with only categorical variable
#' @param names.CV Name of common categorical variable
#' @param names.CV_Cat_prio Name of common categorical variable that have to matche
#' @param print.details Logical Value, print details or not
#' @importFrom dplyr intersect
#' @noRd

statMatch.CCA.possibilities <- function(data_cat_rec,
                                        data_cat_don,
                                        names.CV = NULL,
                                        names.CV_Cat_prio = NULL,
                                        print.details = FALSE) {
  start.time2 <- Sys.time()
  compte <- 0
  compte2 <- 0
  compteur <- 0
  nbr <- 0
  nbr2 <- 1
  iwh <- 0
  row_don <- nrow(data_cat_don)
  row_rec <- nrow(data_cat_rec)
  if (is.null(names.CV)) {
    names.CV <- intersect(colnames(data_cat_rec), colnames(data_cat_don))
  }
  matrix.tot <- matrix(1L, ncol = row_rec, nrow = row_don)
  for (i in 1:row_rec) {
    iwh <- 1
    df.rec_1 <- as.data.frame(matrix(data = data_cat_rec[i, names.CV], ncol = length(names.CV), nrow = nrow(data_cat_don), byrow = T))
    ligne <- which(rowMeans(data_cat_don[, names.CV] == df.rec_1) == 1)
    # if(!identical(ligne,integer(0)))
    if (!length(ligne) == 0) {
      compte <- compte + 1
      # matrix.tot[ligne,i] <- 1L
      matrix.tot[-ligne, i] <- 0L
    } else {
      compte2 <- compte2 + 1
      while (identical(ligne, integer(0))) {
        ligne <- which(rowSums(data_cat_don[, names.CV] == df.rec_1) == (length(names.CV) - iwh))
        iwh <- iwh + 1
        nbr <- max(nbr, iwh)
      }
      # matrix.tot[ligne,i] <- 1L
      matrix.tot[-ligne, i] <- 0L
    }
    # df.possibilities_i <- data_cat_don[ligne,]
    if (!is.null(names.CV_Cat_prio)) {
      df.rec_vec <- as.numeric(data_cat_rec[i, names.CV_Cat_prio])
      df.rec_2 <- matrix(data = df.rec_vec, ncol = length(names.CV_Cat_prio), nrow = nrow(data_cat_don), byrow = T)
      df.don_2 <- data_cat_don[, names.CV_Cat_prio]
      rs2 <- rowSums(df.rec_2 == df.don_2)
      jnew <- which(rs2 != length(names.CV_Cat_prio))
      matrix.tot[jnew, i] <- 0L
      smt <- sum(matrix.tot[, i])
      while (smt == 0) {
        if (nbr2 == 1) {
          compteur <- compteur + 1
        }
        jnew2 <- which(rowSums(data_cat_don[, names.CV] == df.rec_1) == (length(names.CV) - iwh))
        matrix.tot[jnew2, i] <- 1L
        matrix.tot[jnew, i] <- 0L
        smt <- sum(matrix.tot[, i])
        iwh <- iwh + 1
        nbr2 <- nbr2 + 1
      }
      nbr2 <- 1
    }
  }
  if (print.details) {
    message("Il y a eu ", compte2, " lignes ou les categories identiques n'etaient pas disponibles.")
    message("On a du aller jusqu'a ", nbr, " categories utilisees en moins")
    time2 <- Sys.time() - start.time2
    message("Time: ", time2, " min")
    if (!is.null(names.CV_Cat_prio)) {
      message("Il y a eu ", compteur, " lignes ou les categorie les plus proches n'ont pas ete utilise car les categories prioritaires devaient etre utilisees.")
    }
    message("Nombre de colonnes egales a 0 : ", sum(colSums(matrix.tot) == 0))
    message("Nombre de lignes egales a 0 : ", sum(rowSums(matrix.tot) == 0))

    message("Nombre de colonnes egales a 1 : ", sum(colMeans(matrix.tot) == 1))
    message("Nombre de lignes egales a 1 : ", sum(rowMeans(matrix.tot) == 1))
  }
  return(matrix.tot)
}
