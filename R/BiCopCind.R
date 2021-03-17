#' C-index value of a bivariate copula
#'
#' This function computes the theoretical C-index value of a bivariate copula for given parameter values.
#'
#' @export
#' @param family copula family. check \code{BiCopPar2Tau} for detail.
#' @param par Copula parameter.
#' @param par2 Second parameter for the two parameter BB1, BB6, BB7 and BB8 copulas (default: \code{par2 = 0}).
#' @examples
#' # BiCopPar2AUC(3, 1, c0 = 0.5)
BiCopPar2Cind <- function(family, par, par2 = 0){
  tau <- CDVine::BiCopPar2Tau(family, par, par2)
  Cind <- (tau + 1)/2
  Cind
}

#' Parameter of a bivariate copula for a given C-index value
#'
#' This function computes the parameter of a one parameter bivariate copula for a given value of C-index.
#'
#' @export
#' @param family copula family. check \code{BiCopPar2Tau} for detail.
#' @param Cind value (numeric in [0.5, 1])
#' @examples
#' # BiCopCstat2Par(3, 0.8)
BiCopCind2Par <- function(family, Cind){
  tau <- 2 * Cind - 1
  CDVine::BiCopTau2Par(family, tau)
}


