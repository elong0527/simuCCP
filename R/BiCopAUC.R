#' AUC value of a bivariate copula
#'
#' This function computes the theoretical AUC value of a bivariate copula for given parameter values.
#'
#' @export
#' @param family copula family. check \code{BiCopPar2Tau} for detail.
#' @param par Copula parameter.
#' @param par2 Second parameter for the two parameter BB1, BB6, BB7 and BB8 copulas (default: \code{par2 = 0}).
#' @param c0 prevalence; cutoff point of latent uniform parameter
#' @param type type is cumulative or incidence
#' @examples
#' # BiCopPar2AUC(3, 1, c0 = 0.5)
BiCopPar2AUC <- function(family, par, par2 = 0, c0, type = "cumulative"){

  if(! type %in% c("cumulative", "incidence"))
    stop("type should be cumulative or incidence")

  if(type == "cumulative"){
    auc1 <- function(u, family, par, par2 = 0, c0){
      v0 <- rep(c0, length(u))
      BiCopCDF(u, v0, family = family, par = par, par2 = par2)
    }

    term1 <- integral(auc1, 0, 1, family = family, par = par, par2 = par2, c0 = c0)
    auc <- (term1 - c0^2/2)/ c0 / (1 - c0)
  }

  if(type == "incidence"){
    auc5 <- function(u, family, par, par2 = 0, c0){
      v0 <- rep(c0, length(u))
      h <- BiCopHfunc(u, v0, family = family, par = par, par2 = par2)
      h$hfunc2 * (1 - h$hfunc1)
    }
    eta <- 0
    term5 <- pracma::integral(auc5, 0 + eta, 1 - eta, family = family, par = par, c0 = c0)

    auc <- term5 / (1 - c0)

  }

  auc

}




#' Parameter of a bivariate copula for a given Kendall's tau value
#'
#' This function computes the parameter of a one parameter bivariate copula for a given value of AUC.
#'
#' @export
#' @param family copula family. check \code{BiCopPar2Tau} for detail.
#' @param AUC the AUC value
#' @param c0 prevalence; cutoff point of latent uniform parameter
#' @examples
#' # BiCopAUC2Par(3, 0.7274113, c0 = 0.5, lower = 1e-5, upper = 5)
BiCopAUC2Par <- function(family, AUC, c0){
  foo <- function(par){
    BiCopPar2AUC(family = family, par = par, par2 = 0, c0 = c0) - AUC
  }
  # uniroot(foo, lower = lower, upper = upper)$root
  f <- try( multiroot(foo, start = -0.5)$root, silent = T)
  if(class(f) == "try-error") f <- try( multiroot(foo, start = 0.5)$root, silent = T)
  if(class(f) == "try-error") f <- try( multiroot(foo, start = 1.5)$root, silent = T)
  f
}




