#' iAUC value of a bivariate copula
#'
#' This function computes the theoretical AUC value of a bivariate copula for given parameter values.
#'
#' @export
#' @param family copula family. check \code{BiCopPar2Tau} for detail.
#' @param par Copula parameter.
#' @param par2 Second parameter for the two parameter BB1, BB6, BB7 and BB8 copulas (default: \code{par2 = 0}).
#' @param show_AUC display time depednent AUC (default: TRUE)
#' @param c0_length length of cutoff points
#' @param weight weight of time dependent AUC
#' @param q_density marginal quantile function of Time distribution
#' @param d_density marginal PDF of Time distribution
#' @examples
#' BiCopPar2iAUC( family = 3, par = 1)
#' ## Note: the results are not stable
BiCopPar2iAUC <- function(family, par, par2 = 0, show_AUC = T, c0_length = 50, weight = NULL, q_density = NULL, d_density = NULL, type = "cumulative"){

  if(show_AUC){
    c0 <- seq(1e-4,1 - 1e-4, length = c0_length)
    AUC <- sapply(c0, BiCopPar2AUC, family = family, par = par, par2 = par2, type = type)
  }else{
    c0 = NULL
    AUC = NULL
  }

  if( any(is.null(weight), is.null(density) ) ){
    # iAUC = mean(AUC)
    if(type == "cumulative"){
      foo <- function(x){
        sapply(x, BiCopPar2AUC, family = family, par = par, par2 = par2)
      }
      iAUC =  pracma::integral(foo, xmin = 0, xmax = 1, vectorized = F)
      # iAUC = 0
    }
    if(type == "incidence"){
      iAUC = BiCopPar2Cind(family, par)
    }
  }else{
    print("double check how to use the weight in the paper")
  }
  list(AUC = AUC, c0 = c0,  iAUC = iAUC)
}

#' Parameter of a bivariate copula for a given Kendall's tau value
#'
#' This function computes the parameter of a one parameter bivariate copula for a given value of AUC.
#'
#' @export
#' @param family copula family. check \code{BiCopPar2Tau} for detail.
#' @param iAUC the integrated AUC value
#' @param c0 prevalence; cutoff point of latent uniform parameter
#' @examples
#' BiCopiAUC2Par(3, 0.75)
BiCopiAUC2Par <- function(family, iAUC){
  foo <- function(par){
    BiCopPar2iAUC(family = family, par = par, par2 = 0, show_AUC = F)$iAUC - iAUC
  }
  # uniroot(foo, lower = lower, upper = upper)$root
  f <- try( multiroot(foo, start = -0.5)$root, silent = T)
  if(class(f) == "try-error") f <- try( multiroot(foo, start = 0.5)$root, silent = T)
  if(class(f) == "try-error") f <- try( multiroot(foo, start = 1.5)$root, silent = T)
  f
}




