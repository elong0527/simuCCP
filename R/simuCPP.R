#' Simulate Concordance Probability based on bivariate copula with known family and therotical value
#'
#' This function generate multivariate uniform distribution with designed copula family and
#' therotical Concordance probability (C-index, iAUC, AUC)
#'
#'
#' @export
#' @param N sample size
#' @param metric simulate metric (C-index, iAUC, AUC)
#' @param value of therotical value of (C-index, iAUC, AUC)
#' @param family d - 1 dimension of copula family, where d is number of covariates check \code{BiCopPar2Tau} for detail.
#' @param c0 prevalence; cutoff point of latent uniform parameter (Used in AUC)
#' @param fam2 higher order copula family ( dimension: (d^2 - 3d + 2)/2 ) (default: conditional independent)
#' @param par2 parameters of higher order copula family ( dimension: (d^2 - 3d + 2)/2 ) (default: all 0)
#' @examples
#' # C-index
#' U1 = simuCCP(300, metric = "Cind", value = 0.7, family = c(1,3))$data
#' ( cor(U1, method = "kendall") + 1 ) / 2
#'
#' # C-index with a dependent copula
#' U1 = simuCCP(300, metric = "Cind", value = 0.7, family = c(1,3), fam2 = 3, par2 = 1)$data
#' ( cor(U1, method = "kendall") + 1 ) / 2
#'
#' # AUC
#' U1 = simuCCP(300, metric = "AUC", value = 0.7, family = c(1,3))$data
#' library(pROC)
#' apply(U1[,-1], 2, function(x) auc(U1[,1], x) )
#'
#' # AUC with a dependent copula
#' U1 = simuCCP(300, metric = "AUC", value = 0.7, family = c(1,3), fam2 = 3, par2 = 1)$data
#' library(pROC)
#' apply(U1[,-1], 2, function(x) auc(U1[,1], x) )
simuCCP <- function(N, metric = "Cind", value = 0.7, family = c(1,3,4), c0 = 0.5, fam2 = NULL, par2 = NULL){
  if(! metric %in% c("Cind","iAUC","AUC"))
    stop("We only support C-index, iAUC and AUC")

  d <- length(family) + 1
  dd <- d * (d-1) / 2

  par1 <- vector()
  for(i in 1:(d-1) ){
    if(metric == "Cind")  par1[i] <- BiCopCind2Par(family[i], value)
    if(metric == "AUC" )  par1[i] <- BiCopAUC2Par(family[i], value, c0)
    if(metric == "iAUC" ) par1[i] <- BiCopiAUC2Par(family[i], value)
  }

  if(any(is.null(fam2), is.null(par2))){
    fam1 <- c(family, rep(0, dd - d + 1))
    par1 <- c(par1, rep(0, dd - d + 1))
  }else{
    fam1 <- c(family, fam2)
    par1 <- c(par1, par2)
  }
  U1 = CDVineSim(N, fam1, par1, type = 1) # C-Vine

  if(metric == "AUC") U1[,1] <- U1[,1] > c0
  list( data = U1, family = fam1, parameter = par1)
}

#' Simulate Concordance Probability based on bivariate copula with known family and parameters
#'
#' This function generate multivariate uniform distribution with designed copula family and
#' parameters (C-index, iAUC, AUC)
#'
#' @export
#' @param N sample size
#' @param metric simulate metric (C-index, iAUC, AUC)
#' @param family A d*(d-1)/2 integer vector of C-vine pair-copula families (detail is in \code{CDVineSim})
#' @param par	A d*(d-1)/2 vector of pair-copula parameters.
#' @param c0 prevalence; cutoff point of latent uniform parameter (Used in AUC)
simuCCP_para <- function(N, metric, family, par, c0 = 0.5){
  if(! metric %in% c("Cind","iAUC","AUC"))
    stop("We only support C-index, iAUC and AUC")
  U1 = CDVineSim(N, family, par, type = 1) # C-Vine
  if(metric == "AUC") U1[,1] <- U1[,1] > c0

  U1
}


