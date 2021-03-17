context("BiCopPar2AUC")
library(simuCCP)
library(numDeriv)
library(pracma)
library(CDVine)

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
BiCopPar2AUC1 <- function(family, par, par2 = 0, c0, type = "cumulative"){
  #   auc1 <- function(u,v, family, par, par2 = 0, c0){
  #     v0 <- rep(c0, length(u))
  #     (v > c0) * BiCopCDF(u, v0, family = family, par = par, par2 = par2) *
  #                BiCopPDF(u,v, family = family, par = par, par2 = par2)
  #   }
  if(! type %in% c("cumulative", "incidence"))
    stop("type should be cumulative or incidence")

  auc1 <- function(u,v, family, par, par2 = 0, c0){
    v0 <- rep(c0, length(u))
    BiCopCDF(u, v0, family = family, par = par, par2 = par2) *
      BiCopPDF(u,v, family = family, par = par, par2 = par2)
  }

  auc2 <- function(u,vv, family, par, par2 = 0){
    v0 <- rep(vv, length(u))
    BiCopCDF(u, v0, family = family, par = par, par2 = par2) *
      BiCopPDF(u,v0, family = family, par = par, par2 = par2)
  }

  auc3 <- function(u,vv, family, par, par2 = 0, c0){
    foo <- function(vv, uu) BiCopCDF(uu, vv, family = family, par = par, par2 = par2)
    term1 = sapply(u, function(u) grad(foo, x = c0, uu = u) )
    term1 * BiCopPDF(u,vv, family = family, par = par, par2 = par2)
  }

  dd1 <- function(u,v, family, par, par2 = 0){
    BiCopCDF(u, v, family = family, par = par, par2 = par2) *
      BiCopPDF(u,v, family = family, par = par, par2 = par2)
  }

  dd2 <- function(u,v, family, par, par2 = 0, c0 = c0){
    # print(v)
    # print(c0)
    v0 <- pmin(v, c0)
    BiCopCDF(u, v0, family = family, par = par, par2 = par2) *
      BiCopPDF(u,v, family = family, par = par, par2 = par2)
  }


  if(type == "cumulative"){
    pp <- integral2(auc1, xmin = 0, xmax = 1, ymin = c0, ymax = 1, family = family, par = par, par2 = par2, c0 = c0)
    auc = pp$Q / c0 / (1 - c0)
  }
  if(type == "incidence"){

    ## term 1 in the AUC equation (cumulative AUC derivative)
    pp <- function(c0){
      integral2(auc1, xmin = 0, xmax = 1, ymin = c0, ymax = 1, family = family, par = par, par2 = par2, c0 = c0)$Q
    }
    auc00 = numDeriv::grad(pp, x = c0) / (1 - c0)

    ## term 2 in the AUC equation (difference term)
    tmp = pracma::integral(auc2, xmin = 0, xmax = 1, vv = c0, family = family, par = par, par2 = par2)
    auc01 = tmp / (1 - c0)

    ## term 3 in the AUC equation (incidence AUC)
    pp1 <- integral2(auc3, xmin = 0, xmax = 1, ymin = c0, ymax = 1, family = family, par = par, par2 = par2, c0 = c0)
    auc02 = pp1$Q / (1 - c0)

    dd00 = pp(c0)
    dd01 = integral2(dd1, xmin = 0, xmax = 1, ymin = 0, ymax = c0, family = family, par = par, par2 = par2)$Q
    dd02 = integral2(dd2, xmin = 0, xmax = 1, ymin = 0, ymax = 1, family = family, par = par, par2 = par2, c0 = c0)$Q
    # print("this is dd")
    auc = c(auc00, auc01, auc02)
    # auc = c(dd00, dd01, dd02)
    # auc = auc02
    # pp(c0)
  }
  auc
}

test_that("Transfer Copula Parmeter and time dependent AUC", {
  family = 1
  c0 = 0.5
  par = 0.5
  expect_equal(round(BiCopPar2AUC(family = family, par = par, c0 = c0),3),
               round(BiCopPar2AUC1(family = family, par = par, c0 = c0),3) )

  expect_equal(round(BiCopPar2AUC(family = family, par = par, c0 = c0, type = "incidence"),3),
               round(BiCopPar2AUC1(family = family, par = par, c0 = c0, type = "incidence")[3],3) )
})



