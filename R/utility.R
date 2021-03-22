#' Function from CDVine
#' https://github.com/cran/CDVine
#' @noRd
Frank.itau.JJ <- function(tau) {
  a <- 1
  if (tau < 0) {
    a <- -1
    tau <- -tau
  }
  f <- function(x) {
    x/(exp(x) - 1)
  }
  tauF <- function(x) 1 - 4/x + 4/x^2 * integrate(f, lower = 0 + .Machine$double.eps^0.5, upper = x)$value
  v <- uniroot(function(x) tau - tauF(x), lower = 0, upper = 500, tol = .Machine$double.eps^0.5)$root
  return(a * v)
}

Joe.itau.JJ <- function(tau) {
  if (tau < 0) {
    return(1.000001)
  } else {

    tauF <- function(a) {
      # euler=0.5772156649015328606 1+((-2+2*euler+2*log(2)+digamma(1/a)+digamma(1/2*(2+a)/a)+a)/(-2+a))
      1 + 4/a^2 * integrate(function(x) log(x) * x * (1 - x)^(2 * (1 - a)/a), 0, 1)$value
    }


    v <- uniroot(function(x) tau - tauF(x), lower = 1, upper = 500, tol = .Machine$double.eps^0.5)$root
    return(v)
  }
}
