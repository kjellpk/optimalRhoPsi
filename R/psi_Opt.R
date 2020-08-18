#' @title Optimal Psi
#'
#' @description Evaluate Yohai and Zamar's (1997) Optimal Psi function, its
#'              derivative, antiderivative, rho, and weight functions.
#'
#' @param x a numeric vector of points where the function will be evaluated.
#'
#' @param cc a length 4 numeric vector: \code{cc[1]} is the tuning constant a,
#'           \code{cc[2]} is the lower end point of the support interval \code{lower},
#'           \code{cc[3]} is the upper end point of the support interval \code{upper},
#'           and \code{cc[4]}
#'
#' @details The family of optimal robust \eqn{\psi}{psi} functions is given by
#'
#'          \deqn{\psi_{a}(x) = \mathrm{sign}(x) \max\left(|x| - \frac{a}{\phi(x)}, 0\right)}{%
#'                  psi(x, a) = sign(x)*max(abs(x) - a / dnorm(x), 0)}
#'
#'          where \eqn{\phi(x)}{dnorm(x)} is the standard normal density function and
#'          \eqn{a} is the tuning parameter. For \eqn{x > 0}, \eqn{\psi_{a}(x)}{psi(x, a)}
#'          has support on the interval \eqn{(l_{a}, u_{a})}{(l, u)} where \eqn{l} and
#'          \eqn{u} depend on the value of \eqn{a} (see
#'          \code{\link{OptPsiSupportIntervalFromConst}}). Since the support interval
#'          is fixed for a given \eqn{a}, values for \eqn{l} and \eqn{u} are stored in
#'          the parameter vector \code{cc} to avoid unnecessary recomputation. Lastly,
#'          we introduce a scaling parameter \eqn{c}
#'
#'          \deqn{\mathrm{psi_Opt(x, cc)} = \begin{cases} \mathrm{sign}(x) |x/c| - \frac{a}{\phi(x/c)} & |x/c| \in (l_{a}, u_{a}) \\ 0 & \mathrm{Otherwise}}{%
#'                psi(x, a, c) = sign(x/c)*max(abs(x/c) - a / dnorm(x/c), 0)}
#'
#'
#' @export
psi_Opt <- function(x, cc)
{
  x <- x / cc[4]
  y <- numeric(length(x))
  absx <- abs(x)
  idx <- absx > cc[2] & absx < cc[3]
  xidx <- x[idx]
  #' @importFrom stats dnorm
  y[idx] <- xidx - sign(xidx) * cc[1] / dnorm(xidx)
  y
}


#' @describeIn psi_Opt
#'
#' @export
Psi_Opt <- function(x, cc)
  #' @importFrom pracma erfi
  0.5 * x^2 - cc[1] * pi * Re(erfi(x / sqrt(2)))


#' @describeIn psi_Opt
#'
#' @export
rho_Opt <- function(x, cc)
{
  # Since rho functions are even
  x <- abs(x) / cc[4]

  #Allocate vector to store answer
  rho <- numeric(length(x))

  #Case 1: |x| <= lower: rho is already initialized to zeros

  #Case 2: lower < |x| < upper
  idx <- x > cc[2] & x < cc[3]
  rho[idx] <- cc[4] * (Psi_Opt(x[idx], cc) - Psi_Opt(cc[2], cc))

  #Case 3: |x| >= upper
  idx <- x >= cc[3]
  rho[idx] <- cc[4] * (Psi_Opt(cc[3], cc) - Psi_Opt(cc[2], cc))

  rho
}


#' @describeIn psi_Opt
#'
#' @export
psip_Opt <- function(x, cc)
{
  x <- abs(x) / cc[4]
  y <- numeric(length(x))
  idx <- x > cc[2] & x < cc[3]
  xidx <- x[idx]
  #' @importFrom stats dnorm
  y[idx] <- (1.0 - cc[1] * xidx / dnorm(xidx)) / cc[4]
  y
}


#' @describeIn psi_Opt
#'
#' @export
wgt_Opt <- function(x, cc)
{
  y <- numeric(length(x))
  absx <- abs(x) / cc[4]
  idx <- absx > cc[2] & absx < cc[3]
  absxidx <- absx[idx]
  #' @importFrom stats dnorm
  y[idx] <- (1.0 - cc[1] / (absxidx * dnorm(absxidx))) / cc[4]
  y
}










