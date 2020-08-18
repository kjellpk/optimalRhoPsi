#' @export
psi_modOpt <- function(x, cc = c(0.05988905, 1.328913, 2.34452873, 1.00000000))
{
  psi <- numeric(length(x))
  x <- x / cc[4]
  absx <- abs(x)

  idx <- (absx <= 1.0)
  psi[idx] <- x[idx]

  idx <- (absx > 1.0 & absx < cc[3])
  xidx <- x[idx]
  #' @importFrom stats dnorm
  psi[idx] <- cc[2] * (xidx - sign(xidx) * cc[1] / dnorm(xidx))

  psi
}


#' @export
rho_modOpt <- function(x, cc)
{
  x <- abs(x) / cc[4]
  rho <- numeric(length(x))

  #Case 1: |x| <= 1
  idx <- (x <= 1)
  rho[idx] <- cc[4] * (0.5 * x[idx]^2)

  #Case 2: 1 < |x| < upper
  idx <- (x > 1 & x < cc[3])
  rho[idx] <- cc[4] * (0.5 + cc[2] * (Psi_Opt(x[idx], cc) - cc[5]))

  #Case 3: |x| >= upper
  idx <- (x >= cc[3])
  rho[idx] <- cc[4] * cc[6]

  rho
}


#' @export
wgt_modOpt <- function(x, cc)
{
  wgt <- numeric(length(x))
  x <- abs(x) / cc[4]

  idx <- (x <= 1.0)
  wgt[idx] <- 1.0 / cc[4]

  idx <- (x > 1.0 & x < cc[3])
  xidx <- x[idx]
  #' @importFrom stats dnorm
  wgt[idx] <- cc[2] * (1.0 - cc[1] / (xidx * dnorm(xidx))) / cc[4]

  wgt
}


#' @export
psip_modOpt <- function(x, cc)
{
  x <- abs(x) / cc[4];
  psip <- numeric(length(x))

  psip[x <= 1.0] <- 1.0 / cc[4]

  idx <- x > 1.0 & x < cc[3]
  xidx <- x[idx]
  #' @importFrom stats dnorm
  psip[idx] <- cc[2] * (1.0 - cc[1] * xidx / dnorm(xidx)) / cc[4]

  psip
}

