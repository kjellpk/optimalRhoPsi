#' @export
OptPsiSupportIntervalFromConst <- function(a, tol = 1e-8)
{
  a <- as.numeric(a[1])

  if(a > 0.0 && a < 0.2419707) {
    #' @importFrom stats dnorm
    g <- function(x, a) x*dnorm(x) - a
    #' @importFrom stats uniroot
    lower <- uniroot(f = g, interval = c(0.0, 1.0), a = a, check.conv = TRUE, tol = tol)$root
    upper <- uniroot(f = g, interval = c(1.0, 5.0), a = a, extendInt = "downX", check.conv = TRUE, tol = tol)$root
    c(lower, upper)
  }
  else if(a == 0.0) {
    c(0.0, Inf)
  }
  else {
    c(NA_real_, NA_real_)
  }
}


#' @export
computeEfficiencyFromConst_Opt <- function(a, v = 1.0)
{
  support <- OptPsiSupportIntervalFromConst(a)
  cc <- c(a, support, v)

  #' @importFrom stats dnorm
  integrand.top <- function(x, cc)
    psi_Opt(x, cc)^2 * dnorm(x)

  #' @importFrom stats integrate
  nu.top <- 2.0 * integrate(integrand.top, support[1], support[2], cc = cc)$value

  #' @importFrom stats dnorm
  integrand.bottom <- function(x, cc)
    psip_Opt(x, cc) * dnorm(x)

  #' @importFrom stats integrate
  nu.bottom <- (2.0 * integrate(integrand.bottom, support[1], support[2], cc = cc)$value)^2

  nu <- nu.top / nu.bottom
  1.0 / nu
}


#' @export
computeConstFromEfficiency_Opt <- function(eff, tol = 1e-8)
{
  obj <- function(a, e)
    e - computeEfficiencyFromConst_Opt(a)

  #' @importFrom stats uniroot
  uniroot(obj, interval = c(tol, 0.2419706), e = eff, check.conv = TRUE, tol = tol)$root
}


#' @export
computeTuningPsi_Opt <- function(efficiency = 0.95, a = NULL)
{
  if(is.null(a))
    a <- computeConstFromEfficiency_Opt(efficiency)
  support <- OptPsiSupportIntervalFromConst(a)
  cc <- c(a, support, 1.0)
  cc[5] <- Psi_Opt(cc[2], cc)
  cc[6] <- Psi_Opt(cc[3], cc) - cc[5]
  names(cc) <- c("a", "lower", "upper", "c", "Psi_Opt(lower)", "rho(Inf)")
  cc
}

#' @export
computeBreakdownPoint_Opt <- function(cc)
{
  f <- function(x, cc)
    rho_Opt(x, cc = cc) * dnorm(x)

  #' @importFrom stats integrate
  l <- integrate(f, lower = -cc[3], upper = -cc[2], cc = cc)$value
  r <- integrate(f, lower = cc[2], upper = cc[3], cc = cc)$value
  (l + r) / rho_Opt(cc[3], cc = cc)
}


#' @export
computeTuningChi_Opt <- function(cc)
{
  g <- function(v, cc) {
    cc[4] <- v
    computeBreakdownPoint_Opt(cc) - 0.5
  }

  #' @importFrom stats uniroot
  cc[4] <- uniroot(g, c(0.1, 1.0), cc = cc, tol = 1e-8)$root
  cc
}



