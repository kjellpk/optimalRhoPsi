#' @export
computeEfficiencyFromConst_modOpt <- function(a, v = 1.0)
{
  cc <- c(a, dnorm(1.0) / (dnorm(1.0) - a), OptPsiSupportIntervalFromConst(a)[2], v)

  #' @importFrom stats dnorm
  integrand.top <- function(x, cc)
    psi_modOpt(x, cc)^2 * dnorm(x)

  #' @importFrom stats integrate
  nu.top <- 2.0 * integrate(integrand.top, 0.0, cc[3], cc = cc)$value

  #' @importFrom stats dnorm
  integrand.bottom <- function(x, cc)
    psip_modOpt(x, cc) * dnorm(x)

  #' @importFrom stats integrate
  nu.bottom <- (2.0 * integrate(integrand.bottom, 0.0, cc[3], cc = cc)$value)^2

  nu <- nu.top / nu.bottom
  1.0 / nu
}


#' @export
computeConstFromEfficiency_modOpt <- function(eff, interval = c(1e-6, 0.2255))
{
  obj <- function(a, e)
    e - computeEfficiencyFromConst_modOpt(a)

  #' @importFrom stats uniroot
  uniroot(obj, interval = interval, e = eff, check.conv = TRUE, tol = 1e-8)$root
}


#' @export
computeTuningPsi_modOpt <- function(efficiency = 0.95, a = NULL)
{
  if(is.null(a))
    a <- computeConstFromEfficiency_modOpt(efficiency)
  support <- OptPsiSupportIntervalFromConst(a)
  #' @importFrom stats dnorm
  cc <- c(a, dnorm(1) / (dnorm(1) - a), support[2], 1.0)
  cc[5] <- Psi_Opt(1.0, cc)
  cc[6] <- (0.5 + cc[2] * (Psi_Opt(cc[3], cc) - cc[5]))
  names(cc) <- c("a", "normConst", "upper", "c", "Psi_Opt(1)", "rho(Inf)")
  cc
}


#' @export
computeBreakdownPoint_modOpt <- function(cc)
{
  f <- function(x, cc)
    rho_modOpt(x, cc = cc) * dnorm(x)

  #' @importFrom stats integrate
  integrate(f, lower = -cc[3], upper = cc[3], cc = cc)$value / (cc[[4]] * cc[[6]])
}


#' @export
computeTuningChi_modOpt <- function(cc)
{
  g <- function(v, cc) {
    cc[4] <- v
    computeBreakdownPoint_modOpt(cc) - 0.5
  }

  #' @importFrom stats uniroot
  cc[4] <- uniroot(g, c(0.1, 1.0), cc = cc, tol = 1e-8)$root
  cc
}


