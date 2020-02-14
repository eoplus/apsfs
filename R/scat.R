
#' Aerosol CDF Look Up Table
#'
#' Calculates a Look Up Table (LUT) for the aersol cumulative distribution 
#' function to be used with the Monte Carlo calculations.
#'
#' @param ph   Phase function table at single wavelength.
#' @param xout Angle points for calculation. 
#' @param type One of "fixed", "dynamic". See Details.
#'
#' @details This function calculates the cumulative density function of aerosol 
#' phase functions. It can use a fixed integration grid based on xout (must be 
#' very fine resolution!) or an adaptive integration (any xout). Note that the 
#' convention in atmospheric optics is to have the PF not normalized by 4pi 
#' steradian, so it is unitless instead of 1/sr as in ocean optics. The integral 
#' to 1 is:
#' integral_0^2pi integral_0^pi (PF(theta) / 4pi) * sin(theta) * dtheta * dphi = 1
#' such that:
#' integral_0^pi PF(theta) / 2 * sin(theta) * dtheta = 1
#'
#' Due to small numerical errors, the phase functions are renormalized to 1 
#' before returning.
#'
#' @examples
#' psi_out <- c(0, 10^seq(log10(0.00001), log10(180), length.out = 1000))
#' calc_cdf(continental_ph_6sv[, c(1, 8)], psi_out)
#'
#' @export

calc_cdf <- function(ph, xout, type = "fixed") {

  if(type == "fixed") {

    ph <- approx(x = ph[, 1], y = log10(ph[, 2] / 4 / pi), xout = xout) %>%
          as.data.frame()
    ph[, 1] <- ph[, 1] * pi / 180
    ph[, 2] <- 10^ph[, 2]
    tp <- 2 * pi * ph[, 2] * sin(ph[, 1]) 
    tp <- (tp[-1] + tp[-length(tp)]) / 2 # trapezoidal rule
    tp <- diff(ph[, 1]) * tp
    ph[, 3] <- c(0, cumsum(tp))
    ph[, 2] <- ph[, 2] / ph[nrow(ph), 3]
    ph[, 3] <- ph[, 3] / ph[nrow(ph), 3]

    return(ph[, c(1, 3)])

  } else {

    ph[, 1] <- ph[, 1] * pi / 180
    xout    <- xout * pi / 180 

    f <- function(x) {

      10^approx(x = ph[, 1], y = log10(ph[, 2] / 4 / pi), xout = x)$y * sin(x)

    }

    cdf <- numeric(length(xout))

    for(i in 2:length(xout)) {

      cdf[i] <-  2 * pi * integrate(f, lower = 0, upper = xout[i])$value

    }

    cdf <- cdf / cdf[length(xout)]

    return(cbind(xout, cdf))

  }

}

#' Rayleigh scattering phase function
#'
#' @param psi Scattering angle (radians), [0,pi].
#'
#' @export

ray_ph  <- function(psi) {
  f <- (1 - 0.04) / (1 + 0.04)
  (1.0 + f * cos(psi)^2) * sin(psi) * 3 / (3 + f)
}


#' Ancillary scattering functions required in mc.R:
#' @export

scat_aer <- function(ru, aer_cdf) {

  id <- findInterval(x = ru, vec = aer_cdf[, 2], left.open = TRUE, all.inside = TRUE)
  aer_cdf[id, 1]

}

scat_ray <- function(ru) {

  id <- findInterval(x = ru, vec = rayleigh_cdf[, 2], left.open = TRUE, all.inside = TRUE)
  rayleigh_cdf[id, 1]

}

#' Calculate Rayleigh optical depth
#'
#' Calculates the Rayleigh optical depth depending on pressure and height above 
#' sea level.
#'
#' @param atm    Atmosphere profile, see Details.
#' @param lambda Wavelength (nm), [200,1100].
#' @param co2    CO2 concentration (ppm), [0,1E6].
#' @param lat    Geodetic latitude (degrees), [-90,90].
#'
#' @details The atmospheric profile must be a data frame with at least two 
#' columns: Z, with height above sea level (km); and Pressure, with atmospheric
#' pressure (mbar).
#'
#' Note that all arguments should be length = 1 for atm with multiple rows 
#' (expected). For atm with a single row (e.g., sea level), one of the arguments
#' can have length > 1.
#'
#' @return
#' A vector with the same length as the number of rows of atm with the Rayleigh
#' optical depth.
#'
#' @references
#' Bodhaine, B.A.; Wood, N.B.; Dutton, E.G.; Slusser, J.R. 1999. On Rayleigh 
#' optical depth calculations. Journal of Atmospheric and Oceanic Technology 16, 
#' 11, 1854-1861. DOI: 10.1175/1520-0426(1999)016<1854:ORODC>2.0.CO;2
#'
#' @example
#' rayleigh_od(us62, lambda = 400, co2 = 400, lat = 45)
#' rayleigh_od(data.frame(Z = 0, Pressure = 1013.25), lambda = 400:700, 
#'             co2 = 400, lat = 45)
#'
#' @export

rayleigh_od <- function(atm, lambda, co2 = 400, lat = 45) {

 lambda_mum <- lambda * 1E-3
 lambda_cm  <- lambda * 1E-7

 Z <- atm$Z * 1000       # km to m
 P <- atm$Pressure * 1E3 # mbar == 1E3 * dyn / cm^2

 # Mean molecular weight of air (gm / mol):
 co2 <- co2 / 1E6
 ma  <- 15.0556 * co2 + 28.9595

 # Mean refractive index of dry air (at 288.15 K and 101.325 Pa):
 na_450 <- rho::n_air(lambda)
 na     <- 1 + (na_450 - 1) * (1 + 0.54 * (co2 - 0.00045))

 # Scattering cross-section of Reyleigh scattering:
 Ns        <- 2.546899E19 # molecules per cm3 (at 288.15 K and 101.325 Pa)
 depolr_n2 <- 1.034 + 3.17E-4 / lambda_mum^2
 depolr_o2 <- 1.096 + 1.385E-3 / lambda_mum^2 + 1.448E-4 / lambda_mum^4
 depolr_a  <- (78.084 * depolr_n2 + 20.946 * depolr_o2 + 0.934 * 1 + 
              co2 * 100 * 1.15) / (78.084 + 20.946 + 0.934 + co2 * 100)
 sigma     <- depolr_a * (24 * pi^3 * (na^2 - 1)^2) / 
              (lambda_cm^4 * Ns^2 * (na^2 + 2)^2)

 # Gravity (cm / s^2):
 lat <- lat * pi / 180
 g0  <- 980.6160 * (1 - 0.0026373 * cos(2 * lat) + 0.0000059 * cos(2 * lat)^2)
 g   <- g0 - (3.085462E-4 + 2.27E-7 * cos(2 * lat)) * Z +
        (7.254E-11 + 1E-13 * cos(2 * lat)) * Z^2 -
        (1.517E-17 + 6E-20 * cos(2 * lat)) * Z^3

 # Rayleigh optical depth:
 An      <- 6.0221367E23 #6.0221409E23
 tau_ray <- sigma * P * An / ma / g

 return(tau_ray)

}
