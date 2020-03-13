
#' Aerosol CDF Look Up Table
#'
#' Calculates a Look Up Table (LUT) for the aersol cumulative distribution 
#' function to be used with the Monte Carlo calculations.
#'
#' @param ph     Phase function table at single wavelength.
#' @param xout   Angle points for calculation.  See Details.
#' @param method One of "fixed", "adaptative". See Details.
#'
#' @details This function calculates the cumulative distribution function of 
#' aerosol phase functions between 0 and pi. It can use a fixed integration grid 
#' based on xout (must be fine angular grid resolution!) or an adaptive 
#' integration (any xout) based on \code{integrate}. Regardless of the method 
#' used, a fine resolution is recommended to reduce bias on the PSF simulation. 
#' If xout is not specified, the defaault is to use 1000 log10 spaced values 
#' between 0 and pi. If xout is provided, it must include 0 and pi.
#'
#' For the integration, a linear interpolation (\code{approx}) is used in log10
#' space for the aerosol phase function. Method 'fixed' uses fixed integration 
#' steps with quadrature rule.
#'
#' Note that the convention in atmospheric optics is to have the PF not 
#' normalized by 4pi steradian, so it is unitless instead of 1/sr as in ocean 
#' optics. The integral to 1 is:
#'
#' integral_0^2pi integral_0^pi (PF(theta, phi) / 4pi) * sin(theta) * dtheta * dphi = 1,
#'
#' which in isotropic conditions (lack of azimuthal dependence) simplified to:
#'
#' integral_0^pi (PF(theta) / 2) * sin(theta) * dtheta = 1
#'
#' Due to small numerical errors, the phase functions are renormalized to 1 
#' before returning.
#'
#' @return A numeric vector with the aerosol scattering CDF.
#'
#' @examples
#' cdf_lut(continental_ph_6sv[, c(1, 8)])
#'
#' @export

cdf_lut <- function(ph, xout, method = "fixed") {

  if(missing(xout)) {
    xout <- c(0, 10^seq(log10(1E-5), log10(180), length.out = 1E3)) * pi / 180
  } else {
    if(min(xout) != 0 | max(xout) != pi)
      stop("CDF limits should include 0 and pi")
    if(length(xout) < 1E3 & method == "fixed")
      paste("Fixed integration uncertainties will be larger for lower angular",
      "resolution. 1000 log10 spaced intervals is recommended.") %>%
      warning(call. = FALSE)
  }

  if(method == "fixed") {

    ph <- approx(x = ph[, 1], y = log10(ph[, 2] / 4 / pi), xout = xout) %>%
          as.data.frame()
    ph[, 2] <- 10^ph[, 2]
    tmp     <- 2 * pi * ph[, 2] * sin(ph[, 1]) 
    tmp     <- (tmp[-1] + tmp[-length(tmp)]) / 2 # trapezoidal rule
    tmp     <- diff(ph[, 1]) * tmp
    ph$cdf  <- c(0, cumsum(tmp))

    # Renormalize to remove effects of numerical approximations:
    ph[, 2] <- ph[, 2] / ph[nrow(ph), 3]
    ph[, 3] <- ph[, 3] / ph[nrow(ph), 3]

    return(ph[, c(1, 3)])

  } else {

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
 na_450 <- n_air(lambda)
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


#' Dispersion formula for the refractive index of air
#'
#' This function calculates the real part of the refractive index of standard 
#' dry air relative to vaccum. See details.
#'
#' @param lambda Wavelength in vacuum (nm).
#'
#' @details The dispersion model implemented is that by Ciddor (1996), for the 
#' refractive index of standard dry air. It is reference for dry air (0 
#' moisture), 15 ºC, 101.325 Pa and 450 ppm CO2. This equation is recomended by 
#' Zhang et al. (2009) for the convertion of (saline) water refractive index 
#' relative to air to that relative to vacuum.
#'
#' Ciddor (1996) also provides equations to account for different levels of CO2 
#' and water vapor. Those are not implemented in the current version.
#'
#' The valid range is 200 nm to 1100 nm. If the refractive index is requested 
#' outside this range, the function will return a warning.
#'
#' @return A numeric vector with the real part of refractive index of standard 
#' dry air (unitless).
#'
#' @references
#' Ciddor, P. E. 1996. Refractive index of air: new equations for the visible 
#' and near infrared. Applied Optics 35, 9, 1566-73. DOI: 10.1364/AO.35.001566
#'
#' Zhang, X.; Hu, L.; He, M.-X. 2009. Scattering by pure seawater: effect of 
#' salinity. Optics Express 17, 7, 5698-710. DOI: 10.1364/OE.17.005698
#'
#' @examples
#' # Get the refractive index for teh visible range: 
#' n_air(400:700) 
#'
#' @export

n_air <- function(lambda) {

  if(any(lambda < 200) || any(lambda > 1100))
    warning("Requested wavelength outside model domain: [200,1100]")

  l2  <- (lambda * 1E-3)^-2  # Convert to wavlength in microns and square

  k0  <- 5792105
  k1  <- 238.0185
  k2  <- 167917
  k3  <- 57.362

  n <- 1 + ((k0 / (k1 - l2) + k2 / (k3 - l2)) * 1E-8)

  return(n)

}

