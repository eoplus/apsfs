#' US76 Standard Atmosphere
#'
#' A dataset containing model profiles of temperature, pressure, and density of 
#' air from -1 km to 1000 km.
#'
#' @format A data frame with 603 rows and 4 variables:
#' \describe{
#'   \item{Z}{Geometric height from surface (km)}
#'   \item{Temperature}{Temperature (degree K)}
#'   \item{Pressure}{Pressure (mbar)}
#'   \item{Density}{Air density (kg / m3)}
#' }
#' @source Generated with the Public Domain Aeronautical Software (PDAS) code 
#' atmos76.f90 \url{http://www.pdas.com/atmosdownload.html}

"us76"

#' 6SV Standard Aerosol Phase Functions
#'
#' A dataset containing the scattering phase functions of three standard 
#' aerosol models included in 6SV (maritime, continental and urban).
#'
#' @format A data frame with 83 rows and 21 variables, the first being the 
#' angular reference (degrees) and the phase function in 20 wavelengths.
#'
#' @source 6SV code \url{http://6s.ltdri.org/}

"continental_ph_6sv"
"maritime_ph_6sv"
"urban_ph_6sv"

#' 6SV Standard Aerosol Coefficients
#'
#' A dataset containing the optical coefficients of three standard 
#' aerosol models included in 6SV (maritime, continental and urban).
#'
#' @format A data frame with 20 rows and 7 variables, the first being the 
#' wavelength reference (nm) and optical coefficients:
#' \describe{
#'   \item{Nor_Ext_Co}{Normalized extinction coefficients}
#'   \item{Nor_Sca_Co}{Normalized scattering coefficients}
#'   ...
#' }
#'
#' @source 6SV code \url{http://6s.ltdri.org/}

"continental_coef_6sv"
"maritime_coef_6sv"
"urban_coef_6sv"

#' 6SV Rayleigh Optical Thickness Coefficients at surface pressure
#'
#' A dataset containing the optical thickness of Rayleigh scattering at selected
#' wavelengths of the 6SV aerosol data grid. The surface pressure is 1013.25 
#' mbar.
#'
#' @format A data frame with 20 rows and 2 variables, the first being the 
#' wavelength reference (nm) and the second the optical thickness.
#'
#' @source 6SV code \url{http://6s.ltdri.org/}

"tau_ray_6sv"

#' Rayleigh Cumulative Distribution Function
#'
#' A dataset containing the Rayleigh CDF at 0.18 degree resolution.
#'
#' @format A data frame with 1000 rows and 2 variables, the first being the 
#' angular reference (radians) and the second the CDF.
#'
#' @source 6SV code \url{http://6s.ltdri.org/}

"rayleigh_cdf"

