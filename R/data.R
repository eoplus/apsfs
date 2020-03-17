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

#' @rdname continental_ph_6sv
 
"maritime_ph_6sv"

#' @rdname continental_ph_6sv 

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

#' @rdname continental_coef_6sv 

"maritime_coef_6sv"

#' @rdname continental_coef_6sv

"urban_coef_6sv"

#' Rayleigh Cumulative Distribution Function
#'
#' A dataset containing the Rayleigh CDF at 0.18 degree resolution.
#'
#' @format A data frame with 1000 rows and 2 variables, the first being the 
#' angular reference (radians) and the second the CDF.
#'
#' @source 6SV code \url{http://6s.ltdri.org/}

"rayleigh_cdf"

#' Annular simulations
#'
#' Contains the atmospheric point spread functions (APSFs) for nadir vieweing 
#' sensors, simulated for the different standard aerosol models of the 6SV 
#' radiative transfer code: maritime (mar), continental (con) and urban (urb). 
#' Simulations were carried out without Rayleigh contribution, and with an 
#' aerosol scale height of 2 km. The APSFs for Rayleigh alone (ray) are also 
#' available at four pressure levels: 1100 mbar (p1100), 1013.25 mbar (p1013), 
#' 750 mbar (p0750) and 500 mbar (p0500). Sensor was positioned at TOA, with 
#' nadir view angle and infinitesimal field of view. The vertical optical 
#' thickness for all aerosol simulations was 0.5, while for Rayleigh, was the
#' Rayleigh optical thickness at the given surface atmospheric pressure at 450 
#' nm. It is noted, however, that dependence os the normalized APSF on optical 
#' thickness is small. By construction, aerosol APSF have no pressure 
#' dependence. All simulations are for a sensor altitude of 800 km, but altitude
#' dependence is only significant at altitudes lower than 10 km (aircraft, 
#' drones).
#'
#' Simulations run with the annular accumulator geometry.
#'
#' @format A list with the following components:
#' \itemize{
#'   \item{mar:}{ APSF simulation for the maritime aerosol type.}
#'   \item{con:}{ APSF simulation for the continental aerosol type.}
#'   \item{urb:}{ APSF simulation for the urban aerosol type.}
#'   \item{ray:}{ A list with APSFs of Rayleigh scattering for different surface pressures.}
#'   \itemize{
#'     \item{p1100:}{ Rayleigh APSF simulation at 1100 mbar.}
#'     \item{p1013:}{ Rayleigh APSF simulation at 1013.25 mbar.}
#'     \item{p0750:}{ Rayleigh APSF simulation at 750 mbar.}
#'     \item{p0500:}{ Rayleigh APSF simulation at 500 mbar.}
#'   }
#' }
#'
#' @source Generated with the apsfs library.

"asim"

#' Sectorial simulations
#'
#' Contains the atmospheric point spread functions (APSFs) for nadir, 30 and 60 
#' degrees vieweing sensors, simulated for the continental (con) standard aerosol 
#' model of the 6SV radiative transfer code. Simulations were carried out 
#' without Rayleigh contribution, and with an aerosol scale height of 2 km.
#' Sensor was positioned at TOA, with an infinitesimal field of view. The 
#' vertical optical thickness for all aerosol simulations was 0.5. It is noted, 
#' however, that dependence os the normalized APSF on optical thickness is small. 
#' By construction, aerosol APSF have no pressure dependence. All simulations 
#' are for a sensor altitude of 800 km, but altitude dependence is only 
#' significant at altitudes lower than 10 km (aircraft, 
#' drones).
#'
#' Simulations run with the sectorial accumulator geometry.
#'
#' @format A list with the following components:
#' \itemize{
#'   \item{conv00:}{ APSF simulation for the continental aerosol type at 0 degree view angle;}
#'   \item{conv30:}{ APSF simulation for the continental aerosol type at 30 degree view angle;}
#'   \item{conv60:}{ APSF simulation for the continental aerosol type at 60 degree view angle.}
#' }
#'
#' @source Generated with the apsfs library.

"ssim"

#' Grid simulations
#'
#' Contains the atmospheric point spread functions (APSFs) for nadir, 30 and 60 
#' degrees vieweing sensors, simulated for the continental (con) standard aerosol 
#' model of the 6SV radiative transfer code. Simulations were carried out 
#' without Rayleigh contribution, and with an aerosol scale height of 2 km.
#' Sensor was positioned at TOA, with an infinitesimal field of view. The 
#' vertical optical thickness for all aerosol simulations was 0.5. It is noted, 
#' however, that dependence os the normalized APSF on optical thickness is small. 
#' By construction, aerosol APSF have no pressure dependence. All simulations 
#' are for a sensor altitude of 800 km, but altitude dependence is only 
#' significant at altitudes lower than 10 km (aircraft, 
#' drones).
#'
#' Simulations run with the grid accumulator geometry.
#'
#' @format A list with the following components:
#' \itemize{
#'   \item{conv00:}{ APSF simulation for the continental aerosol type at 0 degree view angle;}
#'   \item{conv30:}{ APSF simulation for the continental aerosol type at 30 degree view angle;}
#'   \item{conv60:}{ APSF simulation for the continental aerosol type at 60 degree view angle.}
#' }
#'
#' @source Generated with the apsfs library.

"gsim"

#' Fitted annular simulations
#'
#' Contains the fitted coefficients of the atmospheric point spread functions
#' (APSFs) simulated for the different standard aerosol models of the 6SV 
#' radiative transfer code: maritime (mar), continental (con) and urban (urb). 
#' Simulations were carried out without Rayleigh contribution, and with an 
#' aerosol scale height of 2 km. The APSF for Rayleigh alone (ray) is also 
#' available. Sensor was positioned at TOA, with nadir view angle and 
#' infinitesimal field of view. The fitted model can be expanded into a grid 
#' spatial representation with the function \code{apsfs::predict_grid}. The 
#' vertical optical thickness for all simulations was 0.5, but dependence on 
#' optical thickness is small. By construction, aerosol APSF have no pressure 
#' dependence. The fitted Rayleigh simulations include pressure dependence with 
#' limits of 1100 to 500 mbar. See \code{?apsfs::fit_annular} for further 
#' details.
#'
#' @source Generated with the apsfs library \url{https://github.com/AlexCast/apsfs}.

"fasim"

#' Fitted sectorial simulations
#'
#' Contains the fitted coefficients of the atmospheric point spread functions
#' (APSFs) simulated for the continental standard aerosol model of the 6SV 
#' radiative transfer code. Simulations were carried out without Rayleigh 
#' contribution, and with an aerosol scale height of 2 km. Sensor was positioned
#' at TOA, with three view angles (0, 30 and 60 degrees) and  had an infinitesimal 
#' field of view. The fitted model can be expanded into line, sectors or grid 
#' spatial representation with the functions \code{apsfs::predict_annular}, 
#' \code{apsfs::predict_sectorial} and \code{apsfs::predict_grid}. The vertical 
#' optical thickness for all simulations was 0.5, but dependence on optical 
#' thickness is small. By construction, aerosol APSF have no pressure 
#' dependence. See \code{?apsfs::fit_sectorial} for further details.
#'
#' @source Generated with the apsfs library \url{https://github.com/AlexCast/apsfs}.

"fssim"
