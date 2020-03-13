
#' Get molecular absorption profile
#'
#' Calculates the molecular absorption profile for a given atmosphere. 
#'
#' @param atm    Profile of atmospheric temperature, pressure, water vapor and 
#'               ozone.
#' @param lambda Wavelength (nm) of observation.
#'
#' @details Not implemented. At the moment will return a vector of zeros (zero
#' absorption) with equal length to the atmospheric profile.
#'
#' @return A vector with the molecular absorption profile at the same heights as
#' the atmospheric profile.
#'
#' @export

get_mol_par <- function(atm, lambda) {

  # Molecular absorption not implemented; for now:
  a_mol_z  <- rep(0, nrow(atm))
  tau_mol  <- 0 

  return(list(a_mol_z = a_mol_z, tau_mol = tau_mol))
}

#' Get optical properties of atmospheric layers
#'
#' Calculates the average optical properties per atmosphere layer for a given 
#' atmospheric profile.
#'
#' @param atm       A data.frame with columns "Z" (Km height from surface) and 
#'                  "Pressure" (mbar).
#' @param tau_aer   Aerosol optical thickness (unitless), [0,Inf).
#' @param H_aer     Aerosol scale height (Km), [0,Inf).
#' @param w0_aer    Aerosol single scattering albedo (unitless), [0,1].
#' @param tau_ray_z Rayleigh optical thickness (unitless), [0,Inf).
#' @param a_mol_z   Molecular absorption profile at the same Z and Pressure 
#'                  points of atm.
#' @param nlayers   Number of layers of the output.
#'
#' @details
#' Linear interpolation is used for the atmospheric profile of pressure and 
#' molecular absorption. The atmosphere is returned as nlayers between the 
#' surface (0 km) and the top of the input atmospheric profile.
#'
#' Pseudo-code:
#' 1 - Get break points of atmosphere layers;
#' 2 - Based on the exponential profile, integrate the aerosol extinction from 
#'     TOA to surface.
#' 3 - Calculate average aerosol properties per layer.
#' 4 - Based on the pressure profile, get the Rayleigh optical depth profile.
#' 5 - Calculate average Rayleigh scattering per layer.
#' 6 - Integrate molecular absorption from TOA to surface.
#' 7 - Calculate the average molecular absorption per layer.
#' 8 - Calculate the average extinction per layer and single scattering albedo.
#' 9 - Sum the total optical depth from TOA per layer.
#'
#' @return A data.frame with the average optical properties per layer. Note that 
#' the first two columns are the the break points of the atmosphere layers in 
#' km and optical depth. So the last row will always be NA for all average properties,
#' since the last atmosphere layer is nrow - 1.
#'
#' @examples
#'
#' data(us76)
#' tau_ray_z <- rayleigh_od(atm = us76, lambda = 550, co2 = 400, lat = 45)
#' a_mol_z   <- rep(0, nrow(us76))
#' get_opt_atm_prfl(atm = us62, tau_aer = 0.5, H_aer = 2, w0_aer = 0.89, 
#'   tau_ray_z = tau_ray_z, a_mol_z = a_mol_z)
#'
#' @export

get_opt_atm_prfl <- function(atm, tau_aer, H_aer, w0_aer, tau_ray_z, a_mol_z, 
  nlayers = 1000) {
 
  z     <- sort(atm$Z, decreasing = TRUE)
  brk_z <- rev(c(0, 10^seq(log10(0.01), log10(z[2]), length.out = nlayers - 1), 
    z[1]))
  brk_z[2] <- z[2] # Avoid rounding errors... 

  # Integrate aerosol extinction from TOA to surface:
  c_aer_0 <- tau_aer / H_aer

  tau_aer_fun <- function(x) {
    c_aer_0 * exp(-x / H_aer)
  }

  tau_aer_z <- numeric(length(brk_z))

  for(i in 1:length(brk_z)) {
    tau_aer_z[i] <- integrate(tau_aer_fun, lower = brk_z[i], upper = Inf)$value
  }

  tau_aer_z[length(tau_aer_z)] <- tau_aer

  # Get average aerosol properties in each layer:
  c_aer_z <- diff(tau_aer_z) / abs(diff(brk_z))
  b_aer_z <- c_aer_z * w0_aer
  a_aer_z <- c_aer_z * (1 - w0_aer)

  # Interpolate Rayleigh optical depth:
  if(all(tau_ray_z == 0)) {
    tau_ray_z <- rep(0, length(brk_z))
  } else {
    tau_ray_z <- 10^approx(x = atm$Z, y = log10(tau_ray_z), xout = brk_z)$y
  }

  # Get average Rayleigh properties in each layer:
  b_ray_z <- diff(tau_ray_z) / abs(diff(brk_z))

  # Integrate molecular absorption from TOA to surface:
  tau_mol_z <- numeric(length(brk_z))

  if(sum(a_mol_z) != 0) {
    tau_mol_fun <- function(x) {
      10^approx(x = us62$Z, y = log10(a_mol_z), xout = x)$y
    }

    for(i in 2:length(brk_z)) {
      xout <- seq(brk_z[i], max(brk_z), length.out = (max(brk_z) - brk_z[i]))
      tp   <- tau_mol_fun(xout)
      tau_mol_z[i] <- sum(diff(xout) * (tp[-1] + tp[-length(tp)]) / 2)
    }

  }

  # Get average molecular absorption properties in each layer:
  a_mol_z  <- diff(tau_mol_z) / abs(diff(brk_z))

  # Get total attenuation and total single scattering albedo per layer:
  c_atm_z  <- diff(tau_ray_z + tau_aer_z + tau_mol_z) / abs(diff(brk_z))
  w0_atm_z <- diff(tau_ray_z + tau_aer_z * w0_aer) / abs(diff(brk_z)) / c_atm_z

  # Get total tau from TOA to surface:
  tau_atm_prfl <- tau_ray_z + tau_aer_z + tau_mol_z
  
  tau_atm_prfl <- data.frame(km     = brk_z,                                    # Metric distance break points
                             tau    = tau_atm_prfl,                             # Optical distance break points
                             w0_tot = c(w0_atm_z, NA),                          # average single scattering albedo of the layer
                             b_ray  = c(b_ray_z, NA),                           # average Rayleigh scattering coefficient of the layer
                             b_aer  = c(b_aer_z, NA),                           # average Aerosol scattering coefficient of the layer
                             c_tot  = c(c_atm_z, NA))                           # average total attenuation coefficient of the layer

  attr(tau_atm_prfl, "press") <- max(atm$Pressure)
  tau_atm_prfl
}

