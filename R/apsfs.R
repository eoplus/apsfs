#' \pkg{apsfs}: Atmospheric Point Spread Function Simulation
#'
#' This package contains functions to calculate the atmopsheric point spread 
#' function (PSF) for a given geometry and atmospheric composition and profile 
#' and functions to fit this data to annular or grid models.
#'
#' The PSF is calculated with backward Monte Carlo with the following 
#' simplifications: 
#' \itemize{
#'  \item Plane paralel geometry;
#'  \item Layered medium;
#'  \item Lambertian surface;
#'  \item Molecular absorption not included.
#' }
#'
#' Simulations can be recorded in a annular geometry for symmetric conditions 
#' (Lambertian surfaces and sensor looking at nadir) or grid geometry for 
#' asymmetric conditions (surface BRDF and or zenith view angles away from 
#' nadir).
#'
#' The results can then be fitted to models to provide flexibility for 
#' application. Cumulative annular data is fitted to a two term exponential 
#' function, while grid data is fitted with Zernike polynomials.
#'
#' Below is the list of exported functions in this package:
#' \itemize{
#'  \item get_mol_par: Calculate the molecular absorption and optical depth
#'  \item get_opt_atm_prfl: Calculate optical properties of atmospheric layers 
#'        in high resolution
#'  \item cos2sph: Directional cosines to polar angles
#'  \item mc_psf: Solves the radiative transfer with the Monte Carlo method to 
#'        simulate the PSF
#'  \item dpsf: Calculates the PSF density (per area, in m^2)
#'  \item fit_radial_psf: Fit a cumulative annular model to radial PSF geometry
#'  \item fit_grid_psf: Fit grid PSF geometry with Zernike polynomials
#'  \item calc_cdf: Calculate the cumulative density distribution from 
#'        scattering phase functions
#'  \item ray_ph: Rayleigh scattering phase function
#'  \item scat_aer: Get scattering direction after aerosol interaction
#'  \item scat_ray: Get scattering direction after Rayleigh interaction
#'  \item rayleigh_od: Calculate the Rayleigh optical depth at a given pressure 
#'        and altitude
#'  \item update_cdir: Update directional cosines after scattering
#' }
#'
#' The following data is available:
#' \itemize{
#'  \item us62: US standard atmospheric profile (1962)
#'  \item continental_ph_6sv: Continental aerosol scattering phase function
#'  \item maritime_ph_6sv: Maritime aerosol scattering phase function
#'  \item urban_ph_6sv: Urban aerosol scattering phase function
#'  \item continental_coef_6sv: Continental aerosol optical coefficients
#'  \item maritime_coef_6sv: Maritime aerosol optical coefficients
#'  \item urban_coef_6sv: Urban aerosol optical coefficients
#'  \item tau_ray_6sv: Rayleigh optical depth for US62 and 6SV wavelength grid
#'  \item rayleigh_cdf: Cumulative distribution function of the Rayleigh 
#'        scattering phase function
#' }
#'
#' @author Alexandre Castagna
#' @docType package
#' @name apsfs-package
NULL
