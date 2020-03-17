#' Cumulative PSF
#'
#' Calculates the uni (annular) or bidimensional (sectorial, grid) cumulative 
#' PSF.
#'
#' @param psfi An object produced by the function \code{mc_psf}.
#' @param norm Logical. Should the PSF be normalize to integrate to unity?
#'
#' @details The Monte Carlo simulation returns the area integral PSF on the 
#' intervals defined by the break points of the accumulator. This function sum 
#' those contributions evaluating the cumulative distribution at the break 
#' points.
#'
#' @return A numeric vector if the input simulation is in annular geometry and 
#' a matrix if is in sectorial or grid geometry.
#'
#' @export

cum_psf <- function(psfi, norm = TRUE) {

  z  <- psfi$bin_phtw
  if(norm) z <- z / sum(z)

  if(psfi$metadata$geom == "annular") {
    c(0, cumsum(z))
  } else {
    .cum2d_c(z)
  }

}

.cum2d_c <- function(z) {

  cum <- apply(z, 1, cumsum) %>%
         t() %>%
         apply(., 2, cumsum)

  cbind(rep(0, nrow(z)), cum) %>%
  rbind(rep(0, ncol(z) + 1), .)

}

#' PSF normalized per bin area
#'
#' Normalized the bin integral per area of bin used in the Monte Carlo 
#' accumulation.
#'
#' @param psfm An object produced by the function \code{mc_psf}.
#' @param norm Logical. Should the PSF be normalized to sum to unity?
#'
#' @details PSF are normalized per area of each bin depending on geometry. If 
#' norm is TRUE, PSF is normalize to sum to 1 before calculation and the 
#' integral of the area normalized PSF will be 1.
#' 
#' Note that the bin containing data from the last requested distance to 
#' infinite will have a value of zero per area PSF.
#'
#' @return
#' A vector with the average per area PSF in each bin.
#' 
#' @export

dpsf <- function(psfm, norm = TRUE) {

  if(norm) psfm$bin_phtw <- psfm$bin_phtw / sum(psfm$bin_phtw)
  psfm$bin_brks <- psfm$bin_brks * 1E3

  if(psfm$metadata$geom == "grid") {

    mw    <- matrix(diff(psfm$bin_brks), ncol = ncol(psfm$bin_phtw), 
      nrow = nrow(psfm$bin_phtw))
    psfan <- psfm$bin_phtw / (mw * t(mw))

  } else if(psfm$metadata$geom == "sectorial") {

    psfan <- psfm$bin_phtw / (1 * pi / 180) /
      (psfm$bin_brks[-1]^2 - psfm$bin_brks[-length(psfm$bin_brks)]^2)

  } else if(psfm$metadata$geom == "annular") {

    psfan <- psfm$bin_phtw / pi / 
      (psfm$bin_brks[-1]^2 - psfm$bin_brks[-length(psfm$bin_brks)]^2)

  }

  return(psfan)

}

