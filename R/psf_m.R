
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

  } else if(psfm$metadata$geom == "annular") {

    psfan <- psfm$bin_phtw / pi / 
      (psfm$bin_brks[-1]^2 - psfm$bin_brks[-length(psfm$bin_brks)]^2)

  }

  return(psfan)

}

#' PSF annular cumulative integral model
#'
#' Fits an exponential model to the radius cumulative integral of a PSF in 
#' annular geometry of a given scatter model, possibly with pressure dependency. 
#' 
#' @param psfm   A list of objects created by the function \code{mc_psf}.
#' @param press  Pressure in mbar or NULL. See details.
#' @param norm   Logical. Should the PSF be normalized to sum to unity?
#' @param nstart Number of random start positions for the multidimension 
#'               function optimization.
#'
#' @details
#' The function calculates the (normalized) cumulative integral of the annular 
#' PSF and fits the following model by optimization:
#'
#' F(r) = c1 - (c2 * e^(c3 * r) + c4 * e^(c5 * r) + c6 * e^(c7 * r)),
#'
#' c1 = x1,
#' c2 = x2 / p,
#' c3 = x3 / p,
#' c4 = (c1 - c2) * x4,
#' c5 = x5,
#' c6 = (c1 - c2) * (1 - x4),
#' c7 = x6,
#' p  = press / 1013.25,
#' and x1 to x5 are model parameters. Note that x1 = total diffuse transmittance
#' (or 1, if norm = TRUE) and that c2 + c4 + c6 = c1.
#'
#' If there is no pressure dependency (e.g., aerosol only), 'press' does not 
#' need to be specified. Internally, it will be set to 1013.25 mbar, p will 
#' equal unity and the pressure flag will be set to FALSE in the fit object 
#' metadata. The prediction function \code{predict_annular}, will track 
#' dependency on pressure and model domain.
#'
#' Note that adding multiple models in the input list will only be meaningfull 
#' if they are simulations for the same scatter at different surface pressure 
#' levels. If psfm is a list with more than one component and press is 
#' specified, norm is automatically set to TRUE.
#'
#' The cumulative annular integral PSF is as known as "environmental function".
#' The optimization is made with mean absolute relative error (MARE) as the 
#' error function.
#'
#' Note that 'PSF' in this package refers to the diffusely transmitted photons 
#' only. Photons from direct transmission are simulated and are part of the PSF,
#' but not included in the fit calculations.
#'
#' @return
#' A list with the following components: 
#' "type":         type of the fitted model; 
#' "coefficients": named vector with the 6 coefficients;
#' "mare":         the minimum mean absolute relative error of the model fit; 
#' "mare_seq":     a sorted sequence of MARE of different random starts;
#' "convergence":  convergence code from \code{optim};
#' "press_dep":    logical flag indicating if there was pressure dependency; 
#' "press_rng":    the range of pressure values if fitted with pressure 
#'                 dependency.
#' @examples
#'
#' psfm <- list(
#'  ray_res003_v00_p1100_annular,
#'  ray_res003_v00_p1013_annular,
#'  ray_res003_v00_p0750_annular,
#'  ray_res003_v00_p0500_annular
#' )
#' press <- rep(c(1100, 1013.25, 750, 500), each = length(ray_res003_v00_p0500_annular$bin_mid))
#' opt   <- fit_annular_psf(psfm, press = press, norm = TRUE, nstart = 30)
#' opt
#' plot(opt$mare_seq, xlab = "Iteration", ylab = "MARE")
#'
#' @export

fit_annular_psf <- function(psfm, press = NULL, norm = TRUE, nstart = 10) {

  if(length(psfm) > 1 & any(!is.null(press))) norm <- TRUE
  if(any(is.null(press)) & length(psfm) > 1)
    stop("psfm length > 1 but press not specified.", call. = FALSE)


  if(norm) {
    for(i in 1:length(psfm)) {
      psfm[[i]]$bin_phtw <- psfm[[i]]$bin_phtw / sum(psfm[[i]]$bin_phtw)
    }
  }

  if(is.null(press)) {
    press <- 1
    pdep  <- FALSE
  } else {
    press <- press / 1013.25
    pdep  <- TRUE
  }

  # Build function to be optimized:
  optfun <- function(x, ftot, r, finf, error = T, press) {
    xp1  <- x[1] * press^-1
    xp2  <- x[2] * press^-1
    est  <- finf - (xp1 * exp(xp2 * r) + (finf - xp1) * x[3] * exp(x[4] * r) + 
      (finf - xp1) * (1 - x[3]) * exp(x[5] * r))
    if(error) {
      return(mean(abs(est - ftot) / ftot, na.rm = TRUE))
    } else {
      return(est)
    }
  }

  # Fit model:
  # First poisition will always be zero and will be singular in the MARE 
  # calculation.
  vals <- numeric(nstart + 1)
  x    <- rep(psfm[[1]]$bin_brks[-1], length(psfm))
  y    <- NULL
  for(i in 1:length(psfm)) {
    y <- c(y, cumsum(psfm[[i]][[1]]))
  }
  st   <- c(0.05, -0.14, 0.23, -0.32, -0.05)
  opt  <- optim(st, optfun, method = "Nelder-Mead", ftot = y, r = x, 
    finf = max(y), control = list(maxit = 1E6), press = press)
  vals[1] <- opt$value

  for(i in 1:nstart) {
    st   <- runif(5, 0, 1) * c(1, -1, 1, -1, -1)
    optr <- optim(st, optfun, method = "Nelder-Mead", ftot = y, r = x, 
      finf = max(y), control = list(maxit = 1E6), press = press)
    vals[i + 1] <- optr$value
    if(optr$value < opt$value)
      opt <- optr
  }

  coef <- opt$par

  fit <- list(
    type = "annular",
    coefficients = c(
      c1 = max(y), 
      c2 = coef[1], 
      c3 = coef[2], 
      c4 = coef[3], 
      c5 = coef[4],
      c6 = coef[5]
    ),
    mare = opt$value,
    mare_seq = sort(vals, decreasing = T),
    convergence = opt$convergence,
    press_dep = pdep,
    press_rng = range(press * 1013.25)
  )

  fit

}

#' Predict annular PSF fitted model
#'
#' Solves a fitted model from \code{fit_annular_psf} for the PSF of the radial 
#' cummulative PSF at requested radius.
#'
#' @param r    The radius distances (km) at which the model should be evaluated.
#' @param fit  A model fit from \code{fit_annular_psf}.
#' @param type Type of prediction: 'psf', 'dpsf', or 'cumpsf'. See Details.
#'
#' @details If type = cumpsf, the model fit to the cumulative PSF will be 
#' evaluated at the desired radius points. If type = dpsf, the area derivative 
#' the PSF dPSF/dArea is returned. If type == psf, quadrature is used to on the 
#' average area derivative of the annulus and scaled by the area of the annulus. 
#' Note that in this case, the returned values are for the mid points of the 
#' input vector of radius, so will have a length of length(r) - 1. Default is to 
#' return the PSF.
#'
#' @return A numeric vetor with the PSF, dPSF/dArea or cumulative PSF.
#'
#' @export

predict_annular <- function(r, fit, type = c("psf", "dpsf", "cumpsf"), 
  press = NULL) {

  if(is.null(press) & !fit$press_dep) {
    press <- 1013.25
  } else if(is.null(press) & fit$press_dep) {
    stop("'press' not specified for model fitted with pressure dependency", 
      call. = FALSE)
  }

  if(fit$press_dep) {
    if(any(press < fit$press_rng[1] | press > fit$press_rng[2]))
      paste0("Requested pressure beyond model domain of ", fit$press_rng[1], 
        " to ", fit$press_rng[2]) %>%
      warning(., call. = FALSE)
  }

  press <- press / 1013.25

  fun <- switch(type[1],
                "cumpsf" = .pred_annular_cum,
                "dpsf"   = .pred_annular_den,
                "psf"    = .pred_annular_psf,
                stop("type must be one of 'psf', 'dpsf', or 'cumpsf'", 
                  call. = FALSE)
         )
  fun(sort(r), press, fit)
}

.pred_annular_cum <- function(r, press, fit) {
  cf  <- fit$coefficients
  c1  <- cf[1]
  c2  <- cf[2] * press^-1
  c3  <- cf[3] * press^-1
  c4  <- (c1 - c2) * cf[4]
  c5  <- cf[5]
  c6  <- (c1 - c2) * (1 - cf[4])
  c7  <- cf[6]
  
  res <- c1 - (c2 * exp(c3 * r) + c4 * exp(c5 * r) + c6 * exp(c7 * r))

  res[is.na(res)] <- 0
  res
}

.pred_annular_den <- function(r, press, fit) {
  numDeriv::grad(.pred_annular_cum, x = r, press = press, fit = fit) / 2 / pi / r
}

.pred_annular_psf <- function(r, press, fit) {
  dpsf <- .pred_annular_den(r, press, fit)
  psf  <- pi * diff(r^2) * (dpsf[-1] + dpsf[-length(dpsf)]) / 2
  if(r[1] == 0) psf[1] <- .pred_annular_den(r[2], press, fit) * 2 * pi * r[2]
  psf
}




