
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
#' Fits an exponential model to the radius cumulative integral of a PSF 
#' accumulated in annular geometry. 
#' 
#' @param psfm An object produced by the function \code{mc_psf}.
#' @param norm Logical. Should the PSF be normalized to sum to unity?
#'
#' @details 
#' The code calculates the cummulative integral of the annular PSF and fits the 
#' following model by optimization:
#'
#' F(r) = c1 - (c2 * e^(c3 * r) + c4 * e^(c5 * r) + c6 * e^(c7 * r)),
#'
#' where c2 + c4 + c6 = c1.
#'
#' The cummulative annular integral PSF is as known as "environmental function".
#' The optimization is made with mean absolute relative error (MARE) as the 
#' error function.
#'
#' Note that 'PSF' in this package refers to the diffusely transmitted photons 
#' only. Photons from direct transmission are not included in the caluclations.
#'
#' @return
#' A list with the follwoing components: "type": type of the fitted model, 
#' "coefficients": named vector with the 5 coefficients, and "mare" with the 
#' value of the mean absolute relative error of the model fit.
#'
#' @export

fit_annular_psf <- function(psfm, norm = TRUE) {

  if(norm) psfm$bin_phtw <- psfm$bin_phtw / sum(psfm$bin_phtw)

  # Build function to be optimized:
  optfun <- function(x, ftot, r, finf, error = T) {
    est  <- finf - (x[1] * exp(x[2] * r) + (finf - x[1]) * x[3] * exp(x[4] * r) + 
      (finf - x[1]) * (1 - x[3]) * exp(x[5] * r))
    if(error) {
      return(mean(abs(est - ftot) / ftot))
    } else {
      return(est)
    }
  }

  # Fit model:
  # First poisition will always be zero and will be singular in the MARE 
  # calculation.
  x    <- psfm$bin_brks[-1]
  y    <- cumsum(psfm$bin_phtw)
  st   <- c(0.5, -0.234, 0.5, -3.1, -1)
  opt  <- optim(st, optfun, method = "Nelder-Mead", ftot = y, r = x, 
    finf = sum(psfm$bin_phtw), control = list(maxit = 1E6))

  coef <- opt$par

  fit <- list(
    type = "annular",
    coefficients = c(
      c1 = sum(psfm$bin_phtw), 
      c2 = coef[1], 
      c3 = coef[2], 
      c4 = (sum(psfm$bin_phtw) - coef[1]) * coef[3], 
      c5 = coef[4],
      c6 = (sum(psfm$bin_phtw) - coef[1]) * (1 - coef[3]),
      c7 = coef[5]
    ),
    mare = opt$value,
    convergence = opt$convergence
  )

  fit

}

#' Predict annular PSF fitted model
#'
#' Solves a fitted model from \code{fit_annular_psf} for the PSF of the radial 
#' cummulative PSF at requested radius.
#'
#' @param r   The radius distances (km) at which the model should be evaluated.
#' @param fit A model fit from \code{fit_annular_psf}.
#' @param cum Logical. Should the cummulative PSF be returned? Defaults to FALSE.
#'
#' @details If cum = FALSE, the derivative of the cummulative PSF model fit will 
#' be returned, which is the PSF at the desired radius.
#'
#' @return A numeric vetor with the PSF (or cummulative PSF).
#'
#' @export

predict_annular <- function(r, fit, cum = FALSE) {
  if(cum) {
    .pred_annular_cum(r, fit)
  } else {
    .pred_annular_psf(r, fit)
  }
}

.pred_annular_cum <- function(r, fit) {
  coeff <- fit$coefficients
  coeff[1] - (coeff[2] * exp(coeff[3] * r) + coeff[4] * exp(coeff[5] * r) + 
    coeff[6] * exp(coeff[7] * r))
}

.pred_annular_psf <- function(r, fit) {
  numDeriv::grad(.pred_annular_cum, x = r, fit = fit)
}


#' PSF grid density model
#'
#' under development...
#'
#' @export

fit_grid_dpsf <- function(psfm, norm = TRUE) {


#  xy  <- psfm$bin_mid
#  xy  <- xy[-c(1, length(xy))]
#  xy  <- expand.grid(y = xy, x = xy)[, 2:1]

#  wf <- an_psf(psfm)
#  wf <- wf[, -c(1, ncol(wf))]
#  wf <- wf[-c(1, nrow(wf)), ]

#  df <- cbind(xy, wf = as.vector(wf))
#  coordinates(df) <- ~x+y


#  xo  <- seq(0.03 / 2, psfm$metadata$ext, by = 0.03)
#  xo  <- c(rev(-xo), xo)
#  xo  <- (xo[-1] + xo[-length(xo)]) / 2
#  tp  <- akima::interp(df, z = 'wf', xo = xo, yo = xo, linear = T, nx = length(xo), ny = length(xo))

  require(zernike)

  xy  <- psfm$bin_mid
  xy  <- xy[-c(1, length(xy))]
  xy  <- expand.grid(xy, xy)
  rho <- sqrt(apply(xy^2, 1, sum))
  id  <- which(rho > 1)
  rho <- rho[-id]

  theta <- acos(apply(xy, 1, function(x) { sum(x * c(0, 1)) / sqrt(sum(x^2)) }))
  theta[which(is.na(theta))] <- 0
  theta <- matrix(theta, ncol = ncol(psfm$bin_phtw)-2)
  nr <- (nrow(theta)-1) / 2
  theta[1:nr,] <- (2 * pi - theta[1:nr,])
  dim(theta) <- NULL  
  theta <- theta[-id]

  wf <- psfm$bin_phtw
  wf[rho < 0.1] <- 1E-5
#  wf <- log(-log(wf))
#  wf[is.infinite(wf)] <- -6
  wf <- wf[, -c(1, ncol(wf))]
  wf <- wf[-c(1, nrow(wf)), ]
  dim(wf) <- NULL
  wf <- wf[-id]


  zfit  <- zernike::fitzernikes(wf, rho, theta, phi = 0, maxorder = 40, uselm = TRUE)


  return(zfit)
}

