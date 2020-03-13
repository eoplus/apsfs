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
#' F(r) = c1' - (c2' * e^(c3' * r) + c4' * e^(c5' * r) + c6' * e^(c7' * r)),
#' \describe{
#'   \item{c1' =}{ c1,}
#'   \item{c2' =}{ c2 / p,}
#'   \item{c3' =}{ c3 / p,}
#'   \item{c4' =}{ (c1' - c2') * c4,}
#'   \item{c5' =}{ c5,}
#'   \item{c6' =}{ (c1' - c2') * (1 - c4),}
#'   \item{c7' =}{ c6,}
#'   \item{p   =}{ press / 1013.25,}
#' }
#' and c1 to c6 are fitted parameters. Note that c1 = total diffuse transmittance
#' (or 1, if norm = TRUE) and that c2' + c4' + c6' = c1.
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
#' \describe{
#'   \item{type:}{Type of the fitted model;}
#'   \item{coefficients:}{Named vector with the 6 coefficients;}
#'   \item{mare:}{The minimum mean absolute relative error of the model fit;}
#'   \item{mare_seq:}{A sorted sequence of MARE of different random starts;}
#'   \item{convergence:}{Convergence code from \code{optim};}
#'   \item{press_dep:}{Logical flag indicating if there was pressure dependency;}
#'   \item{press_rng:}{The range of pressure values if fitted with pressure dependency.}
#' }
#'
#' @examples
#' # Fitting continental aerosol annular APSF simulation:
#'
#' data(asim)
#' psfm  <- asim["con"] # Note that input to function must be a list!
#' opt   <- fit_annular(psfm, norm = TRUE, nstart = 30)
#' opt
#' plot(opt$mare_seq, xlab = "Iteration", ylab = "MARE")
#'
#' # Fitting Rayleigh annular APSF simulation with pressure dependence:
#'
#' data(asim)
#' psfm  <- asim[["ray"]] # A list of Rayleigh simulations at different pressures 
#' opt   <- fit_annular(psfm, press = TRUE, norm = TRUE, nstart = 30)
#' opt
#' par(mfcol = c(1, 2))
#' plot(opt$mare_seq, xlab = "Iteration", ylab = "MARE")
#' plot(NA, xlim = c(0, 10), ylim = c(0, 1), xaxs = "i", yaxs = "i", 
#'   xlab = "Radius (km)", ylab = "Normalized f(r)")
#' x <- asim[["ray"]][[1]]$bin_brks
#' cols <- rev(rainbow(1:length(asim[["ray"]]), start = 0, end = 0.8))
#' for(i in 1:length(asim[["ray"]])) {
#'   points(x, cum_psf(asim[["ray"]][[i]], norm = TRUE), col = "grey")
#' }
#' for(i in 1:length(asim[["ray"]])) {
#'   press <- asim[["ray"]][[i]]$metadata$press
#'   lines(x, predict_annular(x, opt, type = "cumpsf", press = press))
#' }
#'
#' @export

fit_annular <- function(psfm, press = FALSE, norm = TRUE, nstart = 10) {

  if(length(psfm) > 1 & press) norm <- TRUE
  if(!press & length(psfm) > 1) {
    "psfm length > 1 but press set to FALSE. Only first psf will be used" %>%
    warning(call. = FALSE)
  }

  if(norm) {
    for(i in 1:length(psfm)) {
      psfm[[i]]$bin_phtw <- psfm[[i]]$bin_phtw / sum(psfm[[i]]$bin_phtw)
    }
  }

  if(press) {
    # Check that all simulations are at the same except for pressure:
    meta <- lapply(psfm, function(x) {x$metadata[c("res", "ext", "snsznt", "snsfov", "snspos")]})
    for(i in 2:length(meta)) {
      for(j in 1:4) {
        if(!identical(meta[[1]][[j]], meta[[i]][[j]])) {
          paste("All simulations included in a given fit with pressure dependence", 
            "must have the same parameters: res, ext, snsznt, snsfov, snspos") %>%
          stop(call. = FALSE)
        }
      }
    }

    n     <- length(psfm[[1]]$bin_mid)
    press <- sapply(psfm, function(x) { x$metadata$press }) %>%
             rep(each = n) %>%
             `/`(., 1013.25)
    pdep  <- TRUE
  } else {
    press <- 1
    pdep  <- FALSE
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
    y <- c(y, cum_psf(psfm[[i]])[-1])
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
#' Solves a fitted model from \code{fit_annular} for the PSF of the radial 
#' cummulative PSF at requested radius.
#'
#' @param r     The radius distances (km) at which the model should be evaluated.
#' @param fit   A model fit from \code{fit_annular}.
#' @param type  Type of prediction: 'psf', 'dpsf', or 'cumpsf'. See Details.
#' @param press Pressure in mbar or NULL. See details.
#'
#' @details If type = cumpsf, the model fit to the cumulative PSF will be 
#' evaluated at the desired radius points. If type = dpsf, the area derivative 
#' the PSF dPSF/dArea is returned. If type == psf, quadrature is used to on the 
#' average area derivative of the annulus and scaled by the area of the annulus. 
#' Note that in this case, r will be sorted and the returned values are for the 
#' mid points of the input vector of radius, so will have a length of 
#' length(r) - 1. Default is to return the PSF.
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

  if(type == "psf") r <- sort(r)
  fun(r = r, press = press, fit = fit)
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
  dpsf <- .pred_annular_den(r = r, press = press, fit = fit)
  psf  <- pi * diff(r^2) * (dpsf[-1] + dpsf[-length(dpsf)]) / 2
  if(r[1] == 0)
    psf[1] <- .pred_annular_den(r = r[2], press = press, fit = fit) * 2 * pi * r[2]
  psf
}

#' Predict PSF grid from fitted model
#'
#' Predicts a PSF grid from annular or grid fitted models.
#'
#' @param f_aer Grid or annular fit to aerosol PSF
#' @param f_aer Grid or annular fit to Rayleigh PSF
#' @param tray  Transmittance due to Rayleigh
#' @param taer  Transmittance due to aerosols
#' @param ext   Maximum spatial extent at geom resolution (km).
#' @param res   Resolution (km), [0,10].
#' @param press Pressure in mbar or NULL.
#'
#' @export

predict_grid <- function(f_aer, f_ray, tray, taer, ext = 0, 
  res = 0.03, press = NULL, view = 0) {

# note: grid model will contain the annular model, since Zernikes will be fit to difference to nadir...

  if(f_aer$type == f_ray$type & f_aer$type == "annular") {
    if(view != 0)
      stop("'view' must be zero with annular fitted models", call. = FALSE)
    .annular_kernel(f_aer, f_ray, tray, taer, ext, res, press)
  } else {
    # NOT IMPLEMENTED YET...
  }

}

.annular_kernel <- function(f_aer, f_ray, tray, taer, ext = 0, res = 0.03, 
  press = NULL) {

  if(ext < (res + res / 2))
    stop("Radial extent has to be equal of larger than res+res/2")

  # Construct the the weighted cumulative function for calculation of 
  # derivative.
  fun_grad <- function(r, press) {
    fa <- predict_annular(r, f_aer, "cumpsf", press)
    fr <- predict_annular(r, f_ray, "cumpsf", press)
    ft <- (tray * fr + taer * fa) / (tray + taer)
    ft[which(is.na(ft))] <- 0
    return(ft)
  }

  # If ex == 0 search the extent that results in 60% of diffuse transmittance
  if(ext == 0) {
    ext <- optim(3, function(x) { abs(fun_grad(x) - 0.6) }, method = "Brent", 
      lower = 0, upper = 100)$par * 1000
  }

  # This section creates a matrix of radial distances to pixel corners
  xvec <- seq(res / 2, ext, by = res)
  xvec <- c(rev(-xvec), xvec)
  yvec <- xvec
  xy   <- expand.grid(xvec, yvec)
  r    <- sqrt(apply(xy^2, 1, sum))
  rm   <- matrix(r, ncol = length(xvec))

  # Calculate the per area weight function at the pixel vertices and calculate
  # the trapezoidal average for each pixel and the weight per pixel given its 
  # area.
  require(numDeriv)
  wm   <- grad(fun_grad, x = rm, press = press) / 2 / pi / rm
  wm   <- wm[-1, ] + wm[-nrow(wm), ]
  wm   <- wm[, -1] + wm[, -ncol(wm)]
  wm   <- (res)^2 * wm / 4
  
  # The weight at the center pixel (target pixel) is more correctly calculated 
  # by evaluating the cumulative function at an r of a circle with equivalent 
  # area.
  rgt <- sqrt(res^2 / pi)
  # id  <- (ncol(wm) + 1) / 2
  # wm[id, id] <- fun_grad(rgt, press = press)
  id <- which(wm == max(wm, na.rm = T))
  wm[id] <- fun_grad(rgt, press = press)

  # Focal function requires an odd number of columns and rows.
  if((nrow(wm) %% 2) == 0) w <- rbind(wm, rep(0, ncol(wm)))
  if((ncol(wm) %% 2) == 0) w <- cbind(wm, rep(0, nrow(wm)))

  return(wm)
}

