#' Rotate azimuth of grid simulations
#'
#' Rotate the grid simulation to the requested viewing azimuth.
#'
#' @param psf   A grid psf object.
#' @param vaz   Viewing azimuth (rad).
#'
#' @export

rotate_grid <- function(psf, vaz) {
  z  <- psf$bin_phtw 
  z  <- raster::raster(z[-c(1, nrow(z)), -c(1, ncol(z))])
  
  vaz <- 90 - (vaz * 180 / pi)
  raster::crs(z) <- "+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=0"
  ncrs   <- paste0("+proj=aeqd +ellps=sphere +lat_0=90 +lon_0=", -vaz)
  res    <- raster::projectRaster(z, res = raster::res(z), crs = ncrs)
  res[is.na(res)] <- 0
  raster::as.matrix(res)
}

#' Predict SAF or PSF grid from fitted model
#'
#' Predicts a SAF or PSF grid from annular or grid fitted models.
#'
#' @param f_aer Grid or annular fit to aerosol SAF or PSF
#' @param f_aer Grid or annular fit to Rayleigh SAF or PSF
#' @param wray  Weight (transmittance or reflectance) due to Rayleigh
#' @param waer  Weight (transmittance or reflectance) due to aerosols
#' @param ext   Maximum spatial extent at geom resolution (km).
#' @param res   Resolution (km), [0,10].
#' @param press Pressure in mbar or NULL. For annular fits only.
#' @param tpred Third predictor level or NULL. For sectorial fits only.
#' @param vaz   Viewing azimuth (rad). For sectorial geometry only.
#'
#' @details
#' The grid for each scatter type are calculated individually, then their 
#' weighted average is calculated based on the weights for Rayleigh and aerosol.
#' The weights are given by the spherical albedo (for SAF) and diffuse upward 
#' transmittance (for PSF).
#'
#' @export

predict_grid <- function(f_aer, f_ray, wray, waer, ext = 0, 
  res = 0.03, press = NULL, tpred = NULL, vaz = pi/2) {

  if(f_aer$type == f_ray$type & f_aer$type == "annular") {
    if(!is.null(tpred))
      stop("'tpred' must be NULL with annular fitted models", call. = FALSE)
    .annular_kernel(f_aer, f_ray, wray, waer, ext, res, press)
  } else if(f_aer$type == f_ray$type & f_aer$type == "sectorial"){
    if(vaz < 0 | vaz > 2*pi)
      stop("'vaz' must be between 0 and 2*pi", call. = FALSE)
    .sectorial_kernel(f_aer, f_ray, wray, waer, ext, res, tpred, vaz = vaz)
  }

}

.annular_kernel <- function(f_aer, f_ray, wray, waer, ext = 0, res = 0.03, 
  press = NULL) {

  if(ext < (res + res / 2))
    stop("Radial extent has to be equal of larger than res+res/2")

  # Construct the the weighted cumulative function for calculation of 
  # derivative.
  fun_grad <- function(r, press) {
    fa <- predict_annular(r, f_aer, "cpsf", press)
    fr <- predict_annular(r, f_ray, "cpsf", press)
    ft <- (wray * fr + waer * fa) / (wray + waer)
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

  wm[wm < 0] <- 0
  return(wm)
}

.sectorial_kernel <- function(f_aer, f_ray, wray, waer, ext = 0, res = 0.03, 
  tpred = NULL, vaz = pi/2) {

  if(ext < (res + res / 2))
    stop("Radial extent has to be equal of larger than res+res/2")

  # Construct the the weighted cumulative function for calculation of 
  # derivative.
  fun_grad <- function(r, a, tpred) {
    fa <- .pred_sectorial_cum_VEC(r, a, tpred, fit = f_aer)
    fr <- .pred_sectorial_cum_VEC(r, a, tpred, fit = f_ray)
    ft <- (wray * fr + waer * fa) / (wray + waer)
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

#  a    <- acos(xy[, 2] / r)
#  am   <- matrix(a, ncol = sqrt(length(a)))
#  am[1:((nrow(am)-1)/2), ] <- 2 * pi - am[1:((nrow(am)-1)/2), ]
  a    <- acos(-xy[, 1] / r)
  am   <- matrix(a, ncol = sqrt(length(a)))
  am[, 1:((ncol(am)-1)/2)] <- 2 * pi - am[, 1:((ncol(am)-1)/2)]
  am <- am - vaz
  am[am > 2*pi] <- am[am > 2*pi] - 2*pi
  am[am < 0]    <- 2*pi + am[am < 0]

  # Calculate the per area weight function at the pixel vertices and calculate
  # the trapezoidal average for each pixel and the weight per pixel given its 
  # area.
  w    <- .get_d2psf_dxdy_VEC(fun_grad, x1 = as.vector(rm), x2 = as.vector(am), 
    tpred = tpred) / r
  wm   <- matrix(w, ncol = length(xvec))
  wm   <- wm[-1, ] + wm[-nrow(wm), ]
  wm   <- wm[, -1] + wm[, -ncol(wm)]
  wm   <- (res)^2 * wm / 4
  
  # The weight at the center pixel (target pixel) is more correctly calculated 
  # by evaluating the cumulative function at an r of a circle with equivalent 
  # area.
  rgt <- sqrt(res^2 / pi)
  id  <- ceiling(ncol(wm) / 2)
  ra  <- expand.grid(rgt, (180:360) * pi / 180)
  # The correct way to calculate the center pixel value would be:
  # wm[id] <- fun_grad(rgt, 2*pi, tpred)
  # but this type of model has large errors at the radius boundaries...
  # Although not proper entirelly justifiable, tests showed that using the 
  # average at all angles result in less error...
  wm[id,id] <- mean(fun_grad(ra[, 1], ra[, 2], tpred))

  wm[wm < 0] <- 0
  return(wm)
}
