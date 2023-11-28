#' Fit sectorial simulation to a bi or three dimensional model
#'
#' Fits sectorial APSF(s) with a combination of Singular Value Decomposition and 
#' polynomials.
#' 
#' @param psfl  A list of psf objects.
#' @param norm  Logical: should the psf be normalized to integrate to unity? 
#'              Forced to TRUE when tpred is specified.
#' @param tpred Optional third predictor. See Details.
#' @param pord  Polynomial order. 
#' @param ncmp  Number of Singular Value Decomposition components to use.
#'
#' @details The procedure for fitting APSF as depending on radius, azimuth and 
#' possibly a third variable follow the method outlined by Mandel (1981). The 
#' APSF is first decomposed with Singular Value Deconmposition and the first 
#' \code{ncmp} components of the u and v matrices are fitted with polynomials of 
#' order \code{pord} to the original independent variables radius and azimuth 
#' when no third predictor is provided.
#'
#' If a third predictor is provided (e.g., view angle), the APSFs are combined 
#' by row (long table format for radius and the third predictor) and the each u
#' component is arranged in matrix format as a function of radius and the third
#' predictor and again decomposed by SDV, with its u and v components fitted by 
#' polynomials as a function of radius and the third predictor. The third 
#' predictor can be specified by \code{tpred} as one of: "pressure", "view",
#' "altitude" or "fov".
#'
#' The maximum value of ncmp is 361 due to the fixed angular resolution of the 
#' sectorial geometry. The total number of parameters used to fit the surface(s) 
#' depends on ncmp and pord. Without a third predictor it will be:
#' ncmp + 2 * ncmp * pord. 
#'
#' If a third predictor is requested, the maximum number of parameters will be:
#' ncmp + ncmp^2 + (2 * ncmp + 1) * pord. The final number will be lower 
#' depending if the polynomial order is larger than the number of independent 
#' levels in tpred.
#'
#' @return A list with model parameters.
#'
#' @seealso \code{predict_sectorial}, \code{predict_annular}
#'
#' @references
#' Mandel, J. 1981. Fitting Curves and Surfaces with Monotonic and Non-Monotonic 
#'   Four Parameter Equations. Journal of research of the National Bureau of 
#'   Standards Vol. 86, 1, 1-25.
#'
#' @examples
#' # Fitting a single sectorial APSF:
#' data(ssim)
#' fit1 <- fit_sectorial(ssim[3])
#' prd1 <- predict_sectorial(r = ssim[[1]]$bin_brks, fit = fit1, type = "cumpsf")
#' sqrt(mean((prd1 - cum_psf(ssim[[3]]))^2, na.rm = T)) # RMSE
#'
#' # Fitting a single sectorial APSFs with view angle dependence:
#' fit2 <- fit_sectorial(ssim, tpred = "view")
#' prd2 <- predict_sectorial(r = ssim[[1]]$bin_brks, tpred = 60 * pi / 180, 
#'   fit = fit2, type = "cumpsf")
#' sqrt(mean((prd2 - cum_psf(ssim[[3]]))^2, na.rm = T)) # RMSE
#'
#' @export

fit_sectorial <- function(psfl, norm = TRUE, tpred = NULL, pord = 10, ncmp = 3) {

  x1  <- psfl[[1]]$bin_brks[-length(psfl[[1]]$bin_brks)]
  x2  <- (0:360) * pi / 180
  
  pr_r <- poly(x1, pord, raw = T)
  pr_a <- poly(x2, pord, raw = T)
  # Add intercept as x^0
  pr_r <- cbind(rep(1, nrow(pr_r)), pr_r)
  pr_a <- cbind(rep(1, nrow(pr_a)), pr_a)
  colnames(pr_r) <- colnames(pr_a) <- paste0("x", 0:pord)

  lmr  <- lma <- list()
  form <- as.formula(paste0("z~0+", paste(colnames(pr_r), collapse = "+")))

  if(is.null(tpred)) {
    # Model fitting based only on radius and azimuth.
    # Since u and v vector will be fitted with polynomials, the last integration
    # position, infinite, is removed.
    y  <- cum_psf(psfl[[1]], norm = norm)[-length(psfl[[1]]$bin_brks), ]
    s  <- svd(y)

    # Fit the first ncmp components of the u and v matrices with the original
    # axis.
    for(i in 1:ncmp) {
      z   <- s$u[, i]
      val <- data.frame(pr_r, z = z)
      lmr[[i]] <- lm(form, data = val) %>%
                  coefficients()
      z   <- s$v[, i]
      val <- data.frame(pr_a, z = z)
      lma[[i]] <- lm(form, data = val) %>%
                  coefficients()
    }

    fit <- list(
      type = "sectorial",
      ncmp = ncmp,
      pord = pord,
      coefficients = list(
        scl  = diag(s$d[1:ncmp]),
        lmr  = matrix(unlist(lmr), nrow = pord+1),
        lma  = matrix(unlist(lma), nrow = pord+1)
      ),
      ext = psfl[[1]]$metadata$ext,
      tpred = NULL
    )

    est <- ((pr_r %*% fit$coefficients$lmr) %*% fit$coefficients$scl) %*% 
      t(pr_a %*% fit$coefficients$lma)

    fit$rmse <- sqrt(mean((est - y)^2, na.rm = T))
    y[y == 0] <- NA
    fit$mare <- mean(abs(est - y) / y, na.rm = T)

  } else {

    # Model fitting based on radius, azimuth and a third variable.
    if(length(psfl) == 1) {
      paste0("For sectorial fits with a third predictor, at least two sectorial",
         "psfs must be priovided") %>%
      stop(call. = FALSE)
    }

    meta <- lapply(psfl, function(x) {x$metadata[c("res", "ext")]})
    for(i in 2:length(meta)) {
      for(j in 1:2) {
        if(!identical(meta[[1]][[j]], meta[[i]][[j]])) {
          paste("All simulations included in a given fit with third predictor", 
            "must have the same resolution (res) and extent (ext)") %>%
          stop(call. = FALSE)
        }
      }
    }

    x3 <- switch(tpred,
                 "pressure" = sapply(psfl, function(x) { x$metadata$press }),
                 "view"     = sapply(psfl, function(x) { x$metadata$snsznt }),
                 "altitude" = sapply(psfl, function(x) { x$metadata$snspos[3] }),
                 "fov"      = sapply(psfl, function(x) { x$metadata$snsfov }),
                 stop(paste("tpred must be one of: 'pressure', 'view',", 
                   "'altitude' or 'fov"), call. = FALSE)
         )

    if(sd(x3) == 0) {
      paste("tpred defined as", tpred, ", but no variability is present along",
        "this dimension") %>%
      stop(call. = FALSE)
    }

    if(ncmp > length(x3)) {
      stop(paste("When tpred is defined, ncmp must be equal or lower than the",
        "number of levels in the input tpred"), call. = FALSE)
    }
   
    # Substitute zero for a small number. Polynomials of 0 will be zero, causing
    # a prediction of 0 (e.g., all cells of a nadir view PSF would be zero...). 
    # The absolute zero is substituted here for a small number.
    id <- which(x3 == 0)
    if(length(id) > 0) x3[id] <- 1E-6
    pr_t <- poly(x3, pord, raw = T)
    # Add intercept as x^0
    pr_t <- cbind(rep(1, nrow(pr_t)), pr_t)
    colnames(pr_t) <- colnames(pr_r)
    lmt  <- list()

    # Long table format for radius and third predictor.
    y <- NULL
    for(i in 1:length(psfl)) {
      y  <- rbind(y, cum_psf(psfl[[i]])[-length(psfl[[i]]$bin_brks), ])
    }
    s1 <- svd(y)

    scl2 <- list()
    for(i in 1:ncmp) {
      # Components of structural variable v only depends on azimuth.
      z   <- s1$v[, i]
      val <- data.frame(pr_a, z = z)
      lma[[i]] <- lm(form, data = val) %>%
                  coefficients()

      # Components of the structural variable u depend on radius and the third
      # predictor, and are rearranged in matrix format and decomposed again. 
      # Those second level structural components u and v are dependent on only 
      # one variable: radius or third predictor.
      tmp <- matrix(s1$u[, i], nrow = length(x1), ncol = length(psfl))
      s2  <- svd(tmp)
      scl2[[i]] <- diag(s2$d[1:ncmp]) 
      lmr[[i]] <- lmt[[i]] <- list()
      for(j in 1:ncmp) {
        z   <- s2$u[, j]
        val <- data.frame(pr_r, z = z)
        lmr[[i]][[j]] <- lm(form, data = val) %>%
                         coefficients()

        z   <- s2$v[, j]
        val <- data.frame(pr_t, z = z)
        lmt[[i]][[j]] <- lm(form, data = val) %>%
                         coefficients()
      }
    }

    fit <- list(
      type = "sectorial",
      ncmp = ncmp,
      pord = pord,
      coefficients = list(
        scl1 = diag(s1$d[1:ncmp]),
        scl2 = scl2,
        lma  = matrix(unlist(lma), nrow = pord+1),
        lmr  = lapply(lmr, function(x) { matrix(unlist(x), nrow = pord+1) }),
        lmt  = lapply(lmt, function(x) { matrix(unlist(x), nrow = pord+1) })
      ),
      ext = psfl[[1]]$metadata$ext,
      tpred = tpred,
      tpred_rng = range(x3)
    )
    for(i in 1:ncmp) {
      fit$coefficients$lmt[[i]] <- na.omit(fit$coefficients$lmt[[i]])
    }

    # for reconstruction, reconstruct first level u components first from second 
    # level u and v components. Then reconstruct first level v and scale.
    cf <- fit$coefficients
    est <- NULL
    for(i in 1:ncmp) {
      est <- (pr_r %*% cf$lmr[[i]] %*% cf$scl2[[i]] %*% 
        t(pr_t[, -attr(cf$lmt[[i]],"na.action")] %*% cf$lmt[[i]])) %>%
        as.vector() %>%
        cbind(est, .)
    }
    est <- (est %*% cf$scl1) %*% t(pr_a %*% cf$lma)
    fit$mare <- mean(abs(est - y) / y, na.rm = T)
    fit$rmse <- sqrt(mean(est - y)^2)

  }

  fit
}

#' Predict sectorial PSF fitted model 
#'
#' Solves a fitted model from \code{fit_sectorial} for the PSF of the radial and 
#' azimuthal cummulative PSF at requested radius and azimuth, possibly including 
#' a third predictor variable.
#'
#' @param r     The radial distances (km) at which the model should be evaluated.
#' @param a     The azimuthal distances (rad) at which the model should be 
#'              evaluated.
#' @param fit   A model fit from \code{fit_sectorial}.
#' @param type  Type of prediction: 'psf', 'dpsf', or 'cumpsf'. See Details.
#' @param tpred The level of the third predictor at which the model should be 
#'              evaluated.
#'
#' @details If type = cumpsf, the model fit to the cumulative PSF will be 
#' evaluated at the desired positions. If type = dpsf, the density PSF is 
#' returned (1/km2). If type == psf, quadrature is used to calculate the average 
#' density PSF of the sector and the average is scaled by the area of the sector. 
#' This should equal the Monte Carlo results.  Note that in this case, r will be 
#' sorted and the returned values are for the mid points of the input vector of 
#' positions, so will have a length of length(r) - 1. Default is to return the 
#' PSF.
#'
#' @return A matrix with the PSF, density PSF or cumulative PSF.
#'
#' @seealso \code{fit_sectorial}, \code{fit_annular}, \code{predict_annular}
#'
#' @export

predict_sectorial <- function(r, a = NULL, fit, type = c("psf", "dpsf", "cpsf"), 
  tpred = NULL) {

  if(is.null(a)) a <- (0:360) * pi / 180

  if(fit$type != "sectorial")
    stop("fit must be of sectorial type", call. = FALSE)

  if(is.null(tpred) & !is.null(fit$tpred)) {
    paste("'tpred' level not specified for model fitted with", fit$tpred, 
      "as a third predictor") %>%
    stop(call. = FALSE)
  }

  if(any(r > fit$ext)) {
    paste("Requested radial distance beyond model domain of", fit$ext) %>%
    warning(call. = FALSE)
  }

  fun <- switch(type[1],
                "psf"  = .pred_sectorial_psf,
                "dpsf" = .pred_sectorial_den,
                "cpsf" = .pred_sectorial_cum,
                stop("type must be one of 'psf', 'dpsf', or 'cpsf'", 
                  call. = FALSE)
         )

  if(type == "psf") r <- sort(r)

  fun(r = r, a = a, tpred = tpred, fit = fit)
}

.pred_sectorial_cum <- function(r, a, tpred, fit) {
  pr_r <- poly(r, fit$pord, raw = T)
  pr_a <- poly(a, fit$pord, raw = T)
  pr_r <- cbind(rep(1, nrow(pr_r)), pr_r)
  pr_a <- cbind(rep(1, nrow(pr_a)), pr_a)

  if(!is.null(fit$tpred)) {
    # When fitting a third predictor with zero value, the zero is changed to 
    # a small number defined as 1E-6.
    if(tpred == 0) tpred <- 1E-6
    pr_t <- poly(tpred, fit$pord, raw = T)
  } else {
    pr_t <- NULL
  }

  cf  <- fit$coefficients
  if(is.null(tpred)) {
    est <- ((pr_r %*% cf$lmr) %*% cf$scl) %*% t(pr_a %*% cf$lma)
  } else {
    est <- NULL
    for(i in 1:fit$ncmp) {
      est <- (pr_r %*% cf$lmr[[i]] %*% cf$scl2[[i]] %*% 
        t(pr_t[, -attr(cf$lmt[[i]],"na.action")] %*% cf$lmt[[i]])) %>%
        as.vector() %>%
        cbind(est, .)
    }
    est <- (est %*% cf$scl1) %*% t(pr_a %*% cf$lma)
    est <- est * fit$tpred_rng[2] / tpred 
  }
  est[is.na(est)] <- 0
  est
}

.pred_sectorial_den <- function(r, a, tpred, fit) {
  .get_d2psf_dxdy(.pred_sectorial_cum, x1 = r, x2 = a, tpred = tpred, 
    fit = fit) / r #/ (pi / 180)
}

.pred_sectorial_psf <- function(r, a, tpred, fit) {
  dpsf <- .pred_sectorial_den(r = r, a = a, tpred = tpred, fit = fit)
  dpsf <- (dpsf[-1,] + dpsf[-nrow(dpsf),]) 
  dpsf <- (dpsf[,-1] + dpsf[,-ncol(dpsf)]) 
  psf  <- (pi / 360) * diff(r^2) * dpsf / 4
  if(r[1] == 0)
    psf[1,] <- diff(.pred_sectorial_cum(r = r[2], a = a, tpred = tpred, 
      fit = fit)[1, ])
  psf
}

# The functions above use matrix notation and require only the coordinates of x
# and the coordinates of y. But that assumes a regular grid. To interpolate the 
# polar grid into a cartesian grid, the radial distances and azimuthal positions
# will not be a regular grid, therefore a vector notation function was written
# below to receive the coordinates {x, y}. This function is called only by the 
# auxiliary function ".sectorial_kernel" for the function "predict_grid".

.pred_sectorial_cum_VEC <- function(r, a, tpred, fit) {
  pr_r <- poly(r, fit$pord, raw = T)
  pr_a <- poly(a, fit$pord, raw = T)
  pr_r <- cbind(rep(1, nrow(pr_r)), pr_r)
  pr_a <- cbind(rep(1, nrow(pr_a)), pr_a)

  if(!is.null(fit$tpred)) {
    # When fitting a third predictor with zero value, the zero is changed to 
    # a small number defined as 1E-6.
    if(tpred == 0) tpred <- 1E-6
    pr_t <- poly(tpred, fit$pord, raw = T)
  } else {
    pr_t <- NULL
  }

  cf  <- fit$coefficients
  if(is.null(tpred)) {
    est <- numeric(length(r))
    for(i in 1:fit$ncmp) {
      est <- est + cf$scl[i,i]  * (pr_r %*% cf$lmr[, i, drop = F]) * 
        (pr_a %*% cf$lma[, i, drop = F])
    }
  } else {
    est <- numeric(length(r))
    for(i in 1:fit$ncmp) {
      for(j in 1:fit$ncmp) {
        est <- est + cf$scl2[[j]][i,i] * (pr_r %*% cf$lmr[[j]][, i, drop = F]) * 
          (pr_t %*% cf$lmt[[j]][, i, drop = F]) * cf$scl[i,i] * 
          (pr_a %*% cf$lma[, i, drop = F])
      }
    }
    est <- est * fit$tpred_rng[2] / tpred 
  }
  est[is.na(est)] <- 0
  est
}


