
#' Monte Carlo Point Spread Function
#'
#' Backward Monte Carlo Code to calculate the elastic scattering Point Spread 
#' Function (integrated in bin area) of a vertically stratified but horizontally 
#' homogeneous atmosphere at monocromatic regime.
#'
#' @param atm     Atmosphere model, see Details. 
#' @param cdf_aer Aerosol scattering CDF LUT (unitless).
#' @param snspos  Sensor position (km; x, y, z), [0,Inf).
#' @param snsfov  Sensor field of view (radians), [0,2*pi].
#' @param snsznt  Sensor view zenith angle (radians), [0,pi].
#' @param geom    One of "annular", "sectorial" or "grid", see Details.
#' @param res     Resolution (km), [0,10]. See Details.
#' @param ext     Maximum spatial extent at geom resolution (km), see Details. 
#' @param np      Number of photons to trace, [0,2^64 - 1].
#' @param mnw     Minimal weight to keep tracing photon package (unitless), 
#'                [0,1].
#' @param press   Atmospheric pressure at surface (mbar). Can be omitted if 
#'                present on the attributes of atm. See Details.
#' @param Rmc     Logical. Should the R version of the MC code be used? Defaults
#'                to FALSE. See Details.
#'
#' @details
#' This backward Monte Carlo radiative transfer code calculates the area 
#' integrals of the Point Spread Function of the atmosphere due to elastic 
#' scattering in a layered atmosphere. Turbulence effects and inelastic 
#' scattering are not included. Additional simplifications of the current 
#  version include scalar solution only (polarization not included), Lambertian 
#' surface reflectance (Bidirectional Reflectance Distribution Function effects 
#' not included). The surface is considered totally absorbing, such that 
#' multiple interaction between surface and atmosphere are not included.
#'
#' The atmosphere profile atm should be a named data frame as the one created by 
#' function \code{get_opt_atm_prfl}. If the profile is create manually or with 
#' user functions, it must contain at least the same variables and column names.
#' It may also include pressure in the attributes (attr(atm, "press") <- value),
#' and in this case the parameter press can be omitted.
#'
#' The initial direction of the photon package is sampled assuming equal 
#' sensibility of the sensor to all directions within its FOV. As an example, 
#' the FOV of the Operational Land Imager (OLI/Landsat 8) is ~ 8.6E-5(radians). 
#' To simulate a sensor with infinitesimal FOV, set snsfov = 0.
#' The FOV is centered in sensor zenith angle. The photon package is initiated 
#' at the upper boundary of the layer where the sensor is located. The sensor 
#' azimuth is set to 90 degrees, such that it is to the right of the center of 
#' the PSF grid.
#'
#' At each iteration, a free optical path is sampled and the photon package is 
#' moved in metric scale. A check is made to verify if photon reached the 
#' surface or if it escaped the system. If it reached the surface, its intensity
#' ("photon package weight") is added to the appropriate accumulator bin. The 
#' accumulator geometry can be: (1) 'annular', which only keep track of the 
#' annular bin in which the photon reached the surface; (2) 'sectorial', wich 
#' divide the annuli in sectors for azimuthal assymetry; or (3) 'grid', which 
#' will keep track of the {x,y} position. The 'annular' accumulator is 
#' appropriate for symmetrical conditions, as nadir view over a Lambertian 
#' surface. Note that regardless of extent of the accumulator, a final bin will 
#' be added between ext and Inf. The resolution of the accumulator will have an
#' impact on the performance, as too high resolution, specially in 'grid' 
#' geometry will require more photons to be traced to reduce statistical 
#' sampling variation. Note that when using 'sectorial' geometry, angular 
#' resolution in azimuth is fixed to 1 degree and break points returned by the 
#' function refer only to the radius break points (angular break points will be
#' 0:360 always).
#'
#' The reduced performance with high resolution grid and the need for high 
#' resolution close to the target pixel while low resolution is sufficient for 
#' larger distances suggests two separate runs. Those grids can be combined
#' later with the function \code{fit_grid_dpsf} to produce a Zernike polynomial
#' fit to the density PSF. This model can then be used to generate a discrete 
#' PSF kernel at any sensor resolution.
#'
#' To reduce the random sampling fluctuation in the sectorial and grid 
#' geometries, the grid is 'folded' on its axis of symmetry (sensor azimuth 
#' line), averaged and then 'unfolded'.
#'
#' Note that the returned values are fractional contribution per accumulation 
#' bin (integrals), and so are not normalized per area. Conversion to per area
#' contribution can be made with function \code{an_psf}.
#'
#' An R version of the code is also provided. It is present mainly for 
#' validation reasons but can be useful in a system without a compiler to
#' compile the C code. It will run many times slower though.
#'
#' @return
#' A list with elements 'dirtw' containing the fractional contribution of the
#' upward directly transmitted photons, 'bin_phtw' containing the fractional 
#' contribution of upward diffusely transmitted photons per accumulation bin
#' (depends on geom), 'bin_brks' containing the break points of the accumulation
#' bins, 'bin_mid' with the mid point of the bins and 'metadata' containing the 
#' input parameters.
#'
#' @examples
#' # Using the C version of the code:
#' data(us76)
#' atm <- rayleigh_od(us76, lambda = 550) %>%
#'        get_opt_atm_prfl(atm = us76, tau_aer = 0.5, H_aer = 2, w0_aer = 0.89, 
#'          tau_ray_z = ., a_mol_z = rep(0, nrow(us76)))
#' cdf_aer <- cdf_lut(continental_ph_6sv[, c(1, 8)])
#' psf_int <- mc_psf(atm     = atm, 
#'                   geom    = 'annular', 
#'                   res     = 0.03, 
#'                   ext     = 10, 
#'                   snspos  = c(0, 0, 800), 
#'                   snsfov  = 0, 
#'                   snsznt  = 0, 
#'                   np      = 1E5, 
#'                   mnw     = 1E-6,
#'                   cdf_aer = cdf_aer)
#' par(mar = c(5, 6, 3, 2))
#' x <- psf_int$bin_brks
#' plot(x, cum_psf(psf_int), xaxs = "i",  yaxs = "i", xlab = "Radius (km)", 
#' ylab = expression(2*pi*integral("r'"*f("r'")*d*"r'", 0, r)), xlim = c(0, 10))
#'
#' @useDynLib apsfs C_mc_psf
#' @export

mc_psf <- function(atm, geom, res, ext, snspos, snsfov, snsznt, np, mnw, 
  cdf_aer, cdf_ray = NULL, press = NULL, Rmc = FALSE) {

  ain <- atm
  atm <- data.frame(km     = atm$km,
                    tau    = atm$tau,
                    tau_u  = c(atm$tau[-length(atm$tau)], NA),
                    tau_d  = c(atm$tau[-1], NA),
                    dkm    = c(abs(diff(atm$km)), NA),
                    dtau   = c(abs(diff(atm$tau)), NA),
                    w0_tot = atm$w0_tot,
                    b_ray  = atm$b_ray,
                    b_aer  = atm$b_aer,
                    b_tot  = atm$b_ray + atm$b_aer,
                    c_tot  = atm$c_tot
  )

  if(is.null(cdf_ray)) {
    cdf_ray <- rayleigh_cdf
  }

  # Check pressure information. It is not necessary on the MC code per se, but 
  # it has to be included on the metadata of the PSF object for further 
  # processing.
  if(is.null(press)) {
    press <- attr(ain, "press")
    if(is.null(press)) {
      "Pressure mst be defined in press is not present as an attribute of atm" %>%
      stop(call. = FALSE)
    }
  }

  if(Rmc) {
    psf <- .Rmc_psf(atm, geom, res, ext, snspos, snsfov, snsznt, np, mnw, 
      cdf_aer)
  } else {
    if(geom == "annular") {
      geom <- 1
    } else if(geom == "sectorial") {
      geom <- 2
    } else if(geom == "grid") {
      geom <- 3
    } else {
      stop("geom must be 'annular', 'sectorial' or 'grid', see Details.", 
        call. = FALSE)
    }

    if(np > (2^64 - 1)) {
      warning("Maximum number of photons to simulate is 2^64 - 1 (18,446,744,073,709,551,615). np set to that limit.", call. = FALSE)
      np <- (2^64 - 1)
    }
    psf <- .Call("C_mc_psf", atm, geom, res, ext, snspos, snsfov, snsznt, 
      np, mnw, cdf_aer, cdf_ray)

    names(psf) <- c("bin_phtw", "dirtw", "bin_brks", "bin_mid")

    # Averaging along axis of symmetry to reduce effect of statistical 
    # fluctuations:
    if(geom == 2) {
      psf$bin_phtw <- matrix(psf$bin_phtw, nrow = length(psf$bin_mid))

      nd  <- ncol(psf$bin_phtw)
      ndh <- nd / 2 

      psf$bin_phtw[, 1:ndh] <- (psf$bin_phtw[, 1:ndh] + 
                                psf$bin_phtw[, nd:(ndh+1)]) / 2
      psf$bin_phtw[, nd:(ndh+1)] <- psf$bin_phtw[, 1:ndh]
    }

    if(geom == 3) {
      psf$bin_phtw <- matrix(psf$bin_phtw, ncol = length(psf$bin_mid))

       # Averaging along axis of symmetry to reduce effect of statistical 
       # fluctuations:
       nd  <- nrow(psf$bin_phtw)
       ndh <- (nd - 1) / 2 

       psf$bin_phtw[1:ndh,] <- (psf$bin_phtw[1:ndh,] + 
         psf$bin_phtw[nd:(ndh+2),]) / 2
       psf$bin_phtw[nd:(ndh+2),] <- psf$bin_phtw[1:ndh,]
    }

    geom <- switch(geom, "annular", "sectorial", "grid")
  }

  psf$metadata <- list(
    press   = press,
    atm     = ain,
    cdf_aer = cdf_aer, 
    snspos  = snspos, 
    snsfov  = snsfov, 
    snsznt  = snsznt, 
    geom    = geom,
    res     = res, 
    ext     = ext, 
    np      = np, 
    mnw     = mnw
  )

  class(psf) <- c("psf", class(psf))

  psf
}

# R version of the MC code. It will be many times slower...

.Rmc_psf <- function(atm, geom, res, ext, snspos, snsfov, snsznt, np, mnw, 
      cdf_aer) {

  if(!is(atm, "data.frame"))
    stop("atm must be a data.frame with suitable column names. See Details.", 
      call. = FALSE)

  if(is.infinite(ext))
    stop(paste0("ext cannot be infinite. An infinite limit is added ",
      "automatically as the final bin break."), call. = FALSE)

  if(snsznt != 0 & geom != 'grid')
    stop(paste0("Results will not be informative for annular accumulation if ",
      "sensor zenith angle is not 0"), call. = FALSE)

  np <- round(np)
  np_blck <- round(np / 10)
  i  <- 1
  pb <- txtProgressBar(min = 0, max = 10, style = 3)

  # Set up Monte Carlo spatial breaks, accumulators and accumulator functions.
  # Accumulators are:
  # (1) dirtw    - directly transmitted photons;
  # (2) bin_phtw - diffuse transmitted photons.
  dirtw      <- 0
  if(geom == "annular") {
    bin_brks <- c(0, seq(res / 2, ext, res), Inf)
    bin_phtw <- numeric(length(bin_brks)-1)
    bin_accm <- .geom_annular
  } else if(geom == "sectorial") {
    bin_brks <- c(0, seq(res / 2, ext, res), Inf)
    bin_phtw <- matrix(0, ncol = 360, nrow = length(bin_brks)-1)
    bin_accm <- .geom_sectorial
  } else if(geom == "grid") {
    bin_brks <- c(seq(res / 2, ext, res), Inf)
    bin_brks <- c(-rev(bin_brks), bin_brks)
    bin_phtw <- matrix(0, ncol = length(bin_brks)-1, nrow = length(bin_brks)-1)
    bin_accm <- .geom_grid_regular
  } else {
    stop("geom must be 'annular' or 'grid', see Details.", call. = FALSE)
  }

  # Set variables to reduce repetitive calculations:
  sod     <- approx(x = atm$km, y = atm$tau, xout = snspos[3])$y
  sid     <- findInterval(x = sod, vec = atm$tau, left.open = TRUE, 
               all.inside = TRUE)
  tau_tot <- max(atm$tau)
  snshfov <- snsfov / 2
  pi2     <- 2 * pi
  km_max  <- max(atm$km)
  atm$km  <- km_max - atm$km
  snspos[2] <- tan(snsznt) * snspos[3]
  snspos[3] <- km_max - snspos[3]
  phi_s <- 270 * pi / 180
  snscdir <- c(
    sin(snsznt) * cos(phi_s), 
    sin(snsznt) * sin(phi_s), 
    cos(snsznt)
  )

  for(n in 1:np) {

    # Find the layer that the sensor is in:
    clid <- sid

    # Generate random initial direction within sensor FOV, assuming equal 
    # sensibility to all incoming angles within the FOV:
    if(snsfov == 0) {
      phi  <- phi_s
      psi  <- snsznt
      cdir <- snscdir
    } else {
      phi <- runif(1) * pi2
      psi <- runif(1) * snshfov
      if(snsznt != 0) {
        cdir <- update_cdir (psi, phi, snscdir)
        sdir <- cos2sph (cdir)
        psi  <- sdir[1]
        phi  <- sdir[2]
      }
    }

    # Initiate photon package:
    # cpos - Curent position (m), (x, y, z);
    # ppos - Previous position (m), (x, y, z);
    # codz - Current optical depth position (unitless); 
    # podz - Previous optical depth position (unitless);
    # sdir - Spherical direction (radians), (psi, phi);
    # cdir - Directional cosines (unitless), (dx/dl, dy/dl, dz/dl);
    # stks - Stokes vector (unitless), (I, Q, U, V);
    # scat - Count of scattering events.
    pht <- list(
      cpos = snspos,
      ppos = snspos,
      codz = sod,
      podz = sod,
      sdir = c(psi, phi),
      cdir = cdir,
      stks = c(1, 0, 0, 0),
      scat = 0
    )

    while(pht$stks[1] > mnw) {
      pht$ppos <- pht$cpos

      # Sample free optical pathlength to transverse before next scattering.
      # If ru == 0, s is Inf and that will cause error surface intersection. A 
      # value of 1E-9 will result in s = 20 which is enough for Earth atmospheres.
      ru       <- runif(1)
      if(ru == 0) ru <- 1E-9
      s        <- -log(ru)

      # If photon was not scattered, it is part of direct transmission:
      if(pht$scat == 0 & (pht$codz + s * pht$cdir[3]) >= tau_tot) {                        
        dirtw <- dirtw + pht$stks[1]
        break
      }

      # Move photon in metric distances through an optically inhomogeneous 
      # layered medium:
      pht <- .movpht_m(pht, s, atm, clid)

      # If photon reached the surface, its weight is summed into the appropriate
      # bin covering the intersection position and photon is terminated.
      if(pht$cpos[3] >= km_max) {

        # Backtrace to surface layer:
        pht$cdir <- -pht$cdir
        bs       <- (km_max - pht$cpos[3]) / pht$cdir[3]
        pht$cpos <-  pht$cpos + (pht$cdir * bs)

        bin_phtw <- bin_accm(pht, bin_phtw, bin_brks, res)
        break
      }

      # Photon scapes the atmosphere:
      if(pht$cpos[3] < 0) {
        break
      }

      # Update photon properties with a new scattering event:
      clid <- findInterval(x = pht$cpos[3], vec = atm$km, 
        left.open = TRUE, all.inside = TRUE)
      pht$stks[1] <- pht$stks[1] * atm[clid, ]$w0_tot
      b_rat <- (atm[clid, ]$b_ray / atm[clid, ]$b_tot)
      if(is.na(b_rat)) b_rat <- 0.5                                             # Sometimes b_rat is NA, probably due to division by 0 in the top atmosphere layer.
      if(runif(1) > b_rat) {
        psi = scat_aer (runif(1), cdf_aer)
      } else {
        psi = scat_ray (runif(1))
      }
      phi <- pi2 * runif(1)
      pht$cdir <- update_cdir (psi, phi, pht$cdir)                              
      pht$scat <- pht$scat + 1
    }
    
    if((n %% np_blck) == 0) {
     
      setTxtProgressBar(pb, i)
      i <- i + 1

    }
  }

  # Normalize the results:
  dirtw    <- dirtw / np
  bin_phtw <- bin_phtw / np

  # Fold the sectors and grids to average and reduce noise:
  if(geom == 'sectorial') {
   nd  <- ncol(bin_phtw)
   ndh <- nd / 2 

   bin_phtw[, 1:ndh] <- (bin_phtw[, 1:ndh] + bin_phtw[, nd:(ndh+2)]) / 2
   bin_phtw[, nd:(ndh+2)] <- bin_phtw[, 1:ndh]
  }

  if(geom == 'grid') {
   nd  <- nrow(bin_phtw)
   ndh <- (nd - 1) / 2 

   bin_phtw[1:ndh,] <- (bin_phtw[1:ndh,] + bin_phtw[nd:(ndh+2),]) / 2
   bin_phtw[nd:(ndh+2),] <- bin_phtw[1:ndh,]
  }

  atm$km  <- km_max - atm$km
  psf <- list(
    bin_phtw  = bin_phtw,
    dirtw     = dirtw,
    bin_brks  = bin_brks,
    bin_mid   = (bin_brks[-1] + bin_brks[-length(bin_brks)]) / 2
  )

  setTxtProgressBar(pb, 10)
  close(pb)

  psf
}

#' Ancillary function to move photon in stratified atmosphere:

.movpht_m <- function(pht, s, atm, cid) {

  pid      <- cid
  pht$podz <- pht$codz
  pht$codz <- pht$podz + (s * pht$cdir[3])
  cid      <- findInterval(pht$codz, atm$tau, left.open = TRUE, 
    all.inside = TRUE)

  if(pid == cid) {

    lkm      <- s / atm$c_tot[cid]
    pht$cpos <- pht$cpos + (lkm * pht$cdir)

  } else {

    lkm <- 0
    dir <- ifelse(pht$cdir[3] > 0, "tau_d", "tau_u")
    lkm <- lkm + (atm[pid, dir] - pht$podz) / pht$cdir[3] / atm$c_tot[pid]
    dir <- ifelse(pht$cdir[3] > 0, "tau_u", "tau_d")
    lkm <- lkm + (pht$codz - atm[cid, dir]) / pht$cdir[3] / atm$c_tot[cid]

    if(abs(cid - pid) > 1) {

      idv <- setdiff(cid:pid, c(cid, pid))
      lkm <- lkm + abs(sum(atm$dkm[idv]) / pht$cdir[3])

    }

    pht$cpos <- pht$cpos + (lkm * pht$cdir)

  }

  return(pht)

}

#' Ancillary functions for accumulators:

.geom_annular <- function(pht, bin_phtw, bin_brks, res) {

  rpht <- sqrt(pht$cpos[1]^2 + pht$cpos[2]^2)

#  lng <- length(bin_brks) - 1
#  id  <- ceiling(rpht / res)
#  id[id < 1]   <- 1
#  id[id > lng] <- lng

  id   <- findInterval(x = rpht, vec = bin_brks, left.open = TRUE, 
    all.inside = TRUE)
  bin_phtw[id] <- bin_phtw[id] + pht$stks[1]

  return(bin_phtw)

}

.geom_sectorial <- function(pht, bin_phtw, bin_brks, res) {

  rpht <- sqrt(pht$cpos[1]^2 + pht$cpos[2]^2)
  azmt <- pht$cpos[2] / rpht
  if(pht$cpos[2] > 0) azmt <- (2 * pi) - azmt

  rid  <- findInterval(x = rpht, vec = bin_brks, left.open = TRUE, 
    all.inside = TRUE)
  aid  <- floor(azmt * 180 / pi) + 1

  bin_phtw[rid[1], aid[2]] <- bin_phtw[rid[1], aid[2]] + pht$stks[1]

  return(bin_phtw)

}

.geom_grid_regular <- function(pht, bin_phtw, bin_brks, res) {

  lng  <- length(bin_brks) - 1
  xyid <- ceiling((pht$cpos[1:2] - (bin_brks[2] - res)) / res)
  xyid[xyid < 1]   <- 1
  xyid[xyid > lng] <- lng

  bin_phtw[xyid[1], xyid[2]] <- bin_phtw[xyid[1], xyid[2]] + pht$stks[1]

  return(bin_phtw)

}

