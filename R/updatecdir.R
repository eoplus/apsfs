#' Update directional cosines
#'
#' Rotates the reference frame of scattering to the reference frame of the 
#' photon package previous direction.
#'
#' @param psi  Polar scattering angle (radians), [0,pi].
#' @param phi  Azimuthal scattering angle (radians), [0,2*pi].
#' @param cdir Directional cosines (unitless; dx/dl, dy/dl, dz/dl).
#'
#' @return A vector of length 3 containing the new travelling directional 
#' cosines (dx/dl, dy/dl, dz/dl) of the photon package.
#'
#' @export

update_cdir <- function(psi, phi, cdir) {

  mus  <- cos(psi)
  sint <- sin(acos(cdir[3]))
  sins <- sin(psi)
  scatmb <- c(0, 0, 0)
  scatma <- matrix(0, ncol = 3, nrow = 3)

  scatmb[1] <- sins * cos(phi)
  scatmb[2] <- sins * sin(phi)
  scatmb[3] <- mus

  if(abs(sint) > 1.0E-012) {
 
    scatma[1, 1] <-  cdir[1] * cdir[3] / sint
    scatma[1, 2] <- -cdir[2] / sint
    scatma[1, 3] <-  cdir[1]
    scatma[2, 1] <-  cdir[2] * cdir[3] / sint
    scatma[2, 2] <-  cdir[1] / sint
    scatma[2, 3] <-  cdir[2]
    scatma[3, 1] <- -sint
    scatma[3, 2] <-  0.0
    scatma[3, 3] <-  cdir[3]

    cdir[1] <- (scatma[1, 1] * scatmb[1]) + (scatma[1, 2] * scatmb[2]) + 
      (scatma[1, 3] * scatmb[3])
    cdir[2] <- (scatma[2, 1] * scatmb[1]) + (scatma[2, 2] * scatmb[2]) + 
      (scatma[2, 3] * scatmb[3])
    cdir[3] <- (scatma[3, 1] * scatmb[1]) + (scatma[3, 2] * scatmb[2]) + 
      (scatma[3, 3] * scatmb[3])

  } else {

    cdir <- sign(cdir[3]) * scatmb

  }
 
  return(cdir)

}
