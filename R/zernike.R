
#' Calculate the Zernike polynomials
#'
#' This function calculates the Zernike polynomials for the unit disk, or for 
#' the product of the two unit disks possibly respecting Helmholtz symetry 
#' relations.
#'
#' @param ord    Maximum radial order.
#' @param rho    Radius, [0, 1].
#' @param phi    Azimuth (radians), [0, 2*pi] or [-pi, pi]
#' @param thetav View angle (radians), [0, pi/2]
#' @param sym    Logical. For Zernike products, if it should respect Helmholtz 
#'               symetry.
#'
#' @details If thetav is missing, the function will return the Zernike 
#' terms up with radial order up to ord and all azimuthal orders. If thetav is 
#' present, the same caluclation will be performed for the product of two 
#' Zernike polynomials for isotropic conditions, i.e., azimuthal dependency is 
#' given only by the relative azimuth between pixel center and view vector 
#' projection on the surface. If sym is true, a more complex model is returned
#' that respects Helmholtz symetry and is defined as:
#' N * [Rkm(rhoi) * Rkn(rhos) + Rkm(rhoi) * Rkn(rhor)] * cos(k * Delta_phi)
#' where N is the normalization factor and R are the radial polynomials of 
#' azimuthal order k and radial orders m and n.
#'
#' Zernike calculations implemented in C. Vectorized use is encouraged, but rho 
#' and phi should have the same length.
#'
#' @return A matrix with Zernike terms per row and positions per column.
#'
#' @export
#' @useDynLib apsfs C_zrnk_s1 C_zrnk_s2 C_zrnk_h2

zrnk_exp <- function(ord, rho, phi, thetav, sym = FALSE) {

  if(missing(thetav)) {
    spec <- .get_spec_comp(ord)
    res  <- .Call("C_zrnk_s1", rho, phi, spec)
  } else {
    spec <- .get_spec_comp(ord, thetav, sym)
print(spec)
    if(!sym) res  <- .Call("C_zrnk_s2", rho, phi, thetav, spec)
    else     res  <- .Call("C_zrnk_h2", rho, phi, thetav, spec)
  }

  res

}

.get_spec_comp <- function(ord, thetav, sym = FALSE) {

  spec <- NULL

  if(missing(thetav)) {

    for(n in 0:ord) {
      for(m in -n:n) {
        if(((n - m) %% 2) == 0) {
          spec <- rbind(spec, c(n, m))
        }
      }
    }
    colnames(spec) <- c("n", "m")

  } else {

    for(n in 0:ord) {
      for(m in 0:n) {
        for(k in -m:m) {
          if(((n - k) %% 2) == 0 & ((m - k) %% 2) == 0) {
            if(k < 0 & sym == 1) {
              next
            }
            spec <- rbind(spec, c(n, m, k))
          }
        }
      }
    }
    colnames(spec) <- c("n", "m", "k")
  }

  spec

}

