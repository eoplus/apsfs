
# Calculates the Radial Zernike Polynomials#'
#
# @param n   Radial order, [0, 100].
# @param m   Azimuthal order, [0, 50].
# @param rho Radius, [0, 1].
#
# @useDynLib apsfs C_rzernike

#.rzernike <- function(n, m, rho) {
#  .Call("C_rzernike", n, m, rho)
#}

#' Calculates the Radial Zernike Polynomials up to a given radial order
#'
#' @param rho    Radius, [0, 1].
#' @param phi    Azimuth, [0, 2pi].
#' @param thetav View angle, [0, pi/2].
#' @param order  Maximum radial order.
#'
#' @useDynLib apsfs C_ozernike
#' @export

ozernike <- function(rho, phi, thetav, ord) {

  spec <- .get_spec_comp(ord)
  res  <- .Call("C_ozernike", rho, phi, thetav, spec)
  res
#  ncol <- (order + 1) * (order + 2) / 2
#  orzp <- matrix(NA, ncol = ncol, nrow = length(rho))
#  id <- 1
#  for(n in 0:order) {
#    for(m in (-n):n) {
#      if((n - abs(m)) %% 2 == 1) {
#        next
#      }
#      if(m < 0) {
#        orzp[, id] <- .rzernike(n, m, rho) * sin(m * phi)
#      } else {
#        orzp[, id] <- .rzernike(n, m, rho) * cos(m * phi)
#      }
#      id <- id + 1
#    }
#  }
#  orzp
}

.get_spec_comp <- function(ord) {
# ord is the maximum radial order

  m <- NULL
  for(i in 0:ord) {
    for(j in 0:i) {
      for(p in 0:j) {
        if(((i - p) %% 2) == 0 & ((j - p) %% 2) == 0)
          m <- rbind(m, c(i, j, p))
      }
    }
  }

  m

}
