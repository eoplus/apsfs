#' Convert form directional cosines to spherical directions
#'
#' Calculates the polar spherical directions from directional cosines.
#'
#' @param cdir A vector containing the directional cosines (dx/dl, dy/dl, dz/dl).
#'
#' @return A vector of length two containing the polar angle (psi) and the 
#' azimuthal angle (phi).
#'
#' @export

cos2sph <- function(cdir) {

    sdir    <- numeric(2)
    sdir[1] <- acos(cdir[3])

    if((abs(cdir[3]) - 1) < 1E-12) {
        sdir[2] <- 2 * pi * runif(1)
    } else {
        sint    <- sin(sdir[1])
        sdir[2] <- acos(round((cdir[1] / sint) * 1.0E012) / 1.0E012)
        if(cdir[2] < 0) {
            sdir[2] <- (2 * pi) - sdir[2]
        }
    }

    return(c(psi = sdir[1], phi = sdir[2]))
}
