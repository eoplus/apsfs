#' Calculates the mixed partial derivative:  d^2 * z / dx * dy for a scalar 
#' valued function (z) of two arguments (x, y). This is the bidimensional 
#' density function, the underlying process integrated by Monte Carlo, 
#' sub-sequentially used to calculate the bidimensional cumulative 
#' distribution function fitted by fit_sectorial or fit_grid.
#' 
#' The code is modified from numDeriv:::genD. The original function is generic 
#' and vectorized in the sense of a vector valued function. It was adapted here 
#' to solve a scalar valued function at several coordinates at the same time for 
#' reasons of computational efficiency. The modifications below make this 
#' function specific to the problem at hand, in particular the way the 
#' predict_sectorial subfunctions are expected to work.
#' 
#' It is based on the numDeriv R package source code version 2016.8-1.1, that is 
#' licensed under GPL-2. This package is also licensed as GPL-2 and a copy of the
#' license is provided in the root directory of the package. Modifications made 
#' on 2020-03-14.

.get_d2psf_dxdy <- function(func, x1, x2, ...) {
   args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
             r=4, v=2) # default
   d <- args$d
   r <- args$r
   v <- args$v
   if (v!=2) stop("The current code assumes v is 2 (the default).")	 

   f0 <- func(x1, x2, ...)  #f0 is the value of the function at x.
   f0 <- as.vector(f0)

   n <- 2  #  number of parameters (theta) FIXED for this application
   h0 <- list(
     x1 = abs(d*x1) + args$eps * (abs(x1) < args$zero.tol),
     x2 = abs(d*x2) + args$eps * (abs(x2) < args$zero.tol)
   )
   D <- matrix(0, length(f0), (n*(n + 3))/2)
   # length(f0) is the dim of the sample space
   # (n*(n + 3))/2 is the number of columns of matrix D.( first
   #	der. & lower triangle of Hessian)
   Daprox <-  matrix(0,length(f0),r) 
   Hdiag  <-  matrix(0,length(f0),n)
   Haprox <-  matrix(0,length(f0),r)

   for(i in 1:n) {   # each parameter  - first deriv. & hessian diagonal
     h <- h0
     for(k in 1:r) { # successively reduce h 
       f1 <- as.vector(func(x1+(i==1)*h[[1]], x2+(i==2)*h[[2]], ...))
       f2 <- as.vector(func(x1-(i==1)*h[[1]], x2-(i==2)*h[[2]], ...))
       if(i == 1) {
         norm <- rep(h[[i]], length(x2))
       } else {
         norm <- rep(h[[i]], each = length(x1))
       }
       Daprox[,k] <- (f1 - f2)  / (2 * norm)     # F'(i) 
       Haprox[,k] <- (f1 - 2 * f0 + f2)/ norm^2  # F''(i,i) hessian diagonal
       h <- lapply(h, function(x) { x / v })     # Reduced h by 1/v.
     }
     for(m in 1:(r - 1)) {
       for(k in 1:(r-m)) {
         Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[,k]) / (4^m - 1)
         Haprox[, k] <- (Haprox[, k + 1] * (4^m) - Haprox[,k]) / (4^m - 1)
       }
       D[, i]     <- Daprox[,1]
       Hdiag[, i] <- Haprox[,1]
     }
   }

   u <- n
   for(i in 1:n) {   # 2nd derivative  - do lower half of hessian only
     for(j in 1:i) {
       u <- u + 1
       if(i==j) {
         D[,u] <- Hdiag[,i]
       } else {
         h <- h0
         for(k in 1:r){  # successively reduce h 

           if(i == 1) {
             norm1 <- rep(h[[i]], length(x2))
           } else {
             norm1 <- rep(h[[i]], each = length(x1))
           }
           if(j == 1) {
             norm2 <- rep(h[[j]], length(x2))
           } else {
             norm2 <- rep(h[[j]], each = length(x1))
           }

           f1 <- as.vector(
             func(x1+(i==1)*h[[1]] + (j==1)*h[[1]], 
                  x2+(i==2)*h[[2]] + (j==2)*h[[2]], ...)
           )
           f2 <- as.vector(
             func(x1-(i==1)*h[[1]] - (j==1)*h[[1]], 
                  x2-(i==2)*h[[2]] - (j==2)*h[[2]], ...)
           )

           Daprox[,k] <- (f1 - 2 * f0 + f2 - Hdiag[,i] * norm1^2 - Hdiag[,j] * 
             norm2^2) / (2*norm1*norm2)  # F''(i,j)  
           h <- lapply(h, function(x) { x / v })     # Reduced h by 1/v.
         }
         for(m in 1:(r - 1))
           for ( k in 1:(r-m))
             Daprox[,k] <- (Daprox[, k + 1] * (4^m) - Daprox[, k]) / (4^m - 1)
             D[,u] <- Daprox[, 1]
           }
         }  
      }
  matrix(D[, 4], ncol = length(x2)) # d2Z / dX dY
}

.get_d2psf_dxdy_VEC <- function(func, x1, x2, ...) {
   args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
             r=4, v=2) # default
   d <- args$d
   r <- args$r
   v <- args$v
   if (v!=2) stop("The current code assumes v is 2 (the default).")	 

   f0 <- func(x1, x2, ...)  #f0 is the value of the function at x.
   f0 <- as.vector(f0)

   n <- 2  #  number of parameters (theta) FIXED for this application
   h0 <- list(
     x1 = abs(d*x1) + args$eps * (abs(x1) < args$zero.tol),
     x2 = abs(d*x2) + args$eps * (abs(x2) < args$zero.tol)
   )
   D <- matrix(0, length(f0), (n*(n + 3))/2)
   # length(f0) is the dim of the sample space
   # (n*(n + 3))/2 is the number of columns of matrix D.( first
   #	der. & lower triangle of Hessian)
   Daprox <-  matrix(0,length(f0),r) 
   Hdiag  <-  matrix(0,length(f0),n)
   Haprox <-  matrix(0,length(f0),r)

   for(i in 1:n) {   # each parameter  - first deriv. & hessian diagonal
     h <- h0
     for(k in 1:r) { # successively reduce h 
       f1 <- as.vector(func(x1+(i==1)*h[[1]], x2+(i==2)*h[[2]], ...))
       f2 <- as.vector(func(x1-(i==1)*h[[1]], x2-(i==2)*h[[2]], ...))
       norm <- h[[i]]
       Daprox[,k] <- (f1 - f2)  / (2 * norm)     # F'(i) 
       Haprox[,k] <- (f1 - 2 * f0 + f2)/ norm^2  # F''(i,i) hessian diagonal
       h <- lapply(h, function(x) { x / v })     # Reduced h by 1/v.
     }
     for(m in 1:(r - 1)) {
       for(k in 1:(r-m)) {
         Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[,k]) / (4^m - 1)
         Haprox[, k] <- (Haprox[, k + 1] * (4^m) - Haprox[,k]) / (4^m - 1)
       }
       D[, i]     <- Daprox[,1]
       Hdiag[, i] <- Haprox[,1]
     }
   }

   u <- n
   for(i in 1:n) {   # 2nd derivative  - do lower half of hessian only
     for(j in 1:i) {
       u <- u + 1
       if(i==j) {
         D[,u] <- Hdiag[,i]
       } else {
         h <- h0
         for(k in 1:r){  # successively reduce h 
           norm1 <- h[[i]]
           norm2 <- h[[j]]

           f1 <- as.vector(
             func(x1+(i==1)*h[[1]] + (j==1)*h[[1]], 
                  x2+(i==2)*h[[2]] + (j==2)*h[[2]], ...)
           )
           f2 <- as.vector(
             func(x1-(i==1)*h[[1]] - (j==1)*h[[1]], 
                  x2-(i==2)*h[[2]] - (j==2)*h[[2]], ...)
           )

           Daprox[,k] <- (f1 - 2 * f0 + f2 - Hdiag[,i] * norm1^2 - Hdiag[,j] * 
             norm2^2) / (2*norm1*norm2)  # F''(i,j)  
           h <- lapply(h, function(x) { x / v })     # Reduced h by 1/v.
         }
         for(m in 1:(r - 1))
           for ( k in 1:(r-m))
             Daprox[,k] <- (Daprox[, k + 1] * (4^m) - Daprox[, k]) / (4^m - 1)
             D[,u] <- Daprox[, 1]
           }
         }  
      }
  D[, 4] # d2Z / dX dY
}

# The function above are based on an input of matrix notation and require only 
# the coordinates of x and the coordinates of y. Those are of different lengths
# and all coordinates are combinations of x and y. This assumes regular grid. To 
# interpolate the polar grid into a cartesian grid, the radial distances and 
# azimuthal positions will not be a regular grid, therefore a vector notation 
# function was written below to receive the coordinates {x, y}, were x and y 
# have the same length. This function is called only by the auxiliary function 
# ".sectorial_kernel" for the function "predict_grid".

.get_d2psf_dxdy_VEC <- function(func, x1, x2, ...) {
   args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
             r=4, v=2) # default
   d <- args$d
   r <- args$r
   v <- args$v
   if (v!=2) stop("The current code assumes v is 2 (the default).")	 

   f0 <- func(x1, x2, ...)  #f0 is the value of the function at x.
   f0 <- as.vector(f0)

   n <- 2  #  number of parameters (theta) FIXED for this application
   h0 <- list(
     x1 = abs(d*x1) + args$eps * (abs(x1) < args$zero.tol),
     x2 = abs(d*x2) + args$eps * (abs(x2) < args$zero.tol)
   )
   D <- matrix(0, length(f0), (n*(n + 3))/2)
   # length(f0) is the dim of the sample space
   # (n*(n + 3))/2 is the number of columns of matrix D.( first
   #	der. & lower triangle of Hessian)
   Daprox <-  matrix(0,length(f0),r) 
   Hdiag  <-  matrix(0,length(f0),n)
   Haprox <-  matrix(0,length(f0),r)

   for(i in 1:n) {   # each parameter  - first deriv. & hessian diagonal
     h <- h0
     for(k in 1:r) { # successively reduce h 
       f1 <- as.vector(func(x1+(i==1)*h[[1]], x2+(i==2)*h[[2]], ...))
       f2 <- as.vector(func(x1-(i==1)*h[[1]], x2-(i==2)*h[[2]], ...))
       norm <- h[[i]]
       Daprox[,k] <- (f1 - f2)  / (2 * norm)     # F'(i) 
       Haprox[,k] <- (f1 - 2 * f0 + f2)/ norm^2  # F''(i,i) hessian diagonal
       h <- lapply(h, function(x) { x / v })     # Reduced h by 1/v.
     }
     for(m in 1:(r - 1)) {
       for(k in 1:(r-m)) {
         Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[,k]) / (4^m - 1)
         Haprox[, k] <- (Haprox[, k + 1] * (4^m) - Haprox[,k]) / (4^m - 1)
       }
       D[, i]     <- Daprox[,1]
       Hdiag[, i] <- Haprox[,1]
     }
   }

   u <- n
   for(i in 1:n) {   # 2nd derivative  - do lower half of hessian only
     for(j in 1:i) {
       u <- u + 1
       if(i==j) {
         D[,u] <- Hdiag[,i]
       } else {
         h <- h0
         for(k in 1:r){  # successively reduce h 
           norm1 <- h[[i]]
           norm2 <- h[[j]]

           f1 <- as.vector(
             func(x1+(i==1)*h[[1]] + (j==1)*h[[1]], 
                  x2+(i==2)*h[[2]] + (j==2)*h[[2]], ...)
           )
           f2 <- as.vector(
             func(x1-(i==1)*h[[1]] - (j==1)*h[[1]], 
                  x2-(i==2)*h[[2]] - (j==2)*h[[2]], ...)
           )

           Daprox[,k] <- (f1 - 2 * f0 + f2 - Hdiag[,i] * norm1^2 - Hdiag[,j] * 
             norm2^2) / (2*norm1*norm2)  # F''(i,j)  
           h <- lapply(h, function(x) { x / v })     # Reduced h by 1/v.
         }
         for(m in 1:(r - 1))
           for ( k in 1:(r-m))
             Daprox[,k] <- (Daprox[, k + 1] * (4^m) - Daprox[, k]) / (4^m - 1)
             D[,u] <- Daprox[, 1]
           }
         }  
      }
  D[, 4] # d2Z / dX dY
}



