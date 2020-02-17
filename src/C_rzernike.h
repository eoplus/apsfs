#ifndef ZERNIKE
#define ZERNIKE

SEXP C_zrnk_s1(SEXP rho, SEXP phi, SEXP spec);
SEXP C_zrnk_s2(SEXP rho, SEXP phi, SEXP thetas, SEXP spec);
SEXP C_zrnk_h2(SEXP rho, SEXP phi, SEXP thetas, SEXP spec);

#endif
