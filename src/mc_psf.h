#ifndef MC_PSF
#define MC_PSF

SEXP C_mc_psf(SEXP atm, 
              SEXP geom, 
              SEXP res, 
              SEXP ext, 
              SEXP snspos,
              SEXP snsfov, 
              SEXP snsznt, 
              SEXP np, 
              SEXP mnw, 
              SEXP cdf_aer,
              SEXP cdf_ray);
              
SEXP C_mc_saf(SEXP atm, 
              SEXP geom, 
              SEXP res, 
              SEXP ext, 
              SEXP snspos,
              SEXP snsznt, 
              SEXP np, 
              SEXP mnw, 
              SEXP cdf_aer,
              SEXP cdf_ray);
              
#endif
