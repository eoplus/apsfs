#ifndef MC_PSF
#define MC_PSF

SEXP mc_psf(SEXP atm, 
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

#endif
