/*******************************************************************************
  Spatially resolved spherical albedo

  alexandre.castagna@eoplus.science

  Version: 1.0
*******************************************************************************/

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include "struct.h"

// mc_saf Function:

SEXP C_mc_saf(
    SEXP atm, 
    SEXP geom, 
    SEXP res, 
    SEXP ext, 
    SEXP snspos,
    SEXP snsznt, 
    SEXP np, 
    SEXP mnw, 
    SEXP cdf_aer,
    SEXP cdf_ray) {

    // R format variables for easier handling of R lists:
    SEXP km_in, km, tau, tau_u, tau_d, dkm, w0_tot, b_ray, b_tot;
    SEXP c_tot, cdfas, psias, cdfrs, psirs, saf, bin_brks_sxp;
    SEXP bin_phtw_sxp, bin_mid_sxp;
    R_xlen_t atm_n;

    // Variable declaration:
    unsigned long long int n;
    int    geomi, n_brks, i, sid, clid, n_aer, n_ray;
    double exti, resi, sod, tau_tot, km_max, PI2, phi, psi;
    double sint, s, mnwi, ru, bs, b_rat;
    double snsposi[3], cdir[3];
    double *bin_brks, *bin_phtw, *psia, *cdfa, *psir, *cdfr;
    void (*bin_accm)(struct str_phtpkg, double*, double*, size_t);
    struct str_phtpkg pht;
    unsigned long long int NP;

    // Randon sampling variables:
    gsl_rng *random;
    struct timeval tv;
    unsigned long int seed;

    // Initiate random number generator
    gettimeofday (&tv, 0);
    seed   = tv.tv_sec + tv.tv_usec;
    gsl_rng_env_setup();
    random = gsl_rng_alloc (gsl_rng_ranlxs0);
    gsl_rng_set (random, seed);

    // Internal copy of selected R arguments:
    geomi = Rf_asInteger(geom);
    resi  = Rf_asReal(res);
    exti  = Rf_asReal(ext);
    mnwi  = Rf_asReal(mnw);
    NP    = (unsigned long long int) Rf_asReal(np);
    snsposi[0] = REAL(snspos)[0];
    snsposi[1] = REAL(snspos)[1];
    snsposi[2] = REAL(snspos)[2];

    // Get atmospheric profile variables:
    km_in  = Rf_protect(VECTOR_ELT(atm, 0));
    km     = Rf_protect(Rf_duplicate(km_in));
    tau    = Rf_protect(VECTOR_ELT(atm, 1));
    tau_u  = Rf_protect(VECTOR_ELT(atm, 2));
    tau_d  = Rf_protect(VECTOR_ELT(atm, 3));
    dkm    = Rf_protect(VECTOR_ELT(atm, 4));
    w0_tot = Rf_protect(VECTOR_ELT(atm, 6));
    b_ray  = Rf_protect(VECTOR_ELT(atm, 7));
    b_tot  = Rf_protect(VECTOR_ELT(atm, 9));
    c_tot  = Rf_protect(VECTOR_ELT(atm, 10));
    atm_n  = XLENGTH(tau);

    // Get scattering variables:
    psias = Rf_protect(VECTOR_ELT(cdf_aer, 0));
    cdfas = Rf_protect(VECTOR_ELT(cdf_aer, 1));
    n_aer = Rf_xlength(cdfas);
    psia = (double*) calloc(n_aer, sizeof(double));
    cdfa = (double*) calloc(n_aer, sizeof(double));
    for(i = 0; i < n_aer; i++) {
        psia[i] = REAL(psias)[i];
        cdfa[i] = REAL(cdfas)[i];
    }

    psirs = Rf_protect(VECTOR_ELT(cdf_ray, 0));
    cdfrs = Rf_protect(VECTOR_ELT(cdf_ray, 1));
    n_ray = Rf_xlength(cdfrs);
    psir = (double*) calloc(n_ray, sizeof(double));
    cdfr = (double*) calloc(n_ray, sizeof(double));
    for(i = 0; i < n_ray; i++) {
        psir[i] = REAL(psirs)[i];
        cdfr[i] = REAL(cdfrs)[i];
    }

    // Set up Monte Carlo spatial breaks, accumulators and accumulator functions.
    // Accumulators are:
    // (1) bin_phtw - reflected photons.

    if(geomi == 1) { // (1) annular

        n_brks = floor((exti - resi) / resi) + 3;
        bin_phtw = (double*) calloc(n_brks - 1, sizeof(double));
        bin_brks = (double*) calloc(n_brks, sizeof(double));
        bin_brks[0] = 0;
        bin_brks[1] = resi / 2;
        bin_brks[n_brks - 1] = INFINITY;

        for(i = 2; i < n_brks - 1; i++) {
            bin_brks[i] = bin_brks[i - 1] + resi;
        }

        bin_accm = &accm_annular;

    } else if(geomi == 2) { // (2) "sectorial"
    
        n_brks = floor((exti - resi) / resi) + 3;
        bin_phtw = (double*) calloc((n_brks - 1) * 360, sizeof(double));
        bin_brks = (double*) calloc(n_brks, sizeof(double));
        bin_brks[0] = 0;
        bin_brks[1] = resi / 2;
        bin_brks[n_brks - 1] = INFINITY;
    
        for(i = 2; i < n_brks - 1; i++) {
            bin_brks[i] = bin_brks[i - 1] + resi;
        }

        bin_accm = &accm_sectorial;

    } else if(geomi == 3) { // (3) "grid"

        n_brks   = 2 * (floor((exti - resi) / resi) + 2);
        bin_phtw = (double*) calloc((n_brks - 1) * (n_brks - 1), sizeof(double));
        bin_brks = (double*) calloc(n_brks, sizeof(double));
        bin_brks[0] = -INFINITY;
        bin_brks[n_brks / 2] = resi / 2;
        bin_brks[n_brks - 1] = INFINITY;
    
        for(i = (n_brks / 2) + 1; i < (n_brks - 1); i++) {
            bin_brks[i] = bin_brks[i - 1] + resi;
        }
        for(i = (n_brks / 2) - 1; i > 0; i--) {
            bin_brks[i] = bin_brks[i + 1] - resi;
        }

        bin_accm = &accm_grid;

    }

    // Set variables to reduce repetitive calculations:
    tau_tot = REAL(tau)[atm_n - 1];
    PI2     = 2.0 * M_PI;
    km_max  = REAL(km_in)[0];
    for(i = 0; i < atm_n; i++) {
        REAL(km)[i] = km_max - REAL(km)[i];
    }
    
    // For spherical albedo, the source is always at the surface, so the last 
    // layer of the atmosphere. 
    snsposi[2] = km_max - snsposi[2];
    sid     = atm_n - 2;
    sod     = tau_tot;

    for(n = 0; n < NP; n++) {

        // Set the current layer to be the layer that the sensor is in:
        clid = sid;

        // Generate random initial direction within sensor FOV, assuming equal 
        // sensibility to all incoming angles within the FOV. For Spherical 
        // Albedo, the sensor is a downwelling plane irradiance.
        phi = gsl_rng_uniform (random) * PI2;
        psi = acos( -sqrt( gsl_rng_uniform (random) ) );
        
        cdir[0] = sin(psi) * cos(phi);
        cdir[1] = sin(psi) * sin(phi);
        cdir[2] = cos(psi);
        
        // Initiate photon package:
        pht.cpos[0] = snsposi[0];
        pht.cpos[1] = snsposi[1];
        pht.cpos[2] = snsposi[2];
        pht.ppos[0] = snsposi[0];
        pht.ppos[1] = snsposi[1];
        pht.ppos[2] = snsposi[2];
        pht.codz = sod;
        pht.podz = sod;
        pht.sdir[0] = psi;
        pht.sdir[1] = phi;
        pht.cdir[0] = cdir[0];
        pht.cdir[1] = cdir[1];
        pht.cdir[2] = cdir[2];
        pht.stks[0] = 1.0;
        pht.stks[1] = 0.0;
        pht.stks[2] = 0.0;
        pht.stks[3] = 0.0;
        pht.scat = 0;

        while(pht.stks[0] > mnwi) {

            pht.ppos[0] = pht.cpos[0];
            pht.ppos[1] = pht.cpos[1];
            pht.ppos[2] = pht.cpos[2];

            // Sample free optical pathlength to transverse before next scattering.
            ru =  gsl_rng_uniform_pos (random);
            s  = -log(ru);

            // Photon escapes the atmosphere:
            if( (pht.codz + s * pht.cdir[2]) <= 0 ) {
                break;
            }

            // Move photon in metric distances through an optically 
            // heterogeneous layered medium:
            mvpht_m (&pht, s, tau, tau_u, tau_d, c_tot, dkm, clid, (int) atm_n);

            // If photon reached the surface, its weight is summed into the appropriate
            // bin covering the intersection position and photon is terminated.
            if(pht.cpos[2] >= km_max) {

                // Backtrace to surface layer:
                pht.cdir[0] = -pht.cdir[0];
                pht.cdir[1] = -pht.cdir[1];
                pht.cdir[2] = -pht.cdir[2];
                bs = (km_max - pht.cpos[2]) / pht.cdir[2];
                pht.cpos[0] += (pht.cdir[0] * bs);
                pht.cpos[1] += (pht.cdir[1] * bs);
                pht.cpos[2] += (pht.cdir[2] * bs);

                (*bin_accm)(pht, bin_phtw, bin_brks, n_brks);
                break;
            }

            // Update photon properties with a new scattering event:
            clid = findInterv(pht.cpos[2], (double*) REAL(km), (int) atm_n);
            pht.stks[0] *= REAL(w0_tot)[clid];
            b_rat = REAL(b_ray)[clid] / REAL(b_tot)[clid];

            if(gsl_rng_uniform (random) > b_rat) {
                psi = scat_fun (gsl_rng_uniform (random), cdfa, psia, n_aer);
            } else {
                psi = scat_fun (gsl_rng_uniform (random), cdfr, psir, n_ray);
            }

            phi = PI2 * gsl_rng_uniform (random);
            update_cdir (psi, phi, &pht.cdir[0]);
            pht.scat += 1;
            //COS2SPH(pht.cdir, pht.sdir)
        }
    }

    // Normalize the results:

    if(geomi == 1) {
        for(i = 0; i < n_brks - 1; i++) {
            bin_phtw[i] /= (double) NP;
        }
    }
    if(geomi == 2) {
        for(i = 0; i < ((n_brks - 1) * 360); i++) {
            bin_phtw[i] /= (double) NP;
        }
    }
    if(geomi == 3) {
        for(i = 0; i < ((n_brks - 1) * (n_brks - 1)); i++) {
            bin_phtw[i] /= (double) NP;
        }
    }

    bin_brks_sxp = Rf_protect(Rf_allocVector(REALSXP, n_brks));
    for(i = 0; i < n_brks; i++) {
        REAL(bin_brks_sxp)[i] = bin_brks[i];
    }
    bin_mid_sxp = Rf_protect(Rf_allocVector(REALSXP, n_brks - 1));
    for(i = 0; i < n_brks - 1; i++) {
        REAL(bin_mid_sxp)[i] = (bin_brks[i] + bin_brks[i + 1]) / 2.0;
    }
    
    if(geomi == 1) {
        bin_phtw_sxp = Rf_protect(Rf_allocVector(REALSXP, n_brks - 1));
    } else if(geomi == 2) {
        bin_phtw_sxp = Rf_protect(Rf_allocVector(REALSXP, (n_brks - 1) * 360));
    } else if(geomi == 3) {
        bin_phtw_sxp = Rf_protect(Rf_allocVector(REALSXP, (n_brks - 1) * (n_brks - 1)));
    }
    for(i = 0; i < (int) Rf_xlength(bin_phtw_sxp); i++) {
        REAL(bin_phtw_sxp)[i] = bin_phtw[i];
    }

    saf = Rf_protect(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(saf, (R_xlen_t) 0, bin_phtw_sxp);
    SET_VECTOR_ELT(saf, (R_xlen_t) 1, bin_brks_sxp);
    SET_VECTOR_ELT(saf, (R_xlen_t) 2, bin_mid_sxp);

    free(bin_phtw);
    free(bin_brks);

    Rf_unprotect(18);

    return saf;
}
