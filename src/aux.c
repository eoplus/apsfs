/*******************************************************************************
  Auxiliary functions

  alexandre.castagna@eoplus.science

  Version: 1.0
*******************************************************************************/

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include "struct.h"

// Function to update the directional cosines after photon scattering
//
// psi = polar scattering angle
// phi = azimuthal scattering angle
// *cdir = pointer to directional cosines

void update_cdir (
    double psi, 
    double phi, 
    double *cdir) {

    double mus = cos(psi);
    double sint = sin(acos(cdir[2]));
    double sins = sin(psi);
    double scatmb[3] = {0.0};
    double scatma[3][3] = {{0.0}, {0.0}, {0.0}};

    scatmb[0] = sins * cos(phi);
    scatmb[1] = sins * sin(phi);
    scatmb[2] = mus;

    if(ABS(sint) > 1.0E-012) {
        scatma[0][0] =  cdir[0] * cdir[2] / sint; 
        scatma[0][1] = -cdir[1] / sint; 
        scatma[0][2] =  cdir[0];
        scatma[1][0] =  cdir[1] * cdir[2] / sint;
        scatma[1][1] =  cdir[0] / sint;
        scatma[1][2] =  cdir[1];
        scatma[2][0] = -sint;
        scatma[2][1] =  0.0;
        scatma[2][2] =  cdir[2]; 
        cdir[0] = (scatma[0][0] * scatmb[0]) + (scatma[0][1] * scatmb[1]) + (scatma[0][2] * scatmb[2]); 
        cdir[1] = (scatma[1][0] * scatmb[0]) + (scatma[1][1] * scatmb[1]) + (scatma[1][2] * scatmb[2]);
        cdir[2] = (scatma[2][0] * scatmb[0]) + (scatma[2][1] * scatmb[1]) + (scatma[2][2] * scatmb[2]); 
    } else {
        cdir[0] = SIGN(cdir[2]) * scatmb[0];
        cdir[1] = SIGN(cdir[2]) * scatmb[1];
        cdir[2] = SIGN(cdir[2]) * scatmb[2];
    }
}

// Function to move photon in metric scale through an vertically stratified 
// medium.
//
// pht = photon package
// s   = free path in optical depths
// atm = atmospheric profile R data frame
// cid = last recorded medium layer position index

void mvpht_m (
    struct str_phtpkg* pht, 
    double s, 
    SEXP tau,
    SEXP tau_u,
    SEXP tau_d,
    SEXP c_tot,
    SEXP dkm, 
    int cid,
    R_xlen_t atm_n) {

    int    pid, st, ed, i;
    double lkm = 0;

    pid = cid;

    pht->podz = pht->codz;
    pht->codz = pht->podz + (s * pht->cdir[2]);

    cid = findInterv(pht->codz, (double*) REAL(tau), (int) atm_n - 1);

    if(pid == cid) {
        lkm = s / REAL(c_tot)[cid];
    } else {
        if(pht->cdir[2] > 0) {
            lkm += (REAL(tau_d)[pid] - pht->podz) / pht->cdir[2] / REAL(c_tot)[pid];
            lkm += (pht->codz - REAL(tau_u)[cid]) / pht->cdir[2] / REAL(c_tot)[cid];
        } else {
            lkm += (REAL(tau_u)[pid] - pht->podz) / pht->cdir[2] / REAL(c_tot)[pid];
            lkm += (pht->codz - REAL(tau_d)[cid]) / pht->cdir[2] / REAL(c_tot)[cid];
        }

        if(ABS(cid - pid) > 1) {
            if(cid < pid) {
                st = cid + 1;
                ed = pid;
            } else {
                st = pid + 1;
                ed = cid;
            }

            for(i = st; i < ed; i++) {
                lkm += ABS(REAL(dkm)[i] / pht->cdir[2]);
            }
        }
    }

    pht->cpos[0] += lkm * pht->cdir[0];
    pht->cpos[1] += lkm * pht->cdir[1];
    pht->cpos[2] += lkm * pht->cdir[2];
}

// Function to find the interval index where a value lies. It is open at the 
// left and closed in the right. Since brks will always be passed with extremes
// being -INFINITE and INFINITE, the value will always be inside the breaks.
//
// val  = scalar value
// brks = vector of breaks
// n    = length of brks

int findInterv (
    double val, 
    double* brks, 
    size_t n) {

    int i;

    for(i = 0; i < n; i++) {
        if(((val - brks[i+1]) * (val - brks[i])) <= 0) {
            return i;
        }
    }

    if(val < 0) {
        return 0;
    } else {
        return n - 1;
    }
}

// Accumulator function for annular geometry
//
// pht      = photon package
// bin_phtw = photon accumulator
// bin_brks = vector of breaks
// n        = length of bin_brks

void accm_annular (
    struct str_phtpkg pht, 
    double* bin_phtw, 
    double* bin_brks,
    size_t n) {

    int    id;
    double rpht = 0;

    rpht = sqrt(pow(pht.cpos[0], 2.0) + pow(pht.cpos[1], 2.0));
    id = findInterv(rpht, bin_brks, n);
    bin_phtw[id] += pht.stks[0];

}

// Accumulator function for sectorial geometry
//
// pht      = photon package
// bin_phtw = photon accumulator
// bin_brks = vector of breaks
// n        = length of bin_brks

void accm_sectorial (
    struct str_phtpkg pht, 
    double* bin_phtw, 
    double* bin_brks, 
    size_t n) {

    int id, rid, aid;
    double rpht = 0;
    double azmt = 0;

    rpht = sqrt(pow(pht.cpos[0], 2.0) + pow(pht.cpos[1], 2.0));
    if(rpht < 1E-12) {
        azmt = 2.0 * M_PI;
    } else {
        azmt = acos(pht.cpos[1] / rpht);
        if(pht.cpos[0] > 0.0) {
            azmt = (2.0 * M_PI) - azmt;
        }
    }

    rid = findInterv(rpht, bin_brks, n);
    aid = (int) floor(azmt * 180.0 / M_PI);
    id  = (n - 1) * aid + rid;
    bin_phtw[id] += pht.stks[0];

}

// Accumulator function for grid geometry
//
// pht      = photon package
// bin_phtw = photon accumulator
// bin_brks = vector of breaks
// n        = length of bin_brks

void accm_grid (
    struct str_phtpkg pht, 
    double* bin_phtw, 
    double* bin_brks, 
    size_t n) {

    int id, rid, cid;

    rid = findInterv(pht.cpos[0], bin_brks, n);
    cid = findInterv(pht.cpos[1], bin_brks, n);
    id  = (n - 1) * cid + rid;
    bin_phtw[id] += pht.stks[0];

}


double scat_fun (
    double ru, 
    double* cdf, 
    double *psi, 
    int npsi) {

    int id = findInterv(ru, cdf, npsi);
 
    return psi[id];
}

// Function to print out photon package status, used for development purposes

void print_pht (
    struct str_phtpkg pht) {
    printf("Photon package:\n");
    printf("    PPOS: %f, %f, %f\n", pht.ppos[0], pht.ppos[1], pht.ppos[2]);
    printf("    CPOS: %f, %f, %f\n", pht.cpos[0], pht.cpos[1], pht.cpos[2]);
    printf("    PODZ: %f\n", pht.podz);
    printf("    CODZ: %f\n", pht.codz);
    printf("    SDIR: %f, %f\n", pht.sdir[0] * 180.0 / M_PI, pht.sdir[1] * 180.0 / M_PI);
    printf("    CDIR: %f, %f, %f\n", pht.cdir[0], pht.cdir[1], pht.cdir[2]);
    printf("    STKS: %f, %f, %f, %f\n", pht.stks[0], pht.stks[1], pht.stks[2], pht.stks[3]);
    printf("    SCAT: %d\n\n", pht.scat);
}
