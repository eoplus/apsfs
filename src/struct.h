#ifndef STRUCT
#define STRUCT

// Macro definitions:

#define ABS(x) \
        (((x) >= 0.0) ? (x) : -(x))

#define NUM_EQU(x, y, pre) \
        (ABS((x) - (y)) <= (pre))

#define SIGN(x) \
        (NUM_EQU((x), 0.0, 1.0E-012) ? 0.0 : (((x) > 0.0) ? 1.0 : -1.0))

#define COS2SPH(cdir, sdir) \
        sdir[0] = acos(cdir[2]); \
        if((abs(cdir[2]) - 1.0) < 1E-12) { \
         sdir[1] = 2.0 * M_PI * gsl_rng_uniform (random); \
        } else { \
  	 sint = sin(sdir[0]); \
         if(cdir[1] >= 0.0) { \
          sdir[1] = acos(round((cdir[0] / sint) * 1.0E+012) / 1.0E+012); \
         } else { \
          sdir[1] = (2.0 * M_PI) - acos(round((cdir[0] / sint) * 1.0E+012) / 1.0E+012); \
         } \
        }

// Structure definitions:

struct str_phtpkg {
    double ppos[3]; // Previous cartesian position.          [0] Xp,  [1] Yp,  [2] Zp.
    double cpos[3]; // Current cartesian position.           [0] Xc,  [1] Yc,  [2] Zc.
    double codz;    // Current optical depth.                [0] Z OD.
    double podz;    // Previous optical depth.               [0] Z OD.
    double sdir[2]; // Spherical direction.                  [0] Theta, [1] Phi.
    double cdir[3]; // Cosines direction.                    [0] muX, [1] muY, [2] muZ.
    double stks[4]; // Stokes vector.                        [0] I, [1] Q, [2] U, [3] V.
    int    scat;    // Flag to keep track if photon was scattered.
 };
 
// Function skeletons:

void   print_pht (struct str_phtpkg pht);
int    findInterv(double val, double* brks, size_t n);
void   accm_annular (struct str_phtpkg pht, double* bin_phtw, double* bin_brks, size_t n);
void   accm_sectorial (struct str_phtpkg pht, double* bin_phtw, double* bin_brks, size_t n);
void   accm_grid (struct str_phtpkg pht, double* bin_phtw, double* bin_brks, size_t n);
void   mvpht_m (struct str_phtpkg* pht, double s, SEXP tau, SEXP tau_u, SEXP tau_d, SEXP c_tot, SEXP dkm, int cid, R_xlen_t atm_n);
void   update_cdir (double psi, double phi, double *cdir);
double scat_fun (double ru, double* cdf, double *psi, int npsi);


#endif
