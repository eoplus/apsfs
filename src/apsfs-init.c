

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "mc_psf.h"

void R_init_apsfs(DllInfo *info) {
  R_RegisterCCallable("apsfs", "mc_psf", (DL_FUNC) &mc_psf);
}
