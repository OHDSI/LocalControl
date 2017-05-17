#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP LocalControl_getMaxDist(SEXP);
extern SEXP LocalControl_newCRLC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LocalControl_newLC(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"LocalControl_getMaxDist", (DL_FUNC) &LocalControl_getMaxDist, 1},
    {"LocalControl_newCRLC",    (DL_FUNC) &LocalControl_newCRLC,    5},
    {"LocalControl_newLC",      (DL_FUNC) &LocalControl_newLC,      3},
    {NULL, NULL, 0}
};

void R_init_LocalControl(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

