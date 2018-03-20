#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void pathgraph_lammax(void *, void *, void *, void *, void *, void *);
extern void pathgraph_prox(void *, void *, void *, void *, void *, void *);
extern void pathgraph_prox2(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"pathgraph_lammax", (DL_FUNC) &pathgraph_lammax, 6},
    {"pathgraph_prox",   (DL_FUNC) &pathgraph_prox,   6},
    {"pathgraph_prox2",  (DL_FUNC) &pathgraph_prox2,  7},
    {NULL, NULL, 0}
};

void R_init_hsm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
