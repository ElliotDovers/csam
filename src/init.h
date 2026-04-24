#include <stdlib.h>
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

  /* Your own native routines */
  SEXP penalised_glm_irls_cpp(
      SEXP, SEXP, SEXP, SEXP,
      SEXP, SEXP, SEXP, SEXP,
      SEXP, SEXP
  );

  const static R_CallMethodDef R_CallDef[] = {

    /* ---- TMB framework entry points ---- */
    TMB_CALLDEFS,

    /* ---- Your own routines ---- */
    CALLDEF(penalised_glm_irls_cpp, 10),

    {NULL, NULL, 0}
  };

  void R_init_csam(DllInfo *dll)
  {
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

#ifdef TMB_CCALLABLES
    TMB_CCALLABLES("csam");
#endif
  }
}

