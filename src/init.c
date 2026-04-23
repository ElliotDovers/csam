#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*
 NOTE:
 Argument counts correspond to TMB >= 1.9.0.
 If TMB updates its .Call() API, this table must be synced.
 */


/* ---- TMB .Call entry points (forward declarations) ---- */
extern SEXP MakeADFunObject(SEXP, SEXP, SEXP, SEXP);
extern SEXP InfoADFunObject(SEXP);
extern SEXP EvalADFunObject(SEXP, SEXP, SEXP);
extern SEXP MakeDoubleFunObject(SEXP, SEXP, SEXP, SEXP);
extern SEXP EvalDoubleFunObject(SEXP, SEXP, SEXP);
extern SEXP getParameterOrder(SEXP, SEXP, SEXP, SEXP);
extern SEXP MakeADGradObject(SEXP, SEXP, SEXP);
extern SEXP MakeADHessObject2(SEXP, SEXP, SEXP, SEXP);
extern SEXP usingAtomics(void);
extern SEXP TMBconfig(SEXP, SEXP);

/* ---- Your future IRLS entry ---- */
extern SEXP penalised_glm_irls_cpp(
    SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP
);


/*
 NOTE:
 Argument counts correspond to TMB >= 1.9.0.
 If TMB updates its .Call() API, this table must be synced.
 */


static const R_CallMethodDef CallEntries[] = {
  {"MakeADFunObject",     (DL_FUNC)&MakeADFunObject,     4},
  {"InfoADFunObject",     (DL_FUNC)&InfoADFunObject,     1},
  {"EvalADFunObject",     (DL_FUNC)&EvalADFunObject,     3},
  {"MakeDoubleFunObject", (DL_FUNC)&MakeDoubleFunObject, 4},
  {"EvalDoubleFunObject", (DL_FUNC)&EvalDoubleFunObject, 3},
  {"getParameterOrder",   (DL_FUNC)&getParameterOrder,   4},
  {"MakeADGradObject",    (DL_FUNC)&MakeADGradObject,    3},
  {"MakeADHessObject2",   (DL_FUNC)&MakeADHessObject2,   4},
  {"usingAtomics",        (DL_FUNC)&usingAtomics,        0},
  {"TMBconfig",           (DL_FUNC)&TMBconfig,           2},

  {"penalised_glm_irls_cpp", (DL_FUNC)&penalised_glm_irls_cpp, 10},
  {NULL, NULL, 0}
};

void R_init_csam(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
