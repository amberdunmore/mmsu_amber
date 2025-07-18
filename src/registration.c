#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void model_initmod_desolve(void *);
extern void model_output_dde(void *);
extern void model_rhs_dde(void *);
extern void model_rhs_desolve(void *);
extern void model_sr_initmod_desolve(void *);
extern void model_sr_output_dde(void *);
extern void model_sr_rhs_dde(void *);
extern void model_sr_rhs_desolve(void *);

/* .Call calls */
extern SEXP model_contents(SEXP);
extern SEXP model_create(SEXP);
extern SEXP model_initial_conditions(SEXP, SEXP);
extern SEXP model_metadata(SEXP);
extern SEXP model_rhs_r(SEXP, SEXP, SEXP);
extern SEXP model_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP model_set_user(SEXP, SEXP);
extern SEXP model_sr_contents(SEXP);
extern SEXP model_sr_create(SEXP);
extern SEXP model_sr_initial_conditions(SEXP, SEXP);
extern SEXP model_sr_metadata(SEXP);
extern SEXP model_sr_rhs_r(SEXP, SEXP, SEXP);
extern SEXP model_sr_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP model_sr_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"model_initmod_desolve",    (DL_FUNC) &model_initmod_desolve,    1},
    {"model_output_dde",         (DL_FUNC) &model_output_dde,         1},
    {"model_rhs_dde",            (DL_FUNC) &model_rhs_dde,            1},
    {"model_rhs_desolve",        (DL_FUNC) &model_rhs_desolve,        1},
    {"model_sr_initmod_desolve", (DL_FUNC) &model_sr_initmod_desolve, 1},
    {"model_sr_output_dde",      (DL_FUNC) &model_sr_output_dde,      1},
    {"model_sr_rhs_dde",         (DL_FUNC) &model_sr_rhs_dde,         1},
    {"model_sr_rhs_desolve",     (DL_FUNC) &model_sr_rhs_desolve,     1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"model_contents",              (DL_FUNC) &model_contents,              1},
    {"model_create",                (DL_FUNC) &model_create,                1},
    {"model_initial_conditions",    (DL_FUNC) &model_initial_conditions,    2},
    {"model_metadata",              (DL_FUNC) &model_metadata,              1},
    {"model_rhs_r",                 (DL_FUNC) &model_rhs_r,                 3},
    {"model_set_initial",           (DL_FUNC) &model_set_initial,           4},
    {"model_set_user",              (DL_FUNC) &model_set_user,              2},
    {"model_sr_contents",           (DL_FUNC) &model_sr_contents,           1},
    {"model_sr_create",             (DL_FUNC) &model_sr_create,             1},
    {"model_sr_initial_conditions", (DL_FUNC) &model_sr_initial_conditions, 2},
    {"model_sr_metadata",           (DL_FUNC) &model_sr_metadata,           1},
    {"model_sr_rhs_r",              (DL_FUNC) &model_sr_rhs_r,              3},
    {"model_sr_set_initial",        (DL_FUNC) &model_sr_set_initial,        4},
    {"model_sr_set_user",           (DL_FUNC) &model_sr_set_user,           2},
    {NULL, NULL, 0}
};

void R_init_mmsu(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
