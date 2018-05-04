//// File Name: init.c
//// File Version: 2.007024
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP _sirt_eigenvaluesDsirt(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_firsteigenvalsirt2(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_parameters_jackknife(SEXP);
extern SEXP _sirt_evm_aux(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_choppin_rowaveraging(SEXP, SEXP, SEXP);
extern SEXP _sirt_evm_comp_matrix_poly(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_gooijer_csn_table(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_isop_tests_C(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_ia_optim_lambda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_ia_optim_nu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_interval_index_C(SEXP, SEXP);
extern SEXP _sirt_rowCumsums2_source(SEXP);
extern SEXP _sirt_rowKSmallest_C(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_rowMaxsCPP_source(SEXP);
extern SEXP _sirt_rowmins2_bundle_C(SEXP, SEXP);
extern SEXP _sirt_calc_copula_itemcluster_C(SEXP, SEXP);
extern SEXP _sirt_md_pattern_csource(SEXP);
extern SEXP _sirt_monoreg_rowwise_Cpp(SEXP, SEXP);
extern SEXP _sirt_mle_pcm_group_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_mle_pcm_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_noharm_compute_dj(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_noharm_compute_ej(SEXP, SEXP);
extern SEXP _sirt_noharm_compute_vjk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_noharm_compute_gammajk(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_noharm_compute_pdfk(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_noharm_estFcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_noharm_estPcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_noharm_estPsicpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_pbivnorm2_rcpp(SEXP, SEXP, SEXP);
extern SEXP _sirt_dmvnorm_2dim_rcpp(SEXP, SEXP, SEXP);
extern SEXP _sirt_polychoric2_estequation(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_polychoric2_itempair(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_tetrachoric2_rcpp_aux(SEXP, SEXP, SEXP);
extern SEXP _sirt_polychoric2_aux_rcpp(SEXP, SEXP, SEXP);
extern SEXP _sirt_probs_pcm_groups_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_probs_pcm_nogroups_C(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_calccounts_pcm_groups_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_rm_arraymult1(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_RM_CALCPOST(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_rm_probraterfct1(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_rm_facets_calcprobs_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_SMIRT_CALCPOST(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_SMIRT_CALCPROB_COMP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_SMIRT_CALCPROB_NONCOMP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_SMIRT_CALCPROB_PARTCOMP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_MML2_RASCHTYPE_COUNTS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_MML2_CALCPOST_V1(SEXP, SEXP, SEXP);
extern SEXP _sirt_MML2_CALCPOST_V2(SEXP, SEXP, SEXP);
extern SEXP _sirt_MML2_CALCPOST_V3(SEXP, SEXP, SEXP);
extern SEXP _sirt_sirt_rcpp_rm_proc_datasets_na_indicators(SEXP, SEXP);
extern SEXP _sirt_sirt_rcpp_rm_proc_expand_dataset(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sirt_sirt_rcpp_xxirt_compute_posterior_expected_counts(SEXP, SEXP);
extern SEXP _sirt_sirt_rcpp_xxirt_compute_likelihood(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_sirt_eigenvaluesDsirt", (DL_FUNC) &_sirt_eigenvaluesDsirt, 4},
    {"_sirt_firsteigenvalsirt2", (DL_FUNC) &_sirt_firsteigenvalsirt2, 4},
    {"_sirt_parameters_jackknife", (DL_FUNC) &_sirt_parameters_jackknife, 1},
    {"_sirt_evm_aux", (DL_FUNC) &_sirt_evm_aux, 6},
    {"_sirt_choppin_rowaveraging", (DL_FUNC) &_sirt_choppin_rowaveraging, 3},
    {"_sirt_evm_comp_matrix_poly", (DL_FUNC) &_sirt_evm_comp_matrix_poly, 9},
    {"_sirt_gooijer_csn_table", (DL_FUNC) &_sirt_gooijer_csn_table, 7},
    {"_sirt_isop_tests_C", (DL_FUNC) &_sirt_isop_tests_C, 5},
    {"_sirt_ia_optim_lambda", (DL_FUNC) &_sirt_ia_optim_lambda, 8},
    {"_sirt_ia_optim_nu", (DL_FUNC) &_sirt_ia_optim_nu, 11},
    {"_sirt_interval_index_C", (DL_FUNC) &_sirt_interval_index_C, 2},
    {"_sirt_rowCumsums2_source", (DL_FUNC) &_sirt_rowCumsums2_source, 1},
    {"_sirt_rowKSmallest_C", (DL_FUNC) &_sirt_rowKSmallest_C, 4},
    {"_sirt_rowMaxsCPP_source", (DL_FUNC) &_sirt_rowMaxsCPP_source, 1},
    {"_sirt_rowmins2_bundle_C", (DL_FUNC) &_sirt_rowmins2_bundle_C, 2},
    {"_sirt_calc_copula_itemcluster_C", (DL_FUNC) &_sirt_calc_copula_itemcluster_C, 2},
    {"_sirt_md_pattern_csource", (DL_FUNC) &_sirt_md_pattern_csource, 1},
    {"_sirt_monoreg_rowwise_Cpp", (DL_FUNC) &_sirt_monoreg_rowwise_Cpp, 2},
    {"_sirt_mle_pcm_group_C", (DL_FUNC) &_sirt_mle_pcm_group_C, 9},
    {"_sirt_mle_pcm_C", (DL_FUNC) &_sirt_mle_pcm_C, 8},
    {"_sirt_noharm_compute_dj", (DL_FUNC) &_sirt_noharm_compute_dj, 4},
    {"_sirt_noharm_compute_ej", (DL_FUNC) &_sirt_noharm_compute_ej, 2},
    {"_sirt_noharm_compute_vjk", (DL_FUNC) &_sirt_noharm_compute_vjk, 7},
    {"_sirt_noharm_compute_gammajk", (DL_FUNC) &_sirt_noharm_compute_gammajk, 5},
    {"_sirt_noharm_compute_pdfk", (DL_FUNC) &_sirt_noharm_compute_pdfk, 4},
    {"_sirt_noharm_estFcpp", (DL_FUNC) &_sirt_noharm_estFcpp, 16},
    {"_sirt_noharm_estPcpp", (DL_FUNC) &_sirt_noharm_estPcpp, 16},
    {"_sirt_noharm_estPsicpp", (DL_FUNC) &_sirt_noharm_estPsicpp, 16},
    {"_sirt_pbivnorm2_rcpp", (DL_FUNC) &_sirt_pbivnorm2_rcpp, 3},
    {"_sirt_dmvnorm_2dim_rcpp", (DL_FUNC) &_sirt_dmvnorm_2dim_rcpp, 3},
    {"_sirt_polychoric2_estequation", (DL_FUNC) &_sirt_polychoric2_estequation, 7},
    {"_sirt_polychoric2_itempair", (DL_FUNC) &_sirt_polychoric2_itempair, 4},
    {"_sirt_tetrachoric2_rcpp_aux", (DL_FUNC) &_sirt_tetrachoric2_rcpp_aux, 3},
    {"_sirt_polychoric2_aux_rcpp", (DL_FUNC) &_sirt_polychoric2_aux_rcpp, 3},
    {"_sirt_probs_pcm_groups_C", (DL_FUNC) &_sirt_probs_pcm_groups_C, 6},
    {"_sirt_probs_pcm_nogroups_C", (DL_FUNC) &_sirt_probs_pcm_nogroups_C, 5},
    {"_sirt_calccounts_pcm_groups_C", (DL_FUNC) &_sirt_calccounts_pcm_groups_C, 7},
    {"_sirt_rm_arraymult1", (DL_FUNC) &_sirt_rm_arraymult1, 4},
    {"_sirt_RM_CALCPOST", (DL_FUNC) &_sirt_RM_CALCPOST, 4},
    {"_sirt_rm_probraterfct1", (DL_FUNC) &_sirt_rm_probraterfct1, 5},
    {"_sirt_rm_facets_calcprobs_cpp", (DL_FUNC) &_sirt_rm_facets_calcprobs_cpp, 12},
    {"_sirt_SMIRT_CALCPOST", (DL_FUNC) &_sirt_SMIRT_CALCPOST, 6},
    {"_sirt_SMIRT_CALCPROB_COMP", (DL_FUNC) &_sirt_SMIRT_CALCPROB_COMP, 6},
    {"_sirt_SMIRT_CALCPROB_NONCOMP", (DL_FUNC) &_sirt_SMIRT_CALCPROB_NONCOMP, 6},
    {"_sirt_SMIRT_CALCPROB_PARTCOMP", (DL_FUNC) &_sirt_SMIRT_CALCPROB_PARTCOMP, 7},
    {"_sirt_MML2_RASCHTYPE_COUNTS", (DL_FUNC) &_sirt_MML2_RASCHTYPE_COUNTS, 6},
    {"_sirt_MML2_CALCPOST_V1", (DL_FUNC) &_sirt_MML2_CALCPOST_V1, 3},
    {"_sirt_MML2_CALCPOST_V2", (DL_FUNC) &_sirt_MML2_CALCPOST_V2, 3},
    {"_sirt_MML2_CALCPOST_V3", (DL_FUNC) &_sirt_MML2_CALCPOST_V3, 3},
    {"_sirt_sirt_rcpp_rm_proc_datasets_na_indicators", (DL_FUNC) &_sirt_sirt_rcpp_rm_proc_datasets_na_indicators, 2},
    {"_sirt_sirt_rcpp_rm_proc_expand_dataset", (DL_FUNC) &_sirt_sirt_rcpp_rm_proc_expand_dataset, 5},
    {"_sirt_sirt_rcpp_xxirt_compute_posterior_expected_counts", (DL_FUNC) &_sirt_sirt_rcpp_xxirt_compute_posterior_expected_counts, 2},
    {"_sirt_sirt_rcpp_xxirt_compute_likelihood", (DL_FUNC) &_sirt_sirt_rcpp_xxirt_compute_likelihood, 5},
    {NULL, NULL, 0}
};

void R_init_sirt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
