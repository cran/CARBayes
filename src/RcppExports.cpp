// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// linpredcompute
NumericVector linpredcompute(NumericMatrix X, const int nsites, const int p, NumericVector beta, NumericVector offset);
RcppExport SEXP _CARBayes_linpredcompute(SEXP XSEXP, SEXP nsitesSEXP, SEXP pSEXP, SEXP betaSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(linpredcompute(X, nsites, p, beta, offset));
    return rcpp_result_gen;
END_RCPP
}
// quadform
double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, NumericVector phi, NumericVector theta, double rho);
RcppExport SEXP _CARBayes_quadform(SEXP WtripletSEXP, SEXP WtripletsumSEXP, SEXP n_tripletSEXP, SEXP nsitesSEXP, SEXP phiSEXP, SEXP thetaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Wtripletsum(WtripletsumSEXP);
    Rcpp::traits::input_parameter< const int >::type n_triplet(n_tripletSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi, theta, rho));
    return rcpp_result_gen;
END_RCPP
}
// binomialindepupdateRW
List binomialindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, const NumericVector failures, const double theta_tune, NumericVector offset);
RcppExport SEXP _CARBayes_binomialindepupdateRW(SEXP nsitesSEXP, SEXP thetaSEXP, SEXP sigma2SEXP, SEXP ySEXP, SEXP failuresSEXP, SEXP theta_tuneSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type failures(failuresSEXP);
    Rcpp::traits::input_parameter< const double >::type theta_tune(theta_tuneSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(binomialindepupdateRW(nsites, theta, sigma2, y, failures, theta_tune, offset));
    return rcpp_result_gen;
END_RCPP
}
// binomialcarupdateRW
List binomialcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, NumericVector phi, double tau2, const NumericVector y, const NumericVector failures, const double phi_tune, double rho, NumericVector offset);
RcppExport SEXP _CARBayes_binomialcarupdateRW(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP WtripletsumSEXP, SEXP nsitesSEXP, SEXP phiSEXP, SEXP tau2SEXP, SEXP ySEXP, SEXP failuresSEXP, SEXP phi_tuneSEXP, SEXP rhoSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Wtripletsum(WtripletsumSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type failures(failuresSEXP);
    Rcpp::traits::input_parameter< const double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(binomialcarupdateRW(Wtriplet, Wbegfin, Wtripletsum, nsites, phi, tau2, y, failures, phi_tune, rho, offset));
    return rcpp_result_gen;
END_RCPP
}
// binomialbetaupdateRW
List binomialbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, NumericVector offset, NumericVector y, NumericVector failures, NumericVector prior_meanbeta, NumericVector prior_varbeta, const int nblock, double beta_tune, List block_list);
RcppExport SEXP _CARBayes_binomialbetaupdateRW(SEXP XSEXP, SEXP nsitesSEXP, SEXP pSEXP, SEXP betaSEXP, SEXP offsetSEXP, SEXP ySEXP, SEXP failuresSEXP, SEXP prior_meanbetaSEXP, SEXP prior_varbetaSEXP, SEXP nblockSEXP, SEXP beta_tuneSEXP, SEXP block_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type failures(failuresSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_meanbeta(prior_meanbetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_varbeta(prior_varbetaSEXP);
    Rcpp::traits::input_parameter< const int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tune(beta_tuneSEXP);
    Rcpp::traits::input_parameter< List >::type block_list(block_listSEXP);
    rcpp_result_gen = Rcpp::wrap(binomialbetaupdateRW(X, nsites, p, beta, offset, y, failures, prior_meanbeta, prior_varbeta, nblock, beta_tune, block_list));
    return rcpp_result_gen;
END_RCPP
}
// binomialbetaupdateMALA
List binomialbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, NumericVector offset, NumericVector y, NumericVector failures, NumericVector trials, NumericVector prior_meanbeta, NumericVector prior_varbeta, const int nblock, double beta_tune, List block_list);
RcppExport SEXP _CARBayes_binomialbetaupdateMALA(SEXP XSEXP, SEXP nsitesSEXP, SEXP pSEXP, SEXP betaSEXP, SEXP offsetSEXP, SEXP ySEXP, SEXP failuresSEXP, SEXP trialsSEXP, SEXP prior_meanbetaSEXP, SEXP prior_varbetaSEXP, SEXP nblockSEXP, SEXP beta_tuneSEXP, SEXP block_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type failures(failuresSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type trials(trialsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_meanbeta(prior_meanbetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_varbeta(prior_varbetaSEXP);
    Rcpp::traits::input_parameter< const int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tune(beta_tuneSEXP);
    Rcpp::traits::input_parameter< List >::type block_list(block_listSEXP);
    rcpp_result_gen = Rcpp::wrap(binomialbetaupdateMALA(X, nsites, p, beta, offset, y, failures, trials, prior_meanbeta, prior_varbeta, nblock, beta_tune, block_list));
    return rcpp_result_gen;
END_RCPP
}
// poissonindepupdateRW
List poissonindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, const double theta_tune, NumericVector offset);
RcppExport SEXP _CARBayes_poissonindepupdateRW(SEXP nsitesSEXP, SEXP thetaSEXP, SEXP sigma2SEXP, SEXP ySEXP, SEXP theta_tuneSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type theta_tune(theta_tuneSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(poissonindepupdateRW(nsites, theta, sigma2, y, theta_tune, offset));
    return rcpp_result_gen;
END_RCPP
}
// poissoncarupdateRW
List poissoncarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, NumericVector phi, double tau2, const NumericVector y, const double phi_tune, double rho, NumericVector offset);
RcppExport SEXP _CARBayes_poissoncarupdateRW(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP WtripletsumSEXP, SEXP nsitesSEXP, SEXP phiSEXP, SEXP tau2SEXP, SEXP ySEXP, SEXP phi_tuneSEXP, SEXP rhoSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Wtripletsum(WtripletsumSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(poissoncarupdateRW(Wtriplet, Wbegfin, Wtripletsum, nsites, phi, tau2, y, phi_tune, rho, offset));
    return rcpp_result_gen;
END_RCPP
}
// poissonbetaupdateMALA
List poissonbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, NumericVector offset, NumericVector y, NumericVector prior_meanbeta, NumericVector prior_varbeta, const int nblock, double beta_tune, List block_list);
RcppExport SEXP _CARBayes_poissonbetaupdateMALA(SEXP XSEXP, SEXP nsitesSEXP, SEXP pSEXP, SEXP betaSEXP, SEXP offsetSEXP, SEXP ySEXP, SEXP prior_meanbetaSEXP, SEXP prior_varbetaSEXP, SEXP nblockSEXP, SEXP beta_tuneSEXP, SEXP block_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_meanbeta(prior_meanbetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_varbeta(prior_varbetaSEXP);
    Rcpp::traits::input_parameter< const int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tune(beta_tuneSEXP);
    Rcpp::traits::input_parameter< List >::type block_list(block_listSEXP);
    rcpp_result_gen = Rcpp::wrap(poissonbetaupdateMALA(X, nsites, p, beta, offset, y, prior_meanbeta, prior_varbeta, nblock, beta_tune, block_list));
    return rcpp_result_gen;
END_RCPP
}
// poissonbetaupdateRW
List poissonbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, NumericVector offset, NumericVector y, NumericVector prior_meanbeta, NumericVector prior_varbeta, const int nblock, double beta_tune, List block_list);
RcppExport SEXP _CARBayes_poissonbetaupdateRW(SEXP XSEXP, SEXP nsitesSEXP, SEXP pSEXP, SEXP betaSEXP, SEXP offsetSEXP, SEXP ySEXP, SEXP prior_meanbetaSEXP, SEXP prior_varbetaSEXP, SEXP nblockSEXP, SEXP beta_tuneSEXP, SEXP block_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_meanbeta(prior_meanbetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_varbeta(prior_varbetaSEXP);
    Rcpp::traits::input_parameter< const int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tune(beta_tuneSEXP);
    Rcpp::traits::input_parameter< List >::type block_list(block_listSEXP);
    rcpp_result_gen = Rcpp::wrap(poissonbetaupdateRW(X, nsites, p, beta, offset, y, prior_meanbeta, prior_varbeta, nblock, beta_tune, block_list));
    return rcpp_result_gen;
END_RCPP
}
// zipcarupdateRW
List zipcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, NumericVector phi, double tau2, const NumericVector y, const double phi_tune, double rho, NumericVector offset, NumericVector poiind);
RcppExport SEXP _CARBayes_zipcarupdateRW(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP WtripletsumSEXP, SEXP nsitesSEXP, SEXP phiSEXP, SEXP tau2SEXP, SEXP ySEXP, SEXP phi_tuneSEXP, SEXP rhoSEXP, SEXP offsetSEXP, SEXP poiindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Wtripletsum(WtripletsumSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type poiind(poiindSEXP);
    rcpp_result_gen = Rcpp::wrap(zipcarupdateRW(Wtriplet, Wbegfin, Wtripletsum, nsites, phi, tau2, y, phi_tune, rho, offset, poiind));
    return rcpp_result_gen;
END_RCPP
}
// zipindepupdateRW
List zipindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, const double theta_tune, NumericVector offset, NumericVector poiind);
RcppExport SEXP _CARBayes_zipindepupdateRW(SEXP nsitesSEXP, SEXP thetaSEXP, SEXP sigma2SEXP, SEXP ySEXP, SEXP theta_tuneSEXP, SEXP offsetSEXP, SEXP poiindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type theta_tune(theta_tuneSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type poiind(poiindSEXP);
    rcpp_result_gen = Rcpp::wrap(zipindepupdateRW(nsites, theta, sigma2, y, theta_tune, offset, poiind));
    return rcpp_result_gen;
END_RCPP
}
// gaussiancarupdate
NumericVector gaussiancarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, NumericVector phi, double tau2, double rho, double nu2, NumericVector offset);
RcppExport SEXP _CARBayes_gaussiancarupdate(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP WtripletsumSEXP, SEXP nsitesSEXP, SEXP phiSEXP, SEXP tau2SEXP, SEXP rhoSEXP, SEXP nu2SEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Wtripletsum(WtripletsumSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type nu2(nu2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussiancarupdate(Wtriplet, Wbegfin, Wtripletsum, nsites, phi, tau2, rho, nu2, offset));
    return rcpp_result_gen;
END_RCPP
}
// binomialmcarupdateRW
List binomialmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, const int nsites, const int nvar, NumericMatrix phi, NumericMatrix Y, NumericMatrix failures, NumericMatrix phioffset, NumericVector denoffset, NumericMatrix Sigmainv, double rho, double phi_tune, NumericMatrix innovations);
RcppExport SEXP _CARBayes_binomialmcarupdateRW(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP nsitesSEXP, SEXP nvarSEXP, SEXP phiSEXP, SEXP YSEXP, SEXP failuresSEXP, SEXP phioffsetSEXP, SEXP denoffsetSEXP, SEXP SigmainvSEXP, SEXP rhoSEXP, SEXP phi_tuneSEXP, SEXP innovationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type failures(failuresSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phioffset(phioffsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type denoffset(denoffsetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Sigmainv(SigmainvSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type innovations(innovationsSEXP);
    rcpp_result_gen = Rcpp::wrap(binomialmcarupdateRW(Wtriplet, Wbegfin, nsites, nvar, phi, Y, failures, phioffset, denoffset, Sigmainv, rho, phi_tune, innovations));
    return rcpp_result_gen;
END_RCPP
}
// poissonmcarupdateRW
List poissonmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, const int nsites, const int nvar, NumericMatrix phi, NumericMatrix Y, NumericMatrix phioffset, NumericVector denoffset, NumericMatrix Sigmainv, double rho, double phi_tune, NumericMatrix innovations);
RcppExport SEXP _CARBayes_poissonmcarupdateRW(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP nsitesSEXP, SEXP nvarSEXP, SEXP phiSEXP, SEXP YSEXP, SEXP phioffsetSEXP, SEXP denoffsetSEXP, SEXP SigmainvSEXP, SEXP rhoSEXP, SEXP phi_tuneSEXP, SEXP innovationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phioffset(phioffsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type denoffset(denoffsetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Sigmainv(SigmainvSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type innovations(innovationsSEXP);
    rcpp_result_gen = Rcpp::wrap(poissonmcarupdateRW(Wtriplet, Wbegfin, nsites, nvar, phi, Y, phioffset, denoffset, Sigmainv, rho, phi_tune, innovations));
    return rcpp_result_gen;
END_RCPP
}
// gaussianmcarupdateRW
List gaussianmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, const int nsites, const int nvar, NumericMatrix phi, NumericMatrix phioffset, NumericVector denoffset, NumericMatrix Sigmainv, double rho, NumericVector nu2, double phi_tune, NumericMatrix innovations);
RcppExport SEXP _CARBayes_gaussianmcarupdateRW(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP nsitesSEXP, SEXP nvarSEXP, SEXP phiSEXP, SEXP phioffsetSEXP, SEXP denoffsetSEXP, SEXP SigmainvSEXP, SEXP rhoSEXP, SEXP nu2SEXP, SEXP phi_tuneSEXP, SEXP innovationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phioffset(phioffsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type denoffset(denoffsetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Sigmainv(SigmainvSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu2(nu2SEXP);
    Rcpp::traits::input_parameter< double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type innovations(innovationsSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussianmcarupdateRW(Wtriplet, Wbegfin, nsites, nvar, phi, phioffset, denoffset, Sigmainv, rho, nu2, phi_tune, innovations));
    return rcpp_result_gen;
END_RCPP
}
// multinomialbetaupdateRW
List multinomialbetaupdateRW(NumericMatrix X, const int nsites, const int J, const int p, const int col, NumericMatrix beta, NumericMatrix offset, NumericMatrix y, NumericVector prior_meanbeta, NumericVector prior_varbeta, const int nblock, double beta_tune, List block_list, NumericVector zeros);
RcppExport SEXP _CARBayes_multinomialbetaupdateRW(SEXP XSEXP, SEXP nsitesSEXP, SEXP JSEXP, SEXP pSEXP, SEXP colSEXP, SEXP betaSEXP, SEXP offsetSEXP, SEXP ySEXP, SEXP prior_meanbetaSEXP, SEXP prior_varbetaSEXP, SEXP nblockSEXP, SEXP beta_tuneSEXP, SEXP block_listSEXP, SEXP zerosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type J(JSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type col(colSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_meanbeta(prior_meanbetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_varbeta(prior_varbetaSEXP);
    Rcpp::traits::input_parameter< const int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tune(beta_tuneSEXP);
    Rcpp::traits::input_parameter< List >::type block_list(block_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zeros(zerosSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomialbetaupdateRW(X, nsites, J, p, col, beta, offset, y, prior_meanbeta, prior_varbeta, nblock, beta_tune, block_list, zeros));
    return rcpp_result_gen;
END_RCPP
}
// multinomialmcarupdateRW
List multinomialmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, const int nsites, const int nvar, NumericMatrix phi, NumericMatrix Y, NumericMatrix phioffset, NumericVector denoffset, NumericMatrix Sigmainv, double rho, double phi_tune, NumericMatrix innovations);
RcppExport SEXP _CARBayes_multinomialmcarupdateRW(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP nsitesSEXP, SEXP nvarSEXP, SEXP phiSEXP, SEXP YSEXP, SEXP phioffsetSEXP, SEXP denoffsetSEXP, SEXP SigmainvSEXP, SEXP rhoSEXP, SEXP phi_tuneSEXP, SEXP innovationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< const int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phioffset(phioffsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type denoffset(denoffsetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Sigmainv(SigmainvSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type innovations(innovationsSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomialmcarupdateRW(Wtriplet, Wbegfin, nsites, nvar, phi, Y, phioffset, denoffset, Sigmainv, rho, phi_tune, innovations));
    return rcpp_result_gen;
END_RCPP
}
// gaussiancarmultilevelupdate
NumericVector gaussiancarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, NumericVector n_individual, const int nsites, NumericVector phi, double tau2, double rho, double nu2, NumericVector offset);
RcppExport SEXP _CARBayes_gaussiancarmultilevelupdate(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP WtripletsumSEXP, SEXP n_individualSEXP, SEXP nsitesSEXP, SEXP phiSEXP, SEXP tau2SEXP, SEXP rhoSEXP, SEXP nu2SEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Wtripletsum(WtripletsumSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_individual(n_individualSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type nu2(nu2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussiancarmultilevelupdate(Wtriplet, Wbegfin, Wtripletsum, n_individual, nsites, phi, tau2, rho, nu2, offset));
    return rcpp_result_gen;
END_RCPP
}
// binomialcarmultilevelupdate
List binomialcarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, List ind_area_list, NumericVector n_individual, const int nsites, NumericVector phi, double tau2, const NumericVector y, const NumericVector failures, const double phi_tune, double rho, NumericVector offset);
RcppExport SEXP _CARBayes_binomialcarmultilevelupdate(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP WtripletsumSEXP, SEXP ind_area_listSEXP, SEXP n_individualSEXP, SEXP nsitesSEXP, SEXP phiSEXP, SEXP tau2SEXP, SEXP ySEXP, SEXP failuresSEXP, SEXP phi_tuneSEXP, SEXP rhoSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Wtripletsum(WtripletsumSEXP);
    Rcpp::traits::input_parameter< List >::type ind_area_list(ind_area_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_individual(n_individualSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type failures(failuresSEXP);
    Rcpp::traits::input_parameter< const double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(binomialcarmultilevelupdate(Wtriplet, Wbegfin, Wtripletsum, ind_area_list, n_individual, nsites, phi, tau2, y, failures, phi_tune, rho, offset));
    return rcpp_result_gen;
END_RCPP
}
// poissoncarmultilevelupdate
List poissoncarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, List ind_area_list, NumericVector n_individual, const int nsites, NumericVector phi, double tau2, const NumericVector y, const double phi_tune, double rho, NumericVector offset);
RcppExport SEXP _CARBayes_poissoncarmultilevelupdate(SEXP WtripletSEXP, SEXP WbegfinSEXP, SEXP WtripletsumSEXP, SEXP ind_area_listSEXP, SEXP n_individualSEXP, SEXP nsitesSEXP, SEXP phiSEXP, SEXP tau2SEXP, SEXP ySEXP, SEXP phi_tuneSEXP, SEXP rhoSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Wtriplet(WtripletSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Wbegfin(WbegfinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Wtripletsum(WtripletsumSEXP);
    Rcpp::traits::input_parameter< List >::type ind_area_list(ind_area_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_individual(n_individualSEXP);
    Rcpp::traits::input_parameter< const int >::type nsites(nsitesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type phi_tune(phi_tuneSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(poissoncarmultilevelupdate(Wtriplet, Wbegfin, Wtripletsum, ind_area_list, n_individual, nsites, phi, tau2, y, phi_tune, rho, offset));
    return rcpp_result_gen;
END_RCPP
}
// basiscomputeinverse
NumericMatrix basiscomputeinverse(NumericMatrix D, int nrows, int ncols, NumericVector Z, int startcol);
RcppExport SEXP _CARBayes_basiscomputeinverse(SEXP DSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP, SEXP ZSEXP, SEXP startcolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type startcol(startcolSEXP);
    rcpp_result_gen = Rcpp::wrap(basiscomputeinverse(D, nrows, ncols, Z, startcol));
    return rcpp_result_gen;
END_RCPP
}
// basiscomputelinear
NumericMatrix basiscomputelinear(NumericMatrix D, int nrows, int ncols, NumericVector Z, int startcol);
RcppExport SEXP _CARBayes_basiscomputelinear(SEXP DSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP, SEXP ZSEXP, SEXP startcolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type startcol(startcolSEXP);
    rcpp_result_gen = Rcpp::wrap(basiscomputelinear(D, nrows, ncols, Z, startcol));
    return rcpp_result_gen;
END_RCPP
}
// basiscomputeexponential
NumericMatrix basiscomputeexponential(NumericMatrix D, int nrows, int ncols, NumericVector Z, int startcol);
RcppExport SEXP _CARBayes_basiscomputeexponential(SEXP DSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP, SEXP ZSEXP, SEXP startcolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type startcol(startcolSEXP);
    rcpp_result_gen = Rcpp::wrap(basiscomputeexponential(D, nrows, ncols, Z, startcol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CARBayes_linpredcompute", (DL_FUNC) &_CARBayes_linpredcompute, 5},
    {"_CARBayes_quadform", (DL_FUNC) &_CARBayes_quadform, 7},
    {"_CARBayes_binomialindepupdateRW", (DL_FUNC) &_CARBayes_binomialindepupdateRW, 7},
    {"_CARBayes_binomialcarupdateRW", (DL_FUNC) &_CARBayes_binomialcarupdateRW, 11},
    {"_CARBayes_binomialbetaupdateRW", (DL_FUNC) &_CARBayes_binomialbetaupdateRW, 12},
    {"_CARBayes_binomialbetaupdateMALA", (DL_FUNC) &_CARBayes_binomialbetaupdateMALA, 13},
    {"_CARBayes_poissonindepupdateRW", (DL_FUNC) &_CARBayes_poissonindepupdateRW, 6},
    {"_CARBayes_poissoncarupdateRW", (DL_FUNC) &_CARBayes_poissoncarupdateRW, 10},
    {"_CARBayes_poissonbetaupdateMALA", (DL_FUNC) &_CARBayes_poissonbetaupdateMALA, 11},
    {"_CARBayes_poissonbetaupdateRW", (DL_FUNC) &_CARBayes_poissonbetaupdateRW, 11},
    {"_CARBayes_zipcarupdateRW", (DL_FUNC) &_CARBayes_zipcarupdateRW, 11},
    {"_CARBayes_zipindepupdateRW", (DL_FUNC) &_CARBayes_zipindepupdateRW, 7},
    {"_CARBayes_gaussiancarupdate", (DL_FUNC) &_CARBayes_gaussiancarupdate, 9},
    {"_CARBayes_binomialmcarupdateRW", (DL_FUNC) &_CARBayes_binomialmcarupdateRW, 13},
    {"_CARBayes_poissonmcarupdateRW", (DL_FUNC) &_CARBayes_poissonmcarupdateRW, 12},
    {"_CARBayes_gaussianmcarupdateRW", (DL_FUNC) &_CARBayes_gaussianmcarupdateRW, 12},
    {"_CARBayes_multinomialbetaupdateRW", (DL_FUNC) &_CARBayes_multinomialbetaupdateRW, 14},
    {"_CARBayes_multinomialmcarupdateRW", (DL_FUNC) &_CARBayes_multinomialmcarupdateRW, 12},
    {"_CARBayes_gaussiancarmultilevelupdate", (DL_FUNC) &_CARBayes_gaussiancarmultilevelupdate, 10},
    {"_CARBayes_binomialcarmultilevelupdate", (DL_FUNC) &_CARBayes_binomialcarmultilevelupdate, 13},
    {"_CARBayes_poissoncarmultilevelupdate", (DL_FUNC) &_CARBayes_poissoncarmultilevelupdate, 12},
    {"_CARBayes_basiscomputeinverse", (DL_FUNC) &_CARBayes_basiscomputeinverse, 5},
    {"_CARBayes_basiscomputelinear", (DL_FUNC) &_CARBayes_basiscomputelinear, 5},
    {"_CARBayes_basiscomputeexponential", (DL_FUNC) &_CARBayes_basiscomputeexponential, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_CARBayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
