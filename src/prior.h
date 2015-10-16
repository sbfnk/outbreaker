
#ifndef __PRIOR_H
#define __PRIOR_H

void filter_logprob(double *in);

double gsl_ran_exponential_pdf_fixed(double x, double mu);

/* double logprior_alpha_i(int i, param *par); */

/* double logprior_kappa_i(int i, param *par); */

double logprior_pi(param *par);

double logprior_phi(param *par);

double logprior_spa1(param *par);

double logprior_spa2(param *par);

double logprior_mu1(param *par);

double logprior_gamma(param *par);

double logprior_trans_mat(param *par,int i);

double logprior_dirichlet_tmat(double in[], double mult, int num_groups);

double logprior_all(param *par, mcmc_param *mcmcPar);

#endif
