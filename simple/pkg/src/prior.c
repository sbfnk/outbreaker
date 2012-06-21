
#include "common.h"
#include "structures.h"
#include "prior.h"


/* p(alpha_i): = 1/(n-1) if alpha_i \neq i, and zero otherwise */
double logprior_alpha_i(int i, param *par){
    double out = (i==vec_int_i(par->alpha,i)) ? 0.0 : (double) 1.0/(par->n-1.0);
    return log(out);
}



/* p(kappa_i) = NB(1, 1-pi) with pi: prop obs cases */
double logprior_kappa_i(int i, param *par){
    double out = gsl_ran_negative_binomial_pdf((unsigned int) vec_int_i(par->kappa,i)-1, 1.0-par->pi, 1.0);
    return log(out);
}



/* p(pi) = beta(4,3) */
double logprior_pi(param *par){
    double out = gsl_ran_beta_pdf(par->pi, 4.0, 3.0);
    return log(out);
}



/* p(mu_1) = Unif(0,1) = 1*/
double logprior_mu1(){
    return 0.0; /* log(1) = 0 */
}



/* p(gamma) = logN(1,1.25) */
double logprior_gamma(param *par){
    double out = gsl_ran_lognormal_pdf(par->gamma, 1.0, 1.25);
    return log(out);
}



double logprior_all(param *par){
    int i;
    double out=0.0;

    /* result is the sum of priors over all parameters and individuals */
    for(i=0;i<par->n;i++){
	out += logprior_alpha_i(i,par);
	out += logprior_kappa_i(i,par);
    }

    out += logprior_pi(par);
    out += logprior_mu1();
    out += logprior_gamma(par);

    return(out);
}





/* int main(){ */
/*     double out; */
/*     /\* CREATE AND INIT PARAMETERS *\/ */
/*     param *par = alloc_param(3); */
    
/*     par->alpha->values[0] = -1; */
/*     par->alpha->values[1] = 0; */
/*     par->alpha->values[2] = 0; */
    
/*     par->kappa->values[0] = 1; */
/*     par->kappa->values[1] = 1; */
/*     par->kappa->values[2] = 1; */

/*     par->mu1 = 0.0001; */
/*     par->gamma = 1.0; */
    
/*     par->pi = 0.5; */

/*     printf("\nParameters\n"); */
/*     print_param(par); */

/*     /\* COMPUTE PRIORS *\/ */
/*     out = logprior_all(par); */
/*     printf("\nReturned value [log (prob)]: %.15f (%.15f)\n", out, exp(out)); */

/*     /\* FREE / RETURN *\/ */
/*     free_param(par); */
/*     return 0; */
/* } */



/*
  gcc instructions

  gcc -o prior matvec.c genclasses.c structures.c prior.c -lgsl -lgslcblas -Wall && ./prior

  valgrind --leak-check=full prior

*/

