#include "common.h"
#include "init.h"
#include "InputOutput.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "alloc.h"
#include "genclasses.h"
#include "matvec.h"
#include "tuneVariances.h"





/**************** nb_data ****************/
nb_data * createNbData(int NbPatients, int T, int NbSequences, int NbColonisedPatients){
    nb_data *nb = (nb_data *) malloc(sizeof(nb_data));
    if(nb == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->NbAdmissions = (int *) calloc(NbPatients, sizeof(int));
    if(nb->NbAdmissions == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->NbPosSwabs = (int *) calloc(NbPatients, sizeof(int));
    if(nb->NbPosSwabs == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->NbNegSwabs = (int *) calloc(NbPatients, sizeof(int));
    if(nb->NbNegSwabs == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->indexColonisedPatients = (int *) calloc(NbColonisedPatients, sizeof(int));
    if(nb->indexColonisedPatients == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->M = (int *) calloc(NbPatients, sizeof(int));
    if(nb->M == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->NbColonisedPatients = NbColonisedPatients;
    nb->NbPatients = NbPatients;
    nb->T = T;
    nb->NbSequences = NbSequences;

    return nb;
}





void freeNbData(nb_data *nb){
    free(nb->NbAdmissions);
    free(nb->NbPosSwabs);
    free(nb->NbNegSwabs);
    free(nb->indexColonisedPatients);
    free(nb->M);
    free(nb);
}





void print_nbData(nb_data *nb){
    int i;
    printf("\nNb of patients: %d, time span 0-%d", nb->NbPatients, nb->T);
    printf("\nNb of colonized patients: %d", nb->NbColonisedPatients);
    printf("\nNb of Admissions:\n");
    for(i=0;i<nb->NbPatients;i++) printf("%d ",nb->NbAdmissions[i]);
    printf("\nNb of positive swabs:\n");
    for(i=0;i<nb->NbPatients;i++) printf("%d ",nb->NbPosSwabs[i]);
    printf("\nNb of negative swabs:\n");
    for(i=0;i<nb->NbPatients;i++) printf("%d ",nb->NbNegSwabs[i]);
    printf("\nIndices of colonized patients:\n");
    for(i=0;i<nb->NbColonisedPatients;i++) printf("%d ",nb->indexColonisedPatients[i]);
    printf("\nNb of sequences in each patient:\n");
    for(i=0;i<nb->NbPatients;i++) printf("%d ",nb->M[i]);
    fflush(stdout);
}



/**************** raw_data ****************/
raw_data *createRawData(nb_data *nb){
    int i;
    raw_data *data = (raw_data *) malloc(sizeof(raw_data));

    /* EPI DATA */
    if(data == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->ward = (int *) calloc(nb->NbPatients, sizeof(int));
    if(data->ward == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    /* data->PatientIndex = (int *) calloc(nb->NbPatients, sizeof(int)); */
    /* if(data->PatientIndex == NULL){ */
    /* 	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n"); */
    /* 	exit(1); */
    /* } */

    /* data->timeSeq = (int *) calloc(nb->NbPatients, sizeof(int)); */
    /* if(data->timeSeq == NULL){ */
    /* 	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n"); */
    /* 	exit(1); */
    /* } */

    data->A = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->A == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->D = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->D == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->P = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->P == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->N = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->N == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->IsInHosp = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->IsInHosp == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    for(i=0;i<nb->NbPatients;i++){
	data->A[i] = gsl_vector_calloc(nb->NbAdmissions[i]);
	data->D[i] = gsl_vector_calloc(nb->NbAdmissions[i]);
	if(nb->NbPosSwabs[i]>0){
		data->P[i] = gsl_vector_calloc(nb->NbPosSwabs[i]);
		/* nb->indexColonisedPatients[nb->NbColonisedPatients] = i; */
		/* nb->NbColonisedPatients++; */
	} else {
	    data->P[i] = NULL;
	}
	if(nb->NbNegSwabs[i]>0){
	    data->N[i] = gsl_vector_calloc(nb->NbNegSwabs[i]);
	} else {
	    data->N[i] = NULL;
	}
	data->IsInHosp[i] = gsl_vector_calloc(nb->T);
    }

    /* simple integers */
    data->NbPatients = nb->NbPatients;
    data->T = nb->T;


    /* GENETIC DATA */
    /* S: indices of sequences collected for each patient */
    data->S = (int **) malloc(nb->NbPatients*sizeof(int *));
    if(data->S == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }
    for(i=0;i<nb->NbPatients;i++){
	data->S[i] = (int *) calloc(nb->M[i], sizeof(int));
	if(data->S[i] == NULL){
	    fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	    exit(1);
	}
    }

    /* Tcollec: collection times for each sequence */
    data->Tcollec = (int *) calloc(nb->NbSequences, sizeof(int));
    if(data->Tcollec == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    /* M: number of sequences collected for each patient */
    data->M = (int *) calloc(nb->NbPatients, sizeof(int));
    if(data->M == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }
    for(i=0;i<nb->NbPatients;i++){
	data->M[i] = nb->M[i];
    }
    data->NbSequences = nb->NbSequences;


    /* RANDOM NUMBER GENERATOR */
    /* time_t t = time(NULL); /\* time in seconds, used to change the seed of the random generator *\/ */
    time_t t = 1; /* time in seconds, used to change the seed of the random generator */
    const gsl_rng_type *typ;
    gsl_rng_env_setup();
    typ=gsl_rng_default;
    gsl_rng * rng=gsl_rng_alloc(typ);
    gsl_rng_set(rng,t); /* changes the seed of the random generator */
    data->rng = rng;

    return data;
}





void freeRawData(raw_data *data){
    int i;

    for(i=0 ; i<data->NbPatients ; i++){
	if(data->A[i] != NULL) gsl_vector_free(data->A[i]);
	if(data->D[i] != NULL) gsl_vector_free(data->D[i]);
	if(data->P[i] != NULL) gsl_vector_free(data->P[i]);
	if(data->N[i] != NULL) gsl_vector_free(data->N[i]);
	if(data->IsInHosp[i] != NULL) gsl_vector_free(data->IsInHosp[i]);
	free(data->S[i]);
    }

    free(data->ward);
    /* free(data->PatientIndex); */
    /* free(data->timeSeq); */
    free(data->A);
    free(data->D);
    free(data->P);
    free(data->N);
    free(data->IsInHosp);
    free(data->S);
    free(data->Tcollec);
    free(data->M);
    gsl_rng_free(data->rng);
    free(data);
}




void print_rawData(raw_data *data){
    int i,j;
    printf("\nNb of patients: %d, time span 0-%d", data->NbPatients, data->T);
    printf("\nWard data:\n");
    for(i=0;i<data->NbPatients;i++) printf("%d ", data->ward[i]);
    printf("\nAdmission dates:\n");
    for(i=0;i<data->NbPatients;i++){
	printf("\nPatient %d:\n",i);
	if(data->A[i]!=NULL) gsl_vector_fprintf(stdout, data->A[i], "%.1f"); else printf("NULL\n");
    }
    printf("\nDischarge dates:\n");
    for(i=0;i<data->NbPatients;i++){
	printf("\nPatient %d:\n",i);
	if(data->D[i]!=NULL) gsl_vector_fprintf(stdout, data->D[i], "%.1f"); else printf("NULL\n");
    }
    printf("\nDates of positive swabs:\n");
    for(i=0;i<data->NbPatients;i++){
	printf("Patient %d:\n",i);
	if(data->P[i]!=NULL) gsl_vector_fprintf(stdout, data->P[i], "%.1f"); else printf("NULL\n");
    }
    printf("\nDates of negative swabs:\n");
    for(i=0;i<data->NbPatients;i++){
	printf("Patient %d:\n",i);
	if(data->N[i]!=NULL) gsl_vector_fprintf(stdout, data->N[i], "%.1f"); else printf("NULL\n");
    }
    printf("\nHospital presence:\n");
    for(i=0;i<data->NbPatients;i++){
	printf("Patient %d:\n",i);
	gsl_vector_fprintf(stdout, data->IsInHosp[i], "%.0f");
    }
    printf("\nIndices of sequences for each patient:\n");
    for(i=0;i<data->NbPatients;i++){
	printf("\nPatient %d\n",i);
	for(j=0;j<data->M[i];j++){
	    printf("%d ", data->S[i][j]);
	}
    }
    printf("\nCollection times:\n");
    for(i=0;i<data->NbSequences;i++){
	printf("%d ", data->Tcollec[i]);
    }
    printf("\nNb of sequences per patient:\n");
    for(i=0;i<data->NbPatients;i++) printf("%d ", data->M[i]);
}







/**************** aug_data ****************/
aug_data *createAugData(int NbPatients, int T){
    aug_data *augData = (aug_data *) malloc(sizeof(aug_data));
    if(augData == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }


    augData->C = (int *) calloc(NbPatients, sizeof(int));
    if(augData->C == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }

    augData->E = (int *) calloc(NbPatients, sizeof(int));
    if(augData->E == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }

    augData->I0 = (int *) calloc(T, sizeof(int));
    if(augData->I0 == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }

    augData->I1 = (int *) calloc(T, sizeof(int));
    if(augData->I1 == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }

    augData->NbPatients = NbPatients;
    augData->T = T;

    return augData;
}





void freeAugData(aug_data *augData){
    free(augData->C);
    free(augData->E);
    free(augData->I0);
    free(augData->I1);
    free(augData);
}




void print_augData(aug_data *augData){
    int i;
    printf("\nNb of patients: %d, time span 0-%d", augData->NbPatients, augData->T);
    printf("\nColonisation times: \n");
    for(i=0;i<augData->NbPatients;i++) printf("%d ", augData->C[i]);
    printf("\nClearance times: \n");
    for(i=0;i<augData->NbPatients;i++) printf("%d ", augData->E[i]);
    printf("\nNb of colonised individuals at each time step, ward 0: \n");
    for(i=0;i<augData->T;i++) printf("%d ", augData->I0[i]);
    printf("\nNb of colonised individuals at each time step, ward 1: \n");
    for(i=0;i<augData->T;i++) printf("%d ", augData->I1[i]);
    fflush(stdout);
}




void copyAugData(aug_data *augDataDest, aug_data *augDataSource){
    augDataDest->NbPatients = augDataSource->NbPatients;
    augDataDest->T = augDataSource->T;
    memcpy(augDataDest->C, augDataSource->C, augDataDest->NbPatients*sizeof(int));
    memcpy(augDataDest->E, augDataSource->E, augDataDest->NbPatients*sizeof(int));
    memcpy(augDataDest->I0, augDataSource->I0, augDataSource->T*sizeof(int));
    memcpy(augDataDest->I1, augDataSource->I1, augDataSource->T*sizeof(int));
}







/***************** param ******************/
parameters *createParam(){
    parameters *param = (parameters *) malloc(sizeof(parameters));
    if(param == NULL){
	fprintf(stderr, "\n[in: alloc.c->createParam]\nNo memory left for creating parameters. Exiting.\n");
	exit(1);
    }

    param->beta = gsl_matrix_calloc(2,2);
    param->betaWardOut = 0.1;
    param->betaOutOut = 0.1;
    param->Se = 0.0;
    param->Pi = 0.0;
    param->mu = 0.0;
    param->sigma = 0.0;
    param->nu1 = 0.0;
    param->nu2 = 0.0;
    param->tau = 0.0;
    param->alpha = 0.0;
    param->weightNaGen = 0.0;

    return param;
}





void freeParam(parameters *param){
    gsl_matrix_free(param->beta);
    free(param);
}





void copyParam(parameters * paramDest, parameters * paramSource){
    gsl_matrix_memcpy (paramDest->beta,paramSource->beta);
    paramDest->betaWardOut = paramSource->betaWardOut;
    paramDest->betaOutOut = paramSource->betaOutOut;

    /* paramDest->Sp = paramSource->Sp; */
    paramDest->Se = paramSource->Se;

    paramDest->Pi = paramSource->Pi;

    paramDest->mu = paramSource->mu;
    paramDest->sigma = paramSource->sigma;

    paramDest->nu1 = paramSource->nu1;
    paramDest->nu2 = paramSource->nu2;

    paramDest->tau = paramSource->tau;
    paramDest->alpha = paramSource->alpha;
    paramDest->weightNaGen = paramSource->weightNaGen;
}






void print_param(parameters *param){
    printf("\nBeta matrix:\n");
    gsl_matrix_fprintf(stdout, param->beta, "%.3f");
    printf("\nBeta ward<-out: %.3f", param->betaWardOut);
    printf("\nBeta out<-out: %.3f", param->betaOutOut);
    printf("\nsensibility of the test: %.3f", param->Se);
    printf("\nprobability of being colonized at first admission: %.3f", param->Pi);
    printf("\nmu/sigma - duration of colonisation: %.3f %.3f", param->mu, param->sigma);
    printf("\nnu1: %.3f   nu2: %.3f", param->nu1, param->nu2);
    printf("\ntau: %.3f", param->tau);
    printf("\nalpha: %.3f", param->alpha);
    printf("\nweightNaGen: %.3f", param->weightNaGen);
    fflush(stdout);
}





/************ MCMC internals **************/
mcmcInternals *createMcmcInternals(){
    mcmcInternals *MCMCSettings = (mcmcInternals *) malloc(sizeof(mcmcInternals));
    if(MCMCSettings == NULL){
	fprintf(stderr, "\n[in: alloc.c->createMcmcInternals]\nNo memory left for creating MCMCSettings. Exiting.\n");
	exit(1);
    }

    MCMCSettings->Sigma_beta = gsl_matrix_calloc(2,2);
    if(MCMCSettings->Sigma_beta == NULL){
	fprintf(stderr, "\n[in: alloc.c->createMcmcInternals]\nNo memory left for creating MCMCSettings. Exiting.\n");
	exit(1);
    }

    MCMCSettings->NbSimul = 10000;
    MCMCSettings->SubSample = 200;
    MCMCSettings->BurnIn = 5000;
    MCMCSettings->Sigma_betaWardOut = 0.1;
    MCMCSettings->Sigma_betaOutOut = 0.1;
    MCMCSettings->Sigma_mu = 0.1;
    MCMCSettings->Sigma_sigma = 0.1;
    MCMCSettings->Sigma_nu1 = 0.1;
    MCMCSettings->Sigma_nu2 = 0.1;
    MCMCSettings->Sigma_tau = 0.1;
    MCMCSettings->Sigma_alpha = 0.1;

    return MCMCSettings;
}





void printStdProp(mcmcInternals *MCMCSettings){
    int i,j;

    for (i=0;i<2;i++){
	for (j=0;j<2;j++){
	    printf("Std proposal for beta_%d,%d: %lg\n",i,j,gsl_matrix_get(MCMCSettings->Sigma_beta,i,j));
	}
    }
    printf("Std proposal for betaWardOut: %lg\n",MCMCSettings->Sigma_betaWardOut);
    printf("Std proposal for betaOutOut: %lg\n",MCMCSettings->Sigma_betaOutOut);
    printf("Std proposal for mu: %lg\n",MCMCSettings->Sigma_mu);
    printf("Std proposal for sigma: %lg\n",MCMCSettings->Sigma_sigma);
    printf("Std proposal for nu1: %lg\n",MCMCSettings->Sigma_nu1);
    printf("Std proposal for nu2: %lg\n",MCMCSettings->Sigma_nu2);
    printf("Std proposal for tau: %lg\n",MCMCSettings->Sigma_tau);
    printf("Std proposal for alpha: %lg\n",MCMCSettings->Sigma_alpha);

    fflush(stdout);
}





void freeMcmcInternals(mcmcInternals *MCMCSettings){
    gsl_matrix_free(MCMCSettings->Sigma_beta);

    free(MCMCSettings);
}







/************ Acceptance ***************/
acceptance *createAcceptance(){
    acceptance *accept = (acceptance *) malloc(sizeof(acceptance));
    if(accept == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAcceptance]\nNo memory left for creating accept. Exiting.\n");
	exit(1);
    }

    accept->PourcAcc_beta = gsl_matrix_calloc(2,2);
    accept->PourcAcc_betaWardOut=0.0;
    accept->PourcAcc_betaOutOut=0.0;
    accept->PourcAcc_mu=0.0;
    accept->PourcAcc_sigma=0.0;
    accept->PourcAcc_nu1=0.0;
    accept->PourcAcc_nu2=0.0;
    accept->PourcAcc_tau=0.0;
    accept->PourcAcc_alpha=0.0;

    return accept;
}





void reInitiateAcceptance(acceptance *accept){
    int i,j;
    for(i=0;i<2;i++){
	for(j=0;j<2;j++){
	    gsl_matrix_set(accept->PourcAcc_beta,i,j,0.0);
	}
    }

    accept->PourcAcc_betaWardOut=0.0;
    accept->PourcAcc_betaOutOut=0.0;
    accept->PourcAcc_mu=0.0;
    accept->PourcAcc_sigma=0.0;
    accept->PourcAcc_nu1=0.0;
    accept->PourcAcc_nu2=0.0;
    accept->PourcAcc_tau=0.0;
    accept->PourcAcc_alpha=0.0;

}





void printAcceptance(acceptance *accept, NbProposals *NbProp){
    int i,j;

    for(i=0;i<2;i++){
	for(j=0;j<2;j++){
	    printf("Prob accept beta_%d,%d\t%lg\n",i,j,gsl_matrix_get(accept->PourcAcc_beta,i,j)/gsl_matrix_get(NbProp->NbProp_beta,i,j));
	    fflush(stdout);
	}
    }

    printf("Prob accept betaWardOut\t%lg\n",accept->PourcAcc_betaWardOut/NbProp->NbProp_betaWardOut);
    fflush(stdout);
    printf("Prob accept betaOutOut\t%lg\n",accept->PourcAcc_betaOutOut/NbProp->NbProp_betaOutOut);
    fflush(stdout);
    printf("Prob accept mu\t%lg\n",accept->PourcAcc_mu/NbProp->NbProp_mu);
    fflush(stdout);
    printf("Prob accept sigma\t%lg\n",accept->PourcAcc_sigma/NbProp->NbProp_sigma);
    fflush(stdout);
    printf("Prob accept nu1\t%lg\n",accept->PourcAcc_nu1/NbProp->NbProp_nu1);
    fflush(stdout);
    printf("Prob accept nu2\t%lg\n",accept->PourcAcc_nu2/NbProp->NbProp_nu2);
    fflush(stdout);
    printf("Prob accept tau\t%lg\n",accept->PourcAcc_tau/NbProp->NbProp_tau);
    fflush(stdout);
    printf("Prob accept alpha\t%lg\n",accept->PourcAcc_alpha/NbProp->NbProp_alpha);
    fflush(stdout);

}





void freeAcceptance(acceptance *accept){
    gsl_matrix_free(accept->PourcAcc_beta);
    free(accept);
}






/************* Is Acceptance OK *************/
isAcceptOK *createIsAcceptOK(){
    isAcceptOK *acceptOK = (isAcceptOK *) malloc(sizeof(isAcceptOK));
    if(acceptOK == NULL){
	fprintf(stderr, "\n[in: alloc.c->createIsAcceptOK]\nNo memory left for creating acceptOK. Exiting.\n");
	exit(1);
    }

    acceptOK->IsAccOK_beta = gsl_matrix_calloc(2,2);

    acceptOK->IsAccOK_betaWardOut=0.0;
    acceptOK->IsAccOK_betaOutOut=0.0;
    acceptOK->IsAccOK_mu=0.0;
    acceptOK->IsAccOK_sigma=0.0;
    acceptOK->IsAccOK_nu1=0.0;
    acceptOK->IsAccOK_nu2=0.0;
    acceptOK->IsAccOK_tau=0.0;
    acceptOK->IsAccOK_alpha=0.0;

    return acceptOK;
}





void freeIsAcceptOK(isAcceptOK *acceptOK){
    gsl_matrix_free(acceptOK->IsAccOK_beta);

    free(acceptOK);
}








/************ NbProposals ***************/
NbProposals *createNbProposals(){
    NbProposals *NbProp = (NbProposals *) malloc(sizeof(NbProposals));
    if(NbProp == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbProposals]\nNo memory left for creating NbProp. Exiting.\n");
	exit(1);
    }

    NbProp->NbProp_beta = gsl_matrix_calloc(2,2);

    NbProp->NbProp_betaWardOut=0.0;
    NbProp->NbProp_betaOutOut=0.0;
    NbProp->NbProp_mu=0.0;
    NbProp->NbProp_sigma=0.0;
    NbProp->NbProp_nu1=0.0;
    NbProp->NbProp_nu2=0.0;
    NbProp->NbProp_tau=0.0;
    NbProp->NbProp_alpha=0.0;

    return NbProp;
}





void reInitiateNbProp(NbProposals * NbProp){
    int i,j;

    for(i=0 ; i<2 ; i++)
	{
	    for(j=0 ; j<2 ; j++)
		{
		    gsl_matrix_set(NbProp->NbProp_beta,i,j,0.0);
		}
	}

    NbProp->NbProp_betaWardOut=0.0;
    NbProp->NbProp_betaOutOut=0.0;
    NbProp->NbProp_mu=0.0;
    NbProp->NbProp_sigma=0.0;
    NbProp->NbProp_nu1=0.0;
    NbProp->NbProp_nu2=0.0;
    NbProp->NbProp_tau=0.0;
    NbProp->NbProp_alpha=0.0;

}





void freeNbProposals(NbProposals *NbProp){
    gsl_matrix_free(NbProp->NbProp_beta);

    free(NbProp);
}







/************** OUTPUT FILES **************/
output_files *createFILES(char *workspace){
    output_files *fich = (output_files *) malloc(sizeof(output_files));
    if(fich == NULL){
	fprintf(stderr, "\n[in: alloc.c->createFILES]\nNo memory left for creating fich. Exiting.\n");
	exit(1);
    }

    char fileName[300];

    strcpy(fileName, workspace);
    strcat(fileName,"LogL.txt");
    fich->LogL = fopen(fileName,"w");
    if ( fich->LogL == NULL )
	{
	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
	    exit(1);
	}

    strcpy(fileName, workspace);
    strcat(fileName,"ColonDates.txt");
    fich->ColonDates = fopen(fileName,"w");
    if ( fich->ColonDates == NULL )
	{
	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
	    exit(1);
	}

    strcpy(fileName, workspace);
    strcat(fileName,"EndColonDates.txt");
    fich->EndColonDates = fopen(fileName,"w");
    if ( fich->EndColonDates == NULL )
	{
	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
	    exit(1);
	}

    strcpy(fileName, workspace);
    strcat(fileName,"Parameters.txt");
    fich->Parameters = fopen(fileName,"w");
    if (fich->Parameters == NULL )
	{
	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
	    exit(1);
	}

    return fich;
}






void freeFILES(output_files *fich){
    fclose(fich->LogL);
    fclose(fich->ColonDates);
    fclose(fich->EndColonDates);
    fclose(fich->Parameters);

    free(fich);
}