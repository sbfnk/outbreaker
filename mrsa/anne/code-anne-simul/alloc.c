#include "common.h"
#include "alloc.h"
#include "matvec.h"

/* Si tu lis ca Thibaut c'est que je masterise git ;-))) */



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
		/* if(data->A[i]!=NULL) gsl_vector_fprintf(stdout, data->A[i], "%.1f"); else printf("NULL\n"); */
		if(data->A[i]!=NULL) print_gsl_vector(data->A[i], "%.1f "); else printf("NULL\n");

	}
	printf("\nDischarge dates:\n");
	for(i=0;i<data->NbPatients;i++){
		printf("\nPatient %d:\n",i);
		if(data->D[i]!=NULL) print_gsl_vector(data->D[i], "%.1f "); else printf("NULL\n");
	}
	printf("\nDates of positive swabs:\n");
	for(i=0;i<data->NbPatients;i++){
		printf("Patient %d:\n",i);
		if(data->P[i]!=NULL) print_gsl_vector(data->P[i], "%.1f "); else printf("NULL\n");
	}
	printf("\nDates of negative swabs:\n");
	for(i=0;i<data->NbPatients;i++){
		printf("Patient %d:\n",i);
		if(data->N[i]!=NULL) print_gsl_vector(data->N[i], "%.1f "); else printf("NULL\n");
	}
	printf("\nHospital presence:\n");
	for(i=0;i<data->NbPatients;i++){
		printf("Patient %d:\n",i);
		print_gsl_vector(data->IsInHosp[i], "%.0f ");
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


void readParameters(char* workspace, parameters * param, hospDurationParam *paramHosp)
{
    FILE * paramInit;
	char val[30];
	char file[200];

	strcpy(file, workspace);
	strcat(file,"param.txt");
	if ((paramInit=fopen(file,"r"))==NULL)
	{
		printf("Cannot read param.txt");
		exit(2);
	}

	fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,0,0));
	printf("%s %g\n",val,gsl_matrix_get(param->beta,0,0));

	fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,0,1));
	printf("%s %g\n",val,gsl_matrix_get(param->beta,0,1));

	fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,1,0));
	printf("%s %g\n",val,gsl_matrix_get(param->beta,1,0));

	fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,1,1));
	printf("%s %g\n",val,gsl_matrix_get(param->beta,1,1));

	fscanf(paramInit,"%s %lf",val,&param->betaWardOut);
	printf("%s %g\n",val,param->betaWardOut);

	fscanf(paramInit,"%s %lf",val,&param->betaOutOut);
	printf("%s %g\n",val,param->betaOutOut);

	/*fscanf(paramInit,"%s %lf",val,&param->Sp);
	printf("%s %g\n",val,param->Sp);*/

	fscanf(paramInit,"%s %lf",val,&param->Se);
	printf("%s %g\n",val,param->Se);

	fscanf(paramInit,"%s %lf",val,&param->Pi);
	printf("%s %g\n",val,param->Pi);

	fscanf(paramInit,"%s %lf",val,&paramHosp->mu);
	printf("%s %g\n",val,paramHosp->mu);

	fscanf(paramInit,"%s %lf",val,&paramHosp->sigma);
	printf("%s %g\n",val,paramHosp->sigma);

	fscanf(paramInit,"%s %lf",val,&param->mu);
	printf("%s %g\n",val,param->mu);

	fscanf(paramInit,"%s %lf",val,&param->sigma);
	printf("%s %g\n",val,param->sigma);

	fscanf(paramInit,"%s %lf",val,&param->nu1);
	printf("%s %g\n",val,param->nu1);

	fscanf(paramInit,"%s %lf",val,&param->nu2);
	printf("%s %g\n",val,param->nu2);

	fscanf(paramInit,"%s %lf",val,&param->tau);
	printf("%s %g\n",val,param->tau);

	fscanf(paramInit,"%s %lf",val,&param->alpha);
	printf("%s %g\n",val,param->alpha);

    fclose(paramInit);

}

/***************** hospDurationParam ******************/
hospDurationParam *createHospDurationParam(){
	hospDurationParam *param = (hospDurationParam *) malloc(sizeof(hospDurationParam));
	if(param == NULL){
		fprintf(stderr, "\n[in: alloc.c->createHospDurationParam]\nNo memory left for creating hospDurationParam. Exiting.\n");
		exit(1);
	}

	param->mu = 0.0;
	param->sigma = 0.0;

	return param;
}


void print_HospDurationParam(hospDurationParam *param){
	printf("\nmu/sigma - duration of hospitalisation: %.3f %.3f", param->mu, param->sigma);
	fflush(stdout);
}



void freeHospDurationParam(hospDurationParam *in){
    free(in);
}