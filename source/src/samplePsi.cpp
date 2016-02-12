/*Determine non-parametric (ranked) correlation
Generate random correlated matrices and determine their psi
*/

#if defined _WIN64 || defined _WIN32
	#include <windows.h>
#else
	#include <pthread.h>
#endif 
#include <deque>
#include <vector>

#include <limits>
#include <iostream>
#include "getPsi.h"

#include <algorithm>    // std::min_element, std::max_element
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <random>
using namespace std;

void * ThreadFunc(void * pUserData);

void getMatrix(double* mat,int b,int n,unsigned int* seed){
	//Generate a matrix with random numbers
	std::default_random_engine generator (*seed);
	std::normal_distribution<double> distribution(0.0,1.0);
	for (int i = 0; i < (b*n); i++){
		mat[i] = distribution(generator);
	}
}

void matmultiply(double* mat1,double* mat2,double* outputmat,int b,int n){
	//mat1 is assumed to be a nxb matrix
	//mat2 is assumed to be a square bxb matrix
	//outputmat is assumed to be a nxb matrix to 
	int rowIdx1=0,colIdx2=0,z=0;
	for (rowIdx1 = 0; rowIdx1 < n; rowIdx1++){
	for (colIdx2 = 0; colIdx2 < b; colIdx2++){
		outputmat[rowIdx1+colIdx2*n] = 0;	
	for (z = 0; z < b; z++){
		outputmat[rowIdx1+colIdx2*n] += mat1[z*n+rowIdx1]*mat2[z+colIdx2*b]	;
	}
				 	
	}}
}




//----------------------------------------------------------------
//Job class
//----------------------------------------------------------------

class CJob
{
    private:
	
    public:
		//Whatever you need inside the worker thread for this job
		unsigned int id;
		unsigned int* seeds;
		unsigned int nDraws;
		unsigned int offset;
		
		CJob(){}

};

//----------------------------------------------------------------
//Jobs class
//----------------------------------------------------------------

class CJobs
{
    private:
    public:
		double* pd_choleski;
		int* pi_missingmat1;
		int* pi_missingmat2;
		int nDraws;
		int maxB1;
		int maxB2;
		int n1;
		int n2;
		int maxn;
		double* pd_result1;
		double* pd_result2;
		int counter;


	
		#if defined _WIN64 || defined _WIN32
			CRITICAL_SECTION criticalSection; //for windows
		#else
	        pthread_mutex_t mutex;
		#endif 

        deque<CJob>     queue;
        
		//Whatever you need inside the worker threads for all jobs

		CJobs (double*& pd_choleski, int*& pi_missingmat1, int*& pi_missingmat2, unsigned int& nDraws,int& maxB1,int& maxB2, int& n1,int& n2,  double* pd_result1, double* pd_result2)
		{
			#if defined _WIN64 || defined _WIN32
				InitializeCriticalSection(&criticalSection); //for windows 
			#else
	            pthread_mutex_init(&mutex,NULL);
			#endif

			this->pd_choleski = pd_choleski;
			this->pi_missingmat1 = pi_missingmat1;
			this->pi_missingmat2 = pi_missingmat2;
			this->nDraws = nDraws;
			this->maxB1 = maxB1;
			this->maxB2 = maxB2;
			this->n1 = n1;
			this->n2 = n2;
			this->pd_result1 = pd_result1;
			this->pd_result2 = pd_result2;
			this->counter =0;
			this->maxn = max(n1,n2);
        }
        
        ~CJobs(void)
        {
			#if defined _WIN64 || defined _WIN32
				DeleteCriticalSection(&criticalSection); //for windows            
			#else
				pthread_mutex_destroy(&mutex);
			#endif
        }
};





//----------------------------------------------------------------
//Worker thread
//----------------------------------------------------------------

//For windows
#if defined _WIN64 || defined _WIN32
long unsigned int WINAPI ThreadFuncWin(void * pUserData)
{
	ThreadFunc(pUserData);
	return 0;
} 
#endif
void * ThreadFunc(void * pUserData) 
{
	CJob  job;
    CJobs * pJobs;	

	pJobs=(CJobs*)pUserData;
 
	double *dataMatrix1 = (double*)malloc(sizeof(double) * (pJobs->maxB1+pJobs->maxB2)*pJobs->maxn);	
	double *dataMatrix2 = (double*)malloc(sizeof(double) * (pJobs->maxB1+pJobs->maxB2)*pJobs->maxn);

	double *dataMatrix3 = (double*)malloc(sizeof(double) * pJobs->maxB1*pJobs->n1);	
	double *dataMatrix4 = (double*)malloc(sizeof(double) * pJobs->maxB2*pJobs->n2);


	for(;;)
    {
        //----------------------------------------------------------------
        //Get job from queue
        //----------------------------------------------------------------
    
		#if defined _WIN64 || defined _WIN32
			EnterCriticalSection(&pJobs->criticalSection);    //for windows
		#else
		    pthread_mutex_lock(&pJobs->mutex);
		#endif	
        if(pJobs->queue.size()==0)
        {
			
			#if defined _WIN64 || defined _WIN32
				LeaveCriticalSection(&pJobs->criticalSection);    //for windows		
			#else
	            pthread_mutex_unlock(&pJobs->mutex);
			#endif
			free(dataMatrix1);
			free(dataMatrix2);
			free(dataMatrix3);
			free(dataMatrix4);
 
			return (void*)pJobs;
        }
        
        job=pJobs->queue.front();
        pJobs->queue.pop_front();
	
		#if defined _WIN64 || defined _WIN32
			LeaveCriticalSection(&pJobs->criticalSection);    //for windows
		#else
	        pthread_mutex_unlock(&pJobs->mutex);
		#endif

		
        //----------------------------------------------------------------
		//Do work on job
        //----------------------------------------------------------------
		unsigned int cycle;		
		for (cycle = 0; cycle< job.nDraws; cycle++){
			//if (bSigintReceived) { break;}
			 	
		
			getMatrix(dataMatrix2,pJobs->maxB1+pJobs->maxB2,pJobs->maxn,job.seeds+cycle);
			matmultiply(dataMatrix2,pJobs->pd_choleski,dataMatrix1,pJobs->maxB1+pJobs->maxB2,pJobs->maxn);

			//cpy to datamatrix3 to datamatrix1 for missingmat1;
			for (int i = 0; i < pJobs->n1; i++) {
				for (int j = 0; j < pJobs->maxB1; j++) {
					if (pJobs->pi_missingmat1[j*pJobs->n1+i]==1) {
						dataMatrix3[ j*pJobs->n1+i] = (dataMatrix1[ j*pJobs->maxn+i]);
					} else {
						dataMatrix3[ j*pJobs->n1+i] = std::numeric_limits<double>::infinity();
					}
				}
			}
			//make reference between datamatrix4 to datamatrix1 for missingmat2;
			for (int i = 0; i < pJobs->n2; i++) {
				for (int j = 0; j < pJobs->maxB2; j++) {
					if (pJobs->pi_missingmat2[j*pJobs->n2+i]==1) {
						dataMatrix4[ j*pJobs->n2+i] = (dataMatrix1[ (pJobs->maxB1+j)*pJobs->maxn+i]);
					} else {
						dataMatrix4[ j*pJobs->n2+i] = std::numeric_limits<double>::infinity();
					}
				}
			}
		
			pJobs->pd_result1[job.offset +cycle]  = getPsi(dataMatrix3, pJobs->n1, pJobs->maxB1);
			pJobs->pd_result2[job.offset +cycle]  = getPsi(dataMatrix4, pJobs->n2, pJobs->maxB2);
			
			
			#if defined _WIN64 || defined _WIN32
				EnterCriticalSection(&pJobs->criticalSection);    //for windows
			#else
			    pthread_mutex_lock(&pJobs->mutex);
			#endif
				pJobs->counter ++;
			#if defined _WIN64 || defined _WIN32
				LeaveCriticalSection(&pJobs->criticalSection);    //for windows
			#else
			    pthread_mutex_unlock(&pJobs->mutex);
			#endif
			
		}

	}
}

  

void startMultithreadedSampling(double* pd_choleski,int* pi_missingmat1,int* pi_missingmat2,unsigned int nDraws,int maxB1,int maxB2, int n1, int n2,unsigned int nCPU,double* pd_result1, double* pd_result2, unsigned int seed) {
	
	unsigned long i;
	unsigned int iThread;
	
	unsigned int* seeds = (unsigned int*)malloc(sizeof(unsigned int) * nDraws);  

	std::default_random_engine generator (seed);
	std::uniform_int_distribution<unsigned int> distribution(0);	
	for (unsigned i = 0 ; i < nDraws; i++){
		seeds[i] =  distribution(generator);
	}

	
	CJobs jobs(pd_choleski, pi_missingmat1, pi_missingmat2, nDraws ,maxB1,maxB2, n1,n2, pd_result1, pd_result2);	
	#if defined _WIN64 || defined _WIN32
		HANDLE *	pThreads;
		 pThreads=(HANDLE*)malloc(nCPU*sizeof(HANDLE));
	#else
		pthread_t *     pThreads;   
		pThreads=(pthread_t*)malloc(nCPU*sizeof(pthread_t));

	#endif

	//----------------------------------------------------------------
	//Prepare jobs
	//----------------------------------------------------------------
	CJob* job = (CJob*)malloc(nCPU*sizeof(CJob));
	int offset = 0;
	for (i = 0; i < nCPU; i++){
		job[i].id=i;
		job[i].seeds = seeds+offset;
		job[i].nDraws = (nDraws-offset)/(nCPU-i);
		job[i].offset = offset;
		offset +=  	job[i].nDraws;
		jobs.queue.push_back(job[i]);
	}
	//----------------------------------------------------------------
	//Create threads
	//----------------------------------------------------------------


	for(iThread=0;iThread<nCPU;iThread++)
	{
		#if defined _WIN64 || defined _WIN32
			pThreads[iThread]=CreateThread(NULL,0, ThreadFuncWin,(void*)&jobs,0, NULL); //for windows
		#else 
			pthread_create(&pThreads[iThread],NULL, ThreadFunc, (void*)&jobs);
		#endif
	}
	//----------------------------------------------------------------
	//Wait for threads to finish
	//----------------------------------------------------------------
	for(iThread=0;iThread<(nCPU);iThread++)
	{
		#if defined _WIN64 || defined _WIN32
		while(WaitForMultipleObjects(nCPU,pThreads,TRUE,100)){	
			/*if (this->verbose){			
				Rprintf("\rEstimating variance: %.1f %%        ",(100.0*jobs.counter)/(*pncycles));
				R_FlushConsole();
				R_ProcessEvents();
			}*/					
		} //for windows
		#else
		while (pthread_tryjoin_np(pThreads[iThread],0)) {
			//nanosleep((const struct timespec[]){{0, 100000000}},(struct timespec[]){{0, 100000000}}); //bug:error: taking address of temporary array
			nanosleep((const struct timespec[]){{0, 100000000}},NULL); 
			/*if (this->verbose){
				Rprintf("\rEstimating variance: %.1f \%           ",(100.0*jobs.counter)/(*pncycles));
				R_FlushConsole();
				R_ProcessEvents();
			}*/
		}
		#endif

	}

	
	/*if (this->verbose and !bSigintReceived){
		Rprintf("\rEstimating variance: %.1f %%   \n",100.0);
		R_FlushConsole();
		R_ProcessEvents();
	}*/

	free(job);
	free(pThreads);
	free(seeds);
}

extern "C" {
	SEXP samplePsi(SEXP rP_choleski,SEXP rP_mat1Missing,SEXP rP_mat2Missing,SEXP rP_nDraws,SEXP rP_nCPU, SEXP rP_seed){
		//rP_choleski: A Cholseki matrix of bmax*bmax
		//rP_bn: A vector of numbers of observations per subject
		//rP_nDraws: An integer for the number of matrices to draw
		rP_choleski = coerceVector(rP_choleski, REALSXP);
		rP_mat1Missing = coerceVector(rP_mat1Missing, LGLSXP);
		rP_mat2Missing = coerceVector(rP_mat2Missing, LGLSXP);
		rP_nDraws = coerceVector(rP_nDraws, INTSXP);
		rP_nCPU = coerceVector(rP_nCPU, INTSXP);
		/*check input parameters*/

		SEXP Rdim1,Rdim2,Rdim3;
		PROTECT(Rdim1=getAttrib(rP_choleski,R_DimSymbol));
		PROTECT(Rdim2=getAttrib(rP_mat1Missing,R_DimSymbol));
		PROTECT(Rdim3=getAttrib(rP_mat2Missing,R_DimSymbol));
		int I=INTEGER(Rdim1)[0];
		int J=INTEGER(Rdim1)[1];
		if (rP_mat1Missing==R_NilValue){error("No missing matrix given");}
		if (Rdim2==R_NilValue){error("No missing matrix given");}
		if (rP_mat2Missing!=R_NilValue and Rdim3==R_NilValue){error("missing matrix 2 is given but not a matrix");}
 
		int n1 = INTEGER(Rdim2)[0];
		int n2 = 0;
		if (rP_mat2Missing!=R_NilValue){n2 = INTEGER(Rdim3)[0];}
	

		int maxB1 = INTEGER(Rdim2)[1];//*std::max_element( INTEGER(rP_bn1),INTEGER(rP_bn1)+n  );
		int maxB2 = 0;
		if (n2>0){ maxB2 = INTEGER(Rdim3)[1];}//*std::max_element( INTEGER(rP_bn2),INTEGER(rP_bn2)+n  ); }
		int maxB = maxB1+maxB2;
		if (I!=J) {error("Input matrix must be a square matrix");}
		if (I != maxB) {error("Cholseki matrix dimensions must satisfy bmax*bmax");}
		
		int nCPU = *INTEGER(rP_nCPU);
		 
		SEXP output1;
		PROTECT(output1 = allocMatrix(REALSXP, *INTEGER(rP_nDraws),2));

		startMultithreadedSampling(REAL(rP_choleski),LOGICAL(rP_mat1Missing),LOGICAL(rP_mat2Missing),*INTEGER(rP_nDraws),maxB1,maxB2, n1, n2, nCPU, REAL(output1),REAL(output1)+*INTEGER(rP_nDraws),INTEGER(rP_seed)[0]); 
		
		UNPROTECT(4);
		return output1;
	}

}


/*
 
library(concordance)
r<-0.9999
chol<-matrix(r,ncol=3,nrow=3);diag(chol)<-1;chol<-chol(chol)
bn<-as.integer(rep(3,200));
i<-.Call("samplePsi",chol,bn,100)
*/	
