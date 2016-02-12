
#if defined _WIN64 || defined _WIN32
	#include <windows.h>
#else
	#include <pthread.h>
#endif 

#include <iostream>
#include <fstream>
#include <signal.h>     /* signal, raise, sig_atomic_t */


#include <errno.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <iostream>
#include <vector>
#include <map>
#include <time.h> 

#include <algorithm>

#include <deque>
#include "getPsi.h"
#include <random>

#define  MAX_LEVELS  100000

using namespace std;

bool bSigintReceived = false;
int p;

void sigintHandler(int signal){
	bSigintReceived = true;
}

void * ThreadFunc_bootstrap(void * pUserData);
	
void countingsort203(unsigned int* aiIn, unsigned int nInItems,unsigned int iLo, unsigned int iHi){
	// maak een precies passende turflijst, begin met 0 voor elk getal 
	unsigned  int nTurfItems= 1 + iHi - iLo;
	unsigned int* aiTurf = (unsigned int*)malloc(sizeof(unsigned int) * nTurfItems);
	unsigned  int j,i;
	for (j= 0; j < nTurfItems; j++){
		aiTurf[j]= 0;
	}
	// turven 
	for (i= 0; i < nInItems; i++){
	aiTurf[aiIn[i] - iLo]++; // 1 erbij voor het getal op positie i
	}
	// turflijst langslopen en aftellen, we hergebruiken de gegeven lijst voor het antwoord  
	i= 0;
	for (j= 0; j < nTurfItems; j++){
		while (aiTurf[j] > 0){
		  aiIn[i++]= iLo + j;
		  aiTurf[j]--;
		}
	}
	free(aiTurf);
	// aiIn bevat nu de getallen op de juiste volgorde  
}


class BootstrapperDataClass: public DataClass
{ 
	private:
	unsigned int* bootstrapSelection;
	double* orgS;

	public:
	BootstrapperDataClass(double* mat, int  nrow1, int ncol1):  DataClass(mat, nrow1,  ncol1){
		this->orgS = new double[ nrow1 *  ncol1];
		this->bootstrapSelection = new unsigned int[ ncol1 ];

		for (int i = 0; i < (nrow1 *  ncol1); i++){
			this->orgS[i] = mat[i];		
		}
	};
	BootstrapperDataClass (const BootstrapperDataClass &obj):  DataClass(obj){
		unsigned int i;
		this->orgS = new double[ obj.ncol * obj.nrow ];
		this->bootstrapSelection = new unsigned int[ obj.ncol ];
		for (i = 0; i < (obj.ncol*obj.nrow); i++){
			this->orgS[i] = obj.orgS[i];		
		}
		for (i = 0; i < obj.ncol; i++){
			this->bootstrapSelection[i] = obj.bootstrapSelection[i];		
		}
	};
	~BootstrapperDataClass(void){
		delete [] orgS; 
		delete [] bootstrapSelection; 
	}

	void doBootstrap(unsigned int seed){
			unsigned long i,j,count;
			//Random selection of observations
			
			std::default_random_engine generator (seed);
			std::uniform_int_distribution<unsigned int> distribution(0,(this->ncol)-1);
					
			for (i = 0; i < this->ncol; i++) { this->bootstrapSelection[i] = distribution(generator); 

			}
			//countingsort203(this->bootstrapSelection, this->ncol,(unsigned int) 0, (unsigned int) (this->ncol-1) ); //Sort from low to high
			count =0;			
			for (i = 0; i < this->ncol; i++) { 
				for (j = 0; j < this->nrow; j++) { 			
					this->sMAT[count++] =  this->orgS[ this->bootstrapSelection[i]* this->nrow+j];//+0.001* rand_r( &seed )/(1.0*RAND_MAX) ;
			}} 
			
			this->preprocess();
			
	};

};



//----------------------------------------------------------------
//Job class
//----------------------------------------------------------------

class CJob_bootstrap
{
    private:
	
    public:
		//Whatever you need inside the worker thread for this job
		unsigned int id;
		unsigned int* seeds;
		unsigned int offset;
		unsigned int nBootstraps;
		
		CJob_bootstrap(unsigned int id,unsigned int* seeds): id(id), seeds(seeds){}
		CJob_bootstrap(void){}

};

//----------------------------------------------------------------
//Jobs class
//----------------------------------------------------------------

class CJobs_bootstrap
{
    private:
    public:
		unsigned long counter;
		double* pdOutput;
		BootstrapperDataClass* dc1 ;
		BootstrapperDataClass* dc2 ;

		double* mat1, *mat2;
		unsigned long nrow1, nrow2, ncol1, ncol2;
		
		#if defined _WIN64 || defined _WIN32
			CRITICAL_SECTION criticalSection; //for windows
		#else
	        pthread_mutex_t mutex;
		#endif 

        deque<CJob_bootstrap>     queue;
        
		//Whatever you need inside the worker threads for all jobs
 
		CJobs_bootstrap(BootstrapperDataClass* dc1, BootstrapperDataClass* dc2, double* pdOutput): pdOutput(pdOutput), dc1(dc1), dc2(dc2)
        {
			#if defined _WIN64 || defined _WIN32
				InitializeCriticalSection(&criticalSection); //for windows 
			#else
	            pthread_mutex_init(&mutex,NULL);
			#endif

        }
        
        ~CJobs_bootstrap(void)
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
long unsigned int WINAPI ThreadFuncWin_bootstrap(void * pUserData)
{
	ThreadFunc_bootstrap(pUserData);
	return 0;
} 
#endif

void * ThreadFunc_bootstrap(void * pUserData) 
{
	
	CJob_bootstrap job;
    CJobs_bootstrap * pJobs;	

	pJobs=(CJobs_bootstrap*)pUserData;

	unsigned int cycle;
 

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
		for (cycle = 0; cycle< job.nBootstraps; cycle++){
			if (bSigintReceived) { break;}
			BootstrapperDataClass dc1 = BootstrapperDataClass(*(pJobs->dc1)); //call copy constructor
			BootstrapperDataClass dc2 = BootstrapperDataClass(*(pJobs->dc2)); //call copy constructor
		
			dc1.doBootstrap(job.seeds[cycle]);
			dc2.doBootstrap(job.seeds[cycle]);

			
			
			pJobs->pdOutput[job.offset+cycle] = dc1.calculatePSI()-dc2.calculatePSI();
				
			/*
			#if defined _WIN64 || defined _WIN32
				EnterCriticalSection(&pJobs->criticalSection);    //for windows
			#else
			    pthread_mutex_lock(&pJobs->mutex);
			#endif
				
			//pJobs->counter ++;
			#if defined _WIN64 || defined _WIN32
				LeaveCriticalSection(&pJobs->criticalSection);    //for windows
			#else
			    pthread_mutex_unlock(&pJobs->mutex);
			#endif*/
			
		}
	}
}

//----------------------------------------------------------------
//BootstrapMain class
//----------------------------------------------------------------
class Bootstrapper
{
    private:
	
		double* mat1;
		double* mat2;
		unsigned int nrow1,nrow2,ncol1,ncol2;		

		void (*prev_handler)(int);

		unsigned int ncycles;
		unsigned int nCPU;


    public:
		Bootstrapper(SEXP& MAT1, SEXP& MAT2, unsigned int ui_nDraws,unsigned int ui_nCPU){
			// extract the "dim" attribute 
			SEXP dim1 = getAttrib( MAT1, R_DimSymbol ) ;
			SEXP dim2 = getAttrib( MAT2, R_DimSymbol ) ;
			this->nrow1 = INTEGER(dim1)[0]; 
			this->ncol1 = INTEGER(dim1)[1]; 
			this->nrow2 = INTEGER(dim2)[0]; 
			this->ncol2 = INTEGER(dim2)[1]; 

			if (nrow1 !=nrow2 or ncol1!=ncol2) { error("Matrices not of comparable size!");}
			this->mat1 = REAL(MAT1);
			this->mat2 = REAL(MAT2);	
			this->ncycles = ui_nDraws;
			this->nCPU = ui_nCPU;	
			
			this->prev_handler =signal (SIGINT, sigintHandler);
		}

		~Bootstrapper(){
			signal (SIGINT, this->prev_handler);
			if (bSigintReceived){
				bSigintReceived = false;	
				raise(SIGINT);
			}
			R_CheckUserInterrupt();
		}
 
		void run(double *pdResult, unsigned int seed){
			unsigned long i;
			unsigned int iThread;
			unsigned int* seeds = (unsigned int*)malloc(sizeof(unsigned int) * this->ncycles);  

			std::default_random_engine generator (seed);
			std::uniform_int_distribution<unsigned int> distribution(0);	



			for (i = 0 ; i < this->ncycles; i++){
				seeds[i] = distribution(generator); //each job has its own seed determined in this master thread.
			}
			
			//----------------------------------------------------------------
			//Begin bootstrap
			//----------------------------------------------------------------
			BootstrapperDataClass dc1(this->mat1, this->nrow1, this->ncol1);
			BootstrapperDataClass dc2(this->mat2, this->nrow2, this->ncol2);
			
			CJobs_bootstrap jobs(&dc1, &dc2, pdResult);	
			#if defined _WIN64 || defined _WIN32
				HANDLE *	pThreads;
				 pThreads=(HANDLE*)malloc(this->nCPU*sizeof(HANDLE));
			#else
				pthread_t *     pThreads;   
				pThreads=(pthread_t*)malloc(this->nCPU*sizeof(pthread_t));

			#endif
		
			//----------------------------------------------------------------
			//Prepare jobs
			//----------------------------------------------------------------
			CJob_bootstrap* job = (CJob_bootstrap*)malloc(this->nCPU*sizeof(CJob_bootstrap));
			unsigned long offset = 0;
			unsigned long nBoostraps = this->ncycles;
			for (i = 0; i < this->nCPU; i++){

				job[i]  =  CJob_bootstrap(i,seeds+offset);
				job[i].offset = offset;
				job[i].nBootstraps = nBoostraps/(this->nCPU-i);
				offset+=job[i].nBootstraps;
				nBoostraps-=job[i].nBootstraps;
				jobs.queue.push_back(job[i]);
			}
			//----------------------------------------------------------------
			//Create threads
			//----------------------------------------------------------------
		
				
			for(iThread=0;iThread < this->nCPU;iThread++)
			{
				#if defined _WIN64 || defined _WIN32
					pThreads[iThread]=CreateThread(NULL,0, ThreadFuncWin_bootstrap,(void*)&jobs,0, NULL); //for windows
				#else 
					pthread_create(&pThreads[iThread],NULL, ThreadFunc_bootstrap, (void*)&jobs);
				#endif
			}
			//----------------------------------------------------------------
			//Wait for threads to finish
			//----------------------------------------------------------------

			for(iThread=0;iThread< this->nCPU;iThread++)
			{
				#if defined _WIN64 || defined _WIN32
					while(WaitForMultipleObjects(this->nCPU,pThreads,TRUE,100)){	
			
					} //for windows
				 #else
					while (pthread_tryjoin_np(pThreads[iThread],0)) {
						nanosleep((const struct timespec[]){{0, 100000000}},NULL); 
				}
				#endif

			}

			
			
			 free(job);
			free(pThreads);
			free(seeds);
	
		}
};
//----------------------------------------------------------------


extern "C"{

SEXP doBootstrap203(SEXP MAT1,SEXP MAT2, SEXP ncycles, SEXP r_nCPU, SEXP R_seed ){	
	SEXP result = PROTECT(allocVector(REALSXP, *INTEGER(ncycles) ));
	Bootstrapper bootstr(MAT1, MAT2, *INTEGER(ncycles), *INTEGER(r_nCPU));
	bootstr.run(REAL(result),*INTEGER(R_seed));
	UNPROTECT(1);

	return(result);
}}


