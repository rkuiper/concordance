#include <iostream>
#include <vector>
#include <map>
#include <time.h> 

#include <algorithm>

#include <deque>
#include <limits>
 
#include <stdlib.h>

#include "getPsi.h"


using namespace std;


typedef struct
{
    long    lPosition;
    double  dValue;
    
}lKEY_dVALUE_PAIR;

typedef struct
{
    long    lPosition;
    long	lValue;
    
}lKEY_lVALUE_PAIR;



struct PREDICATE_K_V_PAIR
{
   bool operator() (lKEY_dVALUE_PAIR pair1,lKEY_dVALUE_PAIR pair2)
    {
	     return pair1.dValue<pair2.dValue;
    }
	bool operator() (lKEY_lVALUE_PAIR pair1,lKEY_lVALUE_PAIR pair2)
    {
		return pair1.lValue<pair2.lValue;
    }
} predicate_k_v_pair;



template <typename T>
inline T  transpose (T mat, int const& nrow,int const& ncol) 
{ 

	int n = nrow*ncol;	
	int idx1=0,idx2 = 0;
	int row=0,col=0;
	T placeholder = static_cast<T>(malloc(n * sizeof(mat[0]))) ;
	for (idx1 = 0; idx1 < n;idx1++){
		row = (idx1)%nrow+1;
		col = (idx1)/nrow+1;
		idx2 = (row-1)*ncol+col-1;
		placeholder[idx2] = mat[idx1];
	}
	for (idx1 = 0; idx1 < n;idx1++){
		mat[idx1] = placeholder[idx1];
	}
	delete [] placeholder;
	return mat; 
} 




DataClass::~DataClass(){
	//cout << "CALLED\n";
	if (this->BN!=NULL){ free(this->BN);}

	if (this->sMAT!=NULL){ free(this->sMAT);}
	if (this->qMAT!=NULL){ free(this->qMAT);}
	if (this->rMAT!=NULL){ free(this->rMAT);}

	this->BN=NULL;
	this->sMAT=NULL;
	this->qMAT=NULL;
	this->rMAT=NULL;
}
		 	
double DataClass::calculatePSI(void ){
	unsigned long *QMAT = this->qMAT;
	unsigned int* bn = this->BN ;
	unsigned int n = this->ncol;
	unsigned int bmax = this->nrow;
	unsigned long T = this->T;

	unsigned long k,j;
	
	//cout << "nrow = " << n <<"; ncol = "<< bmax<<endl;
	//Calculate psi
	long theta;
	double* pdTempPsi = (double*)malloc(sizeof(double) * n);
	for (j = 0; j < n; j++){pdTempPsi[j] = 0;}

	for (j = 0; j < n; j++){

		double denom = 	(bn[j]-1.0)*T*(T-bn[j]);// (b-1)*T*(T-b)
		for (k = 1; k < bn[j]; k++){		
			theta = -2 * k*(k-bn[j]) ; // -2*k*(k-b) 
			pdTempPsi[j]+= theta * (QMAT[j*bmax+(k-1)]);
		}
		pdTempPsi[j]/=denom;
	}
	double psi = 1.0;


	for (j = 0; j < n; j++){psi -= pdTempPsi[j];}

	free(pdTempPsi);
	return(psi);
}


void DataClass::R2Q(){
	//Convert R matrix to Q matrix (i.e. Q_1 = R_2-R_1-1):
	unsigned int i=0,j=0;
	for (j = 0; j < this->ncol; j++){
		for (i = 0; i < (this->BN[j]-1); i++){
			this->qMAT[j*this->nrow+i] = this->rMAT[j*this->nrow+i+1]-this->rMAT[j*this->nrow+i] -1;
		}
	}
}


void DataClass::S2R(){
	//Convert to int
	unsigned long count =0,i=0,j=0;
	unsigned long nInf = 0;
	lKEY_dVALUE_PAIR* pKVP = (lKEY_dVALUE_PAIR*)malloc(sizeof(lKEY_dVALUE_PAIR) * (this->ncol* this->nrow)); //a vector (for each element in the matrix) of key value pairs of doubles 
	
	for (j = 0; j < this->ncol; j++){
		for (i = 0; i < this->nrow; i++){
			pKVP[count].lPosition = count; //Set each key value pair
			pKVP[count++].dValue = *(sMAT+j* this->nrow+i);
			if (*(sMAT+j* this->nrow+i) == std::numeric_limits<double>::infinity() ) {nInf++;}
		}
	}
	sort(pKVP,pKVP+(this->nrow*this->ncol),predicate_k_v_pair);//And sort based on the values

	for (i = 0; i < (this->nrow*this->ncol - nInf); i++){
		this->rMAT[pKVP[i].lPosition] = i+1 ;  //Assign the ranks based on the sorted key value pairs
	}
	for (i = (this->nrow*this->ncol - nInf); i < (this->nrow*this->ncol ); i++){
		this->rMAT[pKVP[i].lPosition] = std::numeric_limits<unsigned long>::max();  //Assign the ranks based on the sorted key value pairs
	}
	free(pKVP);
}

void DataClass::orderPerSubject(){
	//pmat1 is a pointer to double matrix, In the rows are the observers (microarrays), in the columns the subjects (patients)
	unsigned int j;	
	for (j = 0; j < this->ncol; j++){		
		sort(this->sMAT+j* this->nrow, this->sMAT+(j+1)* this->nrow); //sort column j
	} 
}


void DataClass::BN_from_R(){
	unsigned long i,j;
	this->T = 0;
	//Determine BN (= vector of length b indicating the number of known patients with 1..b measurments)  and T (=total number of pivots irrespective of b: t_b +b = T);
	for (j = 0; j < this->ncol; j++){
		i = this->nrow-1;
		while ((*(this->rMAT + j * this->nrow + i))==std::numeric_limits<unsigned long>::max()){ i--;} 
		this->BN[j] = i+1; //Set the number of real values for patient j
		(this->T)+=this->BN[j]; //And the number of pivots
	} 
	//cout << "T = " << this->T << endl;
}

void DataClass::preprocess( void ){
	//Order each column
	this->orderPerSubject();
		
	//Convert to int
	this->S2R();

	//Determine BN (= vector of length b indicating the number of known patients with 1..b measurments)  and T (=total number of pivots irrespective of b: t_b +b = T);
	this->BN_from_R();
	 
	this->R2Q();
}

DataClass::DataClass(double* pmat1, unsigned int n, unsigned int maxB){
	//Assume that all values in MAT1 are ordered within a subject with non finite values at the end (i.e. last columns))
	//Assume subject in the rows and observations in the columns
	//Transpose the matrix such that observation are in the rows and subjects in the columns
	transpose(pmat1,n,maxB);

	this->nrow = maxB;
	this->ncol = n;

	this->BN = (unsigned int*)malloc(sizeof(unsigned int) * this->ncol);

	this->sMAT = (double*)malloc(sizeof(double) * (this->ncol*this->nrow));
	this->qMAT = (unsigned long*)malloc(sizeof(unsigned long) * (this->ncol*this->nrow));
	this->rMAT = (unsigned long*)malloc(sizeof(unsigned long) * (this->ncol*this->nrow));	
	
	for (int i = 0; i < (this->ncol*this->nrow); i++){		this->sMAT[i] = pmat1[i];	}

	this->T =0; //Relation to article, : t_b = T-b
	this->seed =1;
	this->preprocess( );
}
	
DataClass::DataClass (const DataClass &obj) {
	this->nrow = obj.nrow;
	this->ncol = obj.ncol;

	this->T = obj.T; //Relation to article, : t_b = T-b
	this->seed = obj.seed;
	
	this->BN = NULL;
	this->qMAT = NULL;
	this->rMAT = NULL;

	if (obj.BN != NULL) {	
		this->BN = (unsigned int*)malloc(sizeof(unsigned int) * this->ncol);
		for (unsigned int i = 0; i < this->ncol; i++){
			this->BN[i] = obj.BN[i];
		}
	}

	if (obj.qMAT != NULL) {	
		this->qMAT = (unsigned long*)malloc(sizeof(unsigned long) * (this->ncol*this->nrow));
		for (int i = 0; i < (this->ncol*this->nrow); i++){
			this->qMAT [i] = obj.qMAT[i];
		}
	}

	if (obj.rMAT != NULL) {	
		this->rMAT = (unsigned long*)malloc(sizeof(unsigned long) * (this->ncol*this->nrow));
		for (int i = 0; i < (this->ncol*this->nrow); i++){
			this->rMAT [i] = obj.rMAT[i];
		}
	}
	if (obj.sMAT != NULL) {	
		this->sMAT = (double*)malloc(sizeof(double) * (this->ncol*this->nrow));
		for (int i = 0; i < (this->ncol*this->nrow); i++){
			this->sMAT [i] = obj.sMAT[i];
		}
	}
}


double getPsi(double* MAT1, unsigned int n, unsigned int maxB){
	DataClass dc(MAT1,n, maxB);
	return( dc.calculatePSI());

}

extern "C"{
SEXP getPsi202(SEXP MAT1){
	SEXP dim = getAttrib( MAT1, R_DimSymbol ) ;
	int nrow = INTEGER(dim)[0]; 
	int ncol = INTEGER(dim)[1]; 
 
	
	SEXP psi = allocVector(REALSXP,1);
	*REAL(psi) = getPsi(REAL(MAT1),nrow,ncol);
	
	return(psi);
}}
