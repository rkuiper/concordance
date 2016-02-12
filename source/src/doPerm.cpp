#include <R.h>
#include <Rmath.h>
#include <random>


void shuffle(double *array, size_t n, unsigned int seed)
{
 	std::default_random_engine generator (seed);
	std::uniform_int_distribution<unsigned int> distribution(0,n-1);	
    if (n > 1) {    
		for (size_t i = 0; i < n - 1; i++) {
			size_t j = distribution(generator);
		 	double t = array[j];
			array[j] = array[i];
			array[i] = t;
		}
    }
}

extern "C"{

//  quickSort
//
//  This public-domain C implementation by Darel Rex Finley.
//
//  * Returns true if sort was successful, or false if the nested
//    pivots went too deep, in which case your array will have
//    been re-ordered, but probably not sorted correctly.
//
//  * This function assumes it is called with valid parameters.
//
//  * Example calls:
//    quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7



bool quickSort(double *arr, int elements) {

  #define  MAX_LEVELS  1000

  int  piv, beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R ;

  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=arr[L]; if (i==MAX_LEVELS-1) return false;
      while (L<R) {
        while (arr[R]>=piv && L<R) R--; if (L<R) arr[L++]=arr[R];
        while (arr[L]<=piv && L<R) L++; if (L<R) arr[R--]=arr[L]; }
      arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L; }
    else {
      i--; }}
  return true; }


double getOmega(int B, int T, int *nub, int* ub, int* nb){
	double teller =0;
	double noemer = 0;
	for (int idx =0; idx< *nub; idx++){
		double t = (double)T;
		double j = (double)ub[idx];
		double n = (double)nb[idx];
		noemer = noemer + n*sqrt( 45.0*(j-1)*j*t/((j+3)*(t-j-1)));
		if (ub[idx] == B){
			teller = teller + sqrt( 45*(j-1)*j*t/((j+3)*(t-j-1)));
		}
	}
	return(teller/noemer);	
}


void getWeight(long int B, long int T,long int N, double* weightB){
	//Version 1.02
	for (long int idx = 1; idx <= (B-1); idx++){
		//weightB[idx-1] = (double)(2*idx*(B-idx))/(double)(T*(T-1)*(B-1));
		//weightB[idx-1] = -1*(double)(2.0*idx*idx-2*idx*B) /(B*(B-1)*N*T+B*B*(1-B)*N);
		weightB[idx-1] = -1*(double)(2.0*idx*idx-2*idx*B) /((B-1)*T*(T-B));
	}
}

 void doPerm(int *ub,  int *nb, int *nub , int *nperm, int* seed ,double* psi){
	//This function is called from within R
	//ub = unique values of replicates (b)
	//nb = number of observations with b replicates
	//nub = number of unique b
	//seed = a seed (positive integer) for the random number generator (or -1 for the default)
	
	/*if (*seed>-1){
		srand((unsigned int)*seed)	;
	}*/
	
	//Determine the number of elements (T); Relation to article: t_b+b = T
	long int T = 0; 
	long int N = 0;
	double ** weightB=new double*[*nub];
	
	
	for (int idx = 0; idx < *nub; idx++){		
		T=T+nb[idx]*ub[idx];
		N=N+nb[idx];	
	}
	//Determine the weights to be used
	for (int idx = 0; idx < *nub; idx++){		
		weightB[idx]=new double[ub[idx]-1];
		getWeight(ub[idx],T,N,weightB[idx]);	
		//Rprintf("\nweights for b = %i:\n",ub[idx]);		
		/*for (int idx2 = 0; idx2 < (ub[idx]-1); idx2++){		
			Rprintf("\t%.10e,",weightB[idx][idx2]);		
		}*/
	}
	//error("HIER VERDER");

	//Define the matrix as a linear array of size T, and fill it with integers from (1..(T-1))
	//double * R=new double[(T-1)];#R1
	//for (int idx =0; idx < (T-1); idx++){ R[idx]=idx+1;	}#R2
	double * R=new double[T];//#R1
	for (int idx =0; idx < T; idx++){ R[idx]=idx+1;	}//#R2
	
	//start determine psi per object
	int cursorIdx = 0;	
	int perm=0;
	int idx1 =0;
	int idx2 =0;
	int idx3 =0;

	for (perm=0; perm<*nperm;perm++){
		R_CheckUserInterrupt();	
		//shuffle(R,T-1); #R4
		shuffle(R,T,(*seed)+perm); //#R4//Shuffle all elements within R
		cursorIdx = 0;	
		psi[perm] = 0;
		for (idx1 = 0; idx1 < *nub; idx1++){ //loop all unique b
			for (idx2 = 0; idx2 < nb[idx1]; idx2++){ //loop all observations within R with b number of replicates
				quickSort(&(R[cursorIdx]),(ub[idx1])); // sort all replicates within an observation
				for (idx3 = 0; idx3<(ub[idx1]-1);idx3++){ // loop all replicates within an observation
					//psi[perm] = psi[perm] + weightB[idx1][idx3] * (R[cursorIdx+1] - R[cursorIdx]-1);#R3
					psi[perm] = psi[perm] + weightB[idx1][idx3] * (R[cursorIdx+1] - R[cursorIdx]-1);//#R3
					cursorIdx++;
				}
				cursorIdx++;
			}
		}
		psi[perm] = 1-psi[perm];	
	}
	//end	

	for (int idx = 0; idx < *nub; idx++){
		delete [] weightB[idx];
	}

	delete [] weightB;
	delete [] R;
}
}

