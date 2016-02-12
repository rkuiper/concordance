#ifndef __GETPSI_H__
#define __GETPSI_H__

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


double getPsi(double* MAT1, unsigned int n, unsigned int maxB);

extern "C" {
SEXP getPsi202(SEXP MAT1);
}


class DataClass
{
    protected:
		unsigned int seed;
		unsigned int nrow, ncol;
		double * sMAT;
		unsigned long* rMAT;
		unsigned long* qMAT;

		unsigned int* BN;
		unsigned long T;

		void R2Q( void );
		void S2R( void );
		void orderPerSubject( void );
		void BN_from_R( void );
		

 
	 public:
		~DataClass();
		 	
		double calculatePSI(void );

		void preprocess( void );

		DataClass(double* pmat1, unsigned int n, unsigned int maxB);
		DataClass(const DataClass &obj);


};
#endif /*__GETPSI_H__*/
