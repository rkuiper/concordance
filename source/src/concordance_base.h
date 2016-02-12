#ifndef CONCORDANCE_H
#define CONCORDANCE_H

#define EPSILON 1.0e-10
#define SEED    437401448   //Some random number that Matlab gave me


#include <stdint.h>
#include <map>
#include "sparsehash/dense_hash_map"

class CComp
{
    private:
    public:
		bool operator()(const int32_t & x,const int32_t & y) const;
};

class CHashFcn
{
    private:
    public:
		size_t operator()(const char * pKey) const;  
};

class CEqualKey
{
    private:
    public:
		bool operator()(const char * pKey1,const char * pKey2) const;
};




extern int         nStateElements;

extern char *      pPrevKey;

extern int *       pStateLimits;
extern int *       pStateSpan;
extern int *       pState;

extern int		   iMaxLevel;

extern double *    pInv1StateLimits;
extern double *    pInvStateLimits;


extern int32_t iKeyLen;

extern uint64_t MurmurHash(const void * key,int32_t len,uint32_t seed);
extern double LogSum(double dLogX,double dLogY);
extern double Concordance(int iPivot);
extern double Probability(int iElement);





extern google::dense_hash_map<const char *,std::map<double,double>,CHashFcn,CEqualKey> *    pPrevLevel;
extern google::dense_hash_map<const char *,std::map<double,double>,CHashFcn,CEqualKey> *    pLevel;


#endif  /*CONCORDANCE_H*/
