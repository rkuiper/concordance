
#include "concordance_base.h"
#include <R.h>
#include <Rmath.h>

using namespace std;
using google::dense_hash_map;
int         nStateElements;

char *      pPrevKey;

int *       pStateLimits;
int *       pStateSpan;
int *       pState;

int		   iMaxLevel;

double *    pInv1StateLimits;
double *    pInvStateLimits;


int32_t iKeyLen;


//----------------------------------------
bool CComp::operator()(const int32_t & x,const int32_t & y) const
{
    return x > y;
}

//----------------------------------------------------------------

size_t CHashFcn::operator()(const char * pKey) const
{
   return MurmurHash(pKey,iKeyLen,SEED);
}    

//----------------------------------------------------------------
bool CEqualKey::operator()(const char * pKey1,const char * pKey2) const
{
	return (pKey1 == pKey2) || (pKey1 && pKey2 && memcmp(pKey1,pKey2,iKeyLen) == 0);
}




uint64_t MurmurHash(const void * key,int32_t len,uint32_t seed) //MurmurHash2 64 bits version by by Austin Appleby
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int32_t r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len>>3);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h ^= k;
		h *= m; 
	}

	const uint8_t * data2 = (const uint8_t*)data;

	switch(len & 7)
	{
	case 7: h ^= uint64_t(data2[6]) << 48;
	case 6: h ^= uint64_t(data2[5]) << 40;
	case 5: h ^= uint64_t(data2[4]) << 32;
	case 4: h ^= uint64_t(data2[3]) << 24;
	case 3: h ^= uint64_t(data2[2]) << 16;
	case 2: h ^= uint64_t(data2[1]) << 8;
	case 1: h ^= uint64_t(data2[0]);
	        h *= m;
	};
 
	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
} 


//----------------------------------------------------------------
double LogSum(double dLogX,double dLogY)
{
    if(dLogX>=dLogY)
    {
        return dLogX+log(1.0+exp(dLogY-dLogX));
    }
    
    return dLogY+log(1.0+exp(dLogX-dLogY));
}
//----------------------------------------------------------------
double Concordance(int iPivot)
{
    int     iElement;
    
    double  dStateElement;
    double  dSum;
    
    dSum=0.0;
    for(iElement=0;iElement<nStateElements;iElement++)
    {
		if (iElement!=iPivot){
       		dStateElement=double(pState[iElement]);
       		dSum+=2.0*dStateElement*pInv1StateLimits[iElement]*(1.0-dStateElement*pInvStateLimits[iElement])*pStateLimits[iElement]/(iMaxLevel-pStateLimits[iElement]);
		}	
    }
    return dSum;
}
//----------------------------------------------------------------
double Probability(int iElement)
{
    int iSum;
    
    iSum=pStateLimits[iElement]-pState[iElement];
    
    for(iElement++;pStateLimits[iElement-1]==pStateLimits[iElement] && pState[iElement-1]==pState[iElement];iElement++)
    {
        iSum+=pStateLimits[iElement]-pState[iElement];
    }
    
    return log(double(iSum));
}

//----------------------------------------------------------------


google::dense_hash_map<const char *,std::map<double,double>,CHashFcn,CEqualKey> *    pPrevLevel;
google::dense_hash_map<const char *,std::map<double,double>,CHashFcn,CEqualKey> *    pLevel;


