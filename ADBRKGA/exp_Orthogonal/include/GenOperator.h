

#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H
double NrmDcd(chromosome& chrom, bool isForward);
void Crs_BPUC(chromosome& chrom1,chromosome& chrom2,chromosome& temCh);
double IFBD(chromosome& ch);
void LBCA(chromosome& chrom);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SelectChrom(vector<double>& A);
void HrsDcd_CTP(chromosome& chrom);


#endif //CSTCHANGE_CROSSOVER_H
