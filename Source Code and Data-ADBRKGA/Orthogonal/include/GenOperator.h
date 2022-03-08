
#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H

chromosome Crs_BPUC(chromosome& chrom1,chromosome& chrom2);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SelectChrom(vector<double>& A);
void AdpDcd (chromosome&chrom ,double& CurTime,double& TotalTime);

#endif //CSTCHANGE_CROSSOVER_H
