

#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H
void CheckChrom(chromosome& chrom);
double DcdEvl(chromosome& chrom, bool isForward);
double NrmDcd(chromosome& chrom, bool isForward);
void Crs_BPUC(chromosome& chrom1,chromosome& chrom2,chromosome& temCh);
double IFBD(chromosome& ch);
void LBCA(chromosome& chrom);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SelectChrom(vector<double>& A);
void SelectionTournament(int& parent_1, int& parent_2 ,int& NumOfChormPerPop);
void MtnSS_TS(chromosome& a);
void MtnMS_SP(chromosome& a);
void Crossover_CGA(chromosome& pop1, chromosome& pop2);
void Mutation_CGA(chromosome& a);
void GnrTskSchLst_HGA(chromosome& chrom);
void HrsDcd_CTP(chromosome& chrom);
void Crossover_HGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_HGA(chromosome& chrom);
void RscLoadAdjust_HGA(vector<chromosome>& chromosomes);
void Crossover_LWSGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_LWSGA(chromosome& chrom);
chromosome Crossover_NGA(chromosome& pop1,chromosome& pop2, bool& flag);
void Mutation_NGA(chromosome& chrom);

#endif //CSTCHANGE_CROSSOVER_H
