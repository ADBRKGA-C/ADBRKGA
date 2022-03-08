
#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H

chromosome Crs_BPUC(chromosome& chrom1,chromosome& chrom2);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SelectChrom(vector<double>& A);
void SelectionTournament(int& parent_1, int& parent_2 ,int& NumOfChormPerPop);
void MtnSS_TS(chromosome& a);
void MtnMS_SP(chromosome& a);
void Crossover_CGA(chromosome& pop1, chromosome& pop2);
void Mutation_CGA(chromosome& a);
void GnrTskSchLst_HGA(chromosome& chrom);
void Crossover_HGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_HGA(chromosome& chrom);
void RscLoadAdjust_HGA(vector<chromosome>& chromosomes);
void Crossover_LWSGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_LWSGA(chromosome& chrom);
//void selection_Tournament(vector<chromosome>& population, int& parent_1, int& parent_2);
chromosome Crossover_MOELS(chromosome &chrom1, chromosome &chrom2);
void Mutation_MOELS(chromosome &chrom1);
void AdpDcd (chromosome&chrom ,double& CurTime,double& TotalTime);

#endif //CSTCHANGE_CROSSOVER_H
