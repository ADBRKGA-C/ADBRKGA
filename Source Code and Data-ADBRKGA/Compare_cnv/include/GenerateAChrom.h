
#ifndef CSTCHANGE_GENERATEACHROM_H
#include "common.h"
#define CSTCHANGE_GENERATEACHROM_H

void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);
void Calculate_Rank_t(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
vector<int> GnrSS_TS();
vector<int> GnrSS_Lvl();
vector<double> GnrDecimalsByAscend();
double DcdEvl(chromosome& chrom, bool isForward);
double NrmDcd(chromosome& chrom, bool isForward);
double HrsDcd_CTP(chromosome& chrom);
double HrsDcd_EFT(chromosome& ch);
double GnrMS_Evl(chromosome& chrom);
double IFBS(chromosome& ch);
void LBCA_IFBS(chromosome& chrom);
chromosome GnrChr_HEFT_b_ADBRKGA(vector<double> Rank_b);
chromosome GnrChr_HEFT_t_ADBRKGA(vector<double> Rank_t);
chromosome GnrChr_HEFT_b(vector<double> Rank_b);
chromosome GnrChr_HEFT_t(vector<double> Rank_t);
chromosome GnrChr_IHEFT3_b(vector<double> Rank_b);
chromosome GnrChr_IHEFT3_t(vector<double> Rank_t);
chromosome GnrChr_IHEFT3(vector<double> Rank_b, vector<double> Rank_t);
chromosome GnrChr_DHEFT_old(vector<double>& Rank_b);
chromosome GnrChr_Ran();
chromosome GnrChr_Lvl_Ran();
chromosome GnrChr_HEFT_Baseline();
chromosome GnrChr_DIHEFT3(vector<double>& Rank_b);
chromosome GnrChr_DHEFT(vector<double>& Rank_b);
chromosome GnrPrtByRank_Rnd(vector<double>& Rank);
chromosome GnrPrtByRank_EFT(vector<double>& Rank);
double ClcAvrReadyTime(int TskId, chromosome& chrom);
double ClcAvrReadyTime2(int TskId, chromosome& chrom);  //计算传输时间取最大
double IHEFT3(chromosome& ch);
void SeletRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskId, int& RscId, double& FinalStartTime, double& FinalEndTime);
double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime);
vector<int> GenerateAChromOrder(vector<double>& rank);
void RepairMapAndGnrRscAlcLst(chromosome& ch);
int FindNearestRscId(int TaskId, double value);
void RepairPriorityAndGnrSchOrd(chromosome& chrom);
void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime);

#endif //CSTCHANGE_GENERATEACHROM_H
