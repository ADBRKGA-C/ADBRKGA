
#ifndef CSTCHANGE_GENERATEACHROM_H
#include "common.h"
#define CSTCHANGE_GENERATEACHROM_H

void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);
void Calculate_Rank_t(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
vector<double> GnrDecimalsByAscend();
double DcdEvl(chromosome& chrom, bool isForward);
double NrmDcd(chromosome& chrom, bool isForward);
double HrsDcd_CTP(chromosome& chrom);
double HrsDcd_EFT(chromosome& ch);
double IFBS(chromosome& ch);
void LBCA_IFBS(chromosome& chrom);
chromosome GnrChr_HEFT_b_ADBRKGA(vector<double> Rank_b);
chromosome GnrChr_HEFT_t_ADBRKGA(vector<double> Rank_t);
chromosome GnrChr_IHEFT3_b(vector<double> Rank_b);
chromosome GnrChr_IHEFT3_t(vector<double> Rank_t);
chromosome GnrChr_Lvl_Ran();
chromosome GnrChr_DHEFT(vector<double>& Rank_b);
double ClcAvrReadyTime(int TskId, chromosome& chrom);
double IHEFT3(chromosome& ch);
void SeletRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskId, int& RscId, double& FinalStartTime, double& FinalEndTime);
double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime);

#endif //CSTCHANGE_GENERATEACHROM_H
