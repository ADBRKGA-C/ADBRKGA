

#ifndef CSTCHANGE_GENERATEACHROM_H
#define CSTCHANGE_GENERATEACHROM_H

#include "common.h"
void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_t(vector<double>& rankList, vector<double>& w, vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
vector<int> GnrSS_TS();
vector<int> GnrSS_Lvl();
double GnrMS_Evl(chromosome& chrom);
double HrsDcd_EFT(chromosome& ch);
chromosome GnrChr_HEFT_ADBRKGA(vector<double> Rank_b);
chromosome GnrChr_HEFT(vector<double> rnk);
chromosome GnrChrDHEFT(vector<double>& Rank_b);
void GnrChr_Ran(chromosome& chrom);
chromosome GnrChr_HEFT_t(vector<double> rank_t);
chromosome GnrChr_HEFT_b_t(vector<double> rank_b_t);
chromosome GnrChr_HEFT_Baseline();

#endif //CSTCHANGE_GENERATEACHROM_H
