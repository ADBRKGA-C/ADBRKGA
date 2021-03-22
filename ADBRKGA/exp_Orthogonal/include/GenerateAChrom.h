

#ifndef CSTCHANGE_GENERATEACHROM_H
#define CSTCHANGE_GENERATEACHROM_H

#include "common.h"
void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
double HrsDcd_EFT(chromosome& ch);
chromosome GnrChrDHEFT(vector<double>& Rank_b);
chromosome GnrChr_HEFT_ADBRKGA(vector<double> Rank_b);
void GnrChr_Ran(chromosome& chrom);


#endif //CSTCHANGE_GENERATEACHROM_H
