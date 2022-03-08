
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"

double runHEFT(string XmlFile, string RscAlcFile, double& SchTime) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);

    CalculateLevelList();
    vector<double> ww(comConst.NumOfTsk, 0.0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0.0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);

    vector<double> Rank_b(comConst.NumOfTsk, 0.0);
    Calculate_Rank_b(Rank_b,cc,ww);
    chromosome Chrom_HEFT_b = GnrChr_HEFT_b(Rank_b);
//    vector<double> Rank_t(comConst.NumOfTsk, 0.0);
//    Calculate_Rank_t(Rank_t,cc,ww);
//    chromosome Chrom_HEFT_t = GnrChr_HEFT_t(Rank_t);

    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    ClearALL();
    return Chrom_HEFT_b.FitnessValue;
}