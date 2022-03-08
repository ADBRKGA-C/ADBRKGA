//
// Created by xieyi on 2022/1/18.
//

#include "DHEFT.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"

double runDHEFT(string XmlFile, string RscAlcFile, double& SchTime) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    CalculateLevelList();

    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));

    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(Rank_b, cc, ww);

//    chromosome chrom_DHEFT = GnrChr_DHEFT(Rank_b);
//    chromosome chrom_DIHEFT3 = GnrChr_DIHEFT3(Rank_b);
    chromosome chrom_DIHEFT = GnrChr_DIHEFT(Rank_b);

    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    ClearALL();
    return chrom_DIHEFT.FitnessValue;
}

