//
// Created by xieyi on 2022/1/18.
//

#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runHPSO(string XmlFile, string RscAlcFile, double& SchTime, int& iteration){
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);               //read model information
//    TaskCombine();
    ConfigParameter_HPSO();
    CalculateLevelList();
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);                           //calculate the average execution time of tasks
    C_Cal_Average(cc);                           //calculate the average transfer time among tasks
    vector<double> Rank_t(comConst.NumOfTsk,0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    Calculate_Rank_t(Rank_t, cc, ww);            //calcualte the rank_t
    Calculate_Rank_b(Rank_b, cc, ww);            //calcualte the rank_b
    vector<double> Rank_b1(comConst.NumOfTsk, 0);
    double MaxRank_b = 0;
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        if(MaxRank_b < Rank_b[i])
            MaxRank_b = Rank_b[i];
    }
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        Rank_b1[i] = MaxRank_b - Rank_b[i];
    }
//    chromosome Chrom_HEFT_b = GnrPrtByRank_EFT(Rank_b1);
//    population.push_back(Chrom_HEFT_b);
//    chromosome Chrom_HEFT_t = GnrPrtByRank_EFT(Rank_t);
//    population.push_back(Chrom_HEFT_t);
    for(int i = 0; i < Parameter_HPSO.NumOfChromPerPop; ++i){
        chromosome TemChrom1 = GnrPrtByRank_Rnd(Rank_t);
        chromosome TemChrom2 = GnrPrtByRank_Rnd(Rank_b1);
        population.push_back(TemChrom1);
        population.push_back(TemChrom2);
    }
    sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
    vector<chromosome>::iterator iter = population.begin();
    advance(iter,Parameter_HPSO.NumOfChromPerPop);
    population.assign(population.begin(),iter);
    vector<chromosome> Pbest = population;
    chromosome Gbest = population[0];
    while (1) {
        ++iteration;
        double RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        if (RunTime >= SchTime) {
            SchTime = RunTime;
            break;
        }
        #pragma omp parallel for
        for (int i = 0; i < Parameter_HPSO.NumOfChromPerPop; ++i) {
            UpdateParticle(population[i], Pbest[i], Gbest, RunTime, SchTime);
            DcdEvl(population[i], true);
        }
        for (int i = 0; i < Parameter_HPSO.NumOfChromPerPop; ++i) {
            if (population[i].FitnessValue + PrecisionValue < Pbest[i].FitnessValue) {
                Pbest[i] = population[i];
                if (Pbest[i].FitnessValue + PrecisionValue < Gbest.FitnessValue)
                {
                    Gbest = Pbest[i];
                }
            }
        }
    }
    ClearALL();
    return Gbest.FitnessValue;
}





