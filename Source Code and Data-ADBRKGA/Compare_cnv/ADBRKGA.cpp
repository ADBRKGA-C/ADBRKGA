#include <common.h>
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runADBRKGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration, double& EndFitness){
    double CurNum = 0; //double RunTime = 0;
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_ADBRKGA();

    CalculateLevelList();
    vector<double> ww(comConst.NumOfTsk,0.0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0.0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    vector<double> Rank_b(comConst.NumOfTsk,0.0);
    Calculate_Rank_b(Rank_b, cc, ww);
    vector<double> Rank_t(comConst.NumOfTsk,0.0);
    Calculate_Rank_t(Rank_t, cc, ww);

    set<chromosome> TemSubPop;
    chromosome chrom_HEFT_b = GnrChr_HEFT_b_ADBRKGA(Rank_b);
    if (chrom_HEFT_b.FitnessValue < EndFitness + PrecisionValue_ct) {
        SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
        return chrom_HEFT_b.FitnessValue;
    }
    chromosome chrom_IHEFT3_b = GnrChr_IHEFT3_b(Rank_b);
    if (chrom_IHEFT3_b.FitnessValue < EndFitness + PrecisionValue_ct) {
        SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
        return chrom_IHEFT3_b.FitnessValue;
    }
    chromosome chrom_HEFT_t = GnrChr_HEFT_t_ADBRKGA(Rank_t);
    if (chrom_HEFT_t.FitnessValue < EndFitness + PrecisionValue_ct) {
        SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
        return chrom_HEFT_t.FitnessValue;
    }
    chromosome chrom_IHEFT3_t = GnrChr_IHEFT3_t(Rank_t);
    if (chrom_IHEFT3_t.FitnessValue < EndFitness + PrecisionValue_ct) {
        SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
        return chrom_IHEFT3_t.FitnessValue;
    }
    chromosome chrom_DHEFT = GnrChr_DHEFT(Rank_b);
    if (chrom_DHEFT.FitnessValue < EndFitness + PrecisionValue_ct) {
        SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
        return chrom_DHEFT.FitnessValue;
    }
    TemSubPop.insert(chrom_HEFT_b);
    TemSubPop.insert(chrom_IHEFT3_b);
    TemSubPop.insert(chrom_HEFT_t);
    TemSubPop.insert(chrom_IHEFT3_t);
    TemSubPop.insert(chrom_DHEFT);

    int terminationNum = ceil(TrmFactor * ModelScale/sqrt(comConst.NumOfTsk)/Parameter_ADBRKGA.NumOfChromPerPop); double TotalNum = terminationNum;
    int NumOfNoImpGen = 0;
    while (TemSubPop.size() < Parameter_ADBRKGA.NumOfChromPerPop){
        chromosome chrom_Ran = GnrChr_Lvl_Ran();
        AdpDcd(chrom_Ran,CurNum,TotalNum);
        if (chrom_Ran.FitnessValue < EndFitness + PrecisionValue_ct){
            SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
            return  chrom_Ran.FitnessValue;
        }
        TemSubPop.insert(chrom_Ran);
    }
    population.assign(TemSubPop.begin(),TemSubPop.end());

    double bestFitness = population[0].FitnessValue;

    vector<double> A(Parameter_ADBRKGA.NumOfChromPerPop);
    CalSlctProb_Rank(1+1.0/Parameter_ADBRKGA.NumOfChromPerPop, A ,Parameter_ADBRKGA.NumOfChromPerPop); //calculate the cumulative probabilities
    int ImprovementNumber =ceil(Parameter_ADBRKGA.ImprovementRate * Parameter_ADBRKGA.NumOfChromPerPop);
    int ImmigrationNumber = ceil(Parameter_ADBRKGA.ImmigrationRate * Parameter_ADBRKGA.NumOfChromPerPop);

    while (1) {
        ++iteration;
//        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        vector<chromosome> NewPopulation(Parameter_ADBRKGA.NumOfChromPerPop);

        for (int n = 0; n < Parameter_ADBRKGA.NumOfChromPerPop; ++n) {
            int ind1 = SelectChrom(A);
            int ind2 = SelectChrom(A);
            while (ind1 == ind2) {
                ind2 = SelectChrom(A);
            }
            chromosome chrom1 ;
            chromosome chrom2 ;
            if (population[ind1].FitnessValue + PrecisionValue < population[ind2].FitnessValue){
                chrom1 = population[ind1];
                chrom2 = population[ind2];
            } else {
                chrom1 = population[ind2];
                chrom2 = population[ind1];
            }
            NewPopulation[n] = Crs_BPUC(chrom1,chrom2);
            AdpDcd(NewPopulation[n],CurNum,TotalNum);
            if (NewPopulation[n].FitnessValue < EndFitness + PrecisionValue_ct) {
                SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
                return NewPopulation[n].FitnessValue;
            }
        }

        sort(NewPopulation.begin(), NewPopulation.end(), SortPopOnFitValueByAscend);

        for (int n = 0; n < ImprovementNumber; ++n) {
            LBCA_IFBS(NewPopulation[n]);
            if (NewPopulation[n].FitnessValue < EndFitness + PrecisionValue_ct) {
                SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
                return NewPopulation[n].FitnessValue;
            }
        }

        set<chromosome> NxtPop;
        NxtPop.insert(population.begin(),population.end());
        NxtPop.insert(NewPopulation.begin(),NewPopulation.end());
        set<chromosome>::iterator iter = NxtPop.begin();
        advance(iter,(Parameter_ADBRKGA.NumOfChromPerPop - ImmigrationNumber));
        NxtPop.erase(iter,NxtPop.end());
        while (NxtPop.size() < Parameter_ADBRKGA.NumOfChromPerPop) {
            chromosome chrom_Ran = GnrChr_Lvl_Ran();
            AdpDcd(chrom_Ran,CurNum,TotalNum);
            if (chrom_Ran.FitnessValue < EndFitness + PrecisionValue_ct) {
                SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
                return chrom_Ran.FitnessValue;
            }
            NxtPop.insert(chrom_Ran);
        }
        population.assign(NxtPop.begin(),NxtPop.end());
        ++NumOfNoImpGen; CurNum = XY_MAX(CurNum,NumOfNoImpGen);
        if( population[0].FitnessValue + PrecisionValue < bestFitness ){
            bestFitness = population[0].FitnessValue;
            NumOfNoImpGen = 0;
        }
        if( (NumOfNoImpGen == terminationNum) || (bestFitness < EndFitness + PrecisionValue_ct) ){
            SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
            break;
        }
    }
//    ClearALL();
    return population[0].FitnessValue;
}