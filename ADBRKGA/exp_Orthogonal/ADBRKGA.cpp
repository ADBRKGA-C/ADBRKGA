


#include <common.h>
#include "ADBRKGA.h"
#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runADBRKGA(string XmlFile, string RscAlcFile, Orthogonal orthogonal,double& SchTime, int& iteration){
    double RunTime = 0;
    double Time = 0;
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_ADRKGA(orthogonal);
    CalculateLevelList();

    vector<double> Rank_b(comConst.NumOfTsk,0);
    vector<double> Rank_t(comConst.NumOfTsk,0);
    vector<double> Rank_b_t(comConst.NumOfTsk,0);
    vector<double> ww(comConst.NumOfTsk,0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));

    W_Cal_Average(ww);
    C_Cal_Average(cc);//-w
    Calculate_Rank_b(Rank_b, cc, ww);
//    multiset<chromosome> TemSubPop;
    set<chromosome> TemSubPop;
    chromosome chrom_DHEFT = GnrChrDHEFT(Rank_b);
    chromosome chrom_HEFT = GnrChr_HEFT_ADBRKGA(Rank_b);
    TemSubPop.insert(chrom_DHEFT);
    TemSubPop.insert(chrom_HEFT);
    while (TemSubPop.size() < Parameter_ADRKGA.NumOfChormPerPop){
        chromosome chrom;
        IntChr(chrom);
        GnrChr_Ran(chrom);
        Time = (double) (clock() - start) / CLOCKS_PER_SEC;
        AdpDcd(chrom,Time,SchTime);
        TemSubPop.insert(chrom);
    }
    population.assign(TemSubPop.begin(),TemSubPop.end());
    double BestFitness = population[0].FitnessValue;

    vector<double> A(Parameter_ADRKGA.NumOfChormPerPop);
    CalSlctProb_Rank(1+1.0/Parameter_ADRKGA.NumOfChormPerPop, A ,Parameter_ADRKGA.NumOfChormPerPop); //calculate the cumulative probabilities
    while (1) {
        ++iteration;
        //{terminate the algorithm according to running time}
        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        if (RunTime >= SchTime) {
            SchTime = RunTime;
            break;
        }

        if (RunTime < SchTime) {
            vector<chromosome> NewPopulation(Parameter_ADRKGA.NumOfChormPerPop);
            for (int n = 0; n < Parameter_ADRKGA.NumOfChormPerPop; ++n) {
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
                chromosome temCh;
                IntChr(temCh);
                Crs_BPUC(chrom1,chrom2,temCh);
                NewPopulation[n] = temCh;
            }

            #pragma omp parallel for
            for ( int n = 0; n < Parameter_ADRKGA.NumOfChormPerPop; ++n ) {
                Time = (double) (clock() - start) / CLOCKS_PER_SEC;
                AdpDcd(NewPopulation[n],Time,SchTime);
            }
            sort(NewPopulation.begin(), NewPopulation.end(),SortPopOnFitValueByAscend);
            int ImprovementNum =ceil(Parameter_ADRKGA.ImprovementRate*Parameter_ADRKGA.NumOfChormPerPop);
            #pragma omp parallel for
            for (int n = 0; n < ImprovementNum; ++n) {
                IFBD(NewPopulation[n]);                                //SS improvement
                LBCA(NewPopulation[n]);
            }

            int ImmigrationNumber = ceil(Parameter_ADRKGA.ImmigrationRate * Parameter_ADRKGA.NumOfChormPerPop);
//            multiset<chromosome> NxtSubPop;
            set<chromosome> NxtPop;
            NxtPop.insert(population.begin(),population.end());
            NxtPop.insert(NewPopulation.begin(),NewPopulation.end());
//            multiset<chromosome>::iterator iter = NxtSubPop.begin();
            set<chromosome>::iterator iter = NxtPop.begin();
            advance(iter,(Parameter_ADRKGA.NumOfChormPerPop - ImmigrationNumber));
            NxtPop.erase(iter,NxtPop.end());
            while (NxtPop.size() < Parameter_ADRKGA.NumOfChormPerPop) {
                chromosome temCh;
                IntChr(temCh);
                GnrChr_Ran(temCh);
                Time = (double) (clock() - start) / CLOCKS_PER_SEC;
                AdpDcd(temCh,Time,SchTime);
                NxtPop.insert(temCh);
            }
            population.assign(NxtPop.begin(),NxtPop.end());
            //select Top N chromosomes to form next population
            if ( population[0].FitnessValue + PrecisionValue < BestFitness ) {
                BestFitness = population[0].FitnessValue;                   //update the best fitness
            }
        }
    }
    ClearALL();
    return BestFitness;
}