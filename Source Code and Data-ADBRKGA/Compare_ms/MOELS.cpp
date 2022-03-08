//
// Created by qqq on 2021/10/31.
//

#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runMOELS(string XmlFile, string RscAlcFile, double& SchTime, int& iteration){
    clock_t start=clock();
    ReadFile(XmlFile,RscAlcFile);       //read model information-w
    ConfigParameter_MOELS();            //set the parameter values
    CalculateLevelList();               //calculate the levels of tasks

    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double> > cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    vector<double> rank_b(comConst.NumOfTsk, 0);
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(rank_b,cc,ww);
    population.resize(Parameter_MOELS.NumOfChromPerPop);
    //initial population
    for(int n = 0; n < Parameter_MOELS.NumOfChromPerPop; ++n){
        chromosome Chrom;
        IntChr(Chrom);
        Chrom.TskSchLst = GenerateAChromOrder(rank_b);   //此处是基于rank为概率进行选择生成任务调度顺序，每次生成调度顺序不一样
        GnrMS_Evl(Chrom);
        population[n] = Chrom;
    }
    sort(population.begin(), population.end(), SortPopOnFitValueByAscend);

    while(1){
        ++iteration;
        double RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        if ( RunTime >= SchTime ) {
            SchTime = RunTime;
            break;
        }
        vector<chromosome> NewPopulation(Parameter_MOELS.NumOfChromPerPop);
        //crossover
        #pragma omp parallel for
        for (int n = 0; n < Parameter_MOELS.NumOfChromPerPop; ++n){
            int parent1, parent2;
            SelectionTournament(parent1, parent2, Parameter_MOELS.NumOfChromPerPop);
            chromosome chrom1 = population[parent1];
            chromosome chrom2 = population[parent2];
            double a = double(rand()%100)/100;
            if(a <= Parameter_MOELS.CrossoverRate){
                chrom1 = Crossover_MOELS(chrom1, chrom2);
            }
            NewPopulation[n] = chrom1;
        }
        //Mutation
        #pragma omp parallel for
        for(int n = 0; n < Parameter_MOELS.NumOfChromPerPop; ++n){
            double a = double(rand()%100)/100;
            if(a <= Parameter_MOELS.MutationRate){
                Mutation_MOELS(NewPopulation[n]);
            }
        }
        #pragma omp parallel for
        for (int n = 0; n < Parameter_MOELS.NumOfChromPerPop; ++n) {
            GnrMS_Evl(NewPopulation[n]);
        }
        //{generate the next population}
        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end());
        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
        population.resize(Parameter_MOELS.NumOfChromPerPop);
    }
    ClearALL();
    return population[0].FitnessValue;
}


