//
// Created by qqq on 2021/10/31.
//

#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runMOELS(string XmlFile, string RscAlcFile, double& SchTime, int& iteration, double& EndFitness){
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
        Chrom.TskSchLst = GenerateAChromOrder(rank_b);   //generate randomly task schedule order based on rank
        GnrMS_Evl(Chrom);
        if (Chrom.FitnessValue < EndFitness + PrecisionValue_ct) {
            SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
            return Chrom.FitnessValue;
        }
        population[n] = Chrom;
    }
    sort(population.begin(), population.end(), SortPopOnFitValueByAscend);

    double bestFitness = population[0].FitnessValue;
    int terminationNum = ceil(TrmFactor * ModelScale/sqrt(comConst.NumOfTsk)/Parameter_MOELS.NumOfChromPerPop);
    int NumOfNoImpGen = 0;
    while(1){
        ++iteration;
        vector<chromosome> NewPopulation(Parameter_MOELS.NumOfChromPerPop);
        //crossover
        for (int n = 0; n < Parameter_MOELS.NumOfChromPerPop; ++n){
            int parent1, parent2;
            SelectionTournament(parent1, parent2, Parameter_MOELS.NumOfChromPerPop);     //selection
            chromosome chrom1 = population[parent1];
            chromosome chrom2 = population[parent2];
            double a = double(rand()%100)/100;
            if(a <= Parameter_MOELS.CrossoverRate){
                chrom1 = Crossover_MOELS(chrom1, chrom2);
            }
            a = double(rand()%100)/100;
            if(a <= Parameter_MOELS.MutationRate){
                Mutation_MOELS(chrom1);
            }
            GnrMS_Evl(chrom1);
            if (chrom1.FitnessValue < EndFitness + PrecisionValue_ct) {
                SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
                return chrom1.FitnessValue;
            }
            NewPopulation[n] = chrom1;
        }

        //{generate the next population}
        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end());
        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
        population.resize(Parameter_MOELS.NumOfChromPerPop);
        ++NumOfNoImpGen;
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