#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

using namespace std;

double runCGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);  //read model information
    ConfigParameter_CGA();          //set the parameter values
    CalculateLevelList();           //calculate the levels of tasks
    chromosome chrom;
    IntChr(chrom);
    chrom.TskSchLst = GnrSS_TS();
    population.resize(Parameter_CGA.NumOfChromPerPop);
//    #pragma omp parallel for
    for ( int n = 0; n < Parameter_CGA.NumOfChromPerPop - 1; ++n ) {
        int k = rand() % comConst.NumOfTsk + 1;
        chromosome TemChrom = chrom;
        while ( k-- ){
            MtnSS_TS(TemChrom);
        }
        for ( int i = 0; i < comConst.NumOfTsk; ++i ) {
            int size = Tasks[i].ElgRsc.size();
            TemChrom.RscAlcLst[i] = Tasks[i].ElgRsc[rand() % size];
        }
        population[n] = TemChrom;
    }
    chromosome ch_b = GnrChr_HEFT_Baseline();
    population[Parameter_CGA.NumOfChromPerPop - 1] = ch_b;
//    #pragma omp parallel for
    for ( int n = 0; n < Parameter_CGA.NumOfChromPerPop; ++n ) {
        DcdEvl(population[n],true);                          //decoding
    }
    sort(population.begin(), population.end(), SortPopOnFitValueByAscend);  //sorting
    chromosome BestChromosome = population[0];
    vector<double> A(Parameter_CGA.NumOfChromPerPop);
    CalSlctProb_Rank(1+1.0/Parameter_CGA.NumOfChromPerPop, A ,Parameter_CGA.NumOfChromPerPop);     //calculate the cumulative probabilities


    while (1) {
        ++iteration;
        //{terminate the algorithm according to running time}
        double RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        if ( RunTime >= SchTime ) {
            SchTime = RunTime;
            break;
        }
        vector<chromosome> NewPopulation(Parameter_CGA.NumOfChromPerPop);
        //{selection and crossover}
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_CGA.NumOfChromPerPop; n += 2 ) {
            int parent1 = SelectChrom(A);
            int parent2 = parent1;
            while (parent1 == parent2)
                parent2 = SelectChrom(A);
            chromosome chrom1 = population[parent1];
            chromosome chrom2 = population[parent2];
            Crossover_CGA(chrom1, chrom2);
            NewPopulation[n] = chrom1;
            NewPopulation[n+1] = chrom2;
        }
        //{mutation}
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_CGA.NumOfChromPerPop; ++n ) {
            Mutation_CGA(NewPopulation[n]);
        }
        //{decoding}
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_CGA.NumOfChromPerPop; ++n ) {
            DcdEvl(NewPopulation[n],true);
        }
        //{sorting}
        sort(NewPopulation.begin(), NewPopulation.end(), SortPopOnFitValueByAscend);
        //{Elite preservation}
        if ( NewPopulation[0].FitnessValue + PrecisionValue <  BestChromosome.FitnessValue ) {
            BestChromosome = NewPopulation[0];
        } else {
            NewPopulation[Parameter_CGA.NumOfChromPerPop - 1] = BestChromosome;
        }
        population = NewPopulation;
    }
    ClearALL();
    return BestChromosome.FitnessValue;
}
