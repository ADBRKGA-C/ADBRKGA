
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runLWSGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration, double& EndFitness) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);     //read model information-w
    ConfigParameter_LWSGA();           //set the parameter values
    CalculateLevelList();              //calculate the levels of tasks
    population.resize(Parameter_LWSGA.NumOfChromPerPop);
    for ( int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; ++n ) {
        chromosome chrom;
        IntChr(chrom);
        chrom.TskSchLst = GnrSS_Lvl();
        for ( int i = 0; i < comConst.NumOfTsk; ++i ) {
            chrom.RscAlcLst[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()];
        }
        population[n] = chrom;
        DcdEvl(population[n], true);                                             //decoding
        if (population[n].FitnessValue < EndFitness + PrecisionValue_ct) {
            SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
            return population[n].FitnessValue;
        }
    }

    sort(population.begin(), population.end(), SortPopOnFitValueByAscend);       //sorting
    double bestFitness = population[0].FitnessValue;
    int terminationNum = ceil(TrmFactor * ModelScale/sqrt(comConst.NumOfTsk)/Parameter_LWSGA.NumOfChromPerPop);
    int NumOfNoImpGen = 0;
    while (1) {
        ++iteration;
        vector<chromosome> NewPopulation(Parameter_LWSGA.NumOfChromPerPop) ;
        for ( int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; n += 2 ) {
            int parent1 = -1;
            int parent2 = -1;
            SelectionTournament(parent1, parent2,Parameter_LWSGA.NumOfChromPerPop);  //select two chromosomes using the tournament method
            double rand = RandomDouble(0, 1);
            chromosome TemChromosome1 = population[parent1];
            chromosome TemChromosome2 = population[parent2];
            if ( rand < Parameter_LWSGA.CrossoverRate ) {
                Crossover_LWSGA(TemChromosome1, TemChromosome2);                     //crossover
            } else {
                Mutation_LWSGA(TemChromosome1);                                      //mutation
                Mutation_LWSGA(TemChromosome2);
            }
            DcdEvl(TemChromosome1, true);
            if (TemChromosome1.FitnessValue < EndFitness + PrecisionValue_ct) {
                SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
                return TemChromosome1.FitnessValue;
            }
            DcdEvl(TemChromosome2, true);
            if (TemChromosome2.FitnessValue < EndFitness + PrecisionValue_ct) {
                SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
                return TemChromosome2.FitnessValue;
            }
            NewPopulation[n] = TemChromosome1;
            NewPopulation[n+1] = TemChromosome2;
        }

        //{generate the next population}
        population.insert(population.end(),NewPopulation.begin(),NewPopulation.end());
        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
        population.resize(Parameter_LWSGA.NumOfChromPerPop);
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