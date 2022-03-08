#ifndef CSTCHANGE_CLASSDEFINE_H

#include "math.h"
#include <vector>
#include <string>
#include <set>
#include <map>
#include <list>

#define CSTCHANGE_CLASSDEFINE_H
using namespace std;
//{file}
class vfile {
public:
    string FileName;     //file name
    int source;          //the source of file, -1:from the shared server; i: from task i
    double size;         //the size of file
};

//{task}
class Task {
public:
    double length;          //the length of task
    vector<int> ElgRsc;     //the set of resources which are eligible to  perform this task -xy4
    vector<int> parents;    //the set of parent tasks
    vector<int> children;   //the set of child tasks
    vector<vfile> IFile;    //the set of input files
    vector<vfile> OFile;    //the set of output files
};

//{resources}
class Resource {
public:
    vector<int> ElgTsk;     //the set of tasks which it is eligible to perform
    double pc, bw;          //processing capacity, bandwidth

    Resource(int id, double pc, double bw) {
        this->pc = pc;
        this->bw = bw;
    }
};

class chromosome {
public:
    vector<int> RscAlcLst;       //resources allocation, task-to-resource mapping, (Match String)
    vector<int> TskSchLst;       //task scheduling order (Scheduling String)  -xy4
    vector<double> Code_RK;
    vector<double> EndTime;      //the finish time of task
    vector<double> RscAlcPart;
    vector<double> TskSchPart;
    vector<double> VTskSchPart;
    vector<double> VRscAlcPart;
    double FitnessValue;         //Fitness (makespan) -xy4
    //{sort chromosome and remove the same chromosome according to fitness value}
    bool operator<(const chromosome &otherChromosome)const {
        return this->FitnessValue + 1e-6 < otherChromosome.FitnessValue;
    }
//    //{sort chromosome according to fitness value and remove the same chromosome according to code }
//    bool operator<(const chromosome &otherChromosome)const {
//        int flag = -1;
//        for (int i = 0; i < this->TaskOrderList.size(); i++) {
//            if (this->VMAllocationList[i] != otherChromosome.VMAllocationList[i]) {
//                flag = i;
//                break;
//            }
//            if (this->TaskOrderList[i] != otherChromosome.TaskOrderList[i]) {
//                flag = i;
//                break;
//            }
//        }
//        if (flag == -1) {
//            return false;
//        } else {
//            if(fabs(this->FitnessValue-otherChromosome.FitnessValue)<1e-6) {
//                if(this->VMAllocationList[flag] != otherChromosome.VMAllocationList[flag]){
//                    return this->VMAllocationList[flag]<otherChromosome.VMAllocationList[flag];
//                }else {
//                    return this->TaskOrderList[flag]<otherChromosome.TaskOrderList[flag];
//                }
//            }else {
//                return this->FitnessValue < otherChromosome.FitnessValue;
//            }
//        }
//    }
};


class Paramet_CGA {
public:
    int NumOfChromPerPop;
    double CrossoverRate;      //crossover rate(CrossoverRate)
    double MutationRate;       //mutation rate (MutationRate)
};

class Paramet_HGA {
public:
    int NumOfChromPerPop;
    double EliteRate;         //crossover rate(CrossoverRate)
    double MutationRate;      //mutation rate (MutationRate)
};

class Paramet_LWSGA {
public:
    int NumOfChromPerPop;
    double CrossoverRate;     //crossover rate(CrossoverRate)
};

class Paramet_MOELS {
public:
    int NumOfChromPerPop;
    float CrossoverRate;
    float MutationRate;
};

class Paramet_HPSO {
public:
    int NumOfChromPerPop;
    double InertiaWeight;
    double c1;
    double c2;
};

class Paramet_ADBRKGA {
public:
    int NumOfChromPerPop;     //the number of chromosomes in each population
    double alpha;
    double beta;
    double BiasesRate;
    double ImmigrationRate;
    double ImprovementRate;
};

class ComConst {
public:
    int NumOfTsk;             //the number of Tasks
    int NumOfRsc;             //the number of resources
};

extern vector<Task> Tasks;
//extern vector<Task> OriginalTasks;
//extern vector<int> Oid;
//extern vector<int> Nid;
extern vector<vector<int> > TskLstInLvl; //task list (set) in each level
extern vector<int> LevelIdOfTask;        //the level of task
extern vector<vector<double> > ParChildTranFileSizeSum;
extern vector<Resource> Rscs;
extern vector<chromosome> population;
//extern vector<vector<chromosome> > populations;
extern vector<set<int> > Descendants;
extern vector<set<int> > Ancestors;
extern ComConst comConst;
extern double ModelScale;
extern Paramet_CGA Parameter_CGA;
extern Paramet_HGA Parameter_HGA;
extern Paramet_LWSGA Parameter_LWSGA;
extern Paramet_MOELS Parameter_MOELS;
extern Paramet_HPSO Parameter_HPSO;
extern Paramet_ADBRKGA Parameter_ADBRKGA;
extern int TrmFactor;

#endif //CSTCHANGE_CLASSDEFINE_H