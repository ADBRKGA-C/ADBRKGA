#include <cstdlib>
#include <GenerateAChrom.h>
#include <unordered_set>
#include "tools.hpp"
#include "common.h"

using namespace std;

double NrmDcd(chromosome& ch, bool IsFrw) {
    double makespan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    vector<int > upr(comConst.NumOfTsk,0.0);
    vector<int> RTI;
    if(IsFrw){
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].parents.size();
            if (upr[i]==0){
                RTI.push_back(i);
            }
        }
    } else{
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].children.size();
            if (upr[i]==0){
                RTI.push_back(i);
            }
        }
    }

    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double decimal ;
        double tmp = 1;
        int taskIndex = -1;
        for(int j = 0; j < RTI.size(); ++j) {
            decimal = ch.Code_RK[RTI[j]] - floor(ch.Code_RK[RTI[j]]);
            if(decimal < tmp) {
                tmp = decimal;
                taskIndex = RTI[j];
            }
        }
        int RscIndex = floor(ch.Code_RK[taskIndex]);  //obtain the resource (Rsc) allocated to the task
        double ReadyTime = 0;
        double StartTime = 0;
        RTI.erase(find(RTI.begin(), RTI.end(), taskIndex));
        if(IsFrw) {                              //forward-loading
            for (int l = 0; l < Tasks[taskIndex].children.size(); ++l) {
                int childId = Tasks[taskIndex].children[l];
                upr[childId] = upr[childId] - 1;
                if (upr[childId]==0){
                    RTI.push_back(childId);
                }
            }
            if (Tasks[taskIndex].parents.size() != 0) {
                for (int j = 0; j < Tasks[taskIndex].parents.size(); ++j) {
                    int ParentTask = Tasks[taskIndex].parents[j];
                    int ParentRsc = floor(ch.Code_RK[ParentTask]);
                    double TransferTime = 0;
                    if(RscIndex != ParentRsc) {
                        TransferTime = ParChildTranFileSizeSum[ParentTask][taskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw)); // -xy
                    }
                    double sum = ch.EndTime[ParentTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        } else {                                //backward-loading
            for (int l = 0; l < Tasks[taskIndex].parents.size(); ++l) {
                int parentId = Tasks[taskIndex].parents[l];
                upr[parentId] = upr[parentId] - 1;
                if (upr[parentId]==0){
                    RTI.push_back(parentId);
                }
            }
            if (Tasks[taskIndex].children.size() != 0) {
                for (int j = 0; j < Tasks[taskIndex].children.size(); ++j) {
                    int ChildTask = Tasks[taskIndex].children[j];
                    int ChildRsc = floor(ch.Code_RK[ChildTask]);
                    double TransferTime = 0;
                    if(RscIndex != ChildRsc) {
                        TransferTime = ParChildTranFileSizeSum[taskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                    }
                    double sum = ch.EndTime[ChildTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        }
        set<double>::iterator pre  = ITL[RscIndex].begin();
        set<double>::iterator post = ITL[RscIndex].begin();
        ++post;
        double ExecutionTime = Tasks[taskIndex].length / Rscs[RscIndex].pc;
        //{find an idle time-slot in ITL which can finish the task  at the earliest}
        while(post != ITL[RscIndex].end()) {
            if((*post - *pre) >= ExecutionTime && ReadyTime <= (*post)-ExecutionTime) {
                StartTime = XY_MAX(*pre, ReadyTime);
                break;
            } else {
                ++pre;
                ++pre;
                ++post;
                ++post;
            }
        }
        ch.EndTime[taskIndex] = StartTime + ExecutionTime;
        if (makespan < ch.EndTime[taskIndex]) {
            makespan = ch.EndTime[taskIndex];
        }
        //{update ITL}
        if(ITL[RscIndex].find(StartTime) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(StartTime);
        } else {
            ITL[RscIndex].insert(StartTime);
        }

        if(ITL[RscIndex].find(ch.EndTime[taskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[taskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[taskIndex]);
        }

    }
    ch.FitnessValue = makespan;
    return ch.FitnessValue;
}

void HrsDcd_CTP(chromosome& ch) {
    vector<double> w(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<int> ind(comConst.NumOfTsk);
    vector<vector<double>> TransferTime(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    //{calculate the transfer time between tasks when resource(Rsc) allocation has been determined}
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = floor(ch.Code_RK[i]);
        if(Tasks[i].parents.size() !=  0){
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {
                int parent = Tasks[i].parents[j];
                int ParRsc = floor(ch.Code_RK[parent]);
                if(ParRsc != RscIndex){
                    TransferTime[parent][i] = ParChildTranFileSizeSum[parent][i] / VALUE * 8 / XY_MIN(Rscs[RscIndex].bw,Rscs[ParRsc].bw) ;
                }
            }
        }
        w[i] = Tasks[i].length / Rscs[RscIndex].pc;
    }
    Calculate_Rank_b(Rank_b,TransferTime, w);
    IndexSort(ind, Rank_b);
    set<double> st;
    vector<double> dect;
    double a = 0;
    while(st.size() < comConst.NumOfTsk) {
        a = double(rand() % 10000)/10000;
        st.insert(a);
    }
    dect.assign(st.begin(),st.end());
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.Code_RK[TaskIndex] = floor(ch.Code_RK[TaskIndex]) + dect[i];
    }
}

void Crs_BPUC(chromosome& chrom1,chromosome& chrom2,chromosome& temCh){
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double alpha = double(rand() % 100) / 100;
        if (alpha < Parameter_ADRKGA.BiasesRate) {
            temCh.Code_RK[i] = chrom1.Code_RK[i];
        } else {
            temCh.Code_RK[i] = chrom2.Code_RK[i];
        }
    }
}

double IFBD(chromosome& ch) {
    bool IsFrw = false;
    chromosome NewChrom = ch;
    chromosome OldChrom;

    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSort(ind, OldChrom.EndTime);
        set<double> st;
        vector<double> dect;
        double a = 0;
        while(st.size() < comConst.NumOfTsk) {
            a = double(rand()%10000)/10000;
            st.insert(a);
        }
        dect.assign(st.begin(),st.end());
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.Code_RK[ind[comConst.NumOfTsk - 1 - i]] = floor(NewChrom.Code_RK[ind[comConst.NumOfTsk - 1 - i]]) + dect[i];
        }
        NrmDcd(NewChrom, IsFrw);
        IsFrw = !IsFrw;
    } while (NewChrom.FitnessValue + PrecisionValue < OldChrom.FitnessValue);
    if (IsFrw) { //the last is backward
        ch = OldChrom;
    } else {
        ch = NewChrom;
    }
    return ch.FitnessValue;
}

//{ Load Balancing with Communication Reduction Improvement (LBCRI)}
void LBCA(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Id(comConst.NumOfRsc,0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = floor(ch.Code_RK[i]);
        Id[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc;
        TSK[RscIndex].push_back(i);
    }
    vector<int> ind(comConst.NumOfRsc);
    IndexSort(ind, Id);         //sorting according to loads
    int RscWithMinLd = ind[0];          //find out the resource (Rsc) with the lowest load;
    set<int> ST;
    if (abs(Id[RscWithMinLd]) < PrecisionValue) {
        ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
    } else {
        //traverse the tasks allocated to the resource with lowest load and add their parents and children to set ST
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            int TaskIndex = TSK[RscWithMinLd][i];
            ST.insert(Tasks[TaskIndex].children.begin(),Tasks[TaskIndex].children.end());
            ST.insert(Tasks[TaskIndex].parents.begin(),Tasks[TaskIndex].parents.end());
        }
        //delete the tasks which have been allocated the resource with lowest load
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            ST.erase(TSK[RscWithMinLd][i]);
        }
        //delete the tasks which can not be performed by the resource with lowest load
        for (auto iter = ST.begin(); iter != ST.end();) {
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) ==
                Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if(ST.empty()){
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s:ST) {
        t.push_back(pair<int, double>(s, Id[floor(ch.Code_RK[s])]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    double decimal = ch.Code_RK[t[0].first] - floor(ch.Code_RK[t[0].first]);
    ch.Code_RK[t[0].first] = RscWithMinLd + decimal;
    NrmDcd(ch, true);
    IFBD(ch);
    if (OldCh.FitnessValue + PrecisionValue < ch.FitnessValue) {
        ch = OldCh;
    }
}

//{calculate the cumulative probabilities for the population whose chromosome have been sorted}
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A , int& NumOfChormPerPop) {
    for (int n = 0; n < NumOfChormPerPop; ++n) {
        A[n] = pow(RtOfSltPrb,NumOfChormPerPop-1-n) * (RtOfSltPrb - 1) / (pow(RtOfSltPrb, NumOfChormPerPop) - 1);
    }
    for (int n = 1; n < NumOfChormPerPop; ++n){
        A[n] = A[n] + A[n - 1];
    }
}

//{select a chromosome using roulette wheel selection scheme}
int SelectChrom(vector<double>& A) {
    double lambda = RandomDouble(0, 1);
    for (int n = 0; n < A.size(); ++n)
        if (lambda <= A[n])
            return n;
}

