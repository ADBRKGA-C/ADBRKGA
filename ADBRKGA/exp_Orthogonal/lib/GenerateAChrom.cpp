#include <cstdlib>
#include "GenerateAChrom.h"
#include "GenOperator.h"
#include "tools.hpp"

//{calculate the average execution time of tasks}
void W_Cal_Average(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0;
        int RscSize = Tasks[i].ElgRsc.size();
        for (int j = 0; j < RscSize; ++j)
            sum += 1.0 / Rscs[Tasks[i].ElgRsc[j]].pc;
        w[i] = Tasks[i].length * sum / RscSize;
    }
}

//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Tasks[i].parents.size() == 0){
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    } else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

//calculate the rank of tasks based on independent IO using transfer time C[i][j]
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w){
    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size()-1].size(); ++i) {
        int TaskId=TskLstInLvl[TskLstInLvl.size()-1][i];
        RankList[TaskId] = w[TaskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId=TskLstInLvl[i][j];
            double ChildMaxRankc = 0;
            for (int k = 0; k < Tasks[TaskId].children.size(); ++k) {
                int tem = Tasks[TaskId].children[k];
                double CompareObject = RankList[tem] + c[TaskId][tem];
                if(ChildMaxRankc  < CompareObject ){
                    ChildMaxRankc = CompareObject;
                }
            }
            RankList[TaskId] =w[TaskId] + ChildMaxRankc;
        }
    }
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk);
    chrom.Code_RK.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk);
    chrom.EndTime.resize(comConst.NumOfTsk);
}


double HrsDcd_EFT(chromosome& ch) {

    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    vector<int > upr(comConst.NumOfTsk,0.0);
    vector<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double decimal ;
        int RscIndex = -1;
        double tmp = 1;
        int TaskIndex;
        for(int j = 0; j < RTI.size(); ++j) {
            decimal = ch.Code_RK[RTI[j]] - floor(ch.Code_RK[RTI[j]]);
            if(decimal < tmp) {
                tmp = decimal;
                TaskIndex = RTI[j];
            }
        }
        double FinalEndTime = 100000000000;
        double FinalStartTime = 0;
        for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
            double ReadyTime = 0;
            int v = Tasks[TaskIndex].ElgRsc[j];
            if(Tasks[TaskIndex].parents.size() != 0){
                for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) {
                    int ParentIndex = Tasks[TaskIndex].parents[n];
                    int ParentRscIndex = floor(ch.Code_RK[ParentIndex]);
                    double max = ch.EndTime[ParentIndex];
                    if(v != ParentRscIndex){
                        double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                        max += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                    }
                    if (ReadyTime < max){
                        ReadyTime = max;
                    }
                }
            }
            double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc;
            double StartTime = 0;
            double EndTime = 0;
            //{Find an idle time-slot as early as possible from ITL}
            set<double>::iterator pre  = ITL[v].begin();
            set<double>::iterator post = ITL[v].begin();
            ++post;
            while(post != ITL[v].end()) {
                if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
                    StartTime = XY_MAX(*pre, ReadyTime);
                    break;
                } else {
                    ++pre;
                    ++pre;
                    ++post;
                    ++post;
                }
            }
            EndTime = StartTime + ExeTime;
            //{find/record the earliest finish time}
            if (EndTime < FinalEndTime) {
                FinalStartTime = StartTime;
                FinalEndTime = EndTime;
                RscIndex = v;
            }
        }
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.Code_RK[TaskIndex] = RscIndex + tmp;
        RTI.erase(find(RTI.begin(), RTI.end(), TaskIndex));
        for (int l = 0; l < Tasks[TaskIndex].children.size(); ++l) {
            int childId = Tasks[TaskIndex].children[l];
            upr[childId] = upr[childId] - 1;
            if (upr[childId]==0){
                RTI.push_back(childId);
            }
        }
        //{update ITL}
        if(ITL[RscIndex].find(FinalStartTime) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(FinalStartTime);
        } else {
            ITL[RscIndex].insert(FinalStartTime);
        }
        if(ITL[RscIndex].find(ch.EndTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[TaskIndex]);
        }
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

chromosome GnrChrDHEFT(vector<double>& Rank_b){
    chromosome chrom;
    IntChr(chrom);
    set<double> st;
    vector<double> dect;
    double a = 0;
    while(st.size() < comConst.NumOfTsk) {
        a = double(rand() % 10000)/10000;
        st.insert(a);
    }
    dect.assign(st.begin(),st.end());

    vector<int > upr(comConst.NumOfTsk,0.0);
    vector<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }
    vector<double > RtSet(comConst.NumOfTsk,0.0);
    double max_EndTi = -100000;

    vector<set<double> > ITL;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }

    int taskIndex = 0;
    for(int q = 0 ; q < comConst.NumOfTsk ; ++q) {
        double max_RANK = -10000, priorty = 0;
        for (auto task:RTI) {
            priorty = Rank_b[task] + RtSet[task];
            if (priorty > max_RANK) {
                taskIndex = task;
                max_RANK = priorty;
            }
        }
        int vmIndex = -1;

        double finalEndTime = 100000000000, finalStartTime = 0;
        for (int j = 0; j < Tasks[taskIndex].ElgRsc.size(); ++j) {
            int v = Tasks[taskIndex].ElgRsc[j];
            double rt = 0;
            double ReadyTime = 0;
            if(Tasks[taskIndex].parents.size() != 0){
                for (int i = 0; i < Tasks[taskIndex].parents.size(); ++i) {
                    int parID = Tasks[taskIndex].parents[i];
                    int parRsc = floor(chrom.Code_RK[parID]);
                    double trantime = 0;
                    if( v != parRsc){
                        double sum = ParChildTranFileSizeSum[parID][taskIndex];
                        trantime = sum / VALUE / (XY_MIN(Rscs[v].bw, Rscs[parRsc].bw) / 8);
                    }
                    double maxfc = chrom.EndTime[parID] + trantime;
                    if(rt < maxfc){
                        rt = maxfc;
                    }
                }
            }
            ReadyTime = rt;
            double ExeTime = Tasks[taskIndex].length / Rscs[v].pc;
            double startTime = 0, endTime = 0;
            set<double>::iterator pre  = ITL[v].begin();
            set<double>::iterator post = ITL[v].begin();
            ++post;
            while(post != ITL[v].end()) {
                if((*post - *pre) >= ExeTime &&
                   ReadyTime <= (*post)-ExeTime) {;
                    startTime = XY_MAX(*pre, ReadyTime);
                    break;
                } else {
                    ++pre;
                    ++pre;
                    ++post;
                    ++post;
                }
            }
            endTime = startTime + ExeTime;
            if (endTime < finalEndTime) {
                finalStartTime = startTime;
                finalEndTime = endTime;
                vmIndex = v;
            }
        }
        chrom.Code_RK[taskIndex] = dect[q] + vmIndex;
        chrom.EndTime[taskIndex] = finalEndTime;
        max_EndTi = XY_MAX(max_EndTi, finalEndTime);

        if(ITL[vmIndex].find(finalStartTime) != ITL[vmIndex].end()) {
            ITL[vmIndex].erase(finalStartTime);
        } else {
            ITL[vmIndex].insert(finalStartTime);
        }
        if(ITL[vmIndex].find(chrom.EndTime[taskIndex]) != ITL[vmIndex].end()) {
            ITL[vmIndex].erase(chrom.EndTime[taskIndex]);
        } else {
            ITL[vmIndex].insert(chrom.EndTime[taskIndex]);
        }
        RTI.erase(find(RTI.begin(), RTI.end(), taskIndex));
        for (int l = 0; l < Tasks[taskIndex].children.size(); ++l) {
            int childid =Tasks[taskIndex].children[l];
            upr[childid]=upr[childid] -1;
            if (upr[childid]==0){
                RTI.push_back(childid);
                if(Tasks[childid].parents.size() != 0) {
                    int j = Tasks[childid].ElgRsc.size();
                    for (int i = 0; i < Tasks[childid].parents.size(); ++i) {
                        int parent = Tasks[childid].parents[i];
                        int parentRsc = floor(chrom.Code_RK[parent]);
                        double transfertime = 0;
                        for (int k = 0; k < j; ++k) {
                            int v = Tasks[childid].ElgRsc[k];
                            if (v != parentRsc) {
                                double sum = ParChildTranFileSizeSum[parent][childid];
                                transfertime +=
                                        sum / VALUE / (XY_MIN(Rscs[v].bw, Rscs[parentRsc].bw) / 8);
                            }
                        }
                        double max = transfertime / j + chrom.EndTime[parent];
                        if (RtSet[childid] < max) {
                            RtSet[childid] = max;
                        }
                    }
                }
            }
        }
    }
    chrom.FitnessValue = max_EndTi;
    return chrom;
}

chromosome GnrChr_HEFT_ADBRKGA(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    set<double> st;
    vector<double> dect;
    double a = 0;
    while(st.size() < comConst.NumOfTsk) {
        a = double(rand() % 10000)/10000;
        st.insert(a);
    }
    dect.assign(st.begin(),st.end());
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.Code_RK[ind[comConst.NumOfTsk - i - 1]] = dect[i];
    }
    HrsDcd_EFT(TemChrom);
    return TemChrom;
}

void GnrChr_Ran(chromosome& chrom) {
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[i] = double(rand()%10000)/10000 + Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()];
    }
}
