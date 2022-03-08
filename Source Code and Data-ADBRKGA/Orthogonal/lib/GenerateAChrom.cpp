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
            RankList[TaskId] = w[TaskId] + ChildMaxRankc;
        }
    }
}

void Calculate_Rank_t(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w) {
    for(int i =1 ;i < TskLstInLvl.size(); ++i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId = TskLstInLvl[i][j];
            for (int k = 0; k < Tasks[TaskId].parents.size(); ++k) {
                int tem = Tasks[TaskId].parents[k];
                double re = w[tem] + c[tem][TaskId] + RankList[tem];
                if (RankList[TaskId] < re) {
                    RankList[TaskId] = re;
                }
            }
        }
    }
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk);
    chrom.Code_RK.resize(comConst.NumOfTsk);
    chrom.RscAlcPart.resize(comConst.NumOfTsk);
    chrom.TskSchPart.resize(comConst.NumOfTsk);
    chrom.VTskSchPart.resize(comConst.NumOfTsk,0.0);
    chrom.VRscAlcPart.resize(comConst.NumOfTsk,0.0);
    chrom.EndTime.resize(comConst.NumOfTsk);
}

vector<double> GnrDecimalsByAscend() {
    vector<double> decimals(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        decimals[i] = (i + RandomDouble(0,1)) / comConst.NumOfTsk;
    }
    return decimals;
}

// I/O independent
double DcdEvl(chromosome& ch, bool IsFrw) {
    double makespan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i];
        int RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the resource (Rsc) allocated to the task
        double ReadyTime = 0;
        if(IsFrw) {                              //forward-loading
            for (int j = 0; j < Tasks[TaskIndex].parents.size(); ++j) {
                int ParentTask = Tasks[TaskIndex].parents[j];
                int ParentRsc = ch.RscAlcLst[ParentTask];
                double fft = ch.EndTime[ParentTask];
                if(RscIndex != ParentRsc) {
                    fft += ParChildTranFileSizeSum[ParentTask][TaskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw)); // -xy
                }
                if (ReadyTime < fft) {
                    ReadyTime = fft;
                }
            }
        } else {                                //backward-loading
            for (int j = 0; j < Tasks[TaskIndex].children.size(); ++j) {
                int ChildTask = Tasks[TaskIndex].children[j];
                int ChildRsc = ch.RscAlcLst[ChildTask];
                double fft = ch.EndTime[ChildTask];
                if(RscIndex != ChildRsc) {
                    fft += ParChildTranFileSizeSum[TaskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                }
                if (ReadyTime < fft) {
                    ReadyTime = fft;
                }
            }
        }
        double ExecutionTime = Tasks[TaskIndex].length / Rscs[RscIndex].pc;
        double StartTime = FindIdleTimeSlot(ITL[RscIndex],ExecutionTime,ReadyTime); //{find an idle time-slot in ITL which can finish the task  at the earliest}
        ch.EndTime[TaskIndex] = StartTime + ExecutionTime;
        if (makespan < ch.EndTime[TaskIndex]) {
            makespan = ch.EndTime[TaskIndex];
        }
        UpdateITL(ITL[RscIndex],StartTime,ch.EndTime[TaskIndex]);       //{update ITL}
    }
    ch.FitnessValue = makespan;
    return ch.FitnessValue;
}

double NrmDcd(chromosome& ch, bool IsFrw) {
    vector<int > upr(comConst.NumOfTsk,-1);
    list<int> RTI;
    if(IsFrw)
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].parents.size();
            if (upr[i]==0)  RTI.push_back(i);
        }
    else
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].children.size();
            if (upr[i]==0)  RTI.push_back(i);
        }
    //generate resource allocation list and task scheduling order list
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.RscAlcLst[i] = floor(ch.Code_RK[i]);
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            int decimal = ch.Code_RK[*lit] - floor(ch.Code_RK[*lit]);
            if (decimal < tmp) {
                tmp = decimal; pit = lit; //小数部分最小的那个优先调度
            }
        }
        ch.TskSchLst[i] = *pit;
        RTI.erase(pit);
        //更新RTI;
        if (IsFrw)
            for (int l = 0; l < Tasks[ch.TskSchLst[i]].children.size(); ++l) {
                int childId = Tasks[ch.TskSchLst[i]].children[l];
                upr[childId] = upr[childId] - 1;
                if (upr[childId]==0)   RTI.push_back(childId);
            }
        else
            for (int l = 0; l < Tasks[ch.TskSchLst[i]].parents.size(); ++l) {
                int parentId = Tasks[ch.TskSchLst[i]].parents[l];
                upr[parentId] = upr[parentId] - 1;
                if (upr[parentId]==0)  RTI.push_back(parentId);
            }
    }
    DcdEvl(ch, IsFrw);
    return ch.FitnessValue;
}

double HrsDcd_CTP(chromosome& ch) {
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
    IndexSortByValueOnAscend(ind, Rank_b);
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.Code_RK[TaskIndex] = floor(ch.Code_RK[TaskIndex]) + Decimals[i];
    }
    NrmDcd(ch, true);
    return ch.FitnessValue;
}

double IFBS(chromosome& ch) {
    bool IsFrw = false;
    chromosome NewChrom = ch;
    chromosome OldChrom;
    vector<double> Decimals = GnrDecimalsByAscend();
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSortByValueOnAscend(ind, OldChrom.EndTime);
        //vector<double> Decimals = GnrDecimalsByAscend();
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.Code_RK[ind[comConst.NumOfTsk - 1 - i]] = floor(NewChrom.Code_RK[ind[comConst.NumOfTsk - 1 - i]]) + Decimals[i];
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
void LBCA_IFBS(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Ld(comConst.NumOfRsc,0.0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = floor(ch.Code_RK[i]);
        Ld[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc;
        TSK[RscIndex].push_back(i);
    }
    int RscWithMinLd = 0;
    double TemLd = Ld[0];
    for (int j = 1; j < comConst.NumOfRsc; ++j) { //find out the resource (Rsc) with the lowest load -new-xy;
        if (TemLd > Ld[j]) {
            TemLd = Ld[j]; RscWithMinLd = j;
        }
    }
    set<int> ST;
    if (abs(Ld[RscWithMinLd]) < PrecisionValue) {
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
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) ==  Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if(ST.empty()){//-w
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s:ST) {
        t.push_back(pair<int, double>(s, Ld[floor(ch.Code_RK[s])]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    double decimal = ch.Code_RK[t[0].first] - floor(ch.Code_RK[t[0].first]);
    ch.Code_RK[t[0].first] = RscWithMinLd + decimal;
    NrmDcd(ch, true);
    IFBS(ch);
    if (OldCh.FitnessValue + PrecisionValue < ch.FitnessValue) {
        ch = OldCh;
    }
}

double HrsDcd_EFT(chromosome& ch) {
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    vector<int > upr(comConst.NumOfTsk,0.0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0)  RTI.push_back(i);
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1;
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            double decimal = ch.Code_RK[*lit] - floor(ch.Code_RK[*lit]);
            if (decimal < tmp) {
                tmp = decimal; pit = lit; //小数部分最小的那个优先调度
            }
        }
        ch.TskSchLst[i] = *pit;
        RTI.erase(pit);
        double FinalEndTime = InfiniteValue;
        double FinalStartTime = 0;
        SeletRsc_EFT(ch,ITL,ch.TskSchLst[i],RscId,FinalStartTime,FinalEndTime);
        ch.EndTime[ch.TskSchLst[i]] = FinalEndTime;
        ch.Code_RK[ch.TskSchLst[i]] = RscId + tmp;
        ch.RscAlcLst[ch.TskSchLst[i]] = RscId;
        UpdateITL(ITL[RscId],FinalStartTime,FinalEndTime); //{update ITL}
        makespan = XY_MAX(makespan, FinalEndTime);
        for (int l = 0; l < Tasks[ch.TskSchLst[i]].children.size(); ++l) {
            int ChildId = Tasks[ch.TskSchLst[i]].children[l];
            upr[ChildId] = upr[ChildId] - 1;
            if (upr[ChildId]==0)   RTI.push_back(ChildId);
        }
    }
    ch.FitnessValue = makespan;
    return makespan;
}

chromosome GnrChr_HEFT_b_ADBRKGA(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    vector<double> Decimals = GnrDecimalsByAscend();
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.Code_RK[ind[comConst.NumOfTsk - i - 1]] = Decimals[i];
    }
    HrsDcd_EFT(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HEFT_t_ADBRKGA(vector<double> Rank_t) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    vector<double> Decimals = GnrDecimalsByAscend();
    IndexSortByValueOnAscend(ind, Rank_t);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.Code_RK[ind[i]] = Decimals[i];
    }
    HrsDcd_EFT(TemChrom);
    return TemChrom;
}

chromosome GnrChr_IHEFT3_b(vector<double> Rank_b) {
    chromosome Chrom_Rank_b;
    IntChr(Chrom_Rank_b);
    IndexSortByValueOnDescend(Chrom_Rank_b.TskSchLst, Rank_b);
    IHEFT3(Chrom_Rank_b);
    return Chrom_Rank_b;
}

chromosome GnrChr_IHEFT3_t(vector<double> Rank_t) {
    chromosome Chrom_Rank_t;
    IntChr(Chrom_Rank_t);
    IndexSortByValueOnAscend(Chrom_Rank_t.TskSchLst, Rank_t);
    IHEFT3(Chrom_Rank_t);
    return Chrom_Rank_t;
}

double IHEFT3(chromosome& ch) {
    list <int> TemTskSchLst;
    TemTskSchLst.assign(ch.TskSchLst.begin(),ch.TskSchLst.end());
    ch.RscAlcLst.resize(comConst.NumOfTsk,-1);
    ch.TskSchLst.resize(comConst.NumOfTsk,-1);
    vector<set<double>> ITL;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    vector<int> upr(comConst.NumOfTsk,0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
    }
    int IndexCount = 0;
    while (!TemTskSchLst.empty()){
        int CrnTask = TemTskSchLst.front();
        ch.TskSchLst[IndexCount] = CrnTask;
        IndexCount++;
        TemTskSchLst.erase(TemTskSchLst.begin());
        int FinalRscForCrnTask = -1;
        double FinalStartTimeOfCrnTask = 0;
        double FinalEndTimeOfCrnTask = InfiniteValue;
        vector<int> NeedProcessChildTaskSet;
        for (int m = 0; m < Tasks[CrnTask].children.size(); ++m) {
            int childId = Tasks[CrnTask].children[m];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                NeedProcessChildTaskSet.push_back(childId);
            }
        }
        if(NeedProcessChildTaskSet.empty()){
            SeletRsc_EFT(ch,ITL,CrnTask,FinalRscForCrnTask,FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
            ch.EndTime[CrnTask] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTask] = FinalRscForCrnTask;
            ch.FitnessValue = XY_MAX(ch.FitnessValue,ch.EndTime[CrnTask]);
            UpdateITL(ITL[FinalRscForCrnTask],FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
        } else{
            int FinalChildTask = -1 , FinalRscForChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = InfiniteValue;
            vector<set<double> > ITL_CrnTaskScheduled;
            for (int j = 0; j < Tasks[CrnTask].ElgRsc.size(); ++j) {
                double ReadyTimeOfCrnTask = 0;
                int CrnTaskRsc = Tasks[CrnTask].ElgRsc[j];
                ch.RscAlcLst[CrnTask] = CrnTaskRsc;
                for (int i2 = 0; i2 < Tasks[CrnTask].parents.size(); ++i2) {
                    int ParentTask = Tasks[CrnTask].parents[i2];
                    int RscOfParentTask = ch.RscAlcLst[ParentTask];
                    double max = ch.EndTime[ParentTask];
                    if(CrnTaskRsc != RscOfParentTask){
                        max += ParChildTranFileSizeSum[ParentTask][CrnTask] / VALUE * 8 / (XY_MIN(Rscs[CrnTaskRsc].bw,Rscs[RscOfParentTask].bw));
                    }
                    if (ReadyTimeOfCrnTask + PrecisionValue < max){
                        ReadyTimeOfCrnTask = max;
                    }
                }
                double ExeTimeOfCrnTask = Tasks[CrnTask].length / Rscs[CrnTaskRsc].pc;
                double CrnTaskStartTime = FindIdleTimeSlot(ITL[CrnTaskRsc],ExeTimeOfCrnTask,ReadyTimeOfCrnTask);
                double CrnTaskEndTime = CrnTaskStartTime + ExeTimeOfCrnTask;
                vector<set<double> > TemITL = ITL;
                UpdateITL(TemITL[CrnTaskRsc],CrnTaskStartTime,CrnTaskEndTime);
                ch.EndTime[CrnTask] = CrnTaskEndTime;
                for (int i3 = 0; i3 < NeedProcessChildTaskSet.size(); ++i3) {
                    int TemChildTask = NeedProcessChildTaskSet[i3];
                    for (int j1 = 0; j1 < Tasks[TemChildTask].ElgRsc.size(); ++j1) {
                        double TemChildReadyTime = 0;
                        int TemChildRsc = Tasks[TemChildTask].ElgRsc[j1];
                        for (int i4 = 0; i4 < Tasks[TemChildTask].parents.size(); ++i4) {
                            int TemParentTask = Tasks[TemChildTask].parents[i4];
                            int TemParRsc = ch.RscAlcLst[TemParentTask];
                            double max = ch.EndTime[TemParentTask];
                            if (TemChildRsc != TemParRsc) {
                                max += ParChildTranFileSizeSum[TemParentTask][TemChildTask] / VALUE * 8 / (XY_MIN(Rscs[TemChildRsc].bw, Rscs[TemParRsc].bw));
                            }
                            if (TemChildReadyTime + PrecisionValue < max) {
                                TemChildReadyTime = max;
                            }
                        }
                        double TemChildExeTime = Tasks[TemChildTask].length / Rscs[TemChildRsc].pc;
                        double TemChildStartTime = FindIdleTimeSlot(TemITL[TemChildRsc],TemChildExeTime,TemChildReadyTime);
                        double TemChildEndTime = TemChildStartTime + TemChildExeTime;
                        if (FinalEndTimeOfChildTask > TemChildEndTime + PrecisionValue ) {
                            FinalEndTimeOfChildTask = TemChildEndTime;
                            FinalRscForChildTask = TemChildRsc;
                            FinalChildTask = TemChildTask;
                            FinalStartTimeOfChildTask = TemChildStartTime;
                            FinalEndTimeOfCrnTask = CrnTaskEndTime;
                            FinalRscForCrnTask = CrnTaskRsc;
                            ITL_CrnTaskScheduled = TemITL;
                        }
                    }
                }
            }
            ch.EndTime[CrnTask] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTask] = FinalRscForCrnTask;
            ch.TskSchLst[IndexCount] = FinalChildTask;
            IndexCount++;
            ch.RscAlcLst[FinalChildTask] = FinalRscForChildTask;
            ch.EndTime[FinalChildTask] = FinalEndTimeOfChildTask;
            ch.FitnessValue = XY_MAX(ch.FitnessValue,ch.EndTime[FinalChildTask]);
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL[FinalRscForChildTask],FinalStartTimeOfChildTask,ch.EndTime[FinalChildTask]);
            TemTskSchLst.erase(find(TemTskSchLst.begin(),TemTskSchLst.end(),FinalChildTask));
            for (int m = 0; m < Tasks[FinalChildTask].children.size(); ++m) {
                int childId = Tasks[FinalChildTask].children[m];
                upr[childId] = upr[childId] - 1;
            }
        }
    }
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.Code_RK[ch.TskSchLst[i]] = ch.RscAlcLst[ch.TskSchLst[i]] + Decimals[i];
    }
    return ch.FitnessValue;
}

chromosome GnrChr_DHEFT(vector<double>& Rank_b) {
    chromosome chrom;
    IntChr(chrom);
    vector<int> upr(comConst.NumOfTsk,0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }
    vector<set<double>> ITL;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    int IndexCount = 0;
    vector<double > AvrRtSet(comConst.NumOfTsk,0.0);
    while (!RTI.empty()){
        double max_RANK = 0;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            double priority = Rank_b[*lit] + AvrRtSet[*lit];
            if (priority - PrecisionValue > max_RANK) {
                pit = lit; max_RANK = priority;
            }
        }
        int CrnTskId = *pit;
        RTI.erase(pit);
        chrom.TskSchLst[IndexCount] = CrnTskId;
        IndexCount++;
        int FinalRscForCrnTask = -1, NeedProcessChildTask = -1;
        double FinalStartTimeOfCrnTask = 0, FinalEndTimeOfCrnTask = InfiniteValue, TemMaxRnk = -1;
        vector<int> ReadyChildTaskSet;
        for (int m = 0; m < Tasks[CrnTskId].children.size(); ++m) {
            int childId = Tasks[CrnTskId].children[m];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                if (Rank_b[childId] - PrecisionValue > TemMaxRnk) {
                    NeedProcessChildTask = childId;
                    TemMaxRnk = Rank_b[NeedProcessChildTask];
                }
                ReadyChildTaskSet.push_back(childId);
            }
        }
        if(NeedProcessChildTask == -1){ //只需要处理当前任务
            SeletRsc_EFT(chrom,ITL,CrnTskId,FinalRscForCrnTask,FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
            chrom.EndTime[CrnTskId] = FinalEndTimeOfCrnTask;
            chrom.RscAlcLst[CrnTskId] = FinalRscForCrnTask;
            chrom.FitnessValue = XY_MAX(chrom.FitnessValue,chrom.EndTime[CrnTskId]);
            UpdateITL(ITL[FinalRscForCrnTask],FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
        } else{  //需要和可以处理的子任务一起考虑分配资源
            int FinalRscForChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = InfiniteValue;
            vector<set<double>> ITL_CrnTaskScheduled;
            for (int j = 0; j < Tasks[CrnTskId].ElgRsc.size(); ++j) {
                double ReadyTimeOfCrnTask = 0;
                int RscForCrnTsk = Tasks[CrnTskId].ElgRsc[j];
                chrom.RscAlcLst[CrnTskId] = RscForCrnTsk;
                for (int i2 = 0; i2 < Tasks[CrnTskId].parents.size(); ++i2) {
                    int ParentTask = Tasks[CrnTskId].parents[i2];
                    int RscOfParentTask = chrom.RscAlcLst[ParentTask];
                    double ftt = chrom.EndTime[ParentTask];
                    if(RscForCrnTsk != RscOfParentTask){
                        ftt += ParChildTranFileSizeSum[ParentTask][CrnTskId] / VALUE * 8 / (XY_MIN(Rscs[RscForCrnTsk].bw,Rscs[RscOfParentTask].bw));
                    }
                    if (ReadyTimeOfCrnTask + PrecisionValue < ftt){
                        ReadyTimeOfCrnTask = ftt;
                    }
                }
                double ExeTimeOfCrnTask = Tasks[CrnTskId].length / Rscs[RscForCrnTsk].pc;
                double CrnTaskStartTime = FindIdleTimeSlot(ITL[RscForCrnTsk],ExeTimeOfCrnTask,ReadyTimeOfCrnTask);
                double CrnTaskEndTime = CrnTaskStartTime + ExeTimeOfCrnTask;
                vector<set<double> > TemITL = ITL;
                UpdateITL(TemITL[RscForCrnTsk],CrnTaskStartTime,CrnTaskEndTime);
                chrom.EndTime[CrnTskId] = CrnTaskEndTime;
                for (int j1 = 0; j1 < Tasks[NeedProcessChildTask].ElgRsc.size(); ++j1) {
                    double ChildReadyTime = 0;
                    int RscOfChildTsk = Tasks[NeedProcessChildTask].ElgRsc[j1];
                    for (int i4 = 0; i4 < Tasks[NeedProcessChildTask].parents.size(); ++i4) {
                        int TemParentTask = Tasks[NeedProcessChildTask].parents[i4];
                        int RscOfTemParTsk = chrom.RscAlcLst[TemParentTask];
                        double ftt = chrom.EndTime[TemParentTask];
                        if (RscOfChildTsk != RscOfTemParTsk) {
                            ftt += ParChildTranFileSizeSum[TemParentTask][NeedProcessChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscOfChildTsk].bw, Rscs[RscOfTemParTsk].bw));
                        }
                        if (ChildReadyTime + PrecisionValue < ftt) {
                            ChildReadyTime = ftt;
                        }
                    }
                    double ChildExeTime = Tasks[NeedProcessChildTask].length / Rscs[RscOfChildTsk].pc;
                    double ChildStartTime = FindIdleTimeSlot(TemITL[RscOfChildTsk],ChildExeTime,ChildReadyTime);
                    double ChildEndTime = ChildStartTime + ChildExeTime;
                    if (FinalEndTimeOfChildTask > ChildEndTime + PrecisionValue ) {
                        FinalEndTimeOfChildTask = ChildEndTime;
                        FinalRscForChildTask = RscOfChildTsk;
                        FinalStartTimeOfChildTask = ChildStartTime;
                        FinalEndTimeOfCrnTask = CrnTaskEndTime;
                        FinalRscForCrnTask = RscForCrnTsk;
                        ITL_CrnTaskScheduled = TemITL;
                    }
                }
            }
            chrom.EndTime[CrnTskId] = FinalEndTimeOfCrnTask;
            chrom.RscAlcLst[CrnTskId] = FinalRscForCrnTask;
            chrom.TskSchLst[IndexCount] = NeedProcessChildTask;
            IndexCount++;
            chrom.RscAlcLst[NeedProcessChildTask] = FinalRscForChildTask;
            chrom.EndTime[NeedProcessChildTask] = FinalEndTimeOfChildTask;
            chrom.FitnessValue = XY_MAX(chrom.FitnessValue,chrom.EndTime[NeedProcessChildTask]);
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL[FinalRscForChildTask],FinalStartTimeOfChildTask,chrom.EndTime[NeedProcessChildTask]);
            for (int q = 0; q < ReadyChildTaskSet.size(); ++q) { //把所有还没有处理的就绪子任务添加到RTI中；
                int TskId = ReadyChildTaskSet[q];
                if (TskId != NeedProcessChildTask) {
                    RTI.push_back(TskId);
                    AvrRtSet[TskId] = ClcAvrReadyTime(TskId, chrom);
                }
            }
            for (int m = 0; m < Tasks[NeedProcessChildTask].children.size(); ++m) {//把就绪的子任务的子任务添加到RTI中；
                int childId = Tasks[NeedProcessChildTask].children[m];
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    RTI.push_back(childId);
                    AvrRtSet[childId] = ClcAvrReadyTime(childId, chrom);
                }
            }
        }
    }
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[chrom.TskSchLst[i]] = chrom.RscAlcLst[chrom.TskSchLst[i]] + Decimals[i];
    }
    return chrom;
}

double ClcAvrReadyTime(int TskId, chromosome& chrom) {
    int j = Tasks[TskId].ElgRsc.size();
    double AvrRt = 0;
    for (int i = 0; i < Tasks[TskId].parents.size(); ++i) {
        int parent = Tasks[TskId].parents[i];
        int parentRsc = chrom.RscAlcLst[parent];
        double transfertime = 0;
        for (int k = 0; k < j; ++k) {
            int v = Tasks[TskId].ElgRsc[k];
            if (v != parentRsc) {
                transfertime += ParChildTranFileSizeSum[parent][TskId] / VALUE / (XY_MIN(Rscs[v].bw, Rscs[parentRsc].bw) / 8);
            }
        }
        double fft = transfertime / j + chrom.EndTime[parent];
        if (AvrRt < fft) {
            AvrRt = fft;
        }
    }
    return AvrRt;
}

void SeletRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskId, int& RscId, double& FinalStartTime, double& FinalEndTime) {
    for (int j = 0; j < Tasks[TaskId].ElgRsc.size(); ++j) {
        double ReadyTime = 0;
        int RscIdOfCrnTsk = Tasks[TaskId].ElgRsc[j];
        for (int n = 0; n < Tasks[TaskId].parents.size(); ++n) { //calculate the ready time of the task
            int PrnTskId = Tasks[TaskId].parents[n];
            int RscIdOfPrnTsk = ch.RscAlcLst[PrnTskId];
            double fft = ch.EndTime[PrnTskId];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                double TransferData = ParChildTranFileSizeSum[PrnTskId][TaskId];
                fft += TransferData / VALUE * 8 / (XY_MIN(Rscs[RscIdOfCrnTsk].bw,Rscs[RscIdOfPrnTsk].bw));
            }
            if (ReadyTime + PrecisionValue < fft){
                ReadyTime = fft;
            }
        }
        double ExeTime = Tasks[TaskId].length / Rscs[RscIdOfCrnTsk].pc;
        double StartTime = FindIdleTimeSlot(ITL[RscIdOfCrnTsk],ExeTime,ReadyTime); //Find an idle time-slot as early as possible from ITL
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + PrecisionValue < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscId = RscIdOfCrnTsk;
        }
    }
}

double FindIdleTimeSlot(set<double>& ITLofRscId,double& ExeTime,double& ReadyTime){
    set<double>::iterator pre  = ITLofRscId.begin();
    set<double>::iterator post = ITLofRscId.begin();
    ++post;
    while(post != ITLofRscId.end()) {
        if((*post - *pre) > ExeTime - PrecisionValue && ReadyTime - PrecisionValue < (*post)-ExeTime) {
            return  XY_MAX(*pre, ReadyTime);
        } else {
            ++pre; ++pre; ++post; ++post;
        }
    }
}

void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime){
    if(ITLofRscId.find(StartTime) != ITLofRscId.end()) {
        ITLofRscId.erase(StartTime);
    } else {
        ITLofRscId.insert(StartTime);
    }
    if(ITLofRscId.find(EndTime) != ITLofRscId.end()) {
        ITLofRscId.erase(EndTime);
    } else {
        ITLofRscId.insert(EndTime);
    }
}

chromosome GnrChr_Lvl_Ran() {
    chromosome chrom;
    IntChr(chrom);
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()] + (LevelIdOfTask[i] + (rand() % 10000) / 10000.0) / TskLstInLvl.size();
    }
    return chrom;
}

