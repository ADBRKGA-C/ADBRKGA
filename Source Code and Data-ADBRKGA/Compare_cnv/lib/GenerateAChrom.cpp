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

//{generate a topological sort randomly}
vector<int> GnrSS_TS() {
    vector<int> SS;
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }
    while (!RTI.empty()) {
        int RandVec = rand() % RTI.size();
        int v = RTI[RandVec];
        vector<int>::iterator iter = RTI.begin() + RandVec;
        RTI.erase(iter);
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            --upr[Tasks[v].children[i]];
            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
        }
        SS.push_back(v);
    }
    return SS;
}

//{generate a task scheduling order by the levels of tasks from small to large} -xy2
//{Those haveing the same level are ranked arbitrarily among them} -xy2
vector<int> GnrSS_Lvl() {
    vector<int> ch;
    vector<vector<int>> tem = TskLstInLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        random_shuffle(tem[i].begin(), tem[i].end());   //arrange the tasks in each level
        for (int j = 0; j < tem[i].size(); ++j) {
            ch.push_back(tem[i][j]);
        }
    }
    return ch;
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

double GnrMS_Evl(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1;
        int TaskId = ch.TskSchLst[i];
        double FinalEndTime = InfiniteValue;
        double FinalStartTime = 0;
        SeletRsc_EFT(ch,ITL,TaskId,RscId,FinalStartTime,FinalEndTime);
        ch.EndTime[TaskId] = FinalEndTime;
        ch.RscAlcLst[TaskId] = RscId;
        UpdateITL(ITL[RscId],FinalStartTime,FinalEndTime);     //{update ITL}
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
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

chromosome GnrChr_HEFT_b(vector<double> Rank_b) {
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnDescend(TemChrom.TskSchLst, Rank_b);
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HEFT_t(vector<double> Rank_t) {
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(TemChrom.TskSchLst, Rank_t);
    GnrMS_Evl(TemChrom);
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

chromosome GnrChr_IHEFT3(vector<double> Rank_b,vector<double> Rank_t) {
    chromosome Chrom_Rank_b = GnrChr_IHEFT3_b(Rank_b);
    chromosome Chrom_Rank_t = GnrChr_IHEFT3_t(Rank_t);
    if (Chrom_Rank_b.FitnessValue + PrecisionValue < Chrom_Rank_t.FitnessValue){
        return Chrom_Rank_b;
    } else{
        return Chrom_Rank_t;
    }
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

chromosome GnrChr_DIHEFT3(vector<double>& Rank_b) {
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
        int FinalRscForCrnTask = -1;
        double FinalStartTimeOfCrnTask = 0, FinalEndTimeOfCrnTask = InfiniteValue;
        vector<int> NeedProcessChildTaskSet;
        for (int m = 0; m < Tasks[CrnTskId].children.size(); ++m) {
            int childId = Tasks[CrnTskId].children[m];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                NeedProcessChildTaskSet.push_back(childId);
            }
        }
        if(NeedProcessChildTaskSet.empty()){ //只需要处理当前任务
            SeletRsc_EFT(chrom,ITL,CrnTskId,FinalRscForCrnTask,FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
            chrom.EndTime[CrnTskId] = FinalEndTimeOfCrnTask;
            chrom.RscAlcLst[CrnTskId] = FinalRscForCrnTask;
            chrom.FitnessValue = XY_MAX(chrom.FitnessValue,chrom.EndTime[CrnTskId]);
            UpdateITL(ITL[FinalRscForCrnTask],FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
        } else{  //需要和可以处理的子任务一起考虑分配资源
            int FinalChildTask = -1 , FinalRscForChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = InfiniteValue;
            vector<set<double> > ITL_CrnTaskScheduled;
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
                for (int i3 = 0; i3 < NeedProcessChildTaskSet.size(); ++i3) {
                    int TemChildTask = NeedProcessChildTaskSet[i3];
                    for (int j1 = 0; j1 < Tasks[TemChildTask].ElgRsc.size(); ++j1) {
                        double TemChildReadyTime = 0;
                        int RscOfTemChildTsk = Tasks[TemChildTask].ElgRsc[j1];
                        for (int i4 = 0; i4 < Tasks[TemChildTask].parents.size(); ++i4) {
                            int TemParentTask = Tasks[TemChildTask].parents[i4];
                            int RscOfTemParTsk = chrom.RscAlcLst[TemParentTask];
                            double ftt = chrom.EndTime[TemParentTask];
                            if (RscOfTemChildTsk != RscOfTemParTsk) {
                                ftt += ParChildTranFileSizeSum[TemParentTask][TemChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscOfTemChildTsk].bw, Rscs[RscOfTemParTsk].bw));
                            }
                            if (TemChildReadyTime + PrecisionValue < ftt) {
                                TemChildReadyTime = ftt;
                            }
                        }
                        double TemChildExeTime = Tasks[TemChildTask].length / Rscs[RscOfTemChildTsk].pc;
                        double TemChildStartTime = FindIdleTimeSlot(TemITL[RscOfTemChildTsk],TemChildExeTime,TemChildReadyTime);
                        double TemChildEndTime = TemChildStartTime + TemChildExeTime;
                        if (FinalEndTimeOfChildTask > TemChildEndTime + PrecisionValue ) {
                            FinalEndTimeOfChildTask = TemChildEndTime;
                            FinalRscForChildTask = RscOfTemChildTsk;
                            FinalChildTask = TemChildTask;
                            FinalStartTimeOfChildTask = TemChildStartTime;
                            FinalEndTimeOfCrnTask = CrnTaskEndTime;
                            FinalRscForCrnTask = RscForCrnTsk;
                            ITL_CrnTaskScheduled = TemITL;
                        }
                    }
                }
            }
            chrom.EndTime[CrnTskId] = FinalEndTimeOfCrnTask;
            chrom.RscAlcLst[CrnTskId] = FinalRscForCrnTask;
            chrom.TskSchLst[IndexCount] = FinalChildTask;
            IndexCount++;
            chrom.RscAlcLst[FinalChildTask] = FinalRscForChildTask;
            chrom.EndTime[FinalChildTask] = FinalEndTimeOfChildTask;
            chrom.FitnessValue = XY_MAX(chrom.FitnessValue,chrom.EndTime[FinalChildTask]);
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL[FinalRscForChildTask],FinalStartTimeOfChildTask,chrom.EndTime[FinalChildTask]);
            for (int q = 0; q < NeedProcessChildTaskSet.size(); ++q) { //把所有还没有处理的就绪子任务添加到RTI中；
                int TskId = NeedProcessChildTaskSet[q];
                if (TskId != FinalChildTask) {
                    RTI.push_back(TskId);
                    AvrRtSet[TskId] = ClcAvrReadyTime2(TskId, chrom);
                }
            }
            for (int m = 0; m < Tasks[FinalChildTask].children.size(); ++m) {//把就绪的子任务的子任务添加到RTI中；
                int childId = Tasks[FinalChildTask].children[m];
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    RTI.push_back(childId);
                    AvrRtSet[childId] = ClcAvrReadyTime2(childId, chrom);
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

double ClcAvrReadyTime2(int TskId, chromosome& chrom) {
    int j = Tasks[TskId].ElgRsc.size();
    double AvrRt = 0;
    for (int i = 0; i < Tasks[TskId].parents.size(); ++i) {
        int ParTsk = Tasks[TskId].parents[i];
        int RscIdOfParTsk = chrom.RscAlcLst[ParTsk];
        double transfertime = 0;
        for (int k = 0; k < j; ++k) {  //计算传输时间取最大
            int v = Tasks[TskId].ElgRsc[k];
            if (v != RscIdOfParTsk) {
                int TemTrnTime = ParChildTranFileSizeSum[ParTsk][TskId] / VALUE / (XY_MIN(Rscs[v].bw, Rscs[RscIdOfParTsk].bw) / 8);
                if (transfertime < TemTrnTime) transfertime = TemTrnTime;
            }
        }
        double fft = transfertime + chrom.EndTime[ParTsk];
        if (AvrRt < fft) {
            AvrRt = fft;
        }
    }
    return AvrRt;
}

chromosome GnrChr_DHEFT_old(vector<double>& Rank_b){
    chromosome chrom;
    IntChr(chrom);
    vector<double> Decimals = GnrDecimalsByAscend();
    vector<int> upr(comConst.NumOfTsk,0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)   RTI.push_back(i);
    }
    vector<double > AvrRtSet(comConst.NumOfTsk,0.0);
    double MakeSpan = -1;
    vector<set<double> > ITL;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);  a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    for(int q = 0 ; q < comConst.NumOfTsk ; ++q) {
        double max_RANK = 0;
        int CrnTskId = -1, RscIndex = -1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            double priorty = Rank_b[*lit] + AvrRtSet[*lit];
            if (priorty > max_RANK) {
                pit = lit; max_RANK = priorty;
            }
        }
        CrnTskId = *pit;
        RTI.erase(pit);
        double FinalStartTime = 0, FinalEndTime = InfiniteValue;
        SeletRsc_EFT(chrom, ITL, CrnTskId, RscIndex, FinalStartTime, FinalEndTime);
        chrom.Code_RK[CrnTskId] = RscIndex + Decimals[q];
        chrom.TskSchLst[q] = CrnTskId;
        chrom.RscAlcLst[CrnTskId] = RscIndex;
        chrom.EndTime[CrnTskId] = FinalEndTime;
        MakeSpan = XY_MAX(MakeSpan, FinalEndTime);
        UpdateITL(ITL[RscIndex], FinalStartTime, FinalEndTime);
        for (int l = 0; l < Tasks[CrnTskId].children.size(); ++l) {
            int ChildId = Tasks[CrnTskId].children[l];
            upr[ChildId] = upr[ChildId] -1;
            if (upr[ChildId] == 0){
                RTI.push_back(ChildId);
                AvrRtSet[ChildId] = ClcAvrReadyTime(ChildId,chrom);
            }
        }
    }
    chrom.FitnessValue = MakeSpan;
    return chrom;
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

chromosome GnrChr_Ran() {
    chromosome chrom;
    IntChr(chrom);
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()] + (rand() % 10000) / 10000.0;
    }
    return chrom;
}

chromosome GnrChr_Lvl_Ran() {
    chromosome chrom;
    IntChr(chrom);
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.Code_RK[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()] + (LevelIdOfTask[i] + (rand() % 10000) / 10000.0) / TskLstInLvl.size();
    }
    return chrom;
}

// {in the SS, tasks are arranged according to the level form small to large, and the tasks in the same level are arranged in descend of the number of child tasks}
chromosome GnrChr_HEFT_Baseline() {
    chromosome TemChrom;
    IntChr(TemChrom);
    int ScheduleOrder = 0;
    for (int j = 0; j < TskLstInLvl.size(); ++j) {
        if (TskLstInLvl[j].size() < 2) {
            TemChrom.TskSchLst[ScheduleOrder++]=TskLstInLvl[j][0];
            continue;
        }
        vector<int> SonTaskNum;
        for (int i = 0; i < TskLstInLvl[j].size(); ++i)
            SonTaskNum.push_back(Tasks[TskLstInLvl[j][i]].children.size());

        vector<int> ind(TskLstInLvl[j].size());
        IndexSortByValueOnAscend(ind, SonTaskNum);
        for (int i = TskLstInLvl[j].size() - 1; i >= 0; i--) {
            TemChrom.TskSchLst[ScheduleOrder++] = TskLstInLvl[j][ind[i]];
        }
    }
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

int chooseNextTask(vector<double> rank, vector<int> RT){
    vector<double> Pro(RT.size());
    int sum = 0;
    for(int i = 0; i < RT.size(); ++i){
        sum += rank[RT[i]];
    }
    double randNum = rand() % 100 / 100.0;
    double ProbSum = 0;
    int chosenTask = 0;
    for(int i = 0; i < RT.size(); ++i){
        Pro[i] = rank[RT[i]] / sum;
        ProbSum += Pro[i];
        if (ProbSum > randNum) {
            chosenTask = RT[i];
            break;
        }
    }
    return chosenTask;
}

vector<int> GenerateAChromOrder(vector<double>& rank){
    vector<int > Order;
    vector<int > upr(comConst.NumOfTsk, 0);
    vector<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i] == 0) {
            RTI.push_back(i);
        }
    }
    while (!RTI.empty()){
        int SlcTsk = chooseNextTask(rank, RTI);
        RTI.erase(find(RTI.begin(), RTI.end(), SlcTsk));
        for (int l = 0; l < Tasks[SlcTsk].children.size(); ++l) {
            int childId = Tasks[SlcTsk].children[l];
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                RTI.push_back(childId);
            }
        }
        Order.push_back(SlcTsk);
    }
    return Order;
}

chromosome GnrPrtByRank_Rnd(vector<double>& Rank) {
    chromosome chrom;
    IntChr(chrom);
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        chrom.RscAlcPart[i] =RandomDouble2(0,comConst.NumOfRsc-1); //RandomDouble2(0,comConst.NumOfRsc);//rand() % comConst.NumOfRsc + rand() % 1000 / 1000.0 - 0.5;
    }
    RepairMapAndGnrRscAlcLst(chrom); //GnrRscAlcLst(chrom); //
    chrom.TskSchPart = Rank;
    RepairPriorityAndGnrSchOrd(chrom);
    DcdEvl(chrom, true);
    return chrom;
}

chromosome GnrPrtByRank_EFT(vector<double>& Rank) {
    chromosome chrom;
    IntChr(chrom);
    chrom.TskSchPart = Rank;
    RepairPriorityAndGnrSchOrd(chrom);
    GnrMS_Evl(chrom);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.RscAlcPart[i] = chrom.RscAlcLst[i] - 0.5 + (rand() % 10000) / 10000.0;
    }
    return chrom;
}

void RepairMapAndGnrRscAlcLst(chromosome& ch) {
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        int RscId = round(ch.RscAlcPart[i]);
        if(RscId < Tasks[i].ElgRsc[0]) { //超出下限的处理
            ch.RscAlcPart[i] = Tasks[i].ElgRsc[0];
            ch.RscAlcLst[i] = Tasks[i].ElgRsc[0];
            continue;
        }
        if(RscId > Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1]) { //超出上限的处理
            ch.RscAlcPart[i] = Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1];
            ch.RscAlcLst[i] = Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1];
            continue;
        }
        if(find(Tasks[i].ElgRsc.begin(), Tasks[i].ElgRsc.end(), RscId) == Tasks[i].ElgRsc.end()){ //不存在的处理
            if(Tasks[i].ElgRsc.size() == 1) {
                ch.RscAlcPart[i] = Tasks[i].ElgRsc[0];
                ch.RscAlcLst[i] = Tasks[i].ElgRsc[0];
            } else {
                int TemRscId = FindNearestRscId(i, ch.RscAlcPart[i]);
                ch.RscAlcPart[i] = TemRscId;
                ch.RscAlcLst[i] = TemRscId;
            }
            continue;
        }
        ch.RscAlcLst[i] = RscId;
    }
}

int FindNearestRscId(int TaskId, double value ) {
    for (int j = 0; j < Tasks[TaskId].ElgRsc.size()-1; ++j ){
        if (Tasks[TaskId].ElgRsc[j] < value && value < Tasks[TaskId].ElgRsc[j+1] ) {
            if ( Tasks[TaskId].ElgRsc[j+1] - value < value - Tasks[TaskId].ElgRsc[j] ) {
                return Tasks[TaskId].ElgRsc[j+1];
            } else {
                return Tasks[TaskId].ElgRsc[j];
            }
        }
    }
}

void RepairPriorityAndGnrSchOrd(chromosome& chrom) {
    vector<int> V, Q;
    vector<int> N(comConst.NumOfTsk, -1);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        N[i] = round(chrom.TskSchPart[i]);
    }
    vector<int> upr(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0) {
            Q.push_back(i);
        }
    }
    int MaxV = -1;
    while (V.size() != comConst.NumOfTsk) {
        for (int i = 0; i < Q.size(); ++i) {
            int TaskId = Q[i];
            int MaxP = -1;
            for (int i1 = 0; i1 < Tasks[TaskId].parents.size(); ++i1) {
                if (MaxP < N[Tasks[TaskId].parents[i1]]) {
                    MaxP = N[Tasks[TaskId].parents[i1]];
                }
            }
            if (N[TaskId] <= MaxP) {
                N[TaskId] = MaxP + 1;
            }
            for (int i1 = 0; i1 < V.size(); ++i1) {
                if (N[TaskId] == N[V[i1]])  {
                    N[TaskId] = MaxV + 1;
                    MaxV += 1;
                    break;
                }
            }
            MaxV = XY_MAX(N[TaskId], MaxV);
            V.push_back(TaskId);
        }
        vector<int> TemQ;
        for (int i = 0; i < Q.size(); ++i) {
            int taskId = Q[i];
            for (int i2 = 0; i2 < Tasks[taskId].children.size(); ++i2) {
                int childId = Tasks[taskId].children[i2];
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    TemQ.push_back(childId);
                }
            }
        }
        Q = TemQ;
    }
//    chrom.TskSchPart = N;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.TskSchPart[i] = N[i];
    }
    IndexSortByValueOnAscend(chrom.TskSchLst, N);
}

void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime){
    Parameter_HPSO.InertiaWeight = 0.1 * (1-(runtime / SchTime)) + 0.9;
    Parameter_HPSO.c1 = 2 * (1-(runtime / SchTime));
    Parameter_HPSO.c2 = 2 * (runtime / SchTime);
    double r1 = RandomDouble(0,1);
    double r2 = RandomDouble(0,1);
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        ch.VTskSchPart[i] = Parameter_HPSO.InertiaWeight * ch.VTskSchPart[i] + Parameter_HPSO.c1 * r1 * (Pbest.TskSchPart[i] - ch.TskSchPart[i])
                               + Parameter_HPSO.c2 * r2 * (Gbest.TskSchPart[i] - ch.TskSchPart[i]);
        ch.TskSchPart[i] += ch.VTskSchPart[i];

        ch.VRscAlcPart[i] = Parameter_HPSO.InertiaWeight * ch.VRscAlcPart[i] + Parameter_HPSO.c1 * r1 * (Pbest.RscAlcPart[i] - ch.RscAlcPart[i])
                               + Parameter_HPSO.c2 * r2 * (Gbest.RscAlcPart[i] - ch.RscAlcPart[i]);
        ch.RscAlcPart[i] += ch.VRscAlcPart[i];
    }
    RepairMapAndGnrRscAlcLst(ch); //GnrRscAlcLst(ch); //
    RepairPriorityAndGnrSchOrd(ch);
}
