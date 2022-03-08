#include <fstream>
#include <sstream>
#include "config.h"
#include "pugixml.hpp"
#include "common.h"

void DeleteFirstLineInFile(string fileName) {
    vector<string> VecContent;
    string StrLine;
    ifstream iFile(fileName);
    if (!iFile) {
        cout << fileName << "filelist open failed!\n";
        exit(0);
    }
    //Read all content in the document into "VecContent"
    while (iFile) {
        getline(iFile, StrLine);
        VecContent.push_back(StrLine);
    }
    iFile.close();
    VecContent.erase(VecContent.begin()); // delete the first line
    ofstream oFile(fileName);
    if (!oFile) {
        cout << fileName << "filelist open failed!\n";
        exit(0);
    }
    vector<string>::const_iterator iter = VecContent.begin();
    //{Rewrite the contents of "vecContent" into the file}.
    for (; VecContent.end() != iter; ++iter) {
        oFile.write((*iter).c_str(), (*iter).size());
        oFile << '\n';
    }
    oFile.close();
}

int ReadID(string id) {
    int ret = 0;
    for (int i = 0; i < id.length(); ++i) {
        if (id[i] <= '9' && id[i] >= '0')
            ret = ret * 10 + id[i] - '0';
    }
    return ret;
}

//｛Read the model information from the file and store them in the data structures}
void ReadFile(string XmlFile, string RscAlcFile) {
    string FilePath = "../data";
    string XmlPath = FilePath + "/" + XmlFile;
    pugi::xml_document doc;
    int w = doc.load_file((XmlPath).c_str());
    pugi::xml_node root = doc.child("adag");
    for (pugi::xml_node job = root.child("job"); job; job = job.next_sibling("job")) {
        Task task = Task();
        task.length = XY_MAX(fabs(atof(job.attribute("runtime").value())),0.001);// read the length of task
        //{read file}
        for (pugi::xml_node uses = job.child("uses"); uses; uses = uses.next_sibling("uses")) {
            vfile file = vfile();
            file.source = -1;
            if (!strcmp(uses.attribute("link").value(), "input")) {    //read input file
                file.FileName = uses.attribute("file").value();
                file.size = fabs(atof(uses.attribute("size").value()));
                task.IFile.push_back(file);
            } else { //read output file
                file.FileName = uses.attribute("file").value();
                file.size = fabs(atof(uses.attribute("size").value()));
                task.OFile.push_back(file);
            }
        }
        Tasks.push_back(task);
    }
    comConst.NumOfTsk = Tasks.size();

    // read the info of relation between tasks
    for (pugi::xml_node child = root.child("child"); child; child = child.next_sibling("child")) {
        int ChildIndex = ReadID(child.attribute("ref").value());
        for (pugi::xml_node parent = child.child("parent"); parent; parent = parent.next_sibling("parent")) {
            int ParentIndex = ReadID(parent.attribute("ref").value());
            Tasks[ChildIndex].parents.push_back(ParentIndex);
            Tasks[ParentIndex].children.push_back(ChildIndex);
        }
    }

    //{calculate the transfer data size among tasks}
    ParChildTranFileSizeSum.resize(comConst.NumOfTsk);
    for (int k = 0; k < comConst.NumOfTsk; ++k) {
        ParChildTranFileSizeSum[k].resize(comConst.NumOfTsk,0);
    }
    for (int i = 0; i < Tasks.size(); ++i) {
        if (Tasks[i].parents.size() == 0) continue;
        for (int p = 0; p < Tasks[i].IFile.size(); ++p) {               //two loop (for p and j) can be switched in order
            string IName = Tasks[i].IFile[p].FileName;
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {         //Traverse the parent task
                int Parent = Tasks[i].parents[j];
                for (int q = 0; q < Tasks[Parent].OFile.size(); ++q) {  //Traverse the output files of the parent task
                    string OName = Tasks[Parent].OFile[q].FileName;
                    if (IName.compare(OName) == 0) {                    // judge whether two file names are the same; 0: same; -1: not same
                        ParChildTranFileSizeSum[Parent][i] += Tasks[i].IFile[p].size;
                        //If multiple identical files from different parent tasks are transferred to the same child task, the "source" records the last parent task
                        Tasks[i].IFile[p].source = Parent;
                        break;
                    }
                }
            }
        }
    }
    //{Rsc can be added here}
    Resource Rsc_0 = Resource(0, 1*(1-0.24), 20);
    Resource Rsc_1 = Resource(0, 1*(1-0.24), 20);
    Resource Rsc_2 = Resource(0, 2*(1-0.24), 30);
    Resource Rsc_3 = Resource(1, 2*(1-0.24), 30);
    Resource Rsc_4 = Resource(1, 3*(1-0.24), 40);
    Resource Rsc_5 = Resource(1, 3*(1-0.24), 40);
    Rscs.push_back(Rsc_0);
    Rscs.push_back(Rsc_1);
    Rscs.push_back(Rsc_2);
    Rscs.push_back(Rsc_3);
    Rscs.push_back(Rsc_4);
    Rscs.push_back(Rsc_5);
    comConst.NumOfRsc = Rscs.size();
    //read the RscAlc file to task data structure
    //in the RscAlc file, each resource can perform at least one task and each task can be performed by at least one resource
    char line[4096] = {0};
    int TskIndex = 0;
    int RscIndex = 0;
    string RscAlcPath = FilePath + "/" + RscAlcFile;  //RscAlcPath xy4
    ifstream fin(RscAlcPath, ios::in);
    if (!fin) {
        cout << "Error at open Rsc file" << endl;
        exit(0);
    } else {
        while (fin.getline(line, sizeof(line))) {
            stringstream Word(line);
            while (1) {
                Word >> TskIndex;
                if (Word.fail()) break;
                Tasks[TskIndex].ElgRsc.push_back(RscIndex);
                Rscs[RscIndex].ElgTsk.push_back(TskIndex);
            }
            ++RscIndex;
        }
    }
    ModelScale = 0;
    for (int i = 0; i < comConst.NumOfTsk; ++i ){
        ModelScale += Tasks[i].ElgRsc.size();
    }
    fin.close();
}

//task combine
//void TaskCombine() {
//    vector<vector<int>> MT(comConst.NumOfTsk, vector<int>(0));
//    vector<int> flg(comConst.NumOfTsk, 0);
////    vector<int> oid(comConst.NumOfTsk, -1); //记录新的任务编号对应的第一个老任务编号， Oid[δ]表示编号为δ的新任务对应的第一个（合并其它任务的）老任务编号；
////    vector<int> nid(comConst.NumOfTsk, -1); //记录每个老任务编号对应的新任务编号； nid[i']表示编号为i' 的老任务对应的新任务编号；
////    Oid.resize(comConst.NumOfTsk, -1); //记录新的任务编号对应的第一个老任务编号， oid[δ]表示编号为δ的新任务对应的第一个（合并其它任务的）老任务编号；
//    Nid.resize(comConst.NumOfTsk, -1); //记录每个老任务编号对应的新任务编号； Nid[i']表示编号为i' 的老任务对应的新任务编号；
//    vector<Task> TasksAfterTskCombined;
//    vector<int> TaskOrderList;
//    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
//    list<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
//    //{生成一个调度顺序}
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        upr[i] = Tasks[i].parents.size();
//        if (upr[i]==0){
//            RTI.push_back(i);
//        }
//    }
//    while (!RTI.empty()) {
//        int TskId = RTI.front();
//        RTI.erase(RTI.begin());
//        for (int i = 0; i < Tasks[TskId].children.size(); ++i) {
//            --upr[Tasks[TskId].children[i]];
//            if (upr[Tasks[TskId].children[i]] == 0) RTI.push_back(Tasks[TskId].children[i]);
//        }
//        TaskOrderList.push_back(TskId);
//    }
//
//    for (int i = 0; i < comConst.NumOfTsk; i++) {
//        MT[i].push_back(i);                    //每一层的MT中存放的是任务i
//    }
//    //任务合并   从前向后找
//    for (int i : TaskOrderList) {
//        if (flg[i] == 0) {
//            while (Tasks[i].children.size() == 1 && Tasks[Tasks[i].children[0]].parents.size() == 1 && Tasks[i].ElgRsc == Tasks[Tasks[i].children[0]].ElgRsc) {
//                int icl = Tasks[i].children[0];            //记录被合并的任务i的子任务
//                MT[i].push_back(icl);                      //向MT集合中添加可合并任务i的子任务
//                MT[icl].clear();                           //将被合并的任务从所属MT中删除
//                flg[icl] = 1;
//                Tasks[i].length += Tasks[icl].length;      //length = length i+
//                Tasks[i].children = Tasks[icl].children;   //SCi = SCi+ 将子任务中集合的子任务全部插入到父任务集合中去
//                Tasks[i].OFile = Tasks[icl].OFile;         //OFLi = OFLi+ 将子任务中集合的输出文件全部插入到父任务集合中去
//                for (int k = 0; k < Tasks[icl].children.size(); k++) {  //PRi++ = PRi++ - ti+ + ti//遍历每一个子任务的子任务
//                    replace(Tasks[Tasks[icl].children[k]].parents.begin(),Tasks[Tasks[icl].children[k]].parents.end(), icl, i);
//                }
//            }
//        }
//    }
//    //重新编号
//    int N = 0;
//    for (int i = 0; i < comConst.NumOfTsk; i++) {
//        if (!MT[i].empty()) {
//            Oid.push_back(i);   //Oid[N] = i;//记录新的任务编号对应的第一个老任务编号， Oid[n]表示编号为n的新任务对应的第一个（合并其它任务的）老任务编号；
//            for (int j : MT[i]) {
//                Nid[j] = N;                 //记录每个老任务编号对应的新任务编号； Nid[i']表示编号为i'的老任务对应的新任务编号；
//            }
//            N = N + 1;
//        }
//    }
//    //根据新编号更新任务信息
//    for (int n = 0; n < N; n++) {
//        Task task = Task();                     //生成一个新的任务对象
//        for (int i : Tasks[Oid[n]].parents) {
//            task.parents.push_back(Nid[i]);
//        }
//
//        for (int i : Tasks[Oid[n]].children) {
//            task.children.push_back(Nid[i]);
//        }
//        //把新的参数赋值给新生成的对象
//        task.IFile = Tasks[Oid[n]].IFile;
//        task.OFile = Tasks[Oid[n]].OFile;
//        task.length = Tasks[Oid[n]].length;
//        task.ElgRsc = Tasks[Oid[n]].ElgRsc;
//        TasksAfterTskCombined.push_back(task);
//    }
//    //更新任务数量
//    comConst.NumOfTsk = TasksAfterTskCombined.size();
//    //更新父子任务之间传输文件大小
//    ParChildTranFileSizeSum.resize(comConst.NumOfTsk); //建议要保存合并前父子任务之间传输文件大小-xy-2022.01.09
//    for (int k = 0; k < comConst.NumOfTsk; ++k) {
//        ParChildTranFileSizeSum[k].resize(comConst.NumOfTsk,0);
//    }
//    for (int i = 0; i < TasksAfterTskCombined.size(); ++i) {
//        if (TasksAfterTskCombined[i].parents.size() == 0) continue;
//        for (int p = 0; p < TasksAfterTskCombined[i].IFile.size(); ++p) {
//            string iName = TasksAfterTskCombined[i].IFile[p].FileName;
//            int flag = 0;
//            for (int j = 0; j < TasksAfterTskCombined[i].parents.size(); ++j) {
//                if (flag == 1) break;
//                int parent = TasksAfterTskCombined[i].parents[j];
//                for (int q = 0; q < TasksAfterTskCombined[parent].OFile.size(); ++q) {
//                    string oName = TasksAfterTskCombined[parent].OFile[q].FileName;
//                    //判断文件名相同
//                    if (iName.compare(oName) == 0) {
//                        ParChildTranFileSizeSum[parent][i] += TasksAfterTskCombined[i].IFile[p].size;
//                        TasksAfterTskCombined[i].IFile[p].source = parent;
//                        flag = 1;
//                        break;
//                    }
//                }
//            }
//        }
//    }
//    OriginalTasks = Tasks;              //保存合并前的任务关系
//    Tasks = TasksAfterTskCombined;      //更新合并后的任务关系
//}

void ClearALL() {
    Tasks.clear();
//    OriginalTasks.clear();
//    Oid.clear();
//    Nid.clear();
    TskLstInLvl.clear();
    LevelIdOfTask.clear();
    ParChildTranFileSizeSum.clear();
    Rscs.clear();
    population.clear();
    populations.clear();
    Descendants.clear();
    Ancestors.clear();
}

void ConfigParameter_CGA() {
    Parameter_CGA.NumOfChromPerPop = 2 * Tasks.size();
    Parameter_CGA.MutationRate = 0.2;
    Parameter_CGA.CrossoverRate = 0.8;
}

void ConfigParameter_HGA() {
    Parameter_HGA.NumOfChromPerPop = 100;
    Parameter_HGA.MutationRate = 0.02;
    Parameter_HGA.EliteRate = 0.2;
}

void ConfigParameter_LWSGA() {
    Parameter_LWSGA.NumOfChromPerPop = 70;
    Parameter_LWSGA.CrossoverRate = 0.7;
}

void ConfigParameter_MOELS() {
    Parameter_MOELS.NumOfChromPerPop = 50;
    Parameter_MOELS.CrossoverRate = 0.7;
    Parameter_MOELS.MutationRate = 0.7;
}

void ConfigParameter_HPSO() {
    Parameter_HPSO.NumOfChromPerPop = 20;
    Parameter_HPSO.InertiaWeight = 1;
    Parameter_HPSO.c1 = 2;
    Parameter_HPSO.c2 = 2;
}

void ConfigParameter_ADBRKGA() {
    Parameter_ADBRKGA.NumOfChromPerPop = Tasks.size()*1.5;        //set the subpopulation size
    if (Parameter_ADBRKGA.NumOfChromPerPop % 2 == 1) {
        ++Parameter_ADBRKGA.NumOfChromPerPop;
    }
    Parameter_ADBRKGA.ImmigrationRate = 0.07;
    Parameter_ADBRKGA.alpha = 2;
    Parameter_ADBRKGA.beta = 0.7;
    Parameter_ADBRKGA.ImprovementRate= 0.05;
    Parameter_ADBRKGA.BiasesRate= 0.75;
}
