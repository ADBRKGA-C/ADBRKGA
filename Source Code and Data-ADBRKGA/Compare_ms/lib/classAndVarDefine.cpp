#include "classAndVarDefine.h"

vector<Task> Tasks;
vector<Task> OriginalTasks;
vector<int> Oid;
vector<int> Nid;
vector<vector<int>> TskLstInLvl;
vector<int> LevelIdOfTask;
vector<vector<double>> ParChildTranFileSizeSum;
vector<Resource> Rscs;
vector<chromosome> population;
vector<vector<chromosome>> populations;
vector<set<int> > Descendants;
vector<set<int> > Ancestors;
Paramet_ADBRKGA Parameter_ADBRKGA;
Paramet_CGA Parameter_CGA;
Paramet_HGA Parameter_HGA;
Paramet_LWSGA Parameter_LWSGA;
Paramet_HPSO Parameter_HPSO;
Paramet_MOELS Parameter_MOELS;
ComConst comConst;
double ModelScale;
