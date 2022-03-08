#include "classAndVarDefine.h"

vector<Task> Tasks;
vector<vector<int>> TskLstInLvl;
vector<int> LevelIdOfTask;
vector<vector<double>> ParChildTranFileSizeSum;
vector<Resource> Rscs;
vector<chromosome> population;
vector<set<int> > Descendants;
vector<set<int> > Ancestors;
Paramet_ADBRKGA Parameter_ADBRKGA;
ComConst comConst;
double ModelScale;
vector<Orthogonal> orthogonal;
