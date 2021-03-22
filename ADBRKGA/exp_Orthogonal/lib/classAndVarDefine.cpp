#include "classAndVarDefine.h"


vector<Task> Tasks;
vector<vector<int>> TskLstInLvl;
vector<vector<double>> ParChildTranFileSizeSum;
vector<int> LevelIdOfTask;
vector<Resource> Rscs;
vector<chromosome> population;
vector<vector<chromosome>> populations;
Paramet_ADRKGA Parameter_ADRKGA;
ComConst comConst;
double ModelScale;
vector<set<int> > Descendants;
vector<set<int> > Ancestors;
vector<Orthogonal> orthogonal;