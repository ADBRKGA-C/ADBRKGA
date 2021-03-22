
#ifndef CSTCHANGE_TOOLS_HPP
#include "common.h"
#define CSTCHANGE_TOOLS_HPP

void CalculateLevelList();
bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b);
bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b);
void IndexSort(vector<int>& ind, vector<int>& value);
void IndexSort(vector<int>& ind, vector<double>& fitness);
double RandomDouble(int start, int end);
int CalculateParNum(int taskIndex, vector<int>& markList);
int CalculateSonNum(int taskIndex, vector<int>& markList);
void AdpDcd (chromosome&chrom, double& Time, double& SchTime);
#endif //CSTCHANGE_TOOLS_HPP
