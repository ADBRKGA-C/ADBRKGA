
#ifndef CSTCHANGE_TOOLS_HPP
#include "common.h"
#define CSTCHANGE_TOOLS_HPP

void CalculateLevelList();
bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b);
bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b);
void IndexSortByValueOnAscend(vector<int>& ind, vector<int>& value);
void IndexSortByValueOnAscend(vector<int>& ind, vector<double>& fitness);
void IndexSortByValueOnDescend(vector<int>& ind, vector<int>& value);
void IndexSortByValueOnDescend(vector<int>& ind, vector<double>& fitness);
double RandomDouble(int start, int end);
double RandomDouble2(int start, int end);
int CalculateParNum(int taskIndex, vector<int>& markList);
int CalculateSonNum(int taskIndex, vector<int>& markList);

#endif //CSTCHANGE_TOOLS_HPP
