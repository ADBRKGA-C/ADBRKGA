
#ifndef CSTCHANGE_CONFIG_H
#include "common.h"
#define CSTCHANGE_CONFIG_H
void DeleteFirstLineInFile(string fileName);
int ReadID(string id);
void ReadFile(string XmlFile,string RscAlcFile);
void ClearALL();
void ConfigParameter_CGA();
void ConfigParameter_HGA();
void ConfigParameter_LWSGA();
void ConfigParameter_NGA();
void ConfigParameter_ADRKGA(Orthogonal orthogonal);
#endif //CSTCHANGE_CONFIG_H
