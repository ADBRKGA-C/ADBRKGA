#include <cstdlib>
#include <GenerateAChrom.h>
#include <unordered_set>
#include "tools.hpp"
#include "common.h"

using namespace std;

chromosome Crs_BPUC(chromosome& chrom1,chromosome& chrom2){
    chromosome temCh;
    IntChr(temCh);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double alpha = double(rand() % 100) / 100;
        if (alpha < Parameter_ADBRKGA.BiasesRate) {
            temCh.Code_RK[i] = chrom1.Code_RK[i];
        } else {
            temCh.Code_RK[i] = chrom2.Code_RK[i];
        }
    }
    return temCh;
}
//{calculate the cumulative probabilities for the population whose chromosome have been sorted}
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A , int& NumOfChromPerPop) {
    for (int n = 0; n < NumOfChromPerPop; ++n) {
        A[n] = pow(RtOfSltPrb,NumOfChromPerPop-1-n) * (RtOfSltPrb - 1) / (pow(RtOfSltPrb, NumOfChromPerPop) - 1);
    }
    for (int n = 1; n < NumOfChromPerPop; ++n){
        A[n] = A[n] + A[n - 1];
    }
}

//{select a chromosome using roulette wheel selection scheme}
int SelectChrom(vector<double>& A) {
    double lambda = RandomDouble(0, 1);
    for (int n = 0; n < A.size(); ++n)  // -xy
        if (lambda + PrecisionValue <= A[n])
            return n;
}


void AdpDcd (chromosome&chrom ,double& CurTime,double& TotalTime) {
    vector<double> SP(3);
    SP[0] = pow(CurTime/TotalTime,Parameter_ADBRKGA.alpha);
    SP[1] = Parameter_ADBRKGA.beta * (1-pow(CurTime/TotalTime,Parameter_ADBRKGA.alpha));
    SP[2] = (1-Parameter_ADBRKGA.beta) * (1-pow(CurTime/TotalTime,Parameter_ADBRKGA.alpha));
    double RandNum = double (rand()%1000) / 1000;
    if (RandNum < SP[0]){
        NrmDcd(chrom, true);
    }
    if (RandNum >= SP[0] && RandNum < (SP[0]+SP[1])){
        HrsDcd_EFT(chrom);
    }
    if (RandNum >= (SP[0]+SP[1])){
        HrsDcd_CTP(chrom);
    }
}

