#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "HEFT.h"
#include "IHEFT3.h"
#include "DHEFT.h"
#include "HGA.h"
#include "MOELS.h"
#include "LWSGA.h"
#include "CGA.h"
#include "HPSO.h"
#include "ADBRKGA.h"

using namespace std;

int main() {
    srand((int) time(0));
    //{set the fitness (makespan) when algorithm is terminated }
    map<string, double> EndFit;
    EndFit["Montage25_0.4"] = 53.93;
    EndFit["Montage25_0.7"] = 36.86;
    EndFit["Montage25_1.0"] = 36.79;
    EndFit["Montage50_0.4"] = 74.85;
    EndFit["Montage50_0.7"] = 69.69;
    EndFit["Montage50_1.0"] = 69.40;
    EndFit["Montage100_0.4"] =178.08;
    EndFit["Montage100_0.7"] =142.64;
    EndFit["Montage100_1.0"] =135.52;

    EndFit["CyberShake30_0.4"] = 497.77;
    EndFit["CyberShake30_0.7"] = 200.83;
    EndFit["CyberShake30_1.0"] = 200.83;
    EndFit["CyberShake50_0.4"] = 560.90;
    EndFit["CyberShake50_0.7"] = 237.20;
    EndFit["CyberShake50_1.0"] = 223.54;
    EndFit["CyberShake100_0.4"] =440.11;
    EndFit["CyberShake100_0.7"] =312.92;
    EndFit["CyberShake100_1.0"] =303.40;

    EndFit["Epigenomics24_0.4"] = 5753.22;
    EndFit["Epigenomics24_0.7"] = 3059.12;
    EndFit["Epigenomics24_1.0"] = 3052.24;
    EndFit["Epigenomics47_0.4"] = 6378.32;
    EndFit["Epigenomics47_0.7"] = 5714.36;
    EndFit["Epigenomics47_1.0"] = 5674.08;
    EndFit["Epigenomics100_0.4"] =60778.14;
    EndFit["Epigenomics100_0.7"] =46797.51;
    EndFit["Epigenomics100_1.0"] =46384.63;

    EndFit["Ligo30_0.4"] = 1196.74;
    EndFit["Ligo30_0.7"] = 888.56;
    EndFit["Ligo30_1.0"] = 862.44;
    EndFit["Ligo50_0.4"] =1921.77;
    EndFit["Ligo50_0.7"] =1333.52;
    EndFit["Ligo50_1.0"] =1309.08;
    EndFit["Ligo100_0.4"] =3780.67;
    EndFit["Ligo100_0.7"] =2323.59;
    EndFit["Ligo100_1.0"] =2311.38;

    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();
    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
        }
        getline(iFile, StrLine);
        if (StrLine.size() < 1) {
            cout << "Empty input file" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(StrLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";

        TrmFactor = 300;

        double EndFitness = EndFit[Model + NumOfTask + "_" + RscAvlRatio];

        double HGA_SchTime  = 0;
        int HGA_Iteration = 0;
        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration, EndFitness);
        ClearALL();
        double MOELS_SchTime  = 0;
        int MOELS_Iteration = 0;
        double MOELS_Result = runMOELS(XmlFile, RscAlcFile, MOELS_SchTime, MOELS_Iteration, EndFitness);
        ClearALL();
        double LWSGA_SchTime  = 0;
        int LWSGA_Iteration = 0;
        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration, EndFitness);
        ClearALL();
        double CGA_SchTime  = 0;
        int CGA_Iteration = 0;
        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration, EndFitness);
        ClearALL();
        double HPSO_SchTime  = 0;
        int HPSO_Iteration = 0;
        double HPSO_Result = runHPSO(XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration, EndFitness);
        ClearALL();
        double ADBRKGA_SchTime  = 0;
        int ADBRKGA_Iteration = 0;
        double ADBRKGA_Result = runADBRKGA(XmlFile, RscAlcFile, ADBRKGA_SchTime, ADBRKGA_Iteration, EndFitness);
        ClearALL();

        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "
                << MOELS_Result << " " << MOELS_SchTime << " " << MOELS_Iteration << " "
                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " "
                << HPSO_Result << " " << HPSO_SchTime << " " << HPSO_Iteration << " "
                << ADBRKGA_Result << " " << ADBRKGA_SchTime << " " << ADBRKGA_Iteration << " "
                << endl;
        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}
