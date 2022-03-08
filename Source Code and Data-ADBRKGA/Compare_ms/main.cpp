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
    //{set the runtime (termination) time of the algorithm} -xy4
    map<string, double> SchTime;
    SchTime["Montage25_0.4"] = 0.571;
    SchTime["Montage25_0.7"] = 0.868;
    SchTime["Montage25_1.0"] = 1.217;
    SchTime["Montage50_0.4"] = 1.430;
    SchTime["Montage50_0.7"] = 2.370;
    SchTime["Montage50_1.0"] = 3.196;
    SchTime["Montage100_0.4"] = 2.037;
    SchTime["Montage100_0.7"] = 5.764;
    SchTime["Montage100_1.0"] = 9.068;

    SchTime["CyberShake30_0.4"] = 0.614;
    SchTime["CyberShake30_0.7"] = 1.104;
    SchTime["CyberShake30_1.0"] = 1.603;
    SchTime["CyberShake50_0.4"] = 0.965;
    SchTime["CyberShake50_0.7"] = 2.517;
    SchTime["CyberShake50_1.0"] = 3.294;
    SchTime["CyberShake100_0.4"] = 2.012;
    SchTime["CyberShake100_0.7"] = 6.863;
    SchTime["CyberShake100_1.0"] = 9.766;

    SchTime["Epigenomics24_0.4"] = 0.532;
    SchTime["Epigenomics24_0.7"] = 0.888;
    SchTime["Epigenomics24_1.0"] = 1.120;
    SchTime["Epigenomics47_0.4"] = 1.417;
    SchTime["Epigenomics47_0.7"] = 1.956;
    SchTime["Epigenomics47_1.0"] = 2.703;
    SchTime["Epigenomics100_0.4"] = 3.004;
    SchTime["Epigenomics100_0.7"] = 7.400;
    SchTime["Epigenomics100_1.0"] = 13.453;

    SchTime["Ligo30_0.4"] = 0.607;
    SchTime["Ligo30_0.7"] = 0.856;
    SchTime["Ligo30_1.0"] = 1.606;
    SchTime["Ligo50_0.4"] = 1.347;
    SchTime["Ligo50_0.7"] = 2.044;
    SchTime["Ligo50_1.0"] = 3.023;
    SchTime["Ligo100_0.4"] = 2.956;
    SchTime["Ligo100_0.7"] = 5.632;
    SchTime["Ligo100_1.0"] = 7.747;

    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(0);
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
        is >> Model >> NumOfTask >> RscAvlRatio;  //NumOfTask, RscAvlRatio
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";

        double HEFT_SchTime  = 0 ;
        double HEFT_Result = runHEFT(XmlFile, RscAlcFile, HEFT_SchTime);

        double IHEFT3_SchTime  = 0 ;
        double IHEFT3_Result = runIHEFT3(XmlFile, RscAlcFile,IHEFT3_SchTime);

        double DHEFT_SchTime  = 0 ;
        double DHEFT_Result = runDHEFT(XmlFile, RscAlcFile,DHEFT_SchTime);

        double HGA_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HGA_Iteration = 0;
        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);

        double MOELS_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int MOELS_Iteration = 0;
        double MOELS_Result = runMOELS(XmlFile, RscAlcFile,MOELS_SchTime, MOELS_Iteration);

        double LWSGA_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int LWSGA_Iteration = 0;
        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);

        double CGA_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int CGA_Iteration = 0;
        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration);

        double HPSO_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HPSO_Iteration = 0;
        double HPSO_Result = runHPSO(XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration);

        double ADBRKGA_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int ADBRKGA_Iteration = 0;
        double ADBRKGA_Result = runADBRKGA(XmlFile, RscAlcFile, ADBRKGA_SchTime, ADBRKGA_Iteration);

        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << HEFT_Result << " " << HEFT_SchTime << " "
                << IHEFT3_Result << " " << IHEFT3_SchTime << " "
                << DHEFT_Result << " " << DHEFT_SchTime << " "
                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "
                << MOELS_Result << " " << MOELS_SchTime << " " << MOELS_Iteration << " "
                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
                << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " "
                << HPSO_Result << " " << HPSO_SchTime << " " << HPSO_Iteration << " "
                << ADBRKGA_Result  << " " << ADBRKGA_SchTime  << " " << ADBRKGA_Iteration  << " "
                << endl;
        outfile.close();
        DeleteFirstLineInFile("../fileList.txt"); //delete the first line in the file
    } while (1);
    //return 0;
}
