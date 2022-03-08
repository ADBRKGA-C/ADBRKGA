#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
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

    string strLine;
    ifstream iFile("../exp.txt");
    if (!iFile) {
        cout << "filelist open failed!\n";
        exit(1);
    }
    while(getline(iFile,strLine)){
        istringstream is(strLine);
        Orthogonal TemOrthogonal;
        is >> TemOrthogonal.PopSizeFac >> TemOrthogonal.alpha >> TemOrthogonal.beta
           >> TemOrthogonal.BiasesRate >> TemOrthogonal.ImmigrationRate >> TemOrthogonal.ImprovementRate;
        orthogonal.push_back(TemOrthogonal);
    }
    iFile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string strLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
        }
        getline(iFile, strLine);
        if (strLine.size() < 1) {
            cout << "Empty input file(fileList)" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(strLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        int index =0;
        for(Orthogonal TemOrthogonal: orthogonal){
            ++index;
            ofstream outfile("../ResultOutput/result"+to_string(index)+".txt", ios::app);
            if (!outfile) {
                cout << "Open the file failure...\n";
                exit(0);
            }
            outfile.setf(ios::fixed, ios::floatfield);
            outfile.precision(5);
            cout <<endl<< "Parameter" + to_string(index) << " " << Model << " " << NumOfTask << " " << RscAvlRatio << " ";
            for (int times = 0; times < 10; ++times) {
                double ADBRKGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
                int ADBRKGA_Iteration = 0;
                double ADBRKGA_Result = runADBRKGA(XmlFile, RscAlcFile, TemOrthogonal, ADBRKGA_SchTime, ADBRKGA_Iteration);
                outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                        << ADBRKGA_Result << " " << ADBRKGA_SchTime << " " << ADBRKGA_Iteration
                        << endl;
            }
            outfile.close();
        }
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
    //return 0;
}
