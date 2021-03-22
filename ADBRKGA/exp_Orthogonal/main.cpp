#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "ADBRKGA.h"

using namespace std;

int main() {
    srand((int) time(0));
    //{set the runtime (termination) time of the algorithm}
    map<string, double> SchTime;
    SchTime["Montage25_0.4"] = 0.585 ;
    SchTime["Montage25_0.7"] = 0.873 ;
    SchTime["Montage25_1.0"] = 1.131 ;
    SchTime["Montage50_0.4"] = 1.920 ;
    SchTime["Montage50_0.7"] = 3.750 ;
    SchTime["Montage50_1.0"] = 6.074  ;
    SchTime["Montage100_0.4"] =7.024  ;
    SchTime["Montage100_0.7"] =16.787 ;
    SchTime["Montage100_1.0"] =29.335 ;

    SchTime["CyberShake30_0.4"] = 0.544 ;
    SchTime["CyberShake30_0.7"] = 1.421 ;
    SchTime["CyberShake30_1.0"] = 1.769 ;
    SchTime["CyberShake50_0.4"] = 1.200 ;
    SchTime["CyberShake50_0.7"] = 2.889 ;
    SchTime["CyberShake50_1.0"] = 4.591  ;
    SchTime["CyberShake100_0.4"] =5.587  ;
    SchTime["CyberShake100_0.7"] =11.647 ;
    SchTime["CyberShake100_1.0"] =30.125 ;

    SchTime["Epigenomics24_0.4"] = 0.581 ;
    SchTime["Epigenomics24_0.7"] = 0.734 ;
    SchTime["Epigenomics24_1.0"] = 1.274 ;
    SchTime["Epigenomics47_0.4"] = 1.991 ;
    SchTime["Epigenomics47_0.7"] = 2.861 ;
    SchTime["Epigenomics47_1.0"] = 3.511 ;
    SchTime["Epigenomics100_0.4"] =4.886 ;
    SchTime["Epigenomics100_0.7"] =12.696  ;
    SchTime["Epigenomics100_1.0"] =28.482  ;

    SchTime["Ligo30_0.4"] = 0.730 ;
    SchTime["Ligo30_0.7"] = 1.135 ;
    SchTime["Ligo30_1.0"] = 2.066 ;
    SchTime["Ligo50_0.4"] = 1.586 ;
    SchTime["Ligo50_0.7"] = 4.541 ;
    SchTime["Ligo50_1.0"] = 8.385  ;
    SchTime["Ligo100_0.4"] =5.038  ;
    SchTime["Ligo100_0.7"] =9.368  ;
    SchTime["Ligo100_1.0"] =15.334 ;

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
        int index =0 ;
        for(Orthogonal TemOrthogonal: orthogonal){
            ++index;
            ofstream outfile("../print/result"+to_string(index)+".txt", ios::app);
            if (!outfile) {
                cout << "Open the file failure...\n";
                exit(0);
            }
            outfile.setf(ios::fixed, ios::floatfield);
            outfile.precision(5);
            for (int times = 0; times < 10; ++times) {
                cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";
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
    return 0;
}
