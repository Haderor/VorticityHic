#include "Experiment_v02.c"
#include <iostream>

int main(int argc, char** argv) {

     TString inFileName, outFileName;

     if (argc < 5){
             std::cerr << "./main -i INPUTFILE -o OUTPUTFILE" << std::endl;
             return 10;
     }

     for (int i=1; i<argc; i++){
             if (std::string(argv[i]) != "-i" &&
                             std::string(argv[i]) != "-o"){
                     std::cerr << "\n[ERROR]: Unknown parameter: " << i << ": " << argv[i] << std::endl;
                     return 10;
             }else{
                     if (std::string(argv[i]) == "-i" && i!=argc-1){
                             inFileName = argv[++i];
                     }
                     if (std::string(argv[i]) == "-i" && i==argc-1){
                             std::cerr << "\n[ERROR]: Input file name was not specified!" << std::endl;
                             return 20;
                     }
                     if (std::string(argv[i]) == "-o" && i!=argc-1){
                             outFileName = argv[++i];
                     }
                     if (std::string(argv[i]) == "-o" && i==argc-1){
                             std::cerr << "\n[ERROR]: Output file name was not specified!" << std::endl;
                             return 21;
                     }
             }
     }

     TFile f(inFileName.Data(), "READ");
     if (!f.IsOpen()) {
        std::cout << "File was not opened!" << std::endl;
	return 0;
     }
     TTree *t = (TTree *)f.Get("data_tree");
     Experiment h(t, outFileName.Data());
     //h.Show(12);
     h.Loop();
     f.Close();

     return 0;
}


