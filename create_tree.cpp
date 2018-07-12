// Reads events from test.f14 of urqmd
#include <iostream>
#include <fstream>
#include <sstream>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"


using namespace std;

int main(int argc, char** argv) {

	TString inFileName, outFileName;

	if (argc < 5){
		std::cerr << "./convert -i INPUTFILE -o OUTPUTFILE" << std::endl;
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

	ifstream fin(inFileName.Data());	// File with events
	if (!fin.is_open()) {
		cout << "Attached file was not opened!" << endl;
		exit(15);
	}

	string str;
	stringstream ss;
	TFile *MyFile = new TFile(outFileName.Data(), "recreate");	// File to save results
	TTree *t = new TTree("data_tree", "Tree fo save");
	const int count1 = 3, count2 = 13;						// Numbers of lines to skip in test.f14
	const double pi = 3.1415926;
        const int n = 2000;

	int n_particles;
	int time;
        double r0[n];
	double rx[n];
	double ry[n];
	double rz[n];
	double p0[n];
	double px[n];
	double py[n];
	double pz[n];
	double m[n];
	double psiRp;
	double imp;
	int pid[n];
	int ityp[n];
	int i3[n];
	int ichg[n];
	int lcl[n];
	int ncl[n];
	int orr[n];

	// Create branches of the tree
	t->Branch("n_particles", &n_particles, "n_particles/I");
	t->Branch("time", &time, "time/D");
	t->Branch("r0", r0, "r0[n_particles]/D");
	t->Branch("rx", rx, "rx[n_particles]/D");
	t->Branch("ry", ry, "ry[n_particles]/D");
	t->Branch("rz", rz, "rz[n_particles]/D");
	t->Branch("p0", p0, "p0[n_particles]/D");
	t->Branch("px", px, "px[n_particles]/D");
	t->Branch("py", py, "py[n_particles]/D");
	t->Branch("pz", pz, "pz[n_particles]/D");
	t->Branch("m", m, "m[n_particles]/D");
	t->Branch("psiRp", &psiRp, "psiRp/D");
	t->Branch("pid", pid, "pid[n_particles]/I");
	t->Branch("ityp", ityp, "ityp[n_particles]/I");
	t->Branch("i3", i3, "i3[n_particles]/I");
	t->Branch("ichg", ichg, "ichg[n_particles]/I");
	t->Branch("lcl", lcl, "lcl[n_particles]/I");
	t->Branch("ncl", ncl, "ncl[n_particles]/I");
	t->Branch("orr", orr, "orr[n_particles]/I");
	t->Branch("imp", &imp, "imp/D");

        // Data for randomization of angle
        TDatime* ti = new TDatime();
        TRandom3* fi = new TRandom3(ti->GetDate()*ti->GetTime());

	// Loop on events
	while(!fin.eof()) {

		// Skip lines
		for (int j = 0; j < count1; j++) {
			getline(fin, str);
		}

		// Read impact parameter
                ss.str("");
                ss.clear();
		getline(fin, str);
		ss << str;
		ss >> str >> imp;
                for (int j = 0; j < count2; j++) {
                        getline(fin, str);
                }


		if(fin.eof()) {
			break;
		}

		// Read number of particles and time
                ss.str("");
                ss.clear();
		getline(fin, str);
		ss << str;
		ss >> n_particles >> time;
		getline(fin, str);

		// Make random reaction plane angle
		psiRp = 2*pi*fi->Rndm();

		// Loop on particles in this event
		for(int j = 0; j < n_particles; j++) {
			ss.str("");
                        ss.clear();
			getline(fin, str);
			ss << str;
			ss  >> r0[j] >> rx[j] >> ry[j] >> rz[j] >> p0[j] >> px[j] >> py[j] >> pz[j] >> m[j] >> ityp[j] >> i3[j] >> ichg[j] >> lcl[j] >> ncl[j] >> orr[j];
			pid[j] = 1000 * (ichg[j] + 2) + ityp[j];
		}
		for(int j = n_particles; j < n; j++) {
			r0[j] = 0.0; rx[j] = 0.0; ry[j] = 0.0; rz[j] = 0.0; p0[j] = 0.0; px[j] = 0.0; py[j] = 0.0; pz[j] = 0.0; m[j] = 0.0; ityp[j] = 0; i3[j] = 0; ichg[j] = 0; lcl[j] = 0; ncl[j] = 0; orr[j] = 0; pid[j] = 0;
		}
		t->Fill();

	}

	t->Write();
	MyFile->Close();

	return 0;
}
