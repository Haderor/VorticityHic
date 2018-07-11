#include <iostream>
#include <fstream>
#include <sstream>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"


using namespace std;

int main() {
	ifstream fin("test1.f14");	// File with events
	if (!fin.is_open()) {
		cout << "Attached file was not opened!" << endl;
		exit(15);
	}
	string str;
	stringstream ss;
	TFile *MyFile = new TFile("tree.root", "recreate");	// File to save results
	TTree *t = new TTree("data_tree", "Tree fo save");
	const int count1 = 3, count2 = 13;						// Numbers of lines to skip in test.f14
	const double pi = 3.1415926;
        double psiRp, rx1, ry1, px1, py1;
	const int n = 1000;
	int n_particles, time, pid[n], ityp[n], i3[n], ichg[n], lcl[n], ncl[n], orr[n];
	double r0[n], rx[n], ry[n], rz[n], p0[n], px[n], py[n], pz[n], m[n], imp;		// Data of a particle in an event

	// Data for randomization of angle
        TDatime* ti = new TDatime();
	TRandom3* fi = new TRandom3(ti->GetDate()*ti->GetTime());

	// Create branches of the tree
	t->Branch("n_particles", &n_particles);
	t->Branch("time", &time);
	t->Branch("r0", r0, "");
	t->Branch("rx", rx, "");
	t->Branch("ry", ry, "");
	t->Branch("rz", rz, "");
	t->Branch("p0", p0, "");
	t->Branch("px", px, "");
	t->Branch("py", py, "");
	t->Branch("pz", pz, "");
	t->Branch("m", m, "");
	t->Branch("psiRp", &psiRp);
	t->Branch("pid", pid, "");
	t->Branch("ityp", ityp);
	t->Branch("i3", i3);
	t->Branch("ichg", ichg);
	t->Branch("lcl", lcl);
	t->Branch("ncl", ncl);
	t->Branch("orr", orr);
	t->Branch("imp", &imp);

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
