#define Experiment_cxx
#include "Experiment.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

void Experiment::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Experiment.C
//      Root > Experiment t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	TFile *MyFile = new TFile(fname, "recreate");	// File to save results

	// Constants for histograms
        const Double_t binsY[] = {-1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};
        const Double_t binsPt[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0};
        const Int_t NbinsY = 24;
        const Int_t NbinsPt = 11;
	const string cuts = "-1.2 < Eta < 1.2, 0.2 < Pt < 3";
	const Int_t bins_x = 20, bins_y = 20, bins_z = 20;

        Double_t rx1, ry1, px1, py1, pt, fi;    //Temporary storage for coordinates and momentum
        Double_t X, Y, Z;               // Limits for coordinates
        X = Y = Z = 200.0;
        Double_t dx, dy, dz;    // Splitting characteristic
        dx = dy = 20.0;
        dz = dx;
        int nx, ny, nz;         // Number of segmentss
        nx = 2*(X/dx) + 1;
        ny = 2*(Y/dy) + 1;
        nz = 2*(Z/dz) + 1;


	// Histograms
	map <Int_t, TH1D*> arrHistPt;
        map <Int_t, TH1D*> arrHistRap;
	// Flow coeffitients
	map <Int_t, TProfile*> arrProfv1PtSmall;
	map <Int_t, TProfile*> arrProfv2PtSmall;
        map <Int_t, TProfile*> arrProfv1RapSmall;
        map <Int_t, TProfile*> arrProfv2RapSmall;
        map <Int_t, TProfile*> arrProfv1PtMid;
        map <Int_t, TProfile*> arrProfv2PtMid;
        map <Int_t, TProfile*> arrProfv1RapMid;
        map <Int_t, TProfile*> arrProfv2RapMid;
        map <Int_t, TProfile*> arrProfv1PtBig;
        map <Int_t, TProfile*> arrProfv2PtBig;
        map <Int_t, TProfile*> arrProfv1RapBig;
        map <Int_t, TProfile*> arrProfv2RapBig;
	// Dependence v(Ro)
        map <Int_t, TProfile*> arrProfvRo;
	// Momentum, energy and velocity
	map <Int_t, TH3D*> arrHistPxSmall;
        map <Int_t, TH3D*> arrHistPySmall;
        map <Int_t, TH3D*> arrHistPzSmall;
        map <Int_t, TH3D*> arrHistESmall;
        map <Int_t, TH3D*> arrHistVx_2Small;	// Sum of velocities in the cell
        map <Int_t, TH3D*> arrHistVy_2Small;
        map <Int_t, TH3D*> arrHistVz_2Small;

	// Names and PIDs of particles
        vector<string> arrNames = {"All", "Pion", "Kaon", "Proton", "Lambda"};
	vector<Int_t> arrPIDs = {0, 3101, 3106, 3001, 2027};

	if(arrPIDs.size() != arrNames.size()) {
		cerr << "Error in comparison of particle name and PID" << endl;
		return;
	}

	MyFile->mkdir("Small parameter");
        MyFile->mkdir("Mid parameter");
	MyFile->mkdir("Big parameter");

	for (Int_t i = 0; i < arrPIDs.size(); i++) {
			arrHistPt[arrPIDs[i]] = new TH1D(("Pt" + arrNames[i]).c_str(), cuts.c_str(), 500, 0, 5);
			arrHistRap[arrPIDs[i]] =  new TH1D(("Rap" + arrNames[i]).c_str(), cuts.c_str(), 200, -5, 5);

			// Velocity from ro
			arrProfvRo[arrPIDs[i]] = new TProfile(("vRo" + arrNames[i]).c_str(), cuts.c_str(), 500, 0, sqrt(nx*nx/4.0 + ny*ny/4.0));
        		arrProfvRo[arrPIDs[i]]->GetXaxis()->SetTitle("r, cells");

			// v1 and v2 coefficients
			MyFile->cd("Small parameter");
			arrProfv1PtSmall[arrPIDs[i]] = new TProfile(("v1PtSmall" + arrNames[i]).c_str(), "imp 0.0-3.0;pt, GeV/c;v1", NbinsPt, binsPt);
			arrProfv2PtSmall[arrPIDs[i]] = new TProfile(("v2PtSmall" + arrNames[i]).c_str(), "imp 0.0-3.0;pt, GeV/c;v2", NbinsPt, binsPt);
                        arrProfv1RapSmall[arrPIDs[i]] = new TProfile(("v1RapSmall" + arrNames[i]).c_str(), "imp 0.0-3.0;Rap;v1", NbinsY, binsY);
                        arrProfv2RapSmall[arrPIDs[i]] = new TProfile(("v2RapSmall" + arrNames[i]).c_str(), "imp 0.0-3.0;Rap;v2", NbinsY, binsY);
			MyFile->cd("../");

        		MyFile->cd("Mid parameter");
			arrProfv1PtMid[arrPIDs[i]] = new TProfile(("v1PtMid" + arrNames[i]).c_str(), "imp 6.0-7.0;pt, GeV/c;v1", NbinsPt, binsPt);
                        arrProfv2PtMid[arrPIDs[i]] = new TProfile(("v2PtMid" + arrNames[i]).c_str(), "imp 6.0-7.0;pt, GeV/c;v2", NbinsPt, binsPt);
                        arrProfv1RapMid[arrPIDs[i]] = new TProfile(("v1RapMid" + arrNames[i]).c_str(), "imp 6.0-7.0;Rap;v1", NbinsY, binsY);
                        arrProfv2RapMid[arrPIDs[i]] = new TProfile(("v2RapMid" + arrNames[i]).c_str(), "imp 6.0-7.0;Rap;v2", NbinsY, binsY);
                        MyFile->cd("../");

        		MyFile->cd("Big parameter");
			arrProfv1PtBig[arrPIDs[i]] = new TProfile(("v1PtBig" + arrNames[i]).c_str(), "imp 10.0-12.0;pt, GeV/c;v1", NbinsPt, binsPt);
                        arrProfv2PtBig[arrPIDs[i]] = new TProfile(("v2PtBig" + arrNames[i]).c_str(), "imp 10.0-12.0;pt, GeV/c;v2", NbinsPt, binsPt);
                        arrProfv1RapBig[arrPIDs[i]] = new TProfile(("v1RapBig" + arrNames[i]).c_str(), "imp 10.0-12.0;Rap;v1"", NbinsY, binsY);
                        arrProfv2RapBig[arrPIDs[i]] = new TProfile(("v2RapBig" + arrNames[i]).c_str(), "imp 10.0-12.0;Rap;v2", NbinsY, binsY);
                        MyFile->cd("../");

			// Momentum, Energy and velocity
                        arrHistPxSmall[arrPIDs[i]] = new TH3D(("PxSmall" + arrNames[i]).c_str(), (cuts + "; x, fm; y, fm; z, fm").c_str(), bins_x, -X, X, bins_y, -Y, Y, bins_z, -Z, Z);
                        arrHistPySmall[arrPIDs[i]] = new TH3D(("PySmall" + arrNames[i]).c_str(), (cuts + "; x, fm; y, fm; z, fm").c_str(),bins_x, -X, X, bins_y, -Y, Y, bins_z, -Z, Z);
                        arrHistPzSmall[arrPIDs[i]] = new TH3D(("PzSmall" + arrNames[i]).c_str(), (cuts + "; x, fm; y, fm; z, fm").c_str(),bins_x, -X, X, bins_y, -Y, Y, bins_z, -Z, Z);
                        arrHistVx_2Small[arrPIDs[i]] = new TH3D(("Vx_2Small" + arrNames[i]).c_str(), (cuts + "; x, fm; y, fm; z, fm").c_str(),bins_x, -X, X, bins_y, -Y, Y, bins_z, -Z, Z);
                        arrHistVy_2Small[arrPIDs[i]] = new TH3D(("Vy_2Small" + arrNames[i]).c_str(), (cuts + "; x, fm; y, fm; z, fm").c_str(),bins_x, -X, X, bins_y, -Y, Y, bins_z, -Z, Z);
                        arrHistVz_2Small[arrPIDs[i]] = new TH3D(("Vz_2Small" + arrNames[i]).c_str(), (cuts + "; x, fm; y, fm; z, fm").c_str(),bins_x, -X, X, bins_y, -Y, Y, bins_z, -Z, Z);
                        arrHistESmall[arrPIDs[i]] = new TH3D(("ESmall" + arrNames[i]).c_str(), (cuts + "; x, fm; y, fm; z, fm").c_str(),bins_x, -X, X, bins_y, -Y, Y, bins_z, -Z, Z);

	}
	TH1D *HistFi = new TH1D("HistFi", "", 50, -7, 7);
        TH1D *HistFiPsi = new TH1D("HistFiPsi", "", 50, -7, 7);

//	TProfile *vRo = new TProfile("vRo", "", 500, 0, sqrt(nx*nx/4.0 + ny*ny/4.0));       // Velocity from ro
//	vRo->GetXaxis()->SetTitle("r, cells");

	Long64_t nbytes = 0; // nb = 0;
	double t0 = 200.0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);

		// Take only events at time t0
		if (time != t0) {
			continue;
		}

		for (Int_t i = 0; i<n_particles; i++) {
			// Rotate RP on angle psiRp
			rx1=rx[i]; ry1=ry[i];
			rx[i] = rx1*cos(-psiRp)+ry1*sin(-psiRp);
			ry[i] = -rx1*sin(-psiRp)+ry1*cos(-psiRp);
			px1=px[i]; py1=py[i];
            		px[i] = px1*cos(-psiRp)+py1*sin(-psiRp);
            		py[i] = -px1*sin(-psiRp)+py1*cos(-psiRp);
			fi = TMath::ATan2(py[i], px[i]);

			TVector3 p(px[i], py[i], pz[i]);
                        TLorentzVector pL(p, p0[i]);

			// Cuts for particles
			if ((p.Pt() < 0.2) || (3.0 < p.Pt()) || (pL.Eta() < -1.2) || (1.2 < pL.Eta())) {
				continue;
			}
                        HistFi->Fill(TMath::ACos(TMath::Cos(fi)));
                        HistFiPsi->Fill(TMath::ACos(TMath::Cos(fi - psiRp)));

			if((0.0 < imp) && (imp < 3.0)) {
                        	arrHistPxSmall.at(0)->Fill(rx[i], ry[i], rz[i], px[i]);
                        	arrHistPySmall.at(0)->Fill(rx[i], ry[i], rz[i], py[i]);
                        	arrHistPzSmall.at(0)->Fill(rx[i], ry[i], rz[i], pz[i]);
                        	arrHistESmall.at(0)->Fill(rx[i], ry[i], rz[i], p0[i]);
				if (p0[i] != 0) {
                        		arrHistVx_2Small.at(0)->Fill(rx[i], ry[i], rz[i], px[i]/p0[i]);
                        		arrHistVy_2Small.at(0)->Fill(rx[i], ry[i], rz[i], py[i]/p0[i]);
                        		arrHistVz_2Small.at(0)->Fill(rx[i], ry[i], rz[i], pz[i]/p0[i]);
				}
			}

			// Draw histogramms
			if (arrHistPt.count(pid[i])) {
				// Draw histogramms for different types of particles
				arrHistPt.at(pid[i])->Fill(p.Pt());
	                        arrHistRap.at(pid[i])->Fill(pL.Rapidity());
				if ((0.0 < imp) && (imp < 3.0)) {
	                        	arrHistPxSmall.at(pid[i])->Fill(rx[i], ry[i], rz[i], px[i]);
	                        	arrHistPySmall.at(pid[i])->Fill(rx[i], ry[i], rz[i], py[i]);
        	                	arrHistPzSmall.at(pid[i])->Fill(rx[i], ry[i], rz[i], pz[i]);
                	        	arrHistESmall.at(pid[i])->Fill(rx[i], ry[i], rz[i], p0[i]);
					if (p0[i] != 0) {
                                		arrHistVx_2Small.at(pid[i])->Fill(rx[i], ry[i], rz[i], px[i]/p0[i]);
                                		arrHistVy_2Small.at(pid[i])->Fill(rx[i], ry[i], rz[i], py[i]/p0[i]);
                                		arrHistVz_2Small.at(pid[i])->Fill(rx[i], ry[i], rz[i], pz[i]/p0[i]);
					}
				}

				// Draw v1 and v2 for different impact parameters and particles
				if ((0.2 < p.Pt()) && (p.Pt() < 3.0)) {
                                        if (imp < 3.0) {
                                               arrProfv2RapSmall.at(pid[i])->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((6.0 < imp) && (imp < 7.0)) {
                                                arrProfv2RapMid.at(pid[i])->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((10.0 < imp) && (imp < 12.0)) {
                                                arrProfv2RapBig.at(pid[i])->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                }
                                if((-1.0 < pL.Rapidity()) && (pL.Rapidity() < 1.0)) {
                                        if (imp < 3.0) {
                                                arrProfv1RapSmall.at(pid[i])->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                arrProfv2PtSmall.at(pid[i])->Fill(p.Pt(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((6.0 < imp) && (imp < 7.0)) {
                                                arrProfv1RapMid.at(pid[i])->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                arrProfv2PtMid.at(pid[i])->Fill(p.Pt(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((10.0 < imp) && (imp < 12.0)) {
                                                arrProfv1RapBig.at(pid[i])->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                arrProfv2PtBig.at(pid[i])->Fill(p.Pt(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                }
                                if ((-1.0 < pL.Rapidity()) && (pL.Rapidity() < 1.0) && ((pL.Rapidity() < -0.2) || (0.2 < pL.Rapidity()))) {
                                        Short_t sign = 1;
                                        if (pL.Rapidity() < 0) { sign = -1; }
                                        if (imp < 3.0) {
                                                arrProfv1PtSmall.at(pid[i])->Fill(p.Pt(), sign*TMath::Cos(fi - psiRp));
                                        }
                                        else if ((6.0 < imp) && (imp < 7.0)) {
                                                arrProfv1PtMid.at(pid[i])->Fill(p.Pt(), sign*TMath::Cos(fi - psiRp));

                                        }
                                        else if ((10.0 < imp) && (imp < 12.0)) {
                                                arrProfv1PtBig.at(pid[i])->Fill(p.Pt(), sign*TMath::Cos(fi - psiRp));
                                        }
                                }
			}
		}

		//nb = fChain->GetEntry(jentry);   nbytes += nb;
		//if (Cut(ientry) < 0) continue;
	}
/*
	// Dependence v on ro = sqrt(x*x + y*y)
	for (int ix = 0; ix < nx; ix++) {
		for(int jy = 0; jy < ny; jy++) {
			for(int kz = 0; kz < nz; kz++) {
				if (arrV[ix][jy][kz]*arrV[ix][jy][kz] == 0) {
					continue;
				}
				for (Int_t i = 0; i < arrPIDs.size(); i++) {
					arrProfvRo.at(arrPIDs[i])->Fill(sqrt((ix - nx/2.0)*(ix - nx/2.0) + (jy - ny/2.0)*(jy - ny/2.0)), sqrt(arrV[ix][jy][kz]*arrV[ix][jy][kz]));
				}
				vRo->Fill(sqrt((ix - nx/2.0)*(ix - nx/2.0) + (jy - ny/2.0)*(jy - ny/2.0)), sqrt(arrV[ix][jy][kz]*arrV[ix][jy][kz]));
			}
		}
	}
*/
	MyFile->Write();
	MyFile->Close();
}
