#define Experiment_cxx
#include "Experiment.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
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


	// HIstograms
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
	/*// Momentum and energy in a cell
	map <Int_t, vector<pair<TVector3, TVector3> > > arrP;
	map <Int_t, vector<> >
	map <Int_t TH3D*> arrHistPx;
	map <Int_t TH3D*> arrHistPy;
	map <Int_t TH3D*> arrHistPz;
        TH3D *HistPx2 = new TH3D("HistPx2", "", nbins, -200.0, 200.0, nbins, -200.0, 200.0, nbins, -200.0, 200.0);
        TH3D *HistPy2 = new TH3D("HistPy2", "", nbins, -200.0, 200.0, nbins, -200.0, 200.0, nbins, -200.0, 200.0);
        TH3D *HistPz2 = new TH3D("HistPz2", "", nbins, -200.0, 200.0, nbins, -200.0, 200.0, nbins, -200.0, 200.0);*/




	// Probably better to do in vectors to have everything dynamically!!
	vector<Int_t> arrPIDs = {3101, 3106, 3001, 2027};
	string arrNames[] = {"Pion", "Kaon", "Proton", "Lambda"};

	MyFile->mkdir("Small parameter");
        MyFile->mkdir("Mid parameter");
	MyFile->mkdir("Big parameter");


	for (Int_t i = 0; i < arrPIDs.size(); i++) {
			arrHistPt[arrPIDs[i]] = new TH1D(("Pt" + arrNames[i]).c_str(), "-1.2 < Eta < 1.2\n0.2 < Pt < 3", 500, 0, 5);
			arrHistRap[arrPIDs[i]] =  new TH1D(("Rap" + arrNames[i]).c_str(), "-1.2 < Eta < 1.2\n0.2 < Pt < 3", 200, -5, 5); 
			arrProfvRo[arrPIDs[i]] = new TProfile(("vRo" + arrNames[i]).c_str(), "", 500, 0, sqrt(nx*nx/4.0 + ny*ny/4.0));       // Velocity from ro
        		arrProfvRo[arrPIDs[i]]->GetXaxis()->SetTitle("r, cells");

			MyFile->cd("Small parameter");
			arrProfv1PtSmall[arrPIDs[i]] = new TProfile(("v1PtSmall" + arrNames[i]).c_str(), "imp 0.0-3.0", NbinsPt, binsPt);
			arrProfv1PtSmall[arrPIDs[i]]->GetXaxis()->SetTitle("pt, GeV/c");
        		arrProfv1PtSmall[arrPIDs[i]]->GetYaxis()->SetTitle("v1");
			arrProfv2PtSmall[arrPIDs[i]] = new TProfile(("v2PtSmall" + arrNames[i]).c_str(), "imp 0.0-3.0", NbinsPt, binsPt);
                        arrProfv2PtSmall[arrPIDs[i]]->GetXaxis()->SetTitle("pt, GeV/c");
                        arrProfv2PtSmall[arrPIDs[i]]->GetYaxis()->SetTitle("v2");
                        arrProfv1RapSmall[arrPIDs[i]] = new TProfile(("v1RapSmall" + arrNames[i]).c_str(), "imp 0.0-3.0", NbinsY, binsY);
                        arrProfv1RapSmall[arrPIDs[i]]->GetXaxis()->SetTitle("Rap");
                        arrProfv1RapSmall[arrPIDs[i]]->GetYaxis()->SetTitle("v1");
                        arrProfv2RapSmall[arrPIDs[i]] = new TProfile(("v2RapSmall" + arrNames[i]).c_str(), "imp 0.0-3.0", NbinsY, binsY);
                        arrProfv2RapSmall[arrPIDs[i]]->GetXaxis()->SetTitle("Rap");
                        arrProfv2RapSmall[arrPIDs[i]]->GetYaxis()->SetTitle("v2");
			MyFile->cd("../");

        		MyFile->cd("Mid parameter");
			arrProfv1PtMid[arrPIDs[i]] = new TProfile(("v1PtMid" + arrNames[i]).c_str(), "imp 6.0-7.0", NbinsPt, binsPt);
                        arrProfv1PtMid[arrPIDs[i]]->GetXaxis()->SetTitle("pt, GeV/c");
                        arrProfv1PtMid[arrPIDs[i]]->GetYaxis()->SetTitle("v1");
                        arrProfv2PtMid[arrPIDs[i]] = new TProfile(("v2PtMid" + arrNames[i]).c_str(), "imp 6.0-7.0", NbinsPt, binsPt);
                        arrProfv2PtMid[arrPIDs[i]]->GetXaxis()->SetTitle("pt, GeV/c");
                        arrProfv2PtMid[arrPIDs[i]]->GetYaxis()->SetTitle("v2");
                        arrProfv1RapMid[arrPIDs[i]] = new TProfile(("v1RapMid" + arrNames[i]).c_str(), "imp 6.0-7.0", NbinsY, binsY);
                        arrProfv1RapMid[arrPIDs[i]]->GetXaxis()->SetTitle("Rap");
                        arrProfv1RapMid[arrPIDs[i]]->GetYaxis()->SetTitle("v1");
                        arrProfv2RapMid[arrPIDs[i]] = new TProfile(("v2RapMid" + arrNames[i]).c_str(), "imp 6.0-7.0", NbinsY, binsY);
                        arrProfv2RapMid[arrPIDs[i]]->GetXaxis()->SetTitle("Rap");
                        arrProfv2RapMid[arrPIDs[i]]->GetYaxis()->SetTitle("v2");
                        MyFile->cd("../");

        		MyFile->cd("Big parameter");
			arrProfv1PtBig[arrPIDs[i]] = new TProfile(("v1PtBig" + arrNames[i]).c_str(), "imp 10.0-12.0", NbinsPt, binsPt);
                        arrProfv1PtBig[arrPIDs[i]]->GetXaxis()->SetTitle("pt, GeV/c");
                        arrProfv1PtBig[arrPIDs[i]]->GetYaxis()->SetTitle("v1");
                        arrProfv2PtBig[arrPIDs[i]] = new TProfile(("v2PtBig" + arrNames[i]).c_str(), "imp 10.0-12.0", NbinsPt, binsPt);
                        arrProfv2PtBig[arrPIDs[i]]->GetXaxis()->SetTitle("pt, GeV/c");
                        arrProfv2PtBig[arrPIDs[i]]->GetYaxis()->SetTitle("v2");
                        arrProfv1RapBig[arrPIDs[i]] = new TProfile(("v1RapBig" + arrNames[i]).c_str(), "imp 10.0-12.0", NbinsY, binsY);
                        arrProfv1RapBig[arrPIDs[i]]->GetXaxis()->SetTitle("Rap");
                        arrProfv1RapBig[arrPIDs[i]]->GetYaxis()->SetTitle("v1");
                        arrProfv2RapBig[arrPIDs[i]] = new TProfile(("v2RapBig" + arrNames[i]).c_str(), "imp 10.0-12.0", NbinsY, binsY);
                        arrProfv2RapBig[arrPIDs[i]]->GetXaxis()->SetTitle("Rap");
                        arrProfv2RapBig[arrPIDs[i]]->GetYaxis()->SetTitle("v2");
                        MyFile->cd("../");





	}


        TProfile *vRo = new TProfile("vRo", "", 500, 0, sqrt(nx*nx/4.0 + ny*ny/4.0));       // Velocity from ro
	vRo->GetXaxis()->SetTitle("r, cells");

	TVector3 vNull(0, 0, 0);	// Null vector

	TVector3 arrP[nx][ny][nz], arrPPi[nx][ny][nz], arrPPr[nx][ny][nz], arrPKa[nx][ny][nz],  arrPLa[nx][ny][nz];

	// Total energy of particles
	Double_t arrE[nx][ny][nz], arrEPi[nx][ny][nz], arrEPr[nx][ny][nz], arrEKa[nx][ny][nz], arrELa[nx][ny][nz];
	for(Int_t i = 0; i < nx; i++) {
		for(Int_t j = 0; j < ny; j++) {
                	for(Int_t k = 0; k < nz; k++) {
                		arrE[i][j][k] = 0.0;	// Set energy as zero
				arrEPi[i][j][k] = 0.0;
				arrEPr[i][j][k] = 0.0;
				arrEKa[i][j][k] = 0.0;
				arrELa[i][j][k] = 0.0
        		}
        	}
	}



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
			rx[i] = rx1*cos(psiRp)+ry1*sin(psiRp);
			ry[i] = -rx1*sin(psiRp)+ry1*cos(psiRp);
			px1=px[i]; py1=py[i];
            		px[i] = px1*cos(psiRp)+py1*sin(psiRp);
            		py[i] = -px1*sin(psiRp)+py1*cos(psiRp);
			TVector3 p(px[i], py[i], pz[i]);
                        TLorentzVector pL(p, p0[i]);
			// Cuts for particles
			if ((p.Pt() < 0.2) || (3.0 < p.Pt()) || (pL.Eta() < -1.2) || (1.2 < pL.Eta())) {
				continue;
			}

			// Determine cell of particle
			int kx, ky, kz;
			kx = (rx[i] + X)/dx;
			ky = (ry[i] + Y)/dy;
			kz = (rz[i] + Z)/dz;

			// Change energy and momentum
			arrE[kx][ky][kz] += p0[i];
			arrP[kx][ky][kz] += p;
                        if (pid[i] == 3101) {
                                arrEPi[kx][ky][kz] += p0[i];
				arrPPi[kx][ky][kz] += p;
                        }
                        else if (pid[i] == 3101) {
                                arrEPi[kx][ky][kz] += p0[i];
                                arrP[kx][ky][kz] += p;
                        }
                        else if (pid[i] == 3101) {
                                arrEPi[kx][ky][kz] += p0[i];
                                arrP[kx][ky][kz] += p;
                        }
                        else if (pid[i] == 3101) {
                                arrEPi[kx][ky][kz] += p0[i];
                                arrP[kx][ky][kz] += p;
                        }



			// Draw histogramms
			if (arrHistPt.count(pid[i])) {
				// Draw histogramms for different types of particles
				arrHistPt.at(pid[i])->Fill(p.Pt());
	                        arrHistRap.at(pid[i])->Fill(pL.Rapidity());
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
	// Histograms for velicity and vorticity
        TH2D *HistV = new TH2D("HistV", "", nx, 0, nx, ny, 0, ny);
	HistV->GetXaxis()->SetTitle("X, cells");
        HistV->GetYaxis()->SetTitle("Y, cells");
        TH2D *HistVW = new TH2D("HistVW", "", nx, 0, nx, ny, 0, ny);       // v - Hr
        TH2D *HistW = new TH2D("HistW", "", nx, 0, nx, ny, 0, ny);
        HistV->GetXaxis()->SetTitle("X, cells");
        HistV->GetYaxis()->SetTitle("Y, cells");
        TH2D *HistPiV = new TH2D("HistPiV", "", nx, 0, nx, ny, 0, ny);
        TH2D *HistPrV = new TH2D("HistPrV", "", nx, 0, nx, ny, 0, ny);

	// Calculation of velocity
	TVector3 arrV[nx][ny][nz];
	TVector3 arrVPi[nx][ny][nz];
	TVector3 arrVPr[nx][ny][nz];
        for(int i = 0; i < nx; i++) {
                for(int j = 0; j < ny; j++) {
                        for(int k = 0; k < nz; k++) {
				if ((arrE[i][j][k] == 0) && (arrP[i][j][k]*arrP[i][j][k] != 0)) {
					cout << "Error: Energy = 0, P != 0 in cell ";
					cout << "(" << i << "; " << j << "; " << k << "); ";
					cout <<  "put zero vector here" << endl;
					arrV[i][j][k] = vNull;
				}
				else if(arrE[i][j][k] != 0) {
                                	arrV[i][j][k] = arrP[i][j][k]*(1.0/arrE[i][j][k]);
				}

                                if ((arrEPi[i][j][k] == 0) && (arrPPi[i][j][k]*arrPPi[i][j][k] != 0)) {
                                        cout << "Error: Energy = 0, P != 0 in cell ";
                                        cout << "(" << i << "; " << j << "; " << k << "); ";
                                        cout <<  "put zero vector here" << endl;
                                        arrVPi[i][j][k] = vNull;
                                }
                                else if(arrE[i][j][k] != 0) {
                                        arrVPi[i][j][k] = arrPPi[i][j][k]*(1.0/arrEPi[i][j][k]);
                                }
                                if ((arrEPr[i][j][k] == 0) && (arrPPr[i][j][k]*arrP[i][j][k] != 0)) {
                                        cout << "Error: Energy = 0, P != 0 in cell ";
                                        cout << "(" << i << "; " << j << "; " << k << "); ";
                                        cout <<  "put zero vector here" << endl;
                                        arrVPr[i][j][k] = vNull;
                                }
                                else if(arrE[i][j][k] != 0) {
                                        arrVPr[i][j][k] = arrPPr[i][j][k]*(1.0/arrEPr[i][j][k]);
                                }
                        }
                }
        }

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

	// Calculation of vorticity and helicity
	TVector3 arrW[nx][ny][nz];
	double H = 0.0;
        for(int i = 1; i < nx-1; i++) {
                for(int j = 1; j < ny-1; j++) {
			for(int k = 1; k < nz-1; k++) {
				if (arrV[i][j][k]*arrV[i][j][k] != 0)  {
					TVector3 dvi = arrV[i+1][j][k] - arrV[i-1][j][k];
					//cout << (arrV[i+1][j][k]*arrV[i+1][j][k]) << " " << (arrV[i-1][j][k]*arrV[i-1][j][k]) << endl;
                                	TVector3 dvj = arrV[i][j+1][k] - arrV[i][j-1][k];
                                	TVector3 dvk = arrV[i][j][k+1] - arrV[i][j][k-1];
                                        if ((arrV[i+1][j][k]*arrV[i+1][j][k])*(arrV[i-1][j][k]*arrV[i-1][j][k]) == 0) {
                                                dvi = vNull;
                                        }
                                        if ((arrV[i][j+1][k]*arrV[i][j+1][k])*(arrV[i][j-1][k]*arrV[i][j-1][k]) == 0) {
                                                dvj = vNull;
                                        }
                                        if ((arrV[i][j][k+1]*arrV[i][j][k+1])*(arrV[i][j][k-1]*arrV[i][j][k-1]) == 0) {
                                                dvk = vNull;
                                        }
					arrW[i][j][k].SetX((dvi.Z() + dvj.Z() + dvk.Z())/(6*dy) - (dvi.Y() + dvj.Y() + dvk.Y())/(6*dz));
                                	arrW[i][j][k].SetY((dvi.X() + dvj.X() + dvk.X())/(6*dz) - (dvi.Z() + dvj.Z() + dvk.Z())/(6*dx));
                                	arrW[i][j][k].SetZ((dvi.Y() + dvj.Y() + dvk.Y())/(6*dx) - (dvi.X() + dvj.X() + dvk.X())/(6*dy));
				}
				H += arrV[i][j][k]*arrW[i][j][k]*dx*dy*dz;      // Calculation of helicity
			}
                }
        }
	vRo->Write();

	// Fill histograms
	int Nz = 5;
	Double_t Hubble = 0.0025*dx;
	for(int x = 0; x < nx; x++) {
                for(int y = 0; y < ny; y++) {
			HistV->Fill(x, y, sqrt(arrV[x][y][Nz]*arrV[x][y][Nz])/986);
			TVector3 vH(Hubble*(x - nx/2.0), Hubble*(y - ny/2.0), 0);
			HistVW->Fill(x, y, sqrt((arrV[x][y][Nz] - vH)*(arrV[x][y][Nz] - vH)));
			HistW->Fill(x, y, sqrt(arrW[x][y][Nz]*arrW[x][y][Nz]));
                        HistPiV->Fill(x, y, sqrt(arrVPi[x][y][Nz]*arrVPi[x][y][Nz]));
			HistPrV->Fill(x, y, sqrt(arrVPr[x][y][Nz]*arrVPr[x][y][Nz]));
                }
        }

	HistV->Write();
	HistVW->Write();
	HistW->Write();
	MyFile->Write();

	MyFile->Close();

	cout << "Helicity = " << H << endl;	// We can try draw vn vs Helicity (+ narrow area on Pt)
}
