#define Experiment_cxx
#include "Experiment.h"
#include <TStyle.h>
#include <TCanvas.h>
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

	// Histogramms for transverse momentums
	TH1D *HistPi = new TH1D("HistPi", "", 100, 0, 5);	// Pions
	TH1D *HistKa = new TH1D("HistKa", "", 100, 0, 5);	// Kaon
	TH1D *HistPr = new TH1D("HistPr", "", 100, 0, 5);	// Proton
	TH1D *HistLa = new TH1D("HistLa", "", 100, 0, 5);	// Lambda

    	TH1D *HistPsiRp = new TH1D("HistPsiRp", "", 100, 0, 6.3);	// psiRp

    	TH1D *HistFi = new TH1D("HistFi", "", 100, 0, 4);	// Angle fi between P and Px
    	TH1D *HistCFi = new TH1D("HistCFi", "", 100, -1, 1);	// cos(fi)
    	TH1D *HistSFi = new TH1D("HistSFi", "", 100, -1, 1);	// sin(fi)

    	TH1D *HistFiPsi = new TH1D("HistFiPsi", "", 100, 0, 4);	// Fi-Psi
    	TH1D *HistCFiPsi = new TH1D("HistCFiPsi", "", 100, -1, 1);	// cos(Fi-Psi)
    	TH1D *HistSFiPsi = new TH1D("HistSFiPsi", "", 100, -1, 1);	// sin(Fi-Psi)

	TH1D *HistRapPi = new TH1D("HistRapPi", "", 200, -5, 5);	// Pions rapidity
	TH1D *HistRapPr = new TH1D("HistRapPr", "", 200, -5, 5);	// Protons rapidity

	TH1D *HistEtaPi = new TH1D("HistEtaPi", "", 200, -5, 5);	// Pions pseudorapidiyu
	TH1D *HistEtaPr = new TH1D("HistEtaPr", "", 200, -5, 5);	// Protons pseudorapidity

    	TProfile *v1Pi = new TProfile("v1Pi", "", 100, -5, 5); // Pions v1
    	TProfile *v2Pi = new TProfile("v2Pi", "", 100, -5, 5); // Pions v2

    	TProfile *v1Pr = new TProfile("v1Pr", "", 100, -5, 5); // Protons v1
    	TProfile *v2Pr = new TProfile("v2Pr", "", 100, -5, 5); // Protons v2

	Double_t rx1, ry1, px1, py1, pt, fi, rap, eta;	//Temporary storage for coordinates and momentum
	Double_t X, Y, Z;		// Limits for coordinates
	X = Y = Z = 200.0;
	Double_t dx, dy, dz;	// Splitting characteristic
	dx = dy = dz = 20.0;
	int nx, ny, nz;         // Number of segmentss
        nx = 2*(X/dx) + 1;
        ny = 2*(Y/dy) + 1;
        nz = 2*(Z/dz) + 1;

	TVector3 vNull(0, 0, 0);	// Null vector

	TVector3 arrP[nx][ny][nz], arrPPi[nx][ny][nz], arrPPr[nx][ny][nz];
	for(Int_t i = 0; i < nx; i++) {
                for(Int_t j = 0; j < ny; j++) {
                        for(Int_t k = 0; k < nz; k++) {
                                arrP[i][j][k] = vNull;	// Set momentum as zero
				arrPPi[i][j][k] = vNull;
				arrPPr[i][j][k] = vNull;
                        }
                }
        }
	// Total energy of particles
	Double_t arrE[nx][ny][nz], arrEPi[nx][ny][nz], arrEPr[nx][ny][nz];
	for(Int_t i = 0; i < nx; i++) {
		for(Int_t j = 0; j < ny; j++) {
                	for(Int_t k = 0; k < nz; k++) {
                		arrE[i][j][k] = 0.0;	// Set energy as zero
				arrEPi[i][j][k] = 0.0;
				arrEPr[i][j][k] = 0.0;
        		}
        	};
	}



	Long64_t nbytes = 0; // nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		if(time != 200.0) {
			continue;
		}
        	HistPsiRp->Fill(psiRp);
		for (Int_t i = 0; i<n_particles; i++) {
			rx1=rx[i]; ry1=ry[i];
			rx[i] = rx1*cos(psiRp)+ry1*sin(psiRp);
			ry[i] = -rx1*sin(psiRp)+ry1*cos(psiRp);
			px1=px[i]; py1=py[i];
            		px[i] = px1*cos(psiRp)+py1*sin(psiRp);
            		py[i] = -px1*sin(psiRp)+py1*cos(psiRp);
            		fi = TMath::ATan2(py[i], px[i]);
            		HistFi->Fill(TMath::ACos(TMath::Cos(fi)));
            		HistCFi->Fill(TMath::Cos(fi));
            		HistSFi->Fill(TMath::Sin(fi));
            		HistFiPsi->Fill(TMath::ACos(TMath::Cos(fi-psiRp)));
            		HistCFiPsi->Fill(TMath::Cos(fi - psiRp));
           		HistSFiPsi->Fill(TMath::Sin(fi - psiRp));

			// Determine cell of particle
			int kx, ky, kz;
			kx = (rx[i] + X)/dx;
			ky = (ry[i] + Y)/dy;
			kz = (rz[i] + Z)/dz;

			// Change energy and momentum
			arrE[kx][ky][kz] += p0[i];
			TVector3 p(px[i], py[i], pz[i]);
            		TLorentzVector pL(p, p0[i]);
			arrP[kx][ky][kz] += p;
			pt = sqrt(px[i]*px[i] + py[i]*py[i]); //transverse momentum

			// Draw histogramms for different particles
			switch (pid[i]) {
			case 3101:
				arrEPi[kx][ky][kz] += p0[i];
				arrPPi[kx][ky][kz] += p;
				HistRapPi->Fill(pL.Rapidity());
            			HistEtaPi->Fill(pL.Eta());
               			v1Pi->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
               			v2Pi->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
				HistPi->Fill(pt);
				break;
			case 3106:
				HistKa->Fill(pt);
				break;
			case 3001:
               			arrEPr[kx][ky][kz] += p0[i];
               			arrPPr[kx][ky][kz] += p;
               			HistRapPr->Fill(pL.Rapidity());
               			HistEtaPr->Fill(pL.Eta());
               			v1Pr->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
               			v2Pr->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
				HistPr->Fill(pt);
				break;
			case 2027:
				HistLa->Fill(pt);
				break;
			}
		}
		//nb = fChain->GetEntry(jentry);   nbytes += nb;
		//if (Cut(ientry) < 0) continue;
	}
	// Histograms for velicity and vorticity
        TH2D *HistV = new TH2D("HistV", "", nx, 0, nx, ny, 0, ny);
        TH2D *HistW = new TH2D("HistW", "", nx, 0, nx, ny, 0, ny);
        TH2D *HistPiV = new TH2D("HistPiV", "", nx, 0, nx, ny, 0, ny);
        TH2D *HistPrV = new TH2D("HistPrV", "", nx, 0, nx, ny, 0, ny);

	// Calculation of velocity
	TVector3 arrV[nx][ny][nz];
	TVector3 arrVPi[nx][ny][nz];
	TVector3 arrVPr[nx][ny][nz];
        for(int i = 0; i < nx; i++) {
                for(int j = 0; j < ny; j++) {
                        for(int k = 0; k < nz; k++) {
				if((arrE[i][j][k] == 0) && (arrP[i][j][k]*arrP[i][j][k] == 0)) {
					arrV[i][j][k] = vNull;
				}
				else if ((arrE[i][j][k] == 0) && (arrP[i][j][k]*arrP[i][j][k] != 0)) {
					cout << "Error: Energy = 0, P != 0 in cell ";
					cout << "(" << i << "; " << j << "; " << k << "); ";
					cout <<  "put zero vector here" << endl;
					arrV[i][j][k] = vNull;
				}
				else {
                                	arrV[i][j][k] = arrP[i][j][k]*(1.0/arrE[i][j][k]);
				}

                                if((arrEPi[i][j][k] == 0) && (arrPPi[i][j][k]*arrPPi[i][j][k] == 0)) {
                                        arrVPi[i][j][k] = vNull;
                                }
                                else if ((arrEPi[i][j][k] == 0) && (arrPPi[i][j][k]*arrPPi[i][j][k] != 0)) {
                                        cout << "Error: Energy = 0, P != 0 in cell ";
                                        cout << "(" << i << "; " << j << "; " << k << "); ";
                                        cout <<  "put zero vector here" << endl;
                                        arrVPi[i][j][k] = vNull;
                                }
                                else {
                                        arrVPi[i][j][k] = arrPPi[i][j][k]*(1.0/arrEPi[i][j][k]);
                                }
                                if((arrEPr[i][j][k] == 0) && (arrPPr[i][j][k]*arrPPr[i][j][k] == 0)) {
                                        arrVPr[i][j][k] = vNull;
                                }
                                else if ((arrEPr[i][j][k] == 0) && (arrPPr[i][j][k]*arrP[i][j][k] != 0)) {
                                        cout << "Error: Energy = 0, P != 0 in cell ";
                                        cout << "(" << i << "; " << j << "; " << k << "); ";
                                        cout <<  "put zero vector here" << endl;
                                        arrVPr[i][j][k] = vNull;
                                }
                                else {
                                        arrVPr[i][j][k] = arrPPr[i][j][k]*(1.0/arrEPr[i][j][k]);
                                }

                        }
                }
        }

	// Calculation of vorticity and helicity
	TVector3 arrW[nx][ny][nz];
	double H = 0.0;
        for(int i = 0; i < nx; i++) {
                for(int j = 0; j < ny; j++) {
                        for(int k = 0; k < nz; k++) {
                                arrW[i][j][k] = vNull;    // Set vorticity as zero
                               // arrWPi[i][j][k] = 0.0;
                               // arrWPr[i][j][k] = 0.0;
                        }
                }
        }
        for(int i = 1; i < nx-1; i++) {
                for(int j = 1; j < ny-1; j++) {
			for(int k = 1; k < nz-1; k++) {
				if (arrV[i][j][k]*arrV[i][j][k] == 0) {
					arrW[i][j][k] = vNull;
				}
				else {
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
	HistPi->Write();
	HistKa->Write();
	HistPr->Write();
	HistLa->Write();
	HistFi->Write();
	HistRapPi->Write();
	HistEtaPi->Write();
	v1Pi->Write();
	v2Pi->Write();
	HistRapPr->Write();
	HistEtaPr->Write();
	v1Pr->Write();
	v2Pr->Write();
	HistCFiPsi->Write();
	HistSFiPsi->Write();
	HistFiPsi->Write();
	HistCFiPsi->Write();
	HistSFiPsi->Write();
	HistPsiRp->Write();
	MyFile->Write();

	// Draw velocity and vorticity
	TGraph2D *V = new TGraph2D();
	TGraph2D *W = new TGraph2D();
        int N = 0;
	int Nz = 5;
	for(int x = 0; x < nx; x++) {
                for(int y = 0; y < ny; y++) {
                        W->SetPoint(N, x, y, sqrt(arrW[x][y][Nz]*arrW[x][y][Nz]));
			V->SetPoint(N, x, y, sqrt(arrV[x][y][Nz]*arrV[x][y][Nz]));
                        N++;
			HistV->Fill(x, y, sqrt(arrV[x][y][Nz]*arrV[x][y][Nz]));
			HistW->Fill(x, y, sqrt(arrW[x][y][Nz]*arrW[x][y][Nz]));
                        HistPiV->Fill(x, y, sqrt(arrVPi[x][y][Nz]*arrVPi[x][y][Nz]));
			HistPrV->Fill(x, y, sqrt(arrVPr[x][y][Nz]*arrVPr[x][y][Nz]));
                }
        }

	// Draw v versus distance
	int x = nx/2;
	int z = nz/2;
	double arrVvsRo[ny - 0];
	double arrRo[ny - 0];
	for (int y = 0; y < ny; y++) {
		arrVvsRo[y - 0] = sqrt(arrV[x][y][z]*arrV[x][y][z]);
		arrRo[y - 0] = y - ny/2;
	}
	TGraph *gr = new TGraph(ny - 0, arrRo, arrVvsRo);

        gStyle->SetPalette(1);
	V->SetTitle("Magnitude of velocity; X; Y; Velocity");
	W->SetTitle("Magnitude of vorticity; X; Y; Vorticity");
        V->Write();
	W->Write();
	gr->Write();
	HistV->Write();
	HistW->Write();
        HistPiV->Write();
        HistPrV->Write();
	MyFile->Close();
	cout << "Helicity = " << H << endl;
	cout << "Max_N = " << Max_N << endl;

}
