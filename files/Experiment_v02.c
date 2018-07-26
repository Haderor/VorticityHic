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
	TH1D *HistPi = new TH1D("HistPi", "", 500, 0, 5);	// Pions
	TH1D *HistKa = new TH1D("HistKa", "", 500, 0, 5);	// Kaon
	TH1D *HistPr = new TH1D("HistPr", "", 500, 0, 5);	// Proton
        TH1D *HistPrPy = new TH1D("HistPrPy", "", 500, 0, 5);       // Proton Py
        TH1D *HistPrPx = new TH1D("HistPrPx", "", 500, 0, 5);       // Proton Px
        TH1D *HistPrPyBig = new TH1D("HistPrPyBig", "", 500, 0, 5);       // Proton Py big paramter
        TH1D *HistPrPxBig = new TH1D("HistPrPxBig", "", 500, 0, 5);       // Proton Px big parameter

	TH1D *HistLa = new TH1D("HistLa", "", 500, 0, 5);	// Lambda

    	TH1D *HistPsiRp = new TH1D("HistPsiRp", "", 500, 0, 6.3);	// psiRp

    	TH1D *HistFi = new TH1D("HistFi", "", 500, 0, 4);	// Angle fi between Pt and Px
        TH1D *HistFiPr = new TH1D("HistFiPr", "", 500, 0, 4);       // Angle fi between Pt and Px for protons

    	TH1D *HistCFi = new TH1D("HistCFi", "", 500, -1, 1);	// cos(fi)
    	TH1D *HistSFi = new TH1D("HistSFi", "", 500, -1, 1);	// sin(fi)

    	TH1D *HistFiPsi = new TH1D("HistFiPsi", "", 500, 0, 4);		// Fi-Psi
    	TH1D *HistCFiPsi = new TH1D("HistCFiPsi", "", 500, -1, 1);	// cos(Fi-Psi)
    	TH1D *HistSFiPsi = new TH1D("HistSFiPsi", "", 500, -1, 1);	// sin(Fi-Psi)

	TH1D *HistRapPi = new TH1D("HistRapPi", "", 200, -5, 5);	// Pions rapidity
	TH1D *HistRapPr = new TH1D("HistRapPr", "", 200, -5, 5);	// Protons rapidity

	TH1D *HistEtaPi = new TH1D("HistEtaPi", "", 200, -5, 5);	// Pions pseudorapidiyu
	TH1D *HistEtaPr = new TH1D("HistEtaPr", "", 200, -5, 5);	// Protons pseudorapidity

    	TProfile *v1Pi = new TProfile("v1Pi", "", 500, -5, 5); // Pions v1
    	TProfile *v2Pi = new TProfile("v2Pi", "", 500, -5, 5); // Pions v2

	const Double_t binsY[] = {-1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};
	const Double_t binsPt[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0};
	const Int_t NbinsY = 24;
	const Int_t NbinsPt = 11;
        MyFile->mkdir("Small parameter");
        MyFile->cd("Small parameter");
	// v1 and v2 for protons
        TProfile *v1PtPrSmall = new TProfile("v1PtPrSmall", "", NbinsPt, binsPt);
        v1PtPrSmall->GetXaxis()->SetTitle("pt, GeV/c");
        v1PtPrSmall->GetYaxis()->SetTitle("v1");
        TProfile *v2PtPrSmall = new TProfile("v2PtPrSmall", "", NbinsPt, binsPt);
        v2PtPrSmall->GetXaxis()->SetTitle("pt, GeV/c");
        v2PtPrSmall->GetYaxis()->SetTitle("v2");
        TProfile *v1RapPrSmall = new TProfile("v1RapPrSmall", "", NbinsY, binsY);
        v1RapPrSmall->GetXaxis()->SetTitle("rapidity");
        v1RapPrSmall->GetYaxis()->SetTitle("v1");
        TProfile *v2RapPrSmall = new TProfile("v2RapPrSmall", "", NbinsY, binsY);
        v1RapPrSmall->GetXaxis()->SetTitle("rapidity");
        v1RapPrSmall->GetYaxis()->SetTitle("v1");
        // v1 and v2 for pions
        TProfile *v1PtPiSmall = new TProfile("v1PtPiSmall", "", NbinsPt, binsPt);
        v1PtPiSmall->GetXaxis()->SetTitle("pt, GeV/c");
        v1PtPiSmall->GetYaxis()->SetTitle("v1");
        TProfile *v2PtPiSmall = new TProfile("v2PtPiSmall", "", NbinsPt, binsPt);
        v2PtPiSmall->GetXaxis()->SetTitle("pt, GeV/c");
        v2PtPiSmall->GetYaxis()->SetTitle("v2");
        TProfile *v1RapPiSmall = new TProfile("v1RapPiSmall", "", NbinsY, binsY);
        v1RapPiSmall->GetXaxis()->SetTitle("rapidity");
        v1RapPiSmall->GetYaxis()->SetTitle("v1");
        TProfile *v2RapPiSmall = new TProfile("v2RapPiSmall", "", NbinsY, binsY);
        v1RapPiSmall->GetXaxis()->SetTitle("rapidity");
        v1RapPiSmall->GetYaxis()->SetTitle("v1");
	MyFile->cd("../");

        MyFile->mkdir("Mid parameter");
        MyFile->cd("Mid parameter");
	// v1 and v2 for protons
        TProfile *v1RapPrMid = new TProfile("v1RapPrMid", "", NbinsY, binsY);
        v1RapPrMid->GetXaxis()->SetTitle("rapidity");
        v1RapPrMid->GetYaxis()->SetTitle("v1");
        TProfile *v2RapPrMid = new TProfile("v2RapPrMid", "", NbinsY, binsY);
        v2RapPrMid->GetXaxis()->SetTitle("rapidity");
        v2RapPrMid->GetYaxis()->SetTitle("v2");
        TProfile *v1PtPrMid = new TProfile("v1PtPrMid", "", NbinsPt, binsPt);
        v1PtPrMid->GetXaxis()->SetTitle("pt, GeV/c");
        v1PtPrMid->GetYaxis()->SetTitle("v1");
        TProfile *v2PtPrMid = new TProfile("v2PtPrMid", "", NbinsPt, binsPt);
        v2PtPrMid->GetXaxis()->SetTitle("pt, GeV/c");
        v2PtPrMid->GetYaxis()->SetTitle("v2");
	// v1 and v2 for pions
        TProfile *v1RapPiMid = new TProfile("v1RapPiMid", "", NbinsY, binsY);
        v1RapPiMid->GetXaxis()->SetTitle("rapidity");
        v1RapPiMid->GetYaxis()->SetTitle("v1");
        TProfile *v2RapPiMid = new TProfile("v2RapPiMid", "", NbinsY, binsY);
        v2RapPiMid->GetXaxis()->SetTitle("rapidity");
        v2RapPiMid->GetYaxis()->SetTitle("v2");
        TProfile *v1PtPiMid = new TProfile("v1PtPiMid", "", NbinsPt, binsPt);
        v1PtPiMid->GetXaxis()->SetTitle("pt, GeV/c");
        v1PtPiMid->GetYaxis()->SetTitle("v1");
        TProfile *v2PtPiMid = new TProfile("v2PtPiMid", "", NbinsPt, binsPt);
        v2PtPiMid->GetXaxis()->SetTitle("pt, GeV/c");
        v2PtPiMid->GetYaxis()->SetTitle("v2");
        MyFile->cd("../");

        MyFile->mkdir("Big parameter");
        MyFile->cd("Big parameter");
	// v1 and v2 for protons
        TProfile *v1RapPrBig = new TProfile("v1RapPrBig", "", NbinsY, binsY);
        v1RapPrBig->GetXaxis()->SetTitle("rapidity");
        v1RapPrBig->GetYaxis()->SetTitle("v1");
        TProfile *v2RapPrBig = new TProfile("v2RapPrBig", "", NbinsY, binsY);
        v2RapPrBig->GetXaxis()->SetTitle("rapidity");
        v2RapPrBig->GetYaxis()->SetTitle("v2");
        TProfile *v1PtPrBig = new TProfile("v1PtPrBig", "", NbinsPt, binsPt);
        v1PtPrBig->GetXaxis()->SetTitle("pt, GeV/c");
        v1PtPrBig->GetYaxis()->SetTitle("v1");
        TProfile *v2PtPrBig = new TProfile("v2PtPrBig", "", NbinsPt, binsPt);
        v2PtPrBig->GetXaxis()->SetTitle("pt, GeV/c");
        v2PtPrBig->GetYaxis()->SetTitle("v2");
	// v1 and v2 for pions
        TProfile *v1RapPiBig = new TProfile("v1RapPiBig", "", NbinsY, binsY);
        v1RapPiBig->GetXaxis()->SetTitle("rapidity");
        v1RapPiBig->GetYaxis()->SetTitle("v1");
        TProfile *v2RapPiBig = new TProfile("v2RapPiBig", "", NbinsY, binsY);
        v2RapPiBig->GetXaxis()->SetTitle("rapidity");
        v2RapPiBig->GetYaxis()->SetTitle("v2");
        TProfile *v1PtPiBig = new TProfile("v1PtPiBig", "", NbinsPt, binsPt);
        v1PtPiBig->GetXaxis()->SetTitle("pt, GeV/c");
        v1PtPiBig->GetYaxis()->SetTitle("v1");
        TProfile *v2PtPiBig = new TProfile("v2PtPiBig", "", NbinsPt, binsPt);
        v2PtPiBig->GetXaxis()->SetTitle("pt, GeV/c");
        v2PtPiBig->GetYaxis()->SetTitle("v2");
	MyFile->cd("../");

    	TProfile *v1Pr = new TProfile("v1Pr", "", 500, -5, 5); // Protons v1
    	TProfile *v2Pr = new TProfile("v2Pr", "", 500, -5, 5); // Protons v2

	Double_t rx1, ry1, px1, py1, pt, fi, rap, eta;	//Temporary storage for coordinates and momentum
	Double_t X, Y, Z;		// Limits for coordinates
	X = Y = Z = 200.0;
	Double_t dx, dy, dz;	// Splitting characteristic
	dx = dy = 20.0;
	dz = dx;
	int nx, ny, nz;         // Number of segmentss
        nx = 2*(X/dx) + 1;
        ny = 2*(Y/dy) + 1;
        nz = 2*(Z/dz) + 1;

        TProfile *vRo = new TProfile("vRo", "", 500, 0, sqrt(nx*nx/4.0 + ny*ny/4.0));       // Velocity from ro
	vRo->GetXaxis()->SetTitle("r, cells");

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
		psiRp = 0.0;
        	HistPsiRp->Fill(psiRp);
		if (psiRp != 0) {
			cout << "error: psiRp = " << psiRp << endl;
		}
		for (Int_t i = 0; i<n_particles; i++) {
			/*rx1=rx[i]; ry1=ry[i];
			rx[i] = rx1*cos(psiRp)+ry1*sin(psiRp);
			ry[i] = -rx1*sin(psiRp)+ry1*cos(psiRp);
			px1=px[i]; py1=py[i];
            		px[i] = px1*cos(psiRp)+py1*sin(psiRp);
            		py[i] = -px1*sin(psiRp)+py1*cos(psiRp);*/
			TVector3 p(px[i], py[i], pz[i]);
                        TLorentzVector pL(p, p0[i]);

			// Cuts for particles
			if ((p.Pt() < 0.2) || (3.0 < p.Pt()) || (pL.Eta() < -1.2) || (1.2 < pL.Eta())) {
				continue;
			}

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
				// Draw v2 for pions dependence
                                if ((0.2 < pt) && (pt < 3.0)) {
                                        if (imp < 3.0) {
                                                v2RapPiSmall->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((6.0 < imp) && (imp < 7.0)) {
                                                v2RapPiMid->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((10.0 < imp) && (imp < 12.0)) {
                                                v2RapPiBig->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
                                        }
                                }
                                if((-1.0 < pL.Rapidity()) && (pL.Rapidity() < 1.0)) {
                                        if (imp < 3.0) {
                                                v1RapPiSmall->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                v2PtPiSmall->Fill(pt, TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((6.0 < imp) && (imp < 7.0)) {
                                                v1RapPiMid->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                v2PtPiMid->Fill(pt, TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((10.0 < imp) && (imp < 12.0)) {
                                                v1RapPiBig->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                v2PtPiBig->Fill(pt, TMath::Cos(2*(fi - psiRp)));
                                        }
                                }
				if ((-1.0 < pL.Rapidity()) && (pL.Rapidity() < 1.0) && ((pL.Rapidity() < -0.2) || (0.2 < pL.Rapidity()))) {
                                        Short_t sign = 1;
                                        if (pL.Rapidity() < 0) { sign = -1; }
                                        if (imp < 3.0) {
                                                v1PtPiSmall->Fill(pt, sign*TMath::Cos(fi - psiRp));
                                        }
                                        else if ((6.0 < imp) && (imp < 7.0)) {
                                                v1PtPiMid->Fill(pt, sign*TMath::Cos(fi - psiRp));

                                        }
                                        else if ((10.0 < imp) && (imp < 12.0)) {
                                                v1PtPiBig->Fill(pt, sign*TMath::Cos(fi - psiRp));
                                        }
                                }
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
				// Draw v2 for protons dependence in rapidity
				if ((0.2 < pt) && (pt < 3.0)) {
					if (imp < 3.0) {
						v2RapPrSmall->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
					}
					else if ((6.0 < imp) && (imp < 7.0)) {
						v2RapPrMid->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
					}
					else if ((10.0 < imp) && (imp < 12.0)) {
                                        	v2RapPrBig->Fill(pL.Rapidity(), TMath::Cos(2*(fi - psiRp)));
					}
				}
				if((-1.0 < pL.Rapidity()) && (pL.Rapidity() < 1.0)) {
		                        if (imp < 3.0) {
                                                v1RapPrSmall->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                v2PtPrSmall->Fill(pt, TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((6.0 < imp) && (imp < 7.0)) {
                                                v1RapPrMid->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                v2PtPrMid->Fill(pt, TMath::Cos(2*(fi - psiRp)));
                                        }
                                        else if ((10.0 < imp) && (imp < 12.0)) {
                                                v1RapPrBig->Fill(pL.Rapidity(), TMath::Cos(fi - psiRp));
                                                v2PtPrBig->Fill(pt, TMath::Cos(2*(fi - psiRp)));
                                        }

				}
                                if ((-1.0 < pL.Rapidity()) && (pL.Rapidity() < 1.0) && ((pL.Rapidity() < -0.2) || (0.2 < pL.Rapidity()))) {
                                        Short_t sign = 1;
					if (pL.Rapidity() < 0) { sign = -1; }
					if (imp < 3.0) {
                                                v1PtPrSmall->Fill(pt, sign*TMath::Cos(fi - psiRp));
					}
                                        else if ((6.0 < imp) && (imp < 7.0)) {
                                                v1PtPrMid->Fill(pt, sign*TMath::Cos(fi - psiRp));

                                        }
                                        else if ((10.0 < imp) && (imp < 12.0)) {
                                                v1PtPrBig->Fill(pt, sign*TMath::Cos(fi - psiRp));
                                        }
                                }

				HistPr->Fill(pt);
	                        HistFiPr->Fill(TMath::ACos(TMath::Cos(fi)));
                                HistPrPy->Fill(py[i]);
                                HistPrPx->Fill(px[i]);
				if ((10.0 < imp) && (imp < 12.0)) {
                                	HistPrPyBig->Fill(py[i]);
                                	HistPrPxBig->Fill(px[i]);
				}
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
	HistV->GetXaxis()->SetTitle("X, cells");
        HistV->GetYaxis()->SetTitle("Y, cells");
        TH2D *HistVW = new TH2D("HistVW", "", nx, 0, nx, ny, 0, ny);
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

	// Dependence v on ro = sqrt(x*x + y*y)
	for (int ix = 0; ix < nx; ix++) {
		for(int jy = 0; jy < ny; jy++) {
			for(int kz = 0; kz < nz; kz++) {
				if (arrV[ix][jy][kz]*arrV[ix][jy][kz] == 0) {
					continue;
				}
				vRo->Fill(sqrt((ix - nx/2.0)*(ix - nx/2.0) + (jy - ny/2.0)*(jy - ny/2.0)), sqrt(arrV[ix][jy][kz]*arrV[ix][jy][kz]));
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
        HistPiV->Write();
        HistPrV->Write();
	HistPrPy->Write();
	HistPrPx->Write();
        HistPrPyBig->Write();
        HistPrPxBig->Write();
	HistFiPr->Write();
	MyFile->Write();

        MyFile->cd("Small parameter");
	v1PtPrSmall->Write();
        v2PtPrSmall->Write();
        v1RapPrSmall->Write();
        v2RapPrSmall->Write();
        v1PtPiSmall->Write();
        v2PtPiSmall->Write();
        v1RapPiSmall->Write();
        v2RapPiSmall->Write();
	MyFile->cd("../Mid parameter");
        v1PtPrMid->Write();
        v2PtPrMid->Write();
        v1RapPrMid->Write();
        v2RapPrMid->Write();
        v1PtPiMid->Write();
        v2PtPiMid->Write();
        v1RapPiMid->Write();
        v2RapPiMid->Write();
	MyFile->cd("../Big parameter");
	v1RapPrBig->Write();
        v2RapPrBig->Write();
        v1PtPrBig->Write();
        v2PtPrBig->Write();
        v1RapPiBig->Write();
        v2RapPiBig->Write();
        v1PtPiBig->Write();
        v2PtPiBig->Write();
        MyFile->cd("../");
	MyFile->Close();

	cout << "Helicity = " << H << endl;
}