#define Experiment_cxx
#include "Experiment.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
//#include "TVector3.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TMath.h"
//#include "TLorentzVector.h"
#include "TProfile.h"
#include <cmath>
#include <iostream>

using namespace std;

// r is distance between particle and point of field
// rad is radius of particle


/*Double_t FIfunc(Double_t r, Double_t rad) {
	Double_t sigma = rad/3.0;
	if (r<20) {
		//cout << -r*r/(2*sigma*sigma) << endl;
		cout << "    " << TMath::Exp(-r*r/(2*sigma*sigma)) << endl;
	}
	return TMath::Exp(-r*r/(2*sigma*sigma));
}*/

#pragma pack(push, 1)
struct TVector3 {
	Double_t x;
        Double_t y;
        Double_t z;
	TVector3() : x(0.0), y(0.0), z(0.0) {}
	TVector3(Double_t vx, Double_t vy, Double_t vz) : x(vx), y(vy), z(vz) {}
	TVector3 & operator += (const TVector3 & v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
        TVector3 & operator = (const TVector3 & v) {
                x = v.x;
                y = v.y;
                z = v.z;
                return *this;
        }
	Double_t  operator * (const TVector3 & v){
		return x * v.x + y * v.y + z * v.z;
	}
	TVector3 operator * (const Double_t a) {
		return TVector3(a*x, a*y, a*z);
	}
};
#pragma pack(pop)

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

	Double_t rx1, ry1, px1, py1, pt, fi, rap, eta;	//Temporary storage for coordinates and momentum
	Double_t X, Y, Z;		// Limits for coordinates
	X = Y = Z = 200.0;
	Double_t dx, dy, dz;	// Splitting characteristic
	dx = dy = 9.0;
	dz = dx;
	int nx, ny, nz;         // Number of segmentss
        nx = 2*(X/dx) + 1;
        ny = 2*(Y/dy) + 1;
        nz = 2*(Z/dz) + 1;

        TProfile *vRo = new TProfile("vRo", "", 500, 0, sqrt(nx*nx/4.0 + ny*ny/4.0));       // Velocity from ro

	TVector3 vNull(0, 0, 0);	// Null vector

	// Velocity in a cell
	TVector3 arrV[nx][ny][nz];	// Array with velocity
        Double_t arrCount[nx][ny][nz];  // Number of particles in a cell
	for(Int_t i = 0; i < nx; i++) {
		for(Int_t j = 0; j < ny; j++) {
                	for(Int_t k = 0; k < nz; k++) {
                		arrV[i][j][k] = vNull;	// Set velocity as zero
				arrCount[i][j][k] = 0.0;
        		}
        	}
	}

	Long64_t nbytes = 0; // nb = 0;
	const Double_t t0 = 200.0;
	const Double_t rad = 0.88;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		if(time != t0) {
			continue;
		}
		for (Int_t i = 0; i<n_particles; i++) {
			rx1=rx[i]; ry1=ry[i];
			rx[i] = rx1*cos(psiRp) + ry1*sin(psiRp);
			ry[i] = -rx1*sin(psiRp) + ry1*cos(psiRp);
			px1 = px[i]; py1 = py[i];
            		px[i] = px1*cos(psiRp) + py1*sin(psiRp);
            		py[i] = -px1*sin(psiRp) + py1*cos(psiRp);

			// Determine cell of the particle
                        int kx, ky, kz;
                        kx = (rx[i] + X)/dx;
                        ky = (ry[i] + Y)/dy;
                        kz = (rz[i] + Z)/dz;

			TVector3 p(px[i], py[i], pz[i]);
			arrV[kx][ky][kz] += p*(1.0/p0[i]);
			arrCount[kx][ky][kz]++;
			/*Int_t Nz = nz/2.0;
			for (int kx = 0; kx < nx; kx++) {
				for (int ky = 0; ky < nx; ky++) {
					for (int kz = 0; kz < nx; kz++) {
						// Change energy and momentum
                       				Double_t r = sqrt((kx + dx/2.0 - rx[i])*(kx + dx/2.0 - rx[i]) + (ky + dy/2.0 - ry[i])*(ky + dy/2.0 - ry[i]) + (kz + dz/2.0 - rz[i])*(kz + dz/2.0 - rz[i]));
                       				arrFI[kx][ky][kz] += FIfunc(r, rad);
                       				TVector3 p(px[i], py[i], pz[i]);
                       				arrPFI[kx][ky][kz] += p*arrFI[kx][ky][kz]*(1.0/p0[i]);
						if(kz == Nz) {
                        				HistV->Fill(x, y, sqrt(arrV[x][y][Nz]*arrV[x][y][Nz]));
						}
					}
				}
			}*/
		}
		//nb = fChain->GetEntry(jentry);   nbytes += nb;
		//if (Cut(ientry) < 0) continue;
	}

	// Make average in every cell
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
//				if (arrCount[i][j][k] > 0) {
//					arrV[i][j][k] = arrV[i][j][k]*(1.0/arrCount[i][j][k]);
//				}
			}
		}
	}

	// Histograms for velicity and vorticity
        TH2D *HistV = new TH2D("HistV", "", nx, 0, nx, ny, 0, ny);
        TH2D *HistW = new TH2D("HistW", "", nx, 0, nx, ny, 0, ny);

	/*// Calculation of velocity
	TVector3 arrV[nx][ny][nz];
	for (int kx = 0; kx < nx; kx++) {
        	for (int ky = 0; ky < nx; ky++) {
        		for (int kz = 0; kz < nx; kz++) {
				arrV[kx][ky][kz] = arrPFI[kx][ky][kz]*(1.0/arrFI[kx][ky][kz]);
        		}
        	}
        }*/

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
//	TVector3 arrW[nx][ny][nz];
/*	double H = 0.0;
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
        }*/
	vRo->Write();
	MyFile->Write();

	// Draw graph velocity and vorticity
	//TGraph2D *V = new TGraph2D();
	//TGraph2D *W = new TGraph2D();
        //int N = 0;
	const int Nz = 5;
	for(int x = 0; x < nx; x++) {
                for(int y = 0; y < ny; y++) {
                        //W->SetPoint(N, x, y, sqrt(arrW[x][y][Nz]*arrW[x][y][Nz]));
			//V->SetPoint(N, x, y, sqrt(arrV[x][y][Nz]*arrV[x][y][Nz]));
                        //N++;
			HistV->Fill(x, y, sqrt(arrV[x][y][Nz]*arrV[x][y][Nz]));
			//HistW->Fill(x, y, sqrt(arrW[x][y][Nz]*arrW[x][y][Nz]));
                }
        }

/*	// Draw graph v versus distance
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
*/
	//V->SetTitle("Magnitude of velocity; X; Y; Velocity");
	//W->SetTitle("Magnitude of vorticity; X; Y; Vorticity");
        //V->Write();
	//W->Write();
//	gr->Write();
	HistV->Write();
	//HistW->Write();
	MyFile->Close();
	//cout << "Helicity = " << H << endl;
	cout << "Max_N = " << Max_N << endl;

}
