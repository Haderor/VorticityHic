// Use UrQMD withint freeze-out, test.f14 as output

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TDatime.h"
using namespace std;

int main() {
	ifstream fin("test1.f14");	// File with events
	if (!fin.is_open()) {
		cout << "Attached file was not opened!" << endl;
		exit(15);
	}
	string str;
	stringstream ss;

	// Histogramms for transverse momentums
	TH1D *HistPi = new TH1D("HistPi", "", 100, 0, 5);	// Pions
	TH1D *HistKa = new TH1D("HistKa", "", 100, 0, 5);	// Kaon
	TH1D *HistPr = new TH1D("HistPr", "", 100, 0, 5);	// Proton
	TH1D *HistLa = new TH1D("HistLa", "", 100, 0, 5);	// Lambda
	// TH1D *HistpsiRp = new TH1D("HistpsiRp", "", 100, 0, 7);

	TFile *MyFile = new TFile("plots.root", "recreate");	// File to save results
	int count = 17, Max_N = 0;						// Numbers of lines to skip in test.f14
	const double pi = 3.1415926;
	double psiRp, rx1, ry1, px1, py1;
	double r0, rx, ry, rz, p0, px, py, pz, m, pt;		// Data of a particle in an event
        TDatime* ti = new TDatime();
	TRandom3* fi = new TRandom3(ti->GetDate()*ti->GetTime());
	// Partition for space
	double X, Y, Z;		// Limits for coordinates
	X = Y = Z = 200.0;
	double dx, dy, dz;	// Splitting characteristic
	dx = dy = dz = 15.0;
	cout << "dx = " << dx << "; dy = " << dy << "; dz = " << dz << endl;
        int nx, ny, nz;         // Number of segmentss
        nx = 2*(X/dx) + 1;
        ny = 2*(Y/dy) + 1;
        nz = 2*(Z/dz) + 1;

	TVector3 vNull(0, 0, 0);	// Null vector

	TVector3 arrP[nx][ny][nz], arrPPi[nx][ny][nz], arrPPr[nx][ny][nz];
	for(int i = 0; i < nx; i++) {
                for(int j = 0; j < ny; j++) {
                        for(int k = 0; k < nz; k++) {
                                arrP[i][j][k] = vNull;	// Set momentum as zero
				arrPPi[i][j][k] = vNull;
				arrPPr[i][j][k] = vNull;
                        }
                }
        }
	// Total energy of particles
	double arrE[nx][ny][nz], arrEPi[nx][ny][nz], arrEPr[nx][ny][nz];
	for(int i = 0; i < nx; i++) {
		for(int j = 0; j < ny; j++) {
                	for(int k = 0; k < nz; k++) {
                		arrE[i][j][k] = 0.0;	// Set energy as zero
				arrEPi[i][j][k] = 0.0;
				arrEPr[i][j][k] = 0.0;
        		}
        	};
	}
	int pid, ityp, i3, ichg, lcl,  ncl, orr;

	// Loop on file
	while(!fin.eof()) {
		for (int j = 0; j < count; j++) {
			getline(fin, str);
		}
		if(fin.eof()) {
			break;
		}

		// Determine number of particles
		int n_particles, time;
		getline(fin, str);
		ss << str;
		ss >> n_particles >> time;
		if(n_particles > Max_N){
			Max_N = n_particles;
		}

		getline(fin, str);
		ss.str("");
		ss.clear();

		// Loop on particles in this event
		for(int j = 0; j < n_particles; j++) {
			getline(fin, str);
			ss << str;
			ss >> r0 >> rx >> ry >> rz >> p0 >> px >> py >> pz >> m >> ityp >> i3 >> ichg >> lcl >> ncl >> orr;
			pid = 1000 * (ichg + 2) + ityp;

			pt = sqrt(px*px + py*py);

			//Change coordinates with fi(reaction plane)
			psiRp = 2*pi*fi->Rndm();
			rx1=rx; ry1=ry;
			rx = rx1*cos(psiRp)+ry1*sin(psiRp);
			ry = -rx1*sin(psiRp)+ry1*cos(psiRp);
			px1=px; py1=py;
                        px = px1*cos(psiRp)+py1*sin(psiRp);
                        py = -px1*sin(psiRp)+py1*cos(psiRp);
			// HistpsiRp->Fill(psiRp);
			// Determine cell of particle
			int kx, ky, kz;
			kx = (rx + X)/dx;
			ky = (ry + Y)/dy;
			kz = (rz + Z)/dz;

			// Change energy and momentum
			arrE[kx][ky][kz] += p0;
			TVector3 p(px, py, pz);
			arrP[kx][ky][kz] += p;

			// Draw histogramms for different particles
			switch (pid) {
			case 3101:
				arrEPi[kx][ky][kz] += p0;
				arrPPi[kx][ky][kz] += p;
				HistPi->Fill(pt);
				break;
			case 3106:
				HistKa->Fill(pt);
				break;
			case 3001:
                                arrEPr[kx][ky][kz] += p0;
                                arrPPr[kx][ky][kz] += p;
				HistPr->Fill(pt);
				break;
			case 2027:
				HistLa->Fill(pt);
				break;
			}
			ss.str("");
			ss.clear();
		}
	}
	fin.close();

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
	// HistpsiRp->Write();
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
	return 0;
}

