// Adds to result.root histogramms with vx, vy, vz and v for all particles (counted as total momentum / total energy)
// Correction needed in determing edges for histogramm v
#include <vector>
#include <string>
#include "TH3D.h"

void GetVel() {

    // Information about types of particles
    const vector<Int_t> arrPIDs = {0, 3101, 3106, 3001, 2027};
    const vector<string> arrNames = {"All", "Pion", "Kaon", "Proton", "Lambda"};
    if(arrPIDs.size() != arrNames.size()) {
        cout << "Error: numbers of PIDs and particle names are not equal" << endl;
        exit(1);
    }

    // Histogramms in root file
    map <Int_t, TH3D*> arrHistPx;
    map <Int_t, TH3D*> arrHistPy;
    map <Int_t, TH3D*> arrHistPz;
    map <Int_t, TH3D*> arrHistE;

    // Now we don't need Vx_2, but probably we will need it later
    /*map <Int_t, TH3D*> arrHistVx_2;
    map <Int_t, TH3D*> arrHistVy_2;
    map <Int_t, TH3D*> arrHistVz_2;*/

    // Getting histogramms from root file
    TFile *f = new TFile("result.root", "UPDATE");
    for (Int_t i = 0; i < arrPIDs.size(); i++) {
        arrHistPx[arrPIDs[i]] = (TH3D*)f->Get(("PxSmall" + arrNames[i]).c_str());
        arrHistPy[arrPIDs[i]] = (TH3D*)f->Get(("PySmall" + arrNames[i]).c_str());
        arrHistPz[arrPIDs[i]] = (TH3D*)f->Get(("PzSmall" + arrNames[i]).c_str());
        arrHistE[arrPIDs[i]] = (TH3D*)f->Get(("ESmall" + arrNames[i]).c_str());

        /*arrHistVx_2[arrPIDs[i]] = (TH3D*)f->Get(("Vx_2Small" + arrNames[i]).c_str());
        arrHistVy_2[arrPIDs[i]] = (TH3D*)f->Get(("Vy_2Small" + arrNames[i]).c_str());
        arrHistVz_2[arrPIDs[i]] = (TH3D*)f->Get(("Vz_2Small" + arrNames[i]).c_str());*/
    }

    // Determing and saving velocity with 1-st way; i.e. sum P / sum E on all events and all particles
    //---------------------------------------------------------------------------------------------
    // Histogramms for velocities, calculated with 1-st way
    map <Int_t, TH3D*> arrHistVx_1;
    map <Int_t, TH3D*> arrHistVy_1;
    map <Int_t, TH3D*> arrHistVz_1;

    // Get velocity with the 1-st way and save to file
    for(Int_t i = 0; i < arrPIDs.size(); i++) {
        arrHistVx_1[arrPIDs[i]] = new TH3D(*arrHistPx[arrPIDs[i]]);
        if (!arrHistVx_1[arrPIDs[i]]->Divide(arrHistE[arrPIDs[i]])) {
            cout << "Error, division by 0" << endl;
            exit(2);
        }
        arrHistVx_1[arrPIDs[i]]->SetName(("Vx_1Small" + arrNames[i]).c_str());
        arrHistVx_1[arrPIDs[i]]->Write();

        arrHistVy_1[arrPIDs[i]] = new TH3D(*arrHistPy[arrPIDs[i]]);
        if (!arrHistVy_1[arrPIDs[i]]->Divide(arrHistE[arrPIDs[i]])) {
            cout << "Error, division by 0" << endl;
            exit(2);
        }
        arrHistVy_1[arrPIDs[i]]->SetName(("Vy_1Small" + arrNames[i]).c_str());
        arrHistVy_1[arrPIDs[i]]->Write();

        arrHistVz_1[arrPIDs[i]] = new TH3D(*arrHistPz[arrPIDs[i]]);
        if (!arrHistVz_1[arrPIDs[i]]->Divide(arrHistE[arrPIDs[i]])) {
            cout << "Error, division by 0" << endl;
            exit(2);
        }
        arrHistVz_1[arrPIDs[i]]->SetName(("Vz_1Small" + arrNames[i]).c_str());
        arrHistVz_1[arrPIDs[i]]->Write();
    }

    // Get v = |(vx, vy, vz)| and save to file
    for(Int_t i = 0; i < arrPIDs.size(); i++) {
        Int_t NbinsX = arrHistVx_1[arrPIDs[i]]->GetNbinsX();
        Int_t NbinsY = arrHistVx_1[arrPIDs[i]]->GetNbinsY();
        Int_t NbinsZ = arrHistVx_1[arrPIDs[i]]->GetNbinsZ();

        Double_t xmin = arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetBinLowEdge(arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetFirst());
        Double_t xmax = arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetBinUpEdge(arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetLast());
        Double_t ymin = arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetBinLowEdge(arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetFirst());
        Double_t ymax = arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetBinUpEdge(arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetLast());
        Double_t zmin = arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetBinLowEdge(arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetFirst());
        Double_t zmax = arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetBinUpEdge(arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetLast());

        TH3D *arrHistV_1 = new TH3D(("V_1Small" + arrNames[i]).c_str(), "Add cuts here;x, fm;y, fm;z, fm", NbinsX, xmin, xmax, NbinsY, ymin, ymax, NbinsZ, zmin, zmax);
        TProfile *arrProfVRo = new TProfile(("VRo" + arrNames[i]).c_str(), "Add cuts here;Ro, fm", 100, 0, TMath::Sqrt(xmax*xmax + ymax*ymax));
        for(Int_t x = 1; x <= NbinsX; x++) {
            for(Int_t y = 1; y <= NbinsY; y++) {
                for(Int_t z = 1; z <= NbinsZ; z++) {
                    Int_t bin = arrHistV_1->GetBin(x, y, z);
                    Int_t profBin = arrProfVRo->GetBin(TMath::Sqrt(arrHistV_1->GetXaxis()->GetBinCenter(x)*arrHistV_1->GetXaxis()->GetBinCenter(x) + arrHistV_1->GetYaxis()->GetBinCenter(y)*arrHistV_1->GetYaxis()->GetBinCenter(y)));
                    arrHistV_1->SetBinContent(bin, TMath::Sqrt(arrHistVx_1[arrPIDs[i]]->GetBinContent(bin)*arrHistVx_1[arrPIDs[i]]->GetBinContent(bin) + arrHistVy_1[arrPIDs[i]]->GetBinContent(bin)*arrHistVy_1[arrPIDs[i]]->GetBinContent(bin) + arrHistVy_1[arrPIDs[i]]->GetBinContent(bin)*arrHistVy_1[arrPIDs[i]]->GetBinContent(bin)));
                    arrProfVRo->SetBinContent(profBin, arrHistV_1->GetBinContent(bin));
                }
            }
        }
        arrProfVRo->Write();
        arrHistV_1->Write();
    }

    // Get vorticity with the 1-st way and save to file
    for(Int_t i = 0; i < arrPIDs.size(); i++) {
        Int_t NbinsX = arrHistVx_1[arrPIDs[i]]->GetNbinsX();
        Int_t NbinsY = arrHistVx_1[arrPIDs[i]]->GetNbinsY();
        Int_t NbinsZ = arrHistVx_1[arrPIDs[i]]->GetNbinsZ();
        Double_t xmin = arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetBinLowEdge(arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetFirst());
        Double_t xmax = arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetBinUpEdge(arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetLast());
        Double_t ymin = arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetBinLowEdge(arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetFirst());
        Double_t ymax = arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetBinUpEdge(arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetLast());
        Double_t zmin = arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetBinLowEdge(arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetFirst());
        Double_t zmax = arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetBinUpEdge(arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetLast());

        TH3D *arrHistWx_1 = new TH3D(("Wx_1Small" + arrNames[i]).c_str(), "Add cuts here;x, fm;y, fm;z, fm", NbinsX, xmin, xmax, NbinsY, ymin, ymax, NbinsZ, zmin, zmax);
        TH3D *arrHistWy_1 = new TH3D(("Wy_1Small" + arrNames[i]).c_str(), "Add cuts here;x, fm;y, fm;z, fm", NbinsX, xmin, xmax, NbinsY, ymin, ymax, NbinsZ, zmin, zmax);
        TH3D *arrHistWz_1 = new TH3D(("Wz_1Small" + arrNames[i]).c_str(), "Add cuts here;x, fm;y, fm;z, fm", NbinsX, xmin, xmax, NbinsY, ymin, ymax, NbinsZ, zmin, zmax);
        TH3D *arrHistW_1 = new TH3D(("W_1Small" + arrNames[i]).c_str(), "Add cuts here;x, fm;y, fm;z, fm", NbinsX, xmin, xmax, NbinsY, ymin, ymax, NbinsZ, zmin, zmax);
        for(Int_t x = 2; x < NbinsX; x++) {
            for(Int_t y = 2; y < NbinsY; y++) {
                for(Int_t z = 2; z < NbinsZ; z++) {
                    Int_t bin = arrHistVx_1[arrPIDs[i]]->GetBin(x, y, z);
                    Int_t binxUp = arrHistVx_1[arrPIDs[i]]->GetBin(x+1, y, z);
                    Int_t binxLow = arrHistVx_1[arrPIDs[i]]->GetBin(x-1, y, z);
                    Int_t binyUp = arrHistVx_1[arrPIDs[i]]->GetBin(x, y+1, z);
                    Int_t binyLow = arrHistVx_1[arrPIDs[i]]->GetBin(x, y-1, z);
                    Int_t binzUp = arrHistVx_1[arrPIDs[i]]->GetBin(x, y, z+1);
                    Int_t binzLow = arrHistVx_1[arrPIDs[i]]->GetBin(x, y, z-1);

                    Double_t dx = arrHistVx_1[arrPIDs[i]]->GetXaxis()->GetBinWidth(x);
                    Double_t dy = arrHistVx_1[arrPIDs[i]]->GetYaxis()->GetBinWidth(y);
                    Double_t dz = arrHistVx_1[arrPIDs[i]]->GetZaxis()->GetBinWidth(z);

                    // If there is no particles in the cell
                    if(!arrHistE[arrPIDs[i]]->GetBinContent(bin)) {
                        arrHistWx_1->SetBinContent(bin, 0.0);
                        arrHistWy_1->SetBinContent(bin, 0.0);
                        arrHistWz_1->SetBinContent(bin, 0.0);
                        continue;
                    }

                    // Derivatives for curl(v)
                    Double_t dvz_dy = 0.0;
                    Double_t dvy_dz = 0.0;
                    Double_t dvx_dz = 0.0;
                    Double_t dvz_dx = 0.0;
                    Double_t dvy_dx = 0.0;
                    Double_t dvx_dy = 0.0;

                    // dvz_dy
                    if(arrHistE[arrPIDs[i]]->GetBinContent(binyUp) && !arrHistE[arrPIDs[i]]->GetBinContent(binyLow)) {
                        dvz_dy = (arrHistVz_1[arrPIDs[i]]->GetBinContent(binyUp) - arrHistVz_1[arrPIDs[i]]->GetBinContent(bin))/dy;
                    }
                    else if (!arrHistE[arrPIDs[i]]->GetBinContent(binyUp) && arrHistE[arrPIDs[i]]->GetBinContent(binyLow)) {
                        dvz_dy = (arrHistVz_1[arrPIDs[i]]->GetBinContent(bin) - arrHistVz_1[arrPIDs[i]]->GetBinContent(binyLow))/dy;
                    }
                    else if(arrHistE[arrPIDs[i]]->GetBinContent(binyUp) && arrHistE[arrPIDs[i]]->GetBinContent(binyLow)) {
                        dvz_dy = (arrHistVz_1[arrPIDs[i]]->GetBinContent(binyUp) - arrHistVz_1[arrPIDs[i]]->GetBinContent(binyLow))/(2*dy);
                    }

                    // dvy_dz
                    if(arrHistE[arrPIDs[i]]->GetBinContent(binzUp) && !arrHistE[arrPIDs[i]]->GetBinContent(binzLow)) {
                        dvy_dz = (arrHistVy_1[arrPIDs[i]]->GetBinContent(binzUp) - arrHistVy_1[arrPIDs[i]]->GetBinContent(bin))/dz;
                    }
                    else if (!arrHistE[arrPIDs[i]]->GetBinContent(binzUp) && arrHistE[arrPIDs[i]]->GetBinContent(binzLow)) {
                        dvy_dz = (arrHistVy_1[arrPIDs[i]]->GetBinContent(bin) - arrHistVy_1[arrPIDs[i]]->GetBinContent(binzLow))/dz;
                    }
                    else if(arrHistE[arrPIDs[i]]->GetBinContent(binzUp) && arrHistE[arrPIDs[i]]->GetBinContent(binzLow)) {
                        dvy_dz = (arrHistVy_1[arrPIDs[i]]->GetBinContent(binzUp) - arrHistVy_1[arrPIDs[i]]->GetBinContent(binzLow))/(2*dz);
                    }

                    // dvx_dz
                    if(arrHistE[arrPIDs[i]]->GetBinContent(binzUp) && !arrHistE[arrPIDs[i]]->GetBinContent(binzLow)) {
                        dvx_dz = (arrHistVx_1[arrPIDs[i]]->GetBinContent(binzUp) - arrHistVx_1[arrPIDs[i]]->GetBinContent(bin))/dz;
                    }
                    else if (!arrHistE[arrPIDs[i]]->GetBinContent(binzUp) && arrHistE[arrPIDs[i]]->GetBinContent(binzLow)) {
                        dvx_dz = (arrHistVx_1[arrPIDs[i]]->GetBinContent(bin) - arrHistVx_1[arrPIDs[i]]->GetBinContent(binzLow))/dz;
                    }
                    else if(arrHistE[arrPIDs[i]]->GetBinContent(binzUp) && arrHistE[arrPIDs[i]]->GetBinContent(binzLow)) {
                        dvx_dz = (arrHistVx_1[arrPIDs[i]]->GetBinContent(binzUp) - arrHistVx_1[arrPIDs[i]]->GetBinContent(binzLow))/(2*dz);
                    }

                    // dvz_dx
                    if(arrHistE[arrPIDs[i]]->GetBinContent(binxUp) && !arrHistE[arrPIDs[i]]->GetBinContent(binxLow)) {
                        dvz_dx = (arrHistVz_1[arrPIDs[i]]->GetBinContent(binxUp) - arrHistVz_1[arrPIDs[i]]->GetBinContent(bin))/dx;
                    }
                    else if (!arrHistE[arrPIDs[i]]->GetBinContent(binxUp) && arrHistE[arrPIDs[i]]->GetBinContent(binxLow)) {
                        dvz_dx = (arrHistVz_1[arrPIDs[i]]->GetBinContent(bin) - arrHistVz_1[arrPIDs[i]]->GetBinContent(binxLow))/dx;
                    }
                    else if(arrHistE[arrPIDs[i]]->GetBinContent(binxUp) && arrHistE[arrPIDs[i]]->GetBinContent(binxLow)) {
                        dvz_dx = (arrHistVz_1[arrPIDs[i]]->GetBinContent(binxUp) - arrHistVz_1[arrPIDs[i]]->GetBinContent(binxLow))/(2*dx);
                    }

                    // dvy_dx
                    if(arrHistE[arrPIDs[i]]->GetBinContent(binxUp) && !arrHistE[arrPIDs[i]]->GetBinContent(binxLow)) {
                        dvy_dx = (arrHistVy_1[arrPIDs[i]]->GetBinContent(binxUp) - arrHistVy_1[arrPIDs[i]]->GetBinContent(bin))/dx;
                    }
                    else if (!arrHistE[arrPIDs[i]]->GetBinContent(binxUp) && arrHistE[arrPIDs[i]]->GetBinContent(binxLow)) {
                        dvy_dx = (arrHistVy_1[arrPIDs[i]]->GetBinContent(bin) - arrHistVy_1[arrPIDs[i]]->GetBinContent(binxLow))/dx;
                    }
                    else if(arrHistE[arrPIDs[i]]->GetBinContent(binxUp) && arrHistE[arrPIDs[i]]->GetBinContent(binxLow)) {
                        dvy_dx = (arrHistVy_1[arrPIDs[i]]->GetBinContent(binxUp) - arrHistVy_1[arrPIDs[i]]->GetBinContent(binxLow))/(2*dx);
                    }

                    // dvx_dy
                    if(arrHistE[arrPIDs[i]]->GetBinContent(binyUp) && !arrHistE[arrPIDs[i]]->GetBinContent(binyLow)) {
                        dvy_dx = (arrHistVx_1[arrPIDs[i]]->GetBinContent(binyUp) - arrHistVx_1[arrPIDs[i]]->GetBinContent(bin))/dy;
                    }
                    else if (!arrHistE[arrPIDs[i]]->GetBinContent(binyUp) && arrHistE[arrPIDs[i]]->GetBinContent(binyLow)) {
                        dvy_dx = (arrHistVx_1[arrPIDs[i]]->GetBinContent(bin) - arrHistVx_1[arrPIDs[i]]->GetBinContent(binyLow))/dy;
                    }
                    else if(arrHistE[arrPIDs[i]]->GetBinContent(binyUp) && arrHistE[arrPIDs[i]]->GetBinContent(binyLow)) {
                        dvy_dx = (arrHistVx_1[arrPIDs[i]]->GetBinContent(binyUp) - arrHistVx_1[arrPIDs[i]]->GetBinContent(binyLow))/(2*dy);
                    }

                    arrHistWx_1->SetBinContent(bin, dvz_dy - dvy_dz);
                    arrHistWy_1->SetBinContent(bin, dvx_dz - dvz_dx);
                    arrHistWz_1->SetBinContent(bin, dvy_dx - dvx_dy);
                    arrHistW_1->SetBinContent(bin, TMath::Sqrt(arrHistWx_1->GetBinContent(bin)*arrHistWx_1->GetBinContent(bin) + arrHistWy_1->GetBinContent(bin)*arrHistWy_1->GetBinContent(bin) + arrHistWz_1->GetBinContent(bin)*arrHistWz_1->GetBinContent(bin)));

                }
            }
        }
        arrHistWx_1->Write();
        arrHistWy_1->Write();
        arrHistWz_1->Write();
        arrHistW_1->Write();
    }

/*    // Determing and saving velocity with 2-st way
    //---------------------------------------------------------------------------------------------
    // Consider all histogramms have the same parameters!
    const Int_t nBinsX = arrHistPx->GetNbinsX();
    const Int_t nBinsY = arrHistPy->GetNbinsY();
    const Int_t nBinsZ = arrHistPz->GetNbinsZ();

    // Looping on histograms
    for(Int_t i = 0; i < nBinsX; i++) {
        for(Int_t i = 0; i < nBinsX; i++) {
            for(Int_t i = 0; i < nBinsX; i++) {
                
            }
        }
    }
*/
}






















