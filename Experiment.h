//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 10 14:36:44 2018 by ROOT version 5.34/30
// from TTree data_tree/Tree fo save
// found on file: tree.root
//////////////////////////////////////////////////////////

#ifndef Experiment_h
#define Experiment_h
#define Max_N 2000
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Experiment {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           n_particles;
   Double_t        time;
   Double_t        r0[Max_N];   //[n_particles]
   Double_t        rx[Max_N];   //[n_particles]
   Double_t        ry[Max_N];   //[n_particles]
   Double_t        rz[Max_N];   //[n_particles]
   Double_t        p0[Max_N];   //[n_particles]
   Double_t        px[Max_N];   //[n_particles]
   Double_t        py[Max_N];   //[n_particles]
   Double_t        pz[Max_N];   //[n_particles]
   Double_t        m[Max_N];   //[n_particles]
   Double_t        psiRp;
   Int_t           pid[Max_N];   //[n_particles]
   Int_t           ityp[Max_N];   //[n_particles]
   Int_t           i3[Max_N];   //[n_particles]
   Int_t           ichg[Max_N];   //[n_particles]
   Int_t           lcl[Max_N];   //[n_particles]
   Int_t           ncl[Max_N];   //[n_particles]
   Int_t           orr[Max_N];   //[n_particles]
   Double_t        imp;

   // List of branches
   TBranch        *b_n_particles;   //!
   TBranch        *b_time;   //!
   TBranch        *b_r0;   //!
   TBranch        *b_rx;   //!
   TBranch        *b_ry;   //!
   TBranch        *b_rz;   //!
   TBranch        *b_p0;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_m;   //!
   TBranch        *b_psiRp;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_ityp;   //!
   TBranch        *b_i3;   //!
   TBranch        *b_ichg;   //!
   TBranch        *b_lcl;   //!
   TBranch        *b_ncl;   //!
   TBranch        *b_orr;   //!
   TBranch        *b_imp;   //!

   Experiment(TTree *tree=0);
   virtual ~Experiment();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Experiment_cxx
Experiment::Experiment(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tree.root");
      }
      f->GetObject("data_tree",tree);

   }
   Init(tree);
}

Experiment::~Experiment()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Experiment::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Experiment::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Experiment::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("n_particles", &n_particles, &b_n_particles);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("r0", r0, &b_r0);
   fChain->SetBranchAddress("rx", rx, &b_rx);
   fChain->SetBranchAddress("ry", ry, &b_ry);
   fChain->SetBranchAddress("rz", rz, &b_rz);
   fChain->SetBranchAddress("p0", p0, &b_p0);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("psiRp", &psiRp, &b_psiRp);
   fChain->SetBranchAddress("pid", pid, &b_pid);
   fChain->SetBranchAddress("ityp", ityp, &b_ityp);
   fChain->SetBranchAddress("i3", i3, &b_i3);
   fChain->SetBranchAddress("ichg", ichg, &b_ichg);
   fChain->SetBranchAddress("lcl", lcl, &b_lcl);
   fChain->SetBranchAddress("ncl", ncl, &b_ncl);
   fChain->SetBranchAddress("orr", orr, &b_orr);
   fChain->SetBranchAddress("imp", &imp, &b_imp);
   Notify();
}

Bool_t Experiment::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Experiment::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Experiment::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Experiment_cxx
