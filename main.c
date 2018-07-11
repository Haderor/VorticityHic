#include "Experiment.C"
#include <iostream>

int main() {
    
     TFile f("/home/verbarius/Work1/tree.root", "READ");
     if (!f.IsOpen()) {
        std::cout << "lol" << std::endl;
     }
     TTree *t = (TTree *)f.Get("data_tree");
     Experiment h(t);
     // h.Show(0);
     h.Loop();
     f.Close();

return 0;
}


