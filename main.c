#include "Experiment.c"
#include <iostream>

int main() {
    
     TFile f("/home/verbarius/Work1/tree.root", "READ");
     if (!f.IsOpen()) {
        std::cout << "File was not opened!" << std::endl;
     }
     TTree *t = (TTree *)f.Get("data_tree");
     Experiment h(t);
     //h.Show(12);
     h.Loop();
     f.Close();

return 0;
}


