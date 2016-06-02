#ifndef _other_functions_h_
#define _other_functions_h_

#include <TLorentzVector.h>
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

Int_t npf, pf_id[5000];

int firstTrackIndex(int jetindex) {
    // Find index of the first track in this jet
    int result = 0;
    for (result = 0; result < npf; ++result) {
        if (pf_id[result]==jetindex) {
            break; // at this point, result is correct
        }
    }
    return result;
}

bool VecSort(TLorentzVector j1, TLorentzVector j2) { return j1.Pt() > j2.Pt(); }
bool sortJetTuple(std::pair<int,float> i, std::pair<int,float> j) { return std::get<1>(i) > std::get<1>(j); }

#endif
