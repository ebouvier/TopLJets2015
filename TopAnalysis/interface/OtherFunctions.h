#ifndef _other_functions_h_
#define _other_functions_h_

#include <TLorentzVector.h>
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include <vector>

Int_t npf, pf_j[5000];

int firstTrackIndex(int jetindex) {
    // Find index of the first track in this jet
    int result = 0;
    for (result = 0; result < npf; ++result) {
        if (pf_j[result]==jetindex) {
            break; // at this point, result is correct
        }
    }
    return result;
}

int firstTrackIndex(int jetindex, std::vector<std::pair<int,float>> *jets) {
    // Find index of the first track in this jet
    int result = 0;
    for (result = 0; result < (int)jets->size(); ++result) {
        if (std::get<0>(jets->at(result))==jetindex) {
            break; // at this point, result is correct
        }
    }
    return result;
}

bool VecSort(TLorentzVector j1, TLorentzVector j2) { return j1.Pt() > j2.Pt(); }
bool sortJetTuple(std::pair<int,float> i, std::pair<int,float> j) { return std::get<1>(i) > std::get<1>(j); }

#endif
