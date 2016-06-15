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

int firstTrackIndex(int jetindex, std::vector<std::tuple<int,float>> *jets) {
    // Find index of the first track in this jet
    int result = 0;
    for (result = 0; result < (int)jets->size(); ++result) {
        if (std::get<0>(jets->at(result))==jetindex) {
            break; // at this point, result is correct
        }
    }
    return result;
}

int firstTrackIndex(int jetindex, std::vector<std::tuple<int,float,float>> *jets) {
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
//bool sortJetTuple(std::tuple<int,float> i, std::tuple<int,float> j) { return std::get<0>(i) > std::get<0>(j); }
bool sortJetTuple(std::tuple<int,float> i, std::tuple<int,float> j) { return std::get<0>(i) < std::get<0>(j) || ( std::get<0>(i) == std::get<0>(j) && std::get<1>(i) > std::get<1>(j) ); }
bool sortJetCSVTuple(std::tuple<int,float,float> i, std::tuple<int,float,float> j) { return std::get<2>(i) > std::get<2>(j); }

#endif
