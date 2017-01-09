#ifndef _topmass_h_
#define _topmass_h_

#include "TH1.h"
#include "TString.h"
#include "TFile.h"

#include <map>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

void RunTopMass(TString filename,
    TString outname,
    Int_t chargeSelection, 
    TH1F *normH, 
    Bool_t runSysts,
    TString era,
    Bool_t debug);
#endif
