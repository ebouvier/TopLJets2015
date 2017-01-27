#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/TOPMass.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "TopLJets2015/TopAnalysis/interface/OtherFunctions.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;

//
/*
   Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
   {
   return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
   }


*/
//
void RunTopMass(TString filename,
    TString outname,
    Int_t chargeSelection, 
    TH1F *normH, 
    Bool_t runSysts,
    TString era,
    Bool_t debug=false)
{
  if (debug) cout << "in RunTop" << endl;

  bool isTTbar( filename.Contains("_TTJets") );

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  if (ev.isData) runSysts=false;
  bool requireEletriggerOnly(false);
  if (ev.isData && filename.Contains("SingleElectron")) requireEletriggerOnly=true;
  bool requireMutriggerOnly(false);
  if (ev.isData && filename.Contains("SingleMuon"))     requireMutriggerOnly=true;
  bool requireEETriggers(false);
  if (ev.isData && filename.Contains("DoubleEG"))       requireEETriggers=true;
  bool requireMMTriggers(false);
  if (ev.isData && filename.Contains("DoubleMuon"))     requireMMTriggers=true;
  bool requireEMTriggers(false);
  if (ev.isData && filename.Contains("MuonEG"))         requireEMTriggers=true;

  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;

  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if (!ev.isData)
  {
    if (debug) cout << "Starting loading pileup weight" << endl;
    TString puWgtUrl(era+"/pileupWgts.root");
    gSystem->ExpandPathName(puWgtUrl);
    TFile *fIn=TFile::Open(puWgtUrl);
    for (size_t i=0; i<3; i++)
    {
      TString grName("pu_nom");
      if (i==1) grName="pu_down";
      if (i==2) grName="pu_up";
      TGraph *puData=(TGraph *)fIn->Get(grName);
      Float_t totalData=puData->Integral();
      TH1 *tmp=(TH1 *)puTrue->Clone("tmp");
      for (Int_t xbin=1; xbin<=tmp->GetXaxis()->GetNbins(); xbin++)
      {
        Float_t yexp=puTrue->GetBinContent(xbin);
        Double_t xobs,yobs;
        puData->GetPoint(xbin-1,xobs,yobs);
        tmp->SetBinContent(xbin, yexp>0 ? yobs/(totalData*yexp) : 0. );
      }
      TGraph *gr=new TGraph(tmp);
      grName.ReplaceAll("pu","puwgts");
      gr->SetName(grName);
      puWgtGr.push_back( gr );
      tmp->Delete();
    }
    TH1 *tmp=(TH1 *)puTrue->Clone("tmp");
    TGraph *gr=new TGraph(tmp);
    gr->SetName("puwgts_nom");
    puWgtGr.push_back( gr );
    tmp->Delete();
  }
  if (debug) cout << "loading pileup weight DONE" << endl;

  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);


  //B-TAG CALIBRATION
  TString btagUncUrl(era+"/btagSFactors.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  TString btagEffExpUrl(era+"/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  if (!ev.isData)
  {
    BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
    sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
    sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down") ); 
    sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "up") );

    sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );
    sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "down") ); 
    sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "up") );

    TFile *beffIn=TFile::Open(btagEffExpUrl);
    expBtagEffPy8["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
    expBtagEffPy8["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
    expBtagEffPy8["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
    beffIn->Close();

    TString btagExpPostFix("");
    if (isTTbar)
    {
      if (filename.Contains("_herwig")) btagExpPostFix="_herwig";
      if (filename.Contains("_scaleup")) btagExpPostFix="_scaleup";
      if (filename.Contains("_scaledown")) btagExpPostFix="_scaledown";
    }
    btagEffExpUrl.ReplaceAll(".root",btagExpPostFix+".root");
    beffIn=TFile::Open(btagEffExpUrl);
    expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
    expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
    expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
    beffIn->Close();

  }

  //jet energy uncertainties
  TString jecUncUrl(era+"/jecUncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);

  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::vector<TString> lfsVec = {"_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  std::vector<TString> cutVec = {"_topsel", "_jpsicand", "_d0cand"};
  std::vector<TString> wgtVec = {"", "_no_weight"};

  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);

  /* take too much space when running over all samples
  allPlots["norm_topsel"] = new TH1F("norm_topsel",";n#sigmaL/N_{gen}; Events / 10^{-13}", 6e7,0., 6e-6); 
  allPlots["wgt_topsel"] = new TH1F("wgt_topsel",";Weight;Events / 10^{-13}", 1e8,0., 1e-5);
  */
  allPlots["sf_topsel"] = new TH1F("sf_topsel",";Scale factor; 0.01", 260,0., 2.6);

  for (int k = 0; k < (int)wgtVec.size(); k++) {
    TString weight(wgtVec[k]);

    allPlots["yields_all_check"+weight] = new TH1F("yields_all_check"+weight,";Event yields;Events", 7, 0., 7.);
    allPlots["yields_all_check"+weight]->SetOption("bar");
    allPlots["yields_all_check"+weight]->SetBarWidth(0.75);
    allPlots["yields_all_check"+weight]->SetBarOffset(0.125);
    allPlots["yields_all_check"+weight]->GetXaxis()->SetBinLabel(1, "Starting"); 
    allPlots["yields_all_check"+weight]->GetXaxis()->SetBinLabel(2, "Trigger"); 
    allPlots["yields_all_check"+weight]->GetXaxis()->SetBinLabel(3, "Lepton(s)");
    allPlots["yields_all_check"+weight]->GetXaxis()->SetBinLabel(4, "Jets");
    allPlots["yields_all_check"+weight]->GetXaxis()->SetBinLabel(5, "b-tagging");
    allPlots["yields_all_check"+weight]->GetXaxis()->SetBinLabel(6, "MET");
    allPlots["yields_all_check"+weight]->GetXaxis()->SetBinLabel(7, "DY vetos");

    allPlots["yields_m_check"+weight] = new TH1F("yields_m_check"+weight,";Event yields;Events", 7, 0., 7.);
    allPlots["yields_m_check"+weight]->SetOption("bar");
    allPlots["yields_m_check"+weight]->SetBarWidth(0.75);
    allPlots["yields_m_check"+weight]->SetBarOffset(0.125);
    allPlots["yields_m_check"+weight]->GetXaxis()->SetBinLabel(1, "Starting");
    allPlots["yields_m_check"+weight]->GetXaxis()->SetBinLabel(2, "Trigger");
    allPlots["yields_m_check"+weight]->GetXaxis()->SetBinLabel(3, "Exactly 1 tight lepton");
    allPlots["yields_m_check"+weight]->GetXaxis()->SetBinLabel(4, "... which is a #mu");
    allPlots["yields_m_check"+weight]->GetXaxis()->SetBinLabel(5, "No additional veto lepton");
    allPlots["yields_m_check"+weight]->GetXaxis()->SetBinLabel(6, "At least 4 jets with p_{T} > 30 GeV");
    allPlots["yields_m_check"+weight]->GetXaxis()->SetBinLabel(7, "At least 1 b jet");

    allPlots["yields_e_check"+weight] = new TH1F("yields_e_check"+weight,";Event yields;Events", 7, 0., 7.);
    allPlots["yields_e_check"+weight]->SetOption("bar");
    allPlots["yields_e_check"+weight]->SetBarWidth(0.75);
    allPlots["yields_e_check"+weight]->SetBarOffset(0.125);
    allPlots["yields_e_check"+weight]->GetXaxis()->SetBinLabel(1, "Starting");
    allPlots["yields_e_check"+weight]->GetXaxis()->SetBinLabel(2, "Trigger");
    allPlots["yields_e_check"+weight]->GetXaxis()->SetBinLabel(3, "Exactly 1 tight lepton");
    allPlots["yields_e_check"+weight]->GetXaxis()->SetBinLabel(4, "... which is an e");
    allPlots["yields_e_check"+weight]->GetXaxis()->SetBinLabel(5, "No additional veto lepton");
    allPlots["yields_e_check"+weight]->GetXaxis()->SetBinLabel(6, "At least 4 jets with p_{T} > 30 GeV");
    allPlots["yields_e_check"+weight]->GetXaxis()->SetBinLabel(7, "At least 1 b jet");

    allPlots["yields_mm_check"+weight] = new TH1F("yields_mm_check"+weight,";Event yields;Events", 10, 0., 10.);
    allPlots["yields_mm_check"+weight]->SetOption("bar");
    allPlots["yields_mm_check"+weight]->SetBarWidth(0.75);
    allPlots["yields_mm_check"+weight]->SetBarOffset(0.125);
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(1, "Starting");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(2, "Trigger");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(3, "Exactly 2 medium leptons");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(4, "... which are 2 #mu");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(5, "... oppositely charged");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(6, "At least 2 jets with p_{T} > 30 GeV");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(7, "At least 1 b jet");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(8, "MET > 40 GeV");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(9, "Z veto");
    allPlots["yields_mm_check"+weight]->GetXaxis()->SetBinLabel(10, "Low mass resonances");

    allPlots["yields_em_check"+weight] = new TH1F("yields_em_check"+weight,";Event yields;Events", 7, 0., 7.);
    allPlots["yields_em_check"+weight]->SetOption("bar");
    allPlots["yields_em_check"+weight]->SetBarWidth(0.75);
    allPlots["yields_em_check"+weight]->SetBarOffset(0.125);
    allPlots["yields_em_check"+weight]->GetXaxis()->SetBinLabel(1, "Starting");
    allPlots["yields_em_check"+weight]->GetXaxis()->SetBinLabel(2, "Trigger");
    allPlots["yields_em_check"+weight]->GetXaxis()->SetBinLabel(3, "Exactly 2 medium leptons");
    allPlots["yields_em_check"+weight]->GetXaxis()->SetBinLabel(4, "... which are e#mu");
    allPlots["yields_em_check"+weight]->GetXaxis()->SetBinLabel(5, "... oppositely charged");
    allPlots["yields_em_check"+weight]->GetXaxis()->SetBinLabel(6, "At least 2 jets with p_{T} > 30 GeV");
    allPlots["yields_em_check"+weight]->GetXaxis()->SetBinLabel(7, "At least 1 b jet");

    allPlots["yields_ee_check"+weight] = new TH1F("yields_ee_check"+weight,";Event yields;Events", 10, 0., 10.);
    allPlots["yields_ee_check"+weight]->SetOption("bar");
    allPlots["yields_ee_check"+weight]->SetBarWidth(0.75);
    allPlots["yields_ee_check"+weight]->SetBarOffset(0.125);
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(1, "Starting");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(2, "Trigger");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(3, "Exactly 2 medium leptons");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(4, "... which are 2 e");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(5, "... oppositely charged");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(6, "At least 2 jets with p_{T} > 30 GeV");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(7, "At least 1 b jet");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(8, "MET > 40 GeV");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(9, "Z veto");
    allPlots["yields_ee_check"+weight]->GetXaxis()->SetBinLabel(10, "Low mass resonances");

    for (int i = 0; i < (int)lfsVec.size(); i++) {
      TString tag(lfsVec[i]);

      allPlots["lp_n"+tag+"_check"+weight] = new TH1F("lp_n"+tag+"_check"+weight,";Lepton multiplicity;Events", 6, 0., 6.);
      allPlots["lp_pt"+tag+"_check"+weight] = new TH1F("lp_pt"+tag+"_check"+weight,";Lepton p_{T} (GeV);Events / (5 GeV)", 40, 0., 200.);
      allPlots["lp_eta"+tag+"_check"+weight] = new TH1F("lp_eta"+tag+"_check"+weight,";Lepton #eta;Events / 0.1", 50, -2.5, 2.5);
      allPlots["lp_phi"+tag+"_check"+weight] = new TH1F("lp_phi"+tag+"_check"+weight,";Lepton #phi;Events / 0.1", 64, -3.2, 3.2);
      allPlots["lp_iso"+tag+"_check"+weight] = new TH1F("lp_iso"+tag+"_check"+weight,";Lepton isolation;Events / 0.025", 40, 0., 1.);
      allPlots["lp_id"+tag+"_check"+weight] = new TH1F("lp_id"+tag+"_check"+weight,";Lepton ID;Events", 29, -14., 15.);

      allPlots["vlp_n"+tag+"_check"+weight] = new TH1F("vlp_n"+tag+"_check"+weight,";Veto lepton multiplicity;Events", 5, 0., 5.);
      allPlots["vlp_pt"+tag+"_check"+weight] = new TH1F("vlp_pt"+tag+"_check"+weight,";Veto lepton p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
      allPlots["vlp_eta"+tag+"_check"+weight] = new TH1F("vlp_eta"+tag+"_check"+weight,";Veto lepton #eta;Events / 0.2", 25, -2.5, 2.5);
      allPlots["vlp_phi"+tag+"_check"+weight] = new TH1F("vlp_phi"+tag+"_check"+weight,";Veto lepton #phi;Events / 0.2", 32, -3.2, 3.2);
      allPlots["vlp_iso"+tag+"_check"+weight] = new TH1F("vlp_iso"+tag+"_check"+weight,";Veto lepton isolation;Events / 0.01", 30, 0., 0.3);
      allPlots["vlp_id"+tag+"_check"+weight] = new TH1F("vlp_id"+tag+"_check"+weight,";Veto lepton ID;Events", 29, -14., 15.);

      allPlots["mlp_n"+tag+"_check"+weight] = new TH1F("mlp_n"+tag+"_check"+weight,";Medium lepton multiplicity;Events", 4, 0., 4.);
      allPlots["mlp_pt"+tag+"_check"+weight] = new TH1F("mlp_pt"+tag+"_check"+weight,";Medium lepton p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
      allPlots["mlp_eta"+tag+"_check"+weight] = new TH1F("mlp_eta"+tag+"_check"+weight,";Medium lepton #eta;Events / 0.2", 25, -2.5, 2.5);
      allPlots["mlp_phi"+tag+"_check"+weight] = new TH1F("mlp_phi"+tag+"_check"+weight,";Medium lepton #phi;Events / 0.2", 32, -3.2, 3.2);
      allPlots["mlp_iso"+tag+"_check"+weight] = new TH1F("mlp_iso"+tag+"_check"+weight,";Medium lepton isolation;Events / 0.01", 20, 0., 0.2);
      allPlots["mlp_id"+tag+"_check"+weight] = new TH1F("mlp_id"+tag+"_check"+weight,";Medium lepton ID;Events", 29, -14., 15.);

      allPlots["tlp_n"+tag+"_check"+weight] = new TH1F("tlp_n"+tag+"_check"+weight,";Tight lepton multiplicity;Events", 4, 0., 4.);
      allPlots["tlp_pt"+tag+"_check"+weight] = new TH1F("tlp_pt"+tag+"_check"+weight,";Tight lepton p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
      allPlots["tlp_eta"+tag+"_check"+weight] = new TH1F("tlp_eta"+tag+"_check"+weight,";Tight lepton #eta;Events / 0.2", 25, -2.5, 2.5);
      allPlots["tlp_phi"+tag+"_check"+weight] = new TH1F("tlp_phi"+tag+"_check"+weight,";Tight lepton #phi;Events / 0.2", 32, -3.2, 3.2);
      allPlots["tlp_iso"+tag+"_check"+weight] = new TH1F("tlp_iso"+tag+"_check"+weight,";Tight lepton isolation;Events / 0.01", 20, 0., 0.2);
      allPlots["tlp_id"+tag+"_check"+weight] = new TH1F("tlp_id"+tag+"_check"+weight,";Tight lepton ID;Events", 29, -14., 15.);

      allPlots["slp_n"+tag+"_check"+weight] = new TH1F("slp_n"+tag+"_check"+weight,";Selected lepton multiplicity;Events", 4, 0., 4.);
      allPlots["slp_pt"+tag+"_check"+weight] = new TH1F("slp_pt"+tag+"_check"+weight,";Selected lepton p_{T} (GeV);Events / (10 GeV)", 40, 0., 200.);
      allPlots["slp_eta"+tag+"_check"+weight] = new TH1F("slp_eta"+tag+"_check"+weight,";Selected lepton #eta;Events / 0.2", 25, -2.5, 2.5);
      allPlots["slp_phi"+tag+"_check"+weight] = new TH1F("slp_phi"+tag+"_check"+weight,";Selected lepton #phi;Events / 0.2", 32, -3.2, 3.2);

      allPlots["j_n"+tag+"_check"+weight] = new TH1F("j_n"+tag+"_check"+weight,";Jet multiplicity;Events", 10, 0., 10.);
      allPlots["j_pt"+tag+"_check"+weight] = new TH1F("j_pt"+tag+"_check"+weight,";Jet p_{T} (GeV);Events / (10 GeV)", 30, 0., 300.); 
      allPlots["j_eta"+tag+"_check"+weight] = new TH1F("j_eta"+tag+"_check"+weight,";Jet #eta;Events / 0.2", 25, -2.5, 2.5);
      allPlots["j_phi"+tag+"_check"+weight] = new TH1F("j_phi"+tag+"_check"+weight,";Jet #phi;Events / 0.2", 32, -3.2, 3.2);
      allPlots["j_csv"+tag+"_check"+weight] = new TH1F("j_csv"+tag+"_check"+weight,";Jet CSV;Events / 0.05", 20, 0., 1.);

      allPlots["sj_n"+tag+"_check"+weight] = new TH1F("sj_n"+tag+"_check"+weight,";Selected jet multiplicity;Events", 10, 0., 10.);
      allPlots["sj40_n"+tag+"_check"+weight] = new TH1F("sj_n"+tag+"_check"+weight,";Selected jet multiplicity (p_{T} > 40 GeV);Events", 10, 0., 10.);
      allPlots["sj_pt"+tag+"_check"+weight] = new TH1F("sj_pt"+tag+"_check"+weight,";Selected jet p_{T} (GeV);Events / (10 GeV)", 30, 0., 300.); 
      allPlots["sj_eta"+tag+"_check"+weight] = new TH1F("sj_eta"+tag+"_check"+weight,";Selected jet #eta;Events / 0.2", 25, -2.5, 2.5);
      allPlots["sj_phi"+tag+"_check"+weight] = new TH1F("sj_phi"+tag+"_check"+weight,";Selected jet #phi;Events / 0.2", 32, -3.2, 3.2);
      allPlots["sj_csv"+tag+"_check"+weight] = new TH1F("sj_csv"+tag+"_check"+weight,";Selected jet CSV;Events / 0.05", 20, 0., 1.);

      allPlots["bj_n"+tag+"_check"+weight] = new TH1F("bj_n"+tag+"_check"+weight,";b jet multiplicity (CSV > 0.8);Events", 5, 0., 5.);
      allPlots["bj_pt"+tag+"_check"+weight] = new TH1F("bj_pt"+tag+"_check"+weight,";b jet p_{T} (GeV);Events / (10 GeV)", 30, 0., 300.);
      allPlots["bj_eta"+tag+"_check"+weight] = new TH1F("bj_eta"+tag+"_check"+weight,";b jet #eta;Events / 0.2", 25, -2.5, 2.5);
      allPlots["bj_phi"+tag+"_check"+weight] = new TH1F("bj_phi"+tag+"_check"+weight,";b jet #phi;Events / 0.2", 32, -3.2, 3.2);
      allPlots["bj_csv"+tag+"_check"+weight] = new TH1F("bj_csv"+tag+"_check"+weight,";b jet CSV;Events / 0.01", 20, 0.8, 1.);

      allPlots["sbj_n"+tag+"_check"+weight] = new TH1F("sbj_n"+tag+"_check"+weight,";Selected b jet multiplicity (CSV > 0.8);Events", 4, 1., 5.);
      allPlots["sbj_pt"+tag+"_check"+weight] = new TH1F("sbj_pt"+tag+"_check"+weight,";Selected b jet p_{T} (GeV);Events / (10 GeV)", 30, 0., 300.);
      allPlots["sbj_eta"+tag+"_check"+weight] = new TH1F("sbj_eta"+tag+"_check"+weight,";Selected b jet #eta;Events / 0.2", 25, -2.5, 2.5);
      allPlots["sbj_phi"+tag+"_check"+weight] = new TH1F("sbj_phi"+tag+"_check"+weight,";Selected b jet #phi;Events / 0.2", 32, -3.2, 3.2);
      allPlots["sbj_csv"+tag+"_check"+weight] = new TH1F("sbj_csv"+tag+"_check"+weight,";Selected b jet CSV;Events / 0.01", 20, 0.8, 1.);

      allPlots["lp_jets_dR"+tag+"_check"+weight] = new TH1F("lp_jets_dR"+tag+"_check"+weight,";#DeltaR(l, jets);Events / 0.05", 80, 0., 4.);

      allPlots["met"+tag+"_check"+weight] = new TH1F("met"+tag+"_check"+weight,";MET (GeV);Events / (10 GeV)", 20, 0., 200.);

      allPlots["z_mass"+tag+"_check"+weight] = new TH1F("z_mass"+tag+"_check"+weight,";Z boson candidate mass (GeV);Events / (1 GeV)",30, 76., 106.);

      allPlots["jpsi_mass"+tag+"_check"+weight] = new TH1F("jpsi_mass"+tag+"_check"+weight,";M_{#mu^{+}#mu^{-}} (GeV);Events / (25 MeV)", 40, 2.6, 3.6);

      allPlots["d0_mass"+tag+"_check"+weight] = new TH1F("d0_mass"+tag+"_check"+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (15 MeV)", 20, 1.7, 2.0);

      for (int j = 0; j < (int)cutVec.size(); j++) {
        TString cut(cutVec[j]);

        if (tag.Contains("all")) {
          allPlots["decay"+tag+cut+weight] = new TH1F("decay"+tag+cut+weight, ";Decay channel;Events", 5, 0., 5.);
          allPlots["decay"+tag+cut+weight]->SetOption("bar");
          allPlots["decay"+tag+cut+weight]->SetBarWidth(0.75);
          allPlots["decay"+tag+cut+weight]->SetBarOffset(0.125);
          allPlots["decay"+tag+cut+weight]->GetXaxis()->SetBinLabel(1, "e + jets");
          allPlots["decay"+tag+cut+weight]->GetXaxis()->SetBinLabel(2, "#mu + jets");
          allPlots["decay"+tag+cut+weight]->GetXaxis()->SetBinLabel(3, "ee + jets");
          allPlots["decay"+tag+cut+weight]->GetXaxis()->SetBinLabel(4, "e#mu + jets");
          allPlots["decay"+tag+cut+weight]->GetXaxis()->SetBinLabel(5, "#mu#mu + jets");
          allPlots["decay"+tag+cut+weight]->GetYaxis()->SetTitle("Events");
        }

        allPlots["lp_n"+tag+cut+weight] = new TH1F("lp_n"+tag+cut+weight,";Lepton multiplicity;Events", 4, 0., 4.);
        allPlots["lp1_pt"+tag+cut+weight] = new TH1F("lp1_pt"+tag+cut+weight,";Leading lepton p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
        allPlots["lp1_eta"+tag+cut+weight] = new TH1F("lp1_eta"+tag+cut+weight,";Leading lepton #eta;Events / 0.2", 25, -2.5, 2.5);
        allPlots["lp1_phi"+tag+cut+weight] = new TH1F("lp1_phi"+tag+cut+weight,";Leading lepton #phi;Events / 0.2", 32, -3.2, 3.2);
        allPlots["lp2_pt"+tag+cut+weight] = new TH1F("lp2_pt"+tag+cut+weight,";Subleading lepton p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
        allPlots["lp2_eta"+tag+cut+weight] = new TH1F("lp2_eta"+tag+cut+weight,";Subleading lepton #eta;Events / 0.2", 25, -2.5, 2.5);
        allPlots["lp2_phi"+tag+cut+weight] = new TH1F("lp2_phi"+tag+cut+weight,";Subleading lepton #phi;Events / 0.2", 32, -3.2, 3.2);

        allPlots["met"+tag+cut+weight] = new TH1F("met"+tag+cut+weight,";MET (GeV);Events / (10 GeV)", 20, 0., 200.);

        allPlots["dilp_pt"+tag+cut+weight] = new TH1F("dilp_pt"+tag+cut+weight,";Dilepton pair p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
        allPlots["dilp_mass"+tag+cut+weight] = new TH1F("dilp_mass"+tag+cut+weight,";Dilepton pair mass (GeV);Events / (10 GeV)", 20, 0., 200.);
        allPlots["dilp_charge"+tag+cut+weight] = new TH1F("dilp_charge"+tag+cut+weight,";Dilepton pair charge;Events", 5, -2., 2.);

        allPlots["j_n"+tag+cut+weight] = new TH1F("j_n"+tag+cut+weight,";Jet multiplicity;Events", 10, 0., 10.);
        allPlots["j_pt"+tag+cut+weight] = new TH1F("j_pt"+tag+cut+weight,";p_{T}(jet) (GeV);Events / (10 GeV)", 30, 0., 300.); 
        allPlots["j_eta"+tag+cut+weight] = new TH1F("j_eta"+tag+cut+weight,";#eta(jet);Events / 0.2", 25, -2.5, 2.5);
        allPlots["j_phi"+tag+cut+weight] = new TH1F("j_phi"+tag+cut+weight,";#phi(jet);Events / 0.2", 32, -3.2, 3.2);
        allPlots["j_csv"+tag+cut+weight] = new TH1F("j_csv"+tag+cut+weight,";CSV(jet);Events / 0.05", 20, 0., 1.);
        allPlots["j_nch"+tag+cut+weight] = new TH1F("j_nch"+tag+cut+weight,";Track multiplicity per jet; Events", 50, 0., 50.);
        allPlots["j1_pt"+tag+cut+weight] = new TH1F("j1_pt"+tag+cut+weight,";Leading jet p_{T} (GeV);Events / (10 GeV)", 30, 0., 300.); 
        allPlots["j1_eta"+tag+cut+weight] = new TH1F("j1_eta"+tag+cut+weight,";Leading jet #eta;Events / 0.2", 25, -2.5, 2.5);
        allPlots["j1_phi"+tag+cut+weight] = new TH1F("j1_phi"+tag+cut+weight,";Leading jet #phi;Events / 0.2", 32, -3.2, 3.2);
        allPlots["j1_csv"+tag+cut+weight] = new TH1F("j1_csv"+tag+cut+weight,";Leading jet CSV discriminant;Events / 0.05", 20, 0., 1.);
        allPlots["j1_nch"+tag+cut+weight] = new TH1F("j1_nch"+tag+cut+weight,";Track multiplicity per leading jet; Events", 50, 0., 50.);
        allPlots["j2_pt"+tag+cut+weight] = new TH1F("j2_pt"+tag+cut+weight,";Subleading jet p_{T} (GeV);Events / (10 GeV)", 30, 0., 300.); 
        allPlots["j2_eta"+tag+cut+weight] = new TH1F("j2_eta"+tag+cut+weight,";Subleading jet #eta;Events / 0.2", 25, -2.5, 2.5);
        allPlots["j2_phi"+tag+cut+weight] = new TH1F("j2_phi"+tag+cut+weight,";Subleading jet #phi;Events / 0.2", 32, -3.2, 3.2);
        allPlots["j2_csv"+tag+cut+weight] = new TH1F("j2_csv"+tag+cut+weight,";Subleading jet CSV discriminant;Events / 0.05", 20, 0., 1.);
        allPlots["j2_nch"+tag+cut+weight] = new TH1F("j2_nch"+tag+cut+weight,";Track multiplicity per subleading jet; Events", 50, 0., 50.);

        allPlots["lj_n"+tag+cut+weight] = new TH1F("lj_n"+tag+cut+weight,";Light jet multiplicity;Events" ,10,0,10.);
        allPlots["lj_pt"+tag+cut+weight] = new TH1F("lj_pt"+tag+cut+weight,";p_{T}(light jet) (GeV);Events / (10 GeV)", 30, 0., 300.);
        allPlots["lj_eta"+tag+cut+weight] = new TH1F("lj_eta"+tag+cut+weight,";#eta(light jet);Events / 0.2", 25, -2.5, 2.5);
        allPlots["lj_phi"+tag+cut+weight] = new TH1F("lj_phi"+tag+cut+weight,";#phi(light jet);Events / 0.2", 32, -3.2, 3.2);
        allPlots["lj_nch"+tag+cut+weight] = new TH1F("lj_nch"+tag+cut+weight,";Track multiplicity per light jet; Events", 50, 0., 50.);
        allPlots["lj1_pt"+tag+cut+weight] = new TH1F("lj1_pt"+tag+cut+weight,";Leading light jet p_{T} (GeV);Events / (10 GeV)", 30, 0., 300.);
        allPlots["lj1_eta"+tag+cut+weight] = new TH1F("lj1_eta"+tag+cut+weight,";Leading light jet #eta;Events / 0.2", 25, -2.5, 2.5);
        allPlots["lj1_phi"+tag+cut+weight] = new TH1F("lj1_phi"+tag+cut+weight,";Leading light jet #phi;Events / 0.2", 32, -3.2, 3.2);
        allPlots["lj1_nch"+tag+cut+weight] = new TH1F("lj1_nch"+tag+cut+weight,";Track multiplicity per leading light jet; Events", 50, 0., 50.);

        allPlots["bj_n"+tag+cut+weight] = new TH1F("bj_n"+tag+cut+weight,";b-tagged jet multiplicity (CSV > 0.8);Events", 5, 0., 5.);
        allPlots["bj_pt"+tag+cut+weight] = new TH1F("bj_pt"+tag+cut+weight,";p_{T}(b jet) (GeV);Events / (10 GeV)", 30, 0., 300.);
        allPlots["bj_eta"+tag+cut+weight] = new TH1F("bj_eta"+tag+cut+weight,";#eta(b jet);Events / 0.2", 25, -2.5, 2.5);
        allPlots["bj_phi"+tag+cut+weight] = new TH1F("bj_phi"+tag+cut+weight,";#phi(b jet);Events / 0.2", 32, -3.2, 3.2);
        allPlots["bj_csv"+tag+cut+weight] = new TH1F("bj_csv"+tag+cut+weight,";CSV(b jet);Events / 0.01", 20, 0.8, 1.);
        allPlots["bj_nch"+tag+cut+weight] = new TH1F("bj_nch"+tag+cut+weight,";Track multiplicity per b jet; Events", 50, 0., 50.);
        allPlots["bj1_pt"+tag+cut+weight] = new TH1F("bj1_pt"+tag+cut+weight,";Leading b jet p_{T} (GeV);Events / (10 GeV)", 30, 0., 300.);
        allPlots["bj1_eta"+tag+cut+weight] = new TH1F("bj1_eta"+tag+cut+weight,";Leading b jet #eta;Events / 0.2", 25, -2.5, 2.5);
        allPlots["bj1_phi"+tag+cut+weight] = new TH1F("bj1_phi"+tag+cut+weight,";Leading b jet #phi;Events / 0.2", 32, -3.2, 3.2);
        allPlots["bj1_csv"+tag+cut+weight] = new TH1F("bj1_csv"+tag+cut+weight,";Leading b jet CSV discriminant;Events / 0.01", 20, 0.8, 1.);
        allPlots["bj1_nch"+tag+cut+weight] = new TH1F("bj1_nch"+tag+cut+weight,";Track multiplicity per leading b jet; Events", 50, 0., 50.);

        if (cut.Contains("jpsi")) {
          allPlots["jpsi_lp_dR"+tag+cut+weight] = new TH1F("jpsi_lp_dR"+tag+cut+weight,";#DeltaR(#mu^{+}#mu^{-}, l);Events / 0.4", 20, 0., 5.);
          allPlots["jpsi_lp_mass"+tag+cut+weight] = new TH1F("jpsi_lp_mass"+tag+cut+weight,";M_{#mu^{+}#mu^{-}+l} (GeV);Events / (10 GeV)", 25, 0., 250.);
          allPlots["jpsi_n"+tag+cut+weight] = new TH1F("jpsi_n"+tag+cut+weight,";J/#psi multiplicity;Events", 3, 0., 3.);
          allPlots["jpsi_mass"+tag+cut+weight] = new TH1F("jpsi_mass"+tag+cut+weight,";M_{#mu^{+}#mu^{-}} (GeV);Events / (20 MeV)", 20, 2.9, 3.3);
          allPlots["jpsi_ptfrac"+tag+cut+weight] = new TH1F("jpsi_ptfrac"+tag+cut+weight,";p_{T}(#mu^{+}#mu^{-})/p_{T}(jet);Events / 0.04", 25, 0., 1.);
          allPlots["jpsi_pt"+tag+cut+weight] = new TH1F("jpsi_pt"+tag+cut+weight,";p_{T}(#mu^{+}#mu^{-}) (GeV);Events / (5 GeV)", 30, 0., 150.);
          allPlots["jpsi_eta"+tag+cut+weight] = new TH1F("jpsi_eta"+tag+cut+weight,";#eta(#mu^{+}#mu^{-});Events / 0.2", 25, -2.5, 2.5);
          allPlots["jpsi_phi"+tag+cut+weight] = new TH1F("jpsi_phi"+tag+cut+weight,";#phi(#mu^{+}#mu^{-});Events / 0.2", 32, -3.2, 3.2);
          allPlots["mu_jpsi_pt"+tag+cut+weight] = new TH1F("mu_jpsi_pt"+tag+cut+weight,";p_{T}(#mu^{#pm}) (GeV);Events / (5 GeV)", 30, 0., 150.);
          allPlots["mu_jpsi_eta"+tag+cut+weight] = new TH1F("mu_jpsi_eta"+tag+cut+weight,";#eta(#mu^{#pm});Events / 0.2", 25, -2.5, 2.5);
          allPlots["mu_jpsi_phi"+tag+cut+weight] = new TH1F("mu_jpsi_phi"+tag+cut+weight,";#phi(#mu^{#pm});Events / 0.2", 32, -3.2, 3.2);
          allPlots["mu_jpsi_dR"+tag+cut+weight] = new TH1F("mu_jpsi_dR"+tag+cut+weight,";#DeltaR(#mu^{#pm});Events / 0.03", 20, 0., 0.6);
          allPlots["j_jpsi_pt"+tag+cut+weight] = new TH1F("j_jpsi_pt"+tag+cut+weight,";p_{T}(jet with #mu^{+}#mu^{-}) (GeV);Events / (10 GeV)", 30, 0., 300.);
          allPlots["j_jpsi_eta"+tag+cut+weight] = new TH1F("j_jpsi_eta"+tag+cut+weight,";#eta(jet with #mu^{+}#mu^{-});Events / 0.2", 25, -2.5, 2.5);
          allPlots["j_jpsi_phi"+tag+cut+weight] = new TH1F("j_jpsi_phi"+tag+cut+weight,";#phi(jet with #mu^{+}#mu^{-});Events / 0.2", 32, -3.2, 3.2);
          allPlots["j_jpsi_csv"+tag+cut+weight] = new TH1F("j_jpsi_csv"+tag+cut+weight,";CSV(jet with #mu^{+}#mu^{-});Events / 0.05", 20, 0., 1.);
          allPlots["j_jpsi_nch"+tag+cut+weight] = new TH1F("j_jpsi_nch"+tag+cut+weight,";Track multiplicity per jet with #mu^{+}#mu^{-}; Events", 50, 0., 50.);
        }
        if (cut.Contains("d0")) {
          allPlots["d0_n"+tag+cut+weight] = new TH1F("d0_n"+tag+cut+weight,";D^{0} multiplicity;Events", 7, 0., 7.);
          allPlots["d0_mass"+tag+cut+weight] = new TH1F("d0_mass"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (10 MeV)", 11, 1.81, 1.92);
          allPlots["d0_pt"+tag+cut+weight] = new TH1F("d0_pt"+tag+cut+weight,";p_{T}(#kappa^{#pm}#pi^{#mp}) (GeV);Events / (5 GeV)", 30, 0., 150.);
          allPlots["d0_ptfrac"+tag+cut+weight] = new TH1F("d0_ptfrac"+tag+cut+weight,";p_{T}(#kappa^{#pm}#pi^{#mp})/p_{T}(jet);Events / 0.04", 25, 0., 1.);
          allPlots["d0_eta"+tag+cut+weight] = new TH1F("d0_eta"+tag+cut+weight,";#eta(#kappa^{#pm}#pi^{#mp});Events / 0.2", 25, -2.5, 2.5);
          allPlots["d0_phi"+tag+cut+weight] = new TH1F("d0_phi"+tag+cut+weight,";#phi(#kappa^{#pm}#pi^{#mp});Events / 0.2", 32, -3.2, 3.2);
          allPlots["kappa_d0_pt"+tag+cut+weight] = new TH1F("kappa_d0_pt"+tag+cut+weight,";p_{T}(#kappa^{#pm}) (GeV);Events / (5 GeV)", 30, 0., 150.);
          allPlots["kappa_d0_eta"+tag+cut+weight] = new TH1F("kappa_d0_eta"+tag+cut+weight,";#eta(#kappa^{#pm});Events / 0.2", 25, -2.5, 2.5);
          allPlots["kappa_d0_phi"+tag+cut+weight] = new TH1F("kappa_d0_phi"+tag+cut+weight,";#phi(#kappa^{#pm});Events / 0.2", 32, -3.2, 3.2);
          allPlots["pi_d0_pt"+tag+cut+weight] = new TH1F("pi_d0_pt"+tag+cut+weight,";p_{T}(#pi^{#mp}) (GeV);Events / (5 GeV)", 30, 0., 150.);
          allPlots["pi_d0_eta"+tag+cut+weight] = new TH1F("pi_d0_eta"+tag+cut+weight,";#eta(#pi^{#mp});Events / 0.2", 25, -2.5, 2.5);
          allPlots["pi_d0_phi"+tag+cut+weight] = new TH1F("pi_d0_phi"+tag+cut+weight,";#phi(#pi^{#mp});Events / 0.2", 32, -3.2, 3.2);
          allPlots["kappa_pi_d0_dR"+tag+cut+weight] = new TH1F("kappa_pi_d0_dR"+tag+cut+weight,";#DeltaR(#kappa^{#pm}, #pi^{#mp});Events / 0.02", 25, 0., 0.5);
          allPlots["j_d0_pt"+tag+cut+weight] = new TH1F("j_d0_pt"+tag+cut+weight,";p_{T}(jet with #kappa^{#pm}#pi^{#mp}) (GeV);Events / (10 GeV)", 30, 0., 300.);
          allPlots["j_d0_eta"+tag+cut+weight] = new TH1F("j_d0_eta"+tag+cut+weight,";#eta(jet with #kappa^{#pm}#pi^{#mp});Events / 0.2", 25, -2.5, 2.5);
          allPlots["j_d0_phi"+tag+cut+weight] = new TH1F("j_d0_phi"+tag+cut+weight,";#phi(jet with #kappa^{#pm}#pi^{#mp});Events / 0.2", 32, -3.2, 3.2);
          allPlots["j_d0_csv"+tag+cut+weight] = new TH1F("j_d0_csv"+tag+cut+weight,";CSV(jet with #kappa^{#pm}#pi^{#mp});Events / 0.01", 20, 0., 0.8);
          allPlots["j_d0_nch"+tag+cut+weight] = new TH1F("j_d0_nch"+tag+cut+weight,";Track multiplicity per jet with #kappa^{#pm}#pi^{#mp}; Events", 50, 0., 50.);
          allPlots["d0_mass_lp"+tag+cut+weight] = new TH1F("d0_mass_lp"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (10 MeV)", 11, 1.81, 1.92);
          allPlots["d0_mass_mu"+tag+cut+weight] = new TH1F("d0_mass_mu"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (10 MeV)", 11, 1.81, 1.92);
          allPlots["d0_mass_el"+tag+cut+weight] = new TH1F("d0_mass_el"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (10 MeV)", 11, 1.81, 1.92);
          allPlots["ds_mass"+tag+cut+weight] = new TH1F("ds_mass"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}#plus#pi^{#mp}};Events / 10 MeV" ,12,1.95,2.07);
          allPlots["pi_ds_pt"+tag+cut+weight] = new TH1F("pi_ds_pt"+tag+cut+weight,";#pi^{#pm} p_{T} (GeV);Events / (5 GeV)", 10, 0., 50.);
          allPlots["ds_d0_dmass"+tag+cut+weight] = new TH1F("ds_d0_dmass"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}#plus#pi^{#mp}} #minus M_{#kapp^{pm}#pi^{#mp}};Events / (1 MeV)", 30, 0.135, 0.165);
        }
      }
    }
  }


  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  //for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }
  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
  {
    t->GetEntry(iev);
    if (iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

    // nominal event weight for MC (before SF)
    float norm = 1.;
    if (!ev.isData) {
      //float norm( normH ? normH->GetBinContent(1) : 1.0);
      norm =  normH ? normH->GetBinContent(1) : 1.0;
      if (ev.ttbar_nw>0) norm*=ev.ttbar_w[0];
    }

    std::map<TString, int> iCut;
    for(size_t i = 0; i < lfsVec.size(); i++) {
      TString tag(lfsVec[i]);
      iCut[tag] = 0;
      allPlots["yields"+tag+"_check"]->Fill((double)iCut[tag],norm);
      allPlots["yields"+tag+"_check_no_weight"]->Fill((double)iCut[tag],norm);
      ++iCut[tag];
    }

    //check if triggers have fired, see https://indico.cern.ch/event/600194/contributions/2433813/attachments/1397223/2130519/TopTriggers2016_17Jan_v2.pdf
    bool hasEETrigger(((ev.elTrigger>>3)&0x1)!=0); //'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v'
    bool hasMMTrigger(((ev.muTrigger>>4)&0x3)!=0); //'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v'
    bool hasEMTrigger(((ev.elTrigger>>4)&0x3)!=0); //'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_(DZ_)v', 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v, 'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(DZ_)v
    bool hasMuTrigger(((ev.muTrigger>>2)&0x3)!=0); //'HLT_IsoMu24_v', 'HLT_IsoTkMu24_v'
    bool hasEleTrigger((ev.elTrigger & 0x1)!=0);   //'HLT_Ele27_WPTight_Gsf_v'
    if (!ev.isData) {	 
      hasMuTrigger=true;
      hasEleTrigger=true;
      hasEETrigger=true;
      hasMMTrigger=true;
      hasEMTrigger=true;
    } else {
      if (requireMutriggerOnly && !hasMuTrigger) continue;
      if (requireEletriggerOnly && !hasEleTrigger) continue;
      if (requireEETriggers && !hasEETrigger) continue;
      if (requireMMTriggers && !hasMMTrigger) continue;
      if (requireEMTriggers && !hasEMTrigger) continue;
    }
    if (hasMuTrigger && (!ev.isData || requireMutriggerOnly)) {
      allPlots["yields_m_check"]->Fill((double)iCut["_m"],norm);
      allPlots["yields_m_check_no_weight"]->Fill((double)iCut["_m"],norm);
      ++iCut["_m"];
    }
    if (hasEleTrigger && (!ev.isData || requireEletriggerOnly)) {
      allPlots["yields_e_check"]->Fill((double)iCut["_e"],norm);
      allPlots["yields_e_check_no_weight"]->Fill((double)iCut["_e"],norm);
      ++iCut["_e"];
    }
    if (hasMMTrigger && (!ev.isData || requireMMTriggers)) {
      allPlots["yields_mm_check"]->Fill((double)iCut["_mm"],norm);
      allPlots["yields_mm_check_no_weight"]->Fill((double)iCut["_mm"],norm);
      ++iCut["_mm"];
    }
    if (hasEETrigger && (!ev.isData || requireEETriggers)) {
      allPlots["yields_ee_check"]->Fill((double)iCut["_ee"],norm);
      allPlots["yields_ee_check_no_weight"]->Fill((double)iCut["_ee"],norm);
      ++iCut["_ee"];
    }
    if (hasEMTrigger && (!ev.isData || requireEMTriggers)) {
      allPlots["yields_em_check"]->Fill((double)iCut["_em"],norm);
      allPlots["yields_em_check_no_weight"]->Fill((double)iCut["_em"],norm);
      ++iCut["_em"];
    }
    allPlots["yields_all_check"]->Fill((double)iCut["_all"],norm);
    allPlots["yields_all_check_no_weight"]->Fill((double)iCut["_all"],norm);
    ++iCut["_all"];

    //select leptons
    if (debug) cout << "Starting lepton selection -- # all = " << ev.nl << endl;

    std::vector<int> tightLeptons, medLeptons, vetoLeptons;
    if (hasMuTrigger && (!ev.isData || requireMutriggerOnly)) {
      allPlots["lp_n_m_check"]->Fill(ev.nl,norm);
      allPlots["lp_n_m_check_no_weight"]->Fill(ev.nl,norm);
    }
    if (hasEleTrigger && (!ev.isData || requireEletriggerOnly)) {
      allPlots["lp_n_e_check"]->Fill(ev.nl,norm);
      allPlots["lp_n_e_check_no_weight"]->Fill(ev.nl,norm);
    }
    if (hasMMTrigger && (!ev.isData || requireMMTriggers)) {
      allPlots["lp_n_mm_check"]->Fill(ev.nl,norm);
      allPlots["lp_n_mm_check_no_weight"]->Fill(ev.nl,norm);
    }
    if (hasEETrigger && (!ev.isData || requireEETriggers)) {
      allPlots["lp_n_ee_check"]->Fill(ev.nl,norm);
      allPlots["lp_n_ee_check_no_weight"]->Fill(ev.nl,norm);
    }
    if (hasEMTrigger && (!ev.isData || requireEMTriggers)) {
      allPlots["lp_n_em_check"]->Fill(ev.nl,norm);
      allPlots["lp_n_em_check_no_weight"]->Fill(ev.nl,norm);
    }
    allPlots["lp_n_all_check"]->Fill(ev.nl,norm);
    allPlots["lp_n_all_check_no_weight"]->Fill(ev.nl,norm);

    for (int il=0; il<ev.nl; il++) {
      if (hasMuTrigger && (!ev.isData || requireMutriggerOnly)) {
        allPlots["lp_pt_m_check"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_pt_m_check_no_weight"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_eta_m_check"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_eta_m_check_no_weight"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_phi_m_check"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_phi_m_check_no_weight"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_iso_m_check"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_iso_m_check_no_weight"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_id_m_check"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
        allPlots["lp_id_m_check_no_weight"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
      }
      if (hasEleTrigger && (!ev.isData || requireEletriggerOnly)) {
        allPlots["lp_pt_e_check"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_pt_e_check_no_weight"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_eta_e_check"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_eta_e_check_no_weight"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_phi_e_check"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_phi_e_check_no_weight"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_iso_e_check"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_iso_e_check_no_weight"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_id_e_check"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
        allPlots["lp_id_e_check_no_weight"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
      }
      if (hasMMTrigger && (!ev.isData || requireMMTriggers)) {
        allPlots["lp_pt_mm_check"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_pt_mm_check_no_weight"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_eta_mm_check"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_eta_mm_check_no_weight"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_phi_mm_check"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_phi_mm_check_no_weight"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_iso_mm_check"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_iso_mm_check_no_weight"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_id_mm_check"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
        allPlots["lp_id_mm_check_no_weight"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
      }
      if (hasEETrigger && (!ev.isData || requireEETriggers)) {
        allPlots["lp_pt_ee_check"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_pt_ee_check_no_weight"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_eta_ee_check"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_eta_ee_check_no_weight"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_phi_ee_check"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_phi_ee_check_no_weight"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_iso_ee_check"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_iso_ee_check_no_weight"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_id_ee_check"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
        allPlots["lp_id_ee_check_no_weight"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
      }
      if (hasEMTrigger && (!ev.isData || requireEMTriggers)) {
        allPlots["lp_pt_em_check"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_pt_em_check_no_weight"]->Fill(ev.l_pt[il],norm);
        allPlots["lp_eta_em_check"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_eta_em_check_no_weight"]->Fill(ev.l_eta[il],norm);
        allPlots["lp_phi_em_check"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_phi_em_check_no_weight"]->Fill(ev.l_phi[il],norm);
        allPlots["lp_iso_em_check"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_iso_em_check_no_weight"]->Fill(ev.l_relIso[il],norm);
        allPlots["lp_id_em_check"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
        allPlots["lp_id_em_check_no_weight"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
      }
      allPlots["lp_pt_all_check"]->Fill(ev.l_pt[il],norm);
      allPlots["lp_pt_all_check_no_weight"]->Fill(ev.l_pt[il],norm);
      allPlots["lp_eta_all_check"]->Fill(ev.l_eta[il],norm);
      allPlots["lp_eta_all_check_no_weight"]->Fill(ev.l_eta[il],norm);
      allPlots["lp_phi_all_check"]->Fill(ev.l_phi[il],norm);
      allPlots["lp_phi_all_check_no_weight"]->Fill(ev.l_phi[il],norm);
      allPlots["lp_iso_all_check"]->Fill(ev.l_relIso[il],norm);
      allPlots["lp_iso_all_check_no_weight"]->Fill(ev.l_relIso[il],norm);
      allPlots["lp_id_all_check"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);
      allPlots["lp_id_all_check_no_weight"]->Fill(ev.l_id[il]*ev.l_charge[il],norm);

      // FIXME : to be updated with the latest recipes https://indico.cern.ch/event/600194/ 
      /* About the content of l_pid
       - for e : l_pid=(passVetoId | (passTightId<<1) | (passTightIdExceptIso<<2));
       - for mu : l_pid=(isMedium | (isTight<<1));
      */
      if (abs(ev.l_id[il]) == 13) {
        if (ev.l_pt[il] < 15.) continue;
        if (fabs(ev.l_eta[il]) > 2.4) continue;
        if (ev.l_relIso[il] < 0.25) 
          vetoLeptons.push_back(il);
        if (!(ev.l_pid[il] &0x1)) continue;
        if (ev.l_pt[il] < 20.) continue;
        if (fabs(ev.l_eta[il]) > 2.4) continue;
        if (ev.l_relIso[il] < 0.15)
          medLeptons.push_back(il);
        if (!((ev.l_pid[il]>>1)&0x1)) continue;
        if (ev.l_pt[il] < 26.) continue;
        if (fabs(ev.l_eta[il]) > 2.1) continue;
        if (ev.l_relIso[il] < 0.15) 
          tightLeptons.push_back(il);
      } else if (abs(ev.l_id[il]) == 11) {
        if (!(ev.l_pid[il] &0x1)) continue;
        if (ev.l_pt[il] < 15.) continue;
        if (fabs(ev.l_eta[il]) > 2.4) continue;
        if (fabs(ev.l_eta[il]) < 2.4 && ev.l_relIso[il] < 0.15)
          vetoLeptons.push_back(il);
        if (ev.l_pt[il] < 20.) continue;
        if (fabs(ev.l_eta[il]) > 2.4) continue;
        if ((fabs(ev.l_eta[il]) < 1.4442 && ev.l_relIso[il] < 0.0893)
            || (fabs(ev.l_eta[il]) > 1.5660 && fabs(ev.l_eta[il]) < 2.4 && ev.l_relIso[il] < 0.121))
          medLeptons.push_back(il);
        if (!((ev.l_pid[il]>>2)&0x1)) continue;
        if (ev.l_pt[il] < 30.) continue;
        if (fabs(ev.l_eta[il]) > 2.4) continue;
        if ((fabs(ev.l_eta[il]) < 1.4442 && ev.l_relIso[il] < 0.0893)
            || (fabs(ev.l_eta[il]) > 1.5660 && fabs(ev.l_eta[il]) < 2.4 && ev.l_relIso[il] < 0.121))
          tightLeptons.push_back(il);
      }
    }

    if (hasMuTrigger && (!ev.isData || requireMutriggerOnly)) {
      allPlots["vlp_n_m_check"]->Fill(vetoLeptons.size(),norm);
      allPlots["vlp_n_m_check_no_weight"]->Fill(vetoLeptons.size(),norm);
      allPlots["mlp_n_m_check"]->Fill(medLeptons.size(),norm);
      allPlots["mlp_n_m_check_no_weight"]->Fill(medLeptons.size(),norm);
      allPlots["tlp_n_m_check"]->Fill(tightLeptons.size(),norm);
      allPlots["tlp_n_m_check_no_weight"]->Fill(tightLeptons.size(),norm);
    }
    if (hasEleTrigger && (!ev.isData || requireEletriggerOnly)) {
      allPlots["vlp_n_e_check"]->Fill(vetoLeptons.size(),norm);
      allPlots["vlp_n_e_check_no_weight"]->Fill(vetoLeptons.size(),norm);
      allPlots["mlp_n_e_check"]->Fill(medLeptons.size(),norm);
      allPlots["mlp_n_e_check_no_weight"]->Fill(medLeptons.size(),norm);
      allPlots["tlp_n_e_check"]->Fill(tightLeptons.size(),norm);
      allPlots["tlp_n_e_check_no_weight"]->Fill(tightLeptons.size(),norm);
    }
    if (hasMMTrigger && (!ev.isData || requireMMTriggers)) {
      allPlots["vlp_n_mm_check"]->Fill(vetoLeptons.size(),norm);
      allPlots["vlp_n_mm_check_no_weight"]->Fill(vetoLeptons.size(),norm);
      allPlots["mlp_n_mm_check"]->Fill(medLeptons.size(),norm);
      allPlots["mlp_n_mm_check_no_weight"]->Fill(medLeptons.size(),norm);
      allPlots["tlp_n_mm_check"]->Fill(tightLeptons.size(),norm);
      allPlots["tlp_n_mm_check_no_weight"]->Fill(tightLeptons.size(),norm);
    }
    if (hasEETrigger && (!ev.isData || requireEETriggers)) {
      allPlots["vlp_n_ee_check"]->Fill(vetoLeptons.size(),norm);
      allPlots["vlp_n_ee_check_no_weight"]->Fill(vetoLeptons.size(),norm);
      allPlots["mlp_n_ee_check"]->Fill(medLeptons.size(),norm);
      allPlots["mlp_n_ee_check_no_weight"]->Fill(medLeptons.size(),norm);
      allPlots["tlp_n_ee_check"]->Fill(tightLeptons.size(),norm);
      allPlots["tlp_n_ee_check_no_weight"]->Fill(tightLeptons.size(),norm);
    }
    if (hasEMTrigger && (!ev.isData || requireEMTriggers)) {
      allPlots["vlp_n_em_check"]->Fill(vetoLeptons.size(),norm);
      allPlots["vlp_n_em_check_no_weight"]->Fill(vetoLeptons.size(),norm);
      allPlots["mlp_n_em_check"]->Fill(medLeptons.size(),norm);
      allPlots["mlp_n_em_check_no_weight"]->Fill(medLeptons.size(),norm);
      allPlots["tlp_n_em_check"]->Fill(tightLeptons.size(),norm);
      allPlots["tlp_n_em_check_no_weight"]->Fill(tightLeptons.size(),norm);
    }
    allPlots["vlp_n_all_check"]->Fill(vetoLeptons.size(),norm);
    allPlots["vlp_n_all_check_no_weight"]->Fill(vetoLeptons.size(),norm);
    allPlots["mlp_n_all_check"]->Fill(medLeptons.size(),norm);
    allPlots["mlp_n_all_check_no_weight"]->Fill(medLeptons.size(),norm);
    allPlots["tlp_n_all_check"]->Fill(tightLeptons.size(),norm);
    allPlots["tlp_n_all_check_no_weight"]->Fill(tightLeptons.size(),norm);
    for (size_t il = 0; il < vetoLeptons.size(); il++) {
      if (hasMuTrigger && (!ev.isData || requireMutriggerOnly)) {
        allPlots["vlp_pt_m_check"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_pt_m_check_no_weight"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_eta_m_check"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_eta_m_check_no_weight"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_phi_m_check"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_phi_m_check_no_weight"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_iso_m_check"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_iso_m_check_no_weight"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_id_m_check"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
        allPlots["vlp_id_m_check_no_weight"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
      }
      if (hasEleTrigger && (!ev.isData || requireEletriggerOnly)) {
        allPlots["vlp_pt_e_check"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_pt_e_check_no_weight"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_eta_e_check"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_eta_e_check_no_weight"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_phi_e_check"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_phi_e_check_no_weight"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_iso_e_check"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_iso_e_check_no_weight"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_id_e_check"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
        allPlots["vlp_id_e_check_no_weight"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
      }
      if (hasMMTrigger && (!ev.isData || requireMMTriggers)) {
        allPlots["vlp_pt_mm_check"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_pt_mm_check_no_weight"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_eta_mm_check"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_eta_mm_check_no_weight"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_phi_mm_check"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_phi_mm_check_no_weight"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_iso_mm_check"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_iso_mm_check_no_weight"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_id_mm_check"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
        allPlots["vlp_id_mm_check_no_weight"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
      }
      if (hasEETrigger && (!ev.isData || requireEETriggers)) {
        allPlots["vlp_pt_ee_check"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_pt_ee_check_no_weight"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_eta_ee_check"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_eta_ee_check_no_weight"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_phi_ee_check"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_phi_ee_check_no_weight"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_iso_ee_check"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_iso_ee_check_no_weight"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_id_ee_check"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
        allPlots["vlp_id_ee_check_no_weight"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
      }
      if (hasEMTrigger && (!ev.isData || requireEMTriggers)) {
        allPlots["vlp_pt_em_check"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_pt_em_check_no_weight"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
        allPlots["vlp_eta_em_check"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_eta_em_check_no_weight"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
        allPlots["vlp_phi_em_check"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_phi_em_check_no_weight"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
        allPlots["vlp_iso_em_check"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_iso_em_check_no_weight"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
        allPlots["vlp_id_em_check"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
        allPlots["vlp_id_em_check_no_weight"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
      }
      allPlots["vlp_pt_all_check"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
      allPlots["vlp_pt_all_check_no_weight"]->Fill(ev.l_pt[vetoLeptons[il]],norm);
      allPlots["vlp_eta_all_check"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
      allPlots["vlp_eta_all_check_no_weight"]->Fill(ev.l_eta[vetoLeptons[il]],norm);
      allPlots["vlp_phi_all_check"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
      allPlots["vlp_phi_all_check_no_weight"]->Fill(ev.l_phi[vetoLeptons[il]],norm);
      allPlots["vlp_iso_all_check"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
      allPlots["vlp_iso_all_check_no_weight"]->Fill(ev.l_relIso[vetoLeptons[il]],norm);
      allPlots["vlp_id_all_check"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
      allPlots["vlp_id_all_check_no_weight"]->Fill(ev.l_id[vetoLeptons[il]]*ev.l_charge[vetoLeptons[il]],norm);
    }
    for (size_t il = 0; il < medLeptons.size(); il++) {
      if (hasMuTrigger && (!ev.isData || requireMutriggerOnly)) {
        allPlots["mlp_pt_m_check"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_pt_m_check_no_weight"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_eta_m_check"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_eta_m_check_no_weight"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_phi_m_check"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_phi_m_check_no_weight"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_iso_m_check"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_iso_m_check_no_weight"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_id_m_check"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
        allPlots["mlp_id_m_check_no_weight"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
      }
      if (hasEleTrigger && (!ev.isData || requireEletriggerOnly)) {
        allPlots["mlp_pt_e_check"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_pt_e_check_no_weight"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_eta_e_check"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_eta_e_check_no_weight"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_phi_e_check"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_phi_e_check_no_weight"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_iso_e_check"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_iso_e_check_no_weight"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_id_e_check"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
        allPlots["mlp_id_e_check_no_weight"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
      }
      if (hasMMTrigger && (!ev.isData || requireMMTriggers)) {
        allPlots["mlp_pt_mm_check"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_pt_mm_check_no_weight"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_eta_mm_check"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_eta_mm_check_no_weight"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_phi_mm_check"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_phi_mm_check_no_weight"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_iso_mm_check"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_iso_mm_check_no_weight"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_id_mm_check"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
        allPlots["mlp_id_mm_check_no_weight"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
      }
      if (hasEETrigger && (!ev.isData || requireEETriggers)) {
        allPlots["mlp_pt_ee_check"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_pt_ee_check_no_weight"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_eta_ee_check"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_eta_ee_check_no_weight"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_phi_ee_check"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_phi_ee_check_no_weight"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_iso_ee_check"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_iso_ee_check_no_weight"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_id_ee_check"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
        allPlots["mlp_id_ee_check_no_weight"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
      }
      if (hasEMTrigger && (!ev.isData || requireEMTriggers)) {
        allPlots["mlp_pt_em_check"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_pt_em_check_no_weight"]->Fill(ev.l_pt[medLeptons[il]],norm);
        allPlots["mlp_eta_em_check"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_eta_em_check_no_weight"]->Fill(ev.l_eta[medLeptons[il]],norm);
        allPlots["mlp_phi_em_check"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_phi_em_check_no_weight"]->Fill(ev.l_phi[medLeptons[il]],norm);
        allPlots["mlp_iso_em_check"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_iso_em_check_no_weight"]->Fill(ev.l_relIso[medLeptons[il]],norm);
        allPlots["mlp_id_em_check"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
        allPlots["mlp_id_em_check_no_weight"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
      }
      allPlots["mlp_pt_all_check"]->Fill(ev.l_pt[medLeptons[il]],norm);
      allPlots["mlp_pt_all_check_no_weight"]->Fill(ev.l_pt[medLeptons[il]],norm);
      allPlots["mlp_eta_all_check"]->Fill(ev.l_eta[medLeptons[il]],norm);
      allPlots["mlp_eta_all_check_no_weight"]->Fill(ev.l_eta[medLeptons[il]],norm);
      allPlots["mlp_phi_all_check"]->Fill(ev.l_phi[medLeptons[il]],norm);
      allPlots["mlp_phi_all_check_no_weight"]->Fill(ev.l_phi[medLeptons[il]],norm);
      allPlots["mlp_iso_all_check"]->Fill(ev.l_relIso[medLeptons[il]],norm);
      allPlots["mlp_iso_all_check_no_weight"]->Fill(ev.l_relIso[medLeptons[il]],norm);
      allPlots["mlp_id_all_check"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
      allPlots["mlp_id_all_check_no_weight"]->Fill(ev.l_id[medLeptons[il]]*ev.l_charge[medLeptons[il]],norm);
    }
    for (size_t il = 0; il < tightLeptons.size(); il++) {
      if (hasMuTrigger && (!ev.isData || requireMutriggerOnly)) {
        allPlots["tlp_pt_m_check"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_pt_m_check_no_weight"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_eta_m_check"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_eta_m_check_no_weight"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_phi_m_check"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_phi_m_check_no_weight"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_iso_m_check"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_iso_m_check_no_weight"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_id_m_check"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
        allPlots["tlp_id_m_check_no_weight"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
      }
      if (hasEleTrigger && (!ev.isData || requireEletriggerOnly)) {
        allPlots["tlp_pt_e_check"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_pt_e_check_no_weight"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_eta_e_check"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_eta_e_check_no_weight"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_phi_e_check"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_phi_e_check_no_weight"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_iso_e_check"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_iso_e_check_no_weight"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_id_e_check"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
        allPlots["tlp_id_e_check_no_weight"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
      }
      if (hasMMTrigger && (!ev.isData || requireMMTriggers)) {
        allPlots["tlp_pt_mm_check"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_pt_mm_check_no_weight"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_eta_mm_check"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_eta_mm_check_no_weight"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_phi_mm_check"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_phi_mm_check_no_weight"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_iso_mm_check"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_iso_mm_check_no_weight"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_id_mm_check"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
        allPlots["tlp_id_mm_check_no_weight"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
      }
      if (hasEETrigger && (!ev.isData || requireEETriggers)) {
        allPlots["tlp_pt_ee_check"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_pt_ee_check_no_weight"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_eta_ee_check"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_eta_ee_check_no_weight"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_phi_ee_check"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_phi_ee_check_no_weight"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_iso_ee_check"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_iso_ee_check_no_weight"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_id_ee_check"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
        allPlots["tlp_id_ee_check_no_weight"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
      }
      if (hasEMTrigger && (!ev.isData || requireEMTriggers)) {
        allPlots["tlp_pt_em_check"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_pt_em_check_no_weight"]->Fill(ev.l_pt[tightLeptons[il]],norm);
        allPlots["tlp_eta_em_check"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_eta_em_check_no_weight"]->Fill(ev.l_eta[tightLeptons[il]],norm);
        allPlots["tlp_phi_em_check"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_phi_em_check_no_weight"]->Fill(ev.l_phi[tightLeptons[il]],norm);
        allPlots["tlp_iso_em_check"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_iso_em_check_no_weight"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
        allPlots["tlp_id_em_check"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
        allPlots["tlp_id_em_check_no_weight"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
      }
      allPlots["tlp_pt_all_check"]->Fill(ev.l_pt[tightLeptons[il]],norm);
      allPlots["tlp_pt_all_check_no_weight"]->Fill(ev.l_pt[tightLeptons[il]],norm);
      allPlots["tlp_eta_all_check"]->Fill(ev.l_eta[tightLeptons[il]],norm);
      allPlots["tlp_eta_all_check_no_weight"]->Fill(ev.l_eta[tightLeptons[il]],norm);
      allPlots["tlp_phi_all_check"]->Fill(ev.l_phi[tightLeptons[il]],norm);
      allPlots["tlp_phi_all_check_no_weight"]->Fill(ev.l_phi[tightLeptons[il]],norm);
      allPlots["tlp_iso_all_check"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
      allPlots["tlp_iso_all_check_no_weight"]->Fill(ev.l_relIso[tightLeptons[il]],norm);
      allPlots["tlp_id_all_check"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
      allPlots["tlp_id_all_check_no_weight"]->Fill(ev.l_id[tightLeptons[il]]*ev.l_charge[tightLeptons[il]],norm);
    }

    if (debug) cout << "lepton selection DONE -- # tight = " << tightLeptons.size() << ", # med = " << medLeptons.size() << ", #veto = " << vetoLeptons.size() << endl;

    //decide the channel
    if (debug) cout << "Decide channel" << endl;
    TString chTag("");
    std::vector<int> selLeptons;
    if (hasMuTrigger && (!ev.isData || requireMutriggerOnly)) {
      if (tightLeptons.size() == 1) {
        allPlots["yields_m_check"]->Fill((double)iCut["_m"],norm);
        allPlots["yields_m_check_no_weight"]->Fill((double)iCut["_m"],norm);
        ++iCut["_m"];
        if (abs(ev.l_id[tightLeptons[0]]) == 13) {
          allPlots["yields_m_check"]->Fill((double)iCut["_m"],norm);
          allPlots["yields_m_check_no_weight"]->Fill((double)iCut["_m"],norm);
          ++iCut["_m"];
          if (vetoLeptons.size() == 1 && vetoLeptons[0] == tightLeptons[0]) {
            allPlots["yields_m_check"]->Fill((double)iCut["_m"],norm);
            allPlots["yields_m_check_no_weight"]->Fill((double)iCut["_m"],norm);
            ++iCut["_m"];
            chTag = "m";
            selLeptons.push_back(tightLeptons[0]);
          }
        }
      }
    }
    if (hasEleTrigger && (!ev.isData || requireEletriggerOnly)) {
      if (tightLeptons.size() == 1) {
        allPlots["yields_e_check"]->Fill((double)iCut["_e"],norm);
        allPlots["yields_e_check_no_weight"]->Fill((double)iCut["_e"],norm);
        ++iCut["_e"];
        if (abs(ev.l_id[tightLeptons[0]]) == 11) {
          allPlots["yields_e_check"]->Fill((double)iCut["_e"],norm);
          allPlots["yields_e_check_no_weight"]->Fill((double)iCut["_e"],norm);
          ++iCut["_e"];
          if (vetoLeptons.size() == 1 && vetoLeptons[0] == tightLeptons[0]) {
            allPlots["yields_e_check"]->Fill((double)iCut["_e"],norm);
            allPlots["yields_e_check_no_weight"]->Fill((double)iCut["_e"],norm);
            ++iCut["_e"];
            chTag = "e";
            selLeptons.push_back(tightLeptons[0]);
          }
        }
      }
    }
    if (hasMMTrigger && (!ev.isData || requireMMTriggers)) {
      if (medLeptons.size() >= 2) {
        allPlots["yields_mm_check"]->Fill((double)iCut["_mm"],norm);
        allPlots["yields_mm_check_no_weight"]->Fill((double)iCut["_mm"],norm);
        ++iCut["_mm"];
        if (abs(ev.l_id[medLeptons[0]]) == 13 && abs(ev.l_id[medLeptons[1]]) == 13) {
          allPlots["yields_mm_check"]->Fill((double)iCut["_mm"],norm);
          allPlots["yields_mm_check_no_weight"]->Fill((double)iCut["_mm"],norm);
          ++iCut["_mm"];
          if (ev.l_charge[medLeptons[0]]*ev.l_charge[medLeptons[1]] < 0.) {
            allPlots["yields_mm_check"]->Fill((double)iCut["_mm"],norm);
            allPlots["yields_mm_check_no_weight"]->Fill((double)iCut["_mm"],norm);
            ++iCut["_mm"];
            chTag = "mm";
            selLeptons.push_back(medLeptons[0]);
            selLeptons.push_back(medLeptons[1]);
          }
        }
      }
    }
    if (hasEETrigger && (!ev.isData || requireEETriggers)) {
      if (medLeptons.size() >= 2) {
        allPlots["yields_ee_check"]->Fill((double)iCut["_ee"],norm);
        allPlots["yields_ee_check_no_weight"]->Fill((double)iCut["_ee"],norm);
        ++iCut["_ee"];
        if (abs(ev.l_id[medLeptons[0]]) == 11 && abs(ev.l_id[medLeptons[1]]) == 11 && ev.l_pt[medLeptons[0]] > 25.) { // to be consistent with the trigger
          allPlots["yields_ee_check"]->Fill((double)iCut["_ee"],norm);
          allPlots["yields_ee_check_no_weight"]->Fill((double)iCut["_ee"],norm);
          ++iCut["_ee"];      
          if (ev.l_charge[medLeptons[0]]*ev.l_charge[medLeptons[1]] < 0.) {
            allPlots["yields_ee_check"]->Fill((double)iCut["_ee"],norm);
            allPlots["yields_ee_check_no_weight"]->Fill((double)iCut["_ee"],norm);
            ++iCut["_ee"];      
            chTag = "ee";
            selLeptons.push_back(medLeptons[0]);
            selLeptons.push_back(medLeptons[1]);
          }
        }
      }
    }
    if (hasEMTrigger && (!ev.isData || requireEMTriggers)) {
      if (medLeptons.size() >= 2) {
        allPlots["yields_em_check"]->Fill((double)iCut["_em"],norm);
        allPlots["yields_em_check_no_weight"]->Fill((double)iCut["_em"],norm);
        ++iCut["_em"];
        if ((abs(ev.l_id[medLeptons[0]]) == 11 && abs(ev.l_id[medLeptons[1]]) == 13 && ev.l_pt[medLeptons[0]] > 25.)
            || (abs(ev.l_id[medLeptons[0]]) == 13 && abs(ev.l_id[medLeptons[1]]) == 11 && ev.l_pt[medLeptons[1]] > 25.)) { // to be consistent with the trigger
          allPlots["yields_em_check"]->Fill((double)iCut["_em"],norm);
          allPlots["yields_em_check_no_weight"]->Fill((double)iCut["_em"],norm);
          ++iCut["_em"];      
          if (ev.l_charge[medLeptons[0]]*ev.l_charge[medLeptons[1]] < 0.) {
            allPlots["yields_em_check"]->Fill((double)iCut["_em"],norm);
            allPlots["yields_em_check_no_weight"]->Fill((double)iCut["_em"],norm);
            ++iCut["_em"];      
            chTag = "em";
            selLeptons.push_back(medLeptons[0]);
            selLeptons.push_back(medLeptons[1]);
          }
        }
      }
    }

    if (chTag == "") continue;
    chTag = "_"+chTag;
    if (debug) cout << "decide channel DONE -- channel " << chTag << ", # sel = " << selLeptons.size() << endl;

    //save leading lepton kinematics
    std::vector<TLorentzVector> leptons;
    for (size_t il = 0; il < selLeptons.size(); il++) {
      int lepIdx=selLeptons[il];
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);
      leptons.push_back(lp4);
    }
    std::sort(leptons.begin(), leptons.end(), VecSort);

    //event weight
    float wgt(1.0);
    std::vector<float> puWgts(3,1.0),topPtWgts(2,1.0);
    EffCorrection_t lepSelCorrWgt(1.0,0.0), triggerCorrWgt(1.0,0.0);
    if (debug) cout << "Lepton scale factors" << endl;
    if (!ev.isData) {
      //update lepton selection scale factors, if found
      //float lepTriggerSF(1.0),lepSelSF(1.0);
      //FIXME

      //account for pu weights and effect on normalization
      allPlots["puwgtctr"]->Fill(0.,1.0);
      if (debug) cout << "getting puWgts" << endl;
      for (size_t iwgt=0; iwgt<3; iwgt++) {
        puWgts[iwgt]=puWgtGr[iwgt]->Eval(ev.putrue);  
        allPlots["puwgtctr"]->Fill(iwgt+1,puWgts[iwgt]);
      }
      if (debug) cout << "getting puWgts DONE" << endl;
      //trigger/id+iso efficiency corrections
      if (debug) cout << "calling trigger function" << endl;
      std::vector<int> pdgIds;
      for (size_t ilp = 0; ilp < selLeptons.size(); ilp++)
        pdgIds.push_back(ev.l_id[selLeptons[ilp]]);
      triggerCorrWgt=lepEffH.getTriggerCorrection(pdgIds,leptons);
      if (debug) cout << "calling trigger function DONE" << endl;
      // ** selLeptons contains only ev_l position, leptons contains p4 **
      for (size_t il=0; il<selLeptons.size(); il++) {
        EffCorrection_t selSF=lepEffH.getOfflineCorrection(pdgIds[il],leptons[il].Pt(),leptons[il].Eta());
        lepSelCorrWgt.second = sqrt( pow(lepSelCorrWgt.first*selSF.second,2)+pow(lepSelCorrWgt.second*selSF.first,2));
        if (debug) cout << "lepSelCorrWgt=" << lepSelCorrWgt.first << endl;
        if (debug) cout << "selSF=" << selSF.first << endl;
        lepSelCorrWgt.first *= selSF.first;
      }

      //update nominal event weight
      //wgt=lepTriggerSF*lepSelSF*puWeight*norm;
      wgt=triggerCorrWgt.first*lepSelCorrWgt.first*puWgts[0]*norm;
      if (debug) cout << "weight=" << wgt << endl;
      if (debug) cout << "Trigger=" << triggerCorrWgt.first << endl << "Lepton=" << lepSelCorrWgt.first << endl << "PU=" << puWgts[0] << endl << "norm=" << norm  << endl;
      for (size_t il = 0; il < leptons.size(); il++) {
        if (!filename.Contains("_WJets")) continue;
        cout << "pT: " << leptons[il].Pt() << endl;
      }
      //wgt=1.0;
    }
    if (debug) cout << "lepton scale factors DONE -- norm = " << norm << ", wgt = " << wgt << endl;

    allPlots["yields_all_check"]->Fill((double)iCut["_all"],wgt);
    allPlots["yields_all_check_no_weight"]->Fill((double)iCut["_all"],norm);
    ++iCut["_all"];

    allPlots["slp_n"+chTag+"_check"]->Fill(leptons.size(),wgt);
    allPlots["slp_n"+chTag+"_check_no_weight"]->Fill(leptons.size(),norm);
    allPlots["slp_n_all_check"]->Fill(leptons.size(),wgt);
    allPlots["slp_n_all_check_no_weight"]->Fill(leptons.size(),norm);
    for (size_t il = 0; il < leptons.size(); il++) {
      allPlots["slp_pt"+chTag+"_check"]->Fill(leptons[il].Pt(),wgt);
      allPlots["slp_pt"+chTag+"_check_no_weight"]->Fill(leptons[il].Pt(),norm);
      allPlots["slp_eta"+chTag+"_check"]->Fill(leptons[il].Eta(),wgt);
      allPlots["slp_eta"+chTag+"_check_no_weight"]->Fill(leptons[il].Eta(),norm);
      allPlots["slp_phi"+chTag+"_check"]->Fill(leptons[il].Phi(),wgt);
      allPlots["slp_phi"+chTag+"_check_no_weight"]->Fill(leptons[il].Phi(),norm);
      allPlots["slp_pt_all_check"]->Fill(leptons[il].Pt(),wgt);
      allPlots["slp_pt_all_check_no_weight"]->Fill(leptons[il].Pt(),norm);
      allPlots["slp_eta_all_check"]->Fill(leptons[il].Eta(),wgt);
      allPlots["slp_eta_all_check_no_weight"]->Fill(leptons[il].Eta(),norm);
      allPlots["slp_phi_all_check"]->Fill(leptons[il].Phi(),wgt);
      allPlots["slp_phi_all_check_no_weight"]->Fill(leptons[il].Phi(),norm);
    }

    //select jets
    if (debug) cout << "Starting jet selection -- # all = " << ev.nj << endl;

    TLorentzVector jetDiff(0,0,0,0);
    int nbjets(0),ncjets(0),nljets(0);//,leadingJetIdx(-wgt);
    std::vector<int> resolvedJetIdx;
    std::vector<TLorentzVector> resolvedJetP4;
    std::vector<Jet> bJetsVec, lightJetsVec, allJetsVec;
    int nJet40 = 0;
    allPlots["j_n"+chTag+"_check"]->Fill(ev.nj,wgt);
    allPlots["j_n"+chTag+"_check_no_weight"]->Fill(ev.nj,norm);
    allPlots["j_n_all_check"]->Fill(ev.nj,wgt);
    allPlots["j_n_all_check_no_weight"]->Fill(ev.nj,norm);
    for (int k = 0; k < ev.nj; k++) {
      allPlots["j_pt"+chTag+"_check"]->Fill(ev.j_pt[k],wgt);
      allPlots["j_pt"+chTag+"_check_no_weight"]->Fill(ev.j_pt[k],norm);
      allPlots["j_eta"+chTag+"_check"]->Fill(ev.j_eta[k],wgt);
      allPlots["j_eta"+chTag+"_check_no_weight"]->Fill(ev.j_eta[k],norm);
      allPlots["j_phi"+chTag+"_check"]->Fill(ev.j_phi[k],wgt);
      allPlots["j_phi"+chTag+"_check_no_weight"]->Fill(ev.j_phi[k],norm);
      allPlots["j_csv"+chTag+"_check"]->Fill(ev.j_csv[k],wgt);
      allPlots["j_csv"+chTag+"_check_no_weight"]->Fill(ev.j_csv[k],norm);
      allPlots["j_pt_all_check"]->Fill(ev.j_pt[k],wgt);
      allPlots["j_pt_all_check_no_weight"]->Fill(ev.j_pt[k],norm);
      allPlots["j_eta_all_check"]->Fill(ev.j_eta[k],wgt);
      allPlots["j_eta_all_check_no_weight"]->Fill(ev.j_eta[k],norm);
      allPlots["j_phi_all_check"]->Fill(ev.j_phi[k],wgt);
      allPlots["j_phi_all_check_no_weight"]->Fill(ev.j_phi[k],norm);
      allPlots["j_csv_all_check"]->Fill(ev.j_csv[k],wgt);
      allPlots["j_csv_all_check_no_weight"]->Fill(ev.j_csv[k],norm);

      //check kinematics
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

      //cross clean with respect to leptons 
      bool overlapsWithLepton(false);
      for (size_t il = 0; il < leptons.size(); il++) {
        if (jp4.DeltaR(leptons[il]) > 0.4) continue;
        overlapsWithLepton = true;
      }
      if (overlapsWithLepton) continue;
      if (debug) cout << "Overlap with lepton PASSED" << endl;

      //smear jet energy resolution for MC
      //jetDiff -= jp4;
      float genJet_pt(0);
      if (ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
      if (!ev.isData && genJet_pt>0) 
      {
        float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
        jp4 *= jerSmear;
      }
      //jetDiff += jp4;
      resolvedJetIdx.push_back(k);
      resolvedJetP4.push_back(jp4);

      // re-inforce kinematics cuts
      if (jp4.Pt() < 30.) continue;
      if (fabs(jp4.Eta()) > 2.4) continue;
      if (jp4.Pt() > 40.) ++nJet40;

      //b-tag
      if (debug) cout << "Starting b-tagging" << endl;
      float csv = ev.j_csv[k];	  
      bool isBTagged(csv > 0.800);
      if (!ev.isData) {
        float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
        float expEff(1.0), jetBtagSF(1.0);
        if (abs(ev.j_hadflav[k]) == 4) { 
          ncjets++;
          expEff    = expBtagEff["c"]->Eval(jptForBtag); 
          jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
          jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
        } else if (abs(ev.j_hadflav[k]) == 5) { 
          nbjets++;
          expEff    = expBtagEff["b"]->Eval(jptForBtag); 
          jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
          jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
        } else {
          nljets++;
          expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
          jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
          jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
        }

        //updated b-tagging decision with the data/MC scale factor
        myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
      }
      if (debug) cout << "b-tagging DONE -- isBTagged = " << isBTagged << endl;

      //save jet
      Jet tmpj(jp4, csv, k);
      for (int ipf = 0; ipf < ev.npf; ipf++) {
        if (ev.pf_j[ipf] != k) continue;
        if (ev.pf_c[ipf]==0) continue;
        TLorentzVector tkP4(0,0,0,0);
        tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
        pfTrack pftk(tkP4, ev.pf_dxy[ipf], ev.pf_dxyE[ipf], ev.pf_dz[ipf], ev.pf_dzE[ipf], ev.pf_id[ipf]);
        tmpj.addTrack(pftk,ev.pf_id[ipf]);
      }
      tmpj.sortTracksByPt();

      if (isBTagged) bJetsVec.push_back(tmpj);
      else lightJetsVec.push_back(tmpj);
      allJetsVec.push_back(tmpj);
    }
    //sort by Pt
    sort(lightJetsVec.begin(),    lightJetsVec.end(),   sortJetsByPt);
    sort(bJetsVec.begin(),    bJetsVec.end(),   sortJetsByPt);
    sort(allJetsVec.begin(),  allJetsVec.end(), sortJetsByPt);

    if (debug) cout << "jet selection DONE -- # all = " << allJetsVec.size() << ", # light = " << lightJetsVec.size() << ", # b = " << bJetsVec.size() << endl;

    if (selLeptons.size() == 1) {
      if (allJetsVec.size() < 4) continue;
      allPlots["yields"+chTag+"_check"]->Fill((double)iCut[chTag],wgt);
      allPlots["yields"+chTag+"_check_no_weight"]->Fill((double)iCut[chTag],norm);
      ++iCut[chTag];
    }
    if (selLeptons.size() == 2) {
      if (allJetsVec.size() < 2) continue;
      allPlots["yields"+chTag+"_check"]->Fill((double)iCut[chTag],wgt);
      allPlots["yields"+chTag+"_check_no_weight"]->Fill((double)iCut[chTag],norm);
      ++iCut[chTag];
    }
    allPlots["yields_all_check"]->Fill((double)iCut["_all"],wgt);
    allPlots["yields_all_check_no_weight"]->Fill((double)iCut["_all"],norm);
    ++iCut["_all"];

    allPlots["sj_n"+chTag+"_check"]->Fill(allJetsVec.size(),wgt);
    allPlots["sj_n"+chTag+"_check_no_weight"]->Fill(allJetsVec.size(),norm);
    allPlots["sj_n_all_check"]->Fill(allJetsVec.size(),wgt);
    allPlots["sj_n_all_check_no_weight"]->Fill(allJetsVec.size(),norm);
    allPlots["sj40_n"+chTag+"_check"]->Fill(nJet40,wgt);
    allPlots["sj40_n"+chTag+"_check_no_weight"]->Fill(nJet40,norm);
    allPlots["sj40_n_all_check"]->Fill(nJet40,wgt);
    allPlots["sj40_n_all_check_no_weight"]->Fill(nJet40,norm);
    for (size_t ij = 0; ij < allJetsVec.size(); ij++) {
      allPlots["sj_pt"+chTag+"_check"]->Fill(allJetsVec[ij].getVec().Pt(),wgt);
      allPlots["sj_pt"+chTag+"_check_no_weight"]->Fill(allJetsVec[ij].getVec().Pt(),norm);
      allPlots["sj_eta"+chTag+"_check"]->Fill(allJetsVec[ij].getVec().Eta(),wgt);
      allPlots["sj_eta"+chTag+"_check_no_weight"]->Fill(allJetsVec[ij].getVec().Eta(),norm);
      allPlots["sj_phi"+chTag+"_check"]->Fill(allJetsVec[ij].getVec().Phi(),wgt);
      allPlots["sj_phi"+chTag+"_check_no_weight"]->Fill(allJetsVec[ij].getVec().Phi(),norm);
      allPlots["sj_csv"+chTag+"_check"]->Fill(allJetsVec[ij].getCSV(),wgt);
      allPlots["sj_csv"+chTag+"_check_no_weight"]->Fill(allJetsVec[ij].getCSV(),norm);
      allPlots["sj_pt_all_check"]->Fill(allJetsVec[ij].getVec().Pt(),wgt);
      allPlots["sj_pt_all_check_no_weight"]->Fill(allJetsVec[ij].getVec().Pt(),norm);
      allPlots["sj_eta_all_check"]->Fill(allJetsVec[ij].getVec().Eta(),wgt);
      allPlots["sj_eta_all_check_no_weight"]->Fill(allJetsVec[ij].getVec().Eta(),norm);
      allPlots["sj_phi_all_check"]->Fill(allJetsVec[ij].getVec().Phi(),wgt);
      allPlots["sj_phi_all_check_no_weight"]->Fill(allJetsVec[ij].getVec().Phi(),norm);
      allPlots["sj_csv_all_check"]->Fill(allJetsVec[ij].getCSV(),wgt);
      allPlots["sj_csv_all_check_no_weight"]->Fill(allJetsVec[ij].getCSV(),norm);
    }

    allPlots["bj_n"+chTag+"_check"]->Fill(bJetsVec.size(),wgt);
    allPlots["bj_n"+chTag+"_check_no_weight"]->Fill(bJetsVec.size(),norm);
    allPlots["bj_n_all_check"]->Fill(bJetsVec.size(),wgt);
    allPlots["bj_n_all_check_no_weight"]->Fill(bJetsVec.size(),norm);
    for (size_t ij = 0; ij < bJetsVec.size(); ij++) {
      allPlots["bj_pt"+chTag+"_check"]->Fill(bJetsVec[ij].getVec().Pt(),wgt);
      allPlots["bj_pt"+chTag+"_check_no_weight"]->Fill(bJetsVec[ij].getVec().Pt(),norm);
      allPlots["bj_eta"+chTag+"_check"]->Fill(bJetsVec[ij].getVec().Eta(),wgt);
      allPlots["bj_eta"+chTag+"_check_no_weight"]->Fill(bJetsVec[ij].getVec().Eta(),norm);
      allPlots["bj_phi"+chTag+"_check"]->Fill(bJetsVec[ij].getVec().Phi(),wgt);
      allPlots["bj_phi"+chTag+"_check_no_weight"]->Fill(bJetsVec[ij].getVec().Phi(),norm);
      allPlots["bj_csv"+chTag+"_check"]->Fill(bJetsVec[ij].getCSV(),wgt);
      allPlots["bj_csv"+chTag+"_check_no_weight"]->Fill(bJetsVec[ij].getCSV(),norm);
      allPlots["bj_pt_all_check"]->Fill(bJetsVec[ij].getVec().Pt(),wgt);
      allPlots["bj_pt_all_check_no_weight"]->Fill(bJetsVec[ij].getVec().Pt(),norm);
      allPlots["bj_eta_all_check"]->Fill(bJetsVec[ij].getVec().Eta(),wgt);
      allPlots["bj_eta_all_check_no_weight"]->Fill(bJetsVec[ij].getVec().Eta(),norm);
      allPlots["bj_phi_all_check"]->Fill(bJetsVec[ij].getVec().Phi(),wgt);
      allPlots["bj_phi_all_check_no_weight"]->Fill(bJetsVec[ij].getVec().Phi(),norm);
      allPlots["bj_csv_all_check"]->Fill(bJetsVec[ij].getCSV(),wgt);
      allPlots["bj_csv_all_check_no_weight"]->Fill(bJetsVec[ij].getCSV(),norm);
    }

    if (bJetsVec.size() == 0) continue;
    allPlots["yields"+chTag+"_check"]->Fill((double)iCut[chTag],wgt);
    allPlots["yields"+chTag+"_check_no_weight"]->Fill((double)iCut[chTag],norm);
    ++iCut[chTag];
    allPlots["yields_all_check"]->Fill((double)iCut["_all"],wgt);
    allPlots["yields_all_check_no_weight"]->Fill((double)iCut["_all"],norm);
    ++iCut["_all"];

    allPlots["sbj_n"+chTag+"_check"]->Fill(bJetsVec.size(),wgt);
    allPlots["sbj_n"+chTag+"_check_no_weight"]->Fill(bJetsVec.size(),norm);
    allPlots["sbj_n_all_check"]->Fill(bJetsVec.size(),wgt);
    allPlots["sbj_n_all_check_no_weight"]->Fill(bJetsVec.size(),norm);
    for (size_t ij = 0; ij < bJetsVec.size(); ij++) {
      allPlots["sbj_pt"+chTag+"_check"]->Fill(bJetsVec[ij].getVec().Pt(),wgt);
      allPlots["sbj_pt"+chTag+"_check_no_weight"]->Fill(bJetsVec[ij].getVec().Pt(),norm);
      allPlots["sbj_eta"+chTag+"_check"]->Fill(bJetsVec[ij].getVec().Eta(),wgt);
      allPlots["sbj_eta"+chTag+"_check_no_weight"]->Fill(bJetsVec[ij].getVec().Eta(),norm);
      allPlots["sbj_phi"+chTag+"_check"]->Fill(bJetsVec[ij].getVec().Phi(),wgt);
      allPlots["sbj_phi"+chTag+"_check_no_weight"]->Fill(bJetsVec[ij].getVec().Phi(),norm);
      allPlots["sbj_csv"+chTag+"_check"]->Fill(bJetsVec[ij].getCSV(),wgt);
      allPlots["sbj_csv"+chTag+"_check_no_weight"]->Fill(bJetsVec[ij].getCSV(),norm);
      allPlots["sbj_pt_all_check"]->Fill(bJetsVec[ij].getVec().Pt(),wgt);
      allPlots["sbj_pt_all_check_no_weight"]->Fill(bJetsVec[ij].getVec().Pt(),norm);
      allPlots["sbj_eta_all_check"]->Fill(bJetsVec[ij].getVec().Eta(),wgt);
      allPlots["sbj_eta_all_check_no_weight"]->Fill(bJetsVec[ij].getVec().Eta(),norm);
      allPlots["sbj_phi_all_check"]->Fill(bJetsVec[ij].getVec().Phi(),wgt);
      allPlots["sbj_phi_all_check_no_weight"]->Fill(bJetsVec[ij].getVec().Phi(),norm);
      allPlots["sbj_csv_all_check"]->Fill(bJetsVec[ij].getCSV(),wgt);
      allPlots["sbj_csv_all_check_no_weight"]->Fill(bJetsVec[ij].getCSV(),norm);
    }

    for (size_t il = 0; il < leptons.size(); il++) {
      for (size_t ij = 0; ij < bJetsVec.size(); ij++) {
        TLorentzVector jp4 = bJetsVec[ij].getVec();
        allPlots["lp_jets_dR"+chTag+"_check"]->Fill(jp4.DeltaR(leptons[il]),wgt);
        allPlots["lp_jets_dR"+chTag+"_check_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        allPlots["lp_jets_dR_all_check"]->Fill(jp4.DeltaR(leptons[il]),wgt);
        allPlots["lp_jets_dR_all_check_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
      }
      for (size_t ij = 0; ij < lightJetsVec.size(); ij++) {
        TLorentzVector jp4 = lightJetsVec[ij].getVec();
        allPlots["lp_jets_dR"+chTag+"_check"]->Fill(jp4.DeltaR(leptons[il]),wgt);
        allPlots["lp_jets_dR"+chTag+"_check_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        allPlots["lp_jets_dR_all_check"]->Fill(jp4.DeltaR(leptons[il]),wgt);
        allPlots["lp_jets_dR_all_check_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
      }
    }

    //MET and transverse mass
    TLorentzVector met(0,0,0,0);
    met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.); // index 1 for PUPPI MET, but considered as experimental, see slide 8 of https://indico.cern.ch/event/600194/contributions/2423592/attachments/1397230/2130535/2017.01.17_TOP_JetMET-object-review_v3.pdf
    met+=jetDiff;
    met.SetPz(0.); met.SetE(met.Pt());
    allPlots["met"+chTag+"_check"]->Fill(met.Pt(),wgt);
    allPlots["met"+chTag+"_check_no_weight"]->Fill(met.Pt(),norm);
    allPlots["met_all_check"]->Fill(met.Pt(),wgt);
    allPlots["met_all_check_no_weight"]->Fill(met.Pt(),norm);

    // additional cuts in ee and mm candidates
    if (selLeptons.size() == 2 && abs(ev.l_id[selLeptons[0]]) == abs(ev.l_id[selLeptons[1]]) && ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]] < 0) {	  
      if (debug) cout << "MET cut in di-lepton candidates -- MET = " << met.Pt() << " GeV" << endl;
      if (met.Pt() < 40.) continue;
      allPlots["yields"+chTag+"_check"]->Fill((double)iCut[chTag],wgt);
      allPlots["yields"+chTag+"_check_no_weight"]->Fill((double)iCut[chTag],norm);
      ++iCut[chTag];
      if (debug) cout << "MET cut PASSED" << endl;
    }
    allPlots["yields_all_check"]->Fill((double)iCut["_all"],wgt);
    allPlots["yields_all_check_no_weight"]->Fill((double)iCut["_all"],norm);
    ++iCut["_all"];
    if (selLeptons.size() == 2 && abs(ev.l_id[selLeptons[0]]) == abs(ev.l_id[selLeptons[1]]) && ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]] < 0) {	  
      TLorentzVector l1p4, l2p4, dilp4;
      l1p4.SetPtEtaPhiM(ev.l_pt[selLeptons[0]],ev.l_eta[selLeptons[0]],ev.l_phi[selLeptons[0]],ev.l_mass[selLeptons[0]]);
      l2p4.SetPtEtaPhiM(ev.l_pt[selLeptons[1]],ev.l_eta[selLeptons[1]],ev.l_phi[selLeptons[1]],ev.l_mass[selLeptons[1]]);
      dilp4 = l1p4+l2p4;
      if (debug) cout << "DY veto in di-lepton candidates -- M = " << dilp4.M() << " GeV" << endl;
      if (fabs(dilp4.M()-91.) < 15.) { 
        allPlots["z_mass"+chTag+"_check"]->Fill(dilp4.M(),wgt);
        allPlots["z_mass"+chTag+"_check_no_weight"]->Fill(dilp4.M(),norm);
        allPlots["z_mass_all_check"]->Fill(dilp4.M(),wgt);
        allPlots["z_mass_all_check_no_weight"]->Fill(dilp4.M(),norm);
        continue;
      }
      if (debug) cout << "Z veto PASSED" << endl;
      allPlots["yields"+chTag+"_check"]->Fill((double)iCut[chTag],wgt);
      allPlots["yields"+chTag+"_check_no_weight"]->Fill((double)iCut[chTag],norm);
      ++iCut[chTag];
      if (dilp4.M() < 20.) continue;
      if (debug) cout << "Low mass cut PASSED" << endl;
      allPlots["yields"+chTag+"_check"]->Fill((double)iCut[chTag],wgt);
      allPlots["yields"+chTag+"_check_no_weight"]->Fill((double)iCut[chTag],norm);
      ++iCut[chTag];
    }
    allPlots["yields_all_check"]->Fill((double)iCut["_all"],wgt);
    allPlots["yields_all_check_no_weight"]->Fill((double)iCut["_all"],norm);
    ++iCut["_all"];

    if (!ev.isData) {
      /*
      allPlots["norm_topsel"]->Fill(norm);
      allPlots["wgt_topsel"]->Fill(wgt);
      */
      if (fabs(norm) > 1e-15)
        allPlots["sf_topsel"]->Fill(wgt/norm);
    }

    //simple fill
    if (debug) cout << "Starting topsel plots" << endl;

    if (chTag == "_e") {
      allPlots["decay_all_topsel"]->Fill(0.,wgt);
      allPlots["decay_all_topsel_no_weight"]->Fill(0.,norm);
    }
    if (chTag == "_m") {
      allPlots["decay_all_topsel"]->Fill(1.,wgt);
      allPlots["decay_all_topsel_no_weight"]->Fill(1.,norm);
    }
    if (chTag == "_ee") {
      allPlots["decay_all_topsel"]->Fill(2.,wgt);
      allPlots["decay_all_topsel_no_weight"]->Fill(2.,norm);
    }
    if (chTag == "_em") {
      allPlots["decay_all_topsel"]->Fill(3.,wgt);
      allPlots["decay_all_topsel_no_weight"]->Fill(3.,norm);
    }
    if (chTag == "_mm") {
      allPlots["decay_all_topsel"]->Fill(4.,wgt);
      allPlots["decay_all_topsel_no_weight"]->Fill(4.,norm);
    }

    allPlots["lp_n"+chTag+"_topsel"]->Fill(selLeptons.size(),wgt);
    allPlots["lp_n"+chTag+"_topsel_no_weight"]->Fill(selLeptons.size(),norm);
    allPlots["lp_n_all_topsel"]->Fill(selLeptons.size(),wgt);
    allPlots["lp_n_all_topsel_no_weight"]->Fill(selLeptons.size(),norm);
    allPlots["lp1_pt"+chTag+"_topsel"]->Fill(leptons[0].Pt(),wgt);
    allPlots["lp1_pt"+chTag+"_topsel_no_weight"]->Fill(leptons[0].Pt(),norm);
    allPlots["lp1_pt_all_topsel"]->Fill(leptons[0].Pt(),wgt);
    allPlots["lp1_pt_all_topsel_no_weight"]->Fill(leptons[0].Pt(),norm);
    allPlots["lp1_eta"+chTag+"_topsel"]->Fill(leptons[0].Eta(),wgt);
    allPlots["lp1_eta"+chTag+"_topsel_no_weight"]->Fill(leptons[0].Eta(),norm);
    allPlots["lp1_eta_all_topsel"]->Fill(leptons[0].Eta(),wgt);
    allPlots["lp1_eta_all_topsel_no_weight"]->Fill(leptons[0].Eta(),norm);
    allPlots["lp1_phi"+chTag+"_topsel"]->Fill(leptons[0].Phi(),wgt);
    allPlots["lp1_phi"+chTag+"_topsel_no_weight"]->Fill(leptons[0].Phi(),norm);
    allPlots["lp1_phi_all_topsel"]->Fill(leptons[0].Phi(),wgt);
    allPlots["lp1_phi_all_topsel_no_weight"]->Fill(leptons[0].Phi(),norm);
    if (selLeptons.size() == 2) {
      allPlots["lp2_pt"+chTag+"_topsel"]->Fill(leptons[1].Pt(),wgt);
      allPlots["lp2_pt"+chTag+"_topsel_no_weight"]->Fill(leptons[1].Pt(),norm);
      allPlots["lp2_pt_all_topsel"]->Fill(leptons[1].Pt(),wgt);
      allPlots["lp2_pt_all_topsel_no_weight"]->Fill(leptons[1].Pt(),norm);
      allPlots["lp2_eta"+chTag+"_topsel"]->Fill(leptons[1].Eta(),wgt);
      allPlots["lp2_eta"+chTag+"_topsel_no_weight"]->Fill(leptons[1].Eta(),norm);
      allPlots["lp2_eta_all_topsel"]->Fill(leptons[1].Eta(),wgt);
      allPlots["lp2_eta_all_topsel_no_weight"]->Fill(leptons[1].Eta(),norm);
      allPlots["lp2_phi"+chTag+"_topsel"]->Fill(leptons[1].Phi(),wgt);
      allPlots["lp2_phi"+chTag+"_topsel_no_weight"]->Fill(leptons[1].Phi(),norm);
      allPlots["lp2_phi_all_topsel"]->Fill(leptons[1].Phi(),wgt);
      allPlots["lp2_phi_all_topsel_no_weight"]->Fill(leptons[1].Phi(),norm);
    }
    allPlots["met"+chTag+"_topsel"]->Fill(met.Pt(),wgt);
    allPlots["met"+chTag+"_topsel_no_weight"]->Fill(met.Pt(),norm);
    allPlots["met_all_topsel"]->Fill(met.Pt(),wgt);
    allPlots["met_all_topsel_no_weight"]->Fill(met.Pt(),norm);
    if (selLeptons.size() == 2) {
      allPlots["dilp_pt"+chTag+"_topsel"]->Fill((leptons[0]+leptons[1]).Pt(),wgt);
      allPlots["dilp_pt"+chTag+"_topsel_no_weight"]->Fill((leptons[0]+leptons[1]).Pt(),norm);
      allPlots["dilp_pt_all_topsel"]->Fill((leptons[0]+leptons[1]).Pt(),wgt);
      allPlots["dilp_pt_all_topsel_no_weight"]->Fill((leptons[0]+leptons[1]).Pt(),norm);
      allPlots["dilp_mass"+chTag+"_topsel"]->Fill((leptons[0]+leptons[1]).M(),wgt);
      allPlots["dilp_mass"+chTag+"_topsel_no_weight"]->Fill((leptons[0]+leptons[1]).M(),norm);
      allPlots["dilp_mass_all_topsel"]->Fill((leptons[0]+leptons[1]).M(),wgt);
      allPlots["dilp_mass_all_topsel_no_weight"]->Fill((leptons[0]+leptons[1]).M(),norm);
      allPlots["dilp_charge"+chTag+"_topsel"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
      allPlots["dilp_charge"+chTag+"_topsel_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
      allPlots["dilp_charge_all_topsel"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
      allPlots["dilp_charge_all_topsel_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
    }
    allPlots["j_n"+chTag+"_topsel"]->Fill(allJetsVec.size(),wgt);
    allPlots["j_n"+chTag+"_topsel_no_weight"]->Fill(allJetsVec.size(),norm);
    allPlots["j_n_all_topsel"]->Fill(allJetsVec.size(),wgt);
    allPlots["j_n_all_topsel_no_weight"]->Fill(allJetsVec.size(),norm);
    for (size_t iJet = 0; iJet < allJetsVec.size(); iJet++) {
      allPlots["j_pt"+chTag+"_topsel"]->Fill(allJetsVec[iJet].getVec().Pt(),wgt);
      allPlots["j_pt"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[iJet].getVec().Pt(),norm);
      allPlots["j_pt_all_topsel"]->Fill(allJetsVec[iJet].getVec().Pt(),wgt);
      allPlots["j_pt_all_topsel_no_weight"]->Fill(allJetsVec[iJet].getVec().Pt(),norm);
      allPlots["j_eta"+chTag+"_topsel"]->Fill(allJetsVec[iJet].getVec().Eta(),wgt);
      allPlots["j_eta"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[iJet].getVec().Eta(),norm);
      allPlots["j_eta_all_topsel"]->Fill(allJetsVec[iJet].getVec().Eta(),wgt);
      allPlots["j_eta_all_topsel_no_weight"]->Fill(allJetsVec[iJet].getVec().Eta(),norm);
      allPlots["j_phi"+chTag+"_topsel"]->Fill(allJetsVec[iJet].getVec().Phi(),wgt);
      allPlots["j_phi"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[iJet].getVec().Phi(),norm);
      allPlots["j_phi_all_topsel"]->Fill(allJetsVec[iJet].getVec().Phi(),wgt);
      allPlots["j_phi_all_topsel_no_weight"]->Fill(allJetsVec[iJet].getVec().Phi(),norm);
      allPlots["j_csv"+chTag+"_topsel"]->Fill(allJetsVec[iJet].getCSV(),wgt);
      allPlots["j_csv"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[iJet].getCSV(),norm);
      allPlots["j_csv_all_topsel"]->Fill(allJetsVec[iJet].getCSV(),wgt);
      allPlots["j_csv_all_topsel_no_weight"]->Fill(allJetsVec[iJet].getCSV(),norm);
      allPlots["j_nch"+chTag+"_topsel"]->Fill(allJetsVec[iJet].getTracks().size(),wgt);
      allPlots["j_nch"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[iJet].getTracks().size(),norm);
      allPlots["j_nch_all_topsel"]->Fill(allJetsVec[iJet].getTracks().size(),wgt);
      allPlots["j_nch_all_topsel_no_weight"]->Fill(allJetsVec[iJet].getTracks().size(),norm);
    }
    if (allJetsVec.size() > 0) {
      allPlots["j1_pt"+chTag+"_topsel"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      allPlots["j1_pt"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
      allPlots["j1_pt_all_topsel"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      allPlots["j1_pt_all_topsel_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
      allPlots["j1_eta"+chTag+"_topsel"]->Fill(allJetsVec[0].getVec().Eta(),wgt);
      allPlots["j1_eta"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[0].getVec().Eta(),norm);
      allPlots["j1_eta_all_topsel"]->Fill(allJetsVec[0].getVec().Eta(),wgt);
      allPlots["j1_eta_all_topsel_no_weight"]->Fill(allJetsVec[0].getVec().Eta(),norm);
      allPlots["j1_phi"+chTag+"_topsel"]->Fill(allJetsVec[0].getVec().Phi(),wgt);
      allPlots["j1_phi"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[0].getVec().Phi(),norm);
      allPlots["j1_phi_all_topsel"]->Fill(allJetsVec[0].getVec().Phi(),wgt);
      allPlots["j1_phi_all_topsel_no_weight"]->Fill(allJetsVec[0].getVec().Phi(),norm);
      allPlots["j1_csv"+chTag+"_topsel"]->Fill(allJetsVec[0].getCSV(),wgt);
      allPlots["j1_csv"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[0].getCSV(),norm);
      allPlots["j1_csv_all_topsel"]->Fill(allJetsVec[0].getCSV(),wgt);
      allPlots["j1_csv_all_topsel_no_weight"]->Fill(allJetsVec[0].getCSV(),norm);
      allPlots["j1_nch"+chTag+"_topsel"]->Fill(allJetsVec[0].getTracks().size(),wgt);
      allPlots["j1_nch"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[0].getTracks().size(),norm);
      allPlots["j1_nch_all_topsel"]->Fill(allJetsVec[0].getTracks().size(),wgt);
      allPlots["j1_nch_all_topsel_no_weight"]->Fill(allJetsVec[0].getTracks().size(),norm);
      if (allJetsVec.size() > 1) {
        allPlots["j2_pt"+chTag+"_topsel"]->Fill(allJetsVec[1].getVec().Pt(),wgt);
        allPlots["j2_pt"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[1].getVec().Pt(),norm);
        allPlots["j2_pt_all_topsel"]->Fill(allJetsVec[1].getVec().Pt(),wgt);
        allPlots["j2_pt_all_topsel_no_weight"]->Fill(allJetsVec[1].getVec().Pt(),norm);
        allPlots["j2_eta"+chTag+"_topsel"]->Fill(allJetsVec[1].getVec().Eta(),wgt);
        allPlots["j2_eta"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[1].getVec().Eta(),norm);
        allPlots["j2_eta_all_topsel"]->Fill(allJetsVec[1].getVec().Eta(),wgt);
        allPlots["j2_eta_all_topsel_no_weight"]->Fill(allJetsVec[1].getVec().Eta(),norm);
        allPlots["j2_phi"+chTag+"_topsel"]->Fill(allJetsVec[1].getVec().Phi(),wgt);
        allPlots["j2_phi"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[1].getVec().Phi(),norm);
        allPlots["j2_phi_all_topsel"]->Fill(allJetsVec[1].getVec().Phi(),wgt);
        allPlots["j2_phi_all_topsel_no_weight"]->Fill(allJetsVec[1].getVec().Phi(),norm);
        allPlots["j2_csv"+chTag+"_topsel"]->Fill(allJetsVec[1].getCSV(),wgt);
        allPlots["j2_csv"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[1].getCSV(),norm);
        allPlots["j2_csv_all_topsel"]->Fill(allJetsVec[1].getCSV(),wgt);
        allPlots["j2_csv_all_topsel_no_weight"]->Fill(allJetsVec[1].getCSV(),norm);
        allPlots["j2_nch"+chTag+"_topsel"]->Fill(allJetsVec[1].getTracks().size(),wgt);
        allPlots["j2_nch"+chTag+"_topsel_no_weight"]->Fill(allJetsVec[1].getTracks().size(),norm);
        allPlots["j2_nch_all_topsel"]->Fill(allJetsVec[1].getTracks().size(),wgt);
        allPlots["j2_nch_all_topsel_no_weight"]->Fill(allJetsVec[1].getTracks().size(),norm);
      }
    }
    allPlots["lj_n"+chTag+"_topsel"]->Fill(lightJetsVec.size(),wgt);
    allPlots["lj_n"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec.size(),norm);
    allPlots["lj_n_all_topsel"]->Fill(lightJetsVec.size(),wgt);
    allPlots["lj_n_all_topsel_no_weight"]->Fill(lightJetsVec.size(),norm);
    for (size_t iJet = 0; iJet < lightJetsVec.size(); iJet++) {
      allPlots["lj_pt"+chTag+"_topsel"]->Fill(lightJetsVec[iJet].getVec().Pt(),wgt);
      allPlots["lj_pt"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec[iJet].getVec().Pt(),norm);
      allPlots["lj_pt_all_topsel"]->Fill(lightJetsVec[iJet].getVec().Pt(),wgt);
      allPlots["lj_pt_all_topsel_no_weight"]->Fill(lightJetsVec[iJet].getVec().Pt(),norm);
      allPlots["lj_eta"+chTag+"_topsel"]->Fill(lightJetsVec[iJet].getVec().Eta(),wgt);
      allPlots["lj_eta"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec[iJet].getVec().Eta(),norm);
      allPlots["lj_eta_all_topsel"]->Fill(lightJetsVec[iJet].getVec().Eta(),wgt);
      allPlots["lj_eta_all_topsel_no_weight"]->Fill(lightJetsVec[iJet].getVec().Eta(),norm);
      allPlots["lj_phi"+chTag+"_topsel"]->Fill(lightJetsVec[iJet].getVec().Phi(),wgt);
      allPlots["lj_phi"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec[iJet].getVec().Phi(),norm);
      allPlots["lj_phi_all_topsel"]->Fill(lightJetsVec[iJet].getVec().Phi(),wgt);
      allPlots["lj_phi_all_topsel_no_weight"]->Fill(lightJetsVec[iJet].getVec().Phi(),norm);
      allPlots["lj_nch"+chTag+"_topsel"]->Fill(lightJetsVec[iJet].getTracks().size(),wgt);
      allPlots["lj_nch"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec[iJet].getTracks().size(),norm);
      allPlots["lj_nch_all_topsel"]->Fill(lightJetsVec[iJet].getTracks().size(),wgt);
      allPlots["lj_nch_all_topsel_no_weight"]->Fill(lightJetsVec[iJet].getTracks().size(),norm);
    }
    if (lightJetsVec.size() > 0) {
      allPlots["lj1_pt"+chTag+"_topsel"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
      allPlots["lj1_pt"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
      allPlots["lj1_pt_all_topsel"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
      allPlots["lj1_pt_all_topsel_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
      allPlots["lj1_eta"+chTag+"_topsel"]->Fill(lightJetsVec[0].getVec().Eta(),wgt);
      allPlots["lj1_eta"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec[0].getVec().Eta(),norm);
      allPlots["lj1_eta_all_topsel"]->Fill(lightJetsVec[0].getVec().Eta(),wgt);
      allPlots["lj1_eta_all_topsel_no_weight"]->Fill(lightJetsVec[0].getVec().Eta(),norm);
      allPlots["lj1_phi"+chTag+"_topsel"]->Fill(lightJetsVec[0].getVec().Phi(),wgt);
      allPlots["lj1_phi"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec[0].getVec().Phi(),norm);
      allPlots["lj1_phi_all_topsel"]->Fill(lightJetsVec[0].getVec().Phi(),wgt);
      allPlots["lj1_phi_all_topsel_no_weight"]->Fill(lightJetsVec[0].getVec().Phi(),norm);
      allPlots["lj1_nch"+chTag+"_topsel"]->Fill(lightJetsVec[0].getTracks().size(),wgt);
      allPlots["lj1_nch"+chTag+"_topsel_no_weight"]->Fill(lightJetsVec[0].getTracks().size(),norm);
      allPlots["lj1_nch_all_topsel"]->Fill(lightJetsVec[0].getTracks().size(),wgt);
      allPlots["lj1_nch_all_topsel_no_weight"]->Fill(lightJetsVec[0].getTracks().size(),norm);
    }
    allPlots["bj_n"+chTag+"_topsel"]->Fill(bJetsVec.size(),wgt);
    allPlots["bj_n"+chTag+"_topsel_no_weight"]->Fill(bJetsVec.size(),norm);
    allPlots["bj_n_all_topsel"]->Fill(bJetsVec.size(),wgt);
    allPlots["bj_n_all_topsel_no_weight"]->Fill(bJetsVec.size(),norm);
    for (size_t iJet = 0; iJet < bJetsVec.size(); iJet++) {
      allPlots["bj_pt"+chTag+"_topsel"]->Fill(bJetsVec[iJet].getVec().Pt(),wgt);
      allPlots["bj_pt"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[iJet].getVec().Pt(),norm);
      allPlots["bj_pt_all_topsel"]->Fill(bJetsVec[iJet].getVec().Pt(),wgt);
      allPlots["bj_pt_all_topsel_no_weight"]->Fill(bJetsVec[iJet].getVec().Pt(),norm);
      allPlots["bj_eta"+chTag+"_topsel"]->Fill(bJetsVec[iJet].getVec().Eta(),wgt);
      allPlots["bj_eta"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[iJet].getVec().Eta(),norm);
      allPlots["bj_eta_all_topsel"]->Fill(bJetsVec[iJet].getVec().Eta(),wgt);
      allPlots["bj_eta_all_topsel_no_weight"]->Fill(bJetsVec[iJet].getVec().Eta(),norm);
      allPlots["bj_phi"+chTag+"_topsel"]->Fill(bJetsVec[iJet].getVec().Phi(),wgt);
      allPlots["bj_phi"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[iJet].getVec().Phi(),norm);
      allPlots["bj_phi_all_topsel"]->Fill(bJetsVec[iJet].getVec().Phi(),wgt);
      allPlots["bj_phi_all_topsel_no_weight"]->Fill(bJetsVec[iJet].getVec().Phi(),norm);
      allPlots["bj_csv"+chTag+"_topsel"]->Fill(bJetsVec[iJet].getCSV(),wgt);
      allPlots["bj_csv"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[iJet].getCSV(),norm);
      allPlots["bj_csv_all_topsel"]->Fill(bJetsVec[iJet].getCSV(),wgt);
      allPlots["bj_csv_all_topsel_no_weight"]->Fill(bJetsVec[iJet].getCSV(),norm);
      allPlots["bj_nch"+chTag+"_topsel"]->Fill(bJetsVec[iJet].getTracks().size(),wgt);
      allPlots["bj_nch"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[iJet].getTracks().size(),norm);
      allPlots["bj_nch_all_topsel"]->Fill(bJetsVec[iJet].getTracks().size(),wgt);
      allPlots["bj_nch_all_topsel_no_weight"]->Fill(bJetsVec[iJet].getTracks().size(),norm);
    }
    if (bJetsVec.size() > 0) {
      allPlots["bj1_pt"+chTag+"_topsel"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
      allPlots["bj1_pt"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
      allPlots["bj1_pt_all_topsel"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
      allPlots["bj1_pt_all_topsel_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
      allPlots["bj1_eta"+chTag+"_topsel"]->Fill(bJetsVec[0].getVec().Eta(),wgt);
      allPlots["bj1_eta"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[0].getVec().Eta(),norm);
      allPlots["bj1_eta_all_topsel"]->Fill(bJetsVec[0].getVec().Eta(),wgt);
      allPlots["bj1_eta_all_topsel_no_weight"]->Fill(bJetsVec[0].getVec().Eta(),norm);
      allPlots["bj1_phi"+chTag+"_topsel"]->Fill(bJetsVec[0].getVec().Phi(),wgt);
      allPlots["bj1_phi"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[0].getVec().Phi(),norm);
      allPlots["bj1_phi_all_topsel"]->Fill(bJetsVec[0].getVec().Phi(),wgt);
      allPlots["bj1_phi_all_topsel_no_weight"]->Fill(bJetsVec[0].getVec().Phi(),norm);
      allPlots["bj1_csv"+chTag+"_topsel"]->Fill(bJetsVec[0].getCSV(),wgt);
      allPlots["bj1_csv"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[0].getCSV(),norm);
      allPlots["bj1_csv_all_topsel"]->Fill(bJetsVec[0].getCSV(),wgt);
      allPlots["bj1_csv_all_topsel_no_weight"]->Fill(bJetsVec[0].getCSV(),norm);
      allPlots["bj1_nch"+chTag+"_topsel"]->Fill(bJetsVec[0].getTracks().size(),wgt);
      allPlots["bj1_nch"+chTag+"_topsel_no_weight"]->Fill(bJetsVec[0].getTracks().size(),norm);
      allPlots["bj1_nch_all_topsel"]->Fill(bJetsVec[0].getTracks().size(),wgt);
      allPlots["bj1_nch_all_topsel_no_weight"]->Fill(bJetsVec[0].getTracks().size(),norm);
    }

    if (debug) cout << "topsel plots DONE" << endl;

    const float gMassMu(0.1057),gMassK(0.4937),gMassPi(0.1396);    

    // jpsi resonance analysis : use all jets
    if (debug) cout << "Starting J/Psi" << endl;
    int nJpsi = 0;
    for (size_t ij = 0; ij < allJetsVec.size(); ij++) {
      //J/Psi
      std::vector<IdTrack> &tracks = allJetsVec[ij].getTracks(); 

      for (size_t itk1 = 0; itk1 < tracks.size(); itk1++) {
        if (abs(tracks[itk1].second) != 13) continue;
        if (tracks[itk1].first.Pt() < 4.) continue; // FIXME : to be optimized
        for (size_t itk2 = itk1+1; itk2 < tracks.size(); itk2++) {
          if (abs(tracks[itk2].second) != 13) continue;
          if (tracks[itk2].second * tracks[itk1].second > 0) continue;
          if (tracks[itk2].first.Pt() < 4.) continue; // FIXME : to be optimized
          if (nJpsi > 1) continue;

          TLorentzVector mu1P4, mu2P4;
          mu1P4.SetPtEtaPhiM(tracks[itk1].first.Pt(), tracks[itk1].first.Eta(), tracks[itk1].first.Phi(), gMassMu);
          mu2P4.SetPtEtaPhiM(tracks[itk2].first.Pt(), tracks[itk2].first.Eta(), tracks[itk2].first.Phi(), gMassMu);

          float mass12((mu1P4+mu2P4).M());

          if (mass12 > 3.6 || mass12 < 2.6) continue;
          allPlots["jpsi_mass"+chTag+"_check"]->Fill(mass12,wgt);
          allPlots["jpsi_mass"+chTag+"_check_no_weight"]->Fill(mass12,norm);
          allPlots["jpsi_mass_all_check"]->Fill(mass12,wgt);
          allPlots["jpsi_mass_all_check_no_weight"]->Fill(mass12,norm);
          if (mass12 > 3.3 || mass12 < 2.9) continue;
          ++nJpsi;

          if (debug) cout << "starting jpsicand plots" << endl;

          if (chTag == "_e") {
            allPlots["decay_all_jpsicand"]->Fill(0.,wgt);
            allPlots["decay_all_jpsicand_no_weight"]->Fill(0.,norm);
          }
          if (chTag == "_m") {
            allPlots["decay_all_jpsicand"]->Fill(1.,wgt);
            allPlots["decay_all_jpsicand_no_weight"]->Fill(1.,norm);
          }
          if (chTag == "_ee") {
            allPlots["decay_all_jpsicand"]->Fill(2.,wgt);
            allPlots["decay_all_jpsicand_no_weight"]->Fill(2.,norm);
          }
          if (chTag == "_em") {
            allPlots["decay_all_jpsicand"]->Fill(3.,wgt);
            allPlots["decay_all_jpsicand_no_weight"]->Fill(3.,norm);
          }
          if (chTag == "_mm") {
            allPlots["decay_all_jpsicand"]->Fill(4.,wgt);
            allPlots["decay_all_jpsicand_no_weight"]->Fill(4.,norm);
          }

          allPlots["lp_n"+chTag+"_jpsicand"]->Fill(selLeptons.size(),wgt);
          allPlots["lp_n"+chTag+"_jpsicand_no_weight"]->Fill(selLeptons.size(),norm);
          allPlots["lp_n_all_jpsicand"]->Fill(selLeptons.size(),wgt);
          allPlots["lp_n_all_jpsicand_no_weight"]->Fill(selLeptons.size(),norm);
          allPlots["lp1_pt"+chTag+"_jpsicand"]->Fill(leptons[0].Pt(),wgt);
          allPlots["lp1_pt"+chTag+"_jpsicand_no_weight"]->Fill(leptons[0].Pt(),norm);
          allPlots["lp1_pt_all_jpsicand"]->Fill(leptons[0].Pt(),wgt);
          allPlots["lp1_pt_all_jpsicand_no_weight"]->Fill(leptons[0].Pt(),norm);
          allPlots["lp1_eta"+chTag+"_jpsicand"]->Fill(leptons[0].Eta(),wgt);
          allPlots["lp1_eta"+chTag+"_jpsicand_no_weight"]->Fill(leptons[0].Eta(),norm);
          allPlots["lp1_eta_all_jpsicand"]->Fill(leptons[0].Eta(),wgt);
          allPlots["lp1_eta_all_jpsicand_no_weight"]->Fill(leptons[0].Eta(),norm);
          allPlots["lp1_phi"+chTag+"_jpsicand"]->Fill(leptons[0].Phi(),wgt);
          allPlots["lp1_phi"+chTag+"_jpsicand_no_weight"]->Fill(leptons[0].Phi(),norm);
          allPlots["lp1_phi_all_jpsicand"]->Fill(leptons[0].Phi(),wgt);
          allPlots["lp1_phi_all_jpsicand_no_weight"]->Fill(leptons[0].Phi(),norm);          
          if (selLeptons.size() == 2) {
            allPlots["lp2_pt"+chTag+"_jpsicand"]->Fill(leptons[1].Pt(),wgt);
            allPlots["lp2_pt"+chTag+"_jpsicand_no_weight"]->Fill(leptons[1].Pt(),norm);
            allPlots["lp2_pt_all_jpsicand"]->Fill(leptons[1].Pt(),wgt);
            allPlots["lp2_pt_all_jpsicand_no_weight"]->Fill(leptons[1].Pt(),norm);
            allPlots["lp2_eta"+chTag+"_jpsicand"]->Fill(leptons[1].Eta(),wgt);
            allPlots["lp2_eta"+chTag+"_jpsicand_no_weight"]->Fill(leptons[1].Eta(),norm);
            allPlots["lp2_eta_all_jpsicand"]->Fill(leptons[1].Eta(),wgt);
            allPlots["lp2_eta_all_jpsicand_no_weight"]->Fill(leptons[1].Eta(),norm);
            allPlots["lp2_phi"+chTag+"_jpsicand"]->Fill(leptons[1].Phi(),wgt);
            allPlots["lp2_phi"+chTag+"_jpsicand_no_weight"]->Fill(leptons[1].Phi(),norm);
            allPlots["lp2_phi_all_jpsicand"]->Fill(leptons[1].Phi(),wgt);
            allPlots["lp2_phi_all_jpsicand_no_weight"]->Fill(leptons[1].Phi(),norm);            
          }
          allPlots["met"+chTag+"_jpsicand"]->Fill(met.Pt(),wgt);
          allPlots["met"+chTag+"_jpsicand_no_weight"]->Fill(met.Pt(),norm);
          allPlots["met_all_jpsicand"]->Fill(met.Pt(),wgt);
          allPlots["met_all_jpsicand_no_weight"]->Fill(met.Pt(),norm);
          if (selLeptons.size() == 2) {
            allPlots["dilp_pt"+chTag+"_jpsicand"]->Fill((leptons[0]+leptons[1]).Pt(),wgt);
            allPlots["dilp_pt"+chTag+"_jpsicand_no_weight"]->Fill((leptons[0]+leptons[1]).Pt(),norm);
            allPlots["dilp_pt_all_jpsicand"]->Fill((leptons[0]+leptons[1]).Pt(),wgt);
            allPlots["dilp_pt_all_jpsicand_no_weight"]->Fill((leptons[0]+leptons[1]).Pt(),norm);
            allPlots["dilp_mass"+chTag+"_jpsicand"]->Fill((leptons[0]+leptons[1]).M(),wgt);
            allPlots["dilp_mass"+chTag+"_jpsicand_no_weight"]->Fill((leptons[0]+leptons[1]).M(),norm);
            allPlots["dilp_mass_all_jpsicand"]->Fill((leptons[0]+leptons[1]).M(),wgt);
            allPlots["dilp_mass_all_jpsicand_no_weight"]->Fill((leptons[0]+leptons[1]).M(),norm);
            allPlots["dilp_charge"+chTag+"_jpsicand"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
            allPlots["dilp_charge"+chTag+"_jpsicand_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
            allPlots["dilp_charge_all_jpsicand"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
            allPlots["dilp_charge_all_jpsicand_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
          }
          allPlots["j_n"+chTag+"_jpsicand"]->Fill(allJetsVec.size(),wgt);
          allPlots["j_n"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec.size(),norm);
          allPlots["j_n_all_jpsicand"]->Fill(allJetsVec.size(),wgt);
          allPlots["j_n_all_jpsicand_no_weight"]->Fill(allJetsVec.size(),norm);
          for (size_t iJet = 0; iJet < allJetsVec.size(); iJet++) {
            allPlots["j_pt"+chTag+"_jpsicand"]->Fill(allJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["j_pt"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getVec().Pt(),norm);
            allPlots["j_pt_all_jpsicand"]->Fill(allJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["j_pt_all_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getVec().Pt(),norm);
            allPlots["j_eta"+chTag+"_jpsicand"]->Fill(allJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["j_eta"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getVec().Eta(),norm);
            allPlots["j_eta_all_jpsicand"]->Fill(allJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["j_eta_all_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getVec().Eta(),norm);
            allPlots["j_phi"+chTag+"_jpsicand"]->Fill(allJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["j_phi"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getVec().Phi(),norm);
            allPlots["j_phi_all_jpsicand"]->Fill(allJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["j_phi_all_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getVec().Phi(),norm);
            allPlots["j_csv"+chTag+"_jpsicand"]->Fill(allJetsVec[iJet].getCSV(),wgt);
            allPlots["j_csv"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getCSV(),norm);
            allPlots["j_csv_all_jpsicand"]->Fill(allJetsVec[iJet].getCSV(),wgt);
            allPlots["j_csv_all_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getCSV(),norm);
            allPlots["j_nch"+chTag+"_jpsicand"]->Fill(allJetsVec[iJet].getTracks().size(),wgt);
            allPlots["j_nch"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getTracks().size(),norm);
            allPlots["j_nch_all_jpsicand"]->Fill(allJetsVec[iJet].getTracks().size(),wgt);
            allPlots["j_nch_all_jpsicand_no_weight"]->Fill(allJetsVec[iJet].getTracks().size(),norm);
          }
          if (allJetsVec.size() > 0) {
            allPlots["j1_pt"+chTag+"_jpsicand"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
            allPlots["j1_pt"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
            allPlots["j1_pt_all_jpsicand"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
            allPlots["j1_pt_all_jpsicand_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
            allPlots["j1_eta"+chTag+"_jpsicand"]->Fill(allJetsVec[0].getVec().Eta(),wgt);
            allPlots["j1_eta"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[0].getVec().Eta(),norm);
            allPlots["j1_eta_all_jpsicand"]->Fill(allJetsVec[0].getVec().Eta(),wgt);
            allPlots["j1_eta_all_jpsicand_no_weight"]->Fill(allJetsVec[0].getVec().Eta(),norm);
            allPlots["j1_phi"+chTag+"_jpsicand"]->Fill(allJetsVec[0].getVec().Phi(),wgt);
            allPlots["j1_phi"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[0].getVec().Phi(),norm);
            allPlots["j1_phi_all_jpsicand"]->Fill(allJetsVec[0].getVec().Phi(),wgt);
            allPlots["j1_phi_all_jpsicand_no_weight"]->Fill(allJetsVec[0].getVec().Phi(),norm);
            allPlots["j1_csv"+chTag+"_jpsicand"]->Fill(allJetsVec[0].getCSV(),wgt);
            allPlots["j1_csv"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[0].getCSV(),norm);
            allPlots["j1_csv_all_jpsicand"]->Fill(allJetsVec[0].getCSV(),wgt);
            allPlots["j1_csv_all_jpsicand_no_weight"]->Fill(allJetsVec[0].getCSV(),norm);
            allPlots["j1_nch"+chTag+"_jpsicand"]->Fill(allJetsVec[0].getTracks().size(),wgt);
            allPlots["j1_nch"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[0].getTracks().size(),norm);
            allPlots["j1_nch_all_jpsicand"]->Fill(allJetsVec[0].getTracks().size(),wgt);
            allPlots["j1_nch_all_jpsicand_no_weight"]->Fill(allJetsVec[0].getTracks().size(),norm);
            if (allJetsVec.size() > 1) {
              allPlots["j2_pt"+chTag+"_jpsicand"]->Fill(allJetsVec[1].getVec().Pt(),wgt);
              allPlots["j2_pt"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[1].getVec().Pt(),norm);
              allPlots["j2_pt_all_jpsicand"]->Fill(allJetsVec[1].getVec().Pt(),wgt);
              allPlots["j2_pt_all_jpsicand_no_weight"]->Fill(allJetsVec[1].getVec().Pt(),norm);
              allPlots["j2_eta"+chTag+"_jpsicand"]->Fill(allJetsVec[1].getVec().Eta(),wgt);
              allPlots["j2_eta"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[1].getVec().Eta(),norm);
              allPlots["j2_eta_all_jpsicand"]->Fill(allJetsVec[1].getVec().Eta(),wgt);
              allPlots["j2_eta_all_jpsicand_no_weight"]->Fill(allJetsVec[1].getVec().Eta(),norm);
              allPlots["j2_phi"+chTag+"_jpsicand"]->Fill(allJetsVec[1].getVec().Phi(),wgt);
              allPlots["j2_phi"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[1].getVec().Phi(),norm);
              allPlots["j2_phi_all_jpsicand"]->Fill(allJetsVec[1].getVec().Phi(),wgt);
              allPlots["j2_phi_all_jpsicand_no_weight"]->Fill(allJetsVec[1].getVec().Phi(),norm);
              allPlots["j2_csv"+chTag+"_jpsicand"]->Fill(allJetsVec[1].getCSV(),wgt);
              allPlots["j2_csv"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[1].getCSV(),norm);
              allPlots["j2_csv_all_jpsicand"]->Fill(allJetsVec[1].getCSV(),wgt);
              allPlots["j2_csv_all_jpsicand_no_weight"]->Fill(allJetsVec[1].getCSV(),norm);
              allPlots["j2_nch"+chTag+"_jpsicand"]->Fill(allJetsVec[1].getTracks().size(),wgt);
              allPlots["j2_nch"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[1].getTracks().size(),norm);
              allPlots["j2_nch_all_jpsicand"]->Fill(allJetsVec[1].getTracks().size(),wgt);
              allPlots["j2_nch_all_jpsicand_no_weight"]->Fill(allJetsVec[1].getTracks().size(),norm);
            }
          }
          allPlots["lj_n"+chTag+"_jpsicand"]->Fill(lightJetsVec.size(),wgt);
          allPlots["lj_n"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec.size(),norm);
          allPlots["lj_n_all_jpsicand"]->Fill(lightJetsVec.size(),wgt);
          allPlots["lj_n_all_jpsicand_no_weight"]->Fill(lightJetsVec.size(),norm);
          for (size_t iJet = 0; iJet < lightJetsVec.size(); iJet++) {
            allPlots["lj_pt"+chTag+"_jpsicand"]->Fill(lightJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["lj_pt"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Pt(),norm);
            allPlots["lj_pt_all_jpsicand"]->Fill(lightJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["lj_pt_all_jpsicand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Pt(),norm);
            allPlots["lj_eta"+chTag+"_jpsicand"]->Fill(lightJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["lj_eta"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Eta(),norm);
            allPlots["lj_eta_all_jpsicand"]->Fill(lightJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["lj_eta_all_jpsicand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Eta(),norm);
            allPlots["lj_phi"+chTag+"_jpsicand"]->Fill(lightJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["lj_phi"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Phi(),norm);
            allPlots["lj_phi_all_jpsicand"]->Fill(lightJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["lj_phi_all_jpsicand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Phi(),norm);
            allPlots["lj_nch"+chTag+"_jpsicand"]->Fill(lightJetsVec[iJet].getTracks().size(),wgt);
            allPlots["lj_nch"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec[iJet].getTracks().size(),norm);
            allPlots["lj_nch_all_jpsicand"]->Fill(lightJetsVec[iJet].getTracks().size(),wgt);
            allPlots["lj_nch_all_jpsicand_no_weight"]->Fill(lightJetsVec[iJet].getTracks().size(),norm);
          }
          if (lightJetsVec.size() > 0) {
            allPlots["lj1_pt"+chTag+"_jpsicand"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
            allPlots["lj1_pt"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
            allPlots["lj1_pt_all_jpsicand"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
            allPlots["lj1_pt_all_jpsicand_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
            allPlots["lj1_eta"+chTag+"_jpsicand"]->Fill(lightJetsVec[0].getVec().Eta(),wgt);
            allPlots["lj1_eta"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec[0].getVec().Eta(),norm);
            allPlots["lj1_eta_all_jpsicand"]->Fill(lightJetsVec[0].getVec().Eta(),wgt);
            allPlots["lj1_eta_all_jpsicand_no_weight"]->Fill(lightJetsVec[0].getVec().Eta(),norm);
            allPlots["lj1_phi"+chTag+"_jpsicand"]->Fill(lightJetsVec[0].getVec().Phi(),wgt);
            allPlots["lj1_phi"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec[0].getVec().Phi(),norm);
            allPlots["lj1_phi_all_jpsicand"]->Fill(lightJetsVec[0].getVec().Phi(),wgt);
            allPlots["lj1_phi_all_jpsicand_no_weight"]->Fill(lightJetsVec[0].getVec().Phi(),norm);
            allPlots["lj1_nch"+chTag+"_jpsicand"]->Fill(lightJetsVec[0].getTracks().size(),wgt);
            allPlots["lj1_nch"+chTag+"_jpsicand_no_weight"]->Fill(lightJetsVec[0].getTracks().size(),norm);
            allPlots["lj1_nch_all_jpsicand"]->Fill(lightJetsVec[0].getTracks().size(),wgt);
            allPlots["lj1_nch_all_jpsicand_no_weight"]->Fill(lightJetsVec[0].getTracks().size(),norm);
          }
          allPlots["bj_n"+chTag+"_jpsicand"]->Fill(bJetsVec.size(),wgt);
          allPlots["bj_n"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec.size(),norm);
          allPlots["bj_n_all_jpsicand"]->Fill(bJetsVec.size(),wgt);
          allPlots["bj_n_all_jpsicand_no_weight"]->Fill(bJetsVec.size(),norm);
          for (size_t iJet = 0; iJet < bJetsVec.size(); iJet++) {
            allPlots["bj_pt"+chTag+"_jpsicand"]->Fill(bJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["bj_pt"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getVec().Pt(),norm);
            allPlots["bj_pt_all_jpsicand"]->Fill(bJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["bj_pt_all_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getVec().Pt(),norm);
            allPlots["bj_eta"+chTag+"_jpsicand"]->Fill(bJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["bj_eta"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getVec().Eta(),norm);
            allPlots["bj_eta_all_jpsicand"]->Fill(bJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["bj_eta_all_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getVec().Eta(),norm);
            allPlots["bj_phi"+chTag+"_jpsicand"]->Fill(bJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["bj_phi"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getVec().Phi(),norm);
            allPlots["bj_phi_all_jpsicand"]->Fill(bJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["bj_phi_all_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getVec().Phi(),norm);
            allPlots["bj_csv"+chTag+"_jpsicand"]->Fill(bJetsVec[iJet].getCSV(),wgt);
            allPlots["bj_csv"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getCSV(),norm);
            allPlots["bj_csv_all_jpsicand"]->Fill(bJetsVec[iJet].getCSV(),wgt);
            allPlots["bj_csv_all_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getCSV(),norm);
            allPlots["bj_nch"+chTag+"_jpsicand"]->Fill(bJetsVec[iJet].getTracks().size(),wgt);
            allPlots["bj_nch"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getTracks().size(),norm);
            allPlots["bj_nch_all_jpsicand"]->Fill(bJetsVec[iJet].getTracks().size(),wgt);
            allPlots["bj_nch_all_jpsicand_no_weight"]->Fill(bJetsVec[iJet].getTracks().size(),norm);
          }
          if (bJetsVec.size() > 0) {
            allPlots["bj1_pt"+chTag+"_jpsicand"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
            allPlots["bj1_pt"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
            allPlots["bj1_pt_all_jpsicand"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
            allPlots["bj1_pt_all_jpsicand_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
            allPlots["bj1_eta"+chTag+"_jpsicand"]->Fill(bJetsVec[0].getVec().Eta(),wgt);
            allPlots["bj1_eta"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[0].getVec().Eta(),norm);
            allPlots["bj1_eta_all_jpsicand"]->Fill(bJetsVec[0].getVec().Eta(),wgt);
            allPlots["bj1_eta_all_jpsicand_no_weight"]->Fill(bJetsVec[0].getVec().Eta(),norm);
            allPlots["bj1_phi"+chTag+"_jpsicand"]->Fill(bJetsVec[0].getVec().Phi(),wgt);
            allPlots["bj1_phi"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[0].getVec().Phi(),norm);
            allPlots["bj1_phi_all_jpsicand"]->Fill(bJetsVec[0].getVec().Phi(),wgt);
            allPlots["bj1_phi_all_jpsicand_no_weight"]->Fill(bJetsVec[0].getVec().Phi(),norm);
            allPlots["bj1_csv"+chTag+"_jpsicand"]->Fill(bJetsVec[0].getCSV(),wgt);
            allPlots["bj1_csv"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[0].getCSV(),norm);
            allPlots["bj1_csv_all_jpsicand"]->Fill(bJetsVec[0].getCSV(),wgt);
            allPlots["bj1_csv_all_jpsicand_no_weight"]->Fill(bJetsVec[0].getCSV(),norm);
            allPlots["bj1_nch"+chTag+"_jpsicand"]->Fill(bJetsVec[0].getTracks().size(),wgt);
            allPlots["bj1_nch"+chTag+"_jpsicand_no_weight"]->Fill(bJetsVec[0].getTracks().size(),norm);
            allPlots["bj1_nch_all_jpsicand"]->Fill(bJetsVec[0].getTracks().size(),wgt);
            allPlots["bj1_nch_all_jpsicand_no_weight"]->Fill(bJetsVec[0].getTracks().size(),norm);
          }

          allPlots["jpsi_lp_dR"+chTag+"_jpsicand"]->Fill((mu1P4+mu2P4).DeltaR(leptons[0]),wgt);
          allPlots["jpsi_lp_dR"+chTag+"_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).DeltaR(leptons[0]),norm);
          allPlots["jpsi_lp_dR_all_jpsicand"]->Fill((mu1P4+mu2P4).DeltaR(leptons[0]),wgt);
          allPlots["jpsi_lp_dR_all_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).DeltaR(leptons[0]),norm);
          allPlots["jpsi_lp_mass"+chTag+"_jpsicand"]->Fill((mu1P4+mu2P4+leptons[0]).M(),wgt);
          allPlots["jpsi_lp_mass"+chTag+"_jpsicand_no_weight"]->Fill((mu1P4+mu2P4+leptons[0]).M(),norm);
          allPlots["jpsi_lp_mass_all_jpsicand"]->Fill((mu1P4+mu2P4+leptons[0]).M(),wgt);
          allPlots["jpsi_lp_mass_all_jpsicand_no_weight"]->Fill((mu1P4+mu2P4+leptons[0]).M(),norm);

          allPlots["jpsi_mass"+chTag+"_jpsicand"]->Fill(mass12,wgt);
          allPlots["jpsi_mass"+chTag+"_jpsicand_no_weight"]->Fill(mass12,norm);
          allPlots["jpsi_mass_all_jpsicand"]->Fill(mass12,wgt);
          allPlots["jpsi_mass_all_jpsicand_no_weight"]->Fill(mass12,norm);
          allPlots["jpsi_ptfrac"+chTag+"_jpsicand"]->Fill((mu1P4+mu2P4).Pt()/allJetsVec[ij].getVec().Pt(),wgt);
          allPlots["jpsi_ptfrac"+chTag+"_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).Pt()/allJetsVec[ij].getVec().Pt(),norm);
          allPlots["jpsi_ptfrac_all_jpsicand"]->Fill((mu1P4+mu2P4).Pt()/allJetsVec[ij].getVec().Pt(),wgt);
          allPlots["jpsi_ptfrac_all_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).Pt()/allJetsVec[ij].getVec().Pt(),norm);
          allPlots["jpsi_pt"+chTag+"_jpsicand"]->Fill((mu1P4+mu2P4).Pt(),wgt);
          allPlots["jpsi_pt"+chTag+"_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).Pt(),norm);
          allPlots["jpsi_pt_all_jpsicand"]->Fill((mu1P4+mu2P4).Pt(),wgt);
          allPlots["jpsi_pt_all_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).Pt(),norm);
          allPlots["jpsi_eta"+chTag+"_jpsicand"]->Fill((mu1P4+mu2P4).Eta(),wgt);
          allPlots["jpsi_eta"+chTag+"_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).Eta(),norm);
          allPlots["jpsi_eta_all_jpsicand"]->Fill((mu1P4+mu2P4).Eta(),wgt);
          allPlots["jpsi_eta_all_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).Eta(),norm);
          allPlots["jpsi_phi"+chTag+"_jpsicand"]->Fill((mu1P4+mu2P4).Phi(),wgt);
          allPlots["jpsi_phi"+chTag+"_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).Phi(),norm);
          allPlots["jpsi_phi_all_jpsicand"]->Fill((mu1P4+mu2P4).Phi(),wgt);
          allPlots["jpsi_phi_all_jpsicand_no_weight"]->Fill((mu1P4+mu2P4).Phi(),norm);
          allPlots["mu_jpsi_pt"+chTag+"_jpsicand_no_weight"]->Fill(mu1P4.Pt(),norm);
          allPlots["mu_jpsi_pt_all_jpsicand"]->Fill(mu1P4.Pt(),wgt);
          allPlots["mu_jpsi_pt_all_jpsicand_no_weight"]->Fill(mu1P4.Pt(),norm);
          allPlots["mu_jpsi_eta"+chTag+"_jpsicand"]->Fill(mu1P4.Eta(),wgt);
          allPlots["mu_jpsi_eta"+chTag+"_jpsicand_no_weight"]->Fill(mu1P4.Eta(),norm);
          allPlots["mu_jpsi_eta_all_jpsicand"]->Fill(mu1P4.Eta(),wgt);
          allPlots["mu_jpsi_eta_all_jpsicand_no_weight"]->Fill(mu1P4.Eta(),norm);
          allPlots["mu_jpsi_phi"+chTag+"_jpsicand"]->Fill(mu1P4.Phi(),wgt);
          allPlots["mu_jpsi_phi"+chTag+"_jpsicand_no_weight"]->Fill(mu1P4.Phi(),norm);
          allPlots["mu_jpsi_phi_all_jpsicand"]->Fill(mu1P4.Phi(),wgt);
          allPlots["mu_jpsi_phi_all_jpsicand_no_weight"]->Fill(mu1P4.Phi(),norm);
          allPlots["mu_jpsi_pt"+chTag+"_jpsicand_no_weight"]->Fill(mu2P4.Pt(),norm);
          allPlots["mu_jpsi_pt_all_jpsicand"]->Fill(mu2P4.Pt(),wgt);
          allPlots["mu_jpsi_pt_all_jpsicand_no_weight"]->Fill(mu2P4.Pt(),norm);
          allPlots["mu_jpsi_eta"+chTag+"_jpsicand"]->Fill(mu2P4.Eta(),wgt);
          allPlots["mu_jpsi_eta"+chTag+"_jpsicand_no_weight"]->Fill(mu2P4.Eta(),norm);
          allPlots["mu_jpsi_eta_all_jpsicand"]->Fill(mu2P4.Eta(),wgt);
          allPlots["mu_jpsi_eta_all_jpsicand_no_weight"]->Fill(mu2P4.Eta(),norm);
          allPlots["mu_jpsi_phi"+chTag+"_jpsicand"]->Fill(mu2P4.Phi(),wgt);
          allPlots["mu_jpsi_phi"+chTag+"_jpsicand_no_weight"]->Fill(mu2P4.Phi(),norm);
          allPlots["mu_jpsi_phi_all_jpsicand"]->Fill(mu2P4.Phi(),wgt);
          allPlots["mu_jpsi_phi_all_jpsicand_no_weight"]->Fill(mu2P4.Phi(),norm);
          allPlots["mu_jpsi_dR"+chTag+"_jpsicand"]->Fill(mu1P4.DeltaR(mu2P4),wgt);
          allPlots["mu_jpsi_dR"+chTag+"_jpsicand_no_weight"]->Fill(mu1P4.DeltaR(mu2P4),norm);
          allPlots["mu_jpsi_dR_all_jpsicand"]->Fill(mu1P4.DeltaR(mu2P4),wgt);
          allPlots["mu_jpsi_dR_all_jpsicand_no_weight"]->Fill(mu1P4.DeltaR(mu2P4),norm);
          allPlots["j_jpsi_pt"+chTag+"_jpsicand"]->Fill(allJetsVec[ij].getVec().Pt(),wgt);
          allPlots["j_jpsi_pt"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[ij].getVec().Pt(),norm);
          allPlots["j_jpsi_pt_all_jpsicand"]->Fill(allJetsVec[ij].getVec().Pt(),wgt);
          allPlots["j_jpsi_pt_all_jpsicand_no_weight"]->Fill(allJetsVec[ij].getVec().Pt(),norm);
          allPlots["j_jpsi_eta"+chTag+"_jpsicand"]->Fill(allJetsVec[ij].getVec().Eta(),wgt);
          allPlots["j_jpsi_eta"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[ij].getVec().Eta(),norm);
          allPlots["j_jpsi_eta_all_jpsicand"]->Fill(allJetsVec[ij].getVec().Eta(),wgt);
          allPlots["j_jpsi_eta_all_jpsicand_no_weight"]->Fill(allJetsVec[ij].getVec().Eta(),norm);
          allPlots["j_jpsi_phi"+chTag+"_jpsicand"]->Fill(allJetsVec[ij].getVec().Phi(),wgt);
          allPlots["j_jpsi_phi"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[ij].getVec().Phi(),norm);
          allPlots["j_jpsi_phi_all_jpsicand"]->Fill(allJetsVec[ij].getVec().Phi(),wgt);
          allPlots["j_jpsi_phi_all_jpsicand_no_weight"]->Fill(allJetsVec[ij].getVec().Phi(),norm);
          allPlots["j_jpsi_csv"+chTag+"_jpsicand"]->Fill(allJetsVec[ij].getCSV(),wgt);
          allPlots["j_jpsi_csv"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[ij].getCSV(),norm);
          allPlots["j_jpsi_csv_all_jpsicand"]->Fill(allJetsVec[ij].getCSV(),wgt);
          allPlots["j_jpsi_csv_all_jpsicand_no_weight"]->Fill(allJetsVec[ij].getCSV(),norm);
          allPlots["j_jpsi_nch"+chTag+"_jpsicand"]->Fill(allJetsVec[ij].getTracks().size(),wgt);
          allPlots["j_jpsi_nch"+chTag+"_jpsicand_no_weight"]->Fill(allJetsVec[ij].getTracks().size(),norm);
          allPlots["j_jpsi_nch_all_jpsicand"]->Fill(allJetsVec[ij].getTracks().size(),wgt);
          allPlots["j_jpsi_nch_all_jpsicand_no_weight"]->Fill(allJetsVec[ij].getTracks().size(),norm);

        }
        if (debug) cout << "jpsicand plots DONE" << endl;
      }

    }
    if (nJpsi > 0) {
      allPlots["jpsi_n"+chTag+"_jpsicand"]->Fill(nJpsi,wgt);
      allPlots["jpsi_n"+chTag+"_jpsicand_no_weight"]->Fill(nJpsi,norm);
      allPlots["jpsi_n_all_jpsicand"]->Fill(nJpsi,wgt);
      allPlots["jpsi_n_all_jpsicand_no_weight"]->Fill(nJpsi,norm);
    }
    if (debug) cout << "J/Psi DONE -- # jpsicand = " << nJpsi << endl;      

    // D0 resonance analysis : use 2 jets of highest CSV, but more than 0.8
    if (debug) cout << "Starting D0 and D*" << endl;
    const size_t nJets = 2;
    int indJet[nJets];
    double csvJet[nJets];
    for (size_t k = 0; k < nJets; k++) {
      indJet[k] = -1;
      csvJet[k] = 0.;
    }
    for (size_t ij = 0; ij < allJetsVec.size(); ij++) {
      if (allJetsVec[ij].getCSV() < 0.8) continue;
      if (allJetsVec[ij].getCSV() > csvJet[0]) {
        for (size_t k = nJets-1; nJets > 0; k--) {
          csvJet[k] = csvJet[k-1];
          indJet[k] = indJet[k-1];
        }
        csvJet[0] = allJetsVec.size();
        indJet[0] = ij;
      } else if (allJetsVec[ij].getCSV() > csvJet[1]) {
        csvJet[1] = allJetsVec.size();
        indJet[1] = ij;
      }
    }
    if (debug) {
      cout << "Using " << nJets << " jets for D0 and D* reconstruction";
      for (size_t l = 0; l < nJets; l++)
        cout << " -- idx = " << indJet[l] << ", CSV = " << csvJet[l] << ", p_T  = " << allJetsVec[indJet[l]].getVec().Pt() ;
      cout << endl;
    }

    int nD0 = 0;
    for (size_t l = 0; l < nJets; l++) {
      if (indJet[l] < 0) continue;

      std::vector<IdTrack> &tracks = allJetsVec[indJet[l]].getTracks();

      //D0 and D* 
      if (tracks.size() < 3) continue;
      size_t tmax = 4;
      tmax = tracks.size() >= tmax ? tmax : tracks.size();
      for (size_t i = 0; i < tmax; i++) {
        if (fabs(tracks[i].second) != 211) continue;
        if (tracks[i].first.Pt() < 4.) continue; // FIXME : to be optimized
        for (size_t j = 0; j < tmax; j++) {
          if (i == j) continue;
          if (fabs(tracks[j].second) != 211) continue;
          if (tracks[i].second*tracks[j].second > 0) continue;
          if (tracks[j].first.Pt() < 4.) continue; // FIXME : to be optimized

          TLorentzVector p_track1, p_track2;
          p_track1.SetPtEtaPhiM(tracks[i].first.Pt(), tracks[i].first.Eta(), tracks[i].first.Phi(), gMassPi);
          p_track2.SetPtEtaPhiM(tracks[j].first.Pt(), tracks[j].first.Eta(), tracks[j].first.Phi(), gMassK);
          if (debug) cout << i << ": " << tracks[i].first.Pt() << " " << tracks[i].first.Eta() << " " << tracks[i].first.Phi() << " " << gMassPi << endl;
          if (debug) cout << j << ": " << tracks[j].first.Pt() << " " << tracks[j].first.Eta() << " " << tracks[j].first.Phi() << " " << gMassK << endl << endl;

          if ((p_track1+p_track2).Pt() < 10.) continue; // FIXME : to be optimized

          float mass12 = (p_track1+p_track2).M();
          if (debug) cout << mass12 << endl;

          if (mass12 < 1.7 || mass12 > 2.) continue;
          allPlots["d0_mass"+chTag+"_check"]->Fill(mass12,wgt);
          allPlots["d0_mass"+chTag+"_check_no_weight"]->Fill(mass12,norm);
          allPlots["d0_mass_all_check"]->Fill(mass12,wgt);
          allPlots["d0_mass_all_check_no_weight"]->Fill(mass12,norm);
          if (mass12 < 1.81 || mass12 > 1.92) continue;
          ++nD0;
          allPlots["d0_mass"+chTag+"_d0cand"]->Fill(mass12,wgt);
          allPlots["d0_mass"+chTag+"_d0cand_no_weight"]->Fill(mass12,norm);
          allPlots["d0_mass_all_d0cand"]->Fill(mass12,wgt);
          allPlots["d0_mass_all_d0cand_no_weight"]->Fill(mass12,norm);
          allPlots["d0_pt"+chTag+"_d0cand"]->Fill((p_track1+p_track2).Pt(),wgt);
          allPlots["d0_pt"+chTag+"_d0cand_no_weight"]->Fill((p_track1+p_track2).Pt(),norm);
          allPlots["d0_pt_all_d0cand"]->Fill((p_track1+p_track2).Pt(),wgt);
          allPlots["d0_pt_all_d0cand_no_weight"]->Fill((p_track1+p_track2).Pt(),norm);
          allPlots["d0_eta"+chTag+"_d0cand"]->Fill((p_track1+p_track2).Eta(),wgt);
          allPlots["d0_eta"+chTag+"_d0cand_no_weight"]->Fill((p_track1+p_track2).Eta(),norm);
          allPlots["d0_eta_all_d0cand"]->Fill((p_track1+p_track2).Eta(),wgt);
          allPlots["d0_eta_all_d0cand_no_weight"]->Fill((p_track1+p_track2).Eta(),norm);
          allPlots["d0_phi"+chTag+"_d0cand"]->Fill((p_track1+p_track2).Phi(),wgt);
          allPlots["d0_phi"+chTag+"_d0cand_no_weight"]->Fill((p_track1+p_track2).Phi(),norm);
          allPlots["d0_phi_all_d0cand"]->Fill((p_track1+p_track2).Phi(),wgt);
          allPlots["d0_phi_all_d0cand_no_weight"]->Fill((p_track1+p_track2).Phi(),norm);
          allPlots["d0_ptfrac"+chTag+"_d0cand"]->Fill((p_track1+p_track2).Pt()/allJetsVec[indJet[l]].getVec().Pt(),wgt);
          allPlots["d0_ptfrac"+chTag+"_d0cand_no_weight"]->Fill((p_track1+p_track2).Pt()/allJetsVec[indJet[l]].getVec().Pt(),norm);
          allPlots["d0_ptfrac_all_d0cand"]->Fill((p_track1+p_track2).Pt()/allJetsVec[indJet[l]].getVec().Pt(),wgt);
          allPlots["d0_ptfrac_all_d0cand_no_weight"]->Fill((p_track1+p_track2).Pt()/allJetsVec[indJet[l]].getVec().Pt(),norm);
          allPlots["kappa_d0_pt"+chTag+"_d0cand"]->Fill(p_track2.Pt(),wgt);
          allPlots["kappa_d0_pt"+chTag+"_d0cand_no_weight"]->Fill(p_track2.Pt(),norm);
          allPlots["kappa_d0_pt_all_d0cand"]->Fill(p_track2.Pt(),wgt);
          allPlots["kappa_d0_pt_all_d0cand_no_weight"]->Fill(p_track2.Pt(),norm);
          allPlots["kappa_d0_eta"+chTag+"_d0cand"]->Fill(p_track2.Eta(),wgt);
          allPlots["kappa_d0_eta"+chTag+"_d0cand_no_weight"]->Fill(p_track2.Eta(),norm);
          allPlots["kappa_d0_eta_all_d0cand"]->Fill(p_track2.Eta(),wgt);
          allPlots["kappa_d0_eta_all_d0cand_no_weight"]->Fill(p_track2.Eta(),norm);
          allPlots["kappa_d0_phi"+chTag+"_d0cand"]->Fill(p_track2.Phi(),wgt);
          allPlots["kappa_d0_phi"+chTag+"_d0cand_no_weight"]->Fill(p_track2.Phi(),norm);
          allPlots["kappa_d0_phi_all_d0cand"]->Fill(p_track2.Phi(),wgt);
          allPlots["kappa_d0_phi_all_d0cand_no_weight"]->Fill(p_track2.Phi(),norm);
          allPlots["pi_d0_pt"+chTag+"_d0cand"]->Fill(p_track1.Pt(),wgt);
          allPlots["pi_d0_pt"+chTag+"_d0cand_no_weight"]->Fill(p_track1.Pt(),norm);
          allPlots["pi_d0_pt_all_d0cand"]->Fill(p_track1.Pt(),wgt);
          allPlots["pi_d0_pt_all_d0cand_no_weight"]->Fill(p_track1.Pt(),norm);
          allPlots["pi_d0_eta"+chTag+"_d0cand"]->Fill(p_track1.Eta(),wgt);
          allPlots["pi_d0_eta"+chTag+"_d0cand_no_weight"]->Fill(p_track1.Eta(),norm);
          allPlots["pi_d0_eta_all_d0cand"]->Fill(p_track1.Eta(),wgt);
          allPlots["pi_d0_eta_all_d0cand_no_weight"]->Fill(p_track1.Eta(),norm);
          allPlots["pi_d0_phi"+chTag+"_d0cand"]->Fill(p_track1.Phi(),wgt);
          allPlots["pi_d0_phi"+chTag+"_d0cand_no_weight"]->Fill(p_track1.Phi(),norm);
          allPlots["pi_d0_phi_all_d0cand"]->Fill(p_track1.Phi(),wgt);
          allPlots["pi_d0_phi_all_d0cand_no_weight"]->Fill(p_track1.Phi(),norm);
          allPlots["kappa_pi_d0_dR"+chTag+"_d0cand"]->Fill(p_track1.DeltaR(p_track2),wgt);
          allPlots["kappa_pi_d0_dR"+chTag+"_d0cand_no_weight"]->Fill(p_track1.DeltaR(p_track2),norm);
          allPlots["kappa_pi_d0_dR_all_d0cand"]->Fill(p_track1.DeltaR(p_track2),wgt);
          allPlots["kappa_pi_d0_dR_all_d0cand_no_weight"]->Fill(p_track1.DeltaR(p_track2),norm);
          allPlots["j_d0_pt"+chTag+"_d0cand"]->Fill(allJetsVec[indJet[l]].getVec().Pt(),wgt);
          allPlots["j_d0_pt"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getVec().Pt(),norm);
          allPlots["j_d0_pt_all_d0cand"]->Fill(allJetsVec[indJet[l]].getVec().Pt(),wgt);
          allPlots["j_d0_pt_all_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getVec().Pt(),norm);
          allPlots["j_d0_eta"+chTag+"_d0cand"]->Fill(allJetsVec[indJet[l]].getVec().Eta(),wgt);
          allPlots["j_d0_eta"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getVec().Eta(),norm);
          allPlots["j_d0_eta_all_d0cand"]->Fill(allJetsVec[indJet[l]].getVec().Eta(),wgt);
          allPlots["j_d0_eta_all_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getVec().Eta(),norm);
          allPlots["j_d0_phi"+chTag+"_d0cand"]->Fill(allJetsVec[indJet[l]].getVec().Phi(),wgt);
          allPlots["j_d0_phi"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getVec().Phi(),norm);
          allPlots["j_d0_phi_all_d0cand"]->Fill(allJetsVec[indJet[l]].getVec().Phi(),wgt);
          allPlots["j_d0_phi_all_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getVec().Phi(),norm);
          allPlots["j_d0_csv"+chTag+"_d0cand"]->Fill(allJetsVec[indJet[l]].getCSV(),wgt);
          allPlots["j_d0_csv"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getCSV(),norm);
          allPlots["j_d0_csv_all_d0cand"]->Fill(allJetsVec[indJet[l]].getCSV(),wgt);
          allPlots["j_d0_csv_all_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getCSV(),norm);
          allPlots["j_d0_nch"+chTag+"_d0cand"]->Fill(allJetsVec[indJet[l]].getTracks().size(),wgt);
          allPlots["j_d0_nch"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getTracks().size(),norm);
          allPlots["j_d0_nch_all_d0cand"]->Fill(allJetsVec[indJet[l]].getTracks().size(),wgt);
          allPlots["j_d0_nch_all_d0cand_no_weight"]->Fill(allJetsVec[indJet[l]].getTracks().size(),norm);

          if (debug) cout << "starting d0cand plots" << endl;

          if (chTag == "_e") {
            allPlots["decay_all_d0cand"]->Fill(0.,wgt);
            allPlots["decay_all_d0cand_no_weight"]->Fill(0.,norm);
          }
          if (chTag == "_m") {
            allPlots["decay_all_d0cand"]->Fill(1.,wgt);
            allPlots["decay_all_d0cand_no_weight"]->Fill(1.,norm);
          }
          if (chTag == "_ee") {
            allPlots["decay_all_d0cand"]->Fill(2.,wgt);
            allPlots["decay_all_d0cand_no_weight"]->Fill(2.,norm);
          }
          if (chTag == "_em") {
            allPlots["decay_all_d0cand"]->Fill(3.,wgt);
            allPlots["decay_all_d0cand_no_weight"]->Fill(3.,norm);
          }
          if (chTag == "_mm") {
            allPlots["decay_all_d0cand"]->Fill(4.,wgt);
            allPlots["decay_all_d0cand_no_weight"]->Fill(4.,norm);
          }

          allPlots["lp_n"+chTag+"_d0cand"]->Fill(selLeptons.size(),wgt);
          allPlots["lp_n"+chTag+"_d0cand_no_weight"]->Fill(selLeptons.size(),norm);
          allPlots["lp_n_all_d0cand"]->Fill(selLeptons.size(),wgt);
          allPlots["lp_n_all_d0cand_no_weight"]->Fill(selLeptons.size(),norm);
          allPlots["lp1_pt"+chTag+"_d0cand"]->Fill(leptons[0].Pt(),wgt);
          allPlots["lp1_pt"+chTag+"_d0cand_no_weight"]->Fill(leptons[0].Pt(),norm);
          allPlots["lp1_pt_all_d0cand"]->Fill(leptons[0].Pt(),wgt);
          allPlots["lp1_pt_all_d0cand_no_weight"]->Fill(leptons[0].Pt(),norm);
          allPlots["lp1_eta"+chTag+"_d0cand"]->Fill(leptons[0].Eta(),wgt);
          allPlots["lp1_eta"+chTag+"_d0cand_no_weight"]->Fill(leptons[0].Eta(),norm);
          allPlots["lp1_eta_all_d0cand"]->Fill(leptons[0].Eta(),wgt);
          allPlots["lp1_eta_all_d0cand_no_weight"]->Fill(leptons[0].Eta(),norm);
          allPlots["lp1_phi"+chTag+"_d0cand"]->Fill(leptons[0].Phi(),wgt);
          allPlots["lp1_phi"+chTag+"_d0cand_no_weight"]->Fill(leptons[0].Phi(),norm);
          allPlots["lp1_phi_all_d0cand"]->Fill(leptons[0].Phi(),wgt);
          allPlots["lp1_phi_all_d0cand_no_weight"]->Fill(leptons[0].Phi(),norm);          
          if (selLeptons.size() == 2) {
            allPlots["lp2_pt"+chTag+"_d0cand"]->Fill(leptons[1].Pt(),wgt);
            allPlots["lp2_pt"+chTag+"_d0cand_no_weight"]->Fill(leptons[1].Pt(),norm);
            allPlots["lp2_pt_all_d0cand"]->Fill(leptons[1].Pt(),wgt);
            allPlots["lp2_pt_all_d0cand_no_weight"]->Fill(leptons[1].Pt(),norm);
            allPlots["lp2_eta"+chTag+"_d0cand"]->Fill(leptons[1].Eta(),wgt);
            allPlots["lp2_eta"+chTag+"_d0cand_no_weight"]->Fill(leptons[1].Eta(),norm);
            allPlots["lp2_eta_all_d0cand"]->Fill(leptons[1].Eta(),wgt);
            allPlots["lp2_eta_all_d0cand_no_weight"]->Fill(leptons[1].Eta(),norm);
            allPlots["lp2_phi"+chTag+"_d0cand"]->Fill(leptons[1].Phi(),wgt);
            allPlots["lp2_phi"+chTag+"_d0cand_no_weight"]->Fill(leptons[1].Phi(),norm);
            allPlots["lp2_phi_all_d0cand"]->Fill(leptons[1].Phi(),wgt);
            allPlots["lp2_phi_all_d0cand_no_weight"]->Fill(leptons[1].Phi(),norm);            
          }
          allPlots["met"+chTag+"_d0cand"]->Fill(met.Pt(),wgt);
          allPlots["met"+chTag+"_d0cand_no_weight"]->Fill(met.Pt(),norm);
          allPlots["met_all_d0cand"]->Fill(met.Pt(),wgt);
          allPlots["met_all_d0cand_no_weight"]->Fill(met.Pt(),norm);
          if (selLeptons.size() == 2) {
            allPlots["dilp_pt"+chTag+"_d0cand"]->Fill((leptons[0]+leptons[1]).Pt(),wgt);
            allPlots["dilp_pt"+chTag+"_d0cand_no_weight"]->Fill((leptons[0]+leptons[1]).Pt(),norm);
            allPlots["dilp_pt_all_d0cand"]->Fill((leptons[0]+leptons[1]).Pt(),wgt);
            allPlots["dilp_pt_all_d0cand_no_weight"]->Fill((leptons[0]+leptons[1]).Pt(),norm);
            allPlots["dilp_mass"+chTag+"_d0cand"]->Fill((leptons[0]+leptons[1]).M(),wgt);
            allPlots["dilp_mass"+chTag+"_d0cand_no_weight"]->Fill((leptons[0]+leptons[1]).M(),norm);
            allPlots["dilp_mass_all_d0cand"]->Fill((leptons[0]+leptons[1]).M(),wgt);
            allPlots["dilp_mass_all_d0cand_no_weight"]->Fill((leptons[0]+leptons[1]).M(),norm);
            allPlots["dilp_charge"+chTag+"_d0cand"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
            allPlots["dilp_charge"+chTag+"_d0cand_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
            allPlots["dilp_charge_all_d0cand"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
            allPlots["dilp_charge_all_d0cand_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
          }
          allPlots["j_n"+chTag+"_d0cand"]->Fill(allJetsVec.size(),wgt);
          allPlots["j_n"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec.size(),norm);
          allPlots["j_n_all_d0cand"]->Fill(allJetsVec.size(),wgt);
          allPlots["j_n_all_d0cand_no_weight"]->Fill(allJetsVec.size(),norm);
          for (size_t iJet = 0; iJet < allJetsVec.size(); iJet++) {
            allPlots["j_pt"+chTag+"_d0cand"]->Fill(allJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["j_pt"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[iJet].getVec().Pt(),norm);
            allPlots["j_pt_all_d0cand"]->Fill(allJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["j_pt_all_d0cand_no_weight"]->Fill(allJetsVec[iJet].getVec().Pt(),norm);
            allPlots["j_eta"+chTag+"_d0cand"]->Fill(allJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["j_eta"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[iJet].getVec().Eta(),norm);
            allPlots["j_eta_all_d0cand"]->Fill(allJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["j_eta_all_d0cand_no_weight"]->Fill(allJetsVec[iJet].getVec().Eta(),norm);
            allPlots["j_phi"+chTag+"_d0cand"]->Fill(allJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["j_phi"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[iJet].getVec().Phi(),norm);
            allPlots["j_phi_all_d0cand"]->Fill(allJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["j_phi_all_d0cand_no_weight"]->Fill(allJetsVec[iJet].getVec().Phi(),norm);
            allPlots["j_csv"+chTag+"_d0cand"]->Fill(allJetsVec[iJet].getCSV(),wgt);
            allPlots["j_csv"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[iJet].getCSV(),norm);
            allPlots["j_csv_all_d0cand"]->Fill(allJetsVec[iJet].getCSV(),wgt);
            allPlots["j_csv_all_d0cand_no_weight"]->Fill(allJetsVec[iJet].getCSV(),norm);
            allPlots["j_nch"+chTag+"_d0cand"]->Fill(allJetsVec[iJet].getTracks().size(),wgt);
            allPlots["j_nch"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[iJet].getTracks().size(),norm);
            allPlots["j_nch_all_d0cand"]->Fill(allJetsVec[iJet].getTracks().size(),wgt);
            allPlots["j_nch_all_d0cand_no_weight"]->Fill(allJetsVec[iJet].getTracks().size(),norm);
          }
          if (allJetsVec.size() > 0) {
            allPlots["j1_pt"+chTag+"_d0cand"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
            allPlots["j1_pt"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
            allPlots["j1_pt_all_d0cand"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
            allPlots["j1_pt_all_d0cand_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
            allPlots["j1_eta"+chTag+"_d0cand"]->Fill(allJetsVec[0].getVec().Eta(),wgt);
            allPlots["j1_eta"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[0].getVec().Eta(),norm);
            allPlots["j1_eta_all_d0cand"]->Fill(allJetsVec[0].getVec().Eta(),wgt);
            allPlots["j1_eta_all_d0cand_no_weight"]->Fill(allJetsVec[0].getVec().Eta(),norm);
            allPlots["j1_phi"+chTag+"_d0cand"]->Fill(allJetsVec[0].getVec().Phi(),wgt);
            allPlots["j1_phi"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[0].getVec().Phi(),norm);
            allPlots["j1_phi_all_d0cand"]->Fill(allJetsVec[0].getVec().Phi(),wgt);
            allPlots["j1_phi_all_d0cand_no_weight"]->Fill(allJetsVec[0].getVec().Phi(),norm);
            allPlots["j1_csv"+chTag+"_d0cand"]->Fill(allJetsVec[0].getCSV(),wgt);
            allPlots["j1_csv"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[0].getCSV(),norm);
            allPlots["j1_csv_all_d0cand"]->Fill(allJetsVec[0].getCSV(),wgt);
            allPlots["j1_csv_all_d0cand_no_weight"]->Fill(allJetsVec[0].getCSV(),norm);
            allPlots["j1_nch"+chTag+"_d0cand"]->Fill(allJetsVec[0].getTracks().size(),wgt);
            allPlots["j1_nch"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[0].getTracks().size(),norm);
            allPlots["j1_nch_all_d0cand"]->Fill(allJetsVec[0].getTracks().size(),wgt);
            allPlots["j1_nch_all_d0cand_no_weight"]->Fill(allJetsVec[0].getTracks().size(),norm);
            if (allJetsVec.size() > 1) {
              allPlots["j2_pt"+chTag+"_d0cand"]->Fill(allJetsVec[1].getVec().Pt(),wgt);
              allPlots["j2_pt"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[1].getVec().Pt(),norm);
              allPlots["j2_pt_all_d0cand"]->Fill(allJetsVec[1].getVec().Pt(),wgt);
              allPlots["j2_pt_all_d0cand_no_weight"]->Fill(allJetsVec[1].getVec().Pt(),norm);
              allPlots["j2_eta"+chTag+"_d0cand"]->Fill(allJetsVec[1].getVec().Eta(),wgt);
              allPlots["j2_eta"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[1].getVec().Eta(),norm);
              allPlots["j2_eta_all_d0cand"]->Fill(allJetsVec[1].getVec().Eta(),wgt);
              allPlots["j2_eta_all_d0cand_no_weight"]->Fill(allJetsVec[1].getVec().Eta(),norm);
              allPlots["j2_phi"+chTag+"_d0cand"]->Fill(allJetsVec[1].getVec().Phi(),wgt);
              allPlots["j2_phi"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[1].getVec().Phi(),norm);
              allPlots["j2_phi_all_d0cand"]->Fill(allJetsVec[1].getVec().Phi(),wgt);
              allPlots["j2_phi_all_d0cand_no_weight"]->Fill(allJetsVec[1].getVec().Phi(),norm);
              allPlots["j2_csv"+chTag+"_d0cand"]->Fill(allJetsVec[1].getCSV(),wgt);
              allPlots["j2_csv"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[1].getCSV(),norm);
              allPlots["j2_csv_all_d0cand"]->Fill(allJetsVec[1].getCSV(),wgt);
              allPlots["j2_csv_all_d0cand_no_weight"]->Fill(allJetsVec[1].getCSV(),norm);
              allPlots["j2_nch"+chTag+"_d0cand"]->Fill(allJetsVec[1].getTracks().size(),wgt);
              allPlots["j2_nch"+chTag+"_d0cand_no_weight"]->Fill(allJetsVec[1].getTracks().size(),norm);
              allPlots["j2_nch_all_d0cand"]->Fill(allJetsVec[1].getTracks().size(),wgt);
              allPlots["j2_nch_all_d0cand_no_weight"]->Fill(allJetsVec[1].getTracks().size(),norm);
            }
          }
          allPlots["lj_n"+chTag+"_d0cand"]->Fill(lightJetsVec.size(),wgt);
          allPlots["lj_n"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec.size(),norm);
          allPlots["lj_n_all_d0cand"]->Fill(lightJetsVec.size(),wgt);
          allPlots["lj_n_all_d0cand_no_weight"]->Fill(lightJetsVec.size(),norm);
          for (size_t iJet = 0; iJet < lightJetsVec.size(); iJet++) {
            allPlots["lj_pt"+chTag+"_d0cand"]->Fill(lightJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["lj_pt"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Pt(),norm);
            allPlots["lj_pt_all_d0cand"]->Fill(lightJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["lj_pt_all_d0cand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Pt(),norm);
            allPlots["lj_eta"+chTag+"_d0cand"]->Fill(lightJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["lj_eta"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Eta(),norm);
            allPlots["lj_eta_all_d0cand"]->Fill(lightJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["lj_eta_all_d0cand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Eta(),norm);
            allPlots["lj_phi"+chTag+"_d0cand"]->Fill(lightJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["lj_phi"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Phi(),norm);
            allPlots["lj_phi_all_d0cand"]->Fill(lightJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["lj_phi_all_d0cand_no_weight"]->Fill(lightJetsVec[iJet].getVec().Phi(),norm);
            allPlots["lj_nch"+chTag+"_d0cand"]->Fill(lightJetsVec[iJet].getTracks().size(),wgt);
            allPlots["lj_nch"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec[iJet].getTracks().size(),norm);
            allPlots["lj_nch_all_d0cand"]->Fill(lightJetsVec[iJet].getTracks().size(),wgt);
            allPlots["lj_nch_all_d0cand_no_weight"]->Fill(lightJetsVec[iJet].getTracks().size(),norm);
          }
          if (lightJetsVec.size() > 0) {
            allPlots["lj1_pt"+chTag+"_d0cand"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
            allPlots["lj1_pt"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
            allPlots["lj1_pt_all_d0cand"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
            allPlots["lj1_pt_all_d0cand_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
            allPlots["lj1_eta"+chTag+"_d0cand"]->Fill(lightJetsVec[0].getVec().Eta(),wgt);
            allPlots["lj1_eta"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec[0].getVec().Eta(),norm);
            allPlots["lj1_eta_all_d0cand"]->Fill(lightJetsVec[0].getVec().Eta(),wgt);
            allPlots["lj1_eta_all_d0cand_no_weight"]->Fill(lightJetsVec[0].getVec().Eta(),norm);
            allPlots["lj1_phi"+chTag+"_d0cand"]->Fill(lightJetsVec[0].getVec().Phi(),wgt);
            allPlots["lj1_phi"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec[0].getVec().Phi(),norm);
            allPlots["lj1_phi_all_d0cand"]->Fill(lightJetsVec[0].getVec().Phi(),wgt);
            allPlots["lj1_phi_all_d0cand_no_weight"]->Fill(lightJetsVec[0].getVec().Phi(),norm);
            allPlots["lj1_nch"+chTag+"_d0cand"]->Fill(lightJetsVec[0].getTracks().size(),wgt);
            allPlots["lj1_nch"+chTag+"_d0cand_no_weight"]->Fill(lightJetsVec[0].getTracks().size(),norm);
            allPlots["lj1_nch_all_d0cand"]->Fill(lightJetsVec[0].getTracks().size(),wgt);
            allPlots["lj1_nch_all_d0cand_no_weight"]->Fill(lightJetsVec[0].getTracks().size(),norm);
          }
          allPlots["bj_n"+chTag+"_d0cand"]->Fill(bJetsVec.size(),wgt);
          allPlots["bj_n"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec.size(),norm);
          allPlots["bj_n_all_d0cand"]->Fill(bJetsVec.size(),wgt);
          allPlots["bj_n_all_d0cand_no_weight"]->Fill(bJetsVec.size(),norm);
          for (size_t iJet = 0; iJet < bJetsVec.size(); iJet++) {
            allPlots["bj_pt"+chTag+"_d0cand"]->Fill(bJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["bj_pt"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[iJet].getVec().Pt(),norm);
            allPlots["bj_pt_all_d0cand"]->Fill(bJetsVec[iJet].getVec().Pt(),wgt);
            allPlots["bj_pt_all_d0cand_no_weight"]->Fill(bJetsVec[iJet].getVec().Pt(),norm);
            allPlots["bj_eta"+chTag+"_d0cand"]->Fill(bJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["bj_eta"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[iJet].getVec().Eta(),norm);
            allPlots["bj_eta_all_d0cand"]->Fill(bJetsVec[iJet].getVec().Eta(),wgt);
            allPlots["bj_eta_all_d0cand_no_weight"]->Fill(bJetsVec[iJet].getVec().Eta(),norm);
            allPlots["bj_phi"+chTag+"_d0cand"]->Fill(bJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["bj_phi"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[iJet].getVec().Phi(),norm);
            allPlots["bj_phi_all_d0cand"]->Fill(bJetsVec[iJet].getVec().Phi(),wgt);
            allPlots["bj_phi_all_d0cand_no_weight"]->Fill(bJetsVec[iJet].getVec().Phi(),norm);
            allPlots["bj_csv"+chTag+"_d0cand"]->Fill(bJetsVec[iJet].getCSV(),wgt);
            allPlots["bj_csv"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[iJet].getCSV(),norm);
            allPlots["bj_csv_all_d0cand"]->Fill(bJetsVec[iJet].getCSV(),wgt);
            allPlots["bj_csv_all_d0cand_no_weight"]->Fill(bJetsVec[iJet].getCSV(),norm);
            allPlots["bj_nch"+chTag+"_d0cand"]->Fill(bJetsVec[iJet].getTracks().size(),wgt);
            allPlots["bj_nch"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[iJet].getTracks().size(),norm);
            allPlots["bj_nch_all_d0cand"]->Fill(bJetsVec[iJet].getTracks().size(),wgt);
            allPlots["bj_nch_all_d0cand_no_weight"]->Fill(bJetsVec[iJet].getTracks().size(),norm);
          }
          if (bJetsVec.size() > 0) {
            allPlots["bj1_pt"+chTag+"_d0cand"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
            allPlots["bj1_pt"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
            allPlots["bj1_pt_all_d0cand"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
            allPlots["bj1_pt_all_d0cand_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
            allPlots["bj1_eta"+chTag+"_d0cand"]->Fill(bJetsVec[0].getVec().Eta(),wgt);
            allPlots["bj1_eta"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[0].getVec().Eta(),norm);
            allPlots["bj1_eta_all_d0cand"]->Fill(bJetsVec[0].getVec().Eta(),wgt);
            allPlots["bj1_eta_all_d0cand_no_weight"]->Fill(bJetsVec[0].getVec().Eta(),norm);
            allPlots["bj1_phi"+chTag+"_d0cand"]->Fill(bJetsVec[0].getVec().Phi(),wgt);
            allPlots["bj1_phi"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[0].getVec().Phi(),norm);
            allPlots["bj1_phi_all_d0cand"]->Fill(bJetsVec[0].getVec().Phi(),wgt);
            allPlots["bj1_phi_all_d0cand_no_weight"]->Fill(bJetsVec[0].getVec().Phi(),norm);
            allPlots["bj1_csv"+chTag+"_d0cand"]->Fill(bJetsVec[0].getCSV(),wgt);
            allPlots["bj1_csv"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[0].getCSV(),norm);
            allPlots["bj1_csv_all_d0cand"]->Fill(bJetsVec[0].getCSV(),wgt);
            allPlots["bj1_csv_all_d0cand_no_weight"]->Fill(bJetsVec[0].getCSV(),norm);
            allPlots["bj1_nch"+chTag+"_d0cand"]->Fill(bJetsVec[0].getTracks().size(),wgt);
            allPlots["bj1_nch"+chTag+"_d0cand_no_weight"]->Fill(bJetsVec[0].getTracks().size(),norm);
            allPlots["bj1_nch_all_d0cand"]->Fill(bJetsVec[0].getTracks().size(),wgt);
            allPlots["bj1_nch_all_d0cand_no_weight"]->Fill(bJetsVec[0].getTracks().size(),norm);
          }

          //looking for lepton
          if (debug) cout << "third lepton" << endl;
          for (size_t k = 0; k < tracks.size(); k++) {
            if (k == i) continue;
            if (k == j) continue;
            if (debug) cout << "third lepton possible" << endl;

            if (abs(tracks[k].second) != 13 && abs(tracks[k].second) != 11) continue;
            if (debug) cout << "third lepton found" << endl;

            if (tracks[j].second/abs(tracks[j].second) == -tracks[k].second/abs(tracks[k].second)) {
              //Kaon and lepton have same charge
              //correct mass assumption
              if (debug) cout << "correct mass assumption" << endl;
              allPlots["d0_mass_lp"+chTag+"_d0cand"]->Fill(mass12,wgt);
              allPlots["d0_mass_lp"+chTag+"_d0cand_no_weight"]->Fill(mass12,norm);
              allPlots["d0_mass_lp_all_d0cand"]->Fill(mass12,wgt);
              allPlots["d0_mass_lp_all_d0cand_no_weight"]->Fill(mass12,norm);
              if (abs(tracks[k].second) == 13) {
                allPlots["d0_mass_mu"+chTag+"_d0cand"]->Fill(mass12,wgt);
                allPlots["d0_mass_mu"+chTag+"_d0cand_no_weight"]->Fill(mass12,norm);
                allPlots["d0_mass_mu_all_d0cand"]->Fill(mass12,wgt);
                allPlots["d0_mass_mu_all_d0cand_no_weight"]->Fill(mass12,norm);
              }
              if (abs(tracks[k].second) == 11) {
                allPlots["d0_mass_el"+chTag+"_d0cand"]->Fill(mass12,wgt);
                allPlots["d0_mass_el"+chTag+"_d0cand_no_weight"]->Fill(mass12,norm);
                allPlots["d0_mass_el_all_d0cand"]->Fill(mass12,wgt);
                allPlots["d0_mass_el_all_d0cand_no_weight"]->Fill(mass12,norm);
              }
            }
          }
          //looking for pion
          if (debug) cout << "D*->pi+D0" << endl;
          for (size_t k = 0; k < tracks.size(); k++) {
            if (k == i) continue;
            if (k == j) continue;

            if (abs(tracks[k].second) != 211) continue;
            if (debug) cout << "Pion found" << endl;

            TLorentzVector p_track3, p_cand;
            p_track3.SetPtEtaPhiM(tracks[k].first.Pt(), tracks[k].first.Eta(), tracks[k].first.Phi(), gMassPi);
            if (debug) cout << k << ": " << tracks[k].first.Pt() << " " << tracks[k].first.Eta() << " " << tracks[k].first.Phi() << " " << gMassPi << endl;
            allPlots["pi_ds_pt"+chTag+"_d0cand"]->Fill(p_track3.Pt(),wgt);
            allPlots["pi_ds_pt"+chTag+"_d0cand_no_weight"]->Fill(p_track3.Pt(),norm);
            allPlots["pi_ds_pt_all_d0cand"]->Fill(p_track3.Pt(),wgt);
            allPlots["pi_ds_pt_all_d0cand_no_weight"]->Fill(p_track3.Pt(),norm);
            if (tracks[j].second/abs(tracks[j].second) == -tracks[k].second/abs(tracks[k].second)) {
              // Kaon and pion have opposite charges
              // I.e. correct mass assumption
              if (debug) cout << "correct mass assumption" << endl;

              p_cand = p_track1+p_track2+p_track3;
              allPlots["ds_mass"+chTag+"_d0cand"]->Fill(p_cand.M(), wgt);
              allPlots["ds_mass"+chTag+"_d0cand_no_weight"]->Fill(p_cand.M(),norm);
              allPlots["ds_mass_all_d0cand"]->Fill(p_cand.M(), wgt);
              allPlots["ds_mass_all_d0cand_no_weight"]->Fill(p_cand.M(),norm);

              TLorentzVector p_jet = allJetsVec[indJet[l]].getVec();

              float deltam = p_cand.M() - mass12;

              allPlots["ds_d0_dmass"+chTag+"_d0cand"]->Fill(deltam, wgt);
              allPlots["ds_d0_dmass"+chTag+"_d0cand_no_weight"]->Fill(deltam, norm);
              allPlots["ds_d0_dmass_all_d0cand"]->Fill(deltam, wgt);
              allPlots["ds_d0_dmass_all_d0cand_no_weight"]->Fill(deltam, norm);
            }
          }
          if (debug) cout << "d0cand plots DONE" << endl;          
        }
      }
    }
    if (nD0 > 0) {
      allPlots["d0_n"+chTag+"_d0cand"]->Fill(nD0,wgt);
      allPlots["d0_n"+chTag+"_d0cand_no_weight"]->Fill(nD0,norm);
      allPlots["d0_n_all_d0cand"]->Fill(nD0,wgt);
      allPlots["d0_n_all_d0cand_no_weight"]->Fill(nD0,norm);
    }
    if (debug) cout << "D0 and D* DONE" << endl;

  }

  //close input file
  f->Close();

  //save histos to file  
  TString selPrefix("");  
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  fOut->cd();
  if (debug) cout << "writing histograms" << endl;

  for (auto& it : allPlots)  { 
    if (debug) cout << it.second->GetName() << endl;
    if (debug) cout << it.second->GetEntries() << endl;

    //fOut->cd( dir );
    it.second->SetDirectory(fOut); it.second->Write(); 
    fOut->cd();
  }
  if (debug) cout << "writing histograms DONE" << endl;
  if (debug) cout << "closing ROOT file" << endl;
  fOut->Close();
}

