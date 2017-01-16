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

#include "TopLJets2015/TopAnalysis/interface/rochcor2016.h"

#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

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
  debug = true;
  if(debug) cout << "in RunTop" << endl;

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
  if(ev.isData) runSysts=false;
  bool requireEletriggerOnly(false);
  if(ev.isData && filename.Contains("SingleElectron")) requireEletriggerOnly=true;
  bool requireMutriggerOnly(false);
  if(ev.isData && filename.Contains("SingleMuon"))     requireMutriggerOnly=true;
  bool requireEETriggers(false);
  if(ev.isData && filename.Contains("DoubleEG"))       requireEETriggers=true;
  bool requireMMTriggers(false);
  if(ev.isData && filename.Contains("DoubleMuon"))     requireMMTriggers=true;
  bool requireEMTriggers(false);
  if(ev.isData && filename.Contains("MuonEG"))         requireEMTriggers=true;

  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;

  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
  {
    if(debug) cout << "loading pileup weight" << endl;
    TString puWgtUrl(era+"/pileupWgts.root");
    gSystem->ExpandPathName(puWgtUrl);
    TFile *fIn=TFile::Open(puWgtUrl);
    for(size_t i=0; i<3; i++)
    {
      TString grName("pu_nom");
      if(i==1) grName="pu_down";
      if(i==2) grName="pu_up";
      TGraph *puData=(TGraph *)fIn->Get(grName);
      Float_t totalData=puData->Integral();
      TH1 *tmp=(TH1 *)puTrue->Clone("tmp");
      for(Int_t xbin=1; xbin<=tmp->GetXaxis()->GetNbins(); xbin++)
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
  if(debug) cout << "loading pileup weight DONE" << endl;

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
  if(!ev.isData)
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
    if(isTTbar)
    {
      if(filename.Contains("_herwig")) btagExpPostFix="_herwig";
      if(filename.Contains("_scaleup")) btagExpPostFix="_scaleup";
      if(filename.Contains("_scaledown")) btagExpPostFix="_scaledown";
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

  //LIST OF SYSTEMATICS

  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);
  std::vector<TString> lfsVec = { "_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  std::vector<TString> cutVec = { "", "_lep", "_jpsi", "_csv", "_meson" };
  std::vector<TString> wgtVec = { "", "_no_weight" };

  for(int i = 0; i < (int)lfsVec.size(); i++) {
    TString tag(lfsVec[i]);
    for(int k = 0; k < (int)wgtVec.size(); k++) {
      TString weight(wgtVec[k]);
      allPlots["lp_pt"+tag+weight] = new TH1F("lp_pt"+tag+weight,";Leading lepton p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
      allPlots["lp2_pt"+tag+weight] = new TH1F("lp2_pt"+tag+weight,";Subleading lepton p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
      allPlots["met"+tag+weight] = new TH1F("met"+tag+weight,";MET (GeV);Events / (20 GeV)", 10, 0., 200.);
      allPlots["dilp_pt"+tag+weight] = new TH1F("dilp_pt"+tag+weight,";Dilepton pair p_{T} (GeV);Events / (10 GeV)", 20, 0., 200.);
      allPlots["dilp_mass"+tag+weight] = new TH1F("dilp_mass"+tag+weight,";Dilepton pair mass (GeV);Events / (10 GeV)", 20, 0., 200.);
      allPlots["dilp_charge"+tag+weight] = new TH1F("dilp_charge"+tag+weight,";Dilepton pair charge;Events", 5, -2., 2.);
      allPlots["Z_mass"+tag+weight] = new TH1F("Z_mass"+tag+weight,";Z boson candidate mass (GeV);Events / (1 GeV)",30, 81., 111.);
      allPlots["pf_id"+tag+weight] = new TH1F("pf_id"+tag+weight,";Particle PDG ID;Events", 440, -220., 220.);
      allPlots["lp_jets_dR"+tag+weight] = new TH1F("lp_jets_dR"+tag+weight,";#deltR(l, jets);Events / 0.05", 20, 0., 1.);
      for(int j = 0; j < (int)cutVec.size(); j++) {
        TString cut(cutVec[j]);
        allPlots["j_pt"+tag+cut+weight] = new TH1F("j_pt"+tag+cut+weight,";Leading jet p_{T} (GeV);Events / (20 GeV)", 15, 0., 300.); 
        allPlots["lj_pt"+tag+cut+weight] = new TH1F("lj_pt"+tag+cut+weight,";Leading light jet p_{T} (GeV);Events / (20 GeV)", 15, 0., 300.);
        //FIXME : only "", lep and csv
        allPlots["bj_pt"+tag+cut+weight] = new TH1F("bj_pt"+tag+cut+weight,";Leading b jet p_{T} (GeV);Events / (20 GeV)", 15, 0., 300.);
        allPlots["j_csv"+tag+cut+weight] = new TH1F("j_csv"+tag+cut+weight,";Jet CSV discriminant;Events / 0.1", 10, 0., 1.);
        //FIXME : only "", lep, csv, and jpsi
        allPlots["lp_n"+tag+cut+weight] = new TH1F("lp_n"+tag+cut+weight,";Lepton multiplicity;Events", 3, 0., 3.);
        allPlots["j_n"+tag+cut+weight] = new TH1F("j_n"+tag+cut+weight,";Jet multiplicity (p_{T} > 30 GeV);Events", 10, 0., 10.);
        allPlots["lj_n"+tag+cut+weight] = new TH1F("lj_n"+tag+cut+weight,";Light jet multiplicity (p_{T} > 30 GeV);Events" ,10,0,10.);
        allPlots["bj_n"+tag+cut+weight] = new TH1F("bj_n"+tag+cut+weight,";b-tagged jet multiplicity (CSV > 0.8);Events", 4, 1., 5.);
        allPlots["pf_n"+tag+cut+weight] = new TH1F("pf_n"+tag+cut+weight,";Particle multiplicity;Events / 10", 5, 0., 5.);
        //FIXME : only "", lep
        allPlots["jpsi_mass_large"+tag+cut+weight] = new TH1F("jpsi_mass"+tag+cut+weight,";M_{#mu^{+}#mu^{-}} (GeV);Events / (25 MeV)", 40, 2.6, 3.6);
        allPlots["jpsi_mass"+tag+cut+weight] = new TH1F("jpsi_mass"+tag+cut+weight,";M_{#mu^{+}#mu^{-}} (GeV);Events / (20 MeV)", 20, 2.9, 3.3);
        allPlots["jpsi_kaon_mass"+tag+cut+weight] = new TH1F("jpsi_kaon_mass"+tag+cut+weight,";M_{#mu^{+}#mu^{-}#kappa^{#pm}} (GeV);Events / (30 MeV)", 50, 4.5, 6.);
        allPlots["pf_dxy"+tag+cut+weight] = new TH1F("pf_dxy"+tag+cut+weight,";d_{xy} (cm);Events / (20 #mum)", 100, 0., 0.1);
        allPlots["pf_dz"+tag+cut+weight] = new TH1F("pf_dz"+tag+cut+weight,";d_{z} (cm);Events / (20 #mum)", 100, 0., 0.1);
        allPlots["pf_dxyE"+tag+cut+weight] = new TH1F("pf_dxyE"+tag+cut+weight,";#sigma(d_{xy}) (cm);Events / (20 #mum)", 100, 0., 0.1);
        allPlots["pf_dzE"+tag+cut+weight] = new TH1F("pf_dzE"+tag+cut+weight,";#sigma(d_{z}) (cm);Events / (20 #mum)", 100, 0., 0.1);
        allPlots["pf_dxy_sig"+tag+cut+weight] = new TH1F("pf_dxy_significance"+tag+cut+weight,";d_{xy}/#sigmad_{xy};Events / 0.3", 100, 0., 30.);
        allPlots["pf_dz_sig"+tag+cut+weight] = new TH1F("pf_dz_significance"+tag+cut+weight,";d_{z}/#sigma_{z};Events / 0.3", 100, 0., 30.);
        //FIXME : only jpsi
        allPlots["d0_mass_large"+tag+cut+weight] = new TH1F("d0_mass_large"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (15 MeV)", 20, 1.7, 2.0);
        allPlots["d0_mass"+tag+cut+weight] = new TH1F("d0_mass"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (10 MeV)", 11, 1.81, 1.92);
        allPlots["d0_mass_lp"+tag+cut+weight] = new TH1F("d0_mass_lp"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (10 MeV)", 11, 1.81, 1.92);
        allPlots["d0_mass_mu"+tag+cut+weight] = new TH1F("d0_mass_mu"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (10 MeV)", 11, 1.81, 1.92);
        allPlots["d0_mass_el"+tag+cut+weight] = new TH1F("d0_mass_el"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (10 MeV)", 11, 1.81, 1.92);
        allPlots["ds_mass"+tag+cut+weight] = new TH1F("ds_mass"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}#plus#pi^{#mp}};Events / 10 MeV" ,200,0.,2.0);
        allPlots["pi_pt"+tag+cut+weight] = new TH1F("pi_pt"+tag+cut+weight,";#pi^{#pm} p_{T} (GeV);Events / (5 GeV)", 10, 0., 50.);
        allPlots["ds_d0_dmass_loose"+tag+cut+weight] = new TH1F("ds_d0_dmass_loose"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}#plus#pi^{#mp}} #minus M_{#kappa^{#pm}#pi^{#mp}} (GeV);Events / (1 MeV)", 30, 0.14, 0.17);
        allPlots["ds_d0_dmass"+tag+cut+weight] = new TH1F("ds_d0_dmass"+tag+cut+weight,";M_{#kappa^{#pm}#pi^{#mp}#plus#pi^{#mp}} #minus M_{#kapp^{pm}#pi^{#mp}};Events / (1 MeV)", 30, 0.14, 0.17);
        //FIXME : only meson
        allPlots["nevt"+tag+cut+weight] = new TH1F("nevt"+tag+cut+weight,";;Events", 1, 1., 2.);
      }
    }
  }
  allPlots["mu_relIso"] = new TH1F("mu_relIso",";relIso;Events / 0.05", 20,0,1.);
  allPlots["el_relIso"] = new TH1F("el_relIso",";relIso;Events / 0.05", 20,0,1.);


  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  //for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }
  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
  {
    t->GetEntry(iev);
    if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
    allPlots["nevt_all"]->Fill(1,1);

    std::vector<int> tightLeptons,vetoLeptons;
    for(int il=0; il<ev.nl; il++) {
      float relIso = ev.l_relIso[il];
      bool passTightKin((ev.l_pt[il]>20 && fabs(ev.l_eta[il])<2.4 && abs(ev.l_id[il])==13 && relIso<0.25) //muons
          || (ev.l_pt[il]>30 && ((fabs(ev.l_eta[il])<1.479 && relIso<0.0893)  //electrons small eta
              || (fabs(ev.l_eta[il])>1.479 && fabs(ev.l_eta[il])<2.5 && relIso<0.121)) //electrons medium eta
            && abs(ev.l_id[il])==11)); // TOP mu cut for dilep

      bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);

      //Check veto here, but ONLY for lep+jets
      bool passVetoIso(  ev.l_id[il]==13 ? relIso<0.25 : true); 
      bool passVetoKin(  ev.l_pt[il]>10. && fabs(ev.l_eta[il])<2.5); // TOP veto

      if( ev.l_id[il] == 13 )
        allPlots["mu_relIso"]->Fill(relIso,1);
      else if( ev.l_id[il] == 11)
        allPlots["el_relIso"]->Fill(relIso,1);

      if(passTightKin && passTightId) {
        tightLeptons.push_back(il);
      }
      else if(passVetoKin && passVetoIso) vetoLeptons.push_back(il); 
    }
    if(debug) cout << "lepton selection DONE" << endl;

    //check if triggers have fired
    bool hasEETrigger(((ev.elTrigger>>3)&0x1)!=0 || ((ev.elTrigger>>2)&0x1)!=0);//HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
    bool hasMMTrigger(((ev.muTrigger>>4)&0x3)!=0);                              //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v && HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v
    bool hasEMTrigger(((ev.elTrigger>>4)&0x3)!=0);                              //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
    bool hasMuTrigger(((ev.muTrigger>>2)&0x3)!=0);                              //HLT_IsoMu22_v && HLT_IsoTkMu22_v
    bool hasEleTrigger((ev.elTrigger & 0x1)!=0);                                //HLT_Ele27_WPTight_Gsf_v
    if(!ev.isData) {	 
      hasMuTrigger=true;
      hasEleTrigger=true;
      hasEETrigger=true;
      hasMMTrigger=true;
      hasEMTrigger=true;
    } else {
      if(requireMutriggerOnly && !hasMuTrigger) continue;
      if(requireEletriggerOnly && !hasEleTrigger) continue;
      if(requireEETriggers && !hasEETrigger) continue;
      if(requireMMTriggers && !hasMMTrigger) continue;
      if(requireEMTriggers && !hasEMTrigger) continue;
    }

    //decide the channel
    if(debug) cout << "decide channel" << endl;
    TString chTag("");
    std::vector<int> selLeptons;
    bool passTighterKin(false),passIso(false);
    if(tightLeptons.size()==1 ) {
      //** Tighter cuts for lepton + jets **
      if(ev.l_id[tightLeptons[0]]==13) { // muon + jets
        passTighterKin = (ev.l_pt[tightLeptons[0]] > 26. && fabs(ev.l_eta[tightLeptons[0]])<2.1); // TOP mu cut for dilep
        passIso = (ev.l_relIso[tightLeptons[0]] < 0.15); //TOP mu cut for lep+jets
      } else if(ev.l_id[tightLeptons[0]]==11) { // electron + jets
        passTighterKin = (ev.l_pt[tightLeptons[0]] > 30.); //from TOP-15-005
        passIso = (ev.l_relIso[tightLeptons[0]] < 0.15); //TOP mu cut for lep+jets
      }

      if (passTighterKin && passIso) {
        selLeptons.push_back( tightLeptons[0] );
        if(debug) cout << "found 1 tight lepton" << endl;
      }
      //USE VETO HERE
      //no extra isolated leptons
      if(vetoLeptons.size()>0) continue; //veto only on lep+jets
    }
    if(tightLeptons.size()>=2) {
      selLeptons.push_back(tightLeptons[0]);
      selLeptons.push_back(tightLeptons[1]);
      if(debug) cout << "found 2 tight leptons" << endl;
    }
    if(debug) if(selLeptons.size()==0) cout << "NO LEPTONS!!" << endl;
    if(selLeptons.size()==0) continue;
    if(selLeptons.size()==1) {
      if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="e";
      if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger)  chTag="m";
    }
    if(selLeptons.size()==2) {
      if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*11 && hasEETrigger) chTag="ee";
      if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==13*13 && hasMMTrigger) chTag="mm";
      if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*13) {
        if(selLeptons.size()>=2 && (hasEMTrigger)) chTag="em";
        if(selLeptons.size()==1) {
          if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="em";
          if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger) chTag="em";
        }
      }
      if(hasMuTrigger && requireEletriggerOnly) chTag="";
    }
    if(chTag=="") continue;
    chTag = "_"+chTag;
    if(debug) cout << "decide channel DONE" << endl;

    //one good lepton either isolated or in the non-isolated sideband or a Z candidate
    Int_t lepIdx=-1;
    Bool_t isZ(false);//,isZPassingSIP3d(false);
    TLorentzVector l1p4,l2p4,dilp4;
    if(selLeptons.size()==1) 
      lepIdx=selLeptons[0];
    else if(selLeptons.size()==2) {	  
      if(debug) cout << "di-lepton" << endl;
      l1p4.SetPtEtaPhiM(ev.l_pt[selLeptons[0]],ev.l_eta[selLeptons[0]],ev.l_phi[selLeptons[0]],ev.l_mass[selLeptons[0]]);
      l2p4.SetPtEtaPhiM(ev.l_pt[selLeptons[1]],ev.l_eta[selLeptons[1]],ev.l_phi[selLeptons[1]],ev.l_mass[selLeptons[1]]);
      dilp4=l1p4+l2p4;
      if(ev.l_id[selLeptons[0]]==ev.l_id[selLeptons[1]]          && 
          ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]]<0 && 
          fabs(dilp4.M()-91.) < 15.) { 
        isZ=true; 
      }
      lepIdx=selLeptons[0];
      if(debug) cout << "di-lepton DONE" << endl;
    }

    if(lepIdx<0) continue;

    //lepton kinematics
    if(debug) cout << "checking lepton kinematics" << endl;
    TLorentzVector lp4;
    lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);

    if(debug) cout << "checking lepton kinematics DONE" << endl;

    //save lepton kinematics
    std::vector<TLorentzVector> leptons;
    for(size_t il=0; il<selLeptons.size(); il++) {
      int lepIdx=selLeptons[il];
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);
      leptons.push_back(lp4);
    }

    //select jets
    TLorentzVector jetDiff(0,0,0,0);
    int nbjets(0),ncjets(0),nljets(0);//,leadingJetIdx(-wgt);
    std::vector<int> resolvedJetIdx;
    std::vector<TLorentzVector> resolvedJetP4;
    std::vector<Jet> bJetsVec, lightJetsVec, allJetsVec;
    for (int k=0; k<ev.nj;k++) {
      //check kinematics
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

      //cross clean with respect to leptons 
      bool overlapsWithLepton(false);
      for(size_t il=0; il<leptons.size(); il++) {
        if(jp4.DeltaR(leptons[il])>0.4) continue;
        overlapsWithLepton=true;
      }
      if(overlapsWithLepton) continue;
      if(debug) cout << "Overlap with lepton DONE" << endl;

      //smear jet energy resolution for MC
      //jetDiff -= jp4;
      float genJet_pt(0);
      if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
      if(!ev.isData && genJet_pt>0) 
      {
        float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
        jp4 *= jerSmear;
      }
      //jetDiff += jp4;
      resolvedJetIdx.push_back(k);
      resolvedJetP4.push_back(jp4);

      // re-inforce kinematics cuts
      if(jp4.Pt()<30) continue;
      if(fabs(jp4.Eta()) > 2.4) continue;

      //b-tag
      if(debug) cout << "Starting b-tagging" << endl;
      float csv = ev.j_csv[k];	  
      bool isBTagged(csv>0.800);
      if(!ev.isData) {
        float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
        float expEff(1.0), jetBtagSF(1.0);
        if(abs(ev.j_hadflav[k])==4) { 
          ncjets++;
          expEff    = expBtagEff["c"]->Eval(jptForBtag); 
          jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
          jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
        } else if(abs(ev.j_hadflav[k])==5) { 
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
      if(debug) cout << "b-tagging DONE" << endl;

      //save jet
      Jet tmpj(jp4, csv, k);
      for(int ipf = 0; ipf < ev.npf; ipf++) {
        if(ev.pf_j[ipf] != k) continue;
        if(ev.pf_c[ipf]==0) continue;
        TLorentzVector tkP4(0,0,0,0);
        tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
        pfTrack pftk(tkP4, ev.pf_dxy[ipf], ev.pf_dxyE[ipf], ev.pf_dz[ipf], ev.pf_dzE[ipf], ev.pf_id[ipf]);
        tmpj.addTrack(pftk,ev.pf_id[ipf]);
      }
      tmpj.sortTracksByPt();

      if(isBTagged) bJetsVec.push_back(tmpj);
      else lightJetsVec.push_back(tmpj);
      allJetsVec.push_back(tmpj);
    }


    //event weight
    float wgt(1.0);
    float norm(1.0);
    std::vector<float> puWgts(3,1.0),topPtWgts(2,1.0);
    EffCorrection_t lepSelCorrWgt(1.0,0.0), triggerCorrWgt(1.0,0.0);
    if(debug) cout << "Lepton scale factors" << endl;
    if(!ev.isData) {
      //update lepton selection scale factors, if found
      //float lepTriggerSF(1.0),lepSelSF(1.0);
      //FIXME

      //account for pu weights and effect on normalization
      allPlots["puwgtctr"]->Fill(0.,1.0);
      if(debug) cout << "getting puWgts" << endl;
      for(size_t iwgt=0; iwgt<3; iwgt++) {
        puWgts[iwgt]=puWgtGr[iwgt]->Eval(ev.putrue);  
        allPlots["puwgtctr"]->Fill(iwgt+1,puWgts[iwgt]);
      }
      if(debug) cout << "getting puWgts DONE!" << endl;
      //trigger/id+iso efficiency corrections
      if(debug) cout << "calling trigger function" << endl;
      std::vector<int> pdgIds;
      for(size_t ilp = 0; ilp < selLeptons.size(); ilp++)
        pdgIds.push_back(ev.l_id[selLeptons[ilp]]);
      triggerCorrWgt=lepEffH.getTriggerCorrection(pdgIds,leptons);
      if(debug) cout << "calling trigger function DONE!" << endl;
      // ** selLeptons contains only ev_l position, leptons contains p4 **
      for(size_t il=0; il<selLeptons.size(); il++) {
        EffCorrection_t selSF=lepEffH.getOfflineCorrection(pdgIds[il],leptons[il].Pt(),leptons[il].Eta());
        lepSelCorrWgt.second = sqrt( pow(lepSelCorrWgt.first*selSF.second,2)+pow(lepSelCorrWgt.second*selSF.first,2));
        if(debug) cout << "lepSelCorrWgt=" << lepSelCorrWgt.first << endl;
        if(debug) cout << "selSF=" << selSF.first << endl;
        lepSelCorrWgt.first *= selSF.first;
      }

      //update nominal event weight
      //float norm( normH ? normH->GetBinContent(1) : 1.0);
      norm =  normH ? normH->GetBinContent(1) : 1.0;
      //wgt=lepTriggerSF*lepSelSF*puWeight*norm;
      wgt=triggerCorrWgt.first*lepSelCorrWgt.first*puWgts[0]*norm;
      if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
      if(ev.ttbar_nw>0) norm*=ev.ttbar_w[0];
      if(debug) cout << "weight=" << wgt << endl;
      if(debug) cout << "Trigger=" << triggerCorrWgt.first << endl << "Lepton=" << lepSelCorrWgt.first << endl << "PU=" << puWgts[0] << endl << "norm=" << norm  << endl;
      if(filename.Contains("_WJets")) cout << "Trigger=" << triggerCorrWgt.first << endl << "Lepton=" << lepSelCorrWgt.first << endl << "PU=" << puWgts[0] << endl << "norm=" << norm  << endl << "ttbar_w[0]=" << ev.ttbar_w[0] << endl << "wgt=" << wgt << endl;
      for(size_t il = 0; il < leptons.size(); il++) {
        if(!filename.Contains("_WJets")) continue;
        cout << "pT: " << leptons[il].Pt() << endl;
      }
      //wgt=1.0;
    }
    if(debug) cout << "Lepton scale factors DONE!" << endl;

    //sort by Pt
    sort(lightJetsVec.begin(),    lightJetsVec.end(),   sortJetsByPt);
    sort(bJetsVec.begin(),    bJetsVec.end(),   sortJetsByPt);
    sort(allJetsVec.begin(),  allJetsVec.end(), sortJetsByPt);

    for(size_t ij = 0; ij < allJetsVec.size(); ij++) {
      float csv = allJetsVec.at(ij).getCSV();
      allPlots["j_csv"+chTag]->Fill(csv,wgt);
      allPlots["j_csv"+chTag+"_no_weight"]->Fill(csv,norm);
      allPlots["j_csv_all"]->Fill(csv,wgt);
      allPlots["j_csv_all_no_weight"]->Fill(csv,norm);
    }

    if(bJetsVec.size()==0) continue;
    for(size_t il=0; il<leptons.size(); il++) {
      for(size_t ij=0; ij<bJetsVec.size(); ij++) {
        TLorentzVector jp4 = bJetsVec[ij].getVec();
        allPlots["lp_jets_dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
        allPlots["lp_jets_dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        allPlots["lp_jets_dR_all"]->Fill(jp4.DeltaR(leptons[il]),wgt);
        allPlots["lp_jets_dR_all_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
      }
      for(size_t ij=0; ij<lightJetsVec.size(); ij++) {
        TLorentzVector jp4 = lightJetsVec[ij].getVec();
        allPlots["lp_jets_dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
        allPlots["lp_jets_dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        allPlots["lp_jets_dR_all"]->Fill(jp4.DeltaR(leptons[il]),wgt);
        allPlots["lp_jets_dR_all_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
      }
    }


    //MET and transverse mass
    TLorentzVector met(0,0,0,0);
    met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
    met+=jetDiff;
    met.SetPz(0.); met.SetE(met.Pt());

    //simple fill
    bool singleLep(false);
    bool doubleLep(false);
    if(debug) cout << "sorting jets" << endl;
    if(debug) cout << "starting simple plots" << endl;
    allPlots["nevt"+chTag]->Fill(1,wgt);
    allPlots["nevt"+chTag+"_no_weight"]->Fill(1,norm);
    allPlots["nevt_all"]->Fill(1,wgt);
    allPlots["nevt_all_no_weight"]->Fill(1,norm);
    allPlots["j_n"+chTag]->Fill(allJetsVec.size(),wgt);
    allPlots["j_n"+chTag+"_no_weight"]->Fill(allJetsVec.size(),norm);
    allPlots["j_n_all"]->Fill(allJetsVec.size(),wgt);
    allPlots["j_n_all_no_weight"]->Fill(allJetsVec.size(),norm);
    allPlots["lj_n"+chTag]->Fill(lightJetsVec.size(),wgt);
    allPlots["lj_n"+chTag+"_no_weight"]->Fill(lightJetsVec.size(),norm);
    allPlots["lj_n_all"]->Fill(lightJetsVec.size(),wgt);
    allPlots["lj_n_all_no_weight"]->Fill(lightJetsVec.size(),norm);
    allPlots["bj_n"+chTag]->Fill(bJetsVec.size(),wgt);
    allPlots["bj_n"+chTag+"_no_weight"]->Fill(bJetsVec.size(),norm);
    allPlots["bj_n_all"]->Fill(bJetsVec.size(),wgt);
    allPlots["bj_n_all_no_weight"]->Fill(bJetsVec.size(),norm);
    allPlots["lp_n"+chTag]->Fill(selLeptons.size(),wgt);
    allPlots["lp_n"+chTag+"_no_weight"]->Fill(selLeptons.size(),norm);
    allPlots["lp_n_all"]->Fill(selLeptons.size(),wgt);
    allPlots["lp_n_all_no_weight"]->Fill(selLeptons.size(),norm);
    allPlots["pf_n"+chTag]->Fill(ev.npf,wgt);
    allPlots["pf_n"+chTag+"_no_weight"]->Fill(ev.npf,norm);
    allPlots["pf_n_all"]->Fill(ev.npf,wgt);
    allPlots["pf_n_all_no_weight"]->Fill(ev.npf,norm);

    if(allJetsVec.size() > 0) {
      allPlots["j_pt"+chTag]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      allPlots["j_pt"+chTag+"_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
      allPlots["j_pt_all"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      allPlots["j_pt_all_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
    }
    if(lightJetsVec.size() > 0) {
      allPlots["lj_pt"+chTag]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
      allPlots["lj_pt"+chTag+"_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
      allPlots["lj_pt_all"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
      allPlots["lj_pt_all_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
    }
    if(bJetsVec.size() > 0) {
      allPlots["bj_pt"+chTag]->Fill(bJetsVec[0].getVec().Pt(),wgt);
      allPlots["bj_pt"+chTag+"_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
      allPlots["bj_pt_all"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
      allPlots["bj_pt_all_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
    }

    if(selLeptons.size() == 1 && bJetsVec.size() >= 1 && lightJetsVec.size() > 2) {
      singleLep = true;
      std::sort(leptons.begin(), leptons.end(), VecSort);
      allPlots["lp_pt"+chTag]->Fill(leptons[0].Pt(),wgt);
      allPlots["lp_pt"+chTag+"_no_weight"]->Fill(leptons[0].Pt(),norm);
      allPlots["lp_pt_all"]->Fill(leptons[0].Pt(),wgt);
      allPlots["lp_pt_all_no_weight"]->Fill(leptons[0].Pt(),norm);
      allPlots["met"+chTag]->Fill(ev.met_pt[0],wgt);
      allPlots["met"+chTag+"_no_weight"]->Fill(ev.met_pt[0],norm);
      allPlots["met_all"]->Fill(ev.met_pt[0],wgt);
      allPlots["met_all_no_weight"]->Fill(ev.met_pt[0],norm);
      for(int i = 0; i < ev.npf; i++) {
        allPlots["pf_id"+chTag]->Fill(ev.pf_id[i],wgt);
        allPlots["pf_id"+chTag+"_no_weight"]->Fill(ev.pf_id[i],norm);
        allPlots["pf_id_all"]->Fill(ev.pf_id[i],wgt);
        allPlots["pf_id_all_no_weight"]->Fill(ev.pf_id[i],norm);
      }
    }
    else if(selLeptons.size() == 2 && bJetsVec.size() >= 1 && lightJetsVec.size() > 1) {
      if(isZ) {
        allPlots["Z_mass"+chTag]->Fill(dilp4.M(),wgt);
        allPlots["Z_mass"+chTag+"_no_weight"]->Fill(dilp4.M(),norm);
        allPlots["Z_mass_all"]->Fill(dilp4.M(),wgt);
        allPlots["Z_mass_all_no_weight"]->Fill(dilp4.M(),norm);
      }
      if(isZ) continue;
      if(dilp4.M() < 20.) continue;
      if(ev.l_id[selLeptons[0]]==ev.l_id[selLeptons[1]] && met.Pt() < 40) continue; //FIXME
      doubleLep = true;
      allPlots["dilp_pt"+chTag]->Fill(dilp4.Pt(),wgt);
      allPlots["dilp_pt"+chTag+"_no_weight"]->Fill(dilp4.Pt(),norm);
      allPlots["dilp_pt_all"]->Fill(dilp4.Pt(),wgt);
      allPlots["dilp_pt_all_no_weight"]->Fill(dilp4.Pt(),norm);
      allPlots["dilp_mass"+chTag]->Fill(dilp4.M(),wgt);
      allPlots["dilp_mass"+chTag+"_no_weight"]->Fill(dilp4.M(),norm);
      allPlots["dilp_mass_all"]->Fill(dilp4.M(),wgt);
      allPlots["dilp_mass_all_no_weight"]->Fill(dilp4.M(),norm);
      std::sort(leptons.begin(), leptons.end(), VecSort);
      allPlots["lp_pt"+chTag]->Fill(leptons[0].Pt(),wgt);
      allPlots["lp_pt"+chTag+"_no_weight"]->Fill(leptons[0].Pt(),norm);
      allPlots["lp_pt_all"]->Fill(leptons[0].Pt(),wgt);
      allPlots["lp_pt_all_no_weight"]->Fill(leptons[0].Pt(),norm);
      allPlots["lp2_pt"+chTag]->Fill(leptons[1].Pt(),wgt);
      allPlots["lp2_pt"+chTag+"_no_weight"]->Fill(leptons[1].Pt(),norm);
      allPlots["lp2_pt_all"]->Fill(leptons[1].Pt(),wgt);
      allPlots["lp2_pt_all_no_weight"]->Fill(leptons[1].Pt(),norm);
      allPlots["met"+chTag]->Fill(ev.met_pt[0],wgt);
      allPlots["met"+chTag+"_no_weight"]->Fill(ev.met_pt[0],norm);
      allPlots["met_all"]->Fill(ev.met_pt[0],wgt);
      allPlots["met_all_no_weight"]->Fill(ev.met_pt[0],norm);
      allPlots["dilp_charge"+chTag]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
      allPlots["dilp_charge"+chTag+"_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
      allPlots["dilp_charge_all"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
      allPlots["dilp_charge_all_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
      for(int i = 0; i < ev.npf; i++) {
        allPlots["pf_id"+chTag]->Fill(ev.pf_id[i],wgt);
        allPlots["pf_id"+chTag+"_no_weight"]->Fill(ev.pf_id[i],norm);
        allPlots["pf_id_all"]->Fill(ev.pf_id[i],wgt);
        allPlots["pf_id_all_no_weight"]->Fill(ev.pf_id[i],norm);
      }
    }
    if(debug) cout << "simple plots DONE" << endl;


    if(!singleLep && !doubleLep) continue;

    allPlots["nevt"+chTag+"_lep"]->Fill(1,wgt);
    allPlots["nevt"+chTag+"_lep_no_weight"]->Fill(1,norm);
    allPlots["nevt_all_lep"]->Fill(1,wgt);
    allPlots["nevt_all_lep_no_weight"]->Fill(1,norm);
    allPlots["pf_n"+chTag+"_lep"]->Fill(ev.npf,wgt);
    allPlots["pf_n"+chTag+"_lep"+"_no_weight"]->Fill(ev.npf,norm);
    allPlots["pf_n_all_lep"]->Fill(ev.npf,wgt);
    allPlots["pf_n_all_lep_no_weight"]->Fill(ev.npf,norm);

    allPlots["j_n"+chTag+"_lep"]->Fill(allJetsVec.size(),wgt);
    allPlots["j_n"+chTag+"_lep_no_weight"]->Fill(lightJetsVec.size(),norm);
    allPlots["j_n_all_lep"]->Fill(allJetsVec.size(),wgt);
    allPlots["j_n_all_lep_no_weight"]->Fill(lightJetsVec.size(),norm);
    allPlots["lj_n"+chTag+"_lep"]->Fill(lightJetsVec.size(),wgt);
    allPlots["lj_n"+chTag+"_lep_no_weight"]->Fill(lightJetsVec.size(),norm);
    allPlots["lj_n_all_lep"]->Fill(lightJetsVec.size(),wgt);
    allPlots["lj_n_all_lep_no_weight"]->Fill(lightJetsVec.size(),norm);
    allPlots["bj_n"+chTag+"_lep"]->Fill(bJetsVec.size(),wgt);
    allPlots["bj_n"+chTag+"_lep_no_weight"]->Fill(bJetsVec.size(),norm);
    allPlots["bj_n_all_lep"]->Fill(bJetsVec.size(),wgt);
    allPlots["bj_n_all_lep_no_weight"]->Fill(bJetsVec.size(),norm);
    allPlots["lp_n"+chTag+"_lep"]->Fill(selLeptons.size(),wgt);
    allPlots["lp_n"+chTag+"_lep_no_weight"]->Fill(selLeptons.size(),norm);
    allPlots["lp_n_all_lep"]->Fill(selLeptons.size(),wgt);
    allPlots["lp_n_all_lep_no_weight"]->Fill(selLeptons.size(),norm);

    allPlots["j_pt"+chTag+"_lep"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    allPlots["j_pt"+chTag+"_lep_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
    allPlots["j_pt_all_lep"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    allPlots["j_pt_all_lep_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
    allPlots["lj_pt"+chTag+"_lep"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
    allPlots["lj_pt"+chTag+"_lep_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
    allPlots["lj_pt_all_lep"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
    allPlots["lj_pt_all_lep_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
    allPlots["bj_pt"+chTag+"_lep"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
    allPlots["bj_pt"+chTag+"_lep_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
    allPlots["bj_pt_all_lep"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
    allPlots["bj_pt_all_lep_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);

    for(size_t ij = 0; ij < allJetsVec.size(); ij++) {
      float csv = allJetsVec.at(ij).getCSV();
      allPlots["j_csv"+chTag+"_lep"]->Fill(csv,wgt);
      allPlots["j_csv"+chTag+"_lep_no_weight"]->Fill(csv,norm);
      allPlots["j_csv_all_lep"]->Fill(csv,wgt);
      allPlots["j_csv_all_lep_no_weight"]->Fill(csv,norm);
    }

    //charmed resonance analysis : use only jets with CSV>CSVL, up to two per event
    for(size_t ij = 0; ij < bJetsVec.size(); ij++) {

      if(ij > 1) continue;
      if(ij == 0) {
        allPlots["nevt"+chTag+"_csv"]->Fill(1,wgt);
        allPlots["nevt"+chTag+"_csv_no_weight"]->Fill(1,norm);
        allPlots["nevt_all_csv"]->Fill(1,wgt);
        allPlots["nevt_all_csv_no_weight"]->Fill(1,norm);
        allPlots["j_pt"+chTag+"_csv"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
        allPlots["j_pt"+chTag+"_csv_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
        allPlots["j_pt_all_csv"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
        allPlots["j_pt_all_csv_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
        allPlots["lj_pt"+chTag+"_csv"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
        allPlots["lj_pt"+chTag+"_csv_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
        allPlots["lj_pt_all_csv"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
        allPlots["lj_pt_all_csv_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
        allPlots["bj_pt"+chTag+"_csv"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
        allPlots["bj_pt"+chTag+"_csv_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
        allPlots["bj_pt_all_csv"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
        allPlots["bj_pt_all_csv_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
      }

      allPlots["j_csv"+chTag+"_csv"]->Fill(bJetsVec[ij].getCSV(),wgt);
      allPlots["j_csv"+chTag+"_csv_no_weight"]->Fill(bJetsVec[ij].getCSV(),norm);
      allPlots["j_csv_all_csv"]->Fill(bJetsVec[ij].getCSV(),wgt);
      allPlots["j_csv_all_csv_no_weight"]->Fill(bJetsVec[ij].getCSV(),norm);
      std::vector<IdTrack> &tracks = bJetsVec[ij].getTracks();

      //J/Psi
      if(debug) cout << "starting J/Psi" << endl;
      const float gMassMu(0.1057),gMassK(0.4937),gMassPi(0.1396);
      std::vector<pfTrack> pfmuCands,kaonCands;
      for(size_t itk = 0; itk < tracks.size(); itk++) {
        if(abs(tracks[itk].second) == 13) {
          TLorentzVector muP4;
          muP4.SetPtEtaPhiM(tracks[itk].first.Pt(), tracks[itk].first.Eta(), tracks[itk].first.Phi(), gMassMu);
          pfTrack pfmu(muP4, tracks[itk].first.getDxy(), tracks[itk].first.getDxyE(), tracks[itk].first.getDz(), tracks[itk].first.getDzE(), tracks[itk].second);
          pfmuCands.push_back(pfmu);
        }
        if(abs(tracks[itk].second) == 211) {
          TLorentzVector kP4;
          pfTrack pfk(kP4, tracks[itk].first.getDxy(), tracks[itk].first.getDxyE(), tracks[itk].first.getDz(), tracks[itk].first.getDzE(), tracks[itk].second);
          kaonCands.push_back(pfk);
        }
      }

      if(pfmuCands.size()>1) {
        if(pfmuCands[0].getPfid() != -pfmuCands[1].getPfid()) continue;
        float mass12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).M());
        float mass123( kaonCands.size()>0 ? (pfmuCands[0].getVec()+pfmuCands[1].getVec()+kaonCands[0].getVec()).M() : -1);
        allPlots["jpsi_mass"+chTag+"_jpsi"]->Fill(mass12,wgt);
        allPlots["jpsi_mass"+chTag+"_jpsi_no_weight"]->Fill(mass12,norm);
        allPlots["jpsi_mass_all_jpsi"]->Fill(mass12,wgt);
        allPlots["jpsi_mass_all_jpsi_no_weight"]->Fill(mass12,norm);
        allPlots["bj_pt"+chTag+"_jpsi"]->Fill(bJetsVec[ij].getVec().Pt(),wgt);
        allPlots["bj_pt"+chTag+"_jpsi_no_weight"]->Fill(bJetsVec[ij].getVec().Pt(),norm);
        allPlots["bj_pt_all_jpsi"]->Fill(bJetsVec[ij].getVec().Pt(),wgt);
        allPlots["bj_pt_all_jpsi_no_weight"]->Fill(bJetsVec[ij].getVec().Pt(),norm);
        allPlots["j_csv"+chTag+"_jpsi"]->Fill(bJetsVec[ij].getCSV(),wgt);
        allPlots["j_csv"+chTag+"_jpsi_no_weight"]->Fill(bJetsVec[ij].getCSV(),norm);
        allPlots["j_csv_all_jpsi"]->Fill(bJetsVec[ij].getCSV(),wgt);
        allPlots["j_csv_all_jpsi_no_weight"]->Fill(bJetsVec[ij].getCSV(),norm);
        allPlots["nevt"+chTag+"_jpsi"]->Fill(1,wgt);
        allPlots["nevt"+chTag+"_jpsi_no_weight"]->Fill(1,norm);
        allPlots["nevt_all_jpsi"]->Fill(1,wgt);
        allPlots["nevt_all_jpsi_no_weight"]->Fill(1,norm);
        if(mass12<3.0 || mass12>3.2) continue;
        for(int itk = 0; itk < 2; itk++) {

          for(int i = 0; i < 2; i++) {
            allPlots["pf_dxy"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDxy()),wgt);
            allPlots["pf_dz"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDz()),wgt);
            allPlots["pf_dxyE"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDxyE()),wgt);
            allPlots["pf_dzE"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDzE()),wgt);
            allPlots["pf_dz_sig"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
            allPlots["pf_dxy_sig"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDxy())/abs(pfmuCands[i].getDxyE()),wgt);
            allPlots["pf_dz_sig"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
            allPlots["pf_dxy"+chTag+"_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDxy()),norm);
            allPlots["pf_dz"+chTag+"_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDz()),norm);
            allPlots["pf_dxyE"+chTag+"_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDxyE()),norm);
            allPlots["pf_dzE"+chTag+"_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDzE()),norm);
            allPlots["pf_dz_sig"+chTag+"_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),norm);
            allPlots["pf_dxy_sig"+chTag+"_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDxy())/abs(pfmuCands[i].getDxyE()),norm);
            allPlots["pf_dz_sig"+chTag+"_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),norm);
            allPlots["pf_dxy_all_jpsi"]->Fill(abs(pfmuCands[i].getDxy()),wgt);
            allPlots["pf_dz_all_jpsi"]->Fill(abs(pfmuCands[i].getDz()),wgt);
            allPlots["pf_dxyE_all_jpsi"]->Fill(abs(pfmuCands[i].getDxyE()),wgt);
            allPlots["pf_dzE_all_jpsi"]->Fill(abs(pfmuCands[i].getDzE()),wgt);
            allPlots["pf_dz_sig_all_jpsi"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
            allPlots["pf_dxy_sig_all_jpsi"]->Fill(abs(pfmuCands[i].getDxy())/abs(pfmuCands[i].getDxyE()),wgt);
            allPlots["pf_dz_sig_all_jpsi"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
            allPlots["pf_dxy_all_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDxy()),norm);
            allPlots["pf_dz_all_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDz()),norm);
            allPlots["pf_dxyE_all_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDxyE()),norm);
            allPlots["pf_dzE_all_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDzE()),norm);
            allPlots["pf_dz_sig_all_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),norm);
            allPlots["pf_dxy_sig_all_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDxy())/abs(pfmuCands[i].getDxyE()),norm);
            allPlots["pf_dz_sig_all_jpsi_no_weight"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),norm);
          }
        }

        if(filename.Contains("_WJets"))
          cout << endl << mass12 << " " << wgt << endl;
        if(mass123 > 0) {
          allPlots["jpsi_kaon_mass"+chTag+"_jpsi"]->Fill(mass123,wgt);
          allPlots["jpsi_kaon_mass"+chTag+"_jpsi_no_weight"]->Fill(mass123,norm);
          allPlots["jpsi_kaon_mass_all_jpsi"]->Fill(mass123,wgt);
          allPlots["jpsi_kaon_mass_all_jpsi_no_weight"]->Fill(mass123,norm);
        }
      }
      if(debug) cout << "J/Psi DONE" << endl;

      //D0 and D* 
      if(debug) cout << "Starting D0 and D*" << endl;
      int jetindex = allJetsVec[ij].getJetIndex();
      if(tracks.size() < 3) continue;
      size_t tmax = 4;
      tmax = tracks.size() >= tmax ? tmax : tracks.size();
      for(size_t i = 0; i < tmax; i++)
        for(size_t j = 0; j < tmax; j++) {
          if(i == j) continue;

          //opposite sign
          if(tracks[i].second*tracks[j].second != -211*211) continue;

          TLorentzVector p_track1, p_track2;
          p_track1.SetPtEtaPhiM(tracks[i].first.Pt(), tracks[i].first.Eta(), tracks[i].first.Phi(), gMassPi);
          p_track2.SetPtEtaPhiM(tracks[j].first.Pt(), tracks[j].first.Eta(), tracks[j].first.Phi(), gMassK);
          if(debug) cout << i << ": " << tracks[i].first.Pt() << " " << tracks[i].first.Eta() << " " << tracks[i].first.Phi() << " " << gMassPi << endl;
          if(debug) cout << j << ": " << tracks[j].first.Pt() << " " << tracks[j].first.Eta() << " " << tracks[j].first.Phi() << " " << gMassK << endl << endl;
          float mass12 = (p_track1+p_track2).M();
          if(debug) cout << mass12 << endl;

          if (mass12 > 1.7 && mass12 < 2.) {
            allPlots["d0_mass_large"+chTag+"_meson"]->Fill(mass12,wgt);
            allPlots["d0_mass_large"+chTag+"_meson_no_weight"]->Fill(mass12,norm);
            allPlots["d0_mass_large_all_meson"]->Fill(mass12,wgt);
            allPlots["d0_mass_large_all_meson_no_weight"]->Fill(mass12,norm);
            if (mass12 > 1.85 && mass12 < 2.15) {
              allPlots["d0_mass"+chTag+"_meson"]->Fill(mass12,wgt);
              allPlots["d0_mass"+chTag+"_meson_no_weight"]->Fill(mass12,norm);
              allPlots["d0_mass_all_meson"]->Fill(mass12,wgt);
              allPlots["d0_mass_all_meson_no_weight"]->Fill(mass12,norm);
            }
          }

          //looking for lepton
          if(debug) cout << "third lepton" << endl;
          for(size_t k = 0; k < tracks.size(); k++) {
            if(k == i) continue;
            if(k == j) continue;
            if(debug) cout << "third lepton possible" << endl;

            if(abs(tracks[k].second) != 13 && abs(tracks[k].second) != 11) continue;
            if(debug) cout << "third lepton found" << endl;

            if(tracks[j].second/abs(tracks[j].second) == -tracks[k].second/abs(tracks[k].second)) {
              //Kaon and lepton have same charge
              //correct mass assumption
              if (mass12 > 1.85 && mass12 < 2.15) {
                if(debug) cout << "correct mass assumption" << endl;
                allPlots["d0_mass_lp"+chTag+"_meson"]->Fill(mass12,wgt);
                allPlots["d0_mass_lp"+chTag+"_meson_no_weight"]->Fill(mass12,norm);
                allPlots["d0_mass_lp_all_meson"]->Fill(mass12,wgt);
                allPlots["d0_mass_lp_all_no_weight"]->Fill(mass12,norm);
                if(abs(tracks[k].second) == 13) {
                  allPlots["d0_mass_mu"+chTag+"_meson"]->Fill(mass12,wgt);
                  allPlots["d0_mass_mu"+chTag+"_meson_no_weight"]->Fill(mass12,norm);
                  allPlots["d0_mass_mu_all_meson"]->Fill(mass12,wgt);
                  allPlots["d0_mass_mu_all_meson_no_weight"]->Fill(mass12,norm);
                }
                if(abs(tracks[k].second) == 11) {
                  allPlots["d0_mass_el"+chTag+"_meson"]->Fill(mass12,wgt);
                  allPlots["d0_mass_el"+chTag+"_meson_no_weight"]->Fill(mass12,norm);
                  allPlots["d0_mass_el_all_meson"]->Fill(mass12,wgt);
                  allPlots["d0_mass_el_all_meson_no_weight"]->Fill(mass12,norm);
                }
              }
            }
          }
          //looking for pion
          if(debug) cout << "D*->pi+D0" << endl;
          for(size_t k = 0; k < tracks.size(); k++) {
            if(k == i) continue;
            if(k == j) continue;

            if(abs(tracks[k].second) != 211) continue;
            if(debug) cout << "Pion found" << endl;

            TLorentzVector p_track3, p_cand;
            p_track3.SetPtEtaPhiM(tracks[k].first.Pt(), tracks[k].first.Eta(), tracks[k].first.Phi(), gMassPi);
            if(debug) cout << k << ": " << tracks[k].first.Pt() << " " << tracks[k].first.Eta() << " " << tracks[k].first.Phi() << " " << gMassPi << endl;
            allPlots["pi_pt"+chTag+"_meson"]->Fill(p_track3.Pt(),wgt);
            allPlots["pi_pt"+chTag+"_meson_no_weight"]->Fill(p_track3.Pt(),norm);
            allPlots["pi_pt_all_meson"]->Fill(p_track3.Pt(),wgt);
            allPlots["pi_pt_all_meson_no_weight"]->Fill(p_track3.Pt(),norm);
            if( tracks[j].second/abs(tracks[j].second) == -tracks[k].second/abs(tracks[k].second) ) {
              // Kaon and pion have opposite charges
              // I.e. correct mass assumption
              if(debug) cout << "correct mass assumption" << endl;

              p_cand = p_track1+p_track2+p_track3;
              allPlots["ds_mass"+chTag+"_meson"]->Fill(p_cand.M(), wgt);
              allPlots["ds_mass"+chTag+"_meson_no_weight"]->Fill(p_cand.M(),norm);
              allPlots["ds_mass_all"]->Fill(p_cand.M(), wgt);
              allPlots["ds_mass_all_meson_no_weight"]->Fill(p_cand.M(),norm);

              if(abs(mass12-1.864) < 0.10) { // mass window cut
                TLorentzVector p_jet;
                p_jet.SetPtEtaPhiM(ev.j_pt[jetindex], ev.j_eta[jetindex], ev.j_phi[jetindex], 0.);

                float deltam = p_cand.M() - mass12;

                allPlots["ds_d0_dmass_loose"+chTag+"_meson"]->Fill(deltam, wgt);
                allPlots["ds_d0_dmass_loose"+chTag+"_meson_no_weight"]->Fill(deltam,norm);
                allPlots["ds_d0_dmass_loose_all_meson"]->Fill(deltam, wgt);
                allPlots["ds_d0_dmass_loose_all_meson_no_weight"]->Fill(deltam, norm);
                if(abs(mass12-1.864) < 0.05) { // tighter mass window cut
                  if(filename.Contains("_WJets")) {
                    cout << endl << deltam << " " << wgt << endl;
                    cout << "pi1: " << p_track1.Pt() << endl;
                    cout << "K: " << p_track2.Pt() << endl;
                    cout << "pi2: " << p_track3.Pt() << endl;
                  }

                  allPlots["ds_d0_dmass"+chTag+"_meson"]->Fill(deltam, wgt);
                  allPlots["ds_d0_dmass"+chTag+"_meson_no_weight"]->Fill(deltam, norm);
                  allPlots["ds_d0_dmass_all"]->Fill(deltam, wgt);
                  allPlots["ds_d0_dmass_all_meson_no_weight"]->Fill(deltam, norm);
                  if(deltam<0.14 || deltam>0.15) continue;
                  allPlots["nevt"+chTag+"_meson"]->Fill(1,wgt);
                  allPlots["nevt"+chTag+"_meson_no_weight"]->Fill(1,norm);
                  allPlots["nevt_all_meson"]->Fill(1,wgt);
                  allPlots["nevt_all_meson_no_weight"]->Fill(1,norm);
                }
              }
            }
          }
        }
      if(debug) cout << "D0 and D* DONE" << endl;
    }

  }

  //close input file
  f->Close();

  //save histos to file  
  TString selPrefix("");  
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  fOut->cd();
  if(debug) cout << "writing histograms" << endl;

  for (auto& it : allPlots)  { 
    if(debug) cout << it.second->GetName() << endl;
    if(debug) cout << it.second->GetEntries() << endl;

    //fOut->cd( dir );
    it.second->SetDirectory(fOut); it.second->Write(); 
    fOut->cd();
  }
  if(debug) cout << "writing histograms DONE" << endl;
  if(debug) cout << "closing ROOT file" << endl;
  fOut->Close();
}

