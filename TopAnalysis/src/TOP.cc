#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"
#include "TopLJets2015/TopAnalysis/interface/TOPWidth.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "TopLJets2015/TopAnalysis/interface/OtherFunctions.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

using namespace std;

//
/*
Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}


*/
//
void RunTop(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
                 Bool_t debug=false)
{
  if(debug) cout << "in RunTop" << endl;

  bool isTTbar( filename.Contains("_TTJets") );
  //bool debug(false);
  //bool isData( filename.Contains("Data") ? true : false);
  
  //READ TREE FROM FILE
  MiniEvent_t ev;
  //TopWidthEvent_t ev;
  TFile *f = TFile::Open(filename);
  //TTree *t = (TTree*)f->Get("twev");
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev);
  //createTopWidthEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  bool requireEtriggerOnly(false);
  if(ev.isData) runSysts=false;
  if(ev.isData && filename.Contains("SingleElectron")) requireEtriggerOnly=true;
  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;
  
  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
    {
      if(debug) cout << "loading pileup weight" << endl;
      TString puWgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      TGraph *puData=(TGraph *)fIn->Get("pu_nom");
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
      gr->SetName("puwgts_nom");
      puWgtGr.push_back( gr );
      tmp->Delete();
    }
    if(debug) cout << "loading pileup weight DONE" << endl;

  //LEPTON EFFICIENCIES
  TString lepEffUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/leptonEfficiencies.root");
  gSystem->ExpandPathName(lepEffUrl);
  std::map<TString,TH2 *> lepEffH;
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["m_sel"]=(TH2 *)fIn->Get("m_sel");
      lepEffH["m_trig"]=(TH2 *)fIn->Get("m_trig");      
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }

  lepEffUrl="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root";
  gSystem->ExpandPathName(lepEffUrl);
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }

  //B-TAG CALIBRATION
  TString btagUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  TString btagEffExpUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff;
  BTagSFUtil myBTagSFUtil;
  float wgt(1.0);
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
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
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

      //wgt=1.0;
      //float norm( normH ? normH->GetBinContent(1) : 1.0);
      //wgt=norm;//lepTriggerSF[0]*lepSelSF[0]*puWeight[0]*norm;
    }

  //JET ENERGY SCALE: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Summer15_uncertainties
  
  //LIST OF SYSTEMATICS
  
  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;

  std::vector<TString> lfsVec = { "_e", "_ee", "_em", "_mm", "_m" }; 

  for(int i = 0; i < (int)lfsVec.size(); i++) {
    TString tag(lfsVec[i]);
    allPlots["lp_pt"+tag] = new TH1F("lp_pt"+tag,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["l2p_pt"+tag] = new TH1F("l2p_pt"+tag,";Sub-leading Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_pt"+tag] = new TH1F("dilp_pt"+tag,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag] = new TH1F("dilp_m"+tag,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["j_pt"+tag] = new TH1F("j_pt"+tag,";Leading Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["bj_pt"+tag] = new TH1F("bj_pt"+tag,";Leading b Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["nlp"+tag]     = new TH1F("nlp"+tag,";N_{l};Events" ,3,0.,3.);
    allPlots["ndilp"+tag]     = new TH1F("ndilp"+tag,";N_{ll};Events" ,3,0.,3.);
    allPlots["nj"+tag]     = new TH1F("nj"+tag,";N_{jets} (P_{T} > 30 GeV);Events" ,8,2.,10.);
    allPlots["nbj"+tag]     = new TH1F("nbj"+tag,";N_{b-jets} (CSV > 0.8);Events" ,4,1.,5.);
    allPlots["npf"+tag]     = new TH1F("npf"+tag,";N_{pf};Events / 10" ,50,0.,500.);
    allPlots["nstart"+tag]     = new TH1F("nstart"+tag,";N_{start};Events" ,5,0.,5.);
    allPlots["pfid"+tag]     = new TH1F("pfid"+tag,";PFID;Events" ,440,-220.,220.);
/*
    allPlots["massJPsi"+tag]     = new TH1F("massJPsi"+tag,";M_{J/#Psi};Events" ,20,2.,4.);
    allPlots["massD0"+tag]     = new TH1F("massD0"+tag,";M_{jj};Events" ,20,1.,3.);
    allPlots["masslep"+tag]     = new TH1F("masslep"+tag,";M_{K#pi};Events" ,20,0.,10.);
    allPlots["massmu"+tag]     = new TH1F("massmu"+tag,";M_{K#pi};Events" ,20,0.,10.);
    allPlots["masse"+tag]     = new TH1F("masse"+tag,";M_{K#pi};Events" ,20,0.,10.);
    allPlots["massDsmD0loose"+tag]     = new TH1F("massDsmD0loose"+tag,";M_{K#pi};Events" ,20,1.,3.);
    allPlots["massDsmD0"+tag]     = new TH1F("massDsmD0"+tag,";M_{K#pi};Events" ,20,1.,3.);
    allPlots["massDs"+tag]     = new TH1F("massDs"+tag,";M_{D^{*}};Events" ,20,0.,20.);
*/
    allPlots["massJPsi"+tag]     = new TH1F("massJPsi"+tag,";M_{J/#Psi};Events / 0.01 GeV" ,18,2.5,3.4);
    allPlots["massD0"+tag]     = new TH1F("massD0"+tag,";M_{D^{0}};Events / 0.01 GeV" ,30,1.7,2.0);
    allPlots["masslep"+tag]     = new TH1F("masslep"+tag,";M_{K#pi};Events" ,20,0.,20.);
    allPlots["massmu"+tag]     = new TH1F("massmu"+tag,";M_{K#pi};Events" ,20,0.,20.);
    allPlots["masse"+tag]     = new TH1F("masse"+tag,";M_{K#pi};Events" ,20,0.,20.);
    allPlots["massDsmD0loose"+tag]     = new TH1F("massDsmD0loose"+tag,";M_{K#pi};Events / 0.05 GeV" ,6,1.7,2.);
    allPlots["massDsmD0"+tag]     = new TH1F("massDsmD0"+tag,";M_{K#pi};Events / 0.05 GeV" ,6,1.7,2.);
    allPlots["massDs"+tag]     = new TH1F("massDs"+tag,";M_{D^{*}};Events / 0.1 GeV / 0.01" ,30,1.7,2.0);
    allPlots["pi_pt"+tag] = new TH1F("pi_pt"+tag,";#pi^{#pm} P_{T} [GeV];Events / 10 GeV", 15, 0,150);
    allPlots["MET"+tag] = new TH1F("MET"+tag,";MET [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag] = new TH1F("charge"+tag,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["dR"+tag] = new TH1F("dR"+tag,";dR;Events / 0.03", 20,0.4,1.);
    allPlots["lep_pt"+tag] = new TH1F("lep_pt"+tag,";PF lepton P_{T} [GeV];Events / 20 GeV", 15, 0,300);

  }
    allPlots["relIso_m"] = new TH1F("relIso_m",";relIso;Events / 0.05", 20,0,1.);
    allPlots["relIso_e"] = new TH1F("relIso_e",";relIso;Events / 0.05", 20,0,1.);


  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  //for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }
  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //account for pu weights and effect on normalization
      float puWeight(1.0);
      if(!ev.isData) 
	{
	  puWeight=puWgtGr[0]->Eval(ev.putrue);  
          /*
	  allPlots["puwgtctr"]->Fill(0.,1.0);
	  allPlots["puwgtctr"]->Fill(1.,puWeight);
          */
	}

      //select 1 good lepton
      //cout << "entering lepton selection" << endl;
      std::vector<int> tightLeptonsNonIso, vetoLeptons;
      std::vector<int> tightLeptons,looseLeptons;
      for(int il=0; il<ev.nl; il++)
	{
          //cout << "in lepton selection" << endl;
	  bool passTightKin(ev.l_pt[il]>30 && fabs(ev.l_eta[il])<2.1);
	  bool passVetoKin(ev.l_pt[il]>10 && fabs(ev.l_eta[il])<2.5);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
	  float relIso(ev.l_relIso[il]);
	  bool passIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 );
	  bool passNonIso(relIso>0.4);
	  if( ev.l_id[il]==11 && (passIso || relIso<0.4) ) passNonIso=false;
	  bool passVetoIso(  ev.l_id[il]==13 ? relIso<0.25 : true); 
	  //bool passSIP3d(ev.l_ip3dsig[il]<4);
	  //if(channelSelection==21) passSIP3d=true;
	  if( ev.l_id[il] == 13 )
            allPlots["relIso_m"]->Fill(relIso,wgt);
          else if( ev.l_id[il] == 11)
            allPlots["relIso_e"]->Fill(relIso,wgt);

	  if(passTightKin && passTightId)// && passSIP3d)
	    {
	      if(passIso)         tightLeptons.push_back(il);
	      else if(passNonIso) tightLeptonsNonIso.push_back(il);
	    }
	  else if(passVetoKin && passVetoIso) vetoLeptons.push_back(il);
	}
      if(debug) cout << "lepton selection DONE" << endl;

      //check if triggers have fired
      bool hasMuTrigger((ev.muTrigger & 0x3)!=0);
      bool hasEleTrigger((ev.elTrigger & 0x1)!=0);

      //decide the channel
      if(debug) cout << "decide channel" << endl;
      TString chTag("");
      std::vector<int> selLeptons;
      if(tightLeptons.size()==1 && looseLeptons.size()==0)
	{
	  selLeptons.push_back( tightLeptons[0] );
          if(debug) cout << "found 1 tight lepton" << endl;
	}
      if(tightLeptons.size()>=2)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(tightLeptons[1]);
          if(debug) cout << "found 2 tight leptons" << endl;
	}
      if(tightLeptons.size()==1 && looseLeptons.size()>=1)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(looseLeptons[0]);	  
	}
      if(debug) if(selLeptons.size()==0) cout << "NO LEPTONS!!" << endl;
      if(selLeptons.size()==0) continue;
      if(selLeptons.size()==1)
	{
	  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="e";
	  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger)  chTag="m";
	}
      if(selLeptons.size()==2)
	{
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*11 && hasEleTrigger) chTag="ee";
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==13*13 && hasMuTrigger) chTag="mm";
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*13)
	    {
	      if(tightLeptons.size()>=2 && (hasEleTrigger || hasMuTrigger)) chTag="em";
	      if(tightLeptons.size()==1)
		{
		  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="em";
		  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger) chTag="em";
		}
	    }
	  if(hasMuTrigger && requireEtriggerOnly) chTag="";
	}
      if(chTag=="") continue;
      chTag = "_"+chTag;
      if(debug) cout << "decide channel DONE" << endl;

      //one good lepton either isolated or in the non-isolated sideband or a Z candidate
      Int_t lepIdx=-1;
      Bool_t isZ(false);//,isZPassingSIP3d(false);
      TLorentzVector l1p4,l2p4,dilp4;
      if(tightLeptons.size()==1)                                       lepIdx=tightLeptons[0];
      else if (tightLeptons.size()==0 && tightLeptonsNonIso.size()==1) lepIdx=tightLeptonsNonIso[0];
      else if(tightLeptons.size()==2)
	{	  
          if(debug) cout << "di-lepton" << endl;
	  l1p4.SetPtEtaPhiM(ev.l_pt[tightLeptons[0]],ev.l_eta[tightLeptons[0]],ev.l_phi[tightLeptons[0]],ev.l_mass[tightLeptons[0]]);
	  l2p4.SetPtEtaPhiM(ev.l_pt[tightLeptons[1]],ev.l_eta[tightLeptons[1]],ev.l_phi[tightLeptons[1]],ev.l_mass[tightLeptons[1]]);
	  dilp4=l1p4+l2p4;
	  if(ev.l_id[tightLeptons[0]]==ev.l_id[tightLeptons[1]]          && 
	     ev.l_charge[tightLeptons[0]]*ev.l_charge[tightLeptons[1]]<0 && 
	     fabs(dilp4.M()-91)<10) //&& 
	     //dilp4.Pt()>30)
	    { 
	      isZ=true; 
	      //isZPassingSIP3d=(ev.l_ip3dsig[0]<4 && ev.l_ip3dsig[1]<4);
	    }
	  lepIdx=tightLeptons[0];
          if(debug) cout << "di-lepton DONE" << endl;
	}

      if(lepIdx<0) continue;
      
      //no extra isolated leptons
      if(vetoLeptons.size()>0) continue;
      
      //apply trigger requirement
      /*
      if(ev.l_id[lepIdx]==13)
	{
	  if(ev.isData  && (ev.muTrigger & 0x3)==0) continue;
	  if(!ev.isData && (ev.muTrigger & 0x3)==0) continue;
	}
      if(ev.l_id[lepIdx]==11)
	{ 
	  if( ((ev.elTrigger>>0)&0x1)==0 ) continue;
	}

      //select according to the lepton id/charge
      Int_t lid=ev.l_id[lepIdx];
      if(isZ) lid=2100+ev.l_id[lepIdx];
      else if(tightLeptonsNonIso.size()==1) lid=100*ev.l_id[lepIdx];
      if(channelSelection!=0)
	{
	  if(channelSelection==21) { if(!isZ) continue; }
	  else                     { if(lid!=channelSelection) continue; }
	}
      if(chargeSelection!=0  && ev.l_charge[lepIdx]!=chargeSelection) continue;
      */

      //lepton kinematics
      if(debug) cout << "checking lepton kinematics" << endl;
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);

      if(debug) cout << "checking lepton kinematics DONE" << endl;

      //save lepton kinematics
      std::vector<TLorentzVector> leptons;
      for(size_t il=0; il<selLeptons.size(); il++)
	{
	  int lepIdx=selLeptons[il];
	  TLorentzVector lp4;
	  lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);
	  leptons.push_back(lp4);
	}

      //select jets
      Float_t htsum(0);
      TLorentzVector jetDiff(0,0,0,0);
      std::vector<TLorentzVector> bJets,lightJets;
      TLorentzVector visSystem(isZ ? dilp4 : lp4);
      int nbjets(0),ncjets(0),nljets(0);//,leadingJetIdx(-wgt);
      std::vector<int> resolvedJetIdx;
      std::vector<TLorentzVector> resolvedJetP4;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  //jp4=updateJES(jp4,ev.j_rawsf[k],ev.j_area[k],ev.rho,ev.nvtx,jetCorr);

	  //cross clean with respect to leptons 
	  //if(jp4.DeltaR(lp4)<0.5) continue;
	  //if(isZ && jp4.DeltaR(l2p4)<0.5)continue;
          bool overlapsWithLepton(false);
          for(size_t il=0; il<leptons.size(); il++) {
            //allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
            if(jp4.DeltaR(leptons[il])>0.4) continue;
	    overlapsWithLepton=true;
          }
          if(overlapsWithLepton) continue;

	  //smear jet energy resolution for MC
	  //jetDiff -= jp4;
          /*
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
          */
          /*
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }
          */
	  //jetDiff += jp4;
	  resolvedJetIdx.push_back(k);
	  resolvedJetP4.push_back(jp4);

	  //require back-to-back configuration with Z
	  if(isZ && jp4.DeltaPhi(dilp4)<2.7) continue;

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  //if(leadingJetIdx<0) leadingJetIdx=k;
	  htsum += jp4.Pt();
	  if(bJets.size()+lightJets.size()<4) visSystem += jp4;

	  //b-tag
	  float csv = ev.j_csv[k];	  
	  bool isBTagged(csv>0.800);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 
		  ncjets++;
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  nbjets++;
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		}
	      else
		{
		  nljets++;
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		}
	      
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	    }

	  //save jet
	  if(isBTagged) bJets.push_back(jp4);
	  else          lightJets.push_back(jp4);
	}

      
      //event weight
      //float wgt(1.0);
      if(!ev.isData)
	{
	  //update lepton selection scale factors, if found
	  float lepTriggerSF(1.0),lepSelSF(1.0);
	  for(UInt_t il=0; il<leptons.size(); il++)
	    {
	      TString prefix(abs(ev.l_id[il])==11 ? "e" : "m");
	      float minEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmax()-0.01 );
	      float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[il].Eta())),maxEtaForEff),minEtaForEff);
	      Int_t etaBinForEff=lepEffH[prefix+"_sel"]->GetXaxis()->FindBin(etaForEff);
	      
	      float minPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmin() ), maxPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmax()-0.01 );
	      float ptForEff=TMath::Max(TMath::Min(float(leptons[il].Pt()),maxPtForEff),minPtForEff);
	      Int_t ptBinForEff=lepEffH[prefix+"_sel"]->GetYaxis()->FindBin(ptForEff);
		  		  
	      lepSelSF=(lepEffH[prefix+"_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
	      
	      if(il!=0) continue;
	      if(prefix=="m")
		{
		  lepTriggerSF=(lepEffH[prefix+"_trig"]->GetBinContent(etaBinForEff,ptBinForEff));
		
		}
	    }

	  //update nominal event weight
	  float norm( normH ? normH->GetBinContent(1) : 1.0);
	  wgt=lepTriggerSF*lepSelSF*puWeight*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
          //wgt=1.0;
	}

      for(size_t il=0; il<leptons.size(); il++) {
        for(size_t ij=0; ij<bJets.size(); ij++) {
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[ij],ev.j_eta[ij],ev.j_phi[ij],ev.j_mass[ij]);
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
        }
        for(size_t ij=0; ij<lightJets.size(); ij++) {
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[ij],ev.j_eta[ij],ev.j_phi[ij],ev.j_mass[ij]);
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
        }
      }


      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      met+=jetDiff;
      met.SetPz(0.); met.SetE(met.Pt());
      //float mt( computeMT(isZ ? dilp4: lp4,met) );

      //simple fill
      bool singleLep(false);
      bool doubleLep(false);
      if(debug) cout << "starting simple plots" << endl;
      if(selLeptons.size() == 1 && bJets.size() >= 1 && lightJets.size() >= 4) {
        singleLep = true;
        allPlots["nj"+chTag]->Fill(lightJets.size(),wgt);
        allPlots["nbj"+chTag]->Fill(bJets.size(),wgt);
        allPlots["nlp"+chTag]->Fill(selLeptons.size(),wgt);
        allPlots["lp_pt"+chTag]->Fill(lp4.Pt(),wgt);
        if(debug) cout << "sorting jets" << endl;
	std::sort(lightJets.begin(), lightJets.end(), VecSort);
        std::sort(bJets.begin(), bJets.end(), VecSort);
        /*
        for (auto it : lightJets)
          allPlots["j_pt"+chTag]->Fill(it.Pt(),wgt);
        for (auto it : bJets)
          allPlots["bj_pt"+chTag]->Fill(it.Pt(),wgt);
        */
        allPlots["j_pt"+chTag]->Fill(lightJets[0].Pt(),wgt);
        allPlots["bj_pt"+chTag]->Fill(bJets[0].Pt(),wgt);
        allPlots["MET"+chTag]->Fill(ev.met_pt[0],wgt);
      }
      else if(selLeptons.size() == 2 && bJets.size() >= 1 && lightJets.size() > 1) {
        if(isZ) continue;
        if(dilp4.M() < 10) continue; // && ev.l_charge[selLeptons[0]]!=ev.l_charge[selLeptons[1]]) continue;
        if(ev.l_id[selLeptons[0]]==ev.l_id[selLeptons[1]] && met.Pt() < 40) continue;
        doubleLep = true;
        allPlots["nj"+chTag]->Fill(lightJets.size(),wgt);
        allPlots["nbj"+chTag]->Fill(bJets.size(),wgt);
        allPlots["ndilp"+chTag]->Fill(selLeptons.size(),wgt);
        allPlots["dilp_pt"+chTag]->Fill(dilp4.Pt(),wgt);
        allPlots["dilp_m"+chTag]->Fill(dilp4.M(),wgt);
        allPlots["lp_pt"+chTag]->Fill(lp4.Pt(),wgt);
        allPlots["l2p_pt"+chTag]->Fill(l2p4.Pt(),wgt);
        if(debug) cout << "sorting jets" << endl;
	std::sort(lightJets.begin(), lightJets.end(), VecSort);
        std::sort(bJets.begin(), bJets.end(), VecSort);
        /*
        for (auto it : lightJets)
          allPlots["j_pt"+chTag]->Fill(it.Pt(),wgt);
        for (auto it : bJets)
          allPlots["bj_pt"+chTag]->Fill(it.Pt(),wgt);
        */
        allPlots["j_pt"+chTag]->Fill(lightJets[0].Pt(),wgt);
        allPlots["bj_pt"+chTag]->Fill(bJets[0].Pt(),wgt);
        TLorentzVector leadlp41,leadlp42;
        /*
        leadlp41.SetPtEtaPhiM(ev.l_pt[selLeptons[0]],ev.l_eta[selLeptons[0]],ev.l_phi[selLeptons[0]],ev.l_mass[selLeptons[0]]);
        leadlp42.SetPtEtaPhiM(ev.l_pt[selLeptons[1]],ev.l_eta[selLeptons[1]],ev.l_phi[selLeptons[1]],ev.l_mass[selLeptons[1]]);
        allPlots["leadL_pt"+chTag]->Fill(ev.l_pt[selLeptons[0]],wgt);
        allPlots["leadL_pt"+chTag]->Fill(ev.l_pt[selLeptons[1]],wgt);
        */
        allPlots["MET"+chTag]->Fill(ev.met_pt[0],wgt);
        allPlots["charge"+chTag]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
      }
      if(debug) cout << "simple plots DONE" << endl;


      if(!singleLep && !doubleLep) continue;
      
      // Charm resonance stuff:
      float maxcsv[2] = {-wgt, -wgt};
      int maxind[2] = {-1, -1};
      for(int k = 0; k < ev.nj; k++)
        if(ev.j_csv[k] >= maxcsv[1]) {
          if(ev.j_csv[k] >= maxcsv[0]) {
            maxcsv[0] = ev.j_csv[k];
            maxind[0] = k;
          }
          else {
            maxcsv[1] = ev.j_csv[k];
            maxind[1] = k;
          }
        }
      for(auto jetindex : maxind) {
        //int jetindex = maxind[0];
        if(jetindex < 0) continue;

        if(debug) cout << "l or ll" << endl;
        TLorentzVector p_track1, p_track2;
        const float gMassMu = 0.1057;
        int nstart = firstTrackIndex(jetindex);
        allPlots["nstart"+chTag]->Fill(nstart,wgt);
        allPlots["pfid"+chTag]->Fill(ev.pf_id[nstart],wgt);

        //J/Psi
        if(debug) cout << "starting J/Psi" << endl;
        for(int i = 0; i < ev.npf; i++) {
          if(ev.pf_j[i] != jetindex) continue;
          allPlots["lep_pt"+chTag]->Fill(ev.pf_pt[i],wgt);
          if(abs(ev.pf_id[i]) != 13 && abs(ev.pf_id[i]) != 11) continue;
          for(int j = 0; j < ev.npf; j++) {
            if(ev.pf_j[j] != ev.pf_j[i]) continue;
            /*
            if(abs(ev.pf_id[i]) !== abs(ev.pf_id[i])) continue;
            if(ev.pf_id[i]*ev.pf_id[j] > 0) continue; // e^+e^- or mu^+mu^-
            */
            if(ev.pf_id[i] != -ev.pf_id[j]) continue; // e^+e^- or mu^+mu^-
            allPlots["lep_pt"+chTag]->Fill(ev.pf_pt[j],wgt);

            float trackmass = gMassMu;
            if(abs(ev.pf_id[j]*ev.pf_id[j]) == 121) trackmass = 0.;
            p_track1.SetPtEtaPhiM(ev.pf_pt[i], ev.pf_eta[i], ev.pf_phi[i], trackmass);
            p_track2.SetPtEtaPhiM(ev.pf_pt[j], ev.pf_eta[j], ev.pf_phi[j], trackmass);

            float mass12 = (p_track1+p_track2).M();
            if(mass12>2.5 && mass12<3.5)
              allPlots["massJPsi"+chTag]->Fill(mass12,wgt);
          }
        }
        if(debug) cout << "J/Psi DONE" << endl;

        //D0 and D* 
        if(debug) cout << "Starting D0 and D*" << endl;
        //vector<TLorentzVector> jets;
        vector<pair<int,float>> jets;
        for(int i = 0; i < ev.npf; i++) {
          if(ev.pf_j[i] != jetindex) continue;
          jets.push_back(std::make_pair(i,ev.pf_pt[i]));
        }
        std::sort(jets.begin(), jets.end(), sortJetTuple);
        for(int i = 0; i < (int)jets.size(); i++)
          if(std::get<0>(jets.at(i)) == nstart) {
            nstart = i;
            break;
          }
        nstart = firstTrackIndex(jetindex,&jets);
        if(jets.size() < 3) continue;
        if((jets.size() - nstart) < 3) continue;
        for(int i = nstart; i < nstart+3; i++)
        //for(int i = 0; i < 3; i++)
          for(int j = i+1; j < 3; j++) {
            int tk1 = get<0>(jets.at(i));
            int tk2 = get<0>(jets.at(j));
            /*
            int tk1 = i;
            int tk2 = j;
            */
            if(ev.pf_j[tk1] != jetindex) continue;
            if(ev.pf_j[tk2] != jetindex) continue;

            //opposite sign
            if(ev.pf_id[tk1]*ev.pf_id[tk2] != -211*211) continue;

            const float gMassK  = 0.4937;
            const float gMassPi = 0.1396;
          
            p_track1.SetPtEtaPhiM(ev.pf_pt[tk1], ev.pf_eta[tk1], ev.pf_phi[tk1], gMassPi);
            p_track2.SetPtEtaPhiM(ev.pf_pt[tk2], ev.pf_eta[tk2], ev.pf_phi[tk2], gMassK);
            float mass12 = (p_track1+p_track2).M();

            //if (mass12>1.65 && mass12<2.0)
            if (mass12>1.7 && mass12<2.0)
              allPlots["massD0"+chTag]->Fill(mass12,wgt);

            //looking for lepton
            if(debug) cout << "third lepton" << endl;
            //for(int tk3 = 0; tk3 < ev.npf; tk3++) {
            for(int k = 0; k < (int)jets.size(); k++) {
              int tk3 = get<0>(jets.at(k));
              if(ev.pf_j[tk3] != jetindex) continue;
              if(tk3 == tk1) continue;
              if(tk3 == tk2) continue;
              if(debug) cout << "third lepton possible" << endl;
            
              if(abs(ev.pf_id[tk3]) != 13 && abs(ev.pf_id[tk3]) != 11) continue;
              if(debug) cout << "third lepton found" << endl;

              if(ev.pf_id[tk2]/abs(ev.pf_id[tk2]) == -ev.pf_id[tk3]/abs(ev.pf_id[tk3])) {
                //Kaon and lepton have same charge
                //correct mass assumption
                if(debug) cout << "correct mass assumption" << endl;
                allPlots["masslep"+chTag]->Fill(mass12,wgt);

                if(abs(ev.pf_id[tk3]) == 13)
                  allPlots["massmu"+chTag]->Fill(mass12,wgt);
                if(abs(ev.pf_id[tk3]) == 11)
                  allPlots["masse"+chTag]->Fill(mass12,wgt);

              }
            }
            //looking for pion
            if(debug) cout << "D*->pi+D0" << endl;
            //for(int tk3 = 0; tk3 < ev.npf; tk3++) {
            for(int tk3 = 0; tk3 < (int)jets.size(); tk3++) {
              if(ev.pf_j[tk3] != jetindex) continue;
              if(tk3 == tk1) continue;
              if(tk3 == tk2) continue;

              if(abs(ev.pf_id[tk3]) != 211) continue;
              if(debug) cout << "Pion found" << endl;

              TLorentzVector p_track3, p_cand;
              p_track3.SetPtEtaPhiM(ev.pf_pt[tk3], ev.pf_eta[tk3], ev.pf_phi[tk3], gMassPi);
              allPlots["pi_pt"+chTag]->Fill(p_track3.Pt(),wgt);
              if( ev.pf_id[tk2]/abs(ev.pf_id[tk2]) == -ev.pf_id[tk3]/abs(ev.pf_id[tk3]) ) {
                // Kaon and pion have opposite charges
                // I.e. correct mass assumption
                if(debug) cout << "correct mass assumption" << endl;
                
                p_cand = p_track1+p_track2+p_track3;
                allPlots["massDs"+chTag]->Fill(p_cand.M(), wgt);

                if(abs(mass12-1.864) < 0.10) { // mass window cut
                  TLorentzVector p_jet;
                  p_jet.SetPtEtaPhiM(ev.j_pt[jetindex], ev.j_eta[jetindex], ev.j_phi[jetindex], 0.);

                  //float hardpt = std::max(ev.pf_pt[tk3], std::max(ev.pf_pt[tk1], ev.pf_pt[tk2]));
                  //float softpt = std::min(ev.pf_pt[tk3], std::min(ev.pf_pt[tk1], ev.pf_pt[tk2]));
                  float deltam = p_cand.M() - mass12;

                  allPlots["massDsmD0loose"+chTag]->Fill(deltam, wgt);
                  if(abs(mass12-1.864) < 0.05) { // tighter mass window cut
                      //FillCharmTree(413,  jetindex, tk1, gMassPi, tk2, gMassK, tk3, gMassPi);
                      //FillCharmTree(-413, jetindex, deltam, p_cand, p_jet, hardpt, softpt);
                      allPlots["massDsmD0"+chTag]->Fill(deltam, wgt);
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
  if(flavourSplitting!=NOFLAVOURSPLITTING) selPrefix=Form("%d_",flavourSplitting);
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  fOut->cd();
  if(debug) cout << "writing histograms" << endl;
/*
  for(auto it : lfsVec) {
    it.Remove(0,1);
    fOut->mkdir(it);
  }
*/
  for (auto& it : allPlots)  { 
    if(debug) cout << it.second->GetName() << endl;
    if(debug) cout << it.second->GetEntries() << endl;
/*
    TString dir = it.first;
    dir.Remove(0,dir.Last('_')+1);
*/
    //fOut->cd( dir );
    it.second->SetDirectory(fOut); it.second->Write(); 
    fOut->cd();
  }
  if(debug) cout << "writing histograms DONE" << endl;
  if(debug) cout << "closing ROOT file" << endl;
  fOut->Close();
}

