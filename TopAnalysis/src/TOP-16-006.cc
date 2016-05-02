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
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "TopLJets2015/TopAnalysis/interface/OtherFunctions.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

using namespace std;

//
Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}


//
void RunTop16006(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts)
{

  //bool isTTbar( filename.Contains("_TTJets") );
  
  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  if(ev.isData) runSysts=false;
  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;
  
  //PILEUP WEIGHTING

  //LEPTON EFFICIENCIES

  //B-TAG CALIBRATION
  TString btagUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  TString btagEffExpUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff;
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
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
    }

  //JET ENERGY SCALE: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Summer15_uncertainties
  
  //LIST OF SYSTEMATICS
  
  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  allPlots["lp_pt"] = new TH1F("lp_pt",";Lepton P_{T} [GeV];Events", 20, 0,300);
  allPlots["dilp_pt"] = new TH1F("dilp_pt",";Lepton P_{T} [GeV];Events", 20, 0,300);
  allPlots["dilp_m"] = new TH1F("dilp_m",";M_{ll} [GeV];Events", 20, 0,300);
  allPlots["j_pt"] = new TH1F("j_pt",";Jet P_{T} [GeV];Events", 20, 0,300);
  allPlots["bj_pt"] = new TH1F("bj_pt",";b Jet P_{T} [GeV];Events", 20, 0,300);
  allPlots["nlp"]     = new TH1F("nlp",";N_{l};Events" ,3,0.,3.);
  allPlots["ndilp"]     = new TH1F("ndilp",";N_{ll};Events" ,3,0.,3.);
  allPlots["nj"]     = new TH1F("nj",";N_{jets};Events" ,3,0.,5.);
  allPlots["nbj"]     = new TH1F("nbj",";N_{b-jets};Events" ,3,0.,5.);
  allPlots["nstart"]     = new TH1F("nstart",";N_{start};Events" ,3,0.,5.);
  allPlots["pfid"]     = new TH1F("pfid",";PFID;Events" ,440,-220.,220.);
  allPlots["massD0"]     = new TH1F("massD0",";M_{jj};Events" ,20,0.,10.);
  allPlots["masslep"]     = new TH1F("masslep",";M_{K#pi};Events" ,20,0.,10.);
  allPlots["massmu"]     = new TH1F("massmu",";M_{K#pi};Events" ,20,0.,10.);
  allPlots["masse"]     = new TH1F("masse",";M_{K#pi};Events" ,20,0.,10.);
  allPlots["massDsmD0loose"]     = new TH1F("massDsmD0loose",";M_{K#pi};Events" ,20,0.,10.);
  allPlots["massDsmD0"]     = new TH1F("massDsmD0",";M_{K#pi};Events" ,20,0.,10.);


  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //select 1 good lepton
      //cout << "entering lepton selection" << endl;
      std::vector<int> tightLeptonsIso, tightLeptonsNonIso, vetoLeptons;
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
	  bool passSIP3d(ev.l_ip3dsig[il]<4);
	  if(channelSelection==21) passSIP3d=true;

	  if(passTightKin && passTightId && passSIP3d)
	    {
	      if(passIso)         tightLeptonsIso.push_back(il);
	      else if(passNonIso) tightLeptonsNonIso.push_back(il);
	    }
	  else if(passVetoKin && passVetoIso) vetoLeptons.push_back(il);
	}
      //cout << "lepton selection DONE" << endl;

      //one good lepton either isolated or in the non-isolated sideband or a Z candidate
      Int_t lepIdx=-1;
      Bool_t isZ(false);//,isZPassingSIP3d(false);
      TLorentzVector l1p4,l2p4,dilp4;
      if(tightLeptonsIso.size()==1)                                       lepIdx=tightLeptonsIso[0];
      else if (tightLeptonsIso.size()==0 && tightLeptonsNonIso.size()==1) lepIdx=tightLeptonsNonIso[0];
      else if(tightLeptonsIso.size()==2)
	{	  
          //cout << "di-lepton" << endl;
	  l1p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[0]],ev.l_eta[tightLeptonsIso[0]],ev.l_phi[tightLeptonsIso[0]],ev.l_mass[tightLeptonsIso[0]]);
	  l2p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[1]],ev.l_eta[tightLeptonsIso[1]],ev.l_phi[tightLeptonsIso[1]],ev.l_mass[tightLeptonsIso[1]]);
	  dilp4=l1p4+l2p4;
	  if(ev.l_id[tightLeptonsIso[0]]==ev.l_id[tightLeptonsIso[1]]          && 
	     ev.l_charge[tightLeptonsIso[0]]*ev.l_charge[tightLeptonsIso[1]]<0 && 
	     fabs(dilp4.M()-91)<10 && 
	     dilp4.Pt()>30)
	    { 
	      isZ=true; 
	      //isZPassingSIP3d=(ev.l_ip3dsig[0]<4 && ev.l_ip3dsig[1]<4);
	    }
	  lepIdx=tightLeptonsIso[0];
          //cout << "di-lepton DONE" << endl;
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
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);

      //select jets
      Float_t htsum(0);
      TLorentzVector jetDiff(0,0,0,0);
      std::vector<TLorentzVector> bJets,lightJets;
      TLorentzVector visSystem(isZ ? dilp4 : lp4);
      int nbjets(0),ncjets(0),nljets(0),leadingJetIdx(-1);
      std::vector<int> resolvedJetIdx;
      std::vector<TLorentzVector> resolvedJetP4;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  //jp4=updateJES(jp4,ev.j_rawsf[k],ev.j_area[k],ev.rho,ev.nvtx,jetCorr);

	  //cross clean with respect to leptons 
	  if(jp4.DeltaR(lp4)<0.5) continue;
	  if(isZ && jp4.DeltaR(l2p4)<0.5)continue;

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
	  
	  if(leadingJetIdx<0) leadingJetIdx=k;
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

      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      met+=jetDiff;
      met.SetPz(0.); met.SetE(met.Pt());
      //float mt( computeMT(isZ ? dilp4: lp4,met) );

      // Charm resonance stuff:
      float maxcsv(-1.);
      int maxind(-1);
      for(int k = 0; k < ev.nj; k++)
        if(ev.j_csv[k] >= maxcsv) {
          maxcsv = ev.j_csv[k];
          maxind = k;
        }
      int jetindex = maxind;

      //simple fill
      if(tightLeptonsIso.size() == 1 && bJets.size()+lightJets.size() >= 4) {
        allPlots["nj"]->Fill(lightJets.size(),1);
        allPlots["nbj"]->Fill(bJets.size(),1);
        allPlots["nlp"]->Fill(tightLeptonsIso.size(),1);
        allPlots["lp_pt"]->Fill(lp4.Pt(),1);
        for (auto it : lightJets)
          allPlots["j_pt"]->Fill(it.Pt(),1);
        for (auto it : bJets)
          allPlots["bj_pt"]->Fill(it.Pt(),1);
        int nstart = firstTrackIndex(maxind);
        allPlots["nstart"]->Fill(nstart,1);
        allPlots["pfid"]->Fill(ev.pf_id[nstart],1);
        TLorentzVector p_track1, p_track2;
        for(int i = nstart; i < nstart+3; i++)
          for(int j = i+1; j < 3; j++) {
            int tk1 = i;
            int tk2 = j;

            //opposite sign
            if(ev.pf_id[tk1]*ev.pf_id[tk2] != -211*211) continue;

            const float gMassK  = 0.4937;
            const float gMassPi = 0.1396;
            //const float gMassMu = 0.1057;
          
            p_track1.SetPtEtaPhiM(ev.pf_pt[tk1], ev.pf_eta[tk1], ev.pf_phi[tk1], gMassPi);
            p_track2.SetPtEtaPhiM(ev.pf_pt[tk2], ev.pf_eta[tk2], ev.pf_phi[tk2], gMassK);
            float mass12 = (p_track1+p_track2).M();

            if (mass12>1.65 && mass12<2.0)
              allPlots["massD0"]->Fill(mass12,1);

            //looking for lepton
            for(int tk3 = 0; tk3 < ev.npf; tk3++) {
              if(ev.pf_id[tk3] != jetindex) continue;
              if(tk3 == tk1) continue;
              if(tk3 == tk2) continue;
            
              if(abs(ev.pf_id[tk3]) != 13 && abs(ev.pf_id[tk3]) != 11) continue;

              if(ev.pf_id[tk2]/abs(ev.pf_id[tk2]) == -ev.pf_id[tk3]/abs(ev.pf_id[tk3])) {
                //Kaon and lepton have same charge
                //correct mass assumption
                allPlots["masslep"]->Fill(mass12,1);

                if(abs(ev.pf_id[tk3]) == 13)
                  allPlots["massmu"]->Fill(mass12,1);
                if(abs(ev.pf_id[tk3]) == 11)
                  allPlots["masse"]->Fill(mass12,1);

              }
            }
            //looking for pion
            for(int tk3 = 0; tk3 < ev.npf; tk3++) {
              if(ev.pf_id[tk3] != jetindex) continue;
              if(tk3 == tk1) continue;
              if(tk3 == tk2) continue;

              if(abs(ev.pf_id[tk3]) != 211) continue;

              if( ev.pf_id[tk2]/abs(ev.pf_id[tk2]) == -ev.pf_id[tk3]/abs(ev.pf_id[tk3]) ) {
                // Kaon and pion have opposite charges
                // I.e. correct mass assumption
                
                if(abs(mass12-1.864) < 0.10) { // mass window cut
                  TLorentzVector p_track3;
                  p_track3.SetPtEtaPhiM(ev.pf_pt[tk3], ev.pf_eta[tk3], ev.pf_phi[tk3], gMassPi);

                  TLorentzVector p_cand, p_jet;
                  p_cand = p_track1+p_track2+p_track3;
                  p_jet.SetPtEtaPhiM(ev.j_pt[jetindex], ev.j_eta[jetindex], ev.j_phi[jetindex], 0.);

                  //float hardpt = std::max(ev.pf_pt[tk3], std::max(ev.pf_pt[tk1], ev.pf_pt[tk2]));
                  //float softpt = std::min(ev.pf_pt[tk3], std::min(ev.pf_pt[tk1], ev.pf_pt[tk2]));
                  float deltam = (p_track1+p_track2+p_track3).M() - mass12;

                  allPlots["massDsmD0loose"]->Fill(deltam, 1);
                  if(abs(mass12-1.864) < 0.05) { // tighter mass window cut
                      //FillCharmTree(413,  jetindex, tk1, gMassPi, tk2, gMassK, tk3, gMassPi);
                      //FillCharmTree(-413, jetindex, deltam, p_cand, p_jet, hardpt, softpt);
                      allPlots["massDsmD0"]->Fill(deltam, 1);
                  }
                }
              }
            }
        }
      }
      else if(tightLeptonsIso.size() == 2 && bJets.size()+lightJets.size() >= 2) {
        if(isZ) continue;
        if(ev.l_id[tightLeptonsIso[0]]==ev.l_id[tightLeptonsIso[1]] && met.Pt() < 40) continue;
        allPlots["nj"]->Fill(lightJets.size(),1);
        allPlots["nbj"]->Fill(bJets.size(),1);
        allPlots["ndilp"]->Fill(tightLeptonsIso.size(),1);
        allPlots["dilp_pt"]->Fill(dilp4.Pt(),1);
        allPlots["dilp_m"]->Fill(dilp4.M(),1);
        for (auto it : lightJets)
          allPlots["j_pt"]->Fill(it.Pt(),1);
        for (auto it : bJets)
          allPlots["bj_pt"]->Fill(it.Pt(),1);
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
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

