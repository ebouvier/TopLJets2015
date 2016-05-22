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

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"


#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;


//
void RunTopWidth(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
                 Bool_t debug=false)
{

  bool isTTbar( filename.Contains("_TTJets") );

  //prepare output
  TopWidthEvent_t twev;
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("twev","twev");
  createTopWidthEventTree(outT,twev);
  outT->SetDirectory(fOut);

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
  bool requireEtriggerOnly(false);
  twev.isData = ev.isData;
  if(ev.isData && filename.Contains("SingleElectron")) requireEtriggerOnly=true;

  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
    {
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
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  if(!ev.isData)
    {
      BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );

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

  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);

  std::vector<TString> lfsVec = { "E", "EE", "EM", "MM", "M" }; 
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   
    { 
      TString tag(lfsVec[ilfs]);
      allPlots["nvtx_"+tag]  = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events",30,0,30);
      allPlots["njets_"+tag]  = new TH1F("njets_"+tag,";Jet multiplicity;Events",5,0,5);
      for(int i=0; i<2; i++)
	{
	  TString pf(Form("l%d",i));	  
	  allPlots[pf+"pt_"+tag]  = new TH1F(pf+"pt_"+tag,";Lepton p_{t} [GeV];Events",50,0,250);
	  allPlots[pf+"eta_"+tag]  = new TH1F(pf+"eta_"+tag,";Lepton pseudo-rapidity;Events",50,0,2.5);
	}
      for(int i=0; i<6; i++)
	{
	  TString pf(Form("j%d",i));
	  allPlots[pf+"pt_"+tag]  = new TH1F(pf+"pt_"+tag,";Jet transverse momentum [GeV];Events",50,0,250);
	  allPlots[pf+"eta_"+tag]  = new TH1F(pf+"eta_"+tag,";Jet pseudo-rapidity;Events",50,0,4.7);
	}
    }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      resetTopWidthEvent(twev);
      if(iev%10000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //account for pu weights and effect on normalization
      float puWeight(1.0);
      if(!ev.isData) 
	{
	  puWeight=puWgtGr[0]->Eval(ev.putrue);  
	  allPlots["puwgtctr"]->Fill(0.,1.0);
	  allPlots["puwgtctr"]->Fill(1.,puWeight);
	}

      //select 1 good lepton
      std::vector<int> tightLeptons,looseLeptons;
      for(int il=0; il<ev.nl; il++)
	{
	  bool passTightKin(ev.l_pt[il]>30 && fabs(ev.l_eta[il])<2.1);
	  bool passLooseKin(ev.l_pt[il]>10 && fabs(ev.l_eta[il])<2.5);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
	  float relIso(ev.l_relIso[il]);
	  bool passTightIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 );
	  bool passLooseIso(  ev.l_id[il]==13 ? relIso<0.25 : true);
	  if(passTightKin && passTightId && passTightIso) tightLeptons.push_back(il);
	  else if(passLooseKin && passLooseIso)           looseLeptons.push_back(il);
	}
      
      //check if triggers have fired
      bool hasMuTrigger((ev.muTrigger & 0x3)!=0);
      bool hasEleTrigger((ev.elTrigger & 0x1)!=0);
      twev.muTrigger = ev.muTrigger;
      twev.elTrigger = ev.elTrigger;

      //decide the channel
      TString chTag("");
      std::vector<int> selLeptons;
      if(tightLeptons.size()==1 && looseLeptons.size()==0)
	{
	  selLeptons.push_back( tightLeptons[0] );
	}
      if(tightLeptons.size()>=2)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(tightLeptons[1]);
	}
      if(tightLeptons.size()==1 && looseLeptons.size()>=1)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(looseLeptons[0]);	  
	}
      if(selLeptons.size()==0) continue;
      if(selLeptons.size()==1)
	{
	  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="E";
	  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger)  chTag="M";
	}
      if(selLeptons.size()==2)
	{
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*11 && hasEleTrigger) chTag="EE";
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==13*13 && hasMuTrigger) chTag="MM";
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*13)
	    {
	      if(tightLeptons.size()>=2 && (hasEleTrigger || hasMuTrigger)) chTag="EM";
	      if(tightLeptons.size()==1)
		{
		  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="EM";
		  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger) chTag="EM";
		}
	    }
	  if(hasMuTrigger && requireEtriggerOnly) chTag="";
	}
      if(chTag=="") continue;

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
      std::vector<int> genBJetsFlav,genOtherjetsFlav,genBJetsHadFlav,genOtherjetsHadFlav;
      std::vector<TLorentzVector> bJets,otherjets, genBJets, genOtherjets;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  //cross clean with leptons
	  bool overlapsWithLepton(false);
	  for(size_t il=0; il<leptons.size(); il++)
	    {
	      if(jp4.DeltaR(leptons[il])>0.4) continue;
	      overlapsWithLepton=true;
	    }
	  if(overlapsWithLepton) continue;

	  //smear jet energy resolution for MC
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }

	  //b-tag
	  float csv = ev.j_csv[k];	  
          twev.j_csv[k]=ev.j_csv[k];
	  bool isBTagged(csv>0.800);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 	
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		}
	      else
		{
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
		}
	      
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	    }

	  if(jp4.Pt()<30) continue;

	  Int_t hadFlav=ev.j_hadflav[k];
	  Int_t flav=ev.j_flav[k];
	  TLorentzVector gjp4(0,0,0,0);
	  if(ev.j_g[k]>0)
	    {
	      int gidx=ev.j_g[k];
	      gjp4.SetPtEtaPhiM( ev.g_pt[gidx], ev.g_eta[gidx], ev.g_phi[gidx], ev.g_m[gidx] );
	    }

	  //save jets 
	  if(isBTagged && fabs(jp4.Eta()) < 2.4 && bJets.size()<2) 
	    { bJets.push_back(jp4);     genBJets.push_back(gjp4);     genBJetsFlav.push_back(flav);     genBJetsHadFlav.push_back(hadFlav);    }
	  else 
	    { otherjets.push_back(jp4); genOtherjets.push_back(gjp4); genOtherjetsFlav.push_back(flav); genOtherjetsHadFlav.push_back(hadFlav); }
	}

      //2 b-jets are required 
      //notice more can be b-tagged but they are put into the other jets category
      if(bJets.size()!=2) continue;

      //event weight
      float wgt(1.0);
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
          wgt=1.0;
	}
      
      //nominal selection control histograms
      twev.nvtx = ev.nvtx;
      allPlots["nvtx_"+chTag]->Fill(ev.nvtx,wgt);
      allPlots["njets_"+chTag]->Fill(otherjets.size(),wgt);

      twev.nl=leptons.size();
      for(int il=0; il<(int)leptons.size(); il++)
	{
	  TString pf(Form("l%d",il));
	  allPlots[pf+"pt_"+chTag]->Fill(leptons[il].Pt(),wgt);
	  allPlots[pf+"eta_"+chTag]->Fill(fabs(leptons[il].Eta()),wgt);
	  twev.l_pt[il]=leptons[il].Pt();
	  twev.l_eta[il]=leptons[il].Eta();
	  twev.l_phi[il]=leptons[il].Phi();
	  twev.l_m[il]=leptons[il].M();
	  twev.l_id[il]=ev.l_id[ selLeptons[il] ];
	  twev.l_pid[il]=ev.l_pid[ selLeptons[il] ];
	  twev.l_charge[il]=ev.l_charge[ selLeptons[il] ];
	  twev.l_relIso[il]=ev.l_relIso[ selLeptons[il] ];
	  for(Int_t ig=0; ig<ev.ng; ig++)
	    {
	      if(abs(ev.g_id[ig])!=ev.l_id[ selLeptons[il] ]) continue;
	      TLorentzVector glp4;
	      glp4.SetPtEtaPhiM( ev.g_pt[ig], ev.g_eta[ig], ev.g_phi[ig], ev.g_m[ig]);
	      if(glp4.DeltaR( leptons[il] ) > 0.3) continue;
	      twev.gl_id[il]=ev.g_id[ig];
	      twev.gl_pt[il]=ev.g_pt[ig];
	      twev.gl_eta[il]=ev.g_eta[ig];
	      twev.gl_phi[il]=ev.g_phi[ig];
	      twev.gl_m[il]=ev.g_m[ig];
	    }
	}
      
      twev.nj=bJets.size();
      for(int ij=0; ij<(int)bJets.size(); ij++)
	{
	  TString pf(Form("j%d",ij));
	  allPlots[pf+"pt_"+chTag]->Fill(bJets[ij].Pt(),wgt);
	  allPlots[pf+"eta_"+chTag]->Fill(fabs(bJets[ij].Eta()),wgt);
	  twev.j_pt[ij]=bJets[ij].Pt();
	  twev.j_eta[ij]=bJets[ij].Eta();
	  twev.j_phi[ij]=bJets[ij].Phi();
	  twev.j_m[ij]=bJets[ij].M();
	  twev.j_flav[ij]=genBJetsFlav[ij];
	  twev.j_hadflav[ij]=genBJetsHadFlav[ij];
	  twev.gj_pt[ij]=genBJets[ij].Pt();
	  twev.gj_eta[ij]=genBJets[ij].Eta();
	  twev.gj_phi[ij]=genBJets[ij].Phi();
	  twev.gj_m[ij]=genBJets[ij].M();	  
	}
      
      twev.nj+=otherjets.size();
      for(int ij=0; ij<(int)otherjets.size(); ij++)
	{
	  if(ij+2<6)
	    {
	      TString pf(Form("j%d",ij+2));
	      allPlots[pf+"pt_"+chTag]->Fill(otherjets[ij].Pt(),wgt);
	      allPlots[pf+"eta_"+chTag]->Fill(fabs(otherjets[ij].Eta()),wgt);
	    }
	  twev.j_pt[ij+2]=otherjets[ij].Pt();
	  twev.j_eta[ij+2]=otherjets[ij].Eta();
	  twev.j_phi[ij+2]=otherjets[ij].Phi();
	  twev.j_m[ij+2]=otherjets[ij].M();
	  twev.j_flav[ij+2]=genOtherjetsFlav[ij];
	  twev.j_hadflav[ij+2]=genOtherjetsHadFlav[ij];
	  twev.gj_pt[ij+2]=genOtherjets[ij].Pt();
	  twev.gj_eta[ij+2]=genOtherjets[ij].Eta();
	  twev.gj_phi[ij+2]=genOtherjets[ij].Phi();
	  twev.gj_m[ij+2]=genOtherjets[ij].M();	  
	}

      twev.cat=11;
      if(chTag=="M") twev.cat=13;
      if(chTag=="MM") twev.cat=13*13;
      if(chTag=="EM") twev.cat=11*13;
      if(chTag=="EE") twev.cat=11*11;
      twev.nw=1;
      twev.weight[0]=wgt;
      twev.nt=0;
      for(int i = 0; i < 10; i++) {
        twev.met_pt[i]=ev.met_pt[i];
        twev.met_phi[i]=ev.met_phi[i];
      }

      if(ev.ngtop>0)
	{
	  for(int i=0; i<ev.ngtop; i++)
	    {
	      if(abs(ev.gtop_id[i])!=6) continue;
	      twev.t_pt[twev.nt]=ev.gtop_pt[i];
	      twev.t_eta[twev.nt]=ev.gtop_eta[i];
	      twev.t_phi[twev.nt]=ev.gtop_phi[i];
	      twev.t_m[twev.nt]=ev.gtop_m[i];
	      twev.t_id[twev.nt]=ev.gtop_id[i];
	      twev.nt++;
	      if(twev.nt>4) break;
	    }
	}
      for(int i = 0; i < ev.npf; i++) {
        twev.npf = ev.npf;
        twev.pf_j[i] = ev.pf_j[i];
        twev.pf_id[i] = ev.pf_id[i];
        twev.pf_pt[i] = ev.pf_pt[i];
        twev.pf_eta[i] = ev.pf_eta[i];
        twev.pf_phi[i] = ev.pf_phi[i];
        twev.ngpf = ev.ngpf;
        twev.gpf_id[i] = ev.gpf_id[i];
        twev.gpf_pt[i] = ev.gpf_pt[i];
        twev.gpf_eta[i] = ev.gpf_eta[i];
        twev.gpf_phi[i] = ev.gpf_phi[i];
      }
      outT->Fill();
    }
  
  //close input file
  f->Close();

  //save histos to file  
  fOut->cd();
  outT->Write();
  puTrue->Write();
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}



//
void createTopWidthEventTree(TTree *t,TopWidthEvent_t &twev)
{
  //event category
  t->Branch("cat", &twev.cat,"cat/I");
  t->Branch("isData", &twev.isData,"isData/I");

  //event weights
  t->Branch("nw",  &twev.nw, "nw/I");
  t->Branch("weight",  twev.weight, "weight[nw]/F");

  //met
  t->Branch("met_pt",  &twev.met_pt, "met_pt/F");
  t->Branch("met_phi",  &twev.met_phi, "met_phi/F");

  //leptons
  t->Branch("nl",  &twev.nl, "nl/I");
  t->Branch("l_pt",  twev.l_pt ,  "l_pt[nl]/F");
  t->Branch("l_eta", twev.l_eta , "l_eta[nl]/F");
  t->Branch("l_phi", twev.l_phi , "l_phi[nl]/F");
  t->Branch("l_m",   twev.l_m ,   "l_m[nl]/F");
  t->Branch("l_id",   twev.l_id ,   "l_id[nl]/I");
  t->Branch("l_pid",   twev.l_pid ,   "l_pid[nl]/I");
  t->Branch("l_charge",   twev.l_charge ,   "l_charge[nl]/I");
  t->Branch("l_relIso",   twev.l_relIso ,   "l_relIso[nl]/I");
  t->Branch("gl_pt",  twev.gl_pt ,  "gl_pt[nl]/F");
  t->Branch("gl_eta", twev.gl_eta , "gl_eta[nl]/F");
  t->Branch("gl_phi", twev.gl_phi , "gl_phi[nl]/F");
  t->Branch("gl_m",   twev.gl_m ,   "gl_m[nl]/F");
  t->Branch("gl_id",  twev.gl_id ,  "gl_id[nl]/I");

  //jets
  t->Branch("nj",  &twev.nj, "nj/I");
  t->Branch("j_pt",  twev.j_pt ,  "j_pt[nj]/F");
  t->Branch("j_eta", twev.j_eta , "j_eta[nj]/F");
  t->Branch("j_phi", twev.j_phi , "j_phi[nj]/F");
  t->Branch("j_m",   twev.j_m ,   "j_m[nj]/F");
  t->Branch("gj_pt",  twev.gj_pt ,  "gj_pt[nj]/F");
  t->Branch("gj_eta", twev.gj_eta , "gj_eta[nj]/F");
  t->Branch("gj_phi", twev.gj_phi , "gj_phi[nj]/F");
  t->Branch("gj_m",   twev.gj_m ,   "gj_m[nj]/F");
  t->Branch("gj_flav",  twev.gj_flav ,  "gj_flav[nj]/I");
  t->Branch("gj_hadflav",  twev.gj_hadflav ,  "gj_hadflav[nj]/I");
  t->Branch("j_csv",   twev.j_csv ,   "j_csv[nj]/F");
  t->Branch("j_hadflav",  twev.j_hadflav ,  "j_hadflav[nj]/I");
  t->Branch("j_flav",  twev.j_flav ,  "j_flav[nj]/I");

  //mc truth
  t->Branch("nt",  &twev.nt, "nt/I");
  t->Branch("t_pt",  twev.t_pt ,  "t_pt[nt]/F");
  t->Branch("t_eta", twev.t_eta , "t_eta[nt]/F");
  t->Branch("t_phi", twev.t_phi , "t_phi[nt]/F");
  t->Branch("t_m",   twev.t_m ,   "t_m[nt]/F");
  t->Branch("t_id",  twev.t_id ,  "t_id[nt]/I");

  //pf
  t->Branch("npf",        &twev.npf,         "npf/I");
  t->Branch("pf_j",        twev.pf_j,        "pf_j[npf]/I");
  t->Branch("pf_id",       twev.pf_id,       "pf_id[npf]/I");  
  t->Branch("pf_pt",       twev.pf_pt,       "pf_pt[npf]/F");
  t->Branch("pf_eta",      twev.pf_eta,      "pf_eta[npf]/F");
  t->Branch("pf_phi",      twev.pf_phi,      "pf_phi[npf]/F");
  t->Branch("ngpf",        &twev.ngpf,         "ngpf/I");
  t->Branch("gpf_id",       twev.gpf_id,       "gpf_id[npf]/I");  
  t->Branch("gpf_pt",       twev.gpf_pt,       "gpf_pt[npf]/F");
  t->Branch("gpf_eta",      twev.gpf_eta,      "gpf_eta[npf]/F");
  t->Branch("gpf_phi",      twev.gpf_phi,      "gpf_phi[npf]/F");

  //reco level event
  t->Branch("nvtx",      &twev.nvtx,      "nvtx/I");
  t->Branch("muTrigger",        &twev.muTrigger,        "muTrigger/I");
  t->Branch("elTrigger",        &twev.elTrigger,        "elTrigger/I");
}

//
void resetTopWidthEvent(TopWidthEvent_t &twev)
{
  twev.isData=false;
  twev.cat=0;   twev.nw=0;   twev.nl=0;   twev.nj=0;   twev.nt=0;

  for(int i=0; i<10; i++) { twev.weight[i]=0; twev.met_pt[i]=0; twev.met_phi[i]=0; }
  for(int i=0; i<2; i++) { twev.l_pt[i]=0;   twev.l_eta[i]=0;   twev.l_phi[i]=0;   twev.l_m[i]=0; twev.l_id[i]=0; twev.l_pid[i]=0; twev.l_relIso[i]=0; twev.gl_pt[i]=0;   twev.gl_eta[i]=0;   twev.gl_phi[i]=0;   twev.gl_m[i]=0; twev.gl_id[i]=0; twev.l_charge[i]=0; }
  for(int i=0; i<50; i++) { twev.j_pt[i]=0;   twev.j_eta[i]=0;   twev.j_phi[i]=0;   twev.j_m[i]=0; twev.j_hadflav[i]=0; twev.j_flav[i]=0; twev.j_csv[i]=0; }
  for(int i=0; i<50; i++) { twev.gj_pt[i]=0;   twev.gj_eta[i]=0;   twev.gj_phi[i]=0;   twev.gj_m[i]=0 ; twev.gj_flav[i]=0; }
  for(int i=0; i<4; i++) { twev.t_pt[i]=0;   twev.t_eta[i]=0;   twev.t_phi[i]=0;   twev.t_m[i]=0; twev.t_id[i]=0; }
  for(int i =0; i<5000; i++) { twev.npf=0; twev.pf_j[i]=0; twev.pf_id[i]=0; twev.pf_pt[i]=0; twev.pf_eta[i]=0; twev.pf_phi[i]=0; twev.ngpf=0; twev.gpf_pt[i]=0; twev.gpf_eta[i]=0; twev.gpf_phi[i]=0; twev.gpf_id[i]=0; }
}
