#ifndef _topwidth_h_
#define _topwidth_h_

struct TopWidthEvent_t
{
  Bool_t isData;
  Int_t cat,nw,nl,nj,nt;
  Float_t weight[10];
  Float_t l_pt[2],l_eta[2],l_phi[2],l_m[2],l_charge[2],l_relIso[2];
  Int_t l_id[2],l_pid[2];
  Float_t gl_pt[2],gl_eta[2],gl_phi[2],gl_m[2];
  Int_t gl_id[2];
  Float_t j_pt[50],j_eta[50],j_phi[50],j_m[50],j_csv[50],j_flav[50],j_hadflav[50];
  Float_t gj_pt[50],gj_eta[50],gj_phi[50],gj_m[50];
  Int_t gj_flav[50],gj_hadflav[50];
  Float_t t_pt[4],t_eta[4],t_phi[4],t_m[4];
  Int_t t_id[4];
  Float_t met_pt[10],met_phi[10];
  Int_t npf,ngpf,pf_j[5000];
  Int_t pf_id[5000];
  Float_t pf_pt[5000],pf_eta[5000],pf_phi[5000];
  Int_t gpf_id[5000];
  Float_t gpf_pt[5000],gpf_eta[5000],gpf_phi[5000];

  //reco level event
  Int_t muTrigger, elTrigger, nvtx;
};

void createTopWidthEventTree(TTree *t,TopWidthEvent_t &twev);
void resetTopWidthEvent(TopWidthEvent_t &twev);
#endif
