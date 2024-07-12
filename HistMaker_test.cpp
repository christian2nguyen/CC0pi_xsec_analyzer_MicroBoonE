/*Author: Christian
 *
 * Usage: making hist to plot and look at for CC-0pi
 * THis testing script is to be used to check plotting methods
// */
//
#include <fstream> //TODO: Which library is necessary?
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TPad.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include <algorithm>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TTree.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TChain.h"
#include "TLegendEntry.h"
#include <TGraphErrors.h>
#include "TLine.h"
#include "TAttLine.h"
#include "TF1.h"

#include "includes/MnvColors.hh"
#include "includes/NamedCategory.hh"
#include "includes/EventCategory.hh"
#include "includes/HistFolio_slim.hh"
#include "includes/HistFolio_slim.cpp"
//#include "HistUtils.hh"
#include "WeightHandler.hh"
#include "includes/UBTH2Poly.h"
#include "includes/HistUtils_cc0pi.hh"
#include "UniverseMaker.hh"
//#include "HistMaker.hh"
////////////////////////////////////////////////////////
//Categorys For Stacks 
/////////////////////////////////////////////////////////

const std::vector<NamedCategory<EventCategory>>
EventSelectionGroup_categories = {
NamedCategory<EventCategory>({kUnknown},        "Unknown"),
NamedCategory<EventCategory>({kSignalCCQE},     "Signal-(CCQE)"),
NamedCategory<EventCategory>({kSignalCCMEC},    "Signal-(CCMEC)"),
NamedCategory<EventCategory>({kSignalCCRES},    "Signal-(CCRES)"),
NamedCategory<EventCategory>({kSignalOther},    "Signal-(Other)"),
NamedCategory<EventCategory>({kNuMuCCNpi},      "#nu_{#mu}-CCN#pi"),
//NamedCategory<EventCategory>({kNuMuCC0pi0p},    "#nu_{#mu}-CC0#pi0P"),
NamedCategory<EventCategory>({kNuMuCCOther},    "#nu_{#mu}-CCOther"),
NamedCategory<EventCategory>({kNuECC},          "#nu_{e}-CC"),
NamedCategory<EventCategory>({kNC},             "NC"),
NamedCategory<EventCategory>({kOOFV},           "Out FV"),
NamedCategory<EventCategory>({kOther},          "Other")
};


const std::vector<NamedCategory<BDT_Category>>
BTDGroup_categories = {
NamedCategory<BDT_Category>({kBDT_Else},     "BTD Else"),
NamedCategory<BDT_Category>({kBDT_Muon},     "BDT Muon"),
NamedCategory<BDT_Category>({kBDT_Pion},     "BDT Pion"),
NamedCategory<BDT_Category>({kBDT_Proton},   "BDT Proton"),
NamedCategory<BDT_Category>({kBDT_BOGUS},    "BDT Failed Prediction")

};

const std::vector<NamedCategory<FidVol_Category>>
ContainedGroup_categories = {
NamedCategory<FidVol_Category>({kContained},     "Muon Contained"),
NamedCategory<FidVol_Category>({kUnContained},     "Muon UnContained")
};

const std::vector<NamedCategory<Particle_type>>
ParticleGroup_reduced_categories = {
  NamedCategory<Particle_type>({kParticle_N_A},          "N_A"),
  NamedCategory<Particle_type>({kParticle_OTHER},        "Other"),
  NamedCategory<Particle_type>({kParticle_neutral},      "Neutral"),
  NamedCategory<Particle_type>({kMuon},                   "Muon"),
  NamedCategory<Particle_type>({kPion_0_Electron_kGamma},"e^{-/+},#gamma,#pi^{0}"),
  NamedCategory<Particle_type>({kPion_pos_neg},          "Pi^{-/+}"),
  NamedCategory<Particle_type>({kProton},                "Proton")
};


const std::vector<NamedCategory<CCZeroPi_type>>
topology_categories = {
  NamedCategory<CCZeroPi_type>({kCC_0P_ZeroPi},          "CC0Pi"),
  NamedCategory<CCZeroPi_type>({kCC_1P_ZeroPi},          "CC1p0Pi"),
  NamedCategory<CCZeroPi_type>({kCC_2P_ZeroPi},          "CC2p0Pi"),
  NamedCategory<CCZeroPi_type>({kCC_3orGreaterP_ZeroPi}, "CC3>p0Pi")
};


//////Functions//////////////////////////////
void reorderVector(std::vector<int>& vec, int x);
std::vector<int> MakereorderVector(std::vector<int> vec, int GotoFront);



void hist_maker_MC(int input_Int) {
   ///cc0pi_analyzer_1
 std::string base_dir = "/uboone/data/users/gardiner/ntuples-stv/";
//
 ///uboone/data/users/cnguyen/RootFiles_tutorial/stv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.//1_v08_00_00_26_run1_reco2_reco2.root
 
 ///uboone/data/users/cnguyen/RootFiles_tutorial/stv-prodgenie_bnb_int
 //std::string Rootfile1name = "stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
 //std::string Rootfile1name = "rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
 std::string Rootfile2name = "stv-prodgenie_bnb_intrinsic_nue_overlay_run2_v08_00_00_35_run2a_reco2_reco2.root";
 std::string Rootfile3name = "stv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";
 std::string Rootfile4name = "stv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root"; // TRUTH
 
 
 std::string Rootfile1name; // = "rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";
                                          
 if(input_Int ==0 ){ Rootfile1name = "stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";}
 else if(input_Int ==1 ){ Rootfile1name = "stv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";}
   else{Rootfile1name = "stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";}


 
 
 //std::string Rootfile1name = "rustv-run1_neutrinoselection_filt_numu_ALL.root";
 std::string FullRootFileName = base_dir + Rootfile1name;
  std::string FullRootFileName_OutPut = Rootfile1name;
 std::vector<std::string> filenames{FullRootFileName};
 //std::ifstream myfile(infiles);
 //std::copy(std::istream_iterator<std::string>(myfile),
 //           std::istream_iterator<std::string>(),
 //           std::back_inserter(filenames));
 std::cout << "File Count: " << filenames.size() << std::endl;
 
const auto& EventCategory_tool = EventCategoryInterpreter::Instance();
  TH2D* h_matr_pmult = new TH2D("h_matr_pmult", ";true bin; reconstructed bin; ", 5, -1., 4., 5, -1., 4.);
  h_matr_pmult->SetDirectory(0);
  h_matr_pmult->Sumw2();

  std::vector< double > p_edges = {0.1,0.2,0.275,0.335,0.46,0.61,0.785, 0.985, 1.185, 1.385, 1.585 }; //TODO: Rough measurements based on the "resolution" plot.
  
  //std::vector< double > cos_edges = {-1, -0.85, -0.775, -0.7, -0.625, -0.55, -0.475, -0.4, -0.325,
  //   -0.25, -0.175, -0.1, -0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5,
  //    0.575, 0.65, 0.725, 0.8, 0.85, 0.875, 0.9, 0.925, 0.950, 0.975, 1.00};

std::vector< double > cos_edges = {-1, -0.85, -0.775, -0.7, -0.625, -0.55, -0.475, -0.4, -0.325,
     -0.25, -0.175, -0.1, -0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5,
      0.575, 0.65, 0.725, 0.8, 0.85, 0.875, 0.9, 0.925, 0.950, 0.975, 1.00};





  std::vector< double > pn_edges = {0., 0.125, 0.225, 0.325, 0.425, 0.525, 0.65, 0.85 };
  std::vector< double >delta_alphaT_edges = { 0, 0.35, 0.85, 1.35, 1.85, 2.3, 2.7, 2.95};
  std::vector< double >delta_pTx_edges = { -0.6, -0.45, -0.35, -0.25, -0.15, -0.075, 0, 0.075, 0.15, 0.25,  0.35, 0.45, 0.6};
  std::vector< double >delta_pTy_edges = { -0.8, -0.55, -0.39, -0.2125, -0.05, 0.1, 0.225, 0.3375, 0.5};
  std::vector< double >delta_phiT_edges = {0., 0.075, 0.2, 0.35, 0.5, 0.7, 0.9, 1.15, 1.4, 1.65, 1.9, 2.35,2.8};
  std::vector< double >delta_pT_edges = {0, 0.1, 0.2, 0.3, 0.4, 0.525, 0.675, 0.9};
  std::vector< double >openingAngle_edges = {0, 0.52, 0.78, 1.0, 1.15, 1.35, 1.5, 1.65, 1.8, 1.95, 2.1, 2.35, 2.62};
  std::vector< double >openingAngleDeg_edges = {0.0, 15.0, 30.0, 45.0, 55.0, 65.0, 75.0, 85., 95.0, 105., 115., 120.0, 135., 150.0 };
  std::vector< double > VertexX_edges = {-10,0,25,50,100,125,150,175,200,225,250,275,300,325, 350 };
  std::vector< double > VertexY_edges= {-200,-175,-150,-125,-100.,-75,-50,-25,0,25,50,100,125,150,175,200 };
  std::vector< double > VertexZ_edges= {-10,0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000,1050,1100,1150 ,1200};
  std::vector< double > Mult= {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5};
  std::vector< double > Probability_range= {0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.00 };
  ///std::vector< double > p_proton_edges = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.2}; 
   std::vector< double > p_proton_edges = { 0.25, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.2};
                                         
 std::vector< double > cos_edges_proton = {-1.0, -0.5, 0.0, 0.27, 0.45, 0.62, 0.76, 0.86,0.94, 1.0};

  TH1D* Proton_mult  = new TH1D("Proton_mult", "Proton_mult", Mult.size() - 1, Mult.data());
  
  
  
   ////////////////////
    TH1D* h_leadingProton_p             = new TH1D("h_leadingProton_p", ";p_{p}^{reco}", p_proton_edges.size() - 1, p_proton_edges.data());
     h_leadingProton_p->SetDirectory(0); h_leadingProton_p->Sumw2();
  
  
  HistFolio<TH1D, EventCategory> h_leadingProton_p_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_leadingProton_p_EventCategory", p_proton_edges, "h_leadingProton_p_EventCategory; Events");

  
  TH1D* h_costheta_proton             = new TH1D("h_costheta_proton", ";p_{#mu}^{reco}", cos_edges_proton.size() - 1, cos_edges_proton.data());
   h_costheta_proton->SetDirectory(0); h_costheta_proton->Sumw2();
  
  
    HistFolio<TH1D, EventCategory> h_costheta_proton_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_costheta_proton_EventCategory", cos_edges_proton, "h_costheta_proton_EventCategory; Events");
  
  
HistFolio<TH1D, EventCategory> h_NTrack_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_NTrack_EventCategory", Mult, "h_NTrack_EventCategory; Events");

HistFolio<TH1D, Particle_type> h_NTrack_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_NTrack_Particle_categories", Mult, "h_NTrack_Particle_categories; Events");


//HistFolio<TH1D, BDT_Category> h_MuonTrack_BTDGroup_categories =
//HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_MuonTrack_BTDGroup_categories", Mult, "h_MuonTrack_BTDGroup_categories; Events");
//
//HistFolio<TH1D, BDT_Category> h_Proton1Track_BTDGroup_categories =
//HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_Proton1Track_BTDGroup_categories", Mult, "h_Proton1Track_BTDGroup_categories; Events");
//
//HistFolio<TH1D, BDT_Category> h_Proton2Track_BTDGroup_categories =
//HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_Proton2Track_BTDGroup_categories", Mult, "h_Proton2Track_BTDGroup_categories; Events");


HistFolio<TH1D, FidVol_Category> h_NTrack_containment_categories =
HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_NTrack_containment_categories", Mult, "h_NTrack_containment_categories");


HistFolio<TH1D, Particle_type> h_ProtonMupl_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_ProtonMupl_Particle_categories", Mult, "h_ProtonMupl_Particle_categories; Events");


HistFolio<TH1D, EventCategory> Proton_mult_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "Proton_mult_EventCategory", Mult, "Proton_mult_EventCategory; Events");



  /////////////////////////////////ProtonMupl/////////////////////////////////////////////////
  TH1D* h_p             = new TH1D("h_p", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_true        = new TH1D("h_p_true", ";p_{#mu}^{true}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_Resolution  = new TH1D("h_p_Resolution", ";p_{#mu}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_p_Mig  = new TH2D("h_p_Mig", "", p_edges.size() - 1, p_edges.data(), p_edges.size() - 1, p_edges.data());

  HistFolio<TH1D, EventCategory> h_p_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_p_EventCategory", p_edges, "h_p_EventCategory; Events");
 
 HistFolio<TH1D, Particle_type> h_p_Particle_categories =
 HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_p_Particle_categories", p_edges, "h_p_Particle_categories; Events");
 
 HistFolio<TH1D, FidVol_Category> h_p_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_p_containment_categories", p_edges, "h_p_containment_categories");

HistFolio<TH1D, CCZeroPi_type> h_p_topology_categories =
  HistFolio<TH1D, CCZeroPi_type>(topology_categories, "h_p_topogoly_categories", p_edges, "h_p_topogoly_categories");



  h_p->SetDirectory(0); h_p->Sumw2();
  h_p_true->SetDirectory(0); h_p_true->Sumw2();
  h_p_Resolution->SetDirectory(0); h_p_Resolution->Sumw2();
  h_p_Mig->SetDirectory(0); h_p_Mig->Sumw2(); 
    ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_costheta_p             = new TH1D("h_costheta_p", ";p_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data());
  TH1D* h_costheta_p_true        = new TH1D("h_costheta_p_true", ";p_{#mu}^{true}", cos_edges.size() - 1, cos_edges.data());
  TH1D* h_costheta_p_Resolution  = new TH1D("h_costheta_p_Resolution", ";p_{#mu}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_costheta_p_Mig  = new TH2D("h_costheta_p_Mig", "", cos_edges.size() - 1, cos_edges.data(), cos_edges.size() - 1, cos_edges.data());
  h_costheta_p->SetDirectory(0); h_costheta_p->Sumw2();
  h_costheta_p_true->SetDirectory(0); h_costheta_p_true->Sumw2();
  h_costheta_p_Resolution->SetDirectory(0); h_costheta_p_Resolution->Sumw2();
  h_costheta_p_Mig->SetDirectory(0); h_costheta_p_Mig->Sumw2();
  
  
  HistFolio<TH1D, EventCategory> h_costheta_p_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_costheta_p_EventCategory", cos_edges, "h_costheta_p_EventCategory; Events");

  
  HistFolio<TH1D, Particle_type> h_costheta_p_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_costheta_p_Particle_categories", cos_edges, "h_costheta_p_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_costheta_p_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_costheta_p_containment_categories", cos_edges, "h_costheta_p_containment_categories");
  /////////////////////////////////////////////////////////////////////////  
     TH2D* h_p_costheta_p  = new TH2D("h_p_costheta_p", "h_p_costheta_p",  cos_edges.size() - 1, cos_edges.data(), p_edges.size() - 1, p_edges.data());
   h_p_costheta_p->SetDirectory(0);
   h_p_costheta_p->Sumw2();
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_pn             = new TH1D("h_pn", ";pn_{#mu}^{reco}", pn_edges.size() - 1, pn_edges.data());
  TH1D* h_pn_true        = new TH1D("h_pn_true", ";pn_{#mu}^{true}", pn_edges.size() - 1, pn_edges.data());
  TH1D* h_pn_Resolution  = new TH1D("h_pn_Resolution", ";pn_{#mu}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_pn_Mig  = new TH2D("h_pn_Mig", "", pn_edges.size() - 1, pn_edges.data(), pn_edges.size() - 1, pn_edges.data());

  h_pn->SetDirectory(0);            h_pn->Sumw2();
  h_pn_true->SetDirectory(0);       h_pn_true->Sumw2();
  h_pn_Resolution->SetDirectory(0); h_pn_Resolution->Sumw2();
  h_pn_Mig->SetDirectory(0); h_pn_Mig->Sumw2();
  
  HistFolio<TH1D, EventCategory> h_pn_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_pn_EventCategory", pn_edges, "h_pn_EventCategory; Events");


  HistFolio<TH1D, Particle_type> h_pn_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_pn_Particle_categories", pn_edges, "h_pn_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_pn_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_pn_containment_categories", pn_edges, "h_pn_containment_categories");
    
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_alphaT             = new TH1D("h_delta_alphaT", ";#delta_{#alpha T}^{reco}", delta_alphaT_edges.size() - 1, delta_alphaT_edges.data());
  TH1D* h_delta_alphaT_true        = new TH1D("h_delta_alphaT_true", ";#delta_{#alpha T}^{true}", delta_alphaT_edges.size() - 1, delta_alphaT_edges.data());
  TH1D* h_delta_alphaT_Resolution  = new TH1D("h_delta_alphaT_Resolution", ";#delta_{#alpha T}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_delta_alphaT_Mig  = new TH2D("h_delta_alphaT_Mig", "", delta_alphaT_edges.size() - 1, delta_alphaT_edges.data(), delta_alphaT_edges.size() - 1, delta_alphaT_edges.data());

  h_delta_alphaT->SetDirectory(0);            h_delta_alphaT->Sumw2();
  h_delta_alphaT_true->SetDirectory(0);       h_delta_alphaT_true->Sumw2();
  h_delta_alphaT_Resolution->SetDirectory(0); h_delta_alphaT_Resolution->Sumw2();
  h_delta_alphaT_Mig->SetDirectory(0); h_delta_alphaT_Mig->Sumw2();
  
  HistFolio<TH1D, EventCategory> h_delta_alphaT_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_delta_alphaT_EventCategory", delta_alphaT_edges, "h_delta_alphaT_EventCategory; Events");

 
  HistFolio<TH1D, Particle_type> h_delta_alphaT_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_delta_alphaT_Particle_categories", delta_alphaT_edges, "h_delta_alphaT_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_delta_alphaT_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_delta_alphaT_containment_categories", delta_alphaT_edges, "h_delta_alphaT_containment_categories");
    
  
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_pTx             = new TH1D("h_delta_pTx", ";#delta_{#alpha Tx}^{reco}", delta_pTx_edges.size() - 1, delta_pTx_edges.data());
  TH1D* h_delta_pTx_true        = new TH1D("h_delta_pTx_true", ";#delta_{#alpha Tx}^{true}", delta_pTx_edges.size() - 1, delta_pTx_edges.data());
  TH1D* h_delta_pTx_Resolution  = new TH1D("h_delta_pTx_Resolution", ";#delta_{#alpha Tx}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_delta_pTx_Mig  = new TH2D("h_delta_pTx_Mig", "", delta_pTx_edges.size() - 1, delta_pTx_edges.data(), delta_pTx_edges.size() - 1, delta_pTx_edges.data());

  h_delta_pTx->SetDirectory(0);            h_delta_pTx->Sumw2();
  h_delta_pTx_true->SetDirectory(0);       h_delta_pTx_true->Sumw2();
  h_delta_pTx_Resolution->SetDirectory(0); h_delta_pTx_Resolution->Sumw2();
  h_delta_pTx_Mig->SetDirectory(0); h_delta_pTx_Mig->Sumw2();
  
  HistFolio<TH1D, EventCategory> h_delta_pTx_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_delta_pTx_EventCategory", delta_pTx_edges, "h_delta_pTx_EventCategory; Events");

  
  HistFolio<TH1D, Particle_type> h_delta_pTx_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_delta_pTx_Particle_categories", delta_pTx_edges, "h_delta_pTx_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_delta_pTx_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_delta_pTx_containment_categories", delta_pTx_edges, "h_delta_pTx_containment_categories");

  
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_pTy             = new TH1D("h_delta_pTy", ";#delta_{#alpha T}^{reco}", delta_pTy_edges.size() - 1, delta_pTy_edges.data());
  TH1D* h_delta_pTy_true        = new TH1D("h_delta_pTy_true", ";#delta_{#alpha T}^{true}", delta_pTy_edges.size() - 1, delta_pTy_edges.data());
  TH1D* h_delta_pTy_Resolution  = new TH1D("h_delta_pTy_Resolution", ";#delta_{#alpha T}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_delta_pTy_Mig  = new TH2D("h_delta_pTy_Mig", "", delta_pTy_edges.size() - 1, delta_pTy_edges.data(), delta_pTy_edges.size() - 1, delta_pTy_edges.data());
  h_delta_pTy->SetDirectory(0);            h_delta_pTy->Sumw2();
  h_delta_pTy_true->SetDirectory(0);       h_delta_pTy_true->Sumw2();
  h_delta_pTy_Resolution->SetDirectory(0); h_delta_pTy_Resolution->Sumw2();
  h_delta_pTy_Mig->SetDirectory(0); h_delta_pTy_Mig->Sumw2();


  HistFolio<TH1D, EventCategory> h_delta_pTy_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_delta_pTy_EventCategory", delta_pTy_edges, "h_delta_pTy_EventCategory; Events");

   
  HistFolio<TH1D, Particle_type> h_delta_pTy_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_delta_pTy_Particle_categories", delta_pTy_edges, "h_delta_pTy_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_delta_pTy_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_delta_pTy_containment_categories", delta_pTy_edges, "h_delta_pTy_containment_categories");

////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_phiT             = new TH1D("h_delta_phiT", "", delta_phiT_edges.size() - 1, delta_phiT_edges.data());
  TH1D* h_delta_phiT_true        = new TH1D("h_delta_phiT_true", "", delta_phiT_edges.size() - 1, delta_phiT_edges.data());
  TH1D* h_delta_phiT_Resolution  = new TH1D("h_delta_phiT_Resolution", ";", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_delta_phiT_Mig  = new TH2D("h_delta_phiT_Mig", "", delta_phiT_edges.size() - 1, delta_phiT_edges.data(), delta_phiT_edges.size() - 1, delta_phiT_edges.data());

  h_delta_phiT->SetDirectory(0);            h_delta_phiT->Sumw2();
  h_delta_phiT_true->SetDirectory(0);       h_delta_phiT_true->Sumw2();
  h_delta_phiT_Resolution->SetDirectory(0); h_delta_phiT_Resolution->Sumw2();
  h_delta_phiT_Mig->SetDirectory(0); h_delta_phiT_Mig->Sumw2();
  
  
    
  HistFolio<TH1D, EventCategory> h_delta_phiT_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_delta_phiT_EventCategory", delta_phiT_edges, "h_delta_phiT_EventCategory; Events");
  

  HistFolio<TH1D, Particle_type> h_delta_phiT_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_delta_phiT_Particle_categories", delta_phiT_edges, "h_delta_phiT_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_delta_phiT_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_delta_phiT_containment_categories", delta_phiT_edges, "h_delta_phiT_containment_categories");
  
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_openingAngle             = new TH1D("h_openingAngle", "", openingAngle_edges.size() - 1, openingAngle_edges.data());
  TH1D* h_openingAngle_true        = new TH1D("h_openingAngle_true", "", openingAngle_edges.size() - 1, openingAngle_edges.data());
  TH1D* h_openingAngle_Resolution  = new TH1D("h_openingAngle_Resolution", ";", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_openingAngle_Mig  = new TH2D("h_openingAngle_Mig", "", openingAngle_edges.size() - 1, openingAngle_edges.data(), openingAngle_edges.size() - 1, openingAngle_edges.data());
  h_openingAngle->SetDirectory(0);            h_openingAngle->Sumw2();
  h_openingAngle_true->SetDirectory(0);       h_openingAngle_true->Sumw2();
  h_openingAngle_Resolution->SetDirectory(0); h_openingAngle_Resolution->Sumw2();
  h_openingAngle_Mig->SetDirectory(0); h_openingAngle_Mig->Sumw2();
  
  HistFolio<TH1D, EventCategory> h_openingAngle_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_openingAngle_EventCategory", openingAngle_edges, "h_openingAngle_EventCategory; Events");
   
   HistFolio<TH1D, Particle_type> h_openingAngle_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_openingAngle_Particle_categories", openingAngle_edges, "h_openingAngle_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_openingAngle_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_openingAngle_containment_categories", openingAngle_edges, "h_openingAngle_containment_categories");
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_openingAngle_deg             = new TH1D("h_openingAngle_deg", "", openingAngleDeg_edges.size() - 1, openingAngleDeg_edges.data());
  TH1D* h_openingAngle_deg_true        = new TH1D("h_openingAngle_deg_true", "", openingAngleDeg_edges.size() - 1, openingAngleDeg_edges.data());
  TH1D* h_openingAngle_deg_Resolution  = new TH1D("h_openingAngle_deg_Resolution", ";", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_openingAngle_Deg_Mig  = new TH2D("h_openingAngle_Deg_Mig", "", openingAngleDeg_edges.size() - 1, openingAngleDeg_edges.data(), openingAngleDeg_edges.size() - 1, openingAngleDeg_edges.data());
  h_openingAngle_deg->SetDirectory(0);            h_openingAngle_deg->Sumw2();
  h_openingAngle_deg_true->SetDirectory(0);       h_openingAngle_deg_true->Sumw2();
  h_openingAngle_deg_Resolution->SetDirectory(0); h_openingAngle_deg_Resolution->Sumw2();
  h_openingAngle_Deg_Mig->SetDirectory(0); h_openingAngle_Deg_Mig->Sumw2();
  
  
  HistFolio<TH1D, EventCategory> h_openingAngle_Deg_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_openingAngle_Deg_EventCategory", openingAngleDeg_edges, "h_openingAngle_Deg_EventCategory; Events");
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexX             = new TH1D("h_VertexX", "", VertexX_edges.size() - 1, VertexX_edges.data());
  TH1D* h_VertexX_true        = new TH1D("h_VertexX_true", "", VertexX_edges.size() - 1, VertexX_edges.data());
  TH1D* h_VertexX_Resolution  = new TH1D("h_VertexX_Resolution", ";", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_VertexX_Mig  = new TH2D("h_VertexX_Mig", "", VertexX_edges.size() - 1, VertexX_edges.data(), VertexX_edges.size() - 1, VertexX_edges.data());
  h_VertexX->SetDirectory(0);            h_VertexX->Sumw2();
  h_VertexX_true->SetDirectory(0);       h_VertexX_true->Sumw2();
  h_VertexX_Resolution->SetDirectory(0); h_VertexX_Resolution->Sumw2();
  h_VertexX_Mig->SetDirectory(0); h_VertexX_Mig->Sumw2();
  
  
    HistFolio<TH1D, EventCategory> h_VertexX_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_VertexX_EventCategory", VertexX_edges, "h_VertexX_EventCategory; Events");
  
  HistFolio<TH1D, Particle_type> h_VertexX_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_VertexX_Particle_categories", VertexX_edges, "h_VertexX_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_VertexX_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_VertexX_containment_categories", VertexX_edges, "h_VertexX_containment_categories");
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexY             = new TH1D("h_VertexY", "", VertexY_edges.size() - 1, VertexY_edges.data());
  TH1D* h_VertexY_true        = new TH1D("h_VertexY_true", "", VertexY_edges.size() - 1, VertexY_edges.data());
  TH1D* h_VertexY_Resolution  = new TH1D("h_VertexY_Resolution", ";", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_VertexY_Mig         = new TH2D("h_VertexY_Mig", "", VertexY_edges.size() - 1, VertexY_edges.data(), VertexY_edges.size() - 1, VertexY_edges.data());
  h_VertexY->SetDirectory(0);            h_VertexY->Sumw2();
  h_VertexY_true->SetDirectory(0);       h_VertexY_true->Sumw2();
  h_VertexY_Resolution->SetDirectory(0); h_VertexY_Resolution->Sumw2();
  h_VertexY_Mig->SetDirectory(0); h_VertexY_Mig->Sumw2(); 
  
   HistFolio<TH1D, EventCategory> h_VertexY_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_VertexY_EventCategory", VertexY_edges, "h_VertexY_EventCategory; Events");
  
  
  HistFolio<TH1D, Particle_type> h_VertexY_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_VertexY_Particle_categories", VertexY_edges, "h_VertexY_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_VertexY_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_VertexY_containment_categories", VertexY_edges, "h_VertexY_containment_categories");
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexZ             = new TH1D("h_VertexZ", "", VertexZ_edges.size() - 1, VertexZ_edges.data());
  TH1D* h_VertexZ_true        = new TH1D("h_VertexZ_true", "", VertexZ_edges.size() - 1, VertexZ_edges.data());
  TH1D* h_VertexZ_Resolution  = new TH1D("h_VertexZ_Resolution", ";", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_VertexZ_Mig         = new TH2D("h_VertexZ_Mig", "", VertexZ_edges.size() - 1, VertexZ_edges.data(), VertexZ_edges.size() - 1, VertexZ_edges.data());

  h_VertexZ->SetDirectory(0);            h_VertexZ->Sumw2();
  h_VertexZ_true->SetDirectory(0);       h_VertexZ_true->Sumw2();
  h_VertexZ_Resolution->SetDirectory(0); h_VertexZ_Resolution->Sumw2();
  h_VertexZ_Mig->SetDirectory(0); h_VertexZ_Mig->Sumw2();
   HistFolio<TH1D, EventCategory> h_VertexZ_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_VertexZ_EventCategory", VertexY_edges, "h_VertexZ_EventCategory; Events");
  
  
  HistFolio<TH1D, Particle_type> h_VertexZ_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_VertexZ_Particle_categories", VertexY_edges, "h_VertexZ_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_VertexZ_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_VertexZ_containment_categories", VertexY_edges, "h_VertexZ_containment_categories");

  
  
  ////////////////////////////
  //2D vectex
  ///////////////////////////
  TH2D* h_VertexX_Z  = new TH2D("h_VertexX_Z", "", VertexZ_edges.size() - 1, VertexZ_edges.data(), VertexX_edges.size() - 1, VertexX_edges.data());
  h_VertexX_Z->SetDirectory(0);
  h_VertexX_Z->Sumw2();
  TH2D* h_VertexY_Z  = new TH2D("h_VertexY_Z", "", VertexZ_edges.size() - 1, VertexZ_edges.data(), VertexY_edges.size() - 1, VertexY_edges.data());
  h_VertexY_Z->SetDirectory(0);
  h_VertexY_Z->Sumw2();
  TH2D* h_VertexX_Y  = new TH2D("h_VertexX_Y", "", VertexX_edges.size() - 1, VertexX_edges.data(), VertexY_edges.size() - 1, VertexY_edges.data());
 h_VertexX_Y->SetDirectory(0);
 h_VertexX_Y->Sumw2();


  double res_max = -0.01;


  for (auto& infile : filenames) {
    std::cout << "Input File Name: " << infile << std::endl;

    TFile* f = TFile::Open(infile.c_str());
    if (!(f && f->IsOpen() && !f->IsZombie())) {
        std::cout << "Bad file! " << infile << std::endl;
        continue;
    }
    TTree* stv  = (TTree*) f->Get("stv_tree");
    if (!(stv && stv->GetEntries())) {
       std::cout << "Bad file! " << infile << std::endl;
       continue;
    }





  bool mc_signal, sel_signal, sel_muon_contained;
  TVector3 *p3_mu=nullptr, *mc_p3_mu=nullptr;
  TVector3 *p3_lead_p=nullptr, *mc_p3_lead_p=nullptr;
  

  
  std::vector<int> *mc_pdg = nullptr;
  std::vector<int> *pfpdg = nullptr;
  std::vector<int> *backtracked_pdg = nullptr;
  
  std::vector<int> *trk_BDT_PID = nullptr;
  std::vector<int> *trk_BDT_Probility = nullptr;
  
  std::vector<double> *weight_TunedCentralValue_UBGenie = nullptr;
  std::vector<double> *weight_splines_general_Spline = nullptr;
  
  float  spline_weight, tuned_cv_weight; //weight_TunedCentralValue_UBGenie, weight_splines_general_Spline
  float pn, delta_alphaT, delta_pTx, delta_pTy, delta_pL, delta_phiT, delta_pT;
  float reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z;
  float mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z;
  float mc_pn, mc_delta_pL, mc_delta_alphaT, mc_delta_pTx, mc_delta_pTy, mc_delta_p, mc_delta_phiT, mc_delta_pT;
  int mc_interaction,proton_multiplicity,pion_multiplicity,category;  


  //std::vector<float> *trk_BDT_score_0 = nullptr;
  //std::vector<float> *trk_BDT_score_1 = nullptr;
  //std::vector<float> *trk_BDT_score_2 = nullptr;
  //std::vector<float> *trk_BDT_score_3 = nullptr;




  stv->SetBranchStatus("*", 0);
  stv->SetBranchStatus("sel_CCNp0pi", 1);  
  stv->SetBranchStatus("mc_p3_mu", 1);
  stv->SetBranchStatus("p3_mu", 1);
  stv->SetBranchStatus("mc_p3_lead_p", 1);
  stv->SetBranchStatus("p3_lead_p", 1);
  stv->SetBranchStatus("mc_pdg", 1);
  stv->SetBranchStatus("sel_muon_contained", 1);
  stv->SetBranchStatus("spline_weight", 1);
  stv->SetBranchStatus("tuned_cv_weight", 1);
  
  stv->SetBranchStatus("weight_splines_general_Spline", 1);  
  stv->SetBranchStatus("weight_TunedCentralValue_UBGenie", 1);
  stv->SetBranchStatus("pn", 1);
  stv->SetBranchStatus("delta_alphaT",1);
  stv->SetBranchStatus("delta_pTx",1);
  stv->SetBranchStatus("delta_pTy",1);
  stv->SetBranchStatus("delta_pL",1);
  stv->SetBranchStatus("delta_phiT",1);
  stv->SetBranchStatus("delta_pT",1);
  stv->SetBranchStatus("reco_nu_vtx_sce_x",1);
  stv->SetBranchStatus("reco_nu_vtx_sce_y",1);
  stv->SetBranchStatus("reco_nu_vtx_sce_z",1); 

  //stv->SetBranchStatus("proton_multiplicity",1); 
  //stv->SetBranchStatus("pion_multiplicity",1); 
  //stv->SetBranchStatus("trk_BDT_PID",1); 
  //stv->SetBranchStatus("n_pfps",1); 
  stv->SetBranchStatus("pfpdg",1);
  stv->SetBranchStatus("backtracked_pdg",1);
  //stv->SetBranchStatus("trk_BDT_score_0",1);
  //stv->SetBranchStatus("trk_BDT_score_1",1);
  //stv->SetBranchStatus("trk_BDT_score_2",1);
  //stv->SetBranchStatus("trk_BDT_score_3",1);



  stv->SetBranchStatus("mc_pn", 1);
  stv->SetBranchStatus("mc_delta_alphaT",1);
  stv->SetBranchStatus("mc_delta_pTx",1);
  stv->SetBranchStatus("mc_delta_pTy",1);
  stv->SetBranchStatus("mc_delta_pL",1);
  stv->SetBranchStatus("mc_delta_phiT",1);
  stv->SetBranchStatus("mc_delta_pT",1);
  stv->SetBranchStatus("mc_nu_vtx_x",1);
  stv->SetBranchStatus("mc_nu_vtx_y",1);
  stv->SetBranchStatus("mc_nu_vtx_z",1); 
  stv->SetBranchStatus("mc_is_signal",1); 


  stv->SetBranchStatus("mc_interaction",1); 
  stv->SetBranchStatus("category",1);
  stv->SetBranchAddress("category", &category);
  
  
  stv->SetBranchAddress("mc_interaction", &mc_interaction);

  stv->SetBranchAddress("mc_is_signal", &mc_signal);
  stv->SetBranchAddress("sel_CCNp0pi", &sel_signal);
  stv->SetBranchAddress("p3_mu", &p3_mu);
  stv->SetBranchAddress("mc_p3_mu", &mc_p3_mu);
  
  stv->SetBranchAddress("mc_p3_lead_p", &mc_p3_lead_p);
  stv->SetBranchAddress("p3_lead_p", &p3_lead_p);
  
  stv->SetBranchAddress("mc_pdg", &mc_pdg);
  stv->SetBranchAddress("sel_muon_contained", &sel_muon_contained);
  stv->SetBranchAddress("spline_weight", &spline_weight);
  stv->SetBranchAddress("weight_splines_general_Spline", &weight_splines_general_Spline);
  stv->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
  stv->SetBranchAddress("weight_TunedCentralValue_UBGenie", &weight_TunedCentralValue_UBGenie); 
/////////////////////////
//// MC RECO
/////////////////////////
  stv->SetBranchAddress("pn", &pn);
  stv->SetBranchAddress("delta_alphaT",&delta_alphaT);
  stv->SetBranchAddress("delta_pTx",&delta_pTx);
  stv->SetBranchAddress("delta_pTy",&delta_pTy);
  stv->SetBranchAddress("delta_pL",&delta_pL);
  stv->SetBranchAddress("delta_phiT",&delta_phiT);
  stv->SetBranchAddress("delta_pT",&delta_pT);
  stv->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x);
  stv->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y);
  stv->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z); 
  //stv->SetBranchAddress("proton_multiplicity", &proton_multiplicity);
  //stv->SetBranchAddress("pion_multiplicity", &pion_multiplicity);
  //stv->SetBranchAddress("trk_BDT_PID", &trk_BDT_PID);
  ////stv->SetBranchAddress("n_pfps", &n_pfps);
  stv->SetBranchAddress("pfpdg", &pfpdg);
  stv->SetBranchAddress("backtracked_pdg", &backtracked_pdg);


 //stv->SetBranchAddress("trk_BDT_score_0", &trk_BDT_score_0);
 //stv->SetBranchAddress("trk_BDT_score_1", &trk_BDT_score_1);
 //stv->SetBranchAddress("trk_BDT_score_2", &trk_BDT_score_2);
 //stv->SetBranchAddress("trk_BDT_score_3", &trk_BDT_score_3);


/////////////////////////
//// MC TRUE
/////////////////////////
  stv->SetBranchAddress("mc_pn", &mc_pn);
  stv->SetBranchAddress("mc_delta_alphaT",&mc_delta_alphaT);
  stv->SetBranchAddress("mc_delta_pTx",&mc_delta_pTx);
  stv->SetBranchAddress("mc_delta_pTy",&mc_delta_pTy);
  stv->SetBranchAddress("mc_delta_pL",&mc_delta_pL);
  stv->SetBranchAddress("mc_delta_phiT",&mc_delta_phiT);
  stv->SetBranchAddress("mc_delta_pT",&mc_delta_pT);
  stv->SetBranchAddress("mc_nu_vtx_x", &mc_nu_vtx_x);
  stv->SetBranchAddress("mc_nu_vtx_y", &mc_nu_vtx_y);
  stv->SetBranchAddress("mc_nu_vtx_z", &mc_nu_vtx_z); 

  //WeightHandler wh;
  //wh.set_branch_addresses( input_chain_, universe_branch_names );
  //
  //wh.add_branch( input_chain_, "weight_splines_general_Spline", false );
  //wh.add_branch( input_chain_, "weight_TunedCentralValue_UBGenie", false );

  
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout<< "~~~~~~~~~~~~~~~ STarting Loop~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout << " The Number of Entries "<<  stv->GetEntries() << std::endl;
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;


  for (long i=0; i<stv->GetEntries(); i++)
  {
    stv->GetEntry(i);
   // if (sel_signal && mc_signal){
     if ( i % 5 == 0 ) {
      std::cout << "Processing event #" << i << '\n';
    }
      
     if(sel_signal){
   //&& mc_signal
    //std::cout<<"pn = " << pn << std::endl;


  EventCategory TrueEventType = EventCategory_tool.returnEventCategoryType(category);
  FidVol_Category Containmenttype = EventCategory_tool.GetContainmentType(sel_muon_contained,TrueEventType);
  //CCZeroPi_type TopologyType = EventCategory_tool.returnCCZeroPi_type(proton_multiplicity);
std::vector<double> vector_spline  = *weight_splines_general_Spline; 
std::vector<double> vector_Genie  = *weight_TunedCentralValue_UBGenie;



double check1 = weight_splines_general_Spline->front();
double check2 = weight_TunedCentralValue_UBGenie->front();
  double Wgt = spline_weight * tuned_cv_weight;
 
  double Wgt2 = check1* check2;
  double FINAL_Wgt = safe_weight( Wgt );
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<< " spline_weight =  "  << spline_weight << "   tuned_cv_weight = "<< tuned_cv_weight << " mult = " << Wgt<< std::endl;
  std::cout<< " check =  "  << check1 << "   check2 = "<< check2 << " Wgt2 = " << Wgt2<< " Saft wgt = "<<FINAL_Wgt<< std::endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<" "<<std::endl;
   std::cout<<" "<<std::endl;
   
   
   
   
   
   
  double Pmu = p3_mu->Mag();
  double Pmu_true = mc_p3_mu->Mag();
  double Pmu_Resolution = (Pmu_true - Pmu) / Pmu_true;
  
  double CosTheta = p3_mu->CosTheta();
  double CosTheta_true = mc_p3_mu->CosTheta();
  double CosTheta_Resolution =CosTheta - CosTheta_true;
  double P_proton = p3_lead_p->Mag();
  double P_CosTheta = p3_lead_p->CosTheta();
  
  double OpenAngle = TMath::ACos( (p3_mu->X()*p3_lead_p->X() + p3_mu->Y()*p3_lead_p->Y() + p3_mu->Z()*p3_lead_p->Z()) / p3_mu->Mag());
  double OpenAngle_deg = OpenAngle * TMath::RadToDeg();
  double mc_OpenAngle = TMath::ACos( (mc_p3_mu->X()*mc_p3_lead_p->X() + mc_p3_mu->Y()*p3_lead_p->Y() + mc_p3_mu->Z()*mc_p3_lead_p->Z()) / mc_p3_mu->Mag());
  double mc_OpenAngle_deg = mc_OpenAngle* TMath::RadToDeg();
  double OpenAngle_resolution = OpenAngle - mc_OpenAngle;
  int NTracks = proton_multiplicity + 1; // assuming the number of tracks is the Proton Multipleicity + muon track 
  
  
  h_p->Fill(Pmu,FINAL_Wgt); 
  h_p_true->Fill(Pmu_true,FINAL_Wgt); 
  h_p_Resolution->Fill(Pmu_Resolution,FINAL_Wgt); 
  h_p_Mig->Fill(Pmu_true,Pmu,FINAL_Wgt); 
  h_p_EventCategory.GetComponentHist(TrueEventType)->Fill(Pmu,FINAL_Wgt);
  h_p_containment_categories.GetComponentHist(Containmenttype)->Fill(Pmu,FINAL_Wgt);
   
  
  
  h_costheta_p->Fill(CosTheta,FINAL_Wgt);
  h_costheta_p_true->Fill(CosTheta_true,FINAL_Wgt);
  h_costheta_p_Resolution->Fill(CosTheta_Resolution,FINAL_Wgt); 
  h_costheta_p_Mig->Fill(CosTheta_true,CosTheta,FINAL_Wgt); 
  
  
  h_costheta_p_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,FINAL_Wgt);
  h_costheta_p_containment_categories.GetComponentHist(Containmenttype)->Fill(CosTheta,FINAL_Wgt);
  
  
  //h_NTrack_EventCategory.GetComponentHist(TrueEventType)->Fill(NTracks,FINAL_Wgt);
  
  h_p_costheta_p->Fill(CosTheta,Pmu,FINAL_Wgt); 
   
   
  
  
  h_VertexX->Fill(reco_nu_vtx_sce_x,FINAL_Wgt);  
  h_VertexX_true->Fill(mc_nu_vtx_x,FINAL_Wgt);
  h_VertexX_Resolution->Fill(mc_nu_vtx_x-reco_nu_vtx_sce_x,FINAL_Wgt);
  h_VertexX_Mig->Fill(mc_nu_vtx_x,reco_nu_vtx_sce_x,FINAL_Wgt); 
  h_VertexY->Fill(reco_nu_vtx_sce_y,FINAL_Wgt);  
  h_VertexY_true->Fill(mc_nu_vtx_y,FINAL_Wgt);
  h_VertexY_Resolution->Fill(mc_nu_vtx_y-reco_nu_vtx_sce_y,FINAL_Wgt);
  h_VertexY_Mig->Fill(mc_nu_vtx_y,reco_nu_vtx_sce_y,FINAL_Wgt); 
  h_VertexZ->Fill(reco_nu_vtx_sce_z,FINAL_Wgt);  
  h_VertexZ_true->Fill(mc_nu_vtx_z,FINAL_Wgt);
  h_VertexZ_Resolution->Fill(mc_nu_vtx_z-reco_nu_vtx_sce_z,FINAL_Wgt);
  h_VertexZ_Mig->Fill(mc_nu_vtx_z,reco_nu_vtx_sce_z,FINAL_Wgt); 
  h_VertexX_Z->Fill(reco_nu_vtx_sce_z, reco_nu_vtx_sce_x,FINAL_Wgt); 
  h_VertexY_Z->Fill(reco_nu_vtx_sce_z, reco_nu_vtx_sce_y,FINAL_Wgt); 
  h_VertexX_Y->Fill(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y,FINAL_Wgt); 
  
  
  h_openingAngle->Fill(OpenAngle,FINAL_Wgt);  
  h_openingAngle_true->Fill(mc_OpenAngle,FINAL_Wgt);
  h_openingAngle_Resolution->Fill(OpenAngle_resolution,FINAL_Wgt);
  h_openingAngle_Mig->Fill(mc_OpenAngle,OpenAngle,FINAL_Wgt); 
  
  h_openingAngle_deg->Fill(OpenAngle_deg,FINAL_Wgt);  
  h_openingAngle_deg_true->Fill(mc_OpenAngle_deg,FINAL_Wgt); 
  
  h_pn->Fill(pn,FINAL_Wgt); 
  h_pn_true->Fill(mc_pn,FINAL_Wgt);
  h_pn_Resolution->Fill(pn-mc_pn,FINAL_Wgt);
  h_pn_Mig->Fill(mc_pn,pn,FINAL_Wgt); 
  h_pn_EventCategory.GetComponentHist(TrueEventType)->Fill(pn,FINAL_Wgt);
  
  
  
  h_delta_alphaT->Fill(delta_alphaT,FINAL_Wgt); 
  h_delta_alphaT_true->Fill(mc_delta_alphaT,FINAL_Wgt);
  h_delta_alphaT_Resolution->Fill(delta_alphaT-mc_delta_alphaT,FINAL_Wgt);
  h_delta_alphaT_Mig->Fill(mc_delta_alphaT,delta_alphaT,FINAL_Wgt); 
  h_delta_alphaT_EventCategory.GetComponentHist(TrueEventType)->Fill(delta_alphaT,FINAL_Wgt);
  
  
  h_delta_pTx->Fill(delta_pTx,FINAL_Wgt);  
  h_delta_pTx_true->Fill(mc_delta_pTx,FINAL_Wgt);
  h_delta_pTx_Resolution->Fill(delta_pTx-mc_delta_pTx,FINAL_Wgt);
  h_delta_pTx_Mig->Fill(mc_delta_pTx,delta_pTx,FINAL_Wgt); 
  h_delta_pTx_EventCategory.GetComponentHist(TrueEventType)->Fill(delta_pTx,FINAL_Wgt);
  
  
  
  h_delta_pTy->Fill(delta_pTy,FINAL_Wgt);  
  h_delta_pTy_true->Fill(mc_delta_pTy,FINAL_Wgt);
  h_delta_pTy_Resolution->Fill(delta_pTy-mc_delta_pTy,FINAL_Wgt);
  h_delta_pTy_Mig->Fill(mc_delta_pTy,delta_pTy,FINAL_Wgt); 
  h_delta_pTy_EventCategory.GetComponentHist(TrueEventType)->Fill(delta_pTy,FINAL_Wgt);
  
  
  h_delta_phiT->Fill(delta_phiT,FINAL_Wgt);  
  h_delta_phiT_true->Fill(mc_delta_phiT,FINAL_Wgt);
  h_delta_phiT_Resolution->Fill(delta_phiT-mc_delta_phiT,FINAL_Wgt);
  h_delta_phiT_Mig->Fill(mc_delta_phiT,delta_phiT,FINAL_Wgt); 
  h_delta_phiT_EventCategory.GetComponentHist(TrueEventType)->Fill(delta_phiT,FINAL_Wgt);
  
  h_leadingProton_p->Fill(P_proton,FINAL_Wgt);
  h_leadingProton_p_EventCategory.GetComponentHist(TrueEventType)->Fill(P_proton,FINAL_Wgt);
  h_costheta_proton_EventCategory.GetComponentHist(TrueEventType)->Fill(P_CosTheta,FINAL_Wgt);
  h_costheta_proton->Fill(P_CosTheta,FINAL_Wgt);
   
  
  }// ENd of If Statement, that passed Signal Selection 
//////////////////////////////////////////////////////  
}// End of Event Loop
///////////////////////////////////////////////////////
f->Close();







}


/////////////////////////
// Make Root File to look 
/////////////////////////
std::string OutPutRootName = "/uboone/data/users/cnguyen/CC0Pi_Selection/CCNP_0Pion_steven_test/MC_Selected_" + FullRootFileName_OutPut;
std::cout<<" Writing Root File :"<< OutPutRootName << std::endl;

auto outFile = TFile::Open(OutPutRootName.c_str(), "RECREATE");
outFile->cd(); 

h_p->Write();
h_p_true->Write();
h_p_Resolution->Write();
h_p_Mig->Write();

h_costheta_p->Write();
h_costheta_p_true->Write();
h_costheta_p_Resolution->Write();
h_costheta_p_Mig->Write();
h_p_costheta_p->Write();

h_pn->Write();
h_pn_true->Write();
h_pn_Resolution->Write();
h_pn_Mig->Write();

h_p_EventCategory.WriteToFile(*outFile);
h_p_containment_categories.WriteToFile(*outFile);
h_p_topology_categories.WriteToFile(*outFile);

h_NTrack_EventCategory.WriteToFile(*outFile);

h_costheta_p_EventCategory.WriteToFile(*outFile);
h_costheta_p_containment_categories.WriteToFile(*outFile);


h_NTrack_Particle_categories.WriteToFile(*outFile);







h_pn_EventCategory.WriteToFile(*outFile);
h_delta_alphaT_EventCategory.WriteToFile(*outFile);
h_delta_pTx_EventCategory.WriteToFile(*outFile);
h_delta_pTy_EventCategory.WriteToFile(*outFile);
h_delta_phiT_EventCategory.WriteToFile(*outFile);

h_leadingProton_p->Write();
h_leadingProton_p_EventCategory.WriteToFile(*outFile);
h_costheta_proton_EventCategory.WriteToFile(*outFile);
h_costheta_proton->Write();


h_delta_alphaT->Write();
h_delta_alphaT_true->Write();
h_delta_alphaT_Resolution->Write();
h_delta_alphaT_Mig->Write();

h_delta_pTx->Write();
h_delta_pTx_true->Write();
h_delta_pTx_Resolution->Write();
h_delta_pTx_Mig->Write();

h_delta_pTy->Write();
h_delta_pTy_true->Write();
h_delta_pTy_Resolution->Write();
h_delta_pTy_Mig->Write();

h_delta_phiT->Write();
h_delta_phiT_true->Write();
h_delta_phiT_Resolution->Write();
h_delta_phiT_Mig->Write();

h_openingAngle->Write();
h_openingAngle_true->Write();
h_openingAngle_Resolution->Write();
h_openingAngle_Mig->Write();
h_openingAngle_deg->Write();
h_openingAngle_deg_true->Write();

h_VertexX->Write();
h_VertexX_true->Write();
h_VertexX_Resolution->Write();
h_VertexX_Mig->Write();
h_VertexY->Write();
h_VertexY_true->Write();
h_VertexY_Resolution->Write();
h_VertexY_Mig->Write();
h_VertexZ->Write();
h_VertexZ_true->Write();
h_VertexZ_Resolution->Write();
h_VertexZ_Mig->Write();
h_VertexX_Z->Write();
h_VertexY_Z->Write();
h_VertexX_Y->Write();






outFile->Close();

}
int main(int argc, char* argv[]) {

 int myInt = std::atoi(argv[1]);

   hist_maker_MC(myInt);
   return 0;
}


void reorderVector(std::vector<int>& vec, int x) {
    // Find all occurrences of x in the vector
    
    std::vector<int>::iterator it = std::find(vec.begin(), vec.end(), x);
    
    int index = 0;

    while (it != vec.end()) {
        // Swap the found element with the element at index 0
        std::swap(*it, vec[index]);

        // Find the next occurrence of x starting from the next element
        it = std::find(++it, vec.end(), x);
        index++;
    }
}


std::vector<int> MakereorderVector(std::vector<int> vec, int GotoFront){
  std::vector<int>output_vector; 
  for(auto value : vec ){output_vector.push_back(value);}
  reorderVector(output_vector, GotoFront);
  return output_vector;
}

