/*Author: Christian
 *
 * Usage: making hist to plot and look at for CC-0pi
 *
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
#include "includes/UBTH2Poly.h"
#include "includes/HistUtils_cc0pi.hh"
//#include "HistMaker.hh"
#include "UniverseMaker.hh"
#include "ConfigMakerUtils.hh"
////////////////////////////////////////////////////////
//Categorys For Stacks 
/////////////////////////////////////////////////////////

const double DATA_POT = 4.54e+19;

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


const std::vector<NamedCategory<CCZeroPi_type>>
topology_categories = {
  NamedCategory<CCZeroPi_type>({kCC_0P_ZeroPi},          "CC0Pi"),
  NamedCategory<CCZeroPi_type>({kCC_1P_ZeroPi},          "CC1p0Pi"),
  NamedCategory<CCZeroPi_type>({kCC_2P_ZeroPi},          "CC2p0Pi"),
  NamedCategory<CCZeroPi_type>({kCC_3orGreaterP_ZeroPi}, "CC3>p0Pi")
};


const std::vector<NamedCategory<EventCategory>>
EventSelectionGroup_categories = {
NamedCategory<EventCategory>({kUnknown},        "Unknown"),
NamedCategory<EventCategory>({kSignalCCQE},     "Signal-(CCQE)"),
NamedCategory<EventCategory>({kSignalCCMEC},    "Signal-(CCMEC)"),
NamedCategory<EventCategory>({kSignalCCRES},    "Signal-(CCRES)"),
NamedCategory<EventCategory>({kSignalOther},    "Signal-(Other)"),
NamedCategory<EventCategory>({kNuMuCCNpi},      "#nu_{#mu}-CCN#pi"),
NamedCategory<EventCategory>({kNuMuCCOther},    "#nu_{#mu}-CCOther"),
NamedCategory<EventCategory>({kNuECC},          "#nu_{e}-CC"),
NamedCategory<EventCategory>({kNC},             "NC"),
NamedCategory<EventCategory>({kOOFV},           "Out FV"),
NamedCategory<EventCategory>({kOther},          "Other")
};


const std::vector<NamedCategory<FidVol_Category>>
ContainedGroup_categories_MC = {
NamedCategory<FidVol_Category>({kContained},     "Muon Contained"),
NamedCategory<FidVol_Category>({kUnContained},   "Muon UnContained"),
NamedCategory<FidVol_Category>({kContainmentBG},   "BG")
};

const std::vector<NamedCategory<Particle_type>>
ParticleGroup_reduced_categories = {
  NamedCategory<Particle_type>({k_MicroBooNECOSMICs},    "Cosmics"),
  NamedCategory<Particle_type>({kParticle_OTHER},        "Other"),
  NamedCategory<Particle_type>({kParticle_neutral},      "Neutral"),
  NamedCategory<Particle_type>({kKaon},                 "K^{#pm}"),
  NamedCategory<Particle_type>({kMuon},                   "#mu"),
  NamedCategory<Particle_type>({kPion_0_Electron_kGamma}, "e^{#pm}, #gamma, #pi^{0}"),
  NamedCategory<Particle_type>({kPion_pos_neg},          "#pi^{#pm}"),
  NamedCategory<Particle_type>({kProton},                "p")
};


const std::vector<NamedCategory<CCZeroPi_type>>
topology_categories_MC = {
  NamedCategory<CCZeroPi_type>({kCC_0P_ZeroPi},          "CC0#pi"),
  NamedCategory<CCZeroPi_type>({kCC_1P_ZeroPi},          "CC1p0#pi"),
  NamedCategory<CCZeroPi_type>({kCC_2P_ZeroPi},          "CC2p0#pi"),
  NamedCategory<CCZeroPi_type>({kCC_3orGreaterP_ZeroPi}, "CC3>p0#pi"),
  NamedCategory<CCZeroPi_type>({kCC_ZeroPi_BG}, "BG")
};





//////Functions//////////////////////////////
void reorderVector(std::vector<int>& vec, int x);
std::vector<int> MakereorderVector(std::vector<int> vec, int GotoFront);



void hist_maker_DATA(int input_Int) {

 std::string base_dir = "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_2_27_24_PmuCorrection/";
///

 //rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_D_E_all_reco2.rootå
 //rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_G1_beam_good_reco2_1e19.root
 //rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root
 std::string Rootfile1name; // = "rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";
 
 if(input_Int ==0 ){ Rootfile1name = "rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";}
 else if(input_Int ==1 ){ Rootfile1name = "rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root";}
  else if(input_Int ==2 ){ Rootfile1name = "rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_G1_beam_good_reco2_1e19.root";}
   else if(input_Int ==3 ){ Rootfile1name = "rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_D_E_all_reco2.root";}
   else if(input_Int ==4 ){ Rootfile1name = "rustv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_NuWro_reco2_reco2.root";}
   
  else{Rootfile1name = "rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";}


bool isFakeData = false; 
if(input_Int ==4)isFakeData = true; 

 std::string FullRootFileName = base_dir + Rootfile1name;
 
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

  //std::vector< double > p_edges = {0.1,0.2,0.275,0.335,0.46,0.61,0.785, 0.985, 1.185, 1.385, 1.585, 2.0  }; //TODO: Rough measurements based on the "resolution" plot.
  //std::vector< double > cos_edges = {-1, -0.85, -0.775, -0.7, -0.625, -0.55, -0.475, -0.4, -0.325,
  //   -0.25, -0.175, -0.1, -0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5,
  //    0.575, 0.65, 0.725, 0.8, 0.85, 0.875, 0.9, 0.925, 0.950, 0.975, 1.00};

  std::vector <double> cos_edges = GetCCZeroPi_Binning(kBINNING_Costheta);
  std::vector <double> p_edges = GetCCZeroPi_Binning(kBINNING_Pmu);


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
 // std::vector< double > Probability_range= {0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.00 };
  
  std::vector<double> Probability_range = generateBins();

  
  std::vector<double> trk_distance_v_edges = generateBins(4,15, 1 );
std::vector<double> trk_len_v_edges = {0,10,15,20,25,30,40,50,60,70,80,90,100};
   trk_len_v_edges.push_back(125);
   trk_len_v_edges.push_back(150);
   trk_len_v_edges.push_back(175);
   trk_len_v_edges.push_back(200);
   trk_len_v_edges.push_back(250);
  trk_len_v_edges.push_back(300);
  trk_len_v_edges.push_back(350);
  trk_len_v_edges.push_back(400);

  
  
std::vector< double > p_proton_edges = {0.25, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.2}; 
 std::vector< double > cos_edges_proton = {-1.0, -0.5, 0.0, 0.27, 0.45, 0.62, 0.76, 0.86,
      0.94, 1.0};


//std::vector< double > cos_edges2D = {-1.0, 0.0, 0.27, 0.45, 0.62, 0.65, 0.76, 0.86, 0.94, 1.00};
// std::vector< double > p_edges_2D = {0.1, 0.2, 0.275, 0.335, 0.46, 0.61, 0.785, 1.185, 1.585 };
std::vector< double > cos_edges2D = {-1.0, 0.0, 0.27, 0.5, 0.65, 0.8, 0.92, 1.00};
std::vector< double > p_edges_2D = {0.1, 0.24, 0.38, 0.48, 0.785, 1.185, 1.585, 2.0};

  
  
  
TH1D* h_Proton_mult  = new TH1D("h_Proton_mult", "h_Proton_mult", Mult.size() - 1, Mult.data()); 
TH1D* h_Proton_mult_Tracks  = new TH1D("h_Proton_mult_Tracks", "h_Proton_mult_Tracks", Mult.size() - 1, Mult.data()); 


TH1D* h_Else_mult_Tracks  = new TH1D("h_Else_mult_Tracks", "h_Else_mult_Tracks", Mult.size() - 1, Mult.data()); 

TH2D* h_Else_Score_Pmu  = new TH2D("h_Else_Score_Pmu", "", Probability_range.size() - 1, Probability_range.data(), p_edges.size() - 1, p_edges.data());
 TH2D*h_Else_Score_CosTheta  = new TH2D("h_Else_Score_CosTheta", "", Probability_range.size() - 1, Probability_range.data(), cos_edges.size() - 1, cos_edges.data());




TH1D* h_Pion_mult  = new TH1D("h_Pion_mult", "h_Pion_mult", Mult.size() - 1, Mult.data());
TH1D* h_NTracks  = new TH1D("h_NTracks", "h_NTracks", Mult.size() - 1, Mult.data());
  
TH1D* h_CCZeroP_MuonCandidate_BTDPrediction  = new TH1D("h_CCZeroP_MuonCandidate_BTDPrediction", 
"h_CCZeroP_MuonCandidate_BTDPrediction", Probability_range.size() - 1, Probability_range.data());


HistFolio<TH1D, Particle_type> h_CCZeroP_MuonCandidate_BTDPrediction_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_CCZeroP_MuonCandidate_BTDPrediction_Particle_categories", Probability_range, "h_CCZeroP_MuonCandidate_BTDPrediction_Particle_categories; Events");


HistFolio<TH1D, EventCategory> h_CCZeroP_MuonCandidate_BTDPrediction_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_CCZeroP_MuonCandidate_BTDPrediction_EventCategory", Probability_range, "h_CCZeroP_MuonCandidate_BTDPrediction_EventCategory; Events");






TH1D* h_CCZeroP_ProtonCandidate_BTDPrediction  = new TH1D("h_CCZeroP_ProtonCandidate_BTDPrediction",
"h_CCZeroP_ProtonCandidate_BTDPrediction", Probability_range.size() - 1, Probability_range.data());


HistFolio<TH1D, Particle_type> h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories", Probability_range, "h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories; Events");

HistFolio<TH1D, EventCategory> h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory", Probability_range, "h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory; Events");





TH1D* h_CCZeroP_Else_BTDPrediction  = new TH1D("h_CCZeroP_Else_BTDPrediction",
"h_CCZeroP_Else_BTDPrediction", Probability_range.size() - 1, Probability_range.data());

HistFolio<TH1D,  Particle_type> h_CCZeroP_Else_BTDPrediction_Particle_categories =
HistFolio<TH1D,Particle_type>(ParticleGroup_reduced_categories, "h_CCZeroP_Else_BTDPrediction_Particle_categories", Probability_range, "h_CCZeroP_Else_BTDPrediction_Particle_categories; Events");

HistFolio<TH1D, EventCategory> h_CCZeroP_Else_BTDPrediction_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_CCZeroP_Else_BTDPrediction_EventCategory", Probability_range, "hh_CCZeroP_Else_BTDPrediction_EventCategory; Events");





TH1D* h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1  = new TH1D("h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1", 
"h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_BTDGroup_Probability_tobeProton_trk1  = new TH1D("h_CCZeroP_BTDGroup_Probability_tobeProton_trk1", 
"h_CCZeroP_BTDGroup_Probability_tobeProton_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_BTDGroup_Probability_tobePion_trk1  = new TH1D("h_CCZeroP_BTDGroup_Probability_tobePion_trk1",
"h_CCZeroP_BTDGroup_Probability_tobePion_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_BTDGroup_Probability_tobeElse_trk1  = new TH1D("h_CCZeroP_BTDGroup_Probability_tobeElse_trk1",
"h_CCZeroP_BTDGroup_Probability_tobeElse_trk1", Probability_range.size() - 1, Probability_range.data());


TH1D* h_CC1P_BTDGroup_Probability_tobeMuon_tk1  = new TH1D("h_CC1P_BTDGroup_Probability_tobeMuon_tk1",
"h_CC1P_BTDGroup_Probability_tobeMuon_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobeProton_tk1  = new TH1D("h_CC1P_BTDGroup_Probability_tobeProton_tk1",
"h_CC1P_BTDGroup_Probability_tobeProton_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobePion_tk1  = new TH1D("h_CC1P_BTDGroup_Probability_tobePion_tk1",
"h_CC1P_BTDGroup_Probability_tobePion_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobeElse_tk1  = new TH1D("h_CC1P_BTDGroup_Probability_tobeElse_tk1",
"h_CC1P_BTDGroup_Probability_tobeElse_trk1", Probability_range.size() - 1, Probability_range.data());


TH1D* h_CC1P_BTDGroup_Probability_tobeMuon_tk2  = new TH1D("h_CC1P_BTDGroup_Probability_tobeMuon_tk2",
"h_CC1P_BTDGroup_Probability_tobeMuon_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobeProton_tk2  = new TH1D("h_CC1P_BTDGroup_Probability_tobeProton_tk2",
"h_CC1P_BTDGroup_Probability_tobeProton_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobePion_tk2  = new TH1D("h_CC1P_BTDGroup_Probability_tobePion_tk2",
"h_CC1P_BTDGroup_Probability_tobePion_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobeElse_tk2  = new TH1D("h_CC1P_BTDGroup_Probability_tobeElse_tk2",
"h_CC1P_BTDGroup_Probability_tobeElse_trk2", Probability_range.size() - 1, Probability_range.data());

//////////////////////////////
//cc2p
/////////////////////////////
TH1D* h_CC2P_BTDGroup_Probability_tobeMuon_tk1  = new TH1D("h_CC2P_BTDGroup_Probability_tobeMuon_tk1",
"h_CC2P_BTDGroup_Probability_tobeMuon_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeProton_tk1  = new TH1D("h_CC2P_BTDGroup_Probability_tobeProton_tk1",
"h_CC2P_BTDGroup_Probability_tobeProton_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobePion_tk1  = new TH1D("h_CC2P_BTDGroup_Probability_tobePion_tk1",
"h_CC2P_BTDGroup_Probability_tobePion_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeElse_tk1  = new TH1D("h_CC2P_BTDGroup_Probability_tobeElse_tk1",
"h_CC2P_BTDGroup_Probability_tobeElse_trk1", Probability_range.size() - 1, Probability_range.data());

TH1D* h_CC2P_BTDGroup_Probability_tobeMuon_tk2  = new TH1D("h_CC2P_BTDGroup_Probability_tobeMuon_tk2",
"h_CC2P_BTDGroup_Probability_tobeMuon_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeProton_tk2  = new TH1D("h_CC2P_BTDGroup_Probability_tobeProton_tk2",
"h_CC2P_BTDGroup_Probability_tobeProton_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobePion_tk2  = new TH1D("h_CC2P_BTDGroup_Probability_tobePion_tk2",
"h_CC2P_BTDGroup_Probability_tobePion_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeElse_tk2  = new TH1D("h_CC2P_BTDGroup_Probability_tobeElse_tk2",
"h_CC2P_BTDGroup_Probability_tobeElse_trk2", Probability_range.size() - 1, Probability_range.data());


TH1D* h_CC2P_BTDGroup_Probability_tobeMuon_tk3  = new TH1D("h_CC2P_BTDGroup_Probability_tobeMuon_tk3",
"h_CC2P_BTDGroup_Probability_tobeMuon_trk3", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeProton_tk3  = new TH1D("h_CC2P_BTDGroup_Probability_tobeProton_tk3",
"h_CC2P_BTDGroup_Probability_tobeProton_trk3", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobePion_tk3  = new TH1D("h_CC2P_BTDGroup_Probability_tobePion_tk3",
"h_CC2P_BTDGroup_Probability_tobePion_trk3", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeElse_tk3  = new TH1D("h_CC2P_BTDGroup_Probability_tobeElse_tk3",
"h_CC2P_BTDGroup_Probability_tobeElse_trk3", Probability_range.size() - 1, Probability_range.data());

  /////////////////////////////////ProtonMupl/////////////////////////////////////////////////
  TH1D* h_p             = new TH1D("h_p", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_Correction             = new TH1D("h_p_Correction", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  TH1D* h_leadingProton_p             = new TH1D("h_leadingProton_p", ";p_{p}^{reco}", p_proton_edges.size() - 1, p_proton_edges.data());
     h_leadingProton_p->SetDirectory(0); h_leadingProton_p->Sumw2();
  TH1D* h_costheta_proton             = new TH1D("h_costheta_proton", ";p_{#mu}^{reco}", cos_edges_proton.size() - 1, cos_edges_proton.data());
   h_costheta_proton->SetDirectory(0); h_costheta_proton->Sumw2();
  
 HistFolio<TH1D, BDT_Category> h_p_BTDGroup_categories =
 HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_p_BTDGroup_categories", p_edges, "h_p_BTDGroup_categories; Events");
 
 HistFolio<TH1D, FidVol_Category> h_p_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_p_containment_categories", p_edges, "h_p_containment_categories");

HistFolio<TH1D, CCZeroPi_type> h_p_topology_categories =
  HistFolio<TH1D, CCZeroPi_type>(topology_categories, "h_p_topogoly_categories", p_edges, "h_p_topogoly_categories");


HistFolio<TH1D, BDT_Category> h_NTracks_BTDGroup_categories =
HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_NTracks_BTDGroup_categories", Mult, "h_NTracks_EventCategory; Events");


h_p->SetDirectory(0); h_p->Sumw2();
h_p_Correction->SetDirectory(0); h_p_Correction->Sumw2();
 ///////////////////////////////////////
 //CC0P0pi
 /////////////////////////////////////
    TH1D* h_p_CC0P0Pi_Only             = new TH1D("h_p_CC0P0Pi_Only", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
    TH1D* h_p_CCNP0Pi_Only             = new TH1D("h_p_CCNP0Pi_Only", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  h_p_CC0P0Pi_Only->SetDirectory(0); h_p_CC0P0Pi_Only->Sumw2();
  h_p_CCNP0Pi_Only->SetDirectory(0); h_p_CCNP0Pi_Only->Sumw2();
    ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_costheta_p             = new TH1D("h_costheta_p", ";p_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data());
   h_costheta_p->SetDirectory(0); h_costheta_p->Sumw2();
   
  TH1D* h_costheta_p_CC0P0Pi_Only = new TH1D("h_costheta_p_CC0P0Pi_Only",
  ";p_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data());
   h_costheta_p_CC0P0Pi_Only->SetDirectory(0); h_costheta_p_CC0P0Pi_Only->Sumw2();
   
  TH1D* h_costheta_p_CCNP0Pi_Only = new TH1D("h_costheta_p_CCNP0Pi_Only",
  ";p_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data());
   h_costheta_p_CCNP0Pi_Only->SetDirectory(0); h_costheta_p_CCNP0Pi_Only->Sumw2();
  
  
  
 HistFolio<TH1D, BDT_Category> h_costheta_p_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_costheta_p_BTDGroup_categories", cos_edges, "h_costheta_p_BTDGroup_categories; Events");
    
  HistFolio<TH1D, FidVol_Category> h_costheta_p_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_costheta_p_containment_categories", cos_edges, "h_costheta_p_containment_categories");
  /////////////////////////////////////////////////////////////////////////    
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_pn             = new TH1D("h_pn", ";pn_{#mu}^{reco}", pn_edges.size() - 1, pn_edges.data());
  h_pn->SetDirectory(0);            h_pn->Sumw2();
  
  HistFolio<TH1D, BDT_Category> h_pn_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories,
  "h_pn_BTDGroup_categories", pn_edges, "h_pn_EventCategory; Events");
   
  HistFolio<TH1D, FidVol_Category> h_pn_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, 
  "h_pn_containment_categories", pn_edges, "h_pn_containment_categories");
    
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_alphaT             = new TH1D("h_delta_alphaT", ";#delta_{#alpha T}^{reco}", delta_alphaT_edges.size() - 1, delta_alphaT_edges.data());
  
  h_delta_alphaT->SetDirectory(0);            h_delta_alphaT->Sumw2();
  
  HistFolio<TH1D, BDT_Category> h_delta_alphaT_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories,
  "h_delta_alphaT_BTDGroup_categories", delta_alphaT_edges, "h_delta_alphaT_EventCategory; Events");
  
  HistFolio<TH1D, FidVol_Category> h_delta_alphaT_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories,
  "h_delta_alphaT_containment_categories", delta_alphaT_edges, "h_delta_alphaT_containment_categories");
    
  
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_pTx             = new TH1D("h_delta_pTx", ";#delta_{#alpha Tx}^{reco}", delta_pTx_edges.size() - 1, delta_pTx_edges.data());
  h_delta_pTx->SetDirectory(0);            h_delta_pTx->Sumw2();
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_pTy             = new TH1D("h_delta_pTy",
  ";#delta_{#alpha T}^{reco}", delta_pTy_edges.size() - 1, delta_pTy_edges.data());
  h_delta_pTy->SetDirectory(0);            h_delta_pTy->Sumw2();
  

  HistFolio<TH1D, FidVol_Category> h_delta_pTy_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, 
  "h_delta_pTy_containment_categories", delta_pTy_edges, "h_delta_pTy_containment_categories");

////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_phiT             = new TH1D("h_delta_phiT", "", delta_phiT_edges.size() - 1, delta_phiT_edges.data());
  h_delta_phiT->SetDirectory(0);            h_delta_phiT->Sumw2();
   
 
  HistFolio<TH1D, BDT_Category> h_delta_phiT_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories,
  "h_delta_phiT_BTDGroup_categories", delta_phiT_edges, "h_delta_phiT_EventCategory; Events");
  
  HistFolio<TH1D, FidVol_Category> h_delta_phiT_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories,
  "h_delta_phiT_containment_categories", delta_phiT_edges, "h_delta_phiT_containment_categories");
  
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_openingAngle             = new TH1D("h_openingAngle", "", openingAngle_edges.size() - 1, openingAngle_edges.data());
  h_openingAngle->SetDirectory(0);            h_openingAngle->Sumw2();
   
  HistFolio<TH1D, FidVol_Category> h_openingAngle_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories,
  "h_openingAngle_containment_categories", openingAngle_edges, "h_openingAngle_containment_categories");
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_openingAngle_deg             = new TH1D("h_openingAngle_deg", "",
  openingAngleDeg_edges.size() - 1, openingAngleDeg_edges.data());
  h_openingAngle_deg->SetDirectory(0);            h_openingAngle_deg->Sumw2();
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexX             = new TH1D("h_VertexX", "", VertexX_edges.size() - 1, VertexX_edges.data());
  h_VertexX->SetDirectory(0);            h_VertexX->Sumw2();
    
   HistFolio<TH1D, BDT_Category> h_VertexX_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories,
  "h_VertexX_BTDGroup_categories", VertexX_edges, "h_VertexX_EventCategory; Events");
   HistFolio<TH1D, FidVol_Category> h_VertexX_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories,
  "h_VertexX_containment_categories", VertexX_edges, "h_VertexX_containment_categories");
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexY             = new TH1D("h_VertexY", "", VertexY_edges.size() - 1, VertexY_edges.data());
  h_VertexY->SetDirectory(0);            h_VertexY->Sumw2();
    
  HistFolio<TH1D, BDT_Category> h_VertexY_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories,
  "h_VertexY_BTDGroup_categories", VertexY_edges, "h_VertexY_EventCategory; Events");
  
  HistFolio<TH1D, FidVol_Category> h_VertexY_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories,
  "h_VertexY_containment_categories", VertexY_edges, "h_VertexY_containment_categories");
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexZ             = new TH1D("h_VertexZ", "", VertexZ_edges.size() - 1, VertexZ_edges.data());
  
  h_VertexZ->SetDirectory(0);            h_VertexZ->Sumw2();
  
  HistFolio<TH1D, BDT_Category> h_VertexZ_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories,
  "h_VertexZ_BTDGroup_categories", VertexY_edges, "h_VertexZ_EventCategory; Events");
   
  HistFolio<TH1D, FidVol_Category> h_VertexZ_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, 
  "h_VertexZ_containment_categories", VertexY_edges, "h_VertexZ_containment_categories");
//////////////////////////////////////

TH1D* h_topological_score  = new TH1D("h_topological_score", "h_topological_score", Probability_range.size() - 1, Probability_range.data());
h_topological_score->SetDirectory(0);h_topological_score->Sumw2();

TH1D* h_trk_score_v  = new TH1D("h_trk_score_v", "h_trk_score_v", Probability_range.size() - 1, Probability_range.data());
h_trk_score_v->SetDirectory(0);h_trk_score_v->Sumw2();
  
  TH1D* h_trk_len_v  = new TH1D("h_trk_len_v", "h_trk_len_v", trk_len_v_edges.size() - 1, trk_len_v_edges.data());
  h_trk_len_v->SetDirectory(0);h_trk_len_v->Sumw2();
  
  
  TH1D* h_trk_distance_v  = new TH1D("h_trk_distance_v", "h_trk_distance_v", trk_distance_v_edges.size() - 1, trk_distance_v_edges.data());
h_trk_distance_v->SetDirectory(0);h_trk_distance_v->Sumw2();
  
  ////////////////////////////
  //2D vectex
  ///////////////////////////
  TH2D* h_VertexX_Z  = new TH2D("h_VertexX_Z",
  "", VertexZ_edges.size() - 1, VertexZ_edges.data(), VertexX_edges.size() - 1, VertexX_edges.data());
  h_VertexX_Z->SetDirectory(0);
  h_VertexX_Z->Sumw2();
  TH2D* h_VertexY_Z  = new TH2D("h_VertexY_Z", 
  "", VertexZ_edges.size() - 1, VertexZ_edges.data(), VertexY_edges.size() - 1, VertexY_edges.data());
  h_VertexY_Z->SetDirectory(0);
  h_VertexY_Z->Sumw2();
  TH2D* h_VertexX_Y  = new TH2D("h_VertexX_Y", 
  "", VertexX_edges.size() - 1, VertexX_edges.data(), VertexY_edges.size() - 1, VertexY_edges.data());
   h_VertexX_Y->SetDirectory(0);
   h_VertexX_Y->Sumw2();

  TH2D* h_costheta_Pmu  = new TH2D("h_costheta_Pmu", 
  "h_costheta_Pmu", 
  cos_edges2D.size() - 1, cos_edges2D.data(), p_edges_2D.size() - 1, p_edges_2D.data());

   h_costheta_Pmu->SetDirectory(0);
   h_costheta_Pmu->Sumw2();

  double res_max = -0.01;


  std::map<int , BinMap> TH1Poly_binMap_test; 
  
   UBTH2Poly instance;
   gInterpreter->Declare("#include \"includes/UBTH2Poly.h\""); 
  
  
  UBTH2Poly *UBTH2Poly_costheta_Pmu = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "TH2Poly_costheta_Pmu_TRUE");
  UBTH2Poly* h_costheta_Pmu_UBTH2Poly  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly", "h_costheta_Pmu_UBTH2Poly",  MUON_2D_BIN_EDGES, true);
  UBTH2Poly* h_costheta_Pmu_UBTH2Poly_inclusive  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly_inclusive", "h_costheta_Pmu_UBTH2Poly_inclusive",  MUON_2D_BIN_EDGES_inclusive, false);





  for (auto& infile : filenames) {
    std::cout << "File: " << infile << std::endl;

    TFile* f = TFile::Open(infile.c_str());
    
    
    TParameter<Float_t>* myParameter = dynamic_cast<TParameter<Float_t>*>(f->Get("summed_pot"));
    double POT_read = myParameter->GetVal();
    
    if (!(f && f->IsOpen() && !f->IsZombie())) {
        std::cout << "Bad file! " << infile << std::endl;
        continue;
    }
    
    TTree* stv  = (TTree*) f->Get("stv_tree");
    if (!(stv && stv->GetEntries())) {
       std::cout << "Bad file! " << infile << std::endl;
       continue;
    }



   std::cout<<"read POT =  "<< POT_read << std::endl;

  double FakeDataPOTscale = 1.0; 
   //if(isFakeData==true){
   //
   //
   //}

  FakeDataPOTscale = (DATA_POT / POT_read);

  bool mc_signal, sel_signal,sel_signal_CCNP0Pi, sel_muon_contained, sel_BDT_NotBOGUS;
  TVector3 *p3_mu=nullptr, *mc_p3_mu=nullptr;
  TVector3 *p3_lead_p=nullptr, *mc_p3_lead_p=nullptr;
  

  
  std::vector<int> *mc_pdg = nullptr;
  std::vector<int> *pfpdg = nullptr;
  std::vector<int> *backtracked_pdg = nullptr;
  
  
  std::vector<int> *trk_BDT_PID = nullptr;
  std::vector<int> *trk_BDT_Probility = nullptr;
  std::vector<float> *trk_score_v = nullptr;
  
  std::vector<float> *trk_len_v = nullptr;
  std::vector<float> *trk_distance_v = nullptr;
  float spline_weight, tuned_cv_weight, topological_score;
  float pn, delta_alphaT, delta_pTx, delta_pTy, delta_pL, delta_phiT, delta_pT;
  float reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z;
  int proton_multiplicity,pion_multiplicity; 
  int index_muon, index_leadproton,category; 

  float BTD_Muon_Candidate_prediction;
  bool sel_Muon_BDT_SCORE_above_threshold, sel_All_Protons_BDT_SCORE_above_threshold;
  bool sel_BDT_predicts_1plusTrks_tobeProtons;
  bool sel_cosmic_ip_cut_passed;
  std::vector<float> *trk_BDT_score_0 = nullptr;
  std::vector<float> *trk_BDT_score_1 = nullptr;
  std::vector<float> *trk_BDT_score_2 = nullptr;
  std::vector<float> *trk_BDT_score_3 = nullptr;



bool sel_reco_vertex_in_FV , sel_pfp_starts_in_PCV ;
bool sel_has_muon_candidate,  sel_topo_cut_passed  ;
bool sel_no_reco_showers , sel_muon_passed_mom_cuts  ;
bool sel_pions_above_threshold , sel_has_pi_candidate ;


//int pion_multiplicity; 

  stv->SetBranchStatus("*", 0);
  stv->SetBranchStatus("sel_CC0pi", 1);
  stv->SetBranchStatus("sel_CCNp0pi", 1);  
  stv->SetBranchStatus("p3_mu", 1);
  stv->SetBranchStatus("p3_lead_p", 1);
  stv->SetBranchStatus("sel_muon_contained", 1);
  stv->SetBranchStatus("spline_weight", 1);
  stv->SetBranchStatus("tuned_cv_weight", 1);
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



  stv->SetBranchStatus("proton_multiplicity",1); 
  stv->SetBranchStatus("pion_multiplicity",1); 
  stv->SetBranchStatus("trk_BDT_PID",1); 
  //stv->SetBranchStatus("n_pfps",1); 
  stv->SetBranchStatus("pfpdg",1);
  stv->SetBranchStatus("backtracked_pdg",1);
  stv->SetBranchStatus("trk_BDT_score_0",1);
  stv->SetBranchStatus("trk_BDT_score_1",1);
  stv->SetBranchStatus("trk_BDT_score_2",1);
  stv->SetBranchStatus("trk_BDT_score_3",1);

stv->SetBranchStatus("topological_score",1);
stv->SetBranchAddress("topological_score", &topological_score);

stv->SetBranchStatus("trk_score_v",1);
stv->SetBranchAddress("trk_score_v", &trk_score_v);

  stv->SetBranchStatus("trk_len_v",1); 
  stv->SetBranchAddress("trk_len_v", &trk_len_v);
  
  stv->SetBranchStatus("trk_distance_v",1); 
  stv->SetBranchAddress("trk_distance_v", &trk_distance_v);


  stv->SetBranchAddress("sel_CC0pi", &sel_signal);
  stv->SetBranchAddress("sel_CCNp0pi", &sel_signal_CCNP0Pi);
  
  
  stv->SetBranchAddress("p3_mu", &p3_mu);

  stv->SetBranchAddress("p3_lead_p", &p3_lead_p);

  stv->SetBranchAddress("sel_muon_contained", &sel_muon_contained);
  stv->SetBranchAddress("spline_weight", &spline_weight);
  stv->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
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
  stv->SetBranchAddress("proton_multiplicity", &proton_multiplicity);
  stv->SetBranchAddress("pion_multiplicity", &pion_multiplicity);
  stv->SetBranchAddress("trk_BDT_PID", &trk_BDT_PID);
  //stv->SetBranchAddress("n_pfps", &n_pfps);
  stv->SetBranchAddress("pfpdg", &pfpdg);
  stv->SetBranchAddress("backtracked_pdg", &backtracked_pdg);

 stv->SetBranchAddress("trk_BDT_score_0", &trk_BDT_score_0);
 stv->SetBranchAddress("trk_BDT_score_1", &trk_BDT_score_1);
 stv->SetBranchAddress("trk_BDT_score_2", &trk_BDT_score_2);
 stv->SetBranchAddress("trk_BDT_score_3", &trk_BDT_score_3);

stv->SetBranchStatus("category",1);
stv->SetBranchAddress("category", &category);
  
  stv->SetBranchStatus("sel_Muon_BDT_SCORE_above_threshold", 1);
  stv->SetBranchStatus("sel_All_Protons_BDT_SCORE_above_threshold", 1);
  stv->SetBranchStatus("sel_BDT_predicts_1plusTrks_tobeProtons", 1);
  
  
  stv->SetBranchAddress("sel_Muon_BDT_SCORE_above_threshold", &sel_Muon_BDT_SCORE_above_threshold);
  stv->SetBranchAddress("sel_All_Protons_BDT_SCORE_above_threshold", &sel_All_Protons_BDT_SCORE_above_threshold);
  
  stv->SetBranchAddress("sel_BDT_predicts_1plusTrks_tobeProtons", &sel_BDT_predicts_1plusTrks_tobeProtons);
  
  stv->SetBranchStatus("sel_BDT_NotBOGUS",1); 
  stv->SetBranchAddress("sel_BDT_NotBOGUS", &sel_BDT_NotBOGUS);
  
  
  stv->SetBranchStatus("BTD_Muon_Candidate_prediction", 1);
  stv->SetBranchAddress("BTD_Muon_Candidate_prediction", &BTD_Muon_Candidate_prediction);
  
  
  stv->SetBranchStatus("muon_candidate_idx", 1);
  stv->SetBranchStatus("lead_p_candidate_idx", 1);
  stv->SetBranchAddress("muon_candidate_idx", &index_muon);
  stv->SetBranchAddress("lead_p_candidate_idx", &index_leadproton);


  stv->SetBranchStatus("sel_cosmic_ip_cut_passed", 1);
 stv->SetBranchAddress("sel_cosmic_ip_cut_passed", &sel_cosmic_ip_cut_passed);
  

stv->SetBranchStatus("sel_reco_vertex_in_FV",1); 
stv->SetBranchAddress("sel_reco_vertex_in_FV",&sel_reco_vertex_in_FV);
stv->SetBranchStatus("sel_pfp_starts_in_PCV",1); 
stv->SetBranchAddress("sel_pfp_starts_in_PCV",&sel_pfp_starts_in_PCV);
stv->SetBranchStatus("sel_has_muon_candidate",1); 
stv->SetBranchAddress("sel_has_muon_candidate",&sel_has_muon_candidate);
stv->SetBranchStatus("sel_topo_cut_passed ",1); 
stv->SetBranchAddress("sel_topo_cut_passed",&sel_topo_cut_passed);
stv->SetBranchStatus("sel_no_reco_showers",1); 
stv->SetBranchAddress("sel_no_reco_showers",&sel_no_reco_showers);
stv->SetBranchStatus("sel_muon_passed_mom_cuts",1); 
stv->SetBranchAddress("sel_muon_passed_mom_cuts",&sel_muon_passed_mom_cuts);
stv->SetBranchStatus("sel_pions_above_threshold",1); 
stv->SetBranchAddress("sel_pions_above_threshold",&sel_pions_above_threshold);
stv->SetBranchStatus("sel_has_pi_candidate",1); 
stv->SetBranchAddress("sel_has_pi_candidate",&sel_has_pi_candidate);
//stv->SetBranchStatus("pion_multiplicity",1); 
//stv->SetBranchAddress("pion_multiplicity",&pion_multiplicity);


  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout<< "~~~~~~~~~~~~~~~ STarting Loop~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout << " The Number of Entries "<<  stv->GetEntries() << std::endl;
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;


  for (long i=0; i<stv->GetEntries(); i++)
  {
    stv->GetEntry(i);
   // if (sel_signal && mc_signal){
     if ( i % 10000 == 0 ) {
      std::cout << "Processing event #" << i << '\n';
    }
      if (sel_signal && sel_BDT_NotBOGUS ){


//sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV &&
//sel_has_muon_candidate && sel_topo_cut_passed && 
//sel_no_reco_showers && sel_muon_passed_mom_cuts && 
//sel_pions_above_threshold && sel_has_pi_candidate && pion_multiplicity > 0



      //&& sel_BDT_predicts_1plusTrks_tobeProtons
   /* 
   
        bool ProtonGood= true; s
        
     for(int index1 = 0 ; index1 <pfpdg->size(); index1++ ){
      double score_3 = trk_BDT_score_3->at(index1);
     double score_1 = trk_BDT_score_1->at(index1);
      if (index1 != index_muon){
     
     if (score_3 < .28 )
     ProtonGood = false ;
     }     
     }
     
     if(ProtonGood==false) continue; 
   */
//std::cout<<"pn = " << pn << std::endl;

FidVol_Category Containmenttype = EventCategory_tool.GetContainmentType_ForData(sel_muon_contained);
CCZeroPi_type TopologyType = EventCategory_tool.returnCCZeroPi_type_ForData(proton_multiplicity);
EventCategory TrueEventType = EventCategory_tool.returnEventCategoryType(category);


//double Wgt = 1.0;
double Wgt = spline_weight * tuned_cv_weight * FakeDataPOTscale;
//double Wgt = FakeDataPOTscale;
double Pmu = p3_mu->Mag();
//std::cout<<"Wgt"<<Wgt<< std::endl;


double CosTheta = p3_mu->CosTheta();
double P_proton = p3_lead_p->Mag();
double P_CosTheta = p3_lead_p->CosTheta();

double OpenAngle = TMath::ACos( (p3_mu->X()*p3_lead_p->X() + p3_mu->Y()*p3_lead_p->Y() + p3_mu->Z()*p3_lead_p->Z()) / p3_mu->Mag());
double OpenAngle_deg = OpenAngle * TMath::RadToDeg();

int NTracks = proton_multiplicity + 1; // assuming the number of tracks is the Proton Multipleicity + muon track 


h_p->Fill(Pmu,Wgt); 
h_p_Correction->Fill(Pmu,Wgt); 


h_Proton_mult->Fill(proton_multiplicity,Wgt); 
h_Pion_mult->Fill(pion_multiplicity,Wgt); 


h_costheta_p->Fill(CosTheta,Wgt);

h_p_containment_categories.GetComponentHist(Containmenttype)->Fill(Pmu,Wgt);
h_p_topology_categories.GetComponentHist(TopologyType)->Fill(Pmu,Wgt);
h_costheta_p_containment_categories.GetComponentHist(Containmenttype)->Fill(CosTheta,Wgt);

h_costheta_Pmu ->Fill(CosTheta,Pmu,Wgt); 


h_topological_score->Fill(topological_score,Wgt);
 
UBTH2Poly_costheta_Pmu->Fill(CosTheta,Pmu,Wgt); 
h_costheta_Pmu_UBTH2Poly->Fill(CosTheta,Pmu,Wgt);
h_costheta_Pmu_UBTH2Poly_inclusive->Fill(CosTheta,Pmu,Wgt); 
 
 
 
 if(proton_multiplicity ==0){
 h_p_CC0P0Pi_Only->Fill(Pmu,Wgt);
 h_costheta_p_CC0P0Pi_Only->Fill(CosTheta,Wgt);

 
 }
 else{ 
  h_p_CCNP0Pi_Only->Fill(Pmu,Wgt);

 h_costheta_p_CCNP0Pi_Only->Fill(CosTheta,Wgt);
 
 }
 
 
 
 

 
 //std::cout<<"~~~~Before"<<std::endl;
 // for (int i : *trk_BDT_PID) {
 //       std::cout << i << " ";
 //   }
 //   
 //   std::cout<<"~~~~AFter~~"<<std::endl;

std::vector<int> BDTPID_Prediction_Muon_front = MakereorderVector(*trk_BDT_PID, BDT_Category::kBDT_Muon);

  //for (int i : BDTPID_Prediction_Muon_front) {
  //      std::cout << i << " ";
  //  }
int TrackN = 0 ; 
int TrackN_proton = 1; 

for(auto PID : BDTPID_Prediction_Muon_front){


BDT_Category BDT_Type = EventCategory_tool.returnBDTPredictionType(PID);

h_NTracks_BTDGroup_categories.GetComponentHist(BDT_Type)->Fill(TrackN,Wgt);
h_NTracks->Fill(TrackN,Wgt);



 if(PID==BDT_Category::kBDT_Proton){
  h_Proton_mult_Tracks->Fill(TrackN_proton,Wgt);
  TrackN_proton++;
 }
else if (PID==BDT_Category::kBDT_Else){

h_Else_mult_Tracks->Fill(TrackN ,Wgt);

}


TrackN++;
}



for(int index1 = 0 ; index1 <pfpdg->size(); index1++ )
{

    int true_pdg = backtracked_pdg->at(index1);
    Particle_type REduced_Particle =  EventCategory_tool.GetParticlegroup_typeReduced(true_pdg);

 double score_0 = trk_BDT_score_0->at(index1);
 double score_1 = trk_BDT_score_1->at(index1);
 double score_2 = trk_BDT_score_2->at(index1);
 double score_3 = trk_BDT_score_3->at(index1);
int PIDBDT = trk_BDT_PID->at(index1);




 h_trk_score_v->Fill(trk_score_v->at(index1),Wgt);


h_trk_distance_v->Fill(trk_distance_v->at(index1),Wgt);


if(index1 == index_muon){

// Lead Muon Prediction
h_CCZeroP_MuonCandidate_BTDPrediction->Fill(score_1,Wgt);
h_trk_len_v->Fill(trk_len_v->at(index1),Wgt);
}



else{

if(PIDBDT==BDT_Category::kBDT_Proton){
h_CCZeroP_ProtonCandidate_BTDPrediction->Fill(score_3,Wgt);
 h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories.GetComponentHist(REduced_Particle)->Fill(score_3,Wgt);
h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory.GetComponentHist(TrueEventType)->Fill(score_3,Wgt);
}

else if(PIDBDT==BDT_Category::kBDT_Else){

h_CCZeroP_Else_BTDPrediction->Fill(score_0,Wgt);
h_CCZeroP_Else_BTDPrediction_Particle_categories.GetComponentHist(REduced_Particle)->Fill(score_0,Wgt);
h_CCZeroP_Else_BTDPrediction_EventCategory.GetComponentHist(TrueEventType)->Fill(score_0,Wgt);

}

}

 if(index1 == index_leadproton){

//h_CCZeroP_ProtonCandidate_BTDPrediction->Fill(score_3,Wgt);


}



 if(proton_multiplicity ==0){
 
 h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1->Fill(score_1,Wgt);
 h_CCZeroP_BTDGroup_Probability_tobeProton_trk1->Fill(score_3,Wgt);
 h_CCZeroP_BTDGroup_Probability_tobePion_trk1->Fill(score_2,Wgt);
 h_CCZeroP_BTDGroup_Probability_tobeElse_trk1->Fill(score_0,Wgt);
 }

 else if (proton_multiplicity ==1){

 if(index1==1){
 h_CC1P_BTDGroup_Probability_tobeMuon_tk1->Fill(score_1,Wgt);
 h_CC1P_BTDGroup_Probability_tobeProton_tk1->Fill(score_3,Wgt);
 h_CC1P_BTDGroup_Probability_tobePion_tk1->Fill(score_2,Wgt);
 h_CC1P_BTDGroup_Probability_tobeElse_tk1->Fill(score_0,Wgt);
 }
 
 if(index1==2){
 h_CC1P_BTDGroup_Probability_tobeMuon_tk2->Fill(score_1,Wgt);
 h_CC1P_BTDGroup_Probability_tobeProton_tk2->Fill(score_3,Wgt);
 h_CC1P_BTDGroup_Probability_tobePion_tk2->Fill(score_2,Wgt);
 h_CC1P_BTDGroup_Probability_tobeElse_tk2->Fill(score_0,Wgt);
 }

 }// ENd of  1 proton
 
 else if (proton_multiplicity ==2)
 {
  if(index1==1){
  h_CC2P_BTDGroup_Probability_tobeMuon_tk1->Fill(score_1,Wgt);
  h_CC2P_BTDGroup_Probability_tobeProton_tk1->Fill(score_3,Wgt);
  h_CC2P_BTDGroup_Probability_tobePion_tk1->Fill(score_2,Wgt);
  h_CC2P_BTDGroup_Probability_tobeElse_tk1->Fill(score_0,Wgt);
  }

 else if(index1==2){
 h_CC2P_BTDGroup_Probability_tobeMuon_tk2->Fill(score_1,Wgt);
 h_CC2P_BTDGroup_Probability_tobeProton_tk2->Fill(score_3,Wgt);
 h_CC2P_BTDGroup_Probability_tobePion_tk2->Fill(score_2,Wgt);
 h_CC2P_BTDGroup_Probability_tobeElse_tk2->Fill(score_0,Wgt);
 }

 else if(index1==3){
 h_CC2P_BTDGroup_Probability_tobeMuon_tk3->Fill(score_1,Wgt);
 h_CC2P_BTDGroup_Probability_tobeProton_tk3->Fill(score_3,Wgt);
 h_CC2P_BTDGroup_Probability_tobePion_tk3->Fill(score_2,Wgt);
 h_CC2P_BTDGroup_Probability_tobeElse_tk3->Fill(score_0,Wgt);
 }
}// ENd of  2 proton

}// End of For Loop 


h_VertexX->Fill(reco_nu_vtx_sce_x,Wgt);  
h_VertexY->Fill(reco_nu_vtx_sce_y,Wgt);  
h_VertexZ->Fill(reco_nu_vtx_sce_z,Wgt);

h_VertexX_Z->Fill(reco_nu_vtx_sce_z,reco_nu_vtx_sce_x,Wgt); 
h_VertexY_Z->Fill(reco_nu_vtx_sce_z,reco_nu_vtx_sce_y,Wgt); 
h_VertexX_Y->Fill(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,Wgt); 

if(proton_multiplicity !=0){

h_openingAngle->Fill(OpenAngle,Wgt);  
h_openingAngle_deg->Fill(OpenAngle_deg,Wgt);  
h_pn->Fill(pn,Wgt); 
h_delta_alphaT->Fill(delta_alphaT,Wgt); 
h_delta_pTx->Fill(delta_pTx,Wgt);  
h_delta_pTy->Fill(delta_pTy,Wgt);  
h_delta_phiT->Fill(delta_phiT,Wgt);  
h_leadingProton_p->Fill(P_proton,Wgt);
h_costheta_proton->Fill(P_CosTheta,Wgt);

} // End if 1 Proton 





}// ENd of If Statement, that passed Signal Selection 
//////////////////////////////////////////////////////  
}// End of Event Loop
///////////////////////////////////////////////////////


  f->Close();

}


/////////////////////////
// Make Root File to look 
/////////////////////////
std::string OutPulRootName = "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/Data/DATA_Selected_newCC0pi_" + Rootfile1name;
//_addthresholds_WithProtonPIDCUT
std::cout<< "Printing Root File named: "<< OutPulRootName << std::endl;
auto outFile = TFile::Open(OutPulRootName.c_str(), "RECREATE");
outFile->cd(); 

h_p->Write();
h_p_Correction->Write();
h_costheta_p->Write();

h_costheta_Pmu->Write();
h_pn->Write();
h_Proton_mult->Write();
h_Proton_mult_Tracks->Write();
h_Pion_mult->Write();

h_p_containment_categories.WriteToFile(*outFile);
h_p_topology_categories.WriteToFile(*outFile);


h_NTracks_BTDGroup_categories.WriteToFile(*outFile);
h_delta_alphaT->Write();
h_NTracks->Write();

h_delta_pTx->Write();

h_delta_pTy->Write();

h_delta_phiT->Write();

h_openingAngle->Write();
h_openingAngle_deg->Write();


h_VertexX->Write();

h_VertexY->Write();

h_VertexZ->Write();

h_VertexX_Z->Write();
h_VertexY_Z->Write();
h_VertexX_Y->Write();



h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1->Write();
h_CCZeroP_BTDGroup_Probability_tobeProton_trk1->Write();
h_CCZeroP_BTDGroup_Probability_tobePion_trk1->Write();
h_CC1P_BTDGroup_Probability_tobeMuon_tk1->Write();
h_CC1P_BTDGroup_Probability_tobeProton_tk1->Write();
h_CC1P_BTDGroup_Probability_tobePion_tk1->Write();
h_CC1P_BTDGroup_Probability_tobeMuon_tk2->Write();
h_CC1P_BTDGroup_Probability_tobeProton_tk2->Write();
h_CC1P_BTDGroup_Probability_tobePion_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeMuon_tk1->Write();
h_CC2P_BTDGroup_Probability_tobeProton_tk1->Write();
h_CC2P_BTDGroup_Probability_tobePion_tk1->Write();
h_CC2P_BTDGroup_Probability_tobeMuon_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeProton_tk2->Write();
h_CC2P_BTDGroup_Probability_tobePion_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeMuon_tk3->Write();
h_CC2P_BTDGroup_Probability_tobeProton_tk3->Write();
h_CC2P_BTDGroup_Probability_tobePion_tk3->Write();


h_CCZeroP_BTDGroup_Probability_tobeElse_trk1->Write();
h_CC1P_BTDGroup_Probability_tobeElse_tk1->Write();
h_CC1P_BTDGroup_Probability_tobeElse_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeElse_tk1->Write();
h_CC2P_BTDGroup_Probability_tobeElse_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeElse_tk3->Write();

h_CCZeroP_MuonCandidate_BTDPrediction->Write();
h_CCZeroP_ProtonCandidate_BTDPrediction->Write();
h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories.WriteToFile(*outFile);
h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory.WriteToFile(*outFile);
h_CCZeroP_Else_BTDPrediction_Particle_categories.WriteToFile(*outFile);
h_CCZeroP_Else_BTDPrediction_EventCategory.WriteToFile(*outFile);


h_CCZeroP_Else_BTDPrediction->Write();
h_p_CCNP0Pi_Only->Write();
h_p_CC0P0Pi_Only->Write();

h_costheta_p_CC0P0Pi_Only->Write();
h_costheta_p_CCNP0Pi_Only->Write();


h_leadingProton_p->Write();
h_costheta_proton->Write();


h_topological_score->Write();
h_trk_score_v->Write();

h_trk_len_v->Write();
h_trk_distance_v->Write();

char HistName[1024];

PROJECTION_Bin_Map Projection_binMap = GetBin_ProjectionMap();

for(auto BinMap :Projection_binMap ){
  sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_Bin_%i", BinMap.first);
  TH1D* Binning =  h_costheta_Pmu_UBTH2Poly->Projection(HistName, BinMap.second); 
  Binning->Write(HistName);
}

std::vector<int> StartingXbins{1,2,3,4,5,6,7,8,9};
std::vector<int> bin_numbers;

for(auto startingbin: StartingXbins){
//for(auto vec:BinMap.second ){std::cout<<vec<<",";}
  //std::cout<<" "<< std::endl;
  sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_Bin_%i", startingbin);
  int StartBin = startingbin;
  //std::cout<<"StartBin = "<< StartBin<< std::endl;
  TH1D* Bin =  h_costheta_Pmu_UBTH2Poly_inclusive->ProjectionY(HistName, StartBin, bin_numbers);
  Bin->Write(HistName);
}



outFile->Close();


saveUBTH2PolyToTextFile(*UBTH2Poly_costheta_Pmu, "Fake_Data_Rate_StevenBinning.txt");




}

int main(int argc, char* argv[]) {

 if (argc < 2) {
        // No command-line argument provided
        std::cout << "Usage: " << argv[0] << " <integer>" << std::endl;
        return 1; // Return an error code
    }

    // Convert the command-line argument to an integer
    int myInt = std::atoi(argv[1]);




   hist_maker_DATA(myInt);
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

