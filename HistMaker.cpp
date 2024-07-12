/*Author: Christian
 *
 * Usage: making hist to plot and look at for CC-0pi
// *
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
#include "includes/UBTH2Poly.h" 
#include "includes/HistUtils_cc0pi.hh"

#include "WeightHandler.hh"
#include "includes/Pmu_CorrectionFunction.hh"
#include "ConfigMakerUtils.hh"

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
/*NamedCategory<EventCategory>({kNuMuCC0pi0p},    "#nu_{#mu}-CC0#pi0P"),*/
//NamedCategory<EventCategory>({kNuMuBarCC},       "#bar{#nu_{#mu}}-CC"),
NamedCategory<EventCategory>({kNuMuCCOther},    "#nu_{#mu}-CCOther"),
NamedCategory<EventCategory>({kNuECC},          "#nu_{e}-CC"),
NamedCategory<EventCategory>({kNC},             "NC"),
NamedCategory<EventCategory>({kOOFV},           "Out FV"),
NamedCategory<EventCategory>({kOther},          "Other")
};

std::vector<EventCategory> SignalCategory_vector{kUnknown,kSignalCCQE, kSignalCCMEC, kSignalCCRES, kSignalOther, kNuMuCCNpi, kNuMuCCOther, kNuECC, kNC, kOOFV, kOther};


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
topology_categories = {
  NamedCategory<CCZeroPi_type>({kCC_0P_ZeroPi},          "CC0#pi"),
  NamedCategory<CCZeroPi_type>({kCC_1P_ZeroPi},          "CC1p0#pi"),
  NamedCategory<CCZeroPi_type>({kCC_2P_ZeroPi},          "CC2p0#pi"),
  NamedCategory<CCZeroPi_type>({kCC_3orGreaterP_ZeroPi}, "CC3>p0#pi"),
  NamedCategory<CCZeroPi_type>({kCC_ZeroPi_BG}, "BG")
};


//////Functions//////////////////////////////
void reorderVector(std::vector<int>& vec, int x);
void printText(const std::string& text, const std::string& fileName);
std::vector<int> MakereorderVector(std::vector<int> vec, int GotoFront);
int Num_MatchedPDGs(std::vector<int> inputPDGs , int MatchtoPDG);
 Double_t safe_weight_histMaker( Double_t w );
//void saveTH2ToTextFile(const UBTH2Poly& hist, char* fileName) ;



void hist_maker_MC(int input) {
   ///cc0pi_analyzer_1
 std::string base_dir = "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_4_17_2024_PmuCorrection/";
//cc0pi_analyzer_Run1
//cc0pi_analyzer_Run1/rustv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root
 ///uboone/data/users/cnguyen/RootFiles_tutorial/stv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.//1_v08_00_00_26_run1_reco2_reco2.root
 ///exp/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_4_17_2024_PmuCorrection/ru-stv_overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root  
 
 ///uboone/data/users/cnguyen/RootFiles_tutorial/stv-prodgenie_bnb_int
std::string Rootfilename;
if(input == 0){//rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
  Rootfilename = "ru-stv_overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root"; 
}
else if(input == 1){
  Rootfilename = "rustv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_NuWro_reco2_reco2.root";
}
else if(input == 3){
  Rootfilename = "rustv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";
}
 else if(input == 3){
  Rootfilename = "rustv-prodgenie_bnb_intrinsic_nue_overlay_run2_v08_00_00_35_run2a_reco2_reco2.root";
}
 
 
 //std::string Rootfile1name = "rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
 //std::string Rootfile2name = "stv-prodgenie_bnb_intrinsic_nue_overlay_run2_v08_00_00_35_run2a_reco2_reco2.root";
 //std::string Rootfile3name = "rustv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";
//rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
                             
 //std::string Rootfile4name = "stv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root"; // TRUTH
 // = "rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
 //std::string Rootfile1name = "rustv-run1_neutrinoselection_filt_numu_ALL.root";
// std::string FullRootFileName = "/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_2_15_24/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";// base_dir + Rootfile3name;
 
 std::string FullRootFileName = base_dir + Rootfilename;
 std::string OutPutRootName =  "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/MC/MC_Selected_newCC0pi_" + Rootfilename;
 
 ///uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_2_15_24/
 
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

  //std::vector< double > p_edges = {0.1,0.2,0.275,0.335,0.46,0.61,0.785, 0.985, 1.185, 1.385, 1.585, 2.0 }; //TODO: Rough measurements based on the "resolution" plot.
  //std::vector< double > cos_edges = {-1, -0.85, -0.775, -0.7, -0.625, -0.55, -0.475, -0.4, -0.325,
  //   -0.25, -0.175, -0.1, -0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5,
  //    0.575, 0.65, 0.725, 0.8, 0.85, 0.875, 0.9, 0.925, 0.950, 0.975, 1.00};


  std::vector <double> cos_edges = GetCCZeroPi_Binning(kBINNING_Costheta);
  std::vector <double> p_edges = GetCCZeroPi_Binning(kBINNING_Pmu);


//std::vector< double > cos_edges2D = {-1.0, 0.0, 0.27, 0.45, 0.62, 0.65, 0.76, 0.86, 0.94, 1.00};


std::vector< double > cos_edges2D = {-1, -.5, -0.85, -0.775, -0.7, -0.625, -0.55, -0.475, -0.4, -0.325,  -0.25, -0.175, -0.1, -0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85, 0.875, 0.9, 0.925, 0.950, 0.975, 1.00};
 //std::vector< double > p_edges_2D = {0.1, 0.2, 0.275, 0.335, 0.46, 0.61, 0.785, 1.185, 1.585 };
std::vector< double > p_edges_2D = {0.1,0.15, 0.2, 0.25, 0.3, .35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,  0.8, 0.85, 0.9, 0.95,  1.1 ,1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0};



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
  //std::vector< double > Probability_range= {0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.00 };
  
 std::vector< double > p_proton_edges = {0.25, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.2}; 
 std::vector< double > cos_edges_proton = {-1.0, -0.5, 0.0, 0.27, 0.45, 0.62, 0.76, 0.86,0.94, 1.0};
  std::vector< double > p_resolution_edges = {-1, -0.85, -0.8, -0.75, -0.7, -0.65, -0.60, -0.55, -0.50, -.45,-.4,-.35,-.3,-.25,-.2,-.15,-.1,-.05,0.,.05,.1,.15,.20,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1};
  
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
   
  std::vector<double> Cosmic_edges = generateBins(1, 40, 5);
  
  
  TGraph* graph_PonosFunction = new TGraph();


  TH1D* h_Proton_mult  = new TH1D("h_Proton_mult", "h_Proton_mult", Mult.size() - 1, Mult.data());
  TH1D* h_Proton_mult_Purity  = new TH1D("h_Proton_mult_Purity", "h_Proton_mult_Purity", Mult.size() - 1, Mult.data());
  h_Proton_mult_Purity->SetDirectory(0);h_Proton_mult_Purity->Sumw2();
   TH1D* h_CCPiCategory  = new TH1D("h_CCPiCategory", "h_CCPiCategory", Mult.size() - 1, Mult.data());
  
  TH1D* h_Proton_mult_TRUE  = new TH1D("h_Proton_mult_TRUE", "h_Proton_mult_TRUE", Mult.size() - 1, Mult.data());
  TH1D* h_Proton_mult_TRUE_RECO  = new TH1D("h_Proton_mult_TRUE_RECO", "h_Proton_mult_TRUE", Mult.size() - 1, Mult.data());
    h_Proton_mult_TRUE->SetDirectory(0);h_Proton_mult_TRUE->Sumw2();
    h_Proton_mult_TRUE_RECO->SetDirectory(0);h_Proton_mult_TRUE_RECO->Sumw2();
  
  
     
   TH2D* h_Proton_mult_Mig  = new TH2D("h_Proton_mult_Mig", "h_Proton_mult_Mig", Mult.size() - 1, Mult.data(),Mult.size() - 1,Mult.data());
    h_Proton_mult_Mig->SetDirectory(0); h_Proton_mult_Mig->Sumw2();
  
  
  
  
  HistFolio<TH1D, EventCategory> h_Proton_mult_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_Proton_mult_EventCategory", Mult, "h_Proton_mult_EventCategory; Events");

HistFolio<TH1D, EventCategory> h_Proton_mult_Tracks_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_Proton_mult_Tracks_EventCategory", Mult, "h_Proton_mult_Tracks_EventCategory; Events");


HistFolio<TH1D,  Particle_type> h_Proton_mult_Tracks_Particle_categories =
HistFolio<TH1D,Particle_type>(ParticleGroup_reduced_categories, "h_Proton_mult_Tracks_Particle_categories", Mult, "h_Proton_mult_Tracks_Particle_categories; Events");
  
  TH1D* h_Proton_mult_Tracks  = new TH1D("h_Proton_mult_Tracks", "h_Proton_mult_Tracks", Mult.size() - 1, Mult.data());
  TH1D* h_Proton_mult_Tracks_TRUE  = new TH1D("h_Proton_mult_Tracks_TRUE", "h_Proton_mult_Tracks_TRUE", Mult.size() - 1, Mult.data());
  TH1D* h_Proton_mult_Tracks_TRUE_RECO  = new TH1D("h_Proton_mult__TracksTRUE_RECO", "h_Proton_mult_Tracks_TRUE", Mult.size() - 1, Mult.data());
  TH1D* h_Proton_mult_Tracks_Purity  = new TH1D("h_Proton_mult_Tracks_Purity", "h_Proton_mult_Tracks_Purity", Mult.size() - 1, Mult.data());
  h_Proton_mult_Tracks->SetDirectory(0);h_Proton_mult_Tracks->Sumw2();
  h_Proton_mult_Tracks_TRUE->SetDirectory(0);h_Proton_mult_Tracks_TRUE->Sumw2();
  h_Proton_mult_Tracks_TRUE_RECO->SetDirectory(0);h_Proton_mult_Tracks_TRUE_RECO->Sumw2();
  h_Proton_mult_Tracks_Purity->SetDirectory(0);h_Proton_mult_Tracks_Purity->Sumw2();
  
  
  //TH1D* h_Proton_N_Tracks_TRUE  = new TH1D("h_Proton_N_Tracks_TRUE", "h_Proton_N_Tracks_TRUE", Mult.size() - 1, Mult.data());
  //TH1D* h_Proton_N_Tracks_TRUE_RECO  = new TH1D("h_Proton_N_Tracks_TRUE_RECO", "h_Proton_N_Tracks_TRUE_RECOE", Mult.size() - 1, Mult.data());
  
   // h_Proton_N_Tracks->SetDirectory(0);h_Proton_N_Tracks->Sumw2();
  //h_Proton_N_Tracks_TRUE->SetDirectory(0);h_Proton_N_Tracks_TRUE->Sumw2();
  
  
  
  
  TH1D* h_Pion_mult    = new TH1D("h_Pion_mult", "h_Pion_mult", Mult.size() - 1, Mult.data());
  TH1D* Failed_mult  = new TH1D("Failed_mult", "Failed_mult", Mult.size() - 1, Mult.data());  
  
  TH1D* h_NTracks  = new TH1D("h_NTracks", "h_NTracks", Mult.size() - 1, Mult.data());
  TH2D* h_NTracks_Mig  = new TH2D("h_NTracks_Mig", "", Mult.size() - 1, Mult.data(), Mult.size() - 1, Mult.data());

HistFolio<TH1D, EventCategory> h_NTracks_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_NTracks_EventCategory", Mult, "h_NTracks_EventCategory; Events");

HistFolio<TH1D, BDT_Category> h_NTracks_BTDGroup_categories =
HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_NTracks_BTDGroup_categories", Mult, "h_NTracks_EventCategory; Events");

HistFolio<TH1D, Particle_type> h_NTracks_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_NTracks_Particle_categories", Mult, "h_NTracks_Particle_categories; Events");


HistFolio<TH1D,  Particle_type> h_Else_mult_Tracks_Particle_categories =
HistFolio<TH1D,Particle_type>(ParticleGroup_reduced_categories, "h_Else_mult_Tracks_Particle_categories", Mult, "h_Else_mult_Tracks_Particle_categories; Events");

HistFolio<TH1D, EventCategory> h_Else_mult_Tracks_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_Else_mult_Tracks_EventCategory", Mult, "h_Else_mult_Tracks_EventCategory; Events");


TH1D* h_CCZeroP_Else_BTDPrediction  = new TH1D("h_CCZeroP_Else_BTDPrediction",
"h_CCZeroP_Else_BTDPrediction", Probability_range.size() - 1, Probability_range.data());

HistFolio<TH1D,  Particle_type> h_CCZeroP_Else_BTDPrediction_Particle_categories =
HistFolio<TH1D,Particle_type>(ParticleGroup_reduced_categories, "h_CCZeroP_Else_BTDPrediction_Particle_categories", Probability_range, "h_CCZeroP_Else_BTDPrediction_Particle_categories; Events");

HistFolio<TH1D, EventCategory> h_CCZeroP_Else_BTDPrediction_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_CCZeroP_Else_BTDPrediction_EventCategory", Probability_range, "hh_CCZeroP_Else_BTDPrediction_EventCategory; Events");

 TH2D* h_Else_Score_Pmu  = new TH2D("h_Else_Score_Pmu", "", Probability_range.size() - 1, Probability_range.data(), p_edges.size() - 1, p_edges.data());
 TH2D* h_Else_Score_CosTheta  = new TH2D("h_Else_Score_CosTheta", "", Probability_range.size() - 1, Probability_range.data(), cos_edges.size() - 1, cos_edges.data());



TH1D* h_CCZeroP_MuonCandidate_BTDPrediction  = new TH1D("h_CCZeroP_MuonCandidate_BTDPrediction", "h_CCZeroP_MuonCandidate_BTDPrediction", Probability_range.size() - 1, Probability_range.data());

TH1D* h_CCZeroP_MuonCandidate_BTDPrediction_purity  = new TH1D("h_CCZeroP_MuonCandidate_BTDPrediction_purity", "h_CCZeroP_MuonCandidate_BTDPrediction_purity", Probability_range.size() - 1, Probability_range.data());


TH1D* h_CCZeroP_MuonCandidate_BTDPrediction_TRUE  = new TH1D("h_CCZeroP_MuonCandidate_BTDPrediction_TRUE", "h_CCZeroP_MuonCandidate_BTDPrediction_TRUE", Probability_range.size() - 1, Probability_range.data());
h_CCZeroP_MuonCandidate_BTDPrediction_TRUE->SetDirectory(0);h_CCZeroP_MuonCandidate_BTDPrediction_TRUE->Sumw2();

TH1D* h_CCZeroP_MuonCandidate_BTDPrediction_TRUE_RECO  = new TH1D("h_CCZeroP_MuonCandidate_BTDPrediction_TRUE_RECO", "h_CCZeroP_MuonCandidate_BTDPrediction_TRUE_RECO", Probability_range.size() - 1, Probability_range.data());
h_CCZeroP_MuonCandidate_BTDPrediction_TRUE_RECO->SetDirectory(0);h_CCZeroP_MuonCandidate_BTDPrediction_TRUE_RECO->Sumw2();

HistFolio<TH1D, Particle_type> h_CCZeroP_MuonCandidate_BTDPrediction_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_CCZeroP_MuonCandidate_BTDPrediction_Particle_categories", Probability_range, "h_CCZeroP_MuonCandidate_BTDPrediction_Particle_categories; Events");


HistFolio<TH1D, EventCategory> h_CCZeroP_MuonCandidate_BTDPrediction_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_CCZeroP_MuonCandidate_BTDPrediction_EventCategory", Probability_range, "h_CCZeroP_MuonCandidate_BTDPrediction_EventCategory; Events");

HistFolio<TH1D, CCZeroPi_type> h_CCZeroP_MuonCandidate_BTDPrediction_topology_categories =
 HistFolio<TH1D, CCZeroPi_type>(topology_categories, "h_CCZeroP_MuonCandidate_BTDPrediction_topology_categories", Probability_range, "h_CCZeroP_MuonCandidate_BTDPrediction_topology_categories");
 
 HistFolio<TH1D, FidVol_Category> h_CCZeroP_MuonCandidate_BTDPrediction_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_CCZeroP_MuonCandidate_BTDPrediction_containment_categories", Probability_range, "h_CCZeroP_MuonCandidate_BTDPrediction_containment_categories");


TH1D* h_topological_score  = new TH1D("h_topological_score", "h_topological_score", Probability_range.size() - 1, Probability_range.data());

TH1D* h_topological_score_TRUE  = new TH1D("h_topological_score_TRUE", "h_topological_score_TRUE", Probability_range.size() - 1, Probability_range.data());
h_topological_score_TRUE->SetDirectory(0);h_topological_score_TRUE->Sumw2();
TH1D* h_topological_score_TRUE_RECO  = new TH1D("h_topological_score_TRUE_RECO", "h_topological_score_TRUE_RECO", Probability_range.size() - 1, Probability_range.data());
h_topological_score_TRUE_RECO->SetDirectory(0);h_topological_score_TRUE_RECO->Sumw2();


h_topological_score->SetDirectory(0);h_topological_score->Sumw2();
 
HistFolio<TH1D, Particle_type> h_topological_score_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_topological_score_Particle_categories", Probability_range, "h_topological_score_Particle_categories; Events");

HistFolio<TH1D, EventCategory> h_topological_score_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_topological_score_EventCategory", Probability_range, "h_topological_score_EventCategory; Events");



TH1D* h_trk_score_v  = new TH1D("h_trk_score_v", "h_trk_score_v",
Probability_range.size() - 1, Probability_range.data());
h_trk_score_v->SetDirectory(0);h_trk_score_v->Sumw2();
 
 
 TH1D* h_trk_score_v_TRUE  = new TH1D("h_trk_score_v_TRUE", 
 "h_trk_score_v_TRUE", Probability_range.size() - 1, Probability_range.data());
h_trk_score_v_TRUE->SetDirectory(0);h_trk_score_v_TRUE->Sumw2();
 
 
  
 TH1D* h_trk_score_v_TRUE_RECO  = new TH1D("h_trk_score_v_TRUE_RECO", 
 "h_trk_score_v_TRUE_RECO", Probability_range.size() - 1, Probability_range.data());
h_trk_score_v_TRUE_RECO->SetDirectory(0);h_trk_score_v_TRUE_RECO->Sumw2();
 
 
HistFolio<TH1D, Particle_type> h_trk_score_v_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_trk_score_v_Particle_categories", Probability_range, "trk_score_v_Particle_categories; Events");

HistFolio<TH1D, EventCategory> h_trk_score_v_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_trk_score_v_EventCategory", Probability_range, "trk_score_v_EventCategory; Events");



TH1D* h_trk_len_v  = new TH1D("h_trk_len_v", "h_trk_len_v", trk_len_v_edges.size() - 1, trk_len_v_edges.data());
h_trk_len_v->SetDirectory(0);h_trk_len_v->Sumw2();
 
 
TH1D* h_trk_len_v_Purity  = new TH1D("h_trk_len_v_Purity", "h_trk_len_v_Purity", trk_len_v_edges.size() - 1, trk_len_v_edges.data());
h_trk_len_v_Purity->SetDirectory(0);h_trk_len_v_Purity->Sumw2();
 
 
TH1D* h_trk_len_v_TRUE  = new TH1D("h_trk_len_v_TRUE", "h_trk_len_v_TRUE", trk_len_v_edges.size() - 1, trk_len_v_edges.data());
h_trk_len_v_TRUE->SetDirectory(0);h_trk_len_v_TRUE->Sumw2(); 
 
TH1D* h_trk_len_v_TRUE_RECO  = new TH1D("h_trk_len_v_TRUE_RECO", "h_trk_len_v_TRUE", trk_len_v_edges.size() - 1, trk_len_v_edges.data());
h_trk_len_v_TRUE_RECO->SetDirectory(0);h_trk_len_v_TRUE_RECO->Sumw2();  
 
HistFolio<TH1D, Particle_type> h_trk_len_v_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_trk_len_v_Particle_categories", trk_len_v_edges, "trk_score_v_Particle_categories; Events");


HistFolio<TH1D, EventCategory> h_trk_len_v_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_trk_len_v_EventCategory", trk_len_v_edges, "trk_score_v_EventCategory; Events");


TH1D* h_trk_distance_v  = new TH1D("h_trk_distance_v", "h_trk_distance_v", trk_distance_v_edges.size() - 1, trk_distance_v_edges.data());
h_trk_distance_v->SetDirectory(0);h_trk_distance_v->Sumw2();
 
 
 TH1D* h_trk_distance_v_TRUE  = new TH1D("h_trk_distance_v_TRUE", "h_trk_distance_v_TRUE", trk_distance_v_edges.size() - 1, trk_distance_v_edges.data());
  h_trk_distance_v_TRUE->SetDirectory(0);h_trk_distance_v_TRUE->Sumw2();
 
  TH1D* h_trk_distance_v_TRUE_RECO  = new TH1D("h_trk_distance_v_TRUE_RECO", "h_trk_distance_v_TRUE_RECO", trk_distance_v_edges.size() - 1, trk_distance_v_edges.data());
h_trk_distance_v_TRUE_RECO->SetDirectory(0);h_trk_distance_v_TRUE_RECO->Sumw2();
 
 
HistFolio<TH1D, Particle_type> h_trk_distance_v_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_trk_distance_v_Particle_categories", trk_distance_v_edges, "h_trk_distance_v_Particle_categories; Events");

HistFolio<TH1D, EventCategory> h_trk_distance_v_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_trk_distance_v_EventCategory", trk_distance_v_edges, "h_trk_distance_v_EventCategory; Events");



TH1D* h_CCZeroP_ProtonCandidate_BTDPrediction  = new TH1D("h_CCZeroP_ProtonCandidate_BTDPrediction", "h_CCZeroP_ProtonCandidate_BTDPrediction", Probability_range.size() - 1, Probability_range.data());

HistFolio<TH1D, Particle_type> h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories", Probability_range, "h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories; Events");

HistFolio<TH1D, EventCategory> h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory =
HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory", Probability_range, "h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory; Events");



HistFolio<TH1D, CCZeroPi_type> h_CCZeroP_ProtonCandidate_BTDPrediction_topology_categories =
 HistFolio<TH1D, CCZeroPi_type>(topology_categories, "h_CCZeroP_ProtonCandidate_BTDPrediction_topology_categories", Probability_range, "h_CCZeroP_ProtonCandidate_BTDPrediction_topology_categories");
 
 HistFolio<TH1D, FidVol_Category> h_CCZeroP_ProtonCandidate_BTDPrediction_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_CCZeroP_ProtonCandidate_BTDPrediction_containment_categories", Probability_range, "h_CCZeroP_ProtonCandidate_BTDPrediction_containment_categories");





TH1D* h_CCZeroP_ALLProtonCandidate_BTDPrediction  = new TH1D("h_CCZeroP_ALLProtonCandidate_BTDPrediction", "h_CCZeroP_ALLProtonCandidate_BTDPrediction", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE  = new TH1D("h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE", "h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE_RECO  = new TH1D("h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE_RECO", "h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE_RECO", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_ALLProtonCandidate_BTDPrediction_Purity  = new TH1D("h_CCZeroP_ALLProtonCandidate_BTDPrediction_Purity", "h_CCZeroP_ALLProtonCandidate_BTDPrediction_Purity", Probability_range.size() - 1, Probability_range.data());


TH1D* h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1  = new TH1D("h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1", "h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_BTDGroup_Probability_tobeProton_trk1  = new TH1D("h_CCZeroP_BTDGroup_Probability_tobeProton_trk1", "h_CCZeroP_BTDGroup_Probability_tobeProton_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_BTDGroup_Probability_tobePion_trk1  = new TH1D("h_CCZeroP_BTDGroup_Probability_tobePion_trk1", "h_CCZeroP_BTDGroup_Probability_tobePion_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CCZeroP_BTDGroup_Probability_tobeElse_trk1  = new TH1D("h_CCZeroP_BTDGroup_Probability_tobeElse_trk1", "h_CCZeroP_BTDGroup_Probability_tobeElse_trk1", Probability_range.size() - 1, Probability_range.data());


TH1D* h_CC1P_BTDGroup_Probability_tobeMuon_tk1  = new TH1D("h_CC1P_BTDGroup_Probability_tobeMuon_tk1", "h_CC1P_BTDGroup_Probability_tobeMuon_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobeProton_tk1  = new TH1D("h_CC1P_BTDGroup_Probability_tobeProton_tk1", "h_CC1P_BTDGroup_Probability_tobeProton_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobePion_tk1  = new TH1D("h_CC1P_BTDGroup_Probability_tobePion_tk1", "h_CC1P_BTDGroup_Probability_tobePion_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobeElse_tk1  = new TH1D("h_CC1P_BTDGroup_Probability_tobeElse_tk1", "h_CC1P_BTDGroup_Probability_tobeElse_trk1", Probability_range.size() - 1, Probability_range.data());

TH1D* h_CC1P_BTDGroup_Probability_tobeMuon_tk2  = new TH1D("h_CC1P_BTDGroup_Probability_tobeMuon_tk2", "h_CC1P_BTDGroup_Probability_tobeMuon_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobeProton_tk2  = new TH1D("h_CC1P_BTDGroup_Probability_tobeProton_tk2", "h_CC1P_BTDGroup_Probability_tobeProton_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobePion_tk2  = new TH1D("h_CC1P_BTDGroup_Probability_tobePion_tk2", "h_CC1P_BTDGroup_Probability_tobePion_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC1P_BTDGroup_Probability_tobeElse_tk2  = new TH1D("h_CC1P_BTDGroup_Probability_tobeElse_tk2", "h_CC1P_BTDGroup_Probability_tobeElse_trk2", Probability_range.size() - 1, Probability_range.data());
//////////////////////////////
//cc2p
/////////////////////////////
TH1D* h_CC2P_BTDGroup_Probability_tobeMuon_tk1  = new TH1D("h_CC2P_BTDGroup_Probability_tobeMuon_tk1", "h_CC2P_BTDGroup_Probability_tobeMuon_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeProton_tk1  = new TH1D("h_CC2P_BTDGroup_Probability_tobeProton_tk1", "h_CC2P_BTDGroup_Probability_tobeProton_trk1", Probability_range.size() - 1, Probability_range.data());
/////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* h_CC2P_BTDGroup_Probability_tobePion_tk1  = new TH1D("h_CC2P_BTDGroup_Probability_tobePion_tk1", "h_CC2P_BTDGroup_Probability_tobePion_trk1", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeElse_tk1  = new TH1D("h_CC2P_BTDGroup_Probability_tobeElse_tk1", "h_CC2P_BTDGroup_Probability_tobeElse_trk1", Probability_range.size() - 1, Probability_range.data());

TH1D* h_CC2P_BTDGroup_Probability_tobeMuon_tk2  = new TH1D("h_CC2P_BTDGroup_Probability_tobeMuon_tk2", "h_CC2P_BTDGroup_Probability_tobeMuon_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeProton_tk2  = new TH1D("h_CC2P_BTDGroup_Probability_tobeProton_tk2", "h_CC2P_BTDGroup_Probability_tobeProton_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobePion_tk2  = new TH1D("h_CC2P_BTDGroup_Probability_tobePion_tk2", "h_CC2P_BTDGroup_Probability_tobePion_trk2", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeElse_tk2  = new TH1D("h_CC2P_BTDGroup_Probability_tobeElse_tk2", "h_CC2P_BTDGroup_Probability_tobeElse_trk2", Probability_range.size() - 1, Probability_range.data());


TH1D* h_CC2P_BTDGroup_Probability_tobeMuon_tk3  = new TH1D("h_CC2P_BTDGroup_Probability_tobeMuon_tk3", "h_CC2P_BTDGroup_Probability_tobeMuon_trk3", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeProton_tk3  = new TH1D("h_CC2P_BTDGroup_Probability_tobeProton_tk3", "h_CC2P_BTDGroup_Probability_tobeProton_trk3", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobePion_tk3  = new TH1D("h_CC2P_BTDGroup_Probability_tobePion_tk3", "h_CC2P_BTDGroup_Probability_tobePion_trk3", Probability_range.size() - 1, Probability_range.data());
TH1D* h_CC2P_BTDGroup_Probability_tobeElse_tk3  = new TH1D("h_CC2P_BTDGroup_Probability_tobeElse_tk3", "h_CC2P_BTDGroup_Probability_tobeElse_trk3", Probability_range.size() - 1, Probability_range.data());



HistFolio<TH1D, FidVol_Category> h_NTracks_containment_categories =
HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_NTracks_containment_categories", Mult, "h_NTracks_containment_categories");


HistFolio<TH1D, Particle_type> h_ProtonMupl_Particle_categories =
HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_ProtonMupl_Particle_categories", Mult, "h_ProtonMupl_Particle_categories; Events");


HistFolio<TH1D, BDT_Category> h_ProtonMupl_BTDGroup_categories =
HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_ProtonMupl_BTDGroup_categories", Mult, "h_ProtonMupl_BTDGroup_categories; Events");

  /////////////////////////////////ProtonMupl/////////////////////////////////////////////////
  TH1D* h_p             = new TH1D("h_p", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  
  TH1D* h_p_Correction             = new TH1D("h_p_Correction", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  
  TH1D* h_p_true        = new TH1D("h_p_true", ";p_{#mu}^{true}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_Resolution  = new TH1D("h_p_Resolution", ";p_{#mu}^{Resolution}", p_resolution_edges.size() - 1, p_resolution_edges.data());
  TH1D* h_p_MCS_Resolution  = new TH1D("h_p_MCS_Resolution", ";p_{#mu}^{Resolution}", p_resolution_edges.size() - 1, p_resolution_edges.data());
  TH2D* h_p_Mig  = new TH2D("h_p_Mig", "", p_edges.size() - 1, p_edges.data(), p_edges.size() - 1, p_edges.data());
  TH2D* h_p_MCS_Mig  = new TH2D("h_p_MCS_Mig", "", p_edges.size() - 1, p_edges.data(), p_edges.size() - 1, p_edges.data());
  TH2D* h_pionE_Mig  = new TH2D("h_pionE_Mig", "", p_edges.size() - 1, p_edges.data(), p_edges.size() - 1, p_edges.data());
  TH1D* h_p_purity             = new TH1D("h_p_purity", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());


  TH1D* h_p_TRUE             = new TH1D("h_p_TRUE", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_TRUE_RECO             = new TH1D("h_p_TRUE_RECO", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());


 HistFolio<TH1D, EventCategory> h_p_EventCategory =
 HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_p_EventCategory", p_edges, "h_p_EventCategory; Events");

 HistFolio<TH1D, EventCategory> h_p_Correction_EventCategory =
 HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_p_Correction_EventCategory", p_edges, "h_p_Correction_EventCategory; Events");


 HistFolio<TH1D, BDT_Category> h_p_BTDGroup_categories =
 HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_p_BTDGroup_categories", p_edges, "h_p_BTDGroup_categories; Events");

 HistFolio<TH1D, Particle_type> h_p_Particle_categories =
 HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_p_Particle_categories", p_edges, "h_p_Particle_categories; Events");
 
 HistFolio<TH1D, FidVol_Category> h_p_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_p_containment_categories", p_edges, "h_p_containment_categories");

 HistFolio<TH1D, CCZeroPi_type> h_p_topology_categories =
 HistFolio<TH1D, CCZeroPi_type>(topology_categories, "h_p_topology_categories", p_edges, "h_p_topology_categories");
 
 HistFolio<TH1D, FidVol_Category> h_p_Correction_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_p_Correction_containment_categories", p_edges, "h_p_Correction_containment_categories");

 HistFolio<TH1D, CCZeroPi_type> h_p_Correction_topology_categories =
 HistFolio<TH1D, CCZeroPi_type>(topology_categories, "h_p_Correction_topology_categories", p_edges, "h_p_Correction_topology_categories");


  h_p->SetDirectory(0); h_p->Sumw2();
  h_p_true->SetDirectory(0); h_p_true->Sumw2();
  h_p_Resolution->SetDirectory(0); h_p_Resolution->Sumw2();
  h_p_MCS_Resolution->SetDirectory(0);h_p_MCS_Resolution->Sumw2();
  h_p_TRUE->SetDirectory(0); h_p_TRUE->Sumw2();
  h_p_TRUE_RECO->SetDirectory(0); h_p_TRUE_RECO->Sumw2();
  h_p_purity->SetDirectory(0); h_p_purity->Sumw2();
  
  h_p_Mig->SetDirectory(0); h_p_Mig->Sumw2();
  h_p_MCS_Mig->SetDirectory(0); h_p_MCS_Mig->Sumw2();
  h_pionE_Mig->SetDirectory(0); h_pionE_Mig->Sumw2();
  //////////////////////////////////////////////////////////////////////////////////////////
  //CC0Pi_only
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_p_CC0P0Pi_Only             = new TH1D("h_p_CC0P0Pi_Only", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_CC0P0Pi_Only_true        = new TH1D("h_p_CC0P0Pi_Only_true", ";p_{#mu}^{true}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_CC0P0Pi_Only_Resolution  = new TH1D("h_p_CC0P0Pi_Only_Resolution", ";p_{#mu}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_p_CC0P0Pi_Only_Mig  = new TH2D("h_p_CC0P0Pi_Only_Mig", "", p_edges.size() - 1, p_edges.data(), p_edges.size() - 1, p_edges.data());

 HistFolio<TH1D, EventCategory> h_p_CC0P0Pi_Only_EventCategory =
 HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_p_CC0P0Pi_Only_EventCategory", p_edges, "h_p_CC0P0Pi_Only_EventCategory; Events");

 HistFolio<TH1D, BDT_Category> h_p_CC0P0Pi_Only_BTDGroup_categories =
 HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_p_CC0P0Pi_Only_BTDGroup_categories", p_edges, "h_p_CC0P0Pi_Only_BTDGroup_categories; Events");

 HistFolio<TH1D, Particle_type> h_p_CC0P0Pi_Only_Particle_categories =
 HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_p_CC0P0Pi_Only_Particle_categories", p_edges, "h_p_CC0P0Pi_Only_Particle_categories; Events");
 
 HistFolio<TH1D, FidVol_Category> h_p_CC0P0Pi_Only_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_p_CC0P0Pi_Only_containment_categories", p_edges, "h_p_CC0P0Pi_Only_containment_categories");


  h_p_CC0P0Pi_Only->SetDirectory(0); h_p_CC0P0Pi_Only->Sumw2();
  h_p_CC0P0Pi_Only_true->SetDirectory(0); h_p_CC0P0Pi_Only_true->Sumw2();
  h_p_CC0P0Pi_Only_Resolution->SetDirectory(0); h_p_CC0P0Pi_Only_Resolution->Sumw2();
  h_p_CC0P0Pi_Only_Mig->SetDirectory(0); h_p_CC0P0Pi_Only_Mig->Sumw2(); 
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////////////////////////////////////
  //CCNP0Pi_only
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_p_CCNP0Pi_Only             = new TH1D("h_p_CCNP0Pi_Only", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_CCNP0Pi_Only_true        = new TH1D("h_p_CCNP0Pi_Only_true", ";p_{#mu}^{true}", p_edges.size() - 1, p_edges.data());
  TH1D* h_p_CCNP0Pi_Only_Resolution  = new TH1D("h_p_CCNP0Pi_Only_Resolution", ";p_{#mu}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_p_CCNP0Pi_Only_Mig  = new TH2D("h_p_CCNP0Pi_Only_Mig", "", p_edges.size() - 1, p_edges.data(), p_edges.size() - 1, p_edges.data());

 HistFolio<TH1D, EventCategory> h_p_CCNP0Pi_Only_EventCategory =
 HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_p_CCNP0Pi_Only_EventCategory", p_edges, "h_p_CCNP0Pi_Only_EventCategory; Events");

 HistFolio<TH1D, BDT_Category> h_p_CCNP0Pi_Only_BTDGroup_categories =
 HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_p_CCNP0Pi_Only_BTDGroup_categories", p_edges, "h_p_CCNP0Pi_Only_BTDGroup_categories; Events");

 HistFolio<TH1D, Particle_type> h_p_CCNP0Pi_Only_Particle_categories =
 HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_p_CCNP0Pi_Only_Particle_categories", p_edges, "h_p_CCNP0Pi_Only_Particle_categories; Events");
 
 HistFolio<TH1D, FidVol_Category> h_p_CCNP0Pi_Only_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_p_CCNP0Pi_Only_containment_categories", p_edges, "h_p_CCNP0Pi_Only_containment_categories");


  h_p_CCNP0Pi_Only->SetDirectory(0); h_p_CCNP0Pi_Only->Sumw2();
  h_p_CCNP0Pi_Only_true->SetDirectory(0); h_p_CCNP0Pi_Only_true->Sumw2();
  h_p_CCNP0Pi_Only_Resolution->SetDirectory(0); h_p_CCNP0Pi_Only_Resolution->Sumw2();
  h_p_CCNP0Pi_Only_Mig->SetDirectory(0); h_p_CCNP0Pi_Only_Mig->Sumw2(); 
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  //////////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_costheta_p             = new TH1D("h_costheta_p", ";p_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data());
  TH1D* h_costheta_p_true        = new TH1D("h_costheta_p_true", ";p_{#mu}^{true}", cos_edges.size() - 1, cos_edges.data());
  TH1D* h_costheta_p_Resolution  = new TH1D("h_costheta_p_Resolution", ";p_{#mu}^{Resolution}", p_resolution_edges.size() - 1, p_resolution_edges.data());
  
  TH1D* h_costheta_p_MCS_Resolution  = new TH1D("h_costheta_p_MCS_Resolution", ";p_{#mu}^{Resolution}", p_resolution_edges.size() - 1, p_resolution_edges.data());
  TH2D* h_costheta_p_Mig  = new TH2D("h_costheta_p_Mig", "", cos_edges.size() - 1, cos_edges.data(), cos_edges.size() - 1, cos_edges.data());
  TH2D* h_costheta_p_MCS_Mig  = new TH2D("h_costheta_p_MCS_Mig", "", cos_edges.size() - 1, cos_edges.data(), cos_edges.size() - 1, cos_edges.data());

  h_costheta_p->SetDirectory(0); h_costheta_p->Sumw2();
  h_costheta_p_true->SetDirectory(0); h_costheta_p_true->Sumw2();
  h_costheta_p_Resolution->SetDirectory(0); h_costheta_p_Resolution->Sumw2();
  h_costheta_p_MCS_Resolution->SetDirectory(0); h_costheta_p_MCS_Resolution->Sumw2();
  h_costheta_p_Mig->SetDirectory(0); h_costheta_p_Mig->Sumw2();
  h_costheta_p_MCS_Mig->SetDirectory(0); h_costheta_p_MCS_Mig->Sumw2();
  
  
    TH1D* h_costheta_p_TRUE             = new TH1D("h_costheta_p_TRUE", ";Cos(#theta_{#mu}", cos_edges.size() - 1, cos_edges.data());
    TH1D* h_costheta_p_TRUE_RECO             = new TH1D("h_costheta_p_TRUE_RECO", ";Cos(#theta_{#mu}", cos_edges.size() - 1, cos_edges.data());
    TH1D* h_costheta_p_purity             = new TH1D("h_costheta_p_purity", ";Cos(#theta_{#mu}", cos_edges.size() - 1, cos_edges.data());
  
  
    h_costheta_p_TRUE->SetDirectory(0); h_costheta_p_TRUE->Sumw2();
    h_costheta_p_TRUE_RECO->SetDirectory(0); h_costheta_p_TRUE_RECO->Sumw2();
    h_costheta_p_purity->SetDirectory(0); h_costheta_p_purity->Sumw2();
  
  
  HistFolio<TH1D, EventCategory> h_costheta_p_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_costheta_p_EventCategory", cos_edges, "h_costheta_p_EventCategory; Events");

  HistFolio<TH1D, BDT_Category> h_costheta_p_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_costheta_p_BTDGroup_categories", cos_edges, "h_costheta_p_EventCategory; Events");
  
  HistFolio<TH1D, Particle_type> h_costheta_p_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_costheta_p_Particle_categories", cos_edges, "h_costheta_p_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_costheta_p_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_costheta_p_containment_categories", cos_edges, "h_costheta_p_containment_categories");
  
  
  HistFolio<TH1D, CCZeroPi_type> h_costheta_p_topology_categories =
  HistFolio<TH1D, CCZeroPi_type>(topology_categories, "h_costheta_p_topology_categories", cos_edges, "h_costheta_p_topology_categories");
/////////////////////////////////////////////////////////////////////////  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///CC0Pi0P
  //////////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_costheta_p_CC0P0Pi_Only             = new TH1D("h_costheta_p_CC0P0Pi_Only", ";p_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data());
  TH1D* h_costheta_p_CC0P0Pi_Only_true        = new TH1D("h_costheta_p_CC0P0Pi_Only_true", ";p_{#mu}^{true}", cos_edges.size() - 1, cos_edges.data());
  TH1D* h_costheta_p_CC0P0Pi_Only_Resolution  = new TH1D("h_costheta_p_CC0P0Pi_Only_Resolution", ";p_{#mu}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_costheta_p_CC0P0Pi_Only_Mig  = new TH2D("h_costheta_p_CC0P0Pi_Only_Mig", "", cos_edges.size() - 1, cos_edges.data(), cos_edges.size() - 1, cos_edges.data());
  h_costheta_p_CC0P0Pi_Only->SetDirectory(0); h_costheta_p_CC0P0Pi_Only->Sumw2();
  h_costheta_p_CC0P0Pi_Only_true->SetDirectory(0); h_costheta_p_CC0P0Pi_Only_true->Sumw2();
  h_costheta_p_CC0P0Pi_Only_Resolution->SetDirectory(0); h_costheta_p_CC0P0Pi_Only_Resolution->Sumw2();
  h_costheta_p_CC0P0Pi_Only_Mig->SetDirectory(0); h_costheta_p_CC0P0Pi_Only_Mig->Sumw2();
  
  
  HistFolio<TH1D, EventCategory> h_costheta_p_CC0P0Pi_Only_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_costheta_p_CC0P0Pi_Only_EventCategory", cos_edges, "h_costheta_p_CC0P0Pi_Only_EventCategory; Events");

  HistFolio<TH1D, BDT_Category> h_costheta_p_CC0P0Pi_Only_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_costheta_p_CC0P0Pi_Only_BTDGroup_categories", cos_edges, "h_costheta_p_CC0P0Pi_Only_EventCategory; Events");
  
  HistFolio<TH1D, Particle_type> h_costheta_p_CC0P0Pi_Only_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_costheta_p_CC0P0Pi_Only_Particle_categories", cos_edges, "h_costheta_p_CC0P0Pi_Only_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_costheta_p_CC0P0Pi_Only_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_costheta_p_CC0P0Pi_Only_containment_categories", cos_edges, "h_costheta_p_CC0P0Pi_Only_containment_categories");
  
  
/////////////////////////////////////////////////////////////////////////  


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///CC0PiNP
//////////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_costheta_p_CCNP0Pi_Only             = new TH1D("h_costheta_p_CCNP0Pi_Only", ";p_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data());
  TH1D* h_costheta_p_CCNP0Pi_Only_true        = new TH1D("h_costheta_p_CCNP0Pi_Only_true", ";p_{#mu}^{true}", cos_edges.size() - 1, cos_edges.data());
  TH1D* h_costheta_p_CCNP0Pi_Only_Resolution  = new TH1D("h_costheta_p_CCNP0Pi_Only_Resolution", ";p_{#mu}^{Resolution}", cos_edges.size() - 1, cos_edges.data());
  TH2D* h_costheta_p_CCNP0Pi_Only_Mig  = new TH2D("h_costheta_p_CCNP0Pi_Only_Mig", "", cos_edges.size() - 1, cos_edges.data(), cos_edges.size() - 1, cos_edges.data());
  
  h_costheta_p_CCNP0Pi_Only->SetDirectory(0); h_costheta_p_CCNP0Pi_Only->Sumw2();
  h_costheta_p_CCNP0Pi_Only_true->SetDirectory(0); h_costheta_p_CCNP0Pi_Only_true->Sumw2();
  h_costheta_p_CCNP0Pi_Only_Resolution->SetDirectory(0); h_costheta_p_CCNP0Pi_Only_Resolution->Sumw2();
  h_costheta_p_CCNP0Pi_Only_Mig->SetDirectory(0); h_costheta_p_CCNP0Pi_Only_Mig->Sumw2();
  
  
  HistFolio<TH1D, EventCategory> h_costheta_p_CCNP0Pi_Only_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_costheta_p_CCNP0Pi_Only_EventCategory", cos_edges, "h_costheta_p_CCNP0Pi_Only_EventCategory; Events");

  HistFolio<TH1D, BDT_Category> h_costheta_p_CCNP0Pi_Only_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_costheta_p_CCNP0Pi_Only_BTDGroup_categories", cos_edges, "h_costheta_p_CCNP0Pi_Only_EventCategory; Events");
  
  HistFolio<TH1D, Particle_type> h_costheta_p_CCNP0Pi_Only_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_costheta_p_CCNP0Pi_Only_Particle_categories", cos_edges, "h_costheta_p_CCNP0Pi_Only_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_costheta_p_CCNP0Pi_Only_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_costheta_p_CCNP0Pi_Only_containment_categories", cos_edges, "h_costheta_p_CCNP0Pi_Only_containment_categories");
  
  
/////////////////////////////////////////////////////////////////////////  
//////  
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

  HistFolio<TH1D, BDT_Category> h_pn_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_pn_BTDGroup_categories", pn_edges, "h_pn_EventCategory; Events");
  
  HistFolio<TH1D, Particle_type> h_pn_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_pn_Particle_categories", pn_edges, "h_pn_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_pn_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_pn_containment_categories", pn_edges, "h_pn_containment_categories");
    ////////////////////
    TH1D* h_leadingProton_p             = new TH1D("h_leadingProton_p", ";p_{p}^{reco}", p_proton_edges.size() - 1, p_proton_edges.data());
     h_leadingProton_p->SetDirectory(0); h_leadingProton_p->Sumw2();
  
  
  HistFolio<TH1D, EventCategory> h_leadingProton_p_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_leadingProton_p_EventCategory", p_proton_edges, "h_leadingProton_p_EventCategory; Events");

  
  
  TH1D* h_costheta_proton             = new TH1D("h_costheta_proton", ";p_{#mu}^{reco}", cos_edges_proton.size() - 1, cos_edges_proton.data());
   h_costheta_proton->SetDirectory(0); h_costheta_proton->Sumw2();
  
    HistFolio<TH1D, EventCategory> h_costheta_proton_EventCategory =
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_costheta_proton_EventCategory", cos_edges_proton, "h_costheta_proton_EventCategory; Events");
  
  
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

  HistFolio<TH1D, BDT_Category> h_delta_alphaT_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_delta_alphaT_BTDGroup_categories", delta_alphaT_edges, "h_delta_alphaT_EventCategory; Events");
  
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

  HistFolio<TH1D, BDT_Category> h_delta_pTx_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_delta_pTx_BTDGroup_categories", delta_pTx_edges, "h_delta_pTx_EventCategory; Events");
  
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

  HistFolio<TH1D, BDT_Category> h_delta_pTy_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_delta_pTy_BTDGroup_categories", delta_pTy_edges, "h_delta_pTy_EventCategory; Events");
  
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
  
  HistFolio<TH1D, BDT_Category> h_delta_phiT_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_delta_phiT_BTDGroup_categories", delta_phiT_edges, "h_delta_phiT_EventCategory; Events");
  
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
  
  HistFolio<TH1D, BDT_Category> h_openingAngle_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_openingAngle_BTDGroup_categories", openingAngle_edges, "h_openingAngle_EventCategory; Events");
  
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
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_openingAngle_deg_EventCategory", openingAngleDeg_edges, "h_openingAngle_deg_EventCategory; Events");

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
  
  HistFolio<TH1D, BDT_Category> h_VertexX_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_VertexX_BTDGroup_categories", VertexX_edges, "h_VertexX_EventCategory; Events");
  
  HistFolio<TH1D, Particle_type> h_VertexX_leadmuon_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_VertexX_leadmuon_Particle_categories", VertexX_edges, "h_VertexX_leadmuon_Particle_categories; Events");
  
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
  
  HistFolio<TH1D, BDT_Category> h_VertexY_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_VertexY_BTDGroup_categories", VertexY_edges, "h_VertexY_EventCategory; Events");
  
  HistFolio<TH1D, Particle_type> h_VertexY_leadmuon_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_VertexY_leadmuon_Particle_categories", VertexY_edges, "h_VertexY_leadmuon_Particle_categories; Events");
  
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
  HistFolio<TH1D, EventCategory>(EventSelectionGroup_categories, "h_VertexZ_EventCategory", VertexZ_edges, "h_VertexZ_EventCategory; Events");
  
  HistFolio<TH1D, BDT_Category> h_VertexZ_BTDGroup_categories =
  HistFolio<TH1D, BDT_Category>(BTDGroup_categories, "h_VertexZ_BTDGroup_categories", VertexZ_edges, "h_VertexZ_EventCategory; Events");
  
  HistFolio<TH1D, Particle_type> h_VertexZ_leadmuon_Particle_categories =
  HistFolio<TH1D, Particle_type>(ParticleGroup_reduced_categories, "h_VertexZ_leadmuon_Particle_categories", VertexZ_edges, "h_VertexZ_leadmuon_Particle_categories; Events");
  
  HistFolio<TH1D, FidVol_Category> h_VertexZ_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_VertexZ_containment_categories", VertexZ_edges, "h_VertexZ_containment_categories");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  //  2D
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH2D* h_costheta_Pmu  = new TH2D("h_costheta_Pmu", "h_costheta_Pmu", cos_edges2D.size() - 1, cos_edges2D.data(), p_edges_2D.size() - 1, p_edges_2D.data());


  
  
  
  HistFolio<TH2D, EventCategory> h_costheta_Pmu_EventCategory =
  HistFolio<TH2D, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_EventCategory", cos_edges2D, p_edges_2D,"h_costheta_Pmu_EventCategory; Events");


   TH2D* h_costheta_Pmu_TRUE  = new TH2D("h_costheta_Pmu_TRUE", "h_costheta_Pmu_TRUE", cos_edges2D.size() - 1, cos_edges2D.data(), p_edges_2D.size() - 1, p_edges_2D.data());
   h_costheta_Pmu_TRUE->SetDirectory(0); h_costheta_Pmu_TRUE->Sumw2();
   TH2D* h_costheta_Pmu_TRUE_RECO  = new TH2D("h_costheta_Pmu_TRUE_RECO", "h_costheta_Pmu_TRUE_RECO", cos_edges2D.size() - 1, cos_edges2D.data(), p_edges_2D.size() - 1, p_edges_2D.data());
    h_costheta_Pmu_TRUE_RECO->SetDirectory(0); h_costheta_Pmu_TRUE_RECO->Sumw2();
   
   TH2D* h_costheta_Pmu_Purity  = new TH2D("h_costheta_Pmu_Purity", "h_costheta_Pmu_Purity", cos_edges2D.size() - 1, cos_edges2D.data(), p_edges_2D.size() - 1, p_edges_2D.data());
    h_costheta_Pmu_Purity->SetDirectory(0); h_costheta_Pmu_Purity->Sumw2();

  std::map<int , BinMap> TH1Poly_binMap_test;
  std::map<int , BinMap> TH1Poly_binMap_test2;
  std::map<int , BinMap> TH1Poly_binMap_test3; 
  
   UBTH2Poly instance;
   gInterpreter->Declare("#include \"includes/UBTH2Poly.h\""); 
  
  
  UBTH2Poly *TH2Poly_costheta_Pmu_TRUE = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "TH2Poly_costheta_Pmu_TRUE");
  UBTH2Poly *TH2Poly_costheta_Pmu_incluive_TRUE = Make2DHist_inclusive_UB( MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test2, "TH2Poly_costheta_Pmu_incluive_TRUE");

  UBTH2Poly *TH2Poly_costheta_Pmu_TRUE2 = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "TH2Poly_costheta_Pmu_TRUE2");
  UBTH2Poly *TH2Poly_costheta_Pmu_incluive_TRUE2 = Make2DHist_inclusive_UB( MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test2, "TH2Poly_costheta_Pmu_incluive_TRUE2");

 UBTH2Poly *TH2Poly_costheta_Pmu_RECO2 = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "TH2Poly_costheta_Pmu_RECO2");
 UBTH2Poly *TH2Poly_costheta_Pmu_incluive_RECO2 = Make2DHist_inclusive_UB( MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test2, "TH2Poly_costheta_Pmu_incluive_RECO2");


  UBTH2Poly *TH2Poly_costheta_Pmu_TRUE_RECO = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "TH2Poly_costheta_Pmu_TRUE_RECO");
  UBTH2Poly *TH2Poly_costheta_Pmu_incluive_TRUE_RECO = Make2DHist_inclusive_UB( MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test2, "TH2Poly_costheta_Pmu_incluive_TRUE");

  UBTH2Poly *TH2Poly_costheta_Pmu = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "TH2Poly_costheta_Pmu");
  UBTH2Poly *TH2Poly_costheta_Pmu_incluive = Make2DHist_inclusive_UB( MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test2,"TH2Poly_costheta_Pmu_incluive");

  UBTH2Poly *TH2Poly_costheta_Pmu_Purity = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "TH2Poly_costheta_Pmu_Purity");
  UBTH2Poly *TH2Poly_costheta_Pmu_incluive_Purity = Make2DHist_inclusive_UB( MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test2, "TH2Poly_costheta_Pmu_incluive_Purity");

  TH2Poly_costheta_Pmu_TRUE->SetDirectory(0); TH2Poly_costheta_Pmu_TRUE->Sumw2();


HistFolio<UBTH2Poly, EventCategory> h_costheta_Pmu_UBTH2Poly_EventCategory =
HistFolio<UBTH2Poly, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_UBTH2Poly_EventCategory", MUON_2D_BIN_EDGES, true,"h_costheta_Pmu_EventCategory; Events");

HistFolio<UBTH2Poly, EventCategory> h_costheta_Pmu_UBTH2Poly_inclusive_EventCategory =
HistFolio<UBTH2Poly, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_UBTH2Poly_inclusive_EventCategory", MUON_2D_BIN_EDGES_inclusive, false,"h_costheta_Pmu_EventCategory; Events");



HistFolio<UBTH2Poly, EventCategory> h_costheta_Pmu_UBTH2Poly_Signal_EventCategory =
HistFolio<UBTH2Poly, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_UBTH2Poly_Signal_EventCategory", MUON_2D_BIN_EDGES, true,"h_costheta_Pmu_EventCategory; Events");

HistFolio<UBTH2Poly, EventCategory> h_costheta_Pmu_UBTH2Poly_inclusive_Signal_EventCategory =
HistFolio<UBTH2Poly, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_UBTH2Poly_inclusive_Signal_EventCategory", MUON_2D_BIN_EDGES_inclusive, false,"h_costheta_Pmu_EventCategory; Events");




HistFolio<UBTH2Poly, EventCategory> h_costheta_Pmu_UBTH2Poly_Contained_EventCategory =
HistFolio<UBTH2Poly, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_UBTH2Poly_Contained_EventCategory", MUON_2D_BIN_EDGES, true,"h_costheta_Pmu_EventCategory; Events");

HistFolio<UBTH2Poly, EventCategory> h_costheta_Pmu_UBTH2Poly_inclusive_Contained_EventCategory =
HistFolio<UBTH2Poly, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_UBTH2Poly_inclusive_Contained_EventCategory", MUON_2D_BIN_EDGES_inclusive, false,"h_costheta_Pmu_EventCategory; Events");


HistFolio<UBTH2Poly, EventCategory> h_costheta_Pmu_UBTH2Poly_NonContained_EventCategory =
HistFolio<UBTH2Poly, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_UBTH2Poly_NonContained_EventCategory", MUON_2D_BIN_EDGES, true,"h_costheta_Pmu_EventCategory; Events");

HistFolio<UBTH2Poly, EventCategory> h_costheta_Pmu_UBTH2Poly_inclusive_NonContained_EventCategory =
HistFolio<UBTH2Poly, EventCategory>(EventSelectionGroup_categories, "h_costheta_Pmu_UBTH2Poly_inclusive_NonContained_EventCategory", MUON_2D_BIN_EDGES_inclusive, false,"h_costheta_Pmu_EventCategory; Events");





UBTH2Poly* h_costheta_Pmu_UBTH2Poly  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly", "h_costheta_Pmu_UBTH2Poly",  MUON_2D_BIN_EDGES, true);
UBTH2Poly* h_costheta_Pmu_UBTH2Poly_inclusive  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly_inclusive", "h_costheta_Pmu_UBTH2Poly_inclusive",  MUON_2D_BIN_EDGES_inclusive, false);

UBTH2Poly* h_costheta_Pmu_UBTH2Poly_Signal  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly_Signal", "h_costheta_Pmu_UBTH2Poly_Signal",  MUON_2D_BIN_EDGES, true);
UBTH2Poly* h_costheta_Pmu_UBTH2Poly_inclusive_Signal  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly_inclusive_Signal", "h_costheta_Pmu_UBTH2Poly_inclusive_Signal",  MUON_2D_BIN_EDGES_inclusive, false);



UBTH2Poly* h_costheta_Pmu_UBTH2Poly_Contained  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly_Contained", "h_costheta_Pmu_UBTH2Poly_Contained",  MUON_2D_BIN_EDGES, true);
UBTH2Poly* h_costheta_Pmu_UBTH2Poly_inclusive_Contained  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly_inclusive_Contained", "h_costheta_Pmu_UBTH2Poly_inclusive_Contained",  MUON_2D_BIN_EDGES_inclusive, false);

UBTH2Poly* h_costheta_Pmu_UBTH2Poly_NonContained  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly_NonContained", "h_costheta_Pmu_UBTH2Poly_NonContained",  MUON_2D_BIN_EDGES, true);
UBTH2Poly* h_costheta_Pmu_UBTH2Poly_inclusive_NonContained  = new UBTH2Poly("h_costheta_Pmu_UBTH2Poly_inclusive_NonContained", "h_costheta_Pmu_UBTH2Poly_inclusive_NonContained",  MUON_2D_BIN_EDGES_inclusive, false);



  double res_max = -0.01;

TH2Poly * TH2Poly_costheta_Pmu_TRUE_test = new TH2Poly();

 TH1D* h_Cosmic_Cut  = new TH1D("h_Cosmic_Cut", "h_Cosmic_Cut",Cosmic_edges.size() - 1,Cosmic_edges.data() );
  h_Cosmic_Cut->SetDirectory(0);h_Cosmic_Cut->Sumw2();

 TH1D* h_Cosmic_Cut_Purity  = new TH1D("h_Cosmic_Cut_Purity", "h_Cosmic_Cut_Purity",Cosmic_edges.size() - 1,Cosmic_edges.data() );
  h_Cosmic_Cut_Purity->SetDirectory(0);h_Cosmic_Cut_Purity->Sumw2();

 TH1D* h_Cosmic_Cut_TRUE  = new TH1D("h_Cosmic_Cut_TRUE", "h_Cosmic_Cut_TRUE",Cosmic_edges.size() - 1,Cosmic_edges.data() );
  h_Cosmic_Cut_TRUE->SetDirectory(0);h_Cosmic_Cut_TRUE->Sumw2();
 
  TH1D* h_Cosmic_Cut_TRUE_RECO  = new TH1D("h_Cosmic_Cut_TRUE_RECO", "h_Cosmic_Cut_TRUE_RECO", Cosmic_edges.size() - 1,Cosmic_edges.data());
  h_Cosmic_Cut_TRUE_RECO->SetDirectory(0);h_Cosmic_Cut_TRUE_RECO->Sumw2();

//UBTH2Poly  * TH2Poly_costheta_Pmu_TRUE_test = new UBTH2Poly();
	//UBTH2Poly *h2p = new UBTH2Poly();
	
int nbin_mig = NBinsFromMap(MUON_2D_BIN_EDGES);
int nbins_inclusive = NBinsFromMap(MUON_2D_BIN_EDGES_inclusive);
	
	// ADDED ONE BIN TO INCLUDE ZERO , BUT number starts at 1 
	
TH2D* h_matr = new TH2D("h_matr", ";true bin; reconstructed bin; ", nbin_mig+1, 0, nbin_mig+2, nbin_mig+1, 0, nbin_mig+2);
TH2D* h_matr_inclusive = new TH2D("h_matr_inclusive", ";true bin; reconstructed bin; ", nbins_inclusive+1, 0, nbins_inclusive+2, nbins_inclusive+1, 0, nbins_inclusive+2);
	
//TKI_Hists TKI_Track1(EventSelectionGroup_categories,"TKI_Track1",1);
//TKI_Hists TKI_Track2(EventSelectionGroup_categories,"TKI_Track2",2);
//TKI_Hists TKI_Track3(EventSelectionGroup_categories,"TKI_Track3",3);
//TKI_Hists TKI_Track4(EventSelectionGroup_categories,"TKI_Track4",4);	
	
	
	

for (auto it = MUON_2D_BIN_EDGES.begin(); it != MUON_2D_BIN_EDGES.end(); ++it) {
        double currentKey = it->first;
        auto  Xvector = it->second;
        int vector_size = Xvector.size();   
        // Check if there is a next element
        auto nextIt = std::next(it);
        if (nextIt != MUON_2D_BIN_EDGES.end()) {
            double nextKey = nextIt->first;
            //if(nextKey ==2)continue;
            //std::string nextValue = nextIt->second;
             std::cout<<"~~~~~~~~~~~~~~~~~"<< std::endl;
            // Access current key, value, next key, and next value
            std::cout << "Current: Key=" << currentKey  << std::endl;
            std::cout << "Next: Key=" << nextKey  << std::endl;
            std::cout<<"Bins: ";
            for (int k = 0; k <vector_size-1; k++ )
            {std::cout<<Xvector.at(k)<< ", "<< Xvector.at(k+1)<< std::endl;
            double avg_X = .5 * (Xvector.at(k) +Xvector.at(k+1) );
            double avg_Y = .5 * (currentKey + nextKey); 
            
            BinMap BinMapSingle{currentKey,nextKey,Xvector.at(k),Xvector.at(k+1),avg_X,avg_Y};
            int binN = TH2Poly_costheta_Pmu_TRUE_test->AddBin(Xvector.at(k),currentKey, Xvector.at(k+1),nextKey);
            //std::cout<<"BinN = "<< binN<< std::endl;
                TH1Poly_binMap_test3.insert(std::pair<int,BinMap >(binN,BinMapSingle)); 
            }
            std::cout<< " "<< std::endl;
            
            
        } else {
            // If there is no next element, handle accordingly
            std::cout << "Finished last bin  element with Key=" << currentKey << std::endl;
        }
    
    std::cout<<"~~~~~~~~~~~~~~~~~"<< std::endl;
    
    
}


  for (auto& infile : filenames) {
    std::cout << "File: " << infile << std::endl;

   // TFile* f = TFile::Open(infile.c_str());
    
    TFile *f = new TFile(infile.c_str(), "READ");
    
    
    if (!(f && f->IsOpen() && !f->IsZombie())) {
        std::cout << "Bad file! " << infile << std::endl;
        continue;
    }
    TTree* stv  = (TTree*) f->Get("stv_tree");
    
    if (!(stv && stv->GetEntries())) {
       std::cout << "Bad file! " << infile << std::endl;
       continue;
    }





  bool mc_signal, sel_CC0pi, sel_muon_contained;
  bool sel_Muon_BDT_SCORE_above_threshold, sel_All_Protons_BDT_SCORE_above_threshold, sel_BDT_NotBOGUS;
  bool sel_BDT_predicts_1plusTrks_tobeProtons; 
  bool sel_cosmic_ip_cut_passed;
  float BTD_Muon_Candidate_prediction; 
  bool sel_muon_quality_ok; 
  
  TVector3 *p3_mu=nullptr, *mc_p3_mu=nullptr, *p3_mu_mcs=nullptr;
  TVector3 *p3_lead_p=nullptr, *mc_p3_lead_p=nullptr, *p3_pion=nullptr, *mc_p3_pion=nullptr;
  bool sel_signal_CCNP0Pi_only,sel_signal_CC1P0Pi_only,sel_signal_CC2P0Pi_only;
  int run, subrun; 
  
  std::vector<int> *mc_pdg = nullptr;
  std::vector<int> *pfpdg = nullptr;
  std::vector<int> *backtracked_pdg = nullptr;
  
  
  
  std::vector<int> *trk_BDT_PID = nullptr;
  std::vector<int> *trk_BDT_Probility = nullptr;
  std::vector<float> *trk_score_v = nullptr;
  std::vector<float> *trk_len_v = nullptr;
  std::vector<float> *trk_distance_v = nullptr;
  
  
std::vector<float> * P_p_hadronTk= nullptr;
std::vector<float> * P_costheta_hadronTk= nullptr;
std::vector<float> * delta_pT_hadronTk= nullptr;
std::vector<float> * delta_phiT_hadronTk= nullptr;
std::vector<float> * delta_alphaT_hadronTk= nullptr;
std::vector<float> * delta_pL_hadronTk= nullptr;
std::vector<float> * pn_hadronTk= nullptr;
std::vector<float> * delta_pTx_hadronTk= nullptr;
std::vector<float> * delta_pTy_hadronTk= nullptr;
std::vector<float> * deltaPhi3D_hadronTk= nullptr;
std::vector<float> * deltaAlpha3DMu_hadronTk= nullptr;
  
  
std::vector<float> * P_p_hadronTk_mc= nullptr;
std::vector<float> * P_costheta_hadronTk_mc= nullptr;
std::vector<float> * delta_pT_hadronTk_mc= nullptr;
std::vector<float> * delta_phiT_hadronTk_mc= nullptr;
std::vector<float> * delta_alphaT_hadronTk_mc= nullptr;
std::vector<float> * pn_hadronTk_mc= nullptr;
std::vector<float> * delta_pTx_hadronTk_mc= nullptr;
std::vector<float> * delta_pTy_hadronTk_mc= nullptr;  
  
  float  spline_weight, tuned_cv_weight;
  float pn, delta_alphaT, delta_pTx, delta_pTy, delta_pL, delta_phiT, delta_pT;
  float reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z;
  float mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z;
  float mc_pn, mc_delta_pL, mc_delta_alphaT, mc_delta_pTx, mc_delta_pTy, mc_delta_p, mc_delta_phiT, mc_delta_pT;
  float topological_score;
  float muon_candidate_BDT;
  int mc_interaction,proton_multiplicity,pion_multiplicity,category;  
  int index_muon, index_leadproton; 
  int n_pfps;
  //float cosmic_impact_parameter_;
  float CosmicIP;
  std::vector<float> *trk_BDT_score_0 = nullptr;
  std::vector<float> *trk_BDT_score_1 = nullptr;
  std::vector<float> *trk_BDT_score_2 = nullptr;
  std::vector<float> *trk_BDT_score_3 = nullptr;

  std::vector<int> *pfp_generation_v = nullptr;



bool sel_reco_vertex_in_FV , sel_pfp_starts_in_PCV ;
bool sel_has_muon_candidate,  sel_topo_cut_passed  ;
bool sel_no_reco_showers , sel_muon_passed_mom_cuts  ;
bool sel_pions_above_threshold , sel_has_pi_candidate ;


//int pion_multiplicity; 




  std::vector<TVector3> *Proton_RECO_3vector = nullptr;
  std::vector<TVector3> *mc_Proton_RECO_3vector = nullptr;


  stv->SetBranchStatus("*", 0); // I think this turns all the branches to false 
  stv->SetBranchStatus("mc_is_cc0pi_signal", 1);
  stv->SetBranchStatus("sel_CC0pi", 1);
  stv->SetBranchStatus("mc_p3_mu", 1);
  

  stv->SetBranchStatus("p3_mu", 1);
  stv->SetBranchStatus("p3_mu_mcs", 1);
  stv->SetBranchStatus("mc_p3_lead_p", 1);
  stv->SetBranchStatus("mc_p3_lead_pion", 1);
  stv->SetBranchStatus("p3_lead_pion", 1);
  stv->SetBranchStatus("pfpdg",1);
  

  
  stv->SetBranchStatus("muon_candidate_idx", 1);
  stv->SetBranchStatus("lead_p_candidate_idx", 1);
  stv->SetBranchAddress("muon_candidate_idx", &index_muon);
  stv->SetBranchAddress("lead_p_candidate_idx", &index_leadproton);
  
  stv->SetBranchStatus("pfp_generation_v", 1);
  stv->SetBranchAddress("pfp_generation_v", &pfp_generation_v);
  
  
  
  stv->SetBranchStatus("p3_lead_p", 1);
  stv->SetBranchStatus("p3_p_vec", 1);
  stv->SetBranchAddress("p3_p_vec", &Proton_RECO_3vector);
  
  stv->SetBranchStatus("mc_p3_p_vec", 1);
  stv->SetBranchAddress("mc_p3_p_vec", &mc_Proton_RECO_3vector);
  

  
  stv->SetBranchStatus("mc_pdg", 1);
  stv->SetBranchStatus("sel_muon_contained", 1);
  stv->SetBranchStatus("spline_weight", 1);
  stv->SetBranchStatus("tuned_cv_weight", 1);
  stv->SetBranchStatus("CosmicIP", 1);

  
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
  stv->SetBranchStatus("category",1);
  stv->SetBranchAddress("category", &category);


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




  
  
  
  
  
  stv->SetBranchStatus("mc_nu_vtx_y",1);
  stv->SetBranchStatus("mc_nu_vtx_z",1); 

  stv->SetBranchStatus("run",1);
  stv->SetBranchAddress("run", &run);
  stv->SetBranchStatus("subrun",1); 
  stv->SetBranchAddress("subrun", &subrun);
  

  stv->SetBranchStatus("topological_score",1); 
  stv->SetBranchAddress("topological_score", &topological_score);


  stv->SetBranchStatus("trk_score_v",1); 
  stv->SetBranchAddress("trk_score_v", &trk_score_v);
  stv->SetBranchStatus("trk_len_v",1); 
  stv->SetBranchAddress("trk_len_v", &trk_len_v);
  
  stv->SetBranchStatus("trk_distance_v",1); 
  stv->SetBranchAddress("trk_distance_v", &trk_distance_v);
  
//
//stv->SetBranchStatus("P_costheta_hadronTk",1);  
stv->SetBranchStatus("delta_pT_hadronTk",1);
stv->SetBranchStatus("delta_phiT_hadronTk",1);
stv->SetBranchStatus("delta_alphaT_hadronTk",1);
stv->SetBranchStatus("delta_pL_hadronTk",1);
stv->SetBranchStatus("pn_hadronTk",1);
stv->SetBranchStatus("delta_pTx_hadronTk",1);
stv->SetBranchStatus("delta_pTy_hadronTk",1);
stv->SetBranchStatus("deltaPhi3D_hadronTk",1);
stv->SetBranchStatus("deltaAlpha3DMu_hadronTk",1);

//stv->SetBranchAddress("P_costheta_hadronTk",&P_costheta_hadronTk);
//stv->SetBranchAddress("P_costheta_hadronTk",&P_costheta_hadronTk);
stv->SetBranchAddress("delta_pT_hadronTk",&delta_pT_hadronTk);
stv->SetBranchAddress("delta_phiT_hadronTk",&delta_phiT_hadronTk);
stv->SetBranchAddress("delta_alphaT_hadronTk",&delta_alphaT_hadronTk);
stv->SetBranchAddress("delta_pL_hadronTk",&delta_pL_hadronTk);
stv->SetBranchAddress("pn_hadronTk",&pn_hadronTk);
stv->SetBranchAddress("delta_pTx_hadronTk",&delta_pTx_hadronTk);
stv->SetBranchAddress("delta_pTy_hadronTk",&delta_pTy_hadronTk);
stv->SetBranchAddress("deltaPhi3D_hadronTk",&deltaPhi3D_hadronTk);
stv->SetBranchAddress("deltaAlpha3DMu_hadronTk",&deltaAlpha3DMu_hadronTk);

stv->SetBranchStatus("mc_interaction",1); 
stv->SetBranchAddress("mc_interaction", &mc_interaction);
stv->SetBranchAddress("mc_is_cc0pi_signal", &mc_signal);
stv->SetBranchAddress("sel_CC0pi", &sel_CC0pi); 

stv->SetBranchStatus("sel_muon_quality_ok",1);
stv->SetBranchAddress("sel_muon_quality_ok", &sel_muon_quality_ok); 
  
stv->SetBranchStatus("sel_Muon_BDT_SCORE_above_threshold", 1);
stv->SetBranchStatus("sel_All_Protons_BDT_SCORE_above_threshold", 1);
stv->SetBranchStatus("sel_BDT_predicts_1plusTrks_tobeProtons", 1);
  
stv->SetBranchStatus("sel_BDT_NotBOGUS", 1);
stv->SetBranchAddress("sel_Muon_BDT_SCORE_above_threshold", &sel_Muon_BDT_SCORE_above_threshold);
stv->SetBranchAddress("sel_BDT_predicts_1plusTrks_tobeProtons", &sel_BDT_predicts_1plusTrks_tobeProtons);
  
stv->SetBranchAddress("sel_All_Protons_BDT_SCORE_above_threshold", &sel_All_Protons_BDT_SCORE_above_threshold);
stv->SetBranchAddress("sel_BDT_NotBOGUS", &sel_BDT_NotBOGUS);
stv->SetBranchStatus("BTD_Muon_Candidate_prediction", 1);
stv->SetBranchAddress("BTD_Muon_Candidate_prediction", &BTD_Muon_Candidate_prediction);

stv->SetBranchStatus("sel_cosmic_ip_cut_passed", 1);
stv->SetBranchAddress("sel_cosmic_ip_cut_passed", &sel_cosmic_ip_cut_passed);
  
stv->SetBranchAddress("sel_CCNp0pi", &sel_signal_CCNP0Pi_only);
stv->SetBranchAddress("sel_CC1p0pi", &sel_signal_CC1P0Pi_only);
stv->SetBranchAddress("sel_CC2p0pi", &sel_signal_CC2P0Pi_only);

 //sel_CCNp0pi
  stv->SetBranchAddress("p3_mu", &p3_mu);
  stv->SetBranchAddress("p3_mu_mcs", &p3_mu_mcs);
  stv->SetBranchAddress("p3_lead_pion", &p3_pion);
  stv->SetBranchAddress("mc_p3_lead_pion", &mc_p3_pion);
  stv->SetBranchAddress("mc_p3_mu", &mc_p3_mu);
  stv->SetBranchAddress("mc_p3_lead_p", &mc_p3_lead_p);
  stv->SetBranchAddress("p3_lead_p", &p3_lead_p);
  
  stv->SetBranchAddress("mc_pdg", &mc_pdg);
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
 stv->SetBranchAddress("CosmicIP", &CosmicIP);



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
  stv->SetBranchAddress("mc_interaction", &mc_interaction); 
  
//stv->SetBranchStatus("mc_P_p_hadronTk",1);
//stv->SetBranchStatus("mc_P_costheta_hadronTk",1);    
stv->SetBranchStatus("mc_delta_pT_hadronTk",1);
stv->SetBranchStatus("mc_delta_phiT_hadronTk",1);
stv->SetBranchStatus("mc_delta_alphaT_hadronTk",1);
stv->SetBranchStatus("mc_delta_pL_hadronTk",1);
stv->SetBranchStatus("mc_pn_hadronTk",1);
stv->SetBranchStatus("mc_delta_pTx_hadronTk",1);
stv->SetBranchStatus("mc_delta_pTy_hadronTk",1);

//stv->SetBranchAddress("mc_P_p_hadronTk",&P_p_hadronTk_mc);
//stv->SetBranchAddress("mc_P_costheta_hadronTk",&P_costheta_hadronTk_mc);
stv->SetBranchAddress("mc_delta_pT_hadronTk",&delta_pT_hadronTk_mc);
stv->SetBranchAddress("mc_delta_phiT_hadronTk",&delta_phiT_hadronTk_mc);
stv->SetBranchAddress("mc_delta_alphaT_hadronTk",&delta_alphaT_hadronTk_mc);
stv->SetBranchAddress("mc_pn_hadronTk",&pn_hadronTk_mc);
stv->SetBranchAddress("mc_delta_pTx_hadronTk",&delta_pTx_hadronTk_mc);
stv->SetBranchAddress("mc_delta_pTy_hadronTk",&delta_pTy_hadronTk_mc);


//stv->SetBranchStatus("deltaPhi3D_hadronTk_mc",1);
//stv->SetBranchStatus("deltaAlpha3DMu_hadronTk_mc",1);
  
  

  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout<< "~~~~~~~~~~~~~~~ STarting Loop~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout << " The Number of Entries "<<  stv->GetEntries() << std::endl;
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;

  int zz = 1; 
  for (long i=0; i<stv->GetEntries(); i++)
  {
    stv->GetEntry(i);
   // if (sel_signal && mc_signal){
     if ( i % 10000 == 0 ) {
      std::cout << "Processing event #" << i << '\n';
    }

 /*&& BTD_Muon_Candidate_prediction != 9999 
                && !std::isnan(BTD_Muon_Candidate_prediction) 
                && !std::isinf(BTD_Muon_Candidate_prediction) */

//if(sel_All_Protons_BDT_SCORE_above_threshold==false){std::cout<< " this is false"<< std::endl;}

//BTD_Muon_Candidate_prediction != 9999 
if (  sel_CC0pi  ){
// && sel_BDT_NotBOGUS

//sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV &&
//sel_has_muon_candidate && sel_topo_cut_passed && 
//sel_no_reco_showers && sel_muon_passed_mom_cuts && 
//sel_pions_above_threshold && sel_has_pi_candidate && pion_multiplicity > 0
//&& sel_BDT_predicts_1plusTrks_tobeProtons 
  //std::cout<<"sel_CC0pi && sel_BDT_NotBOGUS"<< std::endl;
  //&& sel_All_Protons_BDT_SCORE_above_threshold
  //
     //mc_is_cc0pi_signal && && mc_signal
     //std::cout<<"pn = " << pn << std::endl;
     bool ProtonGood= true; 
      EventCategory TrueEventType = EventCategory_tool.returnEventCategoryType(category);
      bool IsSignal = EventCategory_tool.IsSignal(TrueEventType); 
      FidVol_Category Containmenttype = EventCategory_tool.GetContainmentType(sel_muon_contained,IsSignal);
      CCZeroPi_type TopologyType = EventCategory_tool.returnCCZeroPi_type(proton_multiplicity,IsSignal);
      
      double Wgtt = spline_weight * tuned_cv_weight; // *0.124135
      Double_t Wgt = safe_weight_histMaker( Wgtt );
      double Pmu = p3_mu->Mag();
      double Pmu_true = mc_p3_mu->Mag();
      Double_t Pmu2 = p3_mu->Mag();
      //double Pmu_Correction = Pmu_PANOS_Correction(Pmu, sel_muon_contained);
      
      double Pmu_mcs = p3_mu_mcs->Mag();
     
     
      double Pmu_Resolution = (Pmu_true - Pmu);
      double Pmu_CMS_Resolution = (Pmu_true - Pmu_mcs);
      //std::cout<<"here 1"<< std::endl;
      double CosTheta = p3_mu->CosTheta();
      Double_t CosTheta2 = p3_mu->CosTheta();
      double CosTheta_true = mc_p3_mu->CosTheta();
      double CosTheta_MCS = p3_mu_mcs->CosTheta();
      //Double_t CosTheta_true2 = CosTheta_true;

      double CosTheta_Resolution =CosTheta - CosTheta_true;
      double CosTheta_CMS_Resolution = CosTheta_MCS - CosTheta_true;
      double P_proton = p3_lead_p->Mag();
      double P_CosTheta = p3_lead_p->CosTheta();
      
      double OpenAngle = TMath::ACos( (p3_mu->X()*p3_lead_p->X() + p3_mu->Y()*p3_lead_p->Y() + p3_mu->Z()*p3_lead_p->Z()) / p3_mu->Mag());
      double OpenAngle_deg = OpenAngle * TMath::RadToDeg();
      double mc_OpenAngle = TMath::ACos( (mc_p3_mu->X()*mc_p3_lead_p->X() + mc_p3_mu->Y()*p3_lead_p->Y() + mc_p3_mu->Z()*mc_p3_lead_p->Z()) / mc_p3_mu->Mag());
      double mc_OpenAngle_deg = mc_OpenAngle* TMath::RadToDeg();
      double OpenAngle_resolution = OpenAngle - mc_OpenAngle;
      int NTracks = proton_multiplicity + 1; // assuming the number of tracks is the Proton Multipleicity + muon track 
      float CosmicIP_ = CosmicIP;
      
      double PmuCorrection = Pmu_PANOS_Correction(Pmu, sel_muon_contained, sel_muon_quality_ok);
      double Ponos_Pmu = Pmu+PmuCorrection;
      
      graph_PonosFunction->SetPoint(zz, Pmu, PmuCorrection);
      zz++;      

      h_p->Fill(Pmu,Wgt); 
      h_p_Correction->Fill(Ponos_Pmu,Wgt);
      h_p_true->Fill(Pmu_true,Wgt); 
      h_p_Resolution->Fill(Pmu_Resolution,Wgt); 
      
      h_p_MCS_Resolution->Fill(Pmu_CMS_Resolution,Wgt); 
      h_p_Mig->Fill(Pmu_true,Pmu,Wgt); 
      h_p_MCS_Mig->Fill(Pmu_true,Pmu_mcs,Wgt); 
      h_Proton_mult->Fill(proton_multiplicity,Wgt); 
      h_Cosmic_Cut->Fill(CosmicIP_,Wgt); 
      //std::cout<<" CosmicIP = "<<CosmicIP_<<  std::endl;

      
      
      int TrueProtonN = Num_MatchedPDGs(*backtracked_pdg,2212);
      
      if(TrueProtonN ==proton_multiplicity ){
      h_Proton_mult_Purity->Fill(proton_multiplicity,Wgt);
      
      }
      
      if(sel_BDT_NotBOGUS == false){
      std::string Bogisinfo = std::to_string(i); 
      
         printText(Bogisinfo, "IsBugus_Enties.txt");
      
      
      }
      
      
      //std::cout<<"here 3"<< std::endl;
      
      h_Proton_mult_EventCategory.GetComponentHist(TrueEventType)->Fill(proton_multiplicity,Wgt);
      
      h_Pion_mult->Fill(pion_multiplicity,Wgt); 
            
     h_CCPiCategory->AddBinContent(1,Wgt);
     if(sel_signal_CCNP0Pi_only==true) {h_CCPiCategory->AddBinContent(2,Wgt);}
     if(sel_signal_CC1P0Pi_only ==true){h_CCPiCategory->AddBinContent(3,Wgt);}
     if(sel_signal_CC2P0Pi_only==true) {h_CCPiCategory->AddBinContent(4,Wgt);}

     h_costheta_p->Fill(CosTheta,Wgt);
     h_costheta_p_true->Fill(CosTheta_true,Wgt);
     h_costheta_p_Resolution->Fill(CosTheta_Resolution,Wgt); 
     h_costheta_p_MCS_Resolution->Fill(CosTheta_CMS_Resolution,Wgt);
     
     h_costheta_p_Mig->Fill(CosTheta_true,CosTheta,Wgt); 
     h_costheta_p_MCS_Mig->Fill(CosTheta_true, CosTheta_MCS,Wgt); 
     h_p_EventCategory.GetComponentHist(TrueEventType)->Fill(Pmu,Wgt);
     h_p_containment_categories.GetComponentHist(Containmenttype)->Fill(Pmu,Wgt);
     h_p_topology_categories.GetComponentHist(TopologyType)->Fill(Pmu,Wgt);
     
     
     h_p_Correction_EventCategory.GetComponentHist(TrueEventType)->Fill(Ponos_Pmu,Wgt);
     h_p_Correction_containment_categories.GetComponentHist(Containmenttype)->Fill(Ponos_Pmu,Wgt);
     h_p_Correction_topology_categories.GetComponentHist(TopologyType)->Fill(Ponos_Pmu,Wgt);
     
     h_costheta_p_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Wgt);
     h_costheta_p_containment_categories.GetComponentHist(Containmenttype)->Fill(CosTheta,Wgt);
     h_costheta_p_topology_categories.GetComponentHist(TopologyType)->Fill(CosTheta,Wgt);
     
     h_NTracks_EventCategory.GetComponentHist(TrueEventType)->Fill(NTracks,Wgt);
     h_p_costheta_p->Fill(CosTheta,Pmu,Wgt); 
 
 
    h_costheta_Pmu_UBTH2Poly_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt); 
    h_costheta_Pmu_UBTH2Poly_inclusive_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt); 
 
  
    h_costheta_Pmu_UBTH2Poly->Fill(CosTheta,Pmu,Wgt);
    h_costheta_Pmu_UBTH2Poly_inclusive->Fill(CosTheta,Pmu,Wgt); 

 
 if(sel_muon_contained==true){
    h_costheta_Pmu_UBTH2Poly_Contained_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt); 
    h_costheta_Pmu_UBTH2Poly_inclusive_Contained_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt); 
    h_costheta_Pmu_UBTH2Poly_Contained->Fill(CosTheta,Pmu,Wgt);
    h_costheta_Pmu_UBTH2Poly_inclusive_Contained->Fill(CosTheta,Pmu,Wgt); 
 }
 else if (sel_muon_contained==false){
    h_costheta_Pmu_UBTH2Poly_NonContained_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt); 
    h_costheta_Pmu_UBTH2Poly_inclusive_NonContained_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt); 
    h_costheta_Pmu_UBTH2Poly_NonContained->Fill(CosTheta,Pmu,Wgt);
    h_costheta_Pmu_UBTH2Poly_inclusive_NonContained->Fill(CosTheta,Pmu,Wgt); 
 }
 
if(IsSignal==true){

h_costheta_Pmu_UBTH2Poly_Signal->Fill(CosTheta,Pmu,Wgt);
h_costheta_Pmu_UBTH2Poly_inclusive_Signal->Fill(CosTheta,Pmu,Wgt); 
h_costheta_Pmu_UBTH2Poly_Signal_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt); 
h_costheta_Pmu_UBTH2Poly_inclusive_Signal_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt); 
   
}

   




    if(proton_multiplicity ==0){
    /* IF its NP*/
       h_p_CC0P0Pi_Only->Fill(Pmu,Wgt); 
       h_p_CC0P0Pi_Only_true->Fill(Pmu_true,Wgt); 
       h_p_CC0P0Pi_Only_Resolution->Fill(Pmu_Resolution,Wgt); 
       h_p_CC0P0Pi_Only_Mig->Fill(Pmu_true,Pmu,Wgt);
       h_p_CC0P0Pi_Only_EventCategory.GetComponentHist(TrueEventType)->Fill(Pmu,Wgt);
       h_p_CC0P0Pi_Only_containment_categories.GetComponentHist(Containmenttype)->Fill(Pmu,Wgt);
       h_costheta_p_CC0P0Pi_Only->Fill(CosTheta,Wgt);
       h_costheta_p_CC0P0Pi_Only_true->Fill(CosTheta_true,Wgt);
       h_costheta_p_CC0P0Pi_Only_Resolution->Fill(CosTheta_Resolution,Wgt); 
       h_costheta_p_CC0P0Pi_Only_Mig->Fill(CosTheta_true,CosTheta,Wgt); 
       h_costheta_p_CC0P0Pi_Only_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Wgt);
       h_costheta_p_CC0P0Pi_Only_containment_categories.GetComponentHist(Containmenttype)->Fill(CosTheta,Wgt);
       
       

       }

   else {
         h_costheta_p_CCNP0Pi_Only->Fill(CosTheta,Wgt);
         h_costheta_p_CCNP0Pi_Only_true->Fill(CosTheta_true,Wgt);
         h_costheta_p_CCNP0Pi_Only_Resolution->Fill(CosTheta_Resolution,Wgt); 
         h_costheta_p_CCNP0Pi_Only_Mig->Fill(CosTheta_true,CosTheta,Wgt); 
         h_costheta_p_CCNP0Pi_Only_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Wgt);
         h_costheta_p_CCNP0Pi_Only_containment_categories.GetComponentHist(Containmenttype)->Fill(CosTheta,Wgt);
         h_p_CCNP0Pi_Only->Fill(Pmu,Wgt); 
         h_p_CCNP0Pi_Only_true->Fill(Pmu_true,Wgt); 
         h_p_CCNP0Pi_Only_Resolution->Fill(Pmu_Resolution,Wgt); 
         h_p_CCNP0Pi_Only_Mig->Fill(Pmu_true,Pmu,Wgt);
         h_p_CCNP0Pi_Only_EventCategory.GetComponentHist(TrueEventType)->Fill(Pmu,Wgt);
         h_p_CCNP0Pi_Only_containment_categories.GetComponentHist(Containmenttype)->Fill(Pmu,Wgt);

       }


    std::vector<int> BDTPID_Prediction_Muon_front = MakereorderVector(*trk_BDT_PID, BDT_Category::kBDT_Muon);

   int TrackN = 0;
   int TrackN_proton = 1;




   for(auto PID : BDTPID_Prediction_Muon_front){

    int true_pdg = backtracked_pdg->at(TrackN);
    Particle_type REduced_Particle =  EventCategory_tool.GetParticlegroup_typeReduced(true_pdg);
    BDT_Category BDT_Type = EventCategory_tool.returnBDTPredictionType(PID);

   if(PID==BDT_Category::kBDT_Pion && REduced_Particle == kPion_pos_neg){

   double p3_Pion = p3_pion->Mag();
   double mc_p3_Pion = mc_p3_pion->Mag();
   h_pionE_Mig->Fill(mc_p3_Pion,p3_Pion,Wgt);
 }

  if(PID==BDT_Category::kBDT_Proton){
   h_Proton_mult_Tracks_Particle_categories.GetComponentHist(REduced_Particle)->Fill(TrackN_proton,Wgt);
   h_Proton_mult_Tracks_EventCategory.GetComponentHist(TrueEventType)->Fill(TrackN_proton,Wgt);
   h_Proton_mult_Tracks->Fill(TrackN_proton,Wgt);
   
   if (true_pdg==2212){
        h_Proton_mult_Tracks_Purity->Fill(TrackN_proton,Wgt);
   }
  TrackN_proton++;  
   }


 
  if(PID==BDT_Category::kBDT_Else){
  
   h_Else_mult_Tracks_Particle_categories.GetComponentHist(REduced_Particle)->Fill(TrackN,Wgt);
   h_Else_mult_Tracks_EventCategory.GetComponentHist(TrueEventType)->Fill(TrackN,Wgt);
  //std::cout<< " run = "<< run << " subrun = "<< subrun<< std::endl; 
  std::string interactiontype =  EventCategory_tool.label( TrueEventType );
  
  std::string inputelse = "run,  " + std::to_string(run ) + " , subrun, " + std::to_string(subrun) + " interaction, " + interactiontype;
  
  
   printText(inputelse, "RUB_subrun.txt");
  }
 
 
 
 
 
 h_NTracks_BTDGroup_categories.GetComponentHist(BDT_Type)->Fill(TrackN,Wgt);
 h_NTracks_Particle_categories.GetComponentHist(REduced_Particle)->Fill(TrackN,Wgt);
 
 h_NTracks->Fill(TrackN,Wgt);
 TrackN++;
  }/// End of BDTPID_Prediction_Muon_front
 //////////////////////////////////////////
 bool ThereISaElse = false; 



 for(int index1 = 0 ; index1 <pfpdg->size(); index1++ ){

 double score_0 = trk_BDT_score_0->at(index1);
 double score_1 = trk_BDT_score_1->at(index1);
 double score_2 = trk_BDT_score_2->at(index1);
 double score_3 = trk_BDT_score_3->at(index1);
 int PIDBDT = trk_BDT_PID->at(index1);

  unsigned int generation = pfp_generation_v->at( index1 );
  if ( generation != 2u ) continue;

 int true_pdg = backtracked_pdg->at(index1);
 Particle_type REduced_Particle =  EventCategory_tool.GetParticlegroup_typeReduced(true_pdg);
 
 
 h_trk_distance_v->Fill(trk_distance_v->at(index1),Wgt);
 h_trk_distance_v_Particle_categories.GetComponentHist(REduced_Particle)->Fill(trk_distance_v->at(index1),Wgt);
  h_trk_distance_v_EventCategory.GetComponentHist(TrueEventType)->Fill(trk_distance_v->at(index1),Wgt);


 h_trk_score_v->Fill(trk_score_v->at(index1),Wgt);
 h_trk_score_v_Particle_categories.GetComponentHist(REduced_Particle)->Fill(trk_score_v->at(index1),Wgt);
 h_trk_score_v_EventCategory.GetComponentHist(TrueEventType)->Fill(trk_score_v->at(index1),Wgt);


 if(index1 == index_muon){
  // Lead Muon Prediction
   h_trk_len_v->Fill(trk_len_v->at(index1),Wgt);
 h_trk_len_v_Particle_categories.GetComponentHist(REduced_Particle)->Fill(trk_len_v->at(index1),Wgt);
 h_trk_len_v_EventCategory.GetComponentHist(TrueEventType)->Fill(trk_len_v->at(index1),Wgt);
  h_topological_score->Fill(topological_score,Wgt);
  h_topological_score_Particle_categories.GetComponentHist(REduced_Particle)->Fill(topological_score,Wgt);
  h_topological_score_EventCategory.GetComponentHist(TrueEventType)->Fill(topological_score,Wgt);

  h_CCZeroP_MuonCandidate_BTDPrediction->Fill(score_1,Wgt);
  h_CCZeroP_MuonCandidate_BTDPrediction_Particle_categories.GetComponentHist(REduced_Particle)->Fill(score_1,Wgt);
  h_CCZeroP_MuonCandidate_BTDPrediction_EventCategory.GetComponentHist(TrueEventType)->Fill(score_1,Wgt);
    
  h_CCZeroP_MuonCandidate_BTDPrediction_topology_categories.GetComponentHist(TopologyType)->Fill(score_1,Wgt);
 h_CCZeroP_MuonCandidate_BTDPrediction_containment_categories.GetComponentHist(Containmenttype)->Fill(score_1,Wgt);
  
  h_VertexX_leadmuon_Particle_categories.GetComponentHist(REduced_Particle)->Fill(reco_nu_vtx_sce_x,Wgt);
  h_VertexY_leadmuon_Particle_categories.GetComponentHist(REduced_Particle)->Fill(reco_nu_vtx_sce_y,Wgt);
  h_VertexZ_leadmuon_Particle_categories.GetComponentHist(REduced_Particle)->Fill(reco_nu_vtx_sce_z,Wgt);
  
  
    if(mc_signal== true){
       h_trk_len_v_Purity->Fill(trk_len_v->at(index1),Wgt);
    }
  
  
  if(true_pdg==13){
  h_CCZeroP_MuonCandidate_BTDPrediction_purity->Fill(score_1,Wgt);
  }
 
 
 }
 
 else{

 }

 
 if (index1 != index_muon){
  if(PIDBDT==BDT_Category::kBDT_Else){
    h_CCZeroP_Else_BTDPrediction_Particle_categories.GetComponentHist(REduced_Particle)->Fill(score_0,Wgt);
    h_CCZeroP_Else_BTDPrediction_EventCategory.GetComponentHist(TrueEventType)->Fill(score_0,Wgt);
    
    h_CCZeroP_Else_BTDPrediction->Fill(score_0,Wgt);
    h_Else_Score_Pmu->Fill(score_0,Wgt);
    h_Else_Score_CosTheta->Fill(score_0,Wgt);
    
    ThereISaElse = true; 
    
    //std::cout<<"run :"<< run<< " subrun = "<< subrun << std::endl;
    
    
    
   }
   
  else if(PIDBDT==BDT_Category::kBDT_Proton){
    h_CCZeroP_ALLProtonCandidate_BTDPrediction->Fill(score_3,Wgt);
    h_CCZeroP_ProtonCandidate_BTDPrediction->Fill(score_3,Wgt);
    h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories.GetComponentHist(REduced_Particle)->Fill(score_3,Wgt);
    h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory.GetComponentHist(TrueEventType)->Fill(score_3,Wgt);
    h_CCZeroP_ProtonCandidate_BTDPrediction_topology_categories.GetComponentHist(TopologyType)->Fill(score_3,Wgt);
    h_CCZeroP_ProtonCandidate_BTDPrediction_containment_categories.GetComponentHist(Containmenttype)->Fill(score_3,Wgt);
 
     if(true_pdg == 2212) h_CCZeroP_ALLProtonCandidate_BTDPrediction_Purity->Fill(score_3,Wgt);
  
   }
 
 }// end if Not Muon Candaite
 
 
 
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
  else if(index1==2){
  h_CC1P_BTDGroup_Probability_tobeMuon_tk2->Fill(score_1,Wgt);
  h_CC1P_BTDGroup_Probability_tobeProton_tk2->Fill(score_3,Wgt);
  h_CC1P_BTDGroup_Probability_tobePion_tk2->Fill(score_2,Wgt);
  h_CC1P_BTDGroup_Probability_tobeElse_tk2->Fill(score_0,Wgt);
  }
 }
 else if (proton_multiplicity ==2){
   if(index1==1){
   h_CC2P_BTDGroup_Probability_tobeMuon_tk1->Fill(score_1,Wgt);
   h_CC2P_BTDGroup_Probability_tobeProton_tk1->Fill(score_3,Wgt);
   h_CC2P_BTDGroup_Probability_tobePion_tk1->Fill(score_2,Wgt);
   h_CC2P_BTDGroup_Probability_tobeElse_tk1->Fill(score_0,Wgt);
   }
 
  if(index1==2){
  h_CC2P_BTDGroup_Probability_tobeMuon_tk2->Fill(score_1,Wgt);
  h_CC2P_BTDGroup_Probability_tobeProton_tk2->Fill(score_3,Wgt);
  h_CC2P_BTDGroup_Probability_tobePion_tk2->Fill(score_2,Wgt);
  h_CC2P_BTDGroup_Probability_tobeElse_tk2->Fill(score_0,Wgt);
  }
 
   if(index1==3){
  h_CC2P_BTDGroup_Probability_tobeMuon_tk3->Fill(score_1,Wgt);
  h_CC2P_BTDGroup_Probability_tobeProton_tk3->Fill(score_3,Wgt);
  h_CC2P_BTDGroup_Probability_tobePion_tk3->Fill(score_2,Wgt);
  h_CC2P_BTDGroup_Probability_tobeElse_tk3->Fill(score_0,Wgt);
   }
  }
  
  
  
  
  
  
  
  
 }//// END of pfpdg Loop 
 ////////////////////////////////////////////
  h_costheta_Pmu ->Fill(CosTheta,Pmu,Wgt);  
  h_VertexX->Fill(reco_nu_vtx_sce_x,Wgt);  
  h_VertexX_true->Fill(mc_nu_vtx_x,Wgt);
  h_VertexX_Resolution->Fill(mc_nu_vtx_x-reco_nu_vtx_sce_x,Wgt);
  h_VertexX_Mig->Fill(mc_nu_vtx_x,reco_nu_vtx_sce_x,Wgt); 
  h_VertexY->Fill(reco_nu_vtx_sce_y,Wgt);  
  h_VertexY_true->Fill(mc_nu_vtx_y,Wgt);
  h_VertexY_Resolution->Fill(mc_nu_vtx_y-reco_nu_vtx_sce_y,Wgt);
  h_VertexY_Mig->Fill(mc_nu_vtx_y,reco_nu_vtx_sce_y,Wgt); 
  h_VertexZ->Fill(reco_nu_vtx_sce_z,Wgt);  
  h_VertexZ_true->Fill(mc_nu_vtx_z,Wgt);
  h_VertexZ_Resolution->Fill(mc_nu_vtx_z-reco_nu_vtx_sce_z,Wgt);
  h_VertexZ_Mig->Fill(mc_nu_vtx_z,reco_nu_vtx_sce_z,Wgt); 
  h_VertexX_Z->Fill(reco_nu_vtx_sce_z, reco_nu_vtx_sce_x,Wgt); 
  h_VertexY_Z->Fill(reco_nu_vtx_sce_z, reco_nu_vtx_sce_y,Wgt); 
  h_VertexX_Y->Fill(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y,Wgt); 
  
  h_costheta_Pmu_EventCategory.GetComponentHist(TrueEventType)->Fill(CosTheta,Pmu,Wgt);
  
    
  TH2Poly_costheta_Pmu->Fill(CosTheta2,Pmu2,Wgt); 
  TH2Poly_costheta_Pmu_incluive->Fill(CosTheta2,Pmu2,Wgt); 
  
  if(mc_signal== true){
  h_costheta_Pmu_Purity->Fill(CosTheta,Pmu,Wgt);
  h_Cosmic_Cut_Purity->Fill(CosmicIP_,Wgt); 
  TH2Poly_costheta_Pmu_Purity->Fill(CosTheta2,Pmu2,Wgt);
  TH2Poly_costheta_Pmu_incluive_Purity->Fill(CosTheta2,Pmu2,Wgt);
 
 
  h_p_purity->Fill(Pmu,Wgt);
  h_costheta_p_purity->Fill(CosTheta,Wgt);
  
 }

 h_VertexX_EventCategory.GetComponentHist(TrueEventType)->Fill(reco_nu_vtx_sce_x,Wgt);
 h_VertexY_EventCategory.GetComponentHist(TrueEventType)->Fill(reco_nu_vtx_sce_y,Wgt);
 h_VertexZ_EventCategory.GetComponentHist(TrueEventType)->Fill(reco_nu_vtx_sce_z,Wgt);
  

  // at least 1 proton for this varibels 
  if(proton_multiplicity !=0){
  h_leadingProton_p->Fill(P_proton,Wgt);
  h_costheta_proton->Fill(P_CosTheta,Wgt);
  h_leadingProton_p_EventCategory.GetComponentHist(TrueEventType)->Fill(P_proton,Wgt);
  h_costheta_proton_EventCategory.GetComponentHist(TrueEventType)->Fill(P_CosTheta,Wgt);
  h_openingAngle->Fill(OpenAngle,Wgt);  
  h_openingAngle_true->Fill(mc_OpenAngle,Wgt);
  h_openingAngle_Resolution->Fill(OpenAngle_resolution,Wgt);
  h_openingAngle_Mig->Fill(mc_OpenAngle,OpenAngle,Wgt); 
  
  h_openingAngle_deg->Fill(OpenAngle_deg,Wgt);  
  h_openingAngle_deg_true->Fill(mc_OpenAngle_deg,Wgt); 
  
  h_pn->Fill(pn,Wgt); 
  h_pn_true->Fill(mc_pn,Wgt);
  h_pn_Resolution->Fill(pn-mc_pn,Wgt);
  h_pn_Mig->Fill(mc_pn,pn,Wgt); 
  h_pn_EventCategory.GetComponentHist(TrueEventType)->Fill(pn,Wgt);
  
  h_delta_alphaT->Fill(delta_alphaT,Wgt); 
  h_delta_alphaT_true->Fill(mc_delta_alphaT,Wgt);
  h_delta_alphaT_Resolution->Fill(delta_alphaT-mc_delta_alphaT,Wgt);
  h_delta_alphaT_Mig->Fill(mc_delta_alphaT,delta_alphaT,Wgt); 
  h_delta_alphaT_EventCategory.GetComponentHist(TrueEventType)->Fill(delta_alphaT,Wgt);
  
  h_delta_pTx->Fill(delta_pTx,Wgt);  
  h_delta_pTx_true->Fill(mc_delta_pTx,Wgt);
  h_delta_pTx_Resolution->Fill(delta_pTx-mc_delta_pTx,Wgt);
  h_delta_pTx_Mig->Fill(mc_delta_pTx,delta_pTx,Wgt); 
  h_delta_pTx_EventCategory.GetComponentHist(TrueEventType)->Fill(delta_pTx,Wgt);
  
  h_delta_pTy->Fill(delta_pTy,Wgt);  
  h_delta_pTy_true->Fill(mc_delta_pTy,Wgt);
  h_delta_pTy_Resolution->Fill(delta_pTy-mc_delta_pTy,Wgt);
  h_delta_pTy_Mig->Fill(mc_delta_pTy,delta_pTy,Wgt); 
  h_delta_pTy_EventCategory.GetComponentHist(TrueEventType)->Fill(delta_pTy,Wgt);
  

  h_delta_phiT->Fill(delta_phiT,Wgt);  
  h_delta_phiT_true->Fill(mc_delta_phiT,Wgt);
  h_delta_phiT_Resolution->Fill(delta_phiT-mc_delta_phiT,Wgt);
  h_delta_phiT_Mig->Fill(mc_delta_phiT,delta_phiT,Wgt); 
  h_delta_phiT_EventCategory.GetComponentHist(TrueEventType)->Fill(delta_phiT,Wgt);
  h_openingAngle_Deg_EventCategory.GetComponentHist(TrueEventType)->Fill(OpenAngle_deg,Wgt);
 }

 /*
 std::map<double, TKI_parameters> TKI_proton_map; 

 for(int index1 = 0 ; index1 <pfpdg->size(); index1++ ){

 if(index1 == index_muon)continue;

 TKI_parameters input_TKI = TKI_stuct(
 P_p_hadronTk->at(index1),
 P_p_hadronTk_mc->at(index1),
 P_costheta_hadronTk->at(index1),
 P_costheta_hadronTk_mc->at(index1),
 delta_pT_hadronTk->at(index1),
 delta_pT_hadronTk_mc->at(index1),
 delta_alphaT_hadronTk->at(index1),
 delta_alphaT_hadronTk_mc->at(index1),
 delta_pTx_hadronTk->at(index1),
 delta_pTx_hadronTk_mc->at(index1),
 delta_pTy_hadronTk->at(index1),
 delta_pTy_hadronTk_mc->at(index1),
 delta_phiT_hadronTk->at(index1),
 delta_phiT_hadronTk_mc->at(index1),
 pn_hadronTk->at(index1),
 pn_hadronTk_mc->at(index1),
 TrueEventType);

 TKI_proton_map.insert(std::pair<double,TKI_parameters>(P_p_hadronTk->at(index1),input_TKI));

 }

 int track_type = 1; 

 for(auto TKI_track:TKI_proton_map){
  if(track_type==1){
     TKI_Track1.FillRECOHist(TKI_track.second, Wgt);
  }
  else if(track_type==2){
     TKI_Track2.FillRECOHist(TKI_track.second, Wgt);
  }
  else if(track_type==3){
     TKI_Track3.FillRECOHist(TKI_track.second, Wgt);
  }
  else {
     TKI_Track4.FillRECOHist(TKI_track.second, Wgt);
  }

track_type++;
}
*/






}// ENd of If Statement, that passed Signal Selection for RECO EVENTS 
//////////////////////////////////////////////////////  


if (mc_signal== true){

//std::cout<<"mc_signal== true"<< std::endl;


 double Pmu_true = mc_p3_mu->Mag();
 Double_t Pmu_true2 = mc_p3_mu->Mag();
 
 Double_t Pmu_RECO = p3_mu->Mag();
 Double_t CosTheta_RECO = p3_mu->CosTheta();
 
 
 double CosTheta_true = mc_p3_mu->CosTheta();
 Double_t CosTheta_true2 = mc_p3_mu->CosTheta();
 double Wgtt = spline_weight * tuned_cv_weight;
 Double_t Wgt = safe_weight_histMaker( Wgtt );
 float CosmicIP_ = CosmicIP;
 //std:: cout<<"CosmicIP_ = "<< CosmicIP_<< std::endl;
 int TrueProtonN = Num_MatchedPDGs(*backtracked_pdg,2212);
 int RecoProtonN = Num_MatchedPDGs(*trk_BDT_PID, BDT_Category::kBDT_Proton);
 
 

 
 int TrackN_proton_true = 1;
 EventCategory TrueEventType = EventCategory_tool.returnEventCategoryType(category);
 /*
 std::map<double , TKI_parameters> TKI_proton_map; 

 for(int index1 = 0 ; index1 <pfpdg->size(); index1++ ){

 if(index1 == index_muon)continue;

  TKI_parameters input_TKI = TKI_stuct(
  P_p_hadronTk->at(index1),
  P_p_hadronTk_mc->at(index1),
  P_costheta_hadronTk->at(index1),
  P_costheta_hadronTk_mc->at(index1),
  delta_pT_hadronTk->at(index1),
  delta_pT_hadronTk_mc->at(index1),
  delta_alphaT_hadronTk->at(index1),
  delta_alphaT_hadronTk_mc->at(index1),
  delta_pTx_hadronTk->at(index1),
  delta_pTx_hadronTk_mc->at(index1),
  delta_pTy_hadronTk->at(index1),
  delta_pTy_hadronTk_mc->at(index1),
  delta_phiT_hadronTk->at(index1),
  delta_phiT_hadronTk_mc->at(index1),
  pn_hadronTk->at(index1),
  pn_hadronTk_mc->at(index1),
  TrueEventType);

 TKI_proton_map.insert(std::pair<double, TKI_parameters>(P_p_hadronTk->at(index1),input_TKI));

 }

 int track_type = 1; 

 for(auto TKI_track:TKI_proton_map){
  if(track_type==1){
     TKI_Track1.FillTRUEHist(TKI_track.second, Wgt);
      if (  sel_CC0pi && sel_BDT_NotBOGUS ){
      TKI_Track1.FillTRUE_RECOHist(TKI_track.second, Wgt);
      }
     
     
  }
  else if(track_type==2){
     TKI_Track2.FillTRUEHist(TKI_track.second, Wgt);
      if (  sel_CC0pi && sel_BDT_NotBOGUS ){
      TKI_Track2.FillTRUE_RECOHist(TKI_track.second, Wgt);
      }
  }
  else if(track_type==3){
     TKI_Track3.FillTRUEHist(TKI_track.second, Wgt);
      if (  sel_CC0pi && sel_BDT_NotBOGUS ){
      TKI_Track3.FillTRUE_RECOHist(TKI_track.second, Wgt);
      }
     
  }
  else {
     TKI_Track4.FillTRUEHist(TKI_track.second, Wgt);
    if (  sel_CC0pi && sel_BDT_NotBOGUS ){
      TKI_Track4.FillTRUE_RECOHist(TKI_track.second, Wgt);
      }
  }

 track_type++;
 }

 */

 for(int index1 = 0; index1 < pfpdg->size(); index1++ )
 {
  int true_pdg = backtracked_pdg->at(index1);
  double score_1 = trk_BDT_score_1->at(index1);
  double score_3 = trk_BDT_score_3->at(index1);
  int BDT_Prediction =  trk_BDT_PID->at(index1);

 

 
 
  //if( std::isinf(score_1) || std::isnan(score_1)  || std::isinf(Wgt) || std::isnan(Wgt) )continue; 
   //if( )continue; 
 
   if(index1 == 0  ){
  // Hists that are filled only once in this loop ,
  //I do it this way so all measures would 
  h_Proton_mult_TRUE->Fill(TrueProtonN,Wgt); 
  
  }
 
 
  if(index1 == index_muon){
  //std::cout<<"inside true prediction score_1,Wgt = "<< score_1<< ", "<< Wgt<< std::endl;
    h_CCZeroP_MuonCandidate_BTDPrediction_TRUE->Fill(score_1,Wgt);
    h_costheta_p_TRUE->Fill(CosTheta_true,Wgt); 
    h_costheta_Pmu_TRUE->Fill(CosTheta_true,Pmu_true,Wgt); 
    TH2Poly_costheta_Pmu_TRUE->Fill(CosTheta_true2,Pmu_true2,Wgt); 
    TH2Poly_costheta_Pmu_TRUE_test->Fill(CosTheta_true2,Pmu_true2,Wgt); 
    TH2Poly_costheta_Pmu_incluive_TRUE->Fill(CosTheta_true2,Pmu_true2,Wgt); 
    h_p_TRUE->Fill(Pmu_true,Wgt);
    h_topological_score_TRUE->Fill(topological_score,Wgt);
    h_trk_len_v_TRUE->Fill(trk_len_v->at(index1),Wgt);
  }

  else if(index1 != index_muon){
    h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE->Fill(score_3,Wgt);
     
     if(BDT_Prediction==BDT_Category::kBDT_Proton && true_pdg == 2212){
           h_Proton_mult_Tracks_TRUE->Fill(TrackN_proton_true,Wgt);
           
          if ( sel_CC0pi == true  ){
          h_Proton_mult_Tracks_TRUE_RECO->Fill(TrackN_proton_true,Wgt);
          }
      TrackN_proton_true++;
      }
  }
  
  h_trk_distance_v_TRUE->Fill(trk_distance_v->at(index1),Wgt);
  h_trk_score_v_TRUE->Fill(trk_score_v->at(index1),Wgt);
  h_Cosmic_Cut_TRUE->Fill(CosmicIP_,Wgt);

  /*&& BTD_Muon_Candidate_prediction != 9999 
                && !std::isnan(BTD_Muon_Candidate_prediction) 
                && !std::isinf(BTD_Muon_Candidate_prediction) */
 
  if ( sel_CC0pi    ){
  //&& sel_BDT_predicts_1plusTrks_tobeProtons 
  //&& sel_BDT_NotBOGUS
  // 
  
//sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV &&
//sel_has_muon_candidate && sel_topo_cut_passed && 
//sel_no_reco_showers && sel_muon_passed_mom_cuts && 
//sel_pions_above_threshold && sel_has_pi_candidate && pion_multiplicity > 0
  
   h_trk_score_v_TRUE_RECO->Fill(trk_score_v->at(index1),Wgt); 
   h_trk_distance_v_TRUE_RECO->Fill(trk_distance_v->at(index1),Wgt);
   h_Cosmic_Cut_TRUE_RECO->Fill(CosmicIP_,Wgt);
   
   
   
   
  if(index1 == 0){
      h_Proton_mult_TRUE_RECO->Fill(TrueProtonN,Wgt); 
      h_Proton_mult_Mig->Fill(TrueProtonN,RecoProtonN,Wgt);
      
      
    }

  if(index1 == index_muon){
      h_CCZeroP_MuonCandidate_BTDPrediction_TRUE_RECO->Fill(score_1,Wgt);
      h_costheta_Pmu_TRUE_RECO->Fill(CosTheta_true,Pmu_true,Wgt); 
      TH2Poly_costheta_Pmu_TRUE_RECO->Fill(CosTheta_true2,Pmu_true2,Wgt); 
      TH2Poly_costheta_Pmu_incluive_TRUE_RECO->Fill(CosTheta_true2,Pmu_true2,Wgt); 
      h_costheta_p_TRUE_RECO->Fill(CosTheta_true,Wgt); 
      h_topological_score_TRUE_RECO->Fill(topological_score,Wgt);
      h_p_TRUE_RECO->Fill(Pmu_true,Wgt); 
      h_trk_len_v_TRUE_RECO->Fill(trk_len_v->at(index1),Wgt);
      
      Int_t trueBin1 = TH2Poly_costheta_Pmu_TRUE2->Fill(CosTheta_true,Pmu_true,Wgt); 
      Int_t RECOBin1 = TH2Poly_costheta_Pmu_RECO2->Fill(CosTheta_RECO, Pmu_RECO,Wgt); 
      
      h_matr->Fill(trueBin1,RECOBin1,Wgt);
      
      Int_t trueBininclusive = TH2Poly_costheta_Pmu_incluive_TRUE2->Fill(CosTheta_true, Pmu_true,Wgt); 
      Int_t RECOBininclusive = TH2Poly_costheta_Pmu_incluive_RECO2->Fill(CosTheta_RECO, Pmu_RECO,Wgt); 
      
     h_matr_inclusive->Fill(trueBininclusive,RECOBininclusive,Wgt);


    }
    
   if(index1 != index_muon){
      h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE_RECO->Fill(score_3,Wgt);
    }
    
  }// End Of RECO signal 
 ///////////////////////
}//// End of  Particle Loop ls

///////////////////////
}// End of MC TRUE




}// End of Event Loop
///////////////////////////////////////////////////////
f->Close();
delete f; 

}


/////////////////////////
// Make Root File to look 
/////////////////////////



 //OutPutRootName = 
//addthresholds_WithProtonPIDCUT_
std::cout<<" Writing Root File :"<< OutPutRootName << std::endl;

auto outFile = TFile::Open(OutPutRootName.c_str(), "RECREATE");
outFile->cd(); 

h_p->Write();
h_p_Correction->Write();
h_p_true->Write();
h_p_Resolution->Write();
h_p_MCS_Resolution->Write();
h_p_Mig->Write();
h_p_MCS_Mig->Write();
h_p_purity->Write();
h_costheta_p_purity->Write();
h_p_TRUE->Write();
h_p_TRUE_RECO->Write();
h_costheta_p->Write();
h_costheta_p_true->Write();
h_costheta_p_Resolution->Write();
h_costheta_p_MCS_Resolution->Write();
h_costheta_p_Mig->Write();
h_costheta_p_MCS_Mig->Write();
h_costheta_p_TRUE->Write();
h_costheta_p_TRUE_RECO->Write();
h_p_costheta_p->Write();



h_pn->Write();
h_pn_true->Write();
h_pn_Resolution->Write();
h_pn_Mig->Write();

h_Proton_mult->Write();
h_Proton_mult_TRUE->Write();
h_Proton_mult_TRUE_RECO->Write();
h_Proton_mult_Mig->Write();
h_Proton_mult_Purity->Write(); 

h_Pion_mult->Write();

h_p_EventCategory.WriteToFile(*outFile);
h_p_containment_categories.WriteToFile(*outFile);
h_p_topology_categories.WriteToFile(*outFile);

h_p_Correction_EventCategory.WriteToFile(*outFile);
h_p_Correction_containment_categories.WriteToFile(*outFile);
h_p_Correction_topology_categories.WriteToFile(*outFile);

h_costheta_p_topology_categories.WriteToFile(*outFile);
h_Proton_mult_EventCategory.WriteToFile(*outFile);
h_NTracks_EventCategory.WriteToFile(*outFile);

h_costheta_p_EventCategory.WriteToFile(*outFile);
h_costheta_p_containment_categories.WriteToFile(*outFile);

h_NTracks_BTDGroup_categories.WriteToFile(*outFile);
h_NTracks_Particle_categories.WriteToFile(*outFile);

h_topological_score->Write();
h_topological_score_Particle_categories.WriteToFile(*outFile);
h_topological_score_EventCategory.WriteToFile(*outFile);
h_topological_score_TRUE->Write();
h_topological_score_TRUE_RECO->Write();


h_trk_score_v->Write();
h_trk_score_v_Particle_categories.WriteToFile(*outFile);
h_trk_score_v_EventCategory.WriteToFile(*outFile);
h_trk_score_v_TRUE->Write();
h_trk_score_v_TRUE_RECO->Write();




h_trk_len_v->Write();
h_trk_len_v_Purity->Write();
h_trk_len_v_Particle_categories.WriteToFile(*outFile);
h_trk_len_v_EventCategory.WriteToFile(*outFile);

h_trk_len_v_TRUE->Write();
h_trk_len_v_TRUE_RECO->Write();

h_trk_distance_v->Write();
h_trk_distance_v_Particle_categories.WriteToFile(*outFile);
h_trk_distance_v_EventCategory.WriteToFile(*outFile);
h_trk_distance_v_TRUE->Write();
h_trk_distance_v_TRUE_RECO->Write();


h_Cosmic_Cut_TRUE->Write();
h_Cosmic_Cut_TRUE_RECO->Write();
h_Cosmic_Cut_Purity->Write();
h_Cosmic_Cut->Write();

h_costheta_Pmu_TRUE->Write();
h_costheta_Pmu_TRUE_RECO->Write();
h_costheta_Pmu_Purity->Write();
h_costheta_Pmu_EventCategory.WriteToFile(*outFile);
h_costheta_Pmu->Write();

h_Else_mult_Tracks_Particle_categories.WriteToFile(*outFile);
h_Else_mult_Tracks_EventCategory.WriteToFile(*outFile);

h_NTracks->Write();
h_Proton_mult_Tracks_Particle_categories.WriteToFile(*outFile);
h_Proton_mult_Tracks_EventCategory.WriteToFile(*outFile);
h_Proton_mult_Tracks_TRUE->Write();
h_Proton_mult_Tracks_TRUE_RECO->Write();
h_Proton_mult_Tracks_Purity->Write();
h_Proton_mult_Tracks->Write();
h_pn_EventCategory.WriteToFile(*outFile);
h_delta_alphaT_EventCategory.WriteToFile(*outFile);
h_delta_pTx_EventCategory.WriteToFile(*outFile);
h_delta_pTy_EventCategory.WriteToFile(*outFile);
h_delta_phiT_EventCategory.WriteToFile(*outFile);
h_openingAngle_Deg_EventCategory.WriteToFile(*outFile);
h_VertexX_EventCategory.WriteToFile(*outFile);
h_VertexY_EventCategory.WriteToFile(*outFile);
h_VertexZ_EventCategory.WriteToFile(*outFile);

h_VertexX_leadmuon_Particle_categories.WriteToFile(*outFile);
h_VertexY_leadmuon_Particle_categories.WriteToFile(*outFile);
h_VertexZ_leadmuon_Particle_categories.WriteToFile(*outFile);


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



h_CCZeroP_BTDGroup_Probability_tobeMuon_tk1->Write();
h_CCZeroP_BTDGroup_Probability_tobeProton_trk1->Write();
h_CCZeroP_BTDGroup_Probability_tobePion_trk1->Write();
h_CCZeroP_BTDGroup_Probability_tobeElse_trk1->Write();

h_CC1P_BTDGroup_Probability_tobeMuon_tk1->Write();
h_CC1P_BTDGroup_Probability_tobeProton_tk1->Write();
h_CC1P_BTDGroup_Probability_tobePion_tk1->Write();
h_CC1P_BTDGroup_Probability_tobeElse_tk1->Write();
h_CC1P_BTDGroup_Probability_tobeMuon_tk2->Write();
h_CC1P_BTDGroup_Probability_tobeProton_tk2->Write();
h_CC1P_BTDGroup_Probability_tobePion_tk2->Write();
h_CC1P_BTDGroup_Probability_tobeElse_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeMuon_tk1->Write();
h_CC2P_BTDGroup_Probability_tobeProton_tk1->Write();
h_CC2P_BTDGroup_Probability_tobePion_tk1->Write();
h_CC2P_BTDGroup_Probability_tobeElse_tk1->Write();
h_CC2P_BTDGroup_Probability_tobeMuon_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeProton_tk2->Write();
h_CC2P_BTDGroup_Probability_tobePion_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeElse_tk2->Write();
h_CC2P_BTDGroup_Probability_tobeMuon_tk3->Write();
h_CC2P_BTDGroup_Probability_tobeProton_tk3->Write();
h_CC2P_BTDGroup_Probability_tobePion_tk3->Write();
h_CC2P_BTDGroup_Probability_tobeElse_tk3->Write();


h_CCZeroP_MuonCandidate_BTDPrediction->Write();
h_CCZeroP_MuonCandidate_BTDPrediction_Particle_categories.WriteToFile(*outFile);
h_CCZeroP_MuonCandidate_BTDPrediction_EventCategory.WriteToFile(*outFile);
h_CCZeroP_MuonCandidate_BTDPrediction_purity->Write();
h_CCZeroP_MuonCandidate_BTDPrediction_topology_categories.WriteToFile(*outFile);
h_CCZeroP_MuonCandidate_BTDPrediction_containment_categories.WriteToFile(*outFile);


h_CCZeroP_ProtonCandidate_BTDPrediction->Write();
h_CCZeroP_ProtonCandidate_BTDPrediction_Particle_categories.WriteToFile(*outFile);
h_CCZeroP_ProtonCandidate_BTDPrediction_EventCategory.WriteToFile(*outFile);
h_CCZeroP_ProtonCandidate_BTDPrediction_topology_categories.WriteToFile(*outFile);
h_CCZeroP_ProtonCandidate_BTDPrediction_containment_categories.WriteToFile(*outFile);


h_CCZeroP_Else_BTDPrediction_Particle_categories.WriteToFile(*outFile);
h_CCZeroP_Else_BTDPrediction_EventCategory.WriteToFile(*outFile);
h_CCZeroP_Else_BTDPrediction->Write();
h_CCZeroP_ALLProtonCandidate_BTDPrediction->Write();
h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE->Write();
h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE_RECO->Write();
h_CCZeroP_ALLProtonCandidate_BTDPrediction_Purity->Write();

h_leadingProton_p->Write();
h_costheta_proton->Write();
h_leadingProton_p_EventCategory.WriteToFile(*outFile);
h_costheta_proton_EventCategory.WriteToFile(*outFile);


h_p_CCNP0Pi_Only->Write();
h_p_CCNP0Pi_Only_true->Write();
h_p_CCNP0Pi_Only_Resolution->Write();
h_p_CCNP0Pi_Only_Mig->Write();
h_p_CCNP0Pi_Only_EventCategory.WriteToFile(*outFile);
h_p_CCNP0Pi_Only_containment_categories.WriteToFile(*outFile);
h_costheta_p_CCNP0Pi_Only->Write();
h_costheta_p_CCNP0Pi_Only_true->Write();
h_costheta_p_CCNP0Pi_Only_Resolution->Write();
h_costheta_p_CCNP0Pi_Only_Mig->Write();
h_costheta_p_CCNP0Pi_Only_EventCategory.WriteToFile(*outFile);
h_costheta_p_CCNP0Pi_Only_containment_categories.WriteToFile(*outFile);
h_p_CC0P0Pi_Only->Write();
h_p_CC0P0Pi_Only_true->Write();
h_p_CC0P0Pi_Only_Resolution->Write();
h_p_CC0P0Pi_Only_Mig->Write();
h_p_CC0P0Pi_Only_EventCategory.WriteToFile(*outFile);
h_p_CC0P0Pi_Only_containment_categories.WriteToFile(*outFile);
h_costheta_p_CC0P0Pi_Only->Write();
h_costheta_p_CC0P0Pi_Only_true->Write();
h_costheta_p_CC0P0Pi_Only_Resolution->Write();
h_costheta_p_CC0P0Pi_Only_Mig->Write();
h_costheta_p_CC0P0Pi_Only_EventCategory.WriteToFile(*outFile);
h_costheta_p_CC0P0Pi_Only_containment_categories.WriteToFile(*outFile);
h_pionE_Mig->Write();

h_CCPiCategory->Write();

h_CCZeroP_MuonCandidate_BTDPrediction_TRUE->Write();
h_CCZeroP_MuonCandidate_BTDPrediction_TRUE_RECO->Write();

h_matr->Write();
h_matr_inclusive->Write();
 
//TKI_Track1.WriteAll(*outFile);
//TKI_Track2.WriteAll(*outFile);
//TKI_Track3.WriteAll(*outFile);
//TKI_Track4.WriteAll(*outFile);


//////////////////////////////////////////////////////////////////////


int precision = 2;  // Change this value based on your requirement
gStyle->SetPaintTextFormat(Form(".%df", precision));
float textSize = 0.015;  // Adjust this value based on your requirement
gStyle->SetTextSize(textSize);

gStyle->SetPalette(kTemperatureMap);
float Marksize_eff= .015;


   TCanvas *can = new TCanvas("can");
   gStyle->SetOptStat(0); 
    char text_title_pdf1[2024];
    char text_title_pdf2[2024];
    char text_title_pdf3[2024];
    char text_title_pdf4[2024];

  sprintf(text_title_pdf1, "Make_Plots_2D_Effciency.pdf(","" );
  can -> Print(text_title_pdf1);
  sprintf(text_title_pdf2, "Make_Plots_2D_Effciency.pdf","" );
  sprintf(text_title_pdf3, "Make_Plots_2D_Effciency.pdf)","" );
  sprintf(text_title_pdf4, "Make_Plots_2D_Effciency","" );


  TH2Poly_costheta_Pmu_TRUE_test->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_TRUE_test->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_TRUE_test->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_TRUE_test->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_TRUE_test->GetZaxis()->SetLabelSize (0.017);
  TH2Poly_costheta_Pmu_TRUE_test->SetTitle("Generated MC Selected"); 
  TH2Poly_costheta_Pmu_TRUE_test->Draw("COLZ");
  can -> Print(text_title_pdf2);

  TH2Poly_costheta_Pmu_TRUE_RECO->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_TRUE_RECO->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_TRUE_RECO->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_TRUE_RECO->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_TRUE_RECO->SetTitle("RECO MC Selected"); 
  TH2Poly_costheta_Pmu_TRUE_RECO->GetZaxis()->SetLabelSize (0.017);
  TH2Poly_costheta_Pmu_TRUE_RECO->Draw("COLZ");
   saveUBTH2PolyToTextFile(*TH2Poly_costheta_Pmu_TRUE_RECO, "Num_RECO_Entries.txt");
  can -> Print(text_title_pdf2);
  gStyle->SetPalette(kCool);
  TH2Poly_costheta_Pmu_TRUE_RECO->Divide(TH2Poly_costheta_Pmu_TRUE);
  TH2Poly_costheta_Pmu_TRUE_RECO->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_TRUE_RECO->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_TRUE_RECO->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_TRUE_RECO->GetYaxis()->CenterTitle();
  saveUBTH2PolyToTextFile(*TH2Poly_costheta_Pmu_TRUE_RECO, "Binning_Eff.txt");
  TH2Poly_costheta_Pmu_TRUE_RECO->SetTitle("Efficiency"); 
  TH2Poly_costheta_Pmu_TRUE_RECO->Draw("COLZ");
  can -> Print(text_title_pdf2);
  TH2Poly_costheta_Pmu_TRUE_RECO->SetMarkerSize(.7);
  TH2Poly_costheta_Pmu_TRUE_RECO->Draw("text");
    can -> Print(text_title_pdf2);
    
    
    
/////////////////////////////
 
PROJECTION_Bin_Map Projection_binMap = GetBin_ProjectionMap();
char HistName[1024];

for(auto BinMap :Projection_binMap ){

  sprintf(HistName, "EfficiencyProjections_Bin_%i", BinMap.first);
  TH1D* Binning =  TH2Poly_costheta_Pmu_TRUE_RECO->Projection(HistName, BinMap.second);
  Binning->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  Binning->Write(HistName);

    sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_Bin_%i", BinMap.first);
  TH1D* Binning2 =  h_costheta_Pmu_UBTH2Poly->Projection(HistName, BinMap.second); 
    Binning2->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  Binning2->Write(HistName);
  
   sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_Contained_Bin_%i", BinMap.first);
  TH1D* Binning3 =  h_costheta_Pmu_UBTH2Poly_Contained->Projection(HistName, BinMap.second); 
   Binning3->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  Binning3->Write(HistName);
  
     sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_NonContained_Bin_%i", BinMap.first);
  TH1D* Binning4 =  h_costheta_Pmu_UBTH2Poly_NonContained->Projection(HistName, BinMap.second); 
   Binning4->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  Binning4->Write(HistName);
  
  
  sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_Signal_Bin_%i", BinMap.first);
  TH1D* Binning5 =  h_costheta_Pmu_UBTH2Poly_Signal->Projection(HistName, BinMap.second); 
   Binning5->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  Binning5->Write(HistName);
  
     for(auto x: SignalCategory_vector){
      
     TH1D* Binning_EventType =    h_costheta_Pmu_UBTH2Poly_EventCategory.GetComponentHist(x)->Projection(uniq(), BinMap.second);
     sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_EventCategory_%s_Bin_%i",Binning_EventType->GetTitle(), BinMap.first);
     std::cout<<"Writeing 2D Bin Named: "<< HistName<< std::endl;
     Binning_EventType->Write(HistName);
     
     
      TH1D* Binning_EventType2 =    h_costheta_Pmu_UBTH2Poly_Contained_EventCategory.GetComponentHist(x)->Projection(uniq(), BinMap.second);
     sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_Contained_EventCategory_%s_Bin_%i",Binning_EventType2->GetTitle(), BinMap.first);
     std::cout<<"Writeing 2D Bin Named: "<< HistName<< std::endl;
     Binning_EventType2->Write(HistName);
     
      TH1D* Binning_EventType3 =    h_costheta_Pmu_UBTH2Poly_NonContained_EventCategory.GetComponentHist(x)->Projection(uniq(), BinMap.second);
     sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_NonContained_EventCategory_%s_Bin_%i",Binning_EventType3->GetTitle(), BinMap.first);
     std::cout<<"Writeing 2D Bin Named: "<< HistName<< std::endl;
     Binning_EventType3->Write(HistName);
     
          
      TH1D* Binning_EventType4 =    h_costheta_Pmu_UBTH2Poly_Signal_EventCategory.GetComponentHist(x)->Projection(uniq(), BinMap.second);
     sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_Signal_EventCategory_%s_Bin_%i",Binning_EventType4->GetTitle(), BinMap.first);
     std::cout<<"Writeing 2D Bin Named: "<< HistName<< std::endl;
     Binning_EventType4->Write(HistName);
     
  
  }
  



}


      
     can -> Print(text_title_pdf2);
/////////////////////
gStyle->SetPalette(kTemperatureMap);
  TH2Poly_costheta_Pmu_incluive_TRUE->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_incluive_TRUE->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_TRUE->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_incluive_TRUE->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_TRUE->SetTitle("Generated MC Selected"); 
  TH2Poly_costheta_Pmu_incluive_TRUE->GetZaxis()->SetLabelSize (0.017);
  TH2Poly_costheta_Pmu_incluive_TRUE->Draw("COLZ");
  can -> Print(text_title_pdf2);

  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetZaxis()->SetLabelSize (0.017);
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->SetTitle("RECO MC Selected"); 
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->Draw("COLZ");
  
  can -> Print(text_title_pdf2);
  gStyle->SetPalette(kCool);
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->Divide(TH2Poly_costheta_Pmu_incluive_TRUE);
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->SetTitle("Efficiency"); 
  saveUBTH2PolyToTextFile(*TH2Poly_costheta_Pmu_incluive_TRUE_RECO, "Binning_Eff_inclusive.txt");
  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->Draw("COLZ");
  can -> Print(text_title_pdf2);
TH2Poly_costheta_Pmu_incluive_TRUE_RECO->SetMarkerSize(.7);
TH2Poly_costheta_Pmu_incluive_TRUE_RECO->Draw("text");
can -> Print(text_title_pdf2);
////////////////////////////
/////////////////////
gStyle->SetPalette(kTemperatureMap);
  TH2Poly_costheta_Pmu->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu->SetTitle("RECO Selected"); 
  TH2Poly_costheta_Pmu->Draw("COLZ");
  can -> Print(text_title_pdf2);

  TH2Poly_costheta_Pmu_Purity->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_Purity->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_Purity->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_Purity->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_Purity->SetTitle("RECO + MC Selected"); 
  TH2Poly_costheta_Pmu_Purity->Draw("COLZ");
  can -> Print(text_title_pdf2);
    gStyle->SetPalette(kCool);
  TH2Poly_costheta_Pmu_Purity->Divide(TH2Poly_costheta_Pmu);
  TH2Poly_costheta_Pmu_Purity->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_Purity->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_Purity->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_Purity->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_Purity->GetZaxis()->SetLabelSize (0.017);
  TH2Poly_costheta_Pmu_Purity->SetTitle("Purity"); 
  TH2Poly_costheta_Pmu_Purity->Draw("COLZ");
  can -> Print(text_title_pdf2);
  TH2Poly_costheta_Pmu_Purity->SetMarkerSize(.7);
  TH2Poly_costheta_Pmu_Purity->Draw("text");
  can -> Print(text_title_pdf2);
  
  for(auto BinMap :Projection_binMap ){

  sprintf(HistName, "PurityProjections_Bin_%i", BinMap.first);
  TH1D* Binning =  TH2Poly_costheta_Pmu_Purity->Projection(HistName, BinMap.second);
  Binning->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  Binning->Write(HistName);
  
  }
  
  
////////////////////////////

/////////////////////
gStyle->SetPalette(kViridis);  //kTemperatureMap
  TH2Poly_costheta_Pmu_incluive->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_incluive->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_incluive->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive->SetTitle("RECO Selected"); 
  TH2Poly_costheta_Pmu_incluive->GetZaxis()->SetLabelSize (0.017);
  TH2Poly_costheta_Pmu_incluive->Draw("COLZ");
  can -> Print(text_title_pdf2);

  TH2Poly_costheta_Pmu_incluive_Purity->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_incluive_Purity->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_Purity->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_incluive_Purity->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_Purity->SetTitle("RECO + MC Selected");
  TH2Poly_costheta_Pmu_incluive_Purity->GetZaxis()->SetLabelSize (0.017);
  TH2Poly_costheta_Pmu_incluive_Purity->Draw("COLZ");
  can -> Print(text_title_pdf2);
  TH2Poly_costheta_Pmu_incluive_Purity->Divide(TH2Poly_costheta_Pmu_incluive);
  gStyle->SetPalette(kViridis);
  TH2Poly_costheta_Pmu_incluive_Purity->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
  TH2Poly_costheta_Pmu_incluive_Purity->GetXaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_Purity->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  TH2Poly_costheta_Pmu_incluive_Purity->GetYaxis()->CenterTitle();
  TH2Poly_costheta_Pmu_incluive_Purity->SetTitle("Purity"); 
  TH2Poly_costheta_Pmu_incluive_Purity->Draw("COLZ");
  can -> Print(text_title_pdf2);
  TH2Poly_costheta_Pmu_incluive_Purity->SetMarkerSize(.7);
  TH2Poly_costheta_Pmu_incluive_Purity->Draw("text");
  can -> Print(text_title_pdf2);
////////////////////////////
gStyle->SetPalette(kViridis);
h_matr->GetYaxis()->SetTitle("RECO bins");
h_matr->GetXaxis()->SetTitle("TRUE bins");
h_matr->SetTitle("2D Migration Added little plot"); 
h_matr->GetZaxis()->SetLabelSize (0.017);
h_matr->Draw("COLZ");
can -> Print(text_title_pdf2);

int e; 

PROJECTION_InclusBin_Map Projection_InclusivebinMap = GetInclusvieBin_ProjectionMap(); 



   std::vector<int> bin_numbers;
   std::vector<int> StartingXbins{1,2,3,4,5,6,7,8,9};


for(auto startingbin: StartingXbins){
//for(auto vec:BinMap.second ){std::cout<<vec<<",";}
  //std::cout<<" "<< std::endl;
  sprintf(HistName, "EfficiencyProjections_Inclusive_Bin_%i", startingbin);
  int StartBin = startingbin;
  //std::cout<<"StartBin = "<< StartBin<< std::endl;
  TH1D* Bin =  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->ProjectionY(HistName, StartBin, bin_numbers);
  Bin->Write(HistName);
  
  TH1D* BinPurity =TH2Poly_costheta_Pmu_incluive_Purity->ProjectionY(HistName, StartBin, bin_numbers);
  sprintf(HistName, "PurityProjections_Inclusive_Bin_%i", startingbin);
  BinPurity->Write(HistName);
  
  sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_Bin_%i", startingbin);
  //std::cout<<"StartBin = "<< StartBin<< std::endl;
  TH1D* Bin2 =  h_costheta_Pmu_UBTH2Poly_inclusive->ProjectionY(HistName, StartBin, bin_numbers);
  Bin2->Write(HistName);
  
    sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_Contained_Bin_%i", startingbin);
  //std::cout<<"StartBin = "<< StartBin<< std::endl;
  TH1D* Bin3 =  h_costheta_Pmu_UBTH2Poly_inclusive_Contained->ProjectionY(HistName, StartBin, bin_numbers);
  Bin3->Write(HistName);
  
  sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_NonContained_Bin_%i", startingbin);
  //std::cout<<"StartBin = "<< StartBin<< std::endl;
  TH1D* Bin4 =  h_costheta_Pmu_UBTH2Poly_inclusive_NonContained->ProjectionY(HistName, StartBin, bin_numbers);
  Bin4->Write(HistName);
  
    sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_Signal_Bin_%i", startingbin);
  //std::cout<<"StartBin = "<< StartBin<< std::endl;
  TH1D* Bin5 =  h_costheta_Pmu_UBTH2Poly_inclusive_Signal->ProjectionY(HistName, StartBin, bin_numbers);
  Bin5->Write(HistName);
  
     for(auto x: SignalCategory_vector){
        h_costheta_Pmu_UBTH2Poly_inclusive_EventCategory.GetComponentHist(x)->GetTitle();
         h_costheta_Pmu_UBTH2Poly_inclusive_Contained_EventCategory.GetComponentHist(x)->GetTitle();
         h_costheta_Pmu_UBTH2Poly_inclusive_NonContained_EventCategory.GetComponentHist(x)->GetTitle();
         h_costheta_Pmu_UBTH2Poly_inclusive_Signal_EventCategory.GetComponentHist(x)->GetTitle();
   
        TH1D* Binning_EventType =    h_costheta_Pmu_UBTH2Poly_inclusive_EventCategory.GetComponentHist(x)->ProjectionY(uniq(), StartBin, bin_numbers);
        sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_EventCategory_%s_Bin_%i",Binning_EventType->GetTitle(), StartBin);
        std::cout<<"Writeing 2D Inclusive Bin Named: "<< HistName<< std::endl;
        Binning_EventType->Write(HistName);
        
        
        TH1D* Binning_EventType2 =    h_costheta_Pmu_UBTH2Poly_inclusive_Contained_EventCategory.GetComponentHist(x)->ProjectionY(uniq(), StartBin, bin_numbers);
        sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_Contained_EventCategory_%s_Bin_%i",Binning_EventType2->GetTitle(), StartBin);
        std::cout<<"Writeing 2D Inclusive Bin Named: "<< HistName<< std::endl;
        Binning_EventType2->Write(HistName);
        
        TH1D* Binning_EventType3 =    h_costheta_Pmu_UBTH2Poly_inclusive_NonContained_EventCategory.GetComponentHist(x)->ProjectionY(uniq(), StartBin, bin_numbers);
        sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_NonContained_EventCategory_%s_Bin_%i",Binning_EventType3->GetTitle(), StartBin);
        std::cout<<"Writeing 2D Inclusive Bin Named: "<< HistName<< std::endl;
        Binning_EventType3->Write(HistName);
        
        TH1D* Binning_EventType4 =    h_costheta_Pmu_UBTH2Poly_inclusive_Signal_EventCategory.GetComponentHist(x)->ProjectionY(uniq(), StartBin, bin_numbers);
        sprintf(HistName, "h_costheta_Pmu_UBTH2Poly_inclusive_Signal_EventCategory_%s_Bin_%i",Binning_EventType4->GetTitle(), StartBin);
        std::cout<<"Writeing 2D Inclusive Bin Named: "<< HistName<< std::endl;
        Binning_EventType4->Write(HistName);
        
        
        
  }
  
  
}



 TH2Poly_costheta_Pmu_incluive_TRUE_RECO->SetNBinsX(9);
std::cout<<"Nxbins  = "<< TH2Poly_costheta_Pmu_incluive_TRUE_RECO->GetNBinsX()<< std::endl;

 TH1D*TH2Poly_costheta_Pmu_incluive_TRUE_RECO_BIN1 =  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->ProjectionY("TH2Poly_costheta_Pmu_TRUE_RECO_BIN1", 1, bin_numbers);
 //TH1D* TH2Poly_costheta_Pmu_incluive_TRUE_RECO_test2_BIN1 =  TH2Poly_costheta_Pmu_incluive_TRUE_RECO->Projection("TH2Poly_costheta_Pmu_incluive_TRUE_RECO_test2", Projection_InclusivebinMap[kProj_InclBin1]);

 TH2Poly_costheta_Pmu_incluive_TRUE_RECO_BIN1->GetXaxis()->SetTitle("P_{#mu} [GeV/c]");
 TH2Poly_costheta_Pmu_incluive_TRUE_RECO_BIN1->SetTitle("Bin 1 test");
 TH2Poly_costheta_Pmu_incluive_TRUE_RECO_BIN1->Draw("hist");
 can -> Print(text_title_pdf2);
//

 
 std::cout << " Bin number for projection"<< std::endl;
 for(auto b: bin_numbers){
 std::cout<< "bin N = " <<  b<< std::endl;
 }
 std::cout<<"~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::vector<int>binNumbers{1,2,3,4,5,6,7};
 TH1D* TH2Poly_costheta_Pmu_TRUE_RECO_BIN1 =  TH2Poly_costheta_Pmu_TRUE_RECO->Projection("TTH2Poly_costheta_Pmu_TRUE_RECO_BIN1", binNumbers);
 

       TH2Poly_costheta_Pmu_TRUE_RECO_BIN1->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
     TH2Poly_costheta_Pmu_TRUE_RECO_BIN1->SetTitle("Bin 2 test");
     TH2Poly_costheta_Pmu_TRUE_RECO_BIN1->Draw("hist");
 can -> Print(text_title_pdf2);
 TH2Poly_costheta_Pmu_TRUE_RECO_BIN1->Scale(1.0,"width");
      TH2Poly_costheta_Pmu_TRUE_RECO_BIN1->Draw("hist");
 can -> Print(text_title_pdf2);
 

 
 
 
   std::vector<double> xBinEdges;
    for (int binNumber : binNumbers) {
        // ROOT's bin numbers start from 1, while C++ vectors start from 0
        double xBinLowEdge =TH2Poly_costheta_Pmu_TRUE_RECO->GetXaxis()->GetBinLowEdge(binNumber);
        double xBinUpEdge = TH2Poly_costheta_Pmu_TRUE_RECO->GetXaxis()->GetBinUpEdge(binNumber);
        xBinEdges.push_back(xBinLowEdge);
        xBinEdges.push_back(xBinUpEdge);
    }
  
 
    //can -> Print(text_title_pdf2);
     // "A" for axis, "P" for points
    // Optionally customize the graph appearance
    graph_PonosFunction->SetMarkerStyle(20);  // Set marker style
    graph_PonosFunction->SetMarkerColor(kBlue);  // Set marker color
    graph_PonosFunction->SetTitle("Panos P_{#mu} Correction"); 
    graph_PonosFunction->GetXaxis()->SetTitle("P_{#mu} [GeV]");// Set title
    graph_PonosFunction->GetYaxis()->SetTitle("#delta P_{#mu} [GeV]");// Set title
    graph_PonosFunction->Draw("AP");
    can -> Print(text_title_pdf2);
    // Update the canvas to display the changes
    

 Draw_heatMap(
  h_matr, 
  "TRUE bins",
  "RECO bins",
  "2D Migration (Row Norm)",
  text_title_pdf2,
  1,
  can,
  false );
  
  double Xwindow1 = .10; double Ywindow1 =0.6; double Xwindow2 = 0.45;  double Ywindow2 =.9;
    bool IS_Inclusive = true; 
  
  
   Draw_heatMap_4DMigration(
   h_matr, 
  "TRUE bins",
  "RECO bins",
  "2D Migration (Row Norm)",
  text_title_pdf2,
  1,
  can,
  false,
  !IS_Inclusive,
    MUON_2D_BIN_EDGES, 
  Xwindow1, Ywindow1,
 Xwindow2, Ywindow2);
  
    Draw_heatMap_4DMigration(
   h_matr, 
  "TRUE bins",
  "RECO bins",
  "2D Migration (Col Norm)",
  text_title_pdf2,
  0,
  can,
  false,
  !IS_Inclusive,
    MUON_2D_BIN_EDGES, 
  Xwindow1, Ywindow1,
 Xwindow2, Ywindow2);
  
  
   Draw_heatMap(
  h_matr, 
  "TRUE bins",
  "RECO bins",
  "2D Migration (Col Norm)",
  text_title_pdf2,
  0,
  can,
  false );

h_matr_inclusive->GetYaxis()->SetTitle("RECO bins");
h_matr_inclusive->GetXaxis()->SetTitle("TRUE bins");
h_matr_inclusive->SetTitle("2D Migration (Inclusive Binning)"); 
h_matr_inclusive->GetZaxis()->SetLabelSize (0.017);
h_matr_inclusive->Draw("COLZ");
  
  
   Draw_heatMap(
  h_matr_inclusive, 
  "TRUE bins",
  "RECO bins",
  "2D Migration (inclusive) (Row Norm)",
  text_title_pdf2,
  1,
  can,
  false );
  
  UBTH2Poly *h_UBTH2Poly_binningN = Make2DHist_UB_inclusive_BinNumMap(MUON_2D_BIN_EDGES_inclusive,"h_UBTH2Poly_binningN");
  
  
   Draw_heatMap(
  h_matr_inclusive, 
  "TRUE bins",
  "RECO bins",
  "2D Migration (inclusive)(Col Norm)",
  text_title_pdf2,
  0,
  can,
  false );
  

  
  Draw_heatMap_4DMigration(
   h_matr_inclusive, 
  "TRUE bins",
  "RECO bins",
  "2D Migration (inclusive) (Row Norm)",
  text_title_pdf2,
  1,
  can,
  false,
  IS_Inclusive,
  MUON_2D_BIN_EDGES_inclusive, 
 Xwindow1, Ywindow1,
 Xwindow2, Ywindow2);
  
    Draw_heatMap_4DMigration(
   h_matr_inclusive, 
  "TRUE bins",
  "RECO bins",
  "2D Migration (inclusive) (Col Norm)",
  text_title_pdf2,
  0,
  can,
  false,
  IS_Inclusive,
  MUON_2D_BIN_EDGES_inclusive, 
 Xwindow1, Ywindow1,
 Xwindow2, Ywindow2);
  
  
       Draw_heatMap(
  h_Proton_mult_Mig, 
  "TRUE Proton Multiplicity",
  "RECO Proton Multiplicity",
  "Migration",
  text_title_pdf2,
  3,
  can,
  false );
  
     Draw_heatMap(
  h_Proton_mult_Mig, 
  "TRUE Proton Multiplicity",
  "RECO Proton Multiplicity",
  "Migration (Col Norm)",
  text_title_pdf2,
  0,
  can,
  false );
  
       Draw_heatMap(
  h_Proton_mult_Mig, 
  "TRUE Proton Multiplicity",
  "RECO Proton Multiplicity",
  "Migration (Row Norm)",
  text_title_pdf2,
  1,
  can,
  false );
  
  
  
  
  

 Draw_heatMap(
  h_costheta_Pmu_TRUE, 
  "TRUE Cos(#theta_{#mu})",
  "TRUE P_{#mu} ",
  "True Cuts Events",
  text_title_pdf2,
  3,
  can,
  false );
  
  
   Draw_heatMap(
  h_costheta_Pmu_TRUE_RECO, 
  "TRUE Cos(#theta_{#mu})",
  "TRUE P_{#mu}",
  "True + RECO cuts ",
  text_title_pdf2,
  3,
  can,
  false );
  
  
  h_costheta_Pmu_TRUE_RECO->Divide(h_costheta_Pmu_TRUE);
  h_costheta_Pmu_TRUE_RECO->SetMaximum(.4);
  
     Draw_heatMap_notext(
  h_costheta_Pmu_TRUE_RECO, 
  "TRUE Cos(#theta_{#mu})",
  "TRUE P_{#mu}",
  "2D Efficiency",
  text_title_pdf2,
  1,
  can,
  false );
  
       Draw_heatMap_notext(
  h_costheta_Pmu_Purity, 
  "RECO Cos(#theta_{#mu})",
  "RECO P_{#mu}",
  "Purity numerator",
  text_title_pdf2,
  3,
  can,
  false );
  
Draw_heatMap_notext(
  h_costheta_Pmu, 
  "RECO Cos(#theta_{#mu})",
  "RECO P_{#mu}",
  "Purity denominator",
  text_title_pdf2,
  3,
  can,
  false );
  
  
h_costheta_Pmu_Purity->Divide(h_costheta_Pmu);

  Draw_heatMap_notext(
  h_costheta_Pmu_Purity, 
  "RECO Cos(#theta_{#mu})",
  "RECO P_{#mu}",
  "Purity",
  text_title_pdf2,
  3,
  can,
  false );
  
  
  
  can -> Print(text_title_pdf3);
  can->Close(); 
  outFile->Close();
  
  
}
////////////////////////////////////////////////////////////
// END 
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

 int myInt = std::atoi(argv[1]);
   hist_maker_MC(myInt);
   
   
   return 0;
}

////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////


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

int Num_MatchedPDGs(std::vector<int> inputPDGs , int MatchtoPDG){

int matched = 0 ; 

for(auto pdg: inputPDGs ){
 if(pdg == MatchtoPDG){
 matched++;
 }
}

return matched; 


}




// Event weights that are below MIN_WEIGHT, above MAX_WEIGHT, infinite, or NaN
// are reset to unity by this function. Other weights are returned unaltered.
 Double_t safe_weight_histMaker( Double_t w ) {
 constexpr double MIN_WEIGHT_ = 0.;
 constexpr double MAX_WEIGHT_ = 30.;
 
  if ( std::isfinite(w) && w >= MIN_WEIGHT_ && w <= MAX_WEIGHT_ ) return w;
  else return 1.0;
}

/*
void saveTH2ToTextFile(const UBTH2Poly& hist, char* fileName) {
    // Open the file for writing
    std::ofstream outputFile(fileName);

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open the file " << fileName << " for writing." << std::endl;
        return;
    }

Int_t Nbins = hist.GetNumberOfBins();


    // Iterate over bins and write bin number and content to the file
    for (Int_t bin = 1; bin <= Nbins; ++bin) {
            Double_t binContent = hist.GetBinContent(bin);
            // Write to the file
            outputFile << bin << "," << binContent << std::endl;
        }
    
    
    
    

    // Close the file
    outputFile.close();
}
*/


void printText(const std::string& text, const std::string& fileName) {
    // Create an ofstream object with the specified file name
    
    std::ofstream outputFile(fileName, std::ios_base::app);

//dd
    // Check if the file is open
    if (outputFile.is_open()) {
        // Print to the file
        outputFile << text << std::endl;

        // Close the file
        outputFile.close();
    } else {
        // Print an error message if unable to open the file
        std::cerr << "Error: Unable to open the file '" << fileName << "'." << std::endl;
    }
}

