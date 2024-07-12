/*Author: Christian
 *
 * Usage: making hist to plot and look at for CC-0pi
 *script to test my plotting method 
 */

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
//#include "HistMaker.hh"
////////////////////////////////////////////////////////
//Categorys For Stacks 
/////////////////////////////////////////////////////////

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


//////Functions//////////////////////////////
void reorderVector(std::vector<int>& vec, int x);
std::vector<int> MakereorderVector(std::vector<int> vec, int GotoFront);



void hist_maker_DATA(int input_Int) {

std::string base_dir = "/uboone/data/users/gardiner/ntuples-stv/";


 //rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_D_E_all_reco2.root
 //rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_G1_beam_good_reco2_1e19.root
 //rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root
 std::string Rootfile1name; // = "rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";
                                          
 if(input_Int ==0 ){ Rootfile1name = "stv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";}
 else if(input_Int ==1 ){ Rootfile1name = "stv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root";}
  else if(input_Int ==2 ){ Rootfile1name = "stv-data_bnb_mcc9.1_v08_00_00_25_reco2_G1_beam_good_reco2_1e19.root";}
   else if(input_Int ==3 ){ Rootfile1name = "stv-data_extbnb_mcc9.1_v08_00_00_25_reco2_D_E_all_reco2.root";}
  else{Rootfile1name = "rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";}




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

  std::vector< double > p_edges = {0.1,0.2,0.275,0.335,0.46,0.61,0.785, 0.985, 1.185, 1.385, 1.585 }; //TODO: Rough measurements based on the "resolution" plot.
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
  
//std::vector< double > p_proton_edges = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87,0.93, 1.2}; 
   std::vector< double > p_proton_edges = { 0.25, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.2};
    


     
     
 std::vector< double > cos_edges_proton = {-1.0, -0.5, 0.0, 0.27, 0.45, 0.62, 0.76, 0.86,0.94, 1.0};

  
  /////////////////////////////////ProtonMupl/////////////////////////////////////////////////
  TH1D* h_p             = new TH1D("h_p", ";p_{#mu}^{reco}", p_edges.size() - 1, p_edges.data());
    h_p->SetDirectory(0); h_p->Sumw2();
  TH1D* h_leadingProton_p             = new TH1D("h_leadingProton_p", ";p_{p}^{reco}", p_proton_edges.size() - 1, p_proton_edges.data());
     h_leadingProton_p->SetDirectory(0); h_leadingProton_p->Sumw2();
  TH1D* h_costheta_proton             = new TH1D("h_costheta_proton", ";p_{#mu}^{reco}", cos_edges_proton.size() - 1, cos_edges_proton.data());
   h_costheta_proton->SetDirectory(0); h_costheta_proton->Sumw2();
  
  
 HistFolio<TH1D, FidVol_Category> h_p_containment_categories =
 HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_p_containment_categories", p_edges, "h_p_containment_categories");

HistFolio<TH1D, CCZeroPi_type> h_p_topology_categories =
  HistFolio<TH1D, CCZeroPi_type>(topology_categories, "h_p_topogoly_categories", p_edges, "h_p_topogoly_categories");



  

    ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_costheta_p             = new TH1D("h_costheta_p", ";p_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data());
   h_costheta_p->SetDirectory(0); h_costheta_p->Sumw2();
   
      
  HistFolio<TH1D, FidVol_Category> h_costheta_p_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_costheta_p_containment_categories", cos_edges, "h_costheta_p_containment_categories");
  /////////////////////////////////////////////////////////////////////////  
     TH2D* h_p_costheta_p  = new TH2D("h_p_costheta_p", "h_p_costheta_p",  cos_edges.size() - 1, cos_edges.data(), p_edges.size() - 1, p_edges.data());
   h_p_costheta_p->SetDirectory(0);
   h_p_costheta_p->Sumw2();
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_pn             = new TH1D("h_pn", ";pn_{#mu}^{reco}", pn_edges.size() - 1, pn_edges.data());
  
  h_pn->SetDirectory(0);            h_pn->Sumw2();
  
 
  HistFolio<TH1D, FidVol_Category> h_pn_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_pn_containment_categories", pn_edges, "h_pn_containment_categories");
    
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_alphaT             = new TH1D("h_delta_alphaT", ";#delta_{#alpha T}^{reco}", delta_alphaT_edges.size() - 1, delta_alphaT_edges.data());
  
  h_delta_alphaT->SetDirectory(0);            h_delta_alphaT->Sumw2();
  

  HistFolio<TH1D, FidVol_Category> h_delta_alphaT_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_delta_alphaT_containment_categories", delta_alphaT_edges, "h_delta_alphaT_containment_categories");
    
  
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_pTx             = new TH1D("h_delta_pTx", ";#delta_{#alpha Tx}^{reco}", delta_pTx_edges.size() - 1, delta_pTx_edges.data());
  h_delta_pTx->SetDirectory(0);            h_delta_pTx->Sumw2();
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_pTy             = new TH1D("h_delta_pTy", ";#delta_{#alpha T}^{reco}", delta_pTy_edges.size() - 1, delta_pTy_edges.data());
  h_delta_pTy->SetDirectory(0);            h_delta_pTy->Sumw2();
  

  HistFolio<TH1D, FidVol_Category> h_delta_pTy_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_delta_pTy_containment_categories", delta_pTy_edges, "h_delta_pTy_containment_categories");

////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_delta_phiT             = new TH1D("h_delta_phiT", "", delta_phiT_edges.size() - 1, delta_phiT_edges.data());
  
  h_delta_phiT->SetDirectory(0);            h_delta_phiT->Sumw2();
   
 
  HistFolio<TH1D, FidVol_Category> h_delta_phiT_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_delta_phiT_containment_categories", delta_phiT_edges, "h_delta_phiT_containment_categories");
  
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_openingAngle             = new TH1D("h_openingAngle", "", openingAngle_edges.size() - 1, openingAngle_edges.data());
  h_openingAngle->SetDirectory(0);            h_openingAngle->Sumw2();
   
  HistFolio<TH1D, FidVol_Category> h_openingAngle_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_openingAngle_containment_categories", openingAngle_edges, "h_openingAngle_containment_categories");
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_openingAngle_deg             = new TH1D("h_openingAngle_deg", "", openingAngleDeg_edges.size() - 1, openingAngleDeg_edges.data());
  h_openingAngle_deg->SetDirectory(0);            h_openingAngle_deg->Sumw2();
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexX             = new TH1D("h_VertexX", "", VertexX_edges.size() - 1, VertexX_edges.data());
  h_VertexX->SetDirectory(0);            h_VertexX->Sumw2();
    
  HistFolio<TH1D, FidVol_Category> h_VertexX_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_VertexX_containment_categories", VertexX_edges, "h_VertexX_containment_categories");
  
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexY             = new TH1D("h_VertexY", "", VertexY_edges.size() - 1, VertexY_edges.data());
  h_VertexY->SetDirectory(0);            h_VertexY->Sumw2();
    

  HistFolio<TH1D, FidVol_Category> h_VertexY_containment_categories =
  HistFolio<TH1D, FidVol_Category>(ContainedGroup_categories, "h_VertexY_containment_categories", VertexY_edges, "h_VertexY_containment_categories");
  
////////////////////////////////////////////////////////////////////////////////////////
  TH1D* h_VertexZ             = new TH1D("h_VertexZ", "", VertexZ_edges.size() - 1, VertexZ_edges.data());
  
  h_VertexZ->SetDirectory(0);            h_VertexZ->Sumw2();
  
  
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
    std::cout << "File: " << infile << std::endl;

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
  
  float spline_weight, tuned_cv_weight;
  float pn, delta_alphaT, delta_pTx, delta_pTy, delta_pL, delta_phiT, delta_pT;
  float reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z;
 


  stv->SetBranchStatus("*", 0);
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


  //stv->SetBranchStatus("n_pfps",1); 
  stv->SetBranchStatus("pfpdg",1);


  stv->SetBranchAddress("sel_CCNp0pi", &sel_signal);
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

  //stv->SetBranchAddress("n_pfps", &n_pfps);
  stv->SetBranchAddress("pfpdg", &pfpdg);





  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout<< "~~~~~~~~~~~~~~~ STarting Loop~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout << " The Number of Entries "<<  stv->GetEntries() << std::endl;
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
  std::cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;


  for (long i=0; i<stv->GetEntries(); i++)
  {
    stv->GetEntry(i);
   // if ( sel_signal&& mc_signal){
     if ( i % 10000 == 0 ) {
      std::cout << "Processing event #" << i << '\n';
    }
     
      if (sel_signal){
   
//std::cout<<"pn = " << pn << std::endl;

FidVol_Category Containmenttype = EventCategory_tool.GetContainmentType_ForData(sel_muon_contained);


//double Wgt = 1.0;
double Wgt = spline_weight * tuned_cv_weight;
//std::cout<< "Wgt = " << Wgt << std::endl;
double Pmu = p3_mu->Mag();

double CosTheta = p3_mu->CosTheta();
double P_proton = p3_lead_p->Mag();
double P_CosTheta = p3_lead_p->CosTheta();

double OpenAngle = TMath::ACos( (p3_mu->X()*p3_lead_p->X() + p3_mu->Y()*p3_lead_p->Y() + p3_mu->Z()*p3_lead_p->Z()) / p3_mu->Mag());
double OpenAngle_deg = OpenAngle * TMath::RadToDeg();

 // assuming the number of tracks is the Proton Multipleicity + muon track

h_p->Fill(Pmu,Wgt); 




h_costheta_p->Fill(CosTheta,Wgt);

h_p_containment_categories.GetComponentHist(Containmenttype)->Fill(Pmu,Wgt);

h_costheta_p_containment_categories.GetComponentHist(Containmenttype)->Fill(CosTheta,Wgt);

h_p_costheta_p->Fill(CosTheta,Pmu,Wgt); 
 




h_VertexX->Fill(reco_nu_vtx_sce_x,Wgt);  
h_VertexY->Fill(reco_nu_vtx_sce_y,Wgt);  
h_VertexZ->Fill(reco_nu_vtx_sce_z,Wgt);

h_VertexX_Z->Fill(reco_nu_vtx_sce_z,reco_nu_vtx_sce_x,Wgt); 
h_VertexY_Z->Fill(reco_nu_vtx_sce_z,reco_nu_vtx_sce_y,Wgt); 
h_VertexX_Y->Fill(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,Wgt); 



h_openingAngle->Fill(OpenAngle,Wgt);  
h_openingAngle_deg->Fill(OpenAngle_deg,Wgt);  
h_pn->Fill(pn,Wgt); 
h_delta_alphaT->Fill(delta_alphaT,Wgt); 
h_delta_pTx->Fill(delta_pTx,Wgt);  
h_delta_pTy->Fill(delta_pTy,Wgt);  
h_delta_phiT->Fill(delta_phiT,Wgt);  
h_leadingProton_p->Fill(P_proton,Wgt);
h_costheta_proton->Fill(P_CosTheta,Wgt);


}// ENd of If Statement, that passed Signal Selection 
//////////////////////////////////////////////////////  
}// End of Event Loop
///////////////////////////////////////////////////////


  f->Close();

}


/////////////////////////
// Make Root File to look 
/////////////////////////
std::string OutPulRootName = "/uboone/data/users/cnguyen/CC0Pi_Selection/CCNP_0Pion_steven_test/DATA_Selected_" + Rootfile1name;
std::cout << "OutPutRoot file name = "<< OutPulRootName<<std::endl;


auto outFile = TFile::Open(OutPulRootName.c_str(), "RECREATE");
outFile->cd(); 

h_p->Write();
h_costheta_p->Write();
h_p_costheta_p->Write();
h_pn->Write();


h_p_containment_categories.WriteToFile(*outFile);




h_delta_alphaT->Write();


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

h_leadingProton_p->Write();
h_costheta_proton->Write();

outFile->Close();

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

