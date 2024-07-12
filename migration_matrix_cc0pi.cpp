/*1st Author: P. Englezos <p.englezos@physics.rutgers.edu>
 *2nd Christian 
 * Usage: ./migration_matrix
 * Modified version of Panos the nested loops are SLLLOOWW
 */
// Code to make 2d migration using steven-like binning or inclusive-like binning 



#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include <cassert>
#include <set>
#include <vector>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TPad.h>

#include "includes/EventCategory.hh"
#include "FiducialVolume.hh"
#include "includes/TreeUtils.hh"


#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "includes/PlotUtils.hh"
#include "includes/SliceBinning.hh"
#include "SliceHistogram.hh"
#include "TLatex.h"
#include "TLine.h"
#include "includes/GridCanvas.hh"
#include "TH2Poly.h"
#include "includes/UBTH2Poly.h" 
#include "includes/HistUtils_cc0pi.hh"
//#include "HistUtils.hh"
#include "ConfigMakerUtils.hh"



int main() {


//std::vector<std::vector<double>> MakeVectorofVectorsFromMap(std::map< double, std::vector<double> > inputMap)
//std::vector<double> Vector_Biningkeys = MakeKeyVectorFromMap(std::map< double, std::vector<double> > inputMap)
/*
std::map< double, std::vector<double> > mapcheck{
  {0.1,{-1,0,1}},
 { 0.17, {-1,-0.2,0.4,1}},
  {0.21,{-1,-0.2,0.4,1}},
  {0.24,{-1,-0.1,0.5,1}},
  { 0.27,{-1,-0.1,0.35,0.6,1}},
  {0.30,{-1,-0.4,-0.1,0.1,0.35,0.5,0.7,0.85,1}},
  {0.38,{-1,0,0.5,0.65,0.8,0.92,1}},
  {.48,{-1,0.2,0.5,0.65,0.8,0.88,0.95,1}},
  {.75,{-1,0.5,0.8,0.88,0.95,1}},
  {1.14,{-1,0.85,0.9,0.95,1}},
  {2.5,{-1,1}},
  {INFINITY,{}}


};
*/
    char text_title_pdf1[2024];
    char text_title_pdf2[2024];
    char text_title_pdf3[2024];
    char text_title_pdf4[2024];
 


int nbins = NBinsFromMap(MUON_2D_BIN_EDGES);
int nbins_inclusive = NBinsFromMap(MUON_2D_BIN_EDGES_inclusive);
std::cout <<"nbins = "<< nbins << std::endl;
std::cout <<"nbins_inclusive = "<< nbins_inclusive << std::endl;


  TH2D* h_matr = new TH2D("Muon", ";true bin; reconstructed bin; ", nbins, 0, nbins, nbins, 0, nbins);
  TH2D* h_matr_inclusive = new TH2D("Muon_inclusive", ";true bin; reconstructed bin; ", nbins_inclusive, 0, nbins_inclusive, nbins_inclusive, 0, nbins_inclusive); 
  h_matr->SetDirectory(0);
  h_matr->Sumw2();
  h_matr_inclusive->SetDirectory(0);
  h_matr_inclusive->Sumw2();
  
  

  float tuned_cv_weight,spline_weight;
  bool mc_signal, sel_CC0pi, Muon_BDT_SCORE_above_threshold, sel_Muon_BDT_SCORE_above_threshold;
  bool sel_cosmic_ip_cut_passed, sel_BDT_NotBOGUS, sel_All_Protons_BDT_SCORE_above_threshold;
  TVector3 *p3_mu, *mc_p3_mu;
  //std::vector<float> p_edges{0.1,0.17,0.21,0.24,0.27,0.30,0.38,0.48,0.75,1.14,2.5,INFINITY};
  //std::vector<std::vector<float>> cos_edges{{-1,0,1},{-1,-0.2,0.4,1},{-1,-0.2,0.4,1},{-1,-0.1,0.5,1},{-1,-0.1,0.35,0.6,1},{-1,-0.4,-0.1,0.1,0.35,0.5,0.7,0.85,1},{-1,0,0.5,0.65,0.8,0.92,1},{-1,0.2,0.5,0.65,0.8,0.88,0.95,1},{-1,0.5,0.8,0.88,0.95,1},{-1,0.85,0.9,0.95,1},{-1,1}};    



std::vector<double> p_edges = MakeKeyVectorFromMap(MUON_2D_BIN_EDGES);
std::vector<std::vector<double>> cos_edges = MakeVectorofVectorsFromMap(MUON_2D_BIN_EDGES);


std::vector<double> cos_edges_inclusive = MakeKeyVectorFromMap(MUON_2D_BIN_EDGES_inclusive);
std::vector<std::vector<double>> p_edges_inclusive = MakeVectorofVectorsFromMap(MUON_2D_BIN_EDGES_inclusive);


  TFile* file = TFile::Open("/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_2_2_24/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root");
  TTree* stv = (TTree*) file->Get("stv_tree");

  stv->SetBranchStatus("*", 0);
  stv->SetBranchStatus("mc_is_cc0pi_signal", 1);
  stv->SetBranchStatus("sel_CC0pi", 1);
  stv->SetBranchStatus("mc_p3_mu", 1);
  stv->SetBranchStatus("p3_mu", 1);
  stv->SetBranchStatus("tuned_cv_weight", 1);
  stv->SetBranchStatus("spline_weight", 1);
  stv->SetBranchStatus("sel_Muon_BDT_SCORE_above_threshold", 1);
  stv->SetBranchStatus("sel_cosmic_ip_cut_passed", 1);
  stv->SetBranchStatus("sel_BDT_NotBOGUS", 1);
  stv->SetBranchStatus("sel_All_Protons_BDT_SCORE_above_threshold", 1);



  stv->SetBranchAddress("mc_is_cc0pi_signal", &mc_signal);
  stv->SetBranchAddress("sel_CC0pi", &sel_CC0pi);
  stv->SetBranchAddress("p3_mu", &p3_mu);
  stv->SetBranchAddress("mc_p3_mu", &mc_p3_mu);
  stv->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
  stv->SetBranchAddress("spline_weight", &spline_weight);
  
  stv->SetBranchAddress("sel_Muon_BDT_SCORE_above_threshold" , &sel_Muon_BDT_SCORE_above_threshold); 
  stv->SetBranchAddress("sel_cosmic_ip_cut_passed" , &sel_cosmic_ip_cut_passed);   
  stv->SetBranchAddress("sel_BDT_NotBOGUS" , &sel_BDT_NotBOGUS);   
  stv->SetBranchAddress("sel_All_Protons_BDT_SCORE_above_threshold" , &sel_All_Protons_BDT_SCORE_above_threshold);


  std::map<int , BinMap> TH1Poly_binMap_test;
  std::map<int , BinMap> TH1Poly_binMap_test2;

  std::map<int , BinMap> TH1Poly_binMap_test3;
  std::map<int , BinMap> TH1Poly_binMap_test4;
  
   UBTH2Poly instance;
   gInterpreter->Declare("#include \"includes/UBTH2Poly.h\""); 
  
  std::cout<<"here "<< std::endl;

  UBTH2Poly *TH2Poly_costheta_Pmu_TRUE = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "TH2Poly_costheta_Pmu_TRUE");
  UBTH2Poly *TH2Poly_costheta_Pmu_RECO = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test2, "TH2Poly_costheta_Pmu_RECO");

  UBTH2Poly *TH2Poly_costheta_Pmu_incluive_TRUE = Make2DHist_inclusive_UB( MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test3, "TH2Poly_costheta_Pmu_incluive_TRUE");
  UBTH2Poly *TH2Poly_costheta_Pmu_incluive_RECO = Make2DHist_inclusive_UB( MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test4, "TH2Poly_costheta_Pmu_incluive_RECO");

  std::cout<<"here 2"<< std::endl;
  
  
  
  for (long i=0; i<stv->GetEntries(); i++) {
      stv->GetEntry(i);
      std::cout<<"i = "<< i << std::endl;
     
     if ( i % 10000 == 0 ) {
      std::cout << "Processing event #" << i << '\n';
    }
    
    
      if ( sel_CC0pi && mc_signal && sel_Muon_BDT_SCORE_above_threshold  && 
      sel_cosmic_ip_cut_passed  && sel_BDT_NotBOGUS &&
      sel_All_Protons_BDT_SCORE_above_threshold){ 

         int true_bin = 0;
         
         Double_t Wgt =   tuned_cv_weight; // spline_weight *
         
         Double_t TRUE_p3_mag = mc_p3_mu->Mag();
         Double_t RECO_p3_mag = p3_mu->Mag();
         Double_t TRUE_THETA = mc_p3_mu->CosTheta();
         Double_t RECO_THETA = p3_mu->CosTheta();
         
        Int_t TRUE_bin_ = TH2Poly_costheta_Pmu_TRUE->Fill(TRUE_THETA,TRUE_p3_mag, Wgt);   
        Int_t RECO_bin_ = TH2Poly_costheta_Pmu_RECO->Fill(RECO_THETA,RECO_p3_mag, Wgt); 

        Int_t TRUE_bin2_ = TH2Poly_costheta_Pmu_incluive_TRUE->Fill(TRUE_THETA,TRUE_p3_mag, Wgt);   
        Int_t RECO_bin2_ = TH2Poly_costheta_Pmu_incluive_RECO->Fill(RECO_THETA,RECO_p3_mag, Wgt); 




       h_matr->Fill(TRUE_bin_, RECO_bin_, Wgt);
       h_matr_inclusive->Fill(TRUE_bin2_, RECO_bin2_, Wgt);   

      
          
      }/// End of Selection
  }// ENd of Entries 



std::string OutPutRootName = "/uboone/data/users/cnguyen/CC0Pi_Selection/MC/MC_Selected_2DMigration.root";
std::cout<<" Writing Root File :"<< OutPutRootName << std::endl;

auto outFile = TFile::Open(OutPutRootName.c_str(), "RECREATE");
outFile->cd(); 


h_matr->Write("Migration_BeforeNorm");
h_matr_inclusive->Write("Migration_BeforeNorm_inclusive");


  for (int a=1; a<nbins; a++) {
       
      int bincon=0;      
 
      for (int b=1; b<nbins; b++) {

           bincon += h_matr->GetBinContent(a,b);
            
      }
       
       for (int c=1; c<nbins; c++) {

           h_matr->SetBinContent(a,c,(h_matr->GetBinContent(a,c)/bincon));
       }  

  }  
  
  
  
  
  
  
    for (int a=1; a<nbins_inclusive; a++) {
       
      int bincon=0;      
 
      for (int b=1; b<nbins_inclusive; b++) {

           bincon += h_matr_inclusive->GetBinContent(a,b);
            
      }
       
       for (int c=1; c<nbins; c++) {

           h_matr_inclusive->SetBinContent(a,c,(h_matr->GetBinContent(a,c)/bincon));
       }  

  }  
  
  h_matr->Write("MigrationeNorm");
  h_matr_inclusive->Write("Migration_Norm_inclusive");
  outFile->Close();
  
   TCanvas* can = new TCanvas;
  sprintf(text_title_pdf1, "migration_matrix_cc0pi.pdf(" );
  can -> Print(text_title_pdf1);
  sprintf(text_title_pdf2, "migration_matrix_cc0pi.pdf" );
  sprintf(text_title_pdf3, "migration_matrix_cc0pi.pdf)" );
  sprintf(text_title_pdf4, "migration_matrix_cc0pi.pdf");
  

 
  gStyle->SetPalette(60);
  h_matr->SetStats(0);
  h_matr->Draw("colz");
 can -> Print(text_title_pdf2);


  h_matr_inclusive->SetStats(0);
  h_matr_inclusive->Draw("colz");
 can -> Print(text_title_pdf2);


 can -> Print(text_title_pdf3);
 can->Close(); 

  return 0;
}

