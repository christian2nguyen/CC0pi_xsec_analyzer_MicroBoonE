/*Author: P. Englezos <p.englezos@physics.rutgers.edu>
 *
 * Usage: ./migration_matrix
 *
 */

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
#include "includes/EventCategory.hh"
#include "FiducialVolume.hh"
#include "includes/TreeUtils.hh"
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


int main() {

  TH2D* h_matr = new TH2D("h_matr", ";true bin; reconstructed bin; ", 45, 0, 45, 45, 0, 45);
  h_matr->SetDirectory(0);
  h_matr->Sumw2();

  float tuned_weight;
  bool mc_signal, sel_signal;
  TVector3 *p3_mu, *mc_p3_mu;
  std::vector<float> p_edges{0.1, 0.17, 0.21, 0.24, 0.27, 0.30, 0.38, 0.48, 0.75, 1.14, 2.5, INFINITY};
  std::vector<std::vector<float>> cos_edges{
  {-1,0,1},
  {-1,-0.2,0.4,1},
  {-1,-0.2,0.4,1},
  {-1,-0.1,0.5,1},
  {-1,-0.1,0.35,0.6,1},
  {-1,-0.4,-0.1,0.1,0.35,0.5,0.7,0.85,1},
  {-1,0,0.5,0.65,0.8,0.92,1},
  {-1,0.2,0.5,0.65,0.8,0.88,0.95,1},
  {-1,0.5,0.8,0.88,0.95,1},
  {-1,0.85,0.9,0.95,1},
  {-1,1}
  };    

/*
std::map< double, std::vector<double> > mapcheck{
{0.1,{-1,0,1}},
{0.17,{-1,-0.2,0.4,1}},
{0.21,{-1,-0.2,0.4,1}},
{0.24,{-1,-0.1,0.5,1}},
{0.27,{-1,-0.1,0.35,0.6,1}},
{0.30,{-1,-0.4,-0.1,0.1,0.35,0.5,0.7,0.85,1}},
{0.38,{-1,0,0.5,0.65,0.8,0.92,1}},
{0.48,{-1,0.2,0.5,0.65,0.8,0.88,0.95,1}},
{0.75,{-1,0.5,0.8,0.88,0.95,1}},
{1.14,{-1,1}},
{2.5,{-1,1}},
{INFINITY,{}}
};*/





  TFile* file = TFile::Open("/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_12_26_24/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root");
  TTree* stv = (TTree*) file->Get("stv_tree");

  stv->SetBranchStatus("*", 0);
  stv->SetBranchStatus("mc_is_cc0pi_signal", 1);
  stv->SetBranchStatus("sel_CC0pi", 1);
  stv->SetBranchStatus("mc_p3_mu", 1);
  stv->SetBranchStatus("p3_mu", 1);
  stv->SetBranchStatus("tuned_cv_weight", 1);

  stv->SetBranchAddress("mc_is_cc0pi_signal", &mc_signal);
  stv->SetBranchAddress("sel_CC0pi", &sel_signal);
  stv->SetBranchAddress("p3_mu", &p3_mu);
  stv->SetBranchAddress("mc_p3_mu", &mc_p3_mu);
  stv->SetBranchAddress("tuned_cv_weight", &tuned_weight);

  for (long i=0; i<stv->GetEntries(); i++) {
      stv->GetEntry(i);

      if (sel_signal && mc_signal){ 

         int true_bin = 0;
         
         for (int j=0; j<cos_edges.size(); j++) { 

             for (int k=0; k<cos_edges[j].size()-1; k++) {  

                  if (p_edges[j]<mc_p3_mu->Mag()  && mc_p3_mu->Mag()<p_edges[j+1] && cos_edges[j][k]<mc_p3_mu->CosTheta() && mc_p3_mu->CosTheta()<cos_edges[j][k+1]){ //true_bin
 
                      int reco_bin = 0; 
                      
                      for (int l=0; l<cos_edges.size(); l++) {

                           for (int m=0; m<cos_edges[l].size()-1; m++) {

                                if (p_edges[l]<p3_mu->Mag()  && p3_mu->Mag()<p_edges[l+1] && cos_edges[l][m]<p3_mu->CosTheta() && p3_mu->CosTheta()<cos_edges[l][m+1]){ //reco_bin

                                    h_matr->Fill(true_bin,reco_bin,tuned_weight);
                                 }
                                 reco_bin += 1; 
                           }
                       }
                   }
                   true_bin += 1;
               }
          }
      }
  }

  for (int a=1; a<46; a++) {
       
      int bincon=0;      
 
      for (int b=1; b<46; b++) {

           bincon += h_matr->GetBinContent(a,b);
            
      }
       
       for (int c=1; c<46; c++) {

           h_matr->SetBinContent(a,c,(h_matr->GetBinContent(a,c)/bincon));
       }  

  }  

  TCanvas* c = new TCanvas;
  gStyle->SetPalette(60);
  h_matr->SetStats(0);
  h_matr->Draw("colz");
  c->SaveAs("migration_matrix.pdf");

  return 0;
}

