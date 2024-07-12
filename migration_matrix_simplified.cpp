/*Author: P. Englezos <p.englezos@physics.rutgers.edu>
 *
 * Usage: ./migration_matrix
 *
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

void migration_matrix(std::string infiles) {

 std::vector<std::string> filenames;
 std::ifstream myfile(infiles);
 std::copy(std::istream_iterator<std::string>(myfile),
            std::istream_iterator<std::string>(),
            std::back_inserter(filenames));
 std::cout << "File Count: " << filenames.size() << std::endl;

  TH2D* h_matr_pmult = new TH2D("h_matr_pmult", ";true bin; reconstructed bin; ", 5, -1., 4., 5, -1., 4.);
  h_matr_pmult->SetDirectory(0);
  h_matr_pmult->Sumw2();

  std::vector< double > p_edges = {0.1,0.2,0.275,0.335,0.46,0.61,0.785, 0.985, 1.185, 1.385, 1.585 }; //TODO: Rough measurements based on the "resolution" plot.

  std::vector< double >  cos_edges = {-1.,-0.96,-0.91,-0.87,-0.82,-0.77,-0.72,-0.67,-0.62,-0.58,-0.53,-0.48,-0.44,-0.4,-0.36,-0.32,-0.28,-0.24,-0.20,-0.16,-0.12,-0.08,-0.04,0.,0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96,1.};
 
  TH2D* h_matr_pmu  = new TH2D("h_matr_pmu", ";p_{#mu}^{true} (GeV); p_{#mu}^{reco} (GeV)", p_edges.size() - 1, p_edges.data(), p_edges.size() - 1, p_edges.data()); 
  //TH2D* h_matr_pmu  = new TH2D("h_matr_pmu", ";p_{#mu}^{reco} - p_{#mu}^{true}; inverse", 600, -3., 3., 600, -3., 3.);
  h_matr_pmu->SetDirectory(0);
  h_matr_pmu->Sumw2();

  TH2D* h_matr_costh  = new TH2D("h_matr_costh", ";cos#theta_{#mu}^{true}; cos#theta_{#mu}^{reco}", cos_edges.size() - 1, cos_edges.data(), cos_edges.size() - 1, cos_edges.data()); 
  h_matr_costh->SetDirectory(0);
  h_matr_costh->Sumw2();

  TH2D* h_matr_pres  = new TH2D("h_matr_pres", ";p_{#mu}^{reco} (GeV); p_{#mu}^{reco} - p_{#mu}^{true}", 240, 0.1, 2.5, 600, -3., 3.);
  h_matr_pres->SetDirectory(0);
  h_matr_pres->Sumw2();

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
  std::vector<int> *mc_pdg = nullptr;
  float  spline_weight, tuned_cv_weight;

  stv->SetBranchStatus("*", 0);
  stv->SetBranchStatus("mc_is_cc0pi_signal", 1);
  stv->SetBranchStatus("sel_CC0pi", 1);
  stv->SetBranchStatus("mc_p3_mu", 1);
  stv->SetBranchStatus("p3_mu", 1);
  stv->SetBranchStatus("mc_pdg", 1);
  stv->SetBranchStatus("sel_muon_contained", 1);
  stv->SetBranchStatus("spline_weight", 1);
  stv->SetBranchStatus("tuned_cv_weight", 1);

  stv->SetBranchAddress("mc_is_cc0pi_signal", &mc_signal);
  stv->SetBranchAddress("sel_CC0pi", &sel_signal);
  stv->SetBranchAddress("p3_mu", &p3_mu);
  stv->SetBranchAddress("mc_p3_mu", &mc_p3_mu);
  stv->SetBranchAddress("mc_pdg", &mc_pdg);
  stv->SetBranchAddress("sel_muon_contained", &sel_muon_contained);
  stv->SetBranchAddress("spline_weight", &spline_weight);
  stv->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);

  std::cout << stv->GetEntries() << std::endl;

  for (long i=0; i<stv->GetEntries(); i++) {
      stv->GetEntry(i);

      if (sel_signal && mc_signal){

         //double p_reco_new = p3_mu->Mag() - (0.685944 - 9.35117*1e-6*std::pow(67367*std::sqrt(4.53831*1e21*resol*resol + 7.12411*1e20*resol + 6.05908*1e19) - 4.53831*1e15*resol - 3.56206*1e14,1./3.) + 49474.6/(std::pow(67367*std::sqrt(4.53831*1e21*resol*resol + 7.12411*1e20*resol + 6.05908*1e19) - 4.53831*1e15*resol - 3.56206*1e14,1./3.)));
  
//Best Combo (but it uses resol):      
        //double p_reco_new = p3_mu->Mag() - 0.134734*std::pow(p3_mu->Mag(),3) + 0.27726*std::pow(p3_mu->Mag(),2) - 0.377187*p3_mu->Mag() + 0.250247;
        //if ((resol<0. && p3_mu->Mag()>0.7) || (resol > 0. && p3_mu->Mag() < 1.)) p_reco_new = p3_mu->Mag() + 0.134734*std::pow(p3_mu->Mag(),3) - 0.27726*std::pow(p3_mu->Mag(),2) + 0.377187*p3_mu->Mag() - 0.250247;
        //double resol = p3_mu->Mag() - mc_p3_mu->Mag();
        double p_reco_new = -1.;       
        //double p_reco_new = p3_mu->Mag(); 

       if (!sel_muon_contained){
           p_reco_new = p3_mu->Mag() + (- 0.146795*std::pow(p3_mu->Mag(),3) + 0.318775*std::pow(p3_mu->Mag(),2) - 0.405845*p3_mu->Mag() + 0.241603)*spline_weight*tuned_cv_weight;
        }
        else{
           p_reco_new = p3_mu->Mag() - 0.02901*std::pow(p3_mu->Mag(),3) + 0.0847543*std::pow(p3_mu->Mag(),2) - 0.111323*p3_mu->Mag() + 0.0632452;
        } 
       
        double resol_new = p_reco_new - mc_p3_mu->Mag();                
        if (p3_mu->Mag() > 2. && mc_p3_mu->Mag() > 2.) std::cout<<"true: "<<mc_p3_mu->Mag()<<" reco_old: "<<p3_mu->Mag()<<" reco_new: "<<p_reco_new<<" spline: "<<spline_weight<<" tunes: "<<tuned_cv_weight<<std::endl;
  
         h_matr_pres->Fill(p_reco_new,resol_new, spline_weight*tuned_cv_weight);    

         h_matr_pmu->Fill(mc_p3_mu->Mag(), p_reco_new);
         //h_matr_pmu->Fill(resol,p3_mu->Mag() - p_reco_new);

         bool  true_flag_cos = false;
         
           for (int true_index_cos=0; true_index_cos<cos_edges.size() - 1; true_index_cos++) {
                  
                  if (cos_edges[true_index_cos]<mc_p3_mu->CosTheta()  && mc_p3_mu->CosTheta()<cos_edges[true_index_cos+1]){
                      
                      for (int reco_index_cos=0; reco_index_cos<cos_edges.size() - 1; reco_index_cos++) {
                                 
                                 if (cos_edges[reco_index_cos]<p3_mu->CosTheta()  && p3_mu->CosTheta()<cos_edges[reco_index_cos + 1]){
                                    
                                   // h_matr_costh->Fill(true_index_cos,reco_index_cos);
                                    h_matr_costh->Fill(mc_p3_mu->CosTheta(),p3_mu->CosTheta());
			            true_flag_cos = true;
                                    break;
                                 }
                                 else{
                                      
                                      continue;
                                 
                                 }
                       }
                   }
                   else{
                       
                       continue;
                  
                  }
             
             if (true_flag_cos == true) break;
          
          } 
 
         bool true_flag_pmult = false;  //fix the conditions!
                   
           for (int true_index_pmult=0; true_index_pmult < 5; true_index_pmult++) {
                       
                  if (true_index_pmult<mc_p3_mu->CosTheta()  && mc_p3_mu->CosTheta()<true_index_pmult+1){ //use the pdg variable
                      
                      for (int reco_index_pmult=0; reco_index_pmult < 5; reco_index_pmult++) {
                                 
                                 if (reco_index_pmult<p3_mu->CosTheta()  && p3_mu->CosTheta()<reco_index_pmult + 1){ //use the pdg variable
                                    
                                   // h_matr_pmult->Fill(true_index_pmult,reco_index_pmult);
                                    h_matr_pmult->Fill(mc_p3_mu->CosTheta(),p3_mu->CosTheta());
				    true_flag_pmult = true;
                                    break;
                                 }
                                 else{
                                      
                                      continue;
                                 
                                 }
                       }         
                   }    
                   else{            
                                    
                       continue;    
                                 
                  }              
                                      
             if (true_flag_pmult == true) break;
           
          }        

       }
    }

  f->Close();

  }

/*
 for (int true_pmu=0; true_pmu < p_edges.size() - 1; true_pmu++) {
       
      int column_sum_pmu=0;      
 
      for (int reco_pmu=0; reco_pmu < p_edges.size() - 1; reco_pmu++) {

           column_sum_pmu += h_matr_pmu->GetBinContent(true_pmu + 1, reco_pmu + 1);            
      }
       
       for (int rep_pmu=0; rep_pmu < p_edges.size() - 1; rep_pmu++) {

           h_matr_pmu->SetBinContent(true_pmu + 1, rep_pmu + 1,(h_matr_pmu->GetBinContent(true_pmu + 1, rep_pmu + 1)/column_sum_pmu));
           if (h_matr_pmu->GetBinContent(true_pmu + 1, rep_pmu + 1)<=0.5) h_matr_pmu->SetBinContent(true_pmu + 1, rep_pmu + 1, 0.);
       }  
  }
*/
  for (int true_cos=0; true_cos < cos_edges.size() - 1; true_cos++) {

      int column_sum_cos=0;

      for (int reco_cos=0; reco_cos < cos_edges.size() - 1; reco_cos++) {

           column_sum_cos += h_matr_costh->GetBinContent(true_cos + 1, reco_cos + 1);
      }

       for (int rep_cos=0; rep_cos < cos_edges.size() - 1; rep_cos++) {
           h_matr_costh->SetBinContent(true_cos + 1, rep_cos + 1,(h_matr_costh->GetBinContent(true_cos + 1, rep_cos + 1)/column_sum_cos));
       }
  }

  for (int true_mult=0; true_mult < 5; true_mult++) {

      int column_sum_mult=0;

      for (int reco_mult=0; reco_mult < 5; reco_mult++) {

           column_sum_mult += h_matr_pmult->GetBinContent(true_mult + 1, reco_mult + 1);
      }

       for (int rep_mult=0; rep_mult < 5; rep_mult++) {

           h_matr_pmult->SetBinContent(true_mult + 1, rep_mult + 1,(h_matr_pmult->GetBinContent(true_mult + 1, rep_mult + 1)/column_sum_mult));
      }
  }  

  TCanvas* c1 = new TCanvas;
  gStyle->SetPalette(60);
  h_matr_pmu->SetStats(0);
  h_matr_pmu->Draw("colz");
  h_matr_pmu->SetAxisRange(0., 1.,"Z"); 
  c1->SaveAs("migration_matrix_pmu.pdf");

  TCanvas* c4 = new TCanvas;
  gStyle->SetPalette(60);
  h_matr_pres->SetStats(0);
  h_matr_pres->Draw("colz");
  c4->SaveAs("migration_matrix_resolt_both_preco.pdf");

  /*TCanvas* c5 = new TCanvas;
  gStyle->SetPalette(60);
  h_matr_pres->FitSlicesX(0,1,240,0);
  TH1D *h_matr_pres_1 = (TH1D*) gDirectory->Get("h_matr_pres_1");
  h_matr_pres_1->SetStats(0);
  h_matr_pres_1->Draw();
  c5->SaveAs("migration_matrix_resolt_th1.pdf");*/

  /* TH2D* h_resolve_mean  = new TH2D("h_resolve_mean", ";p_{#mu}^{true}; Mean & STD", 240, 0.1, 2.5, 1450, -1.25, 0.2);
  h_resolve_mean->SetDirectory(0);
  h_resolve_mean->Sumw2("kFALSE");

  TH2D* h_resolve_std_up  = new TH2D("h_resolve_std_up", ";p_{#mu}^{true}; Mean & STD", 240, 0.1, 2.5, 1450, -1.25, 0.2);
  h_resolve_std_up->SetDirectory(0);
  h_resolve_std_up->Sumw2("kFALSE");

  TH2D* h_resolve_std_down  = new TH2D("h_resolve_std_down", ";p_{#mu}^{true}; Mean & STD", 240, 0.1, 2.5, 1450, -1.25, 0.2);
  h_resolve_std_down->SetDirectory(0);
  h_resolve_std_down->Sumw2("kFALSE"); */

  /*std::vector<double>* h_resolve_std = nullptr;
  std::vector<double>* h_resolve_mean = nullptr;
  std::vector<double>* h_resolve_mom  = nullptr;
  std::vector<double>* h_resolve_empty  = nullptr;*/

  double* h_resolve_std = new double[239]; //239
  double* h_resolve_mean = new double [239];
  double* h_resolve_mom = new double[239];

  for (int rep_1d=1; rep_1d < 240; rep_1d++) {  //240
       TH1D* temp_hist = h_matr_pres->ProjectionY("", rep_1d, rep_1d + 1, "e");
       double mean = temp_hist->GetMean();
       double std = temp_hist->GetStdDev();
       double true_mom = rep_1d*0.01+0.1;
       //if (std==0. && mean==0.) mean=-100.; 
       h_resolve_mean[rep_1d - 1] = mean;
       h_resolve_std[rep_1d - 1] = std;
       h_resolve_mom[rep_1d - 1] = true_mom;
      //h_resolve_empty->push_back(0.); 
  }
 
  TCanvas* c5 = new TCanvas; 
  //int reps = (int)h_resolve_mom->size();
  auto gr = new TGraphErrors(239,h_resolve_mom,h_resolve_mean,0,h_resolve_std);
  //auto gr = new TGraph(239,h_resolve_mom,h_resolve_mean); //239
  gr->SetTitle(";p_{#mu}^{reco} (GeV/c);Mean and StD of (p_{#mu}^{reco} - p_{#mu}^{true});");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(7);
  //gr->SetMarkerSize(3);
  auto xaxis = gr->GetXaxis();
  xaxis->SetLimits(0.1,2.5);
  gr->GetHistogram()->SetMinimum(-0.8);
  gr->GetHistogram()->SetMaximum(2.);
  gr->Draw("AP");

  /* ///Fitting
  gr->Fit("pol3");
  TF1 *linefit = gr->GetFunction("pol3");
  std::cout<<"Probability: "<<linefit->GetProb()<<std::endl;
  linefit->SetLineWidth(1);
  linefit->Draw("same"); */

  TLine* l = new TLine(0.1,0.,2.5,0.);
  l->SetLineColor(kRed);
  l->SetLineStyle(kDashed);
  l->Draw("same");

  //h_resolve_mean->SetStats(0);
  //h_resolve_mean->GetYaxis()->SetRange(-1.6,1.);
  //h_resolve_mean->Draw();
 // h_resolve_std_up->Draw("same");
 // h_resolve_std_down->Draw("same");
  c5->SaveAs("migration_matrix_resolt_logistics_muon_both_preco.pdf");

  TCanvas* c2 = new TCanvas;
  gStyle->SetPalette(60);
  h_matr_costh->SetStats(0);
  h_matr_costh->Draw("colz");
 // c2->SaveAs("migration_matrix_cos.pdf"); 

  TCanvas* c3 = new TCanvas;
  gStyle->SetPalette(60);
  h_matr_pmult->SetStats(0);
  h_matr_pmult->Draw("colz");
 // c3->SaveAs("migration_matrix_pmult.pdf");

}
int main(int argc, char* argv[]) {
   migration_matrix(argv[1]);
   return 0;
}

