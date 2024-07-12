/*Author: P. Englezos <p.englezos@physics.rutgers.edu>
 *
 * Usage: ./chi_square_cc0pi files.txt
 *
 */

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <set>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"
#include "includes/EventCategory.hh"
#include "includes/TreeUtils.hh"
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1D.h"
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
#include "TColor.h"
#include "TInterpreter.h"
#include <algorithm>
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "includes/PlotUtils.hh"
#include "TMatrixT.h"
#include "WienerSVDUnfolder.hh"
#include "FiducialVolume.hh"
#include <iomanip>
#include "TMatrixD.h"
#include "MatrixUtils.hh"

void chi_square_all_gens(std::string infiles) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  std::vector< double > nbins = {0.1,0.275,0.375,0.525,0.7,0.9,1.1,2.};

  std::cout<<"Now examining real data."<<std::endl;

  auto* mcc9 = new MCC9SystematicsCalculator("./output_files/test2.root", "systcalc.conf" ); 
  double total_pot = mcc9->total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV( total_pot );
  double num_Ar = num_Ar_targets_in_FV();
  double conv_factor = num_Ar * integ_flux;

  int num_true_signal_bins = 0;
  for ( int t = 0; t < mcc9->true_bins_.size(); ++t ) {
    const auto& tbin = mcc9->true_bins_.at( t );
    if ( tbin.type_ == kSignalTrueBin ) ++num_true_signal_bins;
  }
  
  TH1D* genie_cv_truth = mcc9->cv_universe().hist_true_.get();
  int num_true_bins = genie_cv_truth->GetNbinsX();
  TH1D* unfolded_events = dynamic_cast< TH1D* >(genie_cv_truth->Clone("unfolded_events") );
  unfolded_events->Reset();
  TH1D* unfolded_events_vals = new TH1D("unfolded_events_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
  const auto& fake_data_univ = mcc9->fake_data_universe();
  TH1D* fake_data_truth_hist = nullptr;
  fake_data_truth_hist = fake_data_univ->hist_true_.get(); 

  auto& cov_mat_map = *mcc9->get_covariances().release();
  auto* cov_mat = cov_mat_map.at( "total" ).cov_matrix_.get();

  std::unique_ptr< Unfolder > unfolder (new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv ) ); //Is this the correct unfolding function for our case?

  auto result = unfolder->unfold( *mcc9 );
  for ( int t = 0; t < num_true_bins; ++t ) { //num_true_bins is signal + background bins
      double evts = 0.;
      double error = 0.;
      if ( t < num_true_signal_bins ) { //num_true_signal_bins is just signal bins (as indicated by its name)
         evts = result.unfolded_signal_->operator()( t, 0 );
         error = std::sqrt( std::max(0., result.cov_matrix_->operator()( t, t )) );
      }
      unfolded_events->SetBinContent( t + 1, evts );
      unfolded_events->SetBinError( t + 1, error );
      unfolded_events_vals->SetBinContent( t + 1, evts*1e38/conv_factor/(nbins[t+1] - nbins[t]) );
      unfolded_events_vals->SetBinError( t + 1, error*1e38/conv_factor/(nbins[t+1] - nbins[t]) );
  }

  std::cout<<"true signal bins: "<<num_true_signal_bins<<std::endl;

  TCanvas* c1 = new TCanvas("c1");
  unfolded_events->SetStats( false );
  unfolded_events->SetLineColor( kBlack );
  unfolded_events->SetLineWidth( 3 );
  unfolded_events->GetXaxis()->SetRangeUser( 0, num_true_signal_bins );
  unfolded_events->GetYaxis()->SetRangeUser( 0, 25000 );
  genie_cv_truth->SetStats( false );
  genie_cv_truth->SetLineColor( kRed );
  genie_cv_truth->SetLineWidth( 3 );
  genie_cv_truth->SetLineStyle( 9 );
  unfolded_events->Draw( "e" );
  genie_cv_truth->Draw( "hist same" );
  fake_data_truth_hist->SetStats( false );
  fake_data_truth_hist->SetLineColor( kBlue );
  fake_data_truth_hist->SetLineWidth( 3 );
  fake_data_truth_hist->SetLineWidth( 3 );
  fake_data_truth_hist->SetLineStyle( 2 );
  fake_data_truth_hist->Draw( "hist same" );
  TLegend* lg = new TLegend( 0.15, 0.7, 0.3, 0.85 );
  lg->AddEntry( unfolded_events, "Unfolded Events", "l" );
  lg->AddEntry( genie_cv_truth, "uB tune", "l" );
  lg->AddEntry( fake_data_truth_hist, "Fake Data truth", "l" );
  lg->Draw( "same" );
  c1->SaveAs("overlay_tune_unfolded_truth.pdf"); 

  TMatrixD fake_data_truth( num_true_signal_bins, 1 );
  for ( int b = 0; b < num_true_signal_bins; ++b ) {
      double true_evts = fake_data_truth_hist->GetBinContent( b + 1 );
      fake_data_truth( b, 0 ) = true_evts ;
  }

 TMatrixD genie_cv_truth_vec( num_true_signal_bins, 1 );
  for ( int b = 0; b < num_true_signal_bins; ++b ) {
      double true_evts = genie_cv_truth->GetBinContent( b + 1 );
      genie_cv_truth_vec( b, 0 ) = true_evts ;
  }


  const TMatrixD& A_C = *result.add_smear_matrix_; 
  TMatrixD ac_truth( A_C, TMatrixD::kMult, fake_data_truth );
  fake_data_truth = ac_truth;

 TMatrixD genie_cv_temp( A_C, TMatrixD::kMult, genie_cv_truth_vec );
 genie_cv_truth_vec = genie_cv_temp;

 /* auto inverse_unfolded_cov_mat = invert_matrix( unfolded_cov_mat );
  auto& inv_unfolded_cov_mat = *inverse_unfolded_cov_mat;

  for ( int t = 0; t < unfolded_cov_mat.GetNcols() ; ++t ) {
      for ( int u = 0; u < unfolded_cov_mat.GetNrows() ; ++u ) {
          // std::cout<<" cov_mat: "<<unfolded_cov_mat(u,t)<<std::endl;
     }
  }

  for ( int t = 0; t < inverse_unfolded_cov_mat->GetNcols() ; ++t ) {
      for ( int u = 0; u < inverse_unfolded_cov_mat->GetNrows() ; ++u ) {
          inv_unfolded_cov_mat(u,t) = inv_unfolded_cov_mat(u,t)/1e72*std::pow((nbins[u+1] -  nbins[u])*conv_factor,2);
     }
  }

   int tot_nbins_real = unfolded_events->GetNbinsX();
   TH1D* unfolded_events_num_bins = new TH1D("unfolded_events_num_bins", ";Bin index;Scaled Events", tot_nbins_real, 0, tot_nbins_real);
   int idx_real = 1;
   for (int i=1; i<unfolded_events_num_bins->GetNbinsX()+1; i++) {
      unfolded_events_num_bins->SetBinContent(idx_real, unfolded_events->GetBinContent(i));
      unfolded_events_num_bins->SetBinError(idx_real, unfolded_events->GetBinError(i));
      idx_real++;
    } 
 */

  std::vector<std::string> filenames;
  std::ifstream myfile(infiles);
  std::copy(std::istream_iterator<std::string>(myfile),
            std::istream_iterator<std::string>(),
            std::back_inserter(filenames));
  std::cout << "File Count: " << filenames.size() << std::endl;

  TCanvas* c2 = new TCanvas("c2");
  TLegend *leg=new TLegend(0.7,0.7,0.9,0.9);

  std::cout<<"Now examining generator true data"<<std::endl; 
  
  for (int i = 0; i < filenames.size(); i++) {

  std::cout << "File: " << filenames.at(i) << std::endl;
  TFile* file = TFile::Open((filenames.at(i)).c_str());
  TTree *tree = (TTree*)file->Get("FlatTree_VARS");

  TH1D* h_muons_gen_no_num_bins = new TH1D("h_muons_gen_no_num_bins", "; True Muon Total Momentum (GeV/c); Number of Events", nbins.size() - 1, nbins.data());

  int nfsp, ninitp, mode, pdg[100], pdg_init[100], pdg_nu;
  char cc;
  float  px[100], py[100], pz[100];
  double fSF;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("nfsp", 1);
  tree->SetBranchStatus("cc", 1);
  tree->SetBranchStatus("pdg", 1);
  tree->SetBranchStatus("px", 1);
  tree->SetBranchStatus("py", 1);
  tree->SetBranchStatus("pz", 1);
  tree->SetBranchStatus("fScaleFactor", 1);
  tree->SetBranchStatus("PDGnu", 1);
  tree->SetBranchStatus("Mode", 1);
  tree->SetBranchStatus("pdg_init", 1);
  tree->SetBranchStatus("ninitp", 1);

  tree->SetBranchAddress("nfsp", &nfsp);
  tree->SetBranchAddress("cc", &cc);
  tree->SetBranchAddress("pdg", &pdg);
  tree->SetBranchAddress("px", &px);
  tree->SetBranchAddress("py", &py);
  tree->SetBranchAddress("pz", &pz);
  tree->SetBranchAddress("fScaleFactor", &fSF);
  tree->SetBranchAddress("PDGnu", &pdg_nu);
  tree->SetBranchAddress("Mode", &mode);
  tree->SetBranchAddress("pdg_init", &pdg_init);
   tree->SetBranchAddress("ninitp", &ninitp);

 for (int i_tree=0; i_tree<tree->GetEntries(); i_tree++) {//Loop over the entries.
       tree->GetEntry(i_tree);

       if ( !(int)cc || pdg_nu != 14 ) continue;//Cut which keeps only events with numu CC interaction
       
       bool flag_pi0 = true;
       bool flag_picharged = true;

       for (int j = 0; j < nfsp; j++) {//Loop over the final state particles.
            if (pdg[j] == 111){//If a neutral pion is found, the event is rejected.
               flag_pi0 = false;
               break;
            }
            if (abs(pdg[j]) == 211){//Search for charged pion.
                float picharged_ptot = sqrt(px[j]*px[j] + py[j]*py[j] + pz[j]*pz[j]);
                if (picharged_ptot > 0.07){//If charged pion's total momentum is over 70 MeV/c, the event is rejected.
                    flag_picharged = false;
                    break;
                 }
             }
        }

       if (!flag_pi0 || !flag_picharged) continue;//The event is rejected because one of the above cuts is not satisfied. 

       bool flag_additional_mu = false;
       double muon_fill = -1.;
       int muon_count = 0;
       for (int n = 0; n < nfsp; n++) {
           if (pdg[n] == 13){//Search for muon
           
               if (muon_count > 0){
                  std::cout<<"An additional muon has been found. Event is rejected?"<<std::endl;
                  flag_additional_mu = true;
                  break; 
               }
                  
               double imuon_ptot = sqrt(px[n]*px[n] + py[n]*py[n] + pz[n]*pz[n]); 
               if (imuon_ptot >= 0.1 && imuon_ptot <= 2.){
                  muon_fill = imuon_ptot;
                  muon_count += 1;              
               }
            }
       }
       if (!flag_additional_mu && muon_fill > -1.){
          h_muons_gen_no_num_bins->Fill(muon_fill);
       }
    }

    int num_bins_mat = A_C.GetNcols();
    TMatrixD hist_mat( num_bins_mat, 1 );
    for ( int r = 0; r < num_bins_mat; ++r ) {
       hist_mat( r, 0 ) = h_muons_gen_no_num_bins->GetBinContent( r + 1 );
    }
    TMatrixD hist_mat_transformed( A_C, TMatrixD::EMatrixCreatorsOp2::kMult,hist_mat );
    for ( int r = 0; r < num_bins_mat; ++r ) {
       double val = hist_mat_transformed( r, 0 );
       h_muons_gen_no_num_bins->SetBinContent( r + 1, val*fSF*40.*1e38/(nbins[r+1] - nbins[r]));
    } //fSF = GetEventHistogram()->Integral("width") * double(1E-38) / double(fNEvents); TODO: How is the width handled by fscalefactor?

    /*for (int p=1; p<h_muons_gen_no_num_bins->GetNbinsX() + 1; p++) {
    
        double init_val = h_muons_gen_no_num_bins->GetBinContent(p);
         h_muons_gen_no_num_bins->SetBinContent(p, init_val*fSF*40.*1e38/(nbins[p] - nbins[p-1]));  
    }*/
 
   int tot_nbins = h_muons_gen_no_num_bins->GetNbinsX(); 
   TH1D* h_muons_gen = new TH1D("h_muons_gen", ";Bin index;Scaled Events", tot_nbins, 0, tot_nbins);
   int idx = 1;
   for (int i=1; i<h_muons_gen->GetNbinsX()+1; i++) {
      h_muons_gen->SetBinContent(idx, h_muons_gen_no_num_bins->GetBinContent(i));
      idx++;
   }

   h_muons_gen_no_num_bins->SetStats(0);
   h_muons_gen_no_num_bins->SetLineWidth(2);

   if (i==0){
       h_muons_gen_no_num_bins->GetYaxis()->SetTitle("d#sigma/dp_{muon} (10^{-38} cm^{2}/GeV/Ar)");
       h_muons_gen_no_num_bins->SetAxisRange(0.,30.,"Y");
       h_muons_gen_no_num_bins->SetLineColor(800);
       h_muons_gen_no_num_bins->Draw();
       h_muons_gen_no_num_bins->SetDirectory(0);
   }
   else{
       if (i==1) h_muons_gen_no_num_bins->SetLineColor(7); 
       if (i==2) h_muons_gen_no_num_bins->SetLineColor(2);
       if (i==3) h_muons_gen_no_num_bins->SetLineColor(3);
       h_muons_gen_no_num_bins->Draw("same");
       h_muons_gen_no_num_bins->SetDirectory(0);
    } 
       if (i==0) leg->AddEntry(h_muons_gen_no_num_bins,"NEUT 5.0.4","f");
       if (i==1) leg->AddEntry(h_muons_gen_no_num_bins,"Genie 3.0.6","f");
       if (i==2) leg->AddEntry(h_muons_gen_no_num_bins,"Genie 2.12.10","f");
       if (i==3) leg->AddEntry(h_muons_gen_no_num_bins,"NuWro 19.02","f");
  
  /*  double chi_square = 0;
  for (int k=0; k<inv_unfolded_cov_mat.GetNrows(); k++) {
       for (int j=0; j<inv_unfolded_cov_mat.GetNcols(); j++) {
           chi_square += (unfolded_events_num_bins->GetBinContent(k)-h_muons_gen->GetBinContent(k))*(inv_unfolded_cov_mat(k,j))*(unfolded_events_num_bins->GetBinContent(j)-h_muons_gen->GetBinContent(j));
       }
   }*/
    
  /*double dof = inv_unfolded_cov_mat.GetNrows();
  double p_value = TMath::Prob(chi_square, dof);
  std::cout<<"chi_square is: "<<chi_square<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  leg->AddEntry((TObject*)0, TString::Format("chi_square / d.o.f = %g / %g ", chi_square, dof), "");*/
  }

  unfolded_events_vals->SetStats(0);
  unfolded_events_vals->SetLineWidth(2);
  unfolded_events_vals->SetLineColor(1);
  unfolded_events_vals->Draw("same");
  unfolded_events_vals->SetDirectory(0);
  leg->AddEntry(unfolded_events,"Unfolded signal","f");
  leg->Draw();
  c2->SaveAs("overlay_chi_sqr_CC0pi.pdf");  
}


int main(int argc, char* argv[]) {
   chi_square_all_gens(argv[1]);
   return 0;
}
