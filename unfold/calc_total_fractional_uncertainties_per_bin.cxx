/*Author: P. Englezos <p.englezos@physics.rutgers.edu>
 *
 This script is used for testing if my outputs match with Steven's ones.
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
#include "EventCategory.hh"
#include "TreeUtils.hh"
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
//#include "includes/SliceBinning.hh"
//#include "includes/SliceHistogram.hh"
#include "TMatrixT.h"
//#include "includes/WienerSVDUnfolder.hh"
#include "FiducialVolume.hh"
#include <iomanip>
#include "TMatrixD.h"
#include "MatrixUtils.hh"

int main() {
  using NFT = NtupleFileType;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  constexpr int FONT_STYLE = 62;
  
 
  //auto* syst_ptr = new MCC9SystematicsCalculator("./univmake_genie_with_uncertainties_runs1-3_muon_angle.root", "systcalc.conf" );
  //auto* syst_ptr = new MCC9SystematicsCalculator("./univmake_genie_with_uncertainties_runs1-3_muon_momentum_really.root", "systcalc.conf" );
  //auto* syst_ptr = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_uncertainties_runs_1-3_4.root", "systcalc.conf" );
  //auto* syst_ptr = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_uncertainties_runs_1-3_angle.root", "systcalc.conf" );
  
  //auto* syst_ptr = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_genie_closure_runs_1-3_momentum_old_2.root", "systcalc.conf" );
  auto* syst_ptr = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_genie_closure_runs_1-3_angle_old.root", "systcalc.conf" );
  auto& syst = *syst_ptr;

  //std::vector< double > nbins = {0.1,0.275,0.375,0.525,0.7,0.9,1.1,2.};
  std::vector< double > nbins = {-1.,-0.525,-0.375,-0.225,-0.1,0.025,0.15,0.275,0.45,0.6,0.775,0.925,1.};

  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  int num_reco_bins = syst.reco_bins_.size();
  TH1D* reco_pred_hist = new TH1D( "reco_pred_hist", "; reco bin; events", num_reco_bins, 0., num_reco_bins );
  reco_pred_hist->Sumw2();


  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
  reco_pred_hist->Add( reco_ext_hist );
  reco_pred_hist->Add( syst.cv_universe().hist_reco_.get() );

  auto* fr_unc_hists = new std::map< std::string, TH1D* >();
  auto& frac_uncertainty_hists = *fr_unc_hists;

  const std::vector< std::string > total_cov_mat_keys = { /*"detVar_total",
    "flux", "reint", "xsec_total", "POT", "numTargets",*/ "MCstats"/*, "EXTstats"*/
  };

  int color = 1;
  for ( const auto& key : total_cov_mat_keys ) {

    const auto& temp_results = matrix_map.at( key );

    /*TH1D* temp_hist = new TH1D( ("myfrac_temp_hist_" + key).c_str(),
      "; reco bin; fractional uncertainty", num_reco_bins, 0., num_reco_bins );*/
    TH1D* temp_hist = new TH1D( ("myfrac_temp_hist_" + key).c_str(),
      "; reco bin; fractional uncertainty", nbins.size() - 1, nbins.data() );

    temp_hist->SetStats( false );

   for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
      double err2 = temp_results.cov_matrix_->GetBinContent( rb, rb );
      double err = std::sqrt( std::max(0., err2) );
      double cv = reco_pred_hist->GetBinContent( rb );
      if ( cv > 0. ) err /= cv;
      else err = 0.;

      temp_hist->SetBinContent( rb, err );

      frac_uncertainty_hists[ key ] = temp_hist;
    }

    if ( color <= 9 ) ++color;
    if ( color == 5 ) ++color;
    if ( color >= 10 ) color += 10;
    if (key == "MCstats") color = kOrange + 3;

    temp_hist->SetLineColor( color );
    temp_hist->SetLineWidth( 3 );
    temp_hist->GetYaxis()->SetRangeUser( 0., 1. );
  }
  
  const auto& total_cov_matrix = matrix_map.at( "total" ).cov_matrix_;
  const auto& total_cov_matrix_mcstats = matrix_map.at( "MCstats" ).cov_matrix_;
  for ( size_t a = 1u; a <= num_reco_bins; ++a ) {
    double cov = total_cov_matrix->GetBinContent( a, a );
    double cov_mcstats = total_cov_matrix_mcstats->GetBinContent( a, a );
    std::cout<<" total: "<<cov<<"  mcstats: "<<cov_mcstats<<std::endl;
    reco_pred_hist->SetBinError( a, std::sqrt(std::max(0., cov)) );
  }

  /*TH1D* total_frac_err_hist = new TH1D( "total_frac_err_hist",
    "; reconstructed bin; events", num_reco_bins, 0., num_reco_bins );*/
   TH1D* total_frac_err_hist = new TH1D( "total_frac_err_hist",
    "; cos#theta_{#mu}; events", nbins.size() - 1, nbins.data() );
  for ( size_t a = 1u; a <= num_reco_bins; ++a ) {
    double cv = reco_pred_hist->GetBinContent( a );
    double err = reco_pred_hist->GetBinError( a );
    if ( cv > 0. ) err /= cv;
    else err = 0.;
    total_frac_err_hist->SetBinContent( a, err );
  }

  TCanvas* c2 = new TCanvas;
  TLegend* lg2 = new TLegend( 0.2, 0.62, 0.49, 0.9 );
 
  total_frac_err_hist->SetStats( false );
  total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
  total_frac_err_hist->GetMaximum() * 1.05 );
  total_frac_err_hist->SetLineColor( kBlack );
  total_frac_err_hist->SetLineWidth( 3 );
 
  total_frac_err_hist->GetYaxis()->SetTitle("fractional uncertainty");
  total_frac_err_hist->Draw( "hist" );

  lg2->AddEntry( total_frac_err_hist, "Total", "l" );

  for ( auto& pair : frac_uncertainty_hists ) {
    const auto& name = pair.first;
    TH1D* hist = pair.second;


    //lg2->AddEntry( hist, name.c_str(), "l" );
    if (name == "xsec_total" ){
       lg2->AddEntry( hist, "Neutrino interaction modeling", "l" );
    }
    else if (name == "reint" ){
       lg2->AddEntry( hist, "Reinteractions", "l" );
    }
    else if (name == "numTargets" ){
       lg2->AddEntry( hist, "Number of target nuclei", "l" );
    }
    else if (name == "flux" ){
       lg2->AddEntry( hist, "Neutrino flux", "l" );
    }
    else if (name == "detVar_total" ){
       lg2->AddEntry( hist, "Detector response", "l" );
    }
    else if (name == "POT" ){
       lg2->AddEntry( hist, "Protons on target", "l" );
    }   
    else if (name == "MCstats" ){
       lg2->AddEntry( hist, "Monte Carlo statistics", "l" );
    }
    else if (name == "EXTstats" ){
       lg2->AddEntry( hist, "Beam-off statistics", "l" );
    }
    hist->Draw( "same hist" );

  }

  lg2->Draw( "same" );
  c2->SaveAs("0pi_fract_syst.pdf");

return 0;
}

