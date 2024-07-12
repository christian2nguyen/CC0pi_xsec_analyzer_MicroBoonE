////////////////////////////////////
// Playing this Slice Plots
///////////////////////////////////
/*

// Standard library includes
#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TObjArray.h"


// STV analysis includes
#include "../FilePropertiesManager.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../includes/PlotUtils.hh"
#include "../includes/SliceBinning.hh"
#include "../SliceHistogram.hh"
#include "TLatex.h"
*/

#include "Make_Plots.hh"

using NFT = NtupleFileType;

#define USE_FAKE_DATA ""


// I want to print figures to PDF

char pdf_title[1024];
char textplace[1024];
std::string Pdf_name = "Tutorial_slice_Figures";

void DrawStack(
  TH2D* category_hist_input,
  TH1D* Data,
  TH1D* Total_MC,
  TH1D* extBNB_input,
  SliceHistogram* slice_bnb,
  SliceHistogram* slice_ext,
  SliceHistogram* slice_mc_plus_ext,
  const Slice& slice,
  std::string pdfTitle,
  std::string TotalTitle,
  TCanvas *c1,double ymax
);

/*
void DrawStackMCandData(
    TH1* h_data_input, 
    TH1* h_MC_Total_input,
    TH1* BG_beamOFF_input,
    THStack* Stack_input,
    char *Yaxis_title,
    char *Xaxis_title,
    bool DoBinWidthNorm, 
    char *Title,
    double YMax,
    TCanvas *Canvas
);
*/
//void ScaleHistogramsInStack(THStack* stack, double scaleFactor, char *option);

//void AddHistogramsToLegend(THStack* stack, TLegend* legend);

/*void AddPlotLabel(
        const char* label,
        const double x,
        const double y,
        const double size = 0.05,
        const int color = 1,
        const int font = 62,
        const int align = 22,
        const double angle = 0
      );

void AddHistoTitle(
  const char* title,
  double titleSize,
  int titleFont
);
*/


void tutorial_slice_plots() {

std::cout<<"Starting tutorial_slice_plots"<< std::endl;
  #ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto& fpm = FilePropertiesManager::Instance();
    std::cout<<"Finished:FilePropertiesManager::Instance() "<<std::endl; 
    std::cout<<"trying to apply load_file_properties"<< std::endl;
    fpm.load_file_properties( "nuwro_file_properties_RUN1.txt" ); /*nuwro_file_properties_run1.txt*/
    std::cout<<" passed "<< std::endl;
  #endif

 

/// Taking from my area
//  uboone/data/users/cnguyen/RootFiles_tutorial/tutorial_univmake_output.root
 ///uboone/data/users/cnguyen/CC0Pi_Selection/RUN1_TEST/univmake_RUN1_output.root
 std::cout<<"Starting MCC9SystematicsCalculator"<< std::endl;

  auto* syst_ptr = new MCC9SystematicsCalculator(
    "/uboone/data/users/gardiner/tutorial_univmake_output.root",
    "systcalc.conf" );
  auto& syst = *syst_ptr;
//
std::cout<<"Finished  MCC9SystematicsCalculator "<< std::endl; 

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate thes
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

std::cout<<"Finished getting hiss"<< std::endl;

  #ifdef USE_FAKE_DATA
    // Add the EXT to the "data" when working with fake data
    reco_bnb_hist->Add( reco_ext_hist );
  #endif

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();
std::cout<<"Finished category_hist"<< std::endl;
  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );
std::cout<<"Finished reco_mc_plus_ext_hist"<< std::endl;
  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;
  
  
std::cout<<"Finished matrix_map"<< std::endl;
  auto* sb_ptr = new SliceBinning( "tutorial_reco_slice_config_org.txt" ); /**/
  std::cout<<"Finished sb_ptr"<< std::endl;
  auto& sb = *sb_ptr;
  
  /////////////////////////////////
  // Starting to plot
  //////////////////////////////////

  std::cout<< " Number of Slices to Plot : "<< sb.slices_.size() << std::endl;

  TCanvas *c1 = new TCanvas("c1");
  sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
  c1 -> Print(pdf_title);


  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx )
  {

    std::cout<<"SLICE : "<< sl_idx<<std::endl;

    TLegend* lg1_stacked = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg1_stacked->SetNColumns(2);
    lg1_stacked->SetBorderSize(0);
    Double_t defaultTextSize = lg1_stacked->GetTextSize();
    //std::cout<<"set text size "<< defaultTextSize << std::endl;
    lg1_stacked->SetTextSize(defaultTextSize*1.05); //
    const auto& slice = sb.slices_.at( sl_idx );

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    

     std::cout<<"EXTstats"<< std::endl;
    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );
    std::cout<<"total"<< std::endl;
    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    std::cout<<"BNBstats"<< std::endl;
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
      << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
      << " p-value = " << chi2_result.p_value_ << '\n';


    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    const auto& eci = EventCategoryInterpreter::Instance();
    eci.set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB


    const auto& cat_map = eci.label_map();

    // Go in reverse so that signal ends up on top. Note that this index is
    // one-based to match the ROOT histograms
    int cat_bin_index = cat_map.size();

    lg1_stacked->AddEntry(slice_bnb->hist_.get(), "Data", "pe" );
    lg1_stacked->AddEntry(slice_mc_plus_ext->hist_.get(), "Total MC", "le" );

   //  std::cout<<"stack Loop removed r crbegin() crend() "<< std::endl;
    for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
    {
      EventCategory cat = iter->first;
      TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
        cat_bin_index, cat_bin_index );
      temp_mc_hist->SetDirectory( nullptr );

      SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
        *temp_mc_hist, slice  );

      eci.set_mc_histogram_style( cat, temp_slice_mc->hist_.get() );

      slice_pred_stack->Add( temp_slice_mc->hist_.get() );

      std::string lg1_label = eci.label(cat);
      lg1_stacked->AddEntry(temp_slice_mc->hist_.get(),lg1_label.c_str() , "f" );
      std::string cat_col_prefix = "MC" + std::to_string( cat );

      --cat_bin_index;
    }
    
   lg1_stacked->AddEntry(slice_ext->hist_.get(),"EXT BNB" , "f" );


    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 3 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 0.8 );
    slice_bnb->hist_->SetStats( false );
    double ymax = std::max( slice_bnb->hist_->GetMaximum(),
      slice_mc_plus_ext->hist_->GetMaximum() ) * 1.5;
    slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );
    slice_bnb->hist_->GetYaxis()->SetTitle( "NEvent" );
    slice_bnb->hist_->SetTitleOffset (.96,"Y");

    slice_bnb->hist_->Draw( "e" );

    slice_pred_stack->Draw( "hist same" );

    slice_mc_plus_ext->hist_->SetLineWidth( 3 );
    slice_mc_plus_ext->hist_->Draw( "same hist e" );

    slice_bnb->hist_->Draw( "same e" );

    lg1_stacked->Draw( "same" );
   //  std::cout<<"printing chi values "<< std::endl;
    TLatex* text = new TLatex;
      text->SetNDC();
      text->SetTextSize(0.03);
      text->SetTextColor(kRed);
      sprintf(textplace, "#chi^{2}/ndf = %.2f / %i",chi2_result.chi2_ , chi2_result.num_bins_ );
      text->DrawLatex(0.15, 0.85, textplace);

      sprintf(textplace, "p-value = %.4f",chi2_result.p_value_ );
      text->DrawLatex(0.15, 0.80, textplace);
      
      float NDataEvents =  slice_bnb->hist_.get()->Integral();
      sprintf(textplace, "NDataEvents = %.4f",NDataEvents );
      text->DrawLatex(0.15, 0.75, textplace);
      
      float NMCEvents =  slice_mc_plus_ext->hist_.get()->Integral();
      sprintf(textplace, "NDataEvents = %.4f",NMCEvents );
       text->DrawLatex(0.15, 0.70, textplace);
      sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());
      c1 -> Print(pdf_title);
    //  std::cout<<"FINISHED printing chi values "<< std::endl;

    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    /////////////////////////////////////
    // testing Drawing Function
    ////////////////////////////////////
    TH1D* h_bnb_input = (TH1D*)reco_bnb_hist->Clone("slice_bnb_input");
    TH1D* h_mc_plus_ext_input = (TH1D*)reco_mc_plus_ext_hist->Clone("slice_mc_plus_ext_input");
    TH1D* h_ext_input = (TH1D*)reco_ext_hist->Clone("slice_ext_input");

    DrawStack(
      category_hist,
      h_bnb_input,
      h_mc_plus_ext_input,
      h_ext_input,
       slice_bnb,
       slice_ext,
       slice_mc_plus_ext,
      slice,
      Pdf_name,
      "Test",
      c1,ymax);


    ////////////////////////////////////////
    ///// Finished with First Plot
    ////////////////////////////////////////

   //if(sl_idx==3)continue; 

    //std::cout<<"Starting Making error Plot  "<< std::endl;

    TH1* slice_hist = dynamic_cast< TH1* >(
      slice.hist_->Clone("slice_hist") );

      slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Show fractional uncertainties computed using these covariance matrices
    // in the ROOT plot. All configured fractional uncertainties will be
    // included in the output pgfplots file regardless of whether they appear
    // in this vector.

    //const std::vector< std::string > cov_mat_keys = { "total",
    //  "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets",
    //  "MCstats", "EXTstats", "BNBstats", "xsec_NormCCCOH"
    //};

 //const std::vector< std::string > cov_mat_keys = { "xsec_total",
 //     "xsec_multi", "xsec_unisim", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC",
 //      "POT","MCstats", "EXTstats", "BNBstats"
 //   };
  
    
        const std::vector< std::string > cov_mat_keys = { "total",
      "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets",
      "MCstats", "EXTstats", "BNBstats"
    };
    
  
  //  const std::vector< std::string > cov_mat_keys = { "xsec_unisim",
  //    "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_XSecShape_CCMEC", 
  //    "xsec_NormCCCOH",  "xsec_NormNCCOH", "xsec_RPA_CCQE",
  //     "xsec_ThetaDelta2NRad",  "xsec_Theta_Delta2Npi",
  //     "MCstats", "EXTstats", "BNBstats"
  //  };

   //  std::cout<<"Looping over the various systematic uncertainties "<< std::endl;
    int loop_count=0;
    // Loop over the various systematic uncertainties
    int color = 0;
    for ( const auto& pair : matrix_map )
      {
      //std::cout<< "loop_count = " << loop_count<< std::endl;
      loop_count++;
      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;

      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      //std::cout<<" this loop starts for ( const auto& bin_pair : slice.bin_map_ )"<< std::endl;
      for ( const auto& bin_pair : slice.bin_map_ )
      {
        int global_bin_idx = bin_pair.first;
        double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
        double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }
      //std::cout<<" this loop Ends for ( const auto& bin_pair : slice.bin_map_ )"<< std::endl;
      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend ) continue;

      frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 10;

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );

    }

     //TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.2, 0.7, 0.7, 0.88 );
    lg2->SetNColumns(2);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(.025); //
    //Double_t defaultTextSize_lg2 = lg2->GetTextSize();
    //std::cout<<"set text size "<< defaultTextSize_lg2 << std::endl;
    //lg2->SetTextSize(defaultTextSize_lg2*1.05);

   //  std::cout<<" trying to call frac_uncertainty_hists.at( Total Error ) " << std::endl;
    //auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
    auto* total_frac_err_hist = frac_uncertainty_hists.at( "xsec_unisim" );
    //  std::cout<<" Finshed " << std::endl;
    frac_uncertainty_hists.at( "BNBstats" )->SetLineStyle(2);
    frac_uncertainty_hists.at( "EXTstats" )->SetLineStyle(2);
    frac_uncertainty_hists.at( "MCstats" )->SetLineStyle(2);
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.25 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->GetYaxis()->SetTitle( "Fractional Error" );
    total_frac_err_hist->Draw( "hist" );

    lg2->AddEntry( total_frac_err_hist, "xsec_unisim", "l" );

   //  std::cout<<"starting error leg names"<<std::endl;
    //////////////////////////////////  
    /// Drawing the Sysmatics 
    //////////////////////////////////  
    for ( auto& pair : frac_uncertainty_hists )
    {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "xsec_unisim" ) continue;

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );

      std::cout << name << " frac err in bin #1 = "
      << hist->GetBinContent( 1 )*100. << "%\n";
    }
    //std::cout<<"END error leg names"<<std::endl;
    lg2->Draw( "same" );

    std::cout << "Total frac error in bin #1 = "
      << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";
      //sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());
      c1 -> Print(pdf_title);


  } 
  //////////////////////////////////
  // End of slices Loop 
  //////////////////////////////////


  sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
  c1 -> Print(pdf_title);


}

void DrawStack(
  TH2D* category_hist_input,
  TH1D* Data_input,
  TH1D* Total_MC_input,
  TH1D* extBNB_input,
  SliceHistogram* slice_bnb,
  SliceHistogram* slice_ext,
  SliceHistogram* slice_mc_plus_ext,
  const Slice& slice,
  std::string pdfTitle,
  std::string TotalTitle,
  TCanvas *c1,
  double ymax)
{
  const auto& EventInterp = EventCategoryInterpreter::Instance();
  TLegend* lg1_stacked = new TLegend( 0.35, 0.6, 0.85, 0.85 );
  lg1_stacked->SetNColumns(2);
  lg1_stacked->SetBorderSize(0);
  Double_t defaultTextSize = lg1_stacked->GetTextSize();
  lg1_stacked->SetTextSize(.03); //

  TH2D* category_hist = (TH2D*)category_hist_input->Clone("category_hist");
  TH1D* h_Data =(TH1D*)slice_bnb->hist_.get()->Clone("h_Data");
  lg1_stacked->AddEntry(h_Data, "Data", "pe" );
  TH1D* h_Total_MC =(TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_Total_MC");
  lg1_stacked->AddEntry(h_Total_MC, "Total MC", "le" );
  
  TH1D* h_extBNB =(TH1D*)extBNB_input->Clone("h_extBNB");
  EventInterp.set_ext_histogram_style( h_extBNB );
  
  EventInterp.set_ext_histogram_style( slice_ext->hist_.get() );
  auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );

  THStack* slice_pred_stack = new THStack( "mc+ext", "" );
  slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

  const auto& cat_map = EventInterp.label_map();
  int cat_bin_index = cat_map.size();

  for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
  {
    EventCategory cat = iter->first;
    TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
      cat_bin_index, cat_bin_index );
    temp_mc_hist->SetDirectory( nullptr );
    
    SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
      *temp_mc_hist, slice  );

    EventInterp.set_mc_histogram_style( cat, temp_slice_mc->hist_.get() );

    slice_pred_stack->Add( temp_slice_mc->hist_.get() );

    std::string lg1_label = EventInterp.label(cat);
    lg1_stacked->AddEntry(temp_slice_mc->hist_.get(),lg1_label.c_str() , "f" );
    //std::string cat_col_prefix = "MC" + std::to_string( cat );

    --cat_bin_index;
  }

 bool makeNormWidth = true; 

  DrawStackMCandData(
    h_Data, 
    h_Total_MC,
    h_extBNB,
    slice_pred_stack,
    "NEvents / Bin Width",
    "",
    makeNormWidth, 
    "title",
    ymax,
    c1);


   lg1_stacked->AddEntry(slice_ext->hist_.get(),"EXT BNB" , "f" );
   lg1_stacked->Draw( "same" );
  //  std::cout<<"printing chi values "<< std::endl;
  TLatex* text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextColor(kRed);
    sprintf(textplace, "#chi^{2}/ndf = %.2f / %i",chi2_result.chi2_ , chi2_result.num_bins_ );
    text->DrawLatex(0.15, 0.85, textplace);

    sprintf(textplace, "p-value = %.4f",chi2_result.p_value_ );
    text->DrawLatex(0.15, 0.80, textplace);
    //AddHistoTitle(TotalTitle.c_str(), .035);
  ////////////////////////////////////////
  ///// Finished with First Plot
  ////////////////////////////////////////

  sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  c1 -> Print(pdf_title);

}
//////////////////////////////////////////////////////////////////////////////
/// 
//////////////////////////////////////////////////////////////////////////////
/*
void DrawStackMCandData(
    TH1* h_data_input, 
    TH1* h_MC_Total_input,
    TH1* BG_beamOFF_input,
    THStack* Stack_input,
    char *Yaxis_title,
    char *Xaxis_title,
    bool DoBinWidthNorm, 
    char *Title,
    double YMax,
    TCanvas *Canvas
)
{

 TH1D* h_Data =(TH1D*)h_data_input->Clone("h_Data");
 TH1D* h_Total_MC =(TH1D*)h_MC_Total_input->Clone("h_Total_MC");
 TH1D* h_extBNB =(TH1D*)BG_beamOFF_input->Clone("h_extBNB");
 THStack* Stack = (THStack*)Stack_input->Clone("Stack");

  if(DoBinWidthNorm==true){
    h_Data->Scale(1.0,"width");
    h_Total_MC->Scale(1.0,"width");
    h_extBNB->Scale(1.0,"width");
    ScaleHistogramsInStack(Stack, 1.0, "width" );
    double ymaxnew = std::max( h_Data->GetMaximum(),
    h_Total_MC->GetMaximum() ) * 1.5;
    h_Data->GetYaxis()->SetRangeUser( 0., ymaxnew );
  }

  h_Data->SetLineColor( kBlack );
  h_Data->SetLineWidth( 3 );
  h_Data->SetMarkerStyle( kFullCircle );
  h_Data->SetMarkerSize( 0.8 );
  h_Data->SetStats( false );
  //double ymax = std::max( slice_bnb->hist_->GetMaximum(),
    //slice_mc_plus_ext->hist_->GetMaximum() ) * 1.4;
  if(DoBinWidthNorm==false) h_Data->GetYaxis()->SetRangeUser( 0., YMax );
  h_Data->GetYaxis()->SetTitle( Yaxis_title );
  //h_Data->GetXaxis()->SetTitle( Xaxis_title );
  h_Data->GetYaxis()->SetLabelSize(.024);
  h_Data->SetTitleOffset (1.01,"Y");
  h_Data->Draw( "e" );
  Stack->Draw( "hist same" );
  h_Total_MC->SetLineWidth( 3 );
  h_Total_MC->Draw( "same hist e" );
  h_Data->Draw( "same e" );
  //lg1_stacked->Draw( "same" );
}
*/


//void AddHistogramsToLegend(THStack* stack, TLegend* legend) {
//      
//      //TObjArray* histArray = stack->GetHists();
//      int nVars =stack->GetNhists ()
//
//    for (int i = 0; i < nVars; i++) {
//        TH1* hist = static_cast<TH1*>(stack->At(i));
//			//std::string name = hist->GetTitle();
//        //TString name = histArray->At(i)->GetName();
//        legend->AddEntry(hist, name, "l");
//    }
//}



//void AddPlotLabel(
//        const char* label,
//        const double x,
//        const double y,
//        const double size /*= 0.05*/,
//        const int color /*= 1*/,
//        const int font /*= 62*/,
//        const int align /*= 22*/,
//        const double angle /*= 0*/
//        )
/*
{
    TLatex *latex = new TLatex( x, y, label );
    AddToTmp( latex );

    latex->SetNDC();
    latex->SetTextSize(size);
    latex->SetTextColor(color);
    latex->SetTextFont(font);
    latex->SetTextAlign(align);
    latex->SetTextAngle(angle);
    latex->Draw();
}

void AddHistoTitle(
        const char* title,
        double titleSize,
        int titleFont
        )
{
    AddPlotLabel(title, 0.5, 1-gStyle->GetPadTopMargin()-extra_top_margin, titleSize, 1, titleFont, 21);
}

*/

/*void ScaleHistogramsInStack(THStack* stack, double scaleFactor, char *option ) {
    TList* list = stack->GetHists();
    
    for (int i = 0; i < list->GetSize(); i++) {
        TH1* hist = static_cast<TH1*>(list->At(i));
        if (hist) {
            hist->Scale(scaleFactor, option);
        }
    }
}
*/


int main() {


std::cout<<" Testing Slice Plots "<< std::endl;

  tutorial_slice_plots();
  return 0;
}
