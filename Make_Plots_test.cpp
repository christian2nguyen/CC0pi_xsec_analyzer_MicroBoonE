////////////////////////////////////
// Playing this Slice Plots
///////////////////////////////////


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
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "includes/PlotUtils.hh"
#include "includes/HistUtils_cc0pi.hh"
#include "includes/SliceBinning.hh"
#include "SliceHistogram.hh"
#include "TLatex.h"

//using NFT = NtupleFileType;

#define USE_FAKE_DATA ""

/// Defined Global EventCateogry 
  


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
void DrawStackMCandData_EventCategory(
TFile *TFile_MC,
TFile *TFile_Data,
TFile *TFile_Data_BeamOff,
std::string InputHisBaseName,
char *Yaxis_title,
char *Xaxis_title,
bool DoBinWidthNorm, 
char *Title,
double YMax,
std::string pdfTitle, 
TCanvas *Canvas,
std::vector<NamedCategory<EventCategory>> inputEventCateory_vector,
float POT_scaler_MC,
float BG_Trigger_scaler
);

void Draw_MC_TrackBTDProbability(
    TH1* h_Muon_input, 
    TH1* h_Proton_input,
    TH1* h_Pion_input,
    TH1* h_Else_input,
    TLegend &leg,
    char *Yaxis_title,
    char *Xaxis_title,
    bool DoBinWidthNorm, 
    char *Title,
    double YMax,
    TCanvas *Canvas);

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
void MakePlots(){

///////////////////////////////////////////////////
// Getting the POT Only looking at run 1 for now need to Change this Method  
///////////////////////////////////////////////
   std::cout<<"Starting MakePlots"<< std::endl;
   
   
   TCanvas *can = new TCanvas("can");
    char text_title_pdf1[2024];
    char text_title_pdf2[2024];
    char text_title_pdf3[2024];
    char text_title_pdf4[2024];

sprintf(text_title_pdf1, "Make_Plots_result_testingPlottong_usingSteven_NTUPLES.pdf(","" );
can -> Print(text_title_pdf1);
sprintf(text_title_pdf2, "Make_Plots_result_testingPlottong_usingSteven_NTUPLES.pdf","" );
sprintf(text_title_pdf3, "Make_Plots_result_testingPlottong_usingSteven_NTUPLES.pdf)","" );
sprintf(text_title_pdf4, "Make_Plots_result_testingPlottong_usingSteven_NTUPLES","" );
 std::string text_title_pdf4_string(text_title_pdf4);
   
   auto& fpm = FilePropertiesManager::Instance();
    const auto& EventCategory_tool = EventCategoryInterpreter::Instance();
    std::vector<NamedCategory<EventCategory>> EventCategory_NamedCategory_vector=EventCategory_tool.EventSelectionGroup_categories_;

   
   std::string InputPath_Data = "/uboone/data/users/cnguyen/CC0Pi_Selection/CCNP_0Pion_steven_test/";
std::string Data_file1 =  "DATA_Selected_stv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root";
//std::string Data_file2 = "DATA_Selected_rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_G1_beam_good_reco2_1e19.root";
//std::string Data_file_bkg = "DATA_Selected_rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";
std::string Data_file_bkg = "DATA_Selected_stv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";

std::string InputPath_MC = "/uboone/data/users/cnguyen/CC0Pi_Selection/CCNP_0Pion_steven_test/";
std::string MC_file1 =  "MC_Selected_stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
std::string MC_file2 =  "MC_Selected_stv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";
std::string MC_file3 =  "MC_Selected_stv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";

std::string MC_file4 = "MC_Selected_stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
                        

std::string InputMCFile_name = InputPath_MC + MC_file4;
std::string InputDataFile_name = InputPath_Data + Data_file1;
std::string InputData_BG_File_name = InputPath_Data + Data_file_bkg;
std::cout<< " this far"<< std::endl;
   
   
  std::map< int, double > run_to_bnb_pot_map;
  std::map< int, double > run_to_bnb_trigs_map;
  std::map< int, double > run_to_ext_trigs_map;
  float total_bnb_data_pot_ = 0.;  
  
   const auto& data_norm_map = fpm.data_norm_map();
   
   for(auto value: data_norm_map){
   std::cout<<"Map values 1st =  " << value.first << " second = " << value.second.pot_ << std::endl;
   
   }

std::cout<<"Running Ntuple File Map "<< std::endl;
  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {
    int run = run_and_type_pair.first;
    if (run ==2 || run ==3) continue; 
    std::cout << "run = "<< run << std::endl;
    
    }


  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {
    int run = run_and_type_pair.first;
    
    if (run ==2 || run ==3) continue; 
    const auto& type_map = run_and_type_pair.second;

    const auto& bnb_file_set = type_map.at( NFT::kOnBNB );
    for ( const std::string& bnb_file : bnb_file_set ) {
      const auto& pot_and_trigs = data_norm_map.at( bnb_file );


      if ( !run_to_bnb_pot_map.count(run) ) {
      std::cout<< "on  (run_to_bnb_pot_map)Run: " << run << std::endl;
        run_to_bnb_pot_map[ run ] = 0.;
        run_to_bnb_trigs_map[ run ] = 0.;
      }

      run_to_bnb_pot_map.at( run ) += pot_and_trigs.pot_;
      run_to_bnb_trigs_map.at( run ) += pot_and_trigs.trigger_count_;

    } // BNB data files

    const auto& ext_file_set = type_map.at( NFT::kExtBNB );
    for ( const std::string& ext_file : ext_file_set ) {
      const auto& pot_and_trigs = data_norm_map.at( ext_file );
 
      if ( !run_to_ext_trigs_map.count(run) ) {
        run_to_ext_trigs_map[ run ] = 0.;
      }

      run_to_ext_trigs_map.at( run ) += pot_and_trigs.trigger_count_;

    } // EXT files

  } // runs

  // Now that we have the accumulated POT over all BNB data runs, sum it
  // into a single number. This will be used to normalize the detVar MC
  // samples.
  total_bnb_data_pot_ = 0.;
  for ( const auto& pair : run_to_bnb_pot_map ) {
     int run = pair.first;
    
    if (run ==2 || run ==3) continue; 
    total_bnb_data_pot_ += pair.second;
  }

std::cout <<" Total POT  = " << total_bnb_data_pot_ << std::endl;

std::string file_name_MC = "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
//stv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root
///stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
//


float file_pot = 0;
  TFile temp_mc_file( file_name_MC.c_str(), "read" );
          TParameter<float>* temp_pot = nullptr;
          temp_mc_file.GetObject( "summed_pot", temp_pot );
          if ( !temp_pot ) throw std::runtime_error(
            "Missing POT in MC file!" );
          file_pot = temp_pot->GetVal();

std::cout<<"MC POT = "<< file_pot << std::endl;
std::cout<<"Data POT = "<< total_bnb_data_pot_  << std::endl;
float POT_MC_Scale = 4.54e+19 / 1.30503e+21 ;//total_bnb_data_pot_ / file_pot;



//float POT_MC_Scale =  4.54e19 /  file_pot;
//4.54e19

std::cout<< "POT scale = "<< POT_MC_Scale << std::endl;


float bnb_trigs = run_to_bnb_trigs_map.at( 1 );
float ext_trigs = run_to_ext_trigs_map.at( 1 );
float bd_scale = bnb_trigs / ext_trigs;
std::cout << " bnb_trigs = " << bnb_trigs << std::endl;
std::cout << " ext_trigs = " << ext_trigs << std::endl;
std::cout << " bg scale = " << bd_scale << std::endl;




char inputName[1024];
sprintf(inputName, "%s",InputMCFile_name.c_str());
std::cout<<"InputFileName MC = "<< inputName << std::endl;
TFile *TFile_MC_test = new TFile(inputName);

sprintf(inputName, "%s",InputDataFile_name.c_str());
std::cout<<"InputFileName Data = "<< inputName << std::endl;
TFile *TFile_Data = new TFile(inputName);

sprintf(inputName, "%s",InputData_BG_File_name.c_str());
std::cout<<"InputFileName BG Data = "<< inputName << std::endl;
TFile *TFile_Data_BG = new TFile(inputName);
  //std::vector<TH1D*> check_vector = MakeList_TH1D(*TFile_MC_test, "name",  EventCategory_NamedCategory_vector );


//std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH1D(*TFile_MC_test, "h_p_EventCategory", EventCategory_NamedCategory_vector );

std::vector<std::string> HistNames; 
std::vector<std::string> AxisNames{"RECO P_{#mu} [GeV/c]","RECO Cos(#theta_{#mu})"," RECO p_{p} (leading Proton) [GeV/c]", "RECO Cos(#theta_{P})"};


HistNames.push_back("h_p");
HistNames.push_back("h_costheta_p");

HistNames.push_back("h_leadingProton_p");
HistNames.push_back("h_costheta_proton");

bool DoBinWidthNorm=false;
char inputTitle[1024];

for(int i = 0 ; i < HistNames.size();i++){

sprintf(inputTitle, "%s",AxisNames.at(i).c_str());
std::cout<<"Xaxis Title = "<< inputTitle<<std::endl;

DrawStackMCandData_EventCategory(
TFile_MC_test,
TFile_Data,
TFile_Data_BG,
HistNames.at(i),
"Events ",
inputTitle,
 DoBinWidthNorm, 
"test",
99,
text_title_pdf4_string, 
can,
EventCategory_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
}
//POT_MC_Scale
//for(auto Hist :  h_p_EventCategory_StackMap){
//std::cout<<"checking titles" <<Hist.second->GetTitle()<< std::endl;
//}






  can -> Print(text_title_pdf3);
  
  
  
  
}/////End Of Function 
//////////////////////





void tutorial_slice_plots() {

   std::cout<<"Starting tutorial_slice_plots"<< std::endl;
  #ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto& fpm = FilePropertiesManager::Instance();
    std::cout<<"Finished:FilePropertiesManager::Instance() "<<std::endl; 
    std::cout<<"trying to apply load_file_properties"<< std::endl;
    fpm.load_file_properties( "nuwro_file_properties_RUN1.txt" );
    std::cout<<" passed "<< std::endl;
  #endif

 

/// Taking from my area
//  uboone/data/users/cnguyen/RootFiles_tutorial/tutorial_univmake_output.root
 
 std::cout<<"Starting MCC9SystematicsCalculator"<< std::endl;

  auto* syst_ptr = new MCC9SystematicsCalculator(
    "/uboone/data/users/gardiner/tutorial_univmake_output.root",
    "systcalc.conf" );
  auto& syst = *syst_ptr;

std::cout<<"Finished  MCC9SystematicsCalculator "<< std::endl; 

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  #ifdef USE_FAKE_DATA
    // Add the EXT to the "data" when working with fake data
    reco_bnb_hist->Add( reco_ext_hist );
  #endif

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning( "tutorial_reco_slice_config.txt" );
  auto& sb = *sb_ptr;
  /////////////////////////////////
  // Starting to plot
  //////////////////////////////////

  //std::cout<< " Number of Slices to Plot : "<< sb.slices_.size() << std::endl;

  TCanvas *c1 = new TCanvas("c1");
  sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
  c1 -> Print(pdf_title);


  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx )
  {
    std::cout<<"SLICE : "<< sl_idx<<std::endl;

    TLegend* lg1_TLegend = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg1_TLegend->SetNColumns(2);
    lg1_TLegend->SetBorderSize(0);
    Double_t defaultTextSize = lg1_TLegend->GetTextSize();
    //std::cout<<"set text size "<< defaultTextSize << std::endl;
    lg1_TLegend->SetTextSize(defaultTextSize*1.05); //
    const auto& slice = sb.slices_.at( sl_idx );

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );


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

    lg1_TLegend->AddEntry(slice_bnb->hist_.get(), "Data", "pe" );
    lg1_TLegend->AddEntry(slice_mc_plus_ext->hist_.get(), "Total MC", "le" );

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
      lg1_TLegend->AddEntry(temp_slice_mc->hist_.get(),lg1_label.c_str() , "f" );
      std::string cat_col_prefix = "MC" + std::to_string( cat );

      --cat_bin_index;
    }


   lg1_TLegend->AddEntry(slice_ext->hist_.get(),"EXT BNB" , "f" );


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

    lg1_TLegend->Draw( "same" );
   //  std::cout<<"printing chi values "<< std::endl;
    TLatex* text = new TLatex;
      text->SetNDC();
      text->SetTextSize(0.03);
      text->SetTextColor(kRed);
      sprintf(textplace, "#chi^{2}/ndf = %.2f / %i",chi2_result.chi2_ , chi2_result.num_bins_ );
      text->DrawLatex(0.15, 0.85, textplace);

      sprintf(textplace, "p-value = %.4f",chi2_result.p_value_ );
      text->DrawLatex(0.15, 0.80, textplace);
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
  
  
    const std::vector< std::string > cov_mat_keys = { "xsec_unisim",
      "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_XSecShape_CCMEC", 
      "xsec_NormCCCOH",  "xsec_NormNCCOH", "xsec_RPA_CCQE",
       "xsec_ThetaDelta2NRad",  "xsec_Theta_Delta2Npi",
       "MCstats", "EXTstats", "BNBstats"
    };

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
  TLegend* lg1_TLegend = new TLegend( 0.35, 0.6, 0.85, 0.85 );
  lg1_TLegend->SetNColumns(2);
  lg1_TLegend->SetBorderSize(0);
  Double_t defaultTextSize = lg1_TLegend->GetTextSize();
  lg1_TLegend->SetTextSize(.03); //

  TH2D* category_hist = (TH2D*)category_hist_input->Clone("category_hist");
  TH1D* h_Data =(TH1D*)slice_bnb->hist_.get()->Clone("h_Data");
  lg1_TLegend->AddEntry(h_Data, "Data", "pe" );
  TH1D* h_Total_MC =(TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_Total_MC");
  lg1_TLegend->AddEntry(h_Total_MC, "Total MC", "le" );
  
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
    lg1_TLegend->AddEntry(temp_slice_mc->hist_.get(),lg1_label.c_str() , "f" );
    //std::string cat_col_prefix = "MC" + std::to_string( cat );

    --cat_bin_index;
  }

 bool makeNormWidth = true; 

  DrawStackMCandData(
    h_Data, 
    h_Total_MC,
    slice_pred_stack,
    "NEvents / Bin Width",
    "",
    makeNormWidth, 
    "title",
    ymax,
    c1);


   lg1_TLegend->AddEntry(slice_ext->hist_.get(),"EXT BNB" , "f" );
   lg1_TLegend->Draw( "same" );
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

}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
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
    TCanvas *Canvas)
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
    //double ymaxnew = std::max( h_Data->GetMaximum(),
    //h_Total_MC->GetMaximum() ) * 1.5;
    //h_Data->GetYaxis()->SetRangeUser( 0., ymaxnew );
  }

  h_Data->SetLineColor( kBlack );
  h_Data->SetLineWidth( 3 );
  h_Data->SetMarkerStyle( kFullCircle );
  h_Data->SetMarkerSize( 0.8 );
  h_Data->SetStats( false );
  double ymax = std::max( h_Data->GetMaximum(), h_Total_MC->GetMaximum())* 1.6;
  if(YMax==99){YMax = ymax; }
 // slice_mc_plus_ext->hist_->GetMaximum() ) * 1.4;
  //if(DoBinWidt//hNorm==false) 
  h_Data->GetYaxis()->SetRangeUser( 0., YMax );
  h_Data->GetYaxis()->SetTitle( Yaxis_title );
  h_Data->GetXaxis()->SetTitle( Xaxis_title );
  h_Data->GetYaxis()->SetLabelSize(.024);
  h_Data->GetXaxis()->SetLabelSize(.025);
  h_Data->GetXaxis()->SetTitleSize(0.035);
  h_Data->GetXaxis()->CenterTitle();
  h_Data->SetTitleOffset (1.01,"Y");
  h_Data->SetTitleOffset (.9,"X");
  h_Data->Draw( "e" );
  Stack->Draw( "hist same" );
  h_Total_MC->SetLineWidth( 3 );
  h_Total_MC->Draw( "same hist e" );
  h_Data->Draw( "same e" );

}
*/
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////

void DrawStackMCandData_EventCategory(
TFile *TFile_MC,
TFile *TFile_Data,
TFile *TFile_Data_BeamOff,
std::string InputHisBaseName,
char *Yaxis_title,
char *Xaxis_title,
bool DoBinWidthNorm, 
char *Title,
double YMax,
std::string pdfTitle, 
TCanvas *Canvas,
std::vector<NamedCategory<EventCategory>> inputEventCateory_vector,
float POT_scaler_MC,
float BG_Trigger_scaler)
{
  std::string BaseStackName = InputHisBaseName+ "_EventCategory";
  std::cout << "Inside::DrawStackMCandData_EventCategory"<< std::endl;
  std::cout << "Inside:input X axis title"<< Xaxis_title<<std::endl;
  TLegend* lg1_TLegend = new TLegend( 0.35, 0.6, 0.85, 0.85 );
  lg1_TLegend->SetNColumns(2);
  lg1_TLegend->SetBorderSize(0);
  Double_t defaultTextSize = lg1_TLegend->GetTextSize();
  lg1_TLegend->SetTextSize(.03); //
  
  const auto& DrawStackMCandData_EventCategory_tool = EventCategoryInterpreter::Instance();
  std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH1D(*TFile_MC, BaseStackName, inputEventCateory_vector );
  // Set Hist to Stevens Style 

  
  TH1D* h_Data_BeamOn =  GetTH1DHist(*TFile_Data, InputHisBaseName );
  TH1D* h_Data_BeamOFF =  GetTH1DHist(*TFile_Data_BeamOff, InputHisBaseName );
  TH1D* h_MC_Total =  GetTH1DHist(*TFile_MC, InputHisBaseName );
  int Bins_BeamONData = h_Data_BeamOn->GetNbinsX();
  int Bins_BeamOFFData = h_Data_BeamOFF->GetNbinsX();
  int Bins_MC = h_MC_Total->GetNbinsX();
    
    h_MC_Total->Scale(POT_scaler_MC); 
    //h_MC_Total->Scale(BG_Trigger_scaler); 
    h_Data_BeamOFF->Scale(BG_Trigger_scaler); 
   
   float Amount_BeamOn = h_Data_BeamOn->Integral();
   float Amount_BeamOff = h_Data_BeamOFF->Integral();
   float Amount_MC_Total =  h_MC_Total->Integral();
   float area_scale = Amount_BeamOn/Amount_MC_Total;
   //h_MC_Total->Scale(area_scale); 
   //ScaleHistogramsInStack(THStack* stack, double scaleFactor, char *option )
   
     for(auto &HistInMap : h_p_EventCategory_StackMap){
     HistInMap.second->Scale(POT_scaler_MC); 
    // HistInMap.second->Scale(area_scale ); 
     }
   
   
   
     std::cout<< "NBins (BeamONData)= "<< Bins_BeamONData<< std::endl;
     std::cout<< "NBins (BeamOFFData)= "<< Bins_BeamOFFData<< std::endl;
     std::cout<< "NBins (MC)= "<< Bins_MC<< std::endl;
   

   std::cout<< "Amount_BeamOn Data = "<< Amount_BeamOn<< std::endl;
   std::cout<< "Amount_BeamOff Data = "<< Amount_BeamOff<< std::endl;
   std::cout<< "Amount_BeamOn MC = "<< Amount_MC_Total<< std::endl;
   std::cout<< "Ratio Data to MC BEAM on "<< Amount_BeamOn/Amount_MC_Total<< std::endl;
  /////////////////////////////////////////////
  // POT Scale hist , all but Data Beam Off
  ////////////////////////////////////////////

  //h_MC_Total->Add(h_Data_BeamOFF);

  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_TLegend->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_TLegend->AddEntry(h_MC_Total, "Total MC", "le" );
  EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB
  h_MC_Total->Add( h_Data_BeamOFF );

  for(auto &HistInMap : h_p_EventCategory_StackMap){
    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
    lg1_TLegend->AddEntry(HistInMap.second, HistInMap.second->GetTitle() , "f" );
  }
  
  
   h_MC_Total->GetXaxis()->SetTitle(Xaxis_title);
   h_Data_BeamOn->GetXaxis()->SetTitle( Xaxis_title );
   h_Data_BeamOFF->GetXaxis()->SetTitle( Xaxis_title );
  
  DrawStackMCandData(
      h_Data_BeamOn, 
      h_MC_Total,
      EventCategory_stack,
      Yaxis_title,
      Xaxis_title,
     DoBinWidthNorm, 
      Title,
      YMax,
      Canvas);
  
  
   TLatex* text = new TLatex;
      text->SetNDC();
      text->SetTextSize(0.03);
      text->SetTextColor(kRed);
        
      float NDataEvents =  h_Data_BeamOn->Integral();
      sprintf(textplace, "N Data Events = %.4f",NDataEvents );
      text->DrawLatex(0.1, 0.75, textplace);
      
      float NMCEvents =  h_MC_Total->Integral();
      sprintf(textplace, "N MC Events = %.4f",NMCEvents );
       text->DrawLatex(0.1, 0.70, textplace);
 

   lg1_TLegend->AddEntry(h_Data_BeamOFF,"EXT BNB" , "f" );
   lg1_TLegend->Draw( "same" );

 sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  Canvas -> Print(pdf_title);

}
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void Draw_MC_TrackBTDProbability(
    TH1* h_Muon_input, 
    TH1* h_Proton_input,
    TH1* h_Pion_input,
    TH1* h_Else_input,
    TLegend &leg,
    char *Yaxis_title,
    char *Xaxis_title,
    bool DoBinWidthNorm, 
    char *Title,
    double YMax,
    TCanvas *Canvas)
{

 TH1D* h_Muon =(TH1D*)h_Muon_input->Clone("h_Muon");
 TH1D* h_Proton =(TH1D*)h_Proton_input->Clone("h_Proton");
 TH1D* h_Pion =(TH1D*)h_Pion_input->Clone("h_Pion");
 TH1D* h_Else =(TH1D*)h_Else_input->Clone("h_Else");
 TH1D* h_Combined =(TH1D*)h_Muon_input->Clone("h_Combined");

  h_Combined ->Add(h_Proton);   
  h_Combined ->Add(h_Pion);
  h_Combined ->Add(h_Else);

  if(DoBinWidthNorm==true){
    h_Muon->Scale(1.0,"width");
    h_Proton->Scale(1.0,"width");
    h_Pion->Scale(1.0,"width");
    h_Else->Scale(1.0,"width");
    h_Combined->Scale(1.0,"width");
    
  }
  
  h_Muon->SetLineColor( 2 );
  h_Proton->SetLineColor( 4 );
  h_Pion->SetLineColor( 8 );
  h_Else->SetLineColor( 6 );
  h_Combined->SetLineColor(1 );
   leg.AddEntry(h_Muon, "Muon", "l");
   leg.AddEntry(h_Proton, "Proton", "l");
   leg.AddEntry(h_Pion, "Pion", "l");
   leg.AddEntry(h_Else, "Else", "l");
   leg.AddEntry(h_Combined, "Combined", "l");
   
    double ymax =h_Combined->GetMaximum()* 1.;
  
  if(YMax==99){YMax = ymax; }
 // slice_mc_plus_ext->hist_->GetMaximum() ) * 1.4;
  //if(DoBinWidt//hNorm==false) 
  h_Combined->GetYaxis()->SetRangeUser( 0., YMax );
  h_Combined->GetYaxis()->SetTitle( Yaxis_title );
  h_Combined->GetXaxis()->SetTitle( Xaxis_title );
  h_Combined->GetYaxis()->SetLabelSize(.024);
  h_Combined->SetTitleOffset (1.01,"Y");
  h_Combined->Draw( "Hist" );
  h_Muon->Draw( "hist same" );
  h_Proton->Draw( "same hist" );
  h_Else->Draw( "same hist" );

}// End of Funtion
///////////////////////////

void DrawProbabilityMC(
TFile *TFile_MC,
std::string InputHisBaseName,
std::vector<std::string> InputProName,
char *Yaxis_title,
char *Xaxis_title,
bool DoBinWidthNorm, 
char *Title,
double YMax,
std::string pdfTitle, 
TCanvas *Canvas,
float POT_scaler_MC,
float BG_Trigger_scaler)
{


  TLegend* legend = new TLegend( 0.35, 0.6, 0.85, 0.85 );
  legend->SetNColumns(2);
  legend->SetBorderSize(0);
  Double_t defaultTextSize = legend->GetTextSize();
  legend->SetTextSize(.03); //
  std::string HistName; 
  HistName = InputHisBaseName + InputProName.at(0);
  TH1D* h_Else =  GetTH1DHist(*TFile_MC ,HistName );
  HistName = InputHisBaseName + InputProName.at(1);
  TH1D* h_Muon =  GetTH1DHist(*TFile_MC, HistName );
  HistName = InputHisBaseName + InputProName.at(2);
  TH1D* h_Pion =  GetTH1DHist(*TFile_MC ,HistName );
  HistName = InputHisBaseName + InputProName.at(3);
  TH1D* h_Proton =  GetTH1DHist(*TFile_MC ,HistName );
  /////////////////////////////////////////////
  // POT Scale hist , all but Data Beam Off
  ////////////////////////////////////////////


   Draw_MC_TrackBTDProbability(
    h_Muon, 
    h_Proton,
    h_Pion,
    h_Else,
    *legend,
    "NEvents",
    "Probability Output",
    false, 
    Title,
    YMax,
    Canvas);

   legend->Draw( "same" );

 sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  Canvas -> Print(pdf_title);

}// End of Funtion
///////////////////////////



int main() {


std::cout<<" Testing Slice Plots "<< std::endl;

  //tutorial_slice_plots();
  MakePlots();
  return 0;
}
