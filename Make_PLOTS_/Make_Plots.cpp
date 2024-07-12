////////////////////////////////////
// Playing this Slice Plots
///////////////////////////////////

#include "Make_Plots.hh"
//using NFT = NtupleFileType;

#define USE_FAKE_DATA ""

/// Defined Global EventCateogry 
  


// I want to print figures to PDF

char pdf_title[1024];
char textplace[1024];
std::string Pdf_name = "Tutorial_slice_Figures";
const double DATA_POT = 4.54e+19;

//const double DATA_POT = 3.07997e+20; // Newro

const double MC_POT = 1.30503e+21;

 
 



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

	void Draw2DStackMCandData_EventCategory(
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
float BG_Trigger_scaler);

void DrawStackMCandData_Containment(
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
std::vector<NamedCategory<FidVol_Category>> inputEventCateory_vector,
float POT_scaler_MC,
float BG_Trigger_scaler);

void DrawStackMCandData_Topology(
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
std::vector<NamedCategory<CCZeroPi_type>> inputEventCateory_vector,
float POT_scaler_MC,
float BG_Trigger_scaler);

void DrawStackMCandData_ParticleReduced(
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
std::vector<NamedCategory<Particle_type>> inputEventCateory_vector,
float POT_scaler_MC,
float BG_Trigger_scaler);

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


void DrawFakeData();

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

  sprintf(text_title_pdf1, "Make_Plots_result_Newcc0pi_BDT_Thresholds_WithNoProtonCut.pdf(","" );
  can -> Print(text_title_pdf1);
  sprintf(text_title_pdf2, "Make_Plots_result_Newcc0pi_BDT_Thresholds_WithNoProtonCut.pdf","" );
  sprintf(text_title_pdf3, "Make_Plots_result_Newcc0pi_BDT_Thresholds_WithNoProtonCut.pdf)","" );
  sprintf(text_title_pdf4, "Make_Plots_result_Newcc0pi_BDT_Thresholds_WithNoProtonCut","" );
  std::string text_title_pdf4_string(text_title_pdf4);
  std::string text_title_pdf2_string(text_title_pdf2);
   auto& fpm = FilePropertiesManager::Instance();
    const auto& EventCategory_tool = EventCategoryInterpreter::Instance();
    std::vector<NamedCategory<EventCategory>> EventCategory_NamedCategory_vector= EventCategory_tool.EventSelectionGroup_categories_;


    std::vector<NamedCategory<FidVol_Category>> Containment_NamedCategory_vector= EventCategory_tool.ContainedGroup_categories_;
    std::vector<NamedCategory<CCZeroPi_type>> Topology_NamedCategory_vector= EventCategory_tool.topology_categories_;

    std::vector<NamedCategory<Particle_type>> ReduceParticle_NamedCategory_vector= EventCategory_tool.ParticleGroup_reduced_categories_;

       std::vector<NamedCategory<BDT_Category>> BTDPrectionType_NamedCategory_vector= EventCategory_tool.BTD_Group_categories_;
   
   std::string InputPath_Data = "/uboone/data/users/cnguyen/CC0Pi_Selection/Data/";
//std::string Data_file1 =  "DATA_Selected_newCC0pi_addthresholds_rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root";

////// New Wro
std::string Data_file1 =  "DATA_Selected_newCC0pi_rustv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_NuWro_reco2_reco2.root";
//addthresholds_WithProtonPIDCUT_

std::string Data_file2 = "DATA_Selected_newCC0pi_rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root";
//std::string Data_file_bkg = "DATA_Selected_rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";
std::string Data_file_bkg = "DATA_Selected_newCC0pi_rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";
//addthresholds_WithProtonPIDCUT_
std::string InputPath_MC = "/uboone/data/users/cnguyen/CC0Pi_Selection/MC/";
std::string MC_file1 =  "MC_Selected_newCC0pi_rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
std::string MC_file2 =  "MC_Selected_rustv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";
std::string MC_file3 =  "MC_Selected_rustv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";
                        //MC_Selected_newCC0pi_addthresholds_rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2
std::string MC_file4 = "MC_Selected_newCC0pi_rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
                        //MC_Selected_newCC0pi_addthresholds_rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
                        //_addthresholds_WithProtonPIDCUT
//std::string MC_file4 = "MC_Selected_newCC0pi_addthresholds_rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
     
     
std::string Fake_data_file = " /uboone/data/users/cnguyen/CC0Pi_Selection/Data/DATA_Selected_newCC0pi_addthresholds_rustv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_NuWro_reco2_reco2.root";
     
                    
std::string InputMCFile_name = InputPath_MC + MC_file4;
std::string InputDataFile_name = InputPath_Data + Data_file1;
std::string InputREAL_DataFile_name = InputPath_Data + Data_file2;
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

  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {
    int run = run_and_type_pair.first;
    std::cout << "run = "<< run << std::endl;
    
    }


  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {
    int run = run_and_type_pair.first;
    
    if (run ==2 || run ==3) continue; 
    const auto& type_map = run_and_type_pair.second;

    const auto& bnb_file_set = type_map.at( NFT::kOnBNB );
    for ( const std::string& bnb_file : bnb_file_set ) {
      const auto& pot_and_trigs = data_norm_map.at( bnb_file );
     if (run ==2 || run ==3) continue;
     
      if ( !run_to_bnb_pot_map.count(run) ) {
        run_to_bnb_pot_map[ run ] = 0.;
        run_to_bnb_trigs_map[ run ] = 0.;
      }

      run_to_bnb_pot_map.at( run ) += pot_and_trigs.pot_;
      run_to_bnb_trigs_map.at( run ) += pot_and_trigs.trigger_count_;

    } // BNB data files

    const auto& ext_file_set = type_map.at( NFT::kExtBNB );
    for ( const std::string& ext_file : ext_file_set ) {
    if (run ==2 || run ==3) continue;
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

std::string file_name_MC =  "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
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
//float POT_MC_Scale = total_bnb_data_pot_ / file_pot;
float POT_MC_Scale = (DATA_POT / MC_POT);
//float POT_MC_Scale =  4.54e19 /  file_pot;
//4.54e19



std::cout<< "POT scale = "<< POT_MC_Scale << std::endl;

 


float bnb_trigs = run_to_bnb_trigs_map.at( 1 );
float ext_trigs = run_to_ext_trigs_map.at( 1 );
float bd_scale2 = bnb_trigs / ext_trigs;
//float bd_scale = 36139233. / 65498807.;
//float bd_scale = 68708488. / 65498807.;
//float bd_scale =   36139233. / 68708488.;


float bd_scale = 10080350. / 65498807.;

std::cout << " bnb_trigs = " << bnb_trigs << std::endl;
std::cout << " ext_trigs = " << ext_trigs << std::endl;
std::cout << " bg scale = " << bd_scale << std::endl;
std::cout << " bg scale2 = " << bd_scale2 << std::endl;



char inputName[1024];
sprintf(inputName, "%s",InputMCFile_name.c_str());


  //std::vector<TH1D*> check_vector = MakeList_TH1D(*TFile_MC_test, "name",  EventCategory_NamedCategory_vector );


//std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH1D(*TFile_MC_test, "h_p_EventCategory", EventCategory_NamedCategory_vector );

std::vector<std::string> HistNames; 
std::vector<std::string> AxisNames{
"RECO P_{#mu} [GeV/c]",
"RECO CC0P0Pi P_{#mu} [GeV/c]",
"RECO  CCNP0Pi P_{#mu} [GeV/c]",
"RECO Proton Multiplicity [N]",
"RECO Cos(#theta_{#mu})",
"RECO CC0P0Pi Cos(#theta_{#mu})",
"RECO CCNP0Pi Cos(#theta_{#mu})",
" RECO p_{p} (leading Proton) [GeV/c]", "RECO Cos(#theta_{P})",
"RECO #theta_{opening}",
"#delta #phi_{T}" , 
"pn [GeV/c]", "#delta #alpha_{T}",
"#delta #alpha_{pTx}", 
"#delta #alpha_{pTy}",
"#delta #alpha #phi_{T}",
"Vertex X [cm]",
"Vertex Y [cm]", 
"Vertex Z [cm]"
};


HistNames.push_back("h_p");
HistNames.push_back("h_p_Correction");
//HistNames.push_back("h_p_CC0P0Pi_Only");
HistNames.push_back("h_p_CCNP0Pi_Only");
HistNames.push_back("h_Proton_mult");
HistNames.push_back("h_costheta_p");

HistNames.push_back("h_costheta_p_CC0P0Pi_Only");
HistNames.push_back("h_costheta_p_CCNP0Pi_Only");

HistNames.push_back("h_leadingProton_p");
HistNames.push_back("h_costheta_proton");

HistNames.push_back("h_openingAngle_deg");

 HistNames.push_back("h_delta_phiT");
 HistNames.push_back("h_pn");
 HistNames.push_back("h_delta_alphaT");
 HistNames.push_back("h_delta_pTx");
 HistNames.push_back("h_delta_pTy");
 HistNames.push_back("h_delta_phiT");
 HistNames.push_back("h_VertexX");
 HistNames.push_back("h_VertexY");
 HistNames.push_back("h_VertexZ");

bool DoBinWidthNorm=true;
char inputTitle[1024];

for(int i = 0 ; i < HistNames.size();i++){

sprintf(inputTitle, "%s",AxisNames.at(i).c_str());
std::cout<<"Xaxis Title = "<< inputTitle<<std::endl;

std::cout<<"InputFileName = "<< inputName << std::endl;
sprintf(inputName, "%s",InputMCFile_name.c_str());
TFile *TFile_MC_test = new TFile(inputName);
TFile *TFile_MC_test2 = new TFile(inputName);
sprintf(inputName, "%s",InputDataFile_name.c_str());
TFile *TFile_Data = new TFile(inputName);
TFile *TFile_Data2 = new TFile(inputName);
sprintf(inputName, "%s",InputData_BG_File_name.c_str());
TFile *TFile_Data_BG = new TFile(inputName);
TFile *TFile_Data_BG2 = new TFile(inputName);




DrawStackMCandData_EventCategory(
TFile_MC_test,
TFile_Data,
TFile_Data_BG,
HistNames.at(i),
"Events / [Binwidth]",
inputTitle,
 DoBinWidthNorm, 
"",
99,
text_title_pdf4_string, 
can,
EventCategory_NamedCategory_vector,
POT_MC_Scale,
bd_scale);


if(HistNames.at(i)=="h_p"|| HistNames.at(i)=="h_costheta_p"|| HistNames.at(i)=="h_p_Correction"){


DrawStackMCandData_Containment(
TFile_MC_test2,
TFile_Data2,
TFile_Data_BG2,
HistNames.at(i),
"Events / [Binwidth]",
inputTitle,
 DoBinWidthNorm, 
"test",
99,
text_title_pdf4_string, 
can,
Containment_NamedCategory_vector,
POT_MC_Scale,
bd_scale);


sprintf(inputName, "%s",InputMCFile_name.c_str());
TFile *TFile_MC_test3 = new TFile(inputName);
sprintf(inputName, "%s",InputDataFile_name.c_str());
TFile *TFile_Data3 = new TFile(inputName);
sprintf(inputName, "%s",InputData_BG_File_name.c_str());
TFile *TFile_Data_BG3 = new TFile(inputName);



DrawStackMCandData_Topology(
TFile_MC_test3,
TFile_Data3,
TFile_Data_BG3,
HistNames.at(i),
"Events / [Binwidth]",
inputTitle,
 DoBinWidthNorm, 
"test",
99,
text_title_pdf4_string, 
can,
Topology_NamedCategory_vector,
POT_MC_Scale,
bd_scale);

TFile_MC_test3->Close();
TFile_Data3->Close();
TFile_Data_BG3->Close();

delete TFile_MC_test3;
delete TFile_Data3;
delete TFile_Data_BG3;

} // End of If


TFile_MC_test->Close();
TFile_Data->Close();
TFile_Data_BG->Close();

TFile_MC_test2->Close();
TFile_Data2->Close();
TFile_Data_BG2->Close();

delete TFile_MC_test;
delete TFile_Data;
delete TFile_Data_BG;


delete TFile_MC_test2;
delete TFile_Data2;
delete TFile_Data_BG2;
   // Delete the TFile object (optional, as it will be deleted when it goes out of scope)

}// End of Loop 
/////////////////////////////////////////////////
/////////////////////////////////////////////////


sprintf(inputName, "%s",InputMCFile_name.c_str());
TFile *TFile_MC_test4 = new TFile(inputName);
sprintf(inputName, "%s",InputDataFile_name.c_str());
TFile *TFile_Data4 = new TFile(inputName);
sprintf(inputName, "%s",InputData_BG_File_name.c_str());
TFile *TFile_Data_BG4 = new TFile(inputName);

 DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_NTracks",
  "Events / [Binwidth]",
  "Track N",
   DoBinWidthNorm, 
  "Track's TRUE pdg",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);
  
   DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_Proton_mult_Tracks",
  "NEvents",
  "Proton Track Candidates",
   !DoBinWidthNorm, 
  "BDT Predicted Proton Tracks",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);
  
    
  DrawStackMCandData_EventCategory(
TFile_MC_test4,
TFile_Data4,
TFile_Data_BG4,
"h_Proton_mult_Tracks",
"Events",
"Proton Track Candidates",
 !DoBinWidthNorm, 
"BDT Predicted Proton Tracks",
99,
text_title_pdf4_string, 
can,
EventCategory_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
  
  
   DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_CCZeroP_MuonCandidate_BTDPrediction",
  "Events",
  "BDT Prediction Score",
   !DoBinWidthNorm, 
  "Muon Candidate",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);
  
  
  DrawStackMCandData_EventCategory(
TFile_MC_test4,
TFile_Data4,
TFile_Data_BG4,
"h_CCZeroP_MuonCandidate_BTDPrediction",
"Events",
"BDT Prediction Score",
 !DoBinWidthNorm, 
"Muon Candidate",
99,
text_title_pdf4_string, 
can,
EventCategory_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
  
 DrawStackMCandData_Containment(
TFile_MC_test4,
TFile_Data4,
TFile_Data_BG4,
"h_CCZeroP_MuonCandidate_BTDPrediction",
"Events / [Binwidth]",
"BDT Prediction Score",
 !DoBinWidthNorm, 
"Muon Candidate",
99,
text_title_pdf4_string, 
can,
Containment_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
 
  DrawStackMCandData_Topology(
TFile_MC_test4,
TFile_Data4,
TFile_Data_BG4,
"h_CCZeroP_MuonCandidate_BTDPrediction",
"Events / [Binwidth]",
"BDT Prediction Score",
 !DoBinWidthNorm, 
"Muon Candidate",
99,
text_title_pdf4_string, 
can,
Topology_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
  
  
     DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_CCZeroP_ProtonCandidate_BTDPrediction",
  "Events",
  "BDT Prediction Score",
   !DoBinWidthNorm, 
  "BDT Predicted Proton Candidates",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);

    DrawStackMCandData_EventCategory(
     TFile_MC_test4,
     TFile_Data4,
     TFile_Data_BG4,
     "h_CCZeroP_ProtonCandidate_BTDPrediction",
     "Events",
     "BDT Prediction Score",
      !DoBinWidthNorm, 
     "BDT Predicted Proton Candidates",
     99,
     text_title_pdf4_string, 
     can,
     EventCategory_NamedCategory_vector,
     POT_MC_Scale,
     bd_scale);
  
  
   DrawStackMCandData_Containment(
TFile_MC_test4,
TFile_Data4,
TFile_Data_BG4,
"h_CCZeroP_ProtonCandidate_BTDPrediction",
"Events / [Binwidth]",
"BDT Prediction Score",
 !DoBinWidthNorm, 
  "BDT Predicted Proton Candidates",
99,
text_title_pdf4_string, 
can,
Containment_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
 
  DrawStackMCandData_Topology(
TFile_MC_test4,
TFile_Data4,
TFile_Data_BG4,
"h_CCZeroP_ProtonCandidate_BTDPrediction",
"Events / [Binwidth]",
"BDT Prediction Score",
 !DoBinWidthNorm, 
  "BDT Predicted Proton Candidates",
99,
text_title_pdf4_string, 
can,
Topology_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
  
  
  
     DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_CCZeroP_Else_BTDPrediction",
  "Events",
  "BDT Prediction Score",
   !DoBinWidthNorm, 
  "RECO tracks BDT Predicted Else Category",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);
  
  
  DrawStackMCandData_EventCategory(
TFile_MC_test4,
TFile_Data4,
TFile_Data_BG4,
"h_CCZeroP_Else_BTDPrediction",
"Events",
"BDT Prediction Score",
 !DoBinWidthNorm, 
"RECO tracks BDT Predicted Else Category",
99,
text_title_pdf4_string, 
can,
EventCategory_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
  
  
  
  
  
   
     DrawStackMCandData_ParticleReduced(
  TFile_Data4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_CCZeroP_Else_BTDPrediction",
  "Events",
  "BDT Prediction Score",
   !DoBinWidthNorm, 
  "RECO tracks BDT Predicted Else : NuWro Truth",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  1.0,
  bd_scale);
  
  
  DrawStackMCandData_EventCategory(
TFile_Data4,
TFile_Data4,
TFile_Data_BG4,
"h_CCZeroP_Else_BTDPrediction",
"Events",
"BDT Prediction Score",
 !DoBinWidthNorm, 
"RECO tracks BDT Predicted Else : NuWro Truth",
99,
text_title_pdf4_string, 
can,
EventCategory_NamedCategory_vector,
1.0,
bd_scale);


sprintf(inputName, "%s",InputREAL_DataFile_name.c_str());
TFile *TFile_Data_REAL = new TFile(inputName);

  DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data_REAL,
  TFile_Data_BG4,
  "h_CCZeroP_Else_BTDPrediction",
  "Events",
  "BDT Prediction Score",
   !DoBinWidthNorm, 
  "RECO tracks BDT Predicted Else Category (REAL DATA)",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);
  
  
  DrawStackMCandData_EventCategory(
TFile_MC_test4,
TFile_Data_REAL,
TFile_Data_BG4,
"h_CCZeroP_Else_BTDPrediction",
"Events",
"BDT Prediction Score",
 !DoBinWidthNorm, 
"RECO tracks BDT Predicted Else Category (REAL DATA)",
99,
text_title_pdf4_string, 
can,
EventCategory_NamedCategory_vector,
POT_MC_Scale,
bd_scale);

  
       DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_topological_score",
  "Events",
  "Topological score",
   !DoBinWidthNorm, 
  "Leading Muon Candidate",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);
  
  
          DrawStackMCandData_EventCategory(
         TFile_MC_test4,
         TFile_Data4,
         TFile_Data_BG4,
         "h_topological_score",
         "Events",
         "Topological score",
          !DoBinWidthNorm, 
         "Leading Muon Candidate",
         99,
         text_title_pdf4_string, 
         can,
         EventCategory_NamedCategory_vector,
         POT_MC_Scale,
         bd_scale);
  
  
  
  DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_trk_score_v",
  "Reconstructed particle candidates",
  "Track score",
   !DoBinWidthNorm, 
  "",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);
  
      DrawStackMCandData_EventCategory(
         TFile_MC_test4,
         TFile_Data4,
         TFile_Data_BG4,
         "h_trk_score_v",
         "Events",
         "Track score",
          !DoBinWidthNorm, 
         "Leading Proton Candidate",
         99,
         text_title_pdf4_string, 
         can,
         EventCategory_NamedCategory_vector,
         POT_MC_Scale,
         bd_scale);
  
  
  
  DrawStackMCandData_ParticleReduced(
  TFile_MC_test4,
  TFile_Data4,
  TFile_Data_BG4,
  "h_trk_len_v",
  "Reconstructed particle candidates",
  "TrackLength [cm]",
   !DoBinWidthNorm, 
  "",
  99,
  text_title_pdf4_string, 
  can,
  ReduceParticle_NamedCategory_vector,
  POT_MC_Scale,
  bd_scale);
  
    
      DrawStackMCandData_EventCategory(
         TFile_MC_test4,
         TFile_Data4,
         TFile_Data_BG4,
         "h_trk_len_v",
         "Events",
         "TrackLength [cm]",
          !DoBinWidthNorm, 
         "Reconstructed particle candidatess",
         99,
         text_title_pdf4_string, 
         can,
         EventCategory_NamedCategory_vector,
         POT_MC_Scale,
         bd_scale);
  
  
  
             DrawStackMCandData_ParticleReduced(
             TFile_MC_test4,
             TFile_Data4,
             TFile_Data_BG4,
             "h_trk_distance_v",
             "Reconstructed particle candidates",
             "Track start distance from vertex [cm]",
              !DoBinWidthNorm, 
             "",
             99,
             text_title_pdf4_string, 
             can,
             ReduceParticle_NamedCategory_vector,
             POT_MC_Scale,
             bd_scale);
  
        DrawStackMCandData_EventCategory(
         TFile_MC_test4,
         TFile_Data4,
         TFile_Data_BG4,
         "h_trk_distance_v",
         "Events",
         "Track start distance from vertex [cm]",
          !DoBinWidthNorm, 
         "Reconstructed particle candidates",
         99,
         text_title_pdf4_string, 
         can,
         EventCategory_NamedCategory_vector,
         POT_MC_Scale,
         bd_scale);
  


////////////////////////////////////////
// Make some Mig
///////////////////////////////////////
bool includeOverFlow= false; 
TH2D* h_p_mig =  GetTH2DHist(*TFile_MC_test4, "h_p_Mig" );
TH2D* h_p_MCS_mig =  GetTH2DHist(*TFile_MC_test4, "h_p_MCS_Mig" );


TH2D* h_costheta_p_mig =  GetTH2DHist(*TFile_MC_test4, "h_costheta_p_Mig" );
TH2D* h_costheta_p_MCS_mig =  GetTH2DHist(*TFile_MC_test4, "h_costheta_p_MCS_Mig" );
TH2D* h_pionE_Mig =  GetTH2DHist(*TFile_MC_test4, "h_pionE_Mig" );



Draw_heatMap(
  h_p_mig, 
  "TRUE P_{#mu} [GeV/c]",
  "RECO P_{#mu} [GeV/c]",
  "dE/dX range + MCS [Col Norm]",
  text_title_pdf2,
  0,
  can, includeOverFlow);

Draw_heatMap(
  h_p_mig, 
  "TRUE P_{#mu} [GeV/c]",
  "RECO P_{#mu} [GeV/c]",
  "dE/dX range + MCS [Row Norm]",
  text_title_pdf2,
  1,
  can, includeOverFlow);

Draw_heatMap(
  h_p_MCS_mig, 
  "TRUE P_{#mu} [GeV/c]",
  "RECO P_{#mu} [GeV/c]",
  "MCS ONLY [Col Norm]",
  text_title_pdf2,
  0,
  can, includeOverFlow);

Draw_heatMap(
  h_p_MCS_mig, 
  "TRUE P_{#mu} [GeV/c]",
  "RECO P_{#mu} [GeV/c]",
  "MCS ONLY [Row Norm]",
  text_title_pdf2,
  1,
  can, includeOverFlow);


Draw_heatMap(
  h_costheta_p_mig, 
  "TRUE Cos(#theta_{#mu})",
  "RECO Cos(#theta_{#mu})",
  "dE/dX range + MCS [Col Norm]",
  text_title_pdf2,
  0,
  can, includeOverFlow);

Draw_heatMap(
  h_costheta_p_mig, 
  "TRUE Cos(#theta_{#mu})",
  "RECO Cos(#theta_{#mu}) ",
  "dE/dX range + MCS [Row Norm]",
  text_title_pdf2,
  1,
  can, includeOverFlow);


Draw_heatMap(
  h_costheta_p_MCS_mig, 
  "TRUE Cos(#theta_{#mu})",
  "RECO Cos(#theta_{#mu})",
  "MCS ONLY [Col Norm]",
  text_title_pdf2,
  0,
  can, includeOverFlow);

Draw_heatMap(
  h_costheta_p_MCS_mig, 
  "TRUE Cos(#theta_{#mu})",
  "RECO Cos(#theta_{#mu}) ",
  "MCS ONLY [Row Norm]",
  text_title_pdf2,
  1,
  can, includeOverFlow);

Draw_heatMap(
  h_pionE_Mig, 
  "TRUE P_{#pi} [GeV/c]",
  "RECO P_{#pi} [GeV/c]",
  "dE/dX range + MCS [raw]",
  text_title_pdf2,
  5,
  can, includeOverFlow);

Draw_heatMap(
  h_pionE_Mig, 
  "TRUE P_{#pi} [GeV/c]",
  "RECO P_{#pi} [GeV/c]",
  "dE/dX range + MCS [Col Norm]",
  text_title_pdf2,
  0,
  can, includeOverFlow);
  
  
  Draw_heatMap(
  h_pionE_Mig, 
  "TRUE P_{#pi} [GeV/c]",
  "RECO P_{#pi} [GeV/c]",
  "dE/dX range + MCS [Row Norm]",
  text_title_pdf2,
  1,
  can, includeOverFlow);




TH1D* h_p_Resolution =  GetTH1DHist(*TFile_MC_test4, "h_p_Resolution" );
TH1D* h_p_MCS_Resolution =  GetTH1DHist(*TFile_MC_test4, "h_p_MCS_Resolution" );

TH1D* h_costheta_p_Resolution =  GetTH1DHist(*TFile_MC_test4, "h_costheta_p_Resolution" );
TH1D* h_costheta_p_MCS_Resolution =  GetTH1DHist(*TFile_MC_test4, "h_costheta_p_MCS_Resolution" );

bool NormArea = true; 
bool Setgrid = true; 

Draw_HIST_Resolution(
  h_p_Resolution,
   "Resolution dE/dx + MCS",
   "RECO - TRUE P_{#mu} [GeV/c]",
   "NEvents [Area Norm]",
   text_title_pdf2_string,
   NormArea, 
    Setgrid,
   DoBinWidthNorm,
   -99,
   can);

Draw_HIST_Resolution(
  h_p_MCS_Resolution,
   "Resolution MCS Only",
   "RECO - TRUE P_{#mu} [GeV/c]",
   "NEvents [Area Norm]",
   text_title_pdf2_string,
   NormArea, 
    Setgrid,
   DoBinWidthNorm,
   -99,
   can);


Draw_HIST_Resolution(
  h_costheta_p_Resolution,
   "Resolution dE/dx + MCS",
   "RECO - TRUE Cos(#theta_{#mu})",
   "NEvents [Area Norm]",
   text_title_pdf2_string,
   NormArea, 
    Setgrid,
   DoBinWidthNorm,
   -99,
   can);

Draw_HIST_Resolution(
  h_costheta_p_MCS_Resolution,
   "Resolution MCS Only",
   "RECO - TRUE Cos(#theta_{#mu})",
   "NEvents [Area Norm]",
   text_title_pdf2_string,
   NormArea, 
    Setgrid,
   DoBinWidthNorm,
   -99,
   can);
   
   
   
   
   Draw2DStackMCandData_EventCategory(
TFile_MC_test4,
TFile_Data4,
TFile_Data_BG4,
"h_costheta_Pmu",
"P#mu",
"Cos(#theta_{#mu})",
 DoBinWidthNorm, 
"title",
1000,
text_title_pdf4_string, 
can,
EventCategory_NamedCategory_vector,
POT_MC_Scale,
bd_scale);
   
   


/*
auto h_NTrack_Map = MakeNamedCategoryMap_TH1D(*TFile_Data4, 
"h_NTracks_BTDGroup_categories", BTDPrectionType_NamedCategory_vector );

TH1D* h_NTrack =  GetTH1DHist(*TFile_Data4, "h_NTracks" );

THStack* stack_input;
std::cout<<"going throught loop "<< std::endl;
int hh=0; 
  for(auto HistInMap : h_NTrack_Map){
  std::cout<<"trying to get Index"<< hh<< std::endl;
   // EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    TH1D* inputhist = (TH1D*)HistInMap.second->Clone(HistInMap.second->GetTitle()); 
    stack_input->Add(inputhist); 
    //EventCategory_stack->Add( HistInMap.second );
    //lg1_stacked->AddEntry(HistInMap.second, HistInMap.second->GetTitle() , "f" );
    hh++;
    
  }

std::cout<<"Finished "<< std::endl
;
    DrawStack(
    h_NTrack, 
    stack_input,
   "Nevents",
    "track [num]",
    false, 
    "BTD Predictions of Fake Data",
    99);
    
    */

//can -> Print(text_title_pdf2);

TFile_MC_test4->Close();
TFile_Data4->Close();
TFile_Data_BG4->Close();

delete TFile_MC_test4;
delete TFile_Data4;
delete TFile_Data_BG4;
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

    TLegend* lg1_stacked = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg1_stacked->SetNColumns(2);
    lg1_stacked->SetBorderSize(0);
    Double_t defaultTextSize = lg1_stacked->GetTextSize();
    //std::cout<<"set text size "<< defaultTextSize << std::endl;
    lg1_stacked->SetTextSize(defaultTextSize*1.05); //
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
  TLegend* lg1_stacked = new TLegend( 0.38, 0.58, 0.88, 0.87 );
  lg1_stacked->SetNColumns(2);
  lg1_stacked->SetBorderSize(0);
  Double_t defaultTextSize = lg1_stacked->GetTextSize();
  lg1_stacked->SetTextSize(.03); //
  
  const auto& DrawStackMCandData_EventCategory_tool = EventCategoryInterpreter::Instance();
  std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap_clone =  MakeNamedCategoryMap_TH1D(*TFile_MC, BaseStackName, inputEventCateory_vector );
  // Set Hist to Stevens Style 
  std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap;
    for(auto map : h_p_EventCategory_StackMap_clone){
  TH1D* hist = (TH1D*) map.second->Clone(uniq());
  hist->SetTitle(map.second->GetTitle());
   h_p_EventCategory_StackMap.insert(std::pair<EventCategory, TH1D*>(map.first,hist)); 
  }
  
  
  TH1D* h_Data_BeamOn_1 =  GetTH1DHist(*TFile_Data, InputHisBaseName );
  TH1D* h_Data_BeamOFF_1 =  GetTH1DHist(*TFile_Data_BeamOff, InputHisBaseName );
  TH1D* h_MC_Total_1 =  GetTH1DHist(*TFile_MC, InputHisBaseName );
  
  TH1D* h_Data_BeamOn = (TH1D*) h_Data_BeamOn_1->Clone("h_Data_BeamOn");
  TH1D* h_Data_BeamOFF = (TH1D*) h_Data_BeamOFF_1->Clone("h_Data_BeamOFF");
  TH1D* h_MC_Total = (TH1D*) h_MC_Total_1->Clone("h_MC_Total");
  

  h_MC_Total->Scale(POT_scaler_MC); 
  h_Data_BeamOFF->Scale(BG_Trigger_scaler); 
   
  float Amount_BeamOn = h_Data_BeamOn->Integral();
  float Amount_BeamOff = h_Data_BeamOFF->Integral();
  float Amount_MC_Total =  h_MC_Total->Integral();
   
   std::cout<< "Amount_BeamOn Data = "<< Amount_BeamOn<< std::endl;
   std::cout<< "Amount_BeamOff Data = "<< Amount_BeamOff<< std::endl;
   std::cout<< "Amount_BeamOn MC = "<< Amount_MC_Total<< std::endl;
   std::cout<< "Ratio Data BEAM on  / MC  "<< Amount_BeamOn/Amount_MC_Total<< std::endl;
  /////////////////////////////////////////////
  // POT Scale hist , all but Data Beam Off
  ////////////////////////////////////////////

  h_MC_Total->Add(h_Data_BeamOFF);
  double Total_Area_MC = h_MC_Total->Integral(); 
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
  EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB
  char leg_title[1024];

  for(auto &HistInMap : h_p_EventCategory_StackMap){
    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    double FractionArea = (HistInMap.second->Integral() / Total_Area_MC) * 100 ;
    sprintf(leg_title, " %s (%2.1f %)",HistInMap.second->GetTitle(),FractionArea );
    EventCategory_stack->Add( HistInMap.second );
    
    
    
    lg1_stacked->AddEntry(HistInMap.second, leg_title , "f" );
  }
  
   h_MC_Total->GetXaxis()->SetTitle(Xaxis_title);
   h_Data_BeamOn->GetXaxis()->SetTitle( Xaxis_title );
   h_Data_BeamOFF->GetXaxis()->SetTitle( Xaxis_title );
  
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
  //pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
  
    h_MC_Total->SetTitle(Title);
    h_Data_BeamOn->SetTitle(Title);
  
  DrawStackMCandData(
      h_Data_BeamOn, 
      h_MC_Total,
      h_Data_BeamOFF,
      EventCategory_stack,
      Yaxis_title,
      Xaxis_title,
     DoBinWidthNorm, 
      Title,
      YMax,
      Canvas);
  
  
  DrawFakeData();
     
 
    double FractionArea_BNB = (h_Data_BeamOFF->Integral() / Total_Area_MC) * 100 ;
    sprintf(leg_title, " %s (%2.1f %)", "EXT BNB", FractionArea_BNB );
    
   lg1_stacked->AddEntry(h_Data_BeamOFF, leg_title , "f" );
   lg1_stacked->Draw( "same" );
   
    TLatex* text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextColor(kRed);
        
    float NDataEvents =  h_Data_BeamOn->Integral();
    sprintf(textplace, "N Data Events = %.1f",NDataEvents );
    text->DrawLatex(0.14, 0.75, textplace);
      
      float NMCEvents =  h_MC_Total->Integral();
      sprintf(textplace, "N MC Events = %.1f",NMCEvents );
       text->DrawLatex(0.14, 0.72, textplace);
       Double_t res[Nbins];
      auto CHIsqrResult =  h_Data_BeamOn->Chi2Test(h_MC_Total,"CHI2",res);
   
     sprintf(textplace, "#chi^{2} / ndf = %.2f / %i",CHIsqrResult, Nbins );
     text->DrawLatex(0.14, 0.68, textplace);
     
    Canvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.2);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.1); 
    pad2->SetGridx(); // vertical grid
    //pad2->SetGridy(); // vertical grid
    pad2->Draw();
    pad2->cd();    
     
    TH1D* h_Data_BeamOn_clone =(TH1D*)h_Data_BeamOn->Clone("h_Data_BeamOn_clone");
    TH1D* h_MC_Total_clone =(TH1D*)h_MC_Total->Clone("h_MC_Total_clone");
     h_Data_BeamOn_clone->SetTitle("");
     h_Data_BeamOn_clone->Divide(h_MC_Total_clone);
     h_Data_BeamOn_clone->GetYaxis()->SetTitle("DATA / MC ");
     h_Data_BeamOn_clone->GetXaxis()->SetTitle(Xaxis_title);
     h_Data_BeamOn_clone->GetYaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleSize(0.125);
     h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(0.09);
   
     
     //double txtsize = h_Data_BeamOn_clone->GetYaxis()->GetTitleSize(.1);
      // h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(.1);
     //std::cout<< " txtsize = "<< txtsize << std::endl;
     
     h_Data_BeamOn_clone->GetXaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->SetTitleOffset(.4);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleOffset(.8);
     h_Data_BeamOn_clone->SetMinimum(.5);
     h_Data_BeamOn_clone->SetMaximum(1.5);
     h_Data_BeamOn_clone->Draw("Hist");
     TLine ratio_1;
     ratio_1.SetLineWidth(2);
     ratio_1.SetLineStyle(7);
     ratio_1.SetLineColor(kBlue);
     double line_min = h_Data_BeamOn_clone->GetBinLowEdge(1);
     double line_max = h_Data_BeamOn_clone->GetBinLowEdge(h_Data_BeamOn_clone->GetNbinsX()+1);
     ratio_1.DrawLine(line_min,1,line_max,1);

     
   Canvas->cd();  
     
 sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  Canvas -> Print(pdf_title);

}
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////
	void Draw2DStackMCandData_EventCategory(
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
  char pdf_title[1024];
  sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  TLegend* lg1_stacked = new TLegend( 0.38, 0.63, 0.88, 0.87 );
  lg1_stacked->SetNColumns(2);
  lg1_stacked->SetBorderSize(0);
  Double_t defaultTextSize = lg1_stacked->GetTextSize();
  lg1_stacked->SetTextSize(.03); //
  
  const auto& DrawStackMCandData_EventCategory_tool = EventCategoryInterpreter::Instance();
  std::cout <<"making th2d Stack Map"<< std::endl;
  std::map<EventCategory, TH2D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH2D(*TFile_MC, BaseStackName, inputEventCateory_vector );
  // Set Hist to Stevens Style 
    std::cout <<"Fin th2d Stack Map"<< std::endl;
  
  TH2D* h_Data_BeamOn_1 =  GetTH2DHist(*TFile_Data, InputHisBaseName );
  TH2D* h_Data_BeamOFF_1 =  GetTH2DHist(*TFile_Data_BeamOff, InputHisBaseName );
  TH2D* h_MC_Total_1 =  GetTH2DHist(*TFile_MC, InputHisBaseName );
  
  TH2D* h_Data_BeamOn = (TH2D*) h_Data_BeamOn_1->Clone("h_Data_BeamOn");
  TH2D* h_Data_BeamOFF = (TH2D*) h_Data_BeamOFF_1->Clone("h_Data_BeamOFF");
  TH2D* h_MC_Total = (TH2D*) h_MC_Total_1->Clone("h_MC_Total");
  
  h_MC_Total->SetTitle("");
  h_Data_BeamOn->SetTitle("");
  h_MC_Total->Scale(POT_scaler_MC); 
  h_Data_BeamOFF->Scale(BG_Trigger_scaler); 
   
  float Amount_BeamOn = h_Data_BeamOn->Integral();
  float Amount_BeamOff = h_Data_BeamOFF->Integral();
  float Amount_MC_Total =  h_MC_Total->Integral();
   
   std::cout<< "Amount_BeamOn Data = "<< Amount_BeamOn<< std::endl;
   std::cout<< "Amount_BeamOff Data = "<< Amount_BeamOff<< std::endl;
   std::cout<< "Amount_BeamOn MC = "<< Amount_MC_Total<< std::endl;
   std::cout<< "Ratio Data BEAM on  / MC  "<< Amount_BeamOn/Amount_MC_Total<< std::endl;
  /////////////////////////////////////////////
  // POT Scale hist , all but Data Beam Off
  ////////////////////////////////////////////

  //h_MC_Total->Add(h_Data_BeamOFF);
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  //THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
 // EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB

TObjArray stack_input;  

  for(auto &HistInMap : h_p_EventCategory_StackMap){
    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    TH2D* inputhist = (TH2D*)HistInMap.second->Clone(HistInMap.second->GetTitle()); 
    stack_input.Add(inputhist); 
    //EventCategory_stack->Add( HistInMap.second );
    lg1_stacked->AddEntry(HistInMap.second, HistInMap.second->GetTitle() , "f" );
  }
   h_MC_Total->GetXaxis()->SetTitle(Xaxis_title);
   h_Data_BeamOn->GetXaxis()->SetTitle( Xaxis_title );
   h_Data_BeamOFF->GetXaxis()->SetTitle( Xaxis_title );
  
  
   // TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.0);
   // pad1->SetBottomMargin(0); // Upper and lower plot are joined
  ////pad1->SetGridx();         // Vertical grid
   // pad1->Draw();             // Draw the upper pad: pad1
   // pad1->cd();               // pad1 becomes the current pad
  
  std::vector<double> YMultipliers{1,1,1,1,1,1,1,1,1,1,1,1};
  std::vector<int> fillColors{1,1,1,1,1,1,1,1,1,1,1,1};
  
  std::cout<<"running PlotDataStackedMC2D_ProjY"<< std::endl;
  bool useMult = true; 
  
  PlotDataStackedMC2D_ProjY(
	  h_Data_BeamOn,
	 h_Data_BeamOFF,
	 h_MC_Total,
	stack_input,
 fillColors,
	pdf_title, "title", Xaxis_title,
	Yaxis_title, "NEvents",
	2000, true, useMult,
	 YMultipliers,
	DoBinWidthNorm,
	.02,
	1.0);
  
  
  
  std::cout<<"Finished PlotDataStackedMC2D_ProjY"<< std::endl;
  
  
    PlotDataStackedMC2D_ProjX(
	  h_Data_BeamOn,
	 h_Data_BeamOFF,
	 h_MC_Total,
	stack_input,
 fillColors,
	pdf_title, "title", Xaxis_title,
	Yaxis_title, "NEvents",
	2000, true, useMult,
	 YMultipliers,
	DoBinWidthNorm,
	.02,
	1.0);

   //lg1_stacked->AddEntry(h_Data_BeamOFF,"EXT BNB" , "f" );
   //lg1_stacked->Draw( "same" );
   /*
    TLatex* text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextColor(kRed);
        
    float NDataEvents =  h_Data_BeamOn->Integral();
    sprintf(textplace, "N Data Events = %.1f",NDataEvents );
    text->DrawLatex(0.14, 0.75, textplace);
      
      float NMCEvents =  h_MC_Total->Integral();
      sprintf(textplace, "N MC Events = %.1f",NMCEvents );
       text->DrawLatex(0.14, 0.72, textplace);
       Double_t res[Nbins];
      auto CHIsqrResult =  h_Data_BeamOn->Chi2Test(h_MC_Total,"CHI2",res);
   
     sprintf(textplace, "#chi^{2} / ndf = %.2f / %i",CHIsqrResult, Nbins );
     text->DrawLatex(0.14, 0.68, textplace);
     
    Canvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.2);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.1); 
    pad2->SetGridx(); // vertical grid
    //pad2->SetGridy(); // vertical grid
    pad2->Draw();
    pad2->cd();    
     
    TH1D* h_Data_BeamOn_clone =(TH1D*)h_Data_BeamOn->Clone("h_Data_BeamOn_clone");
    TH1D* h_MC_Total_clone =(TH1D*)h_MC_Total->Clone("h_MC_Total_clone");

    h_Data_BeamOn_clone->Divide(h_MC_Total_clone);
     h_Data_BeamOn_clone->GetYaxis()->SetTitle("DATA / MC ");
     h_Data_BeamOn_clone->GetXaxis()->SetTitle(Xaxis_title);
     h_Data_BeamOn_clone->GetYaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleSize(0.125);
     h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(0.09);
   
     
     //double txtsize = h_Data_BeamOn_clone->GetYaxis()->GetTitleSize(.1);
      // h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(.1);
     //std::cout<< " txtsize = "<< txtsize << std::endl;
     
     h_Data_BeamOn_clone->GetXaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->SetTitleOffset(.4);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleOffset(.8);
     h_Data_BeamOn_clone->SetMinimum(.5);
     h_Data_BeamOn_clone->SetMaximum(1.5);
     h_Data_BeamOn_clone->Draw("Hist");
     TLine ratio_1;
     ratio_1.SetLineWidth(2);
     ratio_1.SetLineStyle(7);
     ratio_1.SetLineColor(kBlue);
     double line_min = h_Data_BeamOn_clone->GetBinLowEdge(1);
     double line_max = h_Data_BeamOn_clone->GetBinLowEdge(h_Data_BeamOn_clone->GetNbinsX()+1);
     ratio_1.DrawLine(line_min,1,line_max,1);

     
   Canvas->cd();  
     */
 
  //Canvas -> Print(pdf_title);
return; 
}
	
////////////////////////////////////////////
//
/////////////////////////////////////////
void DrawStackMCandData_Containment(
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
std::vector<NamedCategory<FidVol_Category>> inputEventCateory_vector,
float POT_scaler_MC,
float BG_Trigger_scaler)
{
  std::string BaseStackName = InputHisBaseName+ "_containment_categories";
  std::cout << "Inside::DrawStackMCandData_containment_"<< std::endl;
  std::cout << "Inside:input X axis title"<< Xaxis_title<<std::endl;
  TLegend* lg1_stacked = new TLegend( 0.38, 0.63, 0.88, 0.87 );
  lg1_stacked->SetNColumns(2);
  lg1_stacked->SetBorderSize(0);
  Double_t defaultTextSize = lg1_stacked->GetTextSize();
  lg1_stacked->SetTextSize(.03); //
  
  const auto& DrawStackMCandData_EventCategory_tool = EventCategoryInterpreter::Instance();
  std::map<FidVol_Category, TH1D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH1D(*TFile_MC, BaseStackName, inputEventCateory_vector );
  // Set Hist to Stevens Style 
  
  
  TH1D* h_Data_BeamOn =  GetTH1DHist(*TFile_Data, InputHisBaseName );
  TH1D* h_Data_BeamOFF =  GetTH1DHist(*TFile_Data_BeamOff, InputHisBaseName );
  TH1D* h_MC_Total =  GetTH1DHist(*TFile_MC, InputHisBaseName );
  h_MC_Total->SetTitle("");
  h_Data_BeamOn->SetTitle("");
  h_MC_Total->Scale(POT_scaler_MC); 
  h_Data_BeamOFF->Scale(BG_Trigger_scaler); 
   
  float Amount_BeamOn = h_Data_BeamOn->Integral();
  float Amount_BeamOff = h_Data_BeamOFF->Integral();
  float Amount_MC_Total =  h_MC_Total->Integral();
   
   std::cout<< "Amount_BeamOn Data = "<< Amount_BeamOn<< std::endl;
   std::cout<< "Amount_BeamOff Data = "<< Amount_BeamOff<< std::endl;
   std::cout<< "Amount_BeamOn MC = "<< Amount_MC_Total<< std::endl;
   std::cout<< "Ratio Data BEAM on  / MC  "<< Amount_BeamOn/Amount_MC_Total<< std::endl;
  /////////////////////////////////////////////
  // POT Scale hist , all but Data Beam Off
  ////////////////////////////////////////////

  h_MC_Total->Add(h_Data_BeamOFF);
    double Total_Area_MC = h_MC_Total->Integral(); 
  
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  std::cout<<"this far : 1"<< std::endl;
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
  EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB

  char leg_title[1024];
  
  for(auto &HistInMap : h_p_EventCategory_StackMap){
    //DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    double FractionArea = (HistInMap.second->Integral() / Total_Area_MC) * 100 ;
    sprintf(leg_title, " %s (%2.1f %)",HistInMap.second->GetTitle(),FractionArea );

   // EventCategory_stack->Add( HistInMap.second );
    lg1_stacked->AddEntry(HistInMap.second, leg_title , "f" );
  }
  
     for (auto it = h_p_EventCategory_StackMap.rbegin(); it != h_p_EventCategory_StackMap.rend(); ++it) {
    // 'it' is a reverse iterator pointing to the current element in reverse order
    auto& HistInMap = *it;

    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
}
  
  
  
  
   h_MC_Total->GetXaxis()->SetTitle(Xaxis_title);
   h_Data_BeamOn->GetXaxis()->SetTitle( Xaxis_title );
   h_Data_BeamOn->SetTitle( Title );
   h_Data_BeamOFF->GetXaxis()->SetTitle( Xaxis_title );
  
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
  //pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
  

  
  DrawStackMCandData(
      h_Data_BeamOn, 
      h_MC_Total,
      h_Data_BeamOFF,
      EventCategory_stack,
      Yaxis_title,
      Xaxis_title,
     DoBinWidthNorm, 
      Title,
      YMax,
      Canvas);
      
      
   double FractionArea_BNB = (h_Data_BeamOFF->Integral() / Total_Area_MC) * 100 ;
    sprintf(leg_title, " %s (%2.1f %)", "EXT BNB", FractionArea_BNB );
   lg1_stacked->AddEntry(h_Data_BeamOFF,leg_title , "f" );
   lg1_stacked->Draw( "same" );
   
   
   DrawFakeData();
   
    TLatex* text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextColor(kRed);
        
    float NDataEvents =  h_Data_BeamOn->Integral();
    sprintf(textplace, "N Data Events = %.1f",NDataEvents );
    text->DrawLatex(0.14, 0.75, textplace);
      
      float NMCEvents =  h_MC_Total->Integral();
      sprintf(textplace, "N MC Events = %.1f",NMCEvents );
       text->DrawLatex(0.14, 0.72, textplace);
       Double_t res[Nbins];
      auto CHIsqrResult =  h_Data_BeamOn->Chi2Test(h_MC_Total,"CHI2",res);
   
     sprintf(textplace, "#chi^{2} / ndf = %.2f / %i",CHIsqrResult, Nbins );
     text->DrawLatex(0.14, 0.68, textplace);
     
    Canvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.2);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.1); 
    pad2->SetGridx(); // vertical grid
    //pad2->SetGridy(); // vertical grid
    pad2->Draw();
    pad2->cd();    
     
    TH1D* h_Data_BeamOn_clone =(TH1D*)h_Data_BeamOn->Clone("h_Data_BeamOn_clone");
    TH1D* h_MC_Total_clone =(TH1D*)h_MC_Total->Clone("h_MC_Total_clone");

     h_Data_BeamOn_clone->Divide(h_MC_Total_clone);
     h_Data_BeamOn_clone->SetTitle("");
     h_Data_BeamOn_clone->GetYaxis()->SetTitle("DATA / MC ");
     h_Data_BeamOn_clone->GetXaxis()->SetTitle(Xaxis_title);
     h_Data_BeamOn_clone->GetYaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleSize(0.125);
     h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(0.09);
   
     
     //double txtsize = h_Data_BeamOn_clone->GetYaxis()->GetTitleSize(.1);
      // h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(.1);
     //std::cout<< " txtsize = "<< txtsize << std::endl;
     
     h_Data_BeamOn_clone->GetXaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->SetTitleOffset(.4);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleOffset(.8);
     h_Data_BeamOn_clone->SetMinimum(.5);
     h_Data_BeamOn_clone->SetMaximum(1.5);
     h_Data_BeamOn_clone->Draw("Hist");
     TLine ratio_1;
     ratio_1.SetLineWidth(2);
     ratio_1.SetLineStyle(7);
     ratio_1.SetLineColor(kBlue);
     double line_min = h_Data_BeamOn_clone->GetBinLowEdge(1);
     double line_max = h_Data_BeamOn_clone->GetBinLowEdge(h_Data_BeamOn_clone->GetNbinsX()+1);
     ratio_1.DrawLine(line_min,1,line_max,1);

     
   Canvas->cd();  
     
 sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  Canvas -> Print(pdf_title);

}
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void DrawStackMCandData_Topology(
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
std::vector<NamedCategory<CCZeroPi_type>> inputEventCateory_vector,
float POT_scaler_MC,
float BG_Trigger_scaler)
{
  std::string BaseStackName = InputHisBaseName+ "_topology_categories";
  
  std::cout << "Inside::DrawStackMCandData_containment_"<< std::endl;
  std::cout << "Inside:input X axis title"<< Xaxis_title<<std::endl;
  TLegend* lg1_stacked = new TLegend( 0.38, 0.63, 0.88, 0.87 );
  lg1_stacked->SetNColumns(2);
  lg1_stacked->SetBorderSize(0);
  Double_t defaultTextSize = lg1_stacked->GetTextSize();
  lg1_stacked->SetTextSize(.03); //
  
  const auto& DrawStackMCandData_EventCategory_tool = EventCategoryInterpreter::Instance();
  std::map<CCZeroPi_type, TH1D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH1D(*TFile_MC, BaseStackName, inputEventCateory_vector );
  // Set Hist to Stevens Style 
  
  
  TH1D* h_Data_BeamOn =  GetTH1DHist(*TFile_Data, InputHisBaseName );
  TH1D* h_Data_BeamOFF =  GetTH1DHist(*TFile_Data_BeamOff, InputHisBaseName );
  TH1D* h_MC_Total =  GetTH1DHist(*TFile_MC, InputHisBaseName );
  h_MC_Total->SetTitle("");
  h_Data_BeamOn->SetTitle("");
  h_MC_Total->Scale(POT_scaler_MC); 
  h_Data_BeamOFF->Scale(BG_Trigger_scaler); 
   
  float Amount_BeamOn = h_Data_BeamOn->Integral();
  float Amount_BeamOff = h_Data_BeamOFF->Integral();
  float Amount_MC_Total =  h_MC_Total->Integral();
   
   std::cout<< "Amount_BeamOn Data = "<< Amount_BeamOn<< std::endl;
   std::cout<< "Amount_BeamOff Data = "<< Amount_BeamOff<< std::endl;
   std::cout<< "Amount_BeamOn MC = "<< Amount_MC_Total<< std::endl;
   std::cout<< "Ratio Data BEAM on  / MC  "<< Amount_BeamOn/Amount_MC_Total<< std::endl;
  /////////////////////////////////////////////
  // POT Scale hist , all but Data Beam Off
  ////////////////////////////////////////////

  h_MC_Total->Add(h_Data_BeamOFF);
  double Total_Area_MC = h_MC_Total->Integral(); 
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
  EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB
char leg_title[1024];

  for(auto &HistInMap : h_p_EventCategory_StackMap){
    //DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    //EventCategory_stack->Add( HistInMap.second );
    double FractionArea = (HistInMap.second->Integral() / Total_Area_MC) * 100 ;
    sprintf(leg_title, " %s (%2.1f %)",HistInMap.second->GetTitle(),FractionArea );
    lg1_stacked->AddEntry(HistInMap.second, leg_title , "f" );
  }
  
  
  for (auto it = h_p_EventCategory_StackMap.rbegin(); it != h_p_EventCategory_StackMap.rend(); ++it) {
    // 'it' is a reverse iterator pointing to the current element in reverse order
    auto& HistInMap = *it;
    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
}
  
  
  
   h_MC_Total->GetXaxis()->SetTitle(Xaxis_title);
   h_Data_BeamOn->GetXaxis()->SetTitle( Xaxis_title );
   h_Data_BeamOn->SetTitle( Title );
   h_Data_BeamOFF->GetXaxis()->SetTitle( Xaxis_title );
  
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
  //pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
  
  
  DrawStackMCandData(
      h_Data_BeamOn, 
      h_MC_Total,
      h_Data_BeamOFF,
      EventCategory_stack,
      Yaxis_title,
      Xaxis_title,
     DoBinWidthNorm, 
      Title,
      YMax,
      Canvas);
  
      double FractionArea_BNB = (h_Data_BeamOFF->Integral() / Total_Area_MC) * 100 ;
    sprintf(leg_title, " %s (%2.1f %)", "EXT BNB", FractionArea_BNB );
  
   lg1_stacked->AddEntry(h_Data_BeamOFF, leg_title , "f" );
   lg1_stacked->Draw( "same" );
   DrawFakeData();
    TLatex* text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextColor(kRed);
        
    float NDataEvents =  h_Data_BeamOn->Integral();
    sprintf(textplace, "N Data Events = %.1f",NDataEvents );
    text->DrawLatex(0.14, 0.75, textplace);
      
      float NMCEvents =  h_MC_Total->Integral();
      sprintf(textplace, "N MC Events = %.1f",NMCEvents );
       text->DrawLatex(0.14, 0.72, textplace);
       Double_t res[Nbins];
      auto CHIsqrResult =  h_Data_BeamOn->Chi2Test(h_MC_Total,"CHI2",res);
   
     sprintf(textplace, "#chi^{2} / ndf = %.2f / %i",CHIsqrResult, Nbins );
     text->DrawLatex(0.14, 0.68, textplace);
     
    Canvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.2);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.1); 
    pad2->SetGridx(); // vertical grid
    //pad2->SetGridy(); // vertical grid
    pad2->Draw();
    pad2->cd();    
     
    TH1D* h_Data_BeamOn_clone =(TH1D*)h_Data_BeamOn->Clone("h_Data_BeamOn_clone");
    TH1D* h_MC_Total_clone =(TH1D*)h_MC_Total->Clone("h_MC_Total_clone");

     h_Data_BeamOn_clone->Divide(h_MC_Total_clone);
     h_Data_BeamOn_clone->SetTitle("");
     h_Data_BeamOn_clone->GetYaxis()->SetTitle("DATA / MC ");
     h_Data_BeamOn_clone->GetXaxis()->SetTitle(Xaxis_title);
     h_Data_BeamOn_clone->GetYaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleSize(0.125);
     h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(0.09);
   
     
     //double txtsize = h_Data_BeamOn_clone->GetYaxis()->GetTitleSize(.1);
      // h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(.1);
     //std::cout<< " txtsize = "<< txtsize << std::endl;
     
     h_Data_BeamOn_clone->GetXaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->SetTitleOffset(.4);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleOffset(.8);
     h_Data_BeamOn_clone->SetMinimum(.5);
     h_Data_BeamOn_clone->SetMaximum(1.5);
     h_Data_BeamOn_clone->Draw("Hist");
     TLine ratio_1;
     ratio_1.SetLineWidth(2);
     ratio_1.SetLineStyle(7);
     ratio_1.SetLineColor(kBlue);
     double line_min = h_Data_BeamOn_clone->GetBinLowEdge(1);
     double line_max = h_Data_BeamOn_clone->GetBinLowEdge(h_Data_BeamOn_clone->GetNbinsX()+1);
     ratio_1.DrawLine(line_min,1,line_max,1);

     
   Canvas->cd();  
     
 sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  Canvas -> Print(pdf_title);

}
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void DrawStackMCandData_ParticleReduced(
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
std::vector<NamedCategory<Particle_type>> inputEventCateory_vector,
float POT_scaler_MC,
float BG_Trigger_scaler)
{
  std::string BaseStackName = InputHisBaseName+ "_Particle_categories";
  
  std::cout << "Inside::DrawStackMCandData_Particle_categories_"<< std::endl;
  std::cout << "Inside:input X axis title"<< Xaxis_title<<std::endl;
  TLegend* lg1_stacked = new TLegend( 0.38, 0.58, 0.88, 0.88 );
  lg1_stacked->SetNColumns(2);
  lg1_stacked->SetBorderSize(0);
  Double_t defaultTextSize = lg1_stacked->GetTextSize();
  lg1_stacked->SetTextSize(.03); //
  
  std::string MicroBooneTitle = MicroBooNE_LegendTitle();
  lg1_stacked->SetHeader(MicroBooneTitle.c_str());
  
  const auto& DrawStackMCandData_EventCategory_tool = EventCategoryInterpreter::Instance();
  std::map<Particle_type, TH1D*> h_p_EventCategory_StackMap_clone =  MakeNamedCategoryMap_TH1D(*TFile_MC, BaseStackName, inputEventCateory_vector );
  // Set Hist to Stevens Style 
  std::map<Particle_type, TH1D*> h_p_EventCategory_StackMap;
  
  for(auto map : h_p_EventCategory_StackMap_clone){
  TH1D* hist = (TH1D*) map.second->Clone(uniq());
  hist->SetTitle(map.second->GetTitle());
   h_p_EventCategory_StackMap.insert(std::pair<Particle_type, TH1D*>(map.first,hist)); 
  }
  
  
  
  std::cout<<"file Name 1" <<std::endl;
  TH1D* h_Data_BeamOn_1 =  GetTH1DHist(*TFile_Data, InputHisBaseName );
  std::cout<<"file Name 2" <<std::endl;
  TH1D* h_Data_BeamOFF_1 =  GetTH1DHist(*TFile_Data_BeamOff, InputHisBaseName );
    std::cout<<"file Name 3" <<std::endl;
  TH1D* h_MC_Total_1 =  GetTH1DHist(*TFile_MC, InputHisBaseName );
  
  TH1D* h_Data_BeamOn = (TH1D*) h_Data_BeamOn_1->Clone("TH1D* h_Data_BeamOn");
  TH1D* h_Data_BeamOFF = (TH1D*) h_Data_BeamOFF_1->Clone("TH1D* h_Data_BeamOFF");
  TH1D* h_MC_Total = (TH1D*) h_MC_Total_1->Clone("TH1D* h_MC_Total");
  
  
  h_MC_Total->Scale(POT_scaler_MC);
  h_Data_BeamOFF->Scale(BG_Trigger_scaler); 
  
  h_MC_Total->GetXaxis()->SetTitle(Xaxis_title);
  h_Data_BeamOn->GetXaxis()->SetTitle( Xaxis_title );
  h_Data_BeamOFF->GetXaxis()->SetTitle( Xaxis_title );
     
  float Amount_BeamOn = h_Data_BeamOn->Integral();
  float Amount_BeamOff = h_Data_BeamOFF->Integral();
  float Amount_MC_Total =  h_MC_Total->Integral();
   
   std::cout<< "Amount_BeamOn Data = "<< Amount_BeamOn<< std::endl;
   std::cout<< "Amount_BeamOff Data = "<< Amount_BeamOff<< std::endl;
   std::cout<< "Amount_BeamOn MC = "<< Amount_MC_Total<< std::endl;
   std::cout<< "Ratio Data BEAM on  / MC  "<< Amount_BeamOn/Amount_MC_Total<< std::endl;
  /////////////////////////////////////////////
  // POT Scale hist , all but Data Beam Off
  ////////////////////////////////////////////

  h_MC_Total->Add(h_Data_BeamOFF);
  double Total_Area_MC = h_MC_Total->Integral(); 
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  
  TH1D* h_Data_BeamOn_clone =(TH1D*)h_Data_BeamOn->Clone("h_Data_BeamOn_clone");
  TH1D* h_MC_Total_clone =(TH1D*)h_MC_Total->Clone("h_MC_Total_clone");
  
  
  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
  EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB
char leg_title[1024];

  for(auto &HistInMap : h_p_EventCategory_StackMap){
    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
    double FractionArea = (HistInMap.second->Integral() / Total_Area_MC) * 100 ;
    sprintf(leg_title, " %s (%2.1f %)",HistInMap.second->GetTitle(),FractionArea );
    lg1_stacked->AddEntry(HistInMap.second, leg_title, "f" );
  }

  
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
  //pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
  
  
   h_MC_Total->SetTitle(Title);
   h_Data_BeamOn->SetTitle(Title);
  
  
  DrawStackMCandData(
      h_Data_BeamOn, 
      h_MC_Total,
      h_Data_BeamOFF,
      EventCategory_stack,
      Yaxis_title,
      Xaxis_title,
     DoBinWidthNorm, 
      Title,
      YMax,
      Canvas);
      
      double FractionArea_BNB = (h_Data_BeamOFF->Integral() / Total_Area_MC) * 100 ;
    sprintf(leg_title, " %s (%2.1f %)", "EXT BNB", FractionArea_BNB );
  
   lg1_stacked->AddEntry(h_Data_BeamOFF,leg_title , "f" );
   lg1_stacked->Draw( "same" );
   DrawFakeData();
    TLatex* text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextColor(kRed);
        
    float NDataEvents =  h_Data_BeamOn->Integral();
    sprintf(textplace, "N Data Events = %.1f",NDataEvents );
    text->DrawLatex(0.14, 0.75, textplace);
      
      float NMCEvents =  h_MC_Total->Integral();
      sprintf(textplace, "N MC Events = %.1f",NMCEvents );
       text->DrawLatex(0.14, 0.72, textplace);
       Double_t res[Nbins];
      auto CHIsqrResult =  h_Data_BeamOn->Chi2Test(h_MC_Total,"CHI2",res);
   
     sprintf(textplace, "#chi^{2} / ndf = %.2f / %i",CHIsqrResult, Nbins );
     text->DrawLatex(0.14, 0.68, textplace);
     
    Canvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.2);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.1); 
    pad2->SetGridx(); // vertical grid
    pad2->SetGridy(); // vertical grid
    pad2->Draw();
    pad2->cd();    
     
    gStyle->SetOptStat(0);
    h_Data_BeamOn_clone->SetTitle("");
     
     h_Data_BeamOn_clone->Divide(h_MC_Total_clone);
     h_Data_BeamOn_clone->GetYaxis()->SetTitle("DATA / MC ");
     h_Data_BeamOn_clone->GetXaxis()->SetTitle(Xaxis_title);
     h_Data_BeamOn_clone->GetYaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetLabelSize(.1);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleSize(0.125);
     h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(0.09);
   
     
     //double txtsize = h_Data_BeamOn_clone->GetYaxis()->GetTitleSize(.1);
      // h_Data_BeamOn_clone->GetYaxis()->SetTitleSize(.1);
     //std::cout<< " txtsize = "<< txtsize << std::endl;
     
     h_Data_BeamOn_clone->GetXaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->CenterTitle();
     h_Data_BeamOn_clone->GetYaxis()->SetTitleOffset(.4);
     h_Data_BeamOn_clone->GetXaxis()->SetTitleOffset(.8);
     h_Data_BeamOn_clone->SetMinimum(.5);
     h_Data_BeamOn_clone->SetMaximum(1.5);
     h_Data_BeamOn_clone->Draw("Hist");
     TLine ratio_1;
     ratio_1.SetLineWidth(2);
     ratio_1.SetLineStyle(7);
     ratio_1.SetLineColor(kBlue);
     double line_min = h_Data_BeamOn_clone->GetBinLowEdge(1);
     double line_max = h_Data_BeamOn_clone->GetBinLowEdge(h_Data_BeamOn_clone->GetNbinsX()+1);
     ratio_1.DrawLine(line_min,1,line_max,1);

     
   Canvas->cd();  
     
 sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  Canvas -> Print(pdf_title);

}
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
std::string MicroBooNE_LegendTitle(){
return get_legend_title( DATA_POT );
}// End of Funtion
///////////////////////////
void DrawFakeData(){
    TLatex* textfake = new TLatex;
    textfake->SetNDC();
    textfake->SetTextSize(0.033);
    textfake->SetTextColor(kBlue);
    textfake->DrawLatex(0.14, 0.8, "USING FAKE DATA (NuWro)");

}



int main() {


std::cout<<" Testing Slice Plots "<< std::endl;

  //tutorial_slice_plots();
  MakePlots();
  return 0;
}
