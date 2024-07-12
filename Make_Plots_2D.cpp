////////////////////////////////////
// Playing this Slice Plots
///////////////////////////////////

#include "Make_Plots.hh"
#include <TROOT.h>
//using NFT = NtupleFileType;

//#include "HistUtils.hh"

#define USE_FAKE_DATA ""

/// Defined Global EventCateogry 
  
//void DrawBinningInfo(std::map<int , BinMap> TH1Poly_binMap_input);
//void DrawBinningNum(std::map<int , BinMap> TH1Poly_binMap_input);
// I want to print figures to PDF

char pdf_title[1024];
char textplace[1024];
std::string Pdf_name = "Tutorial_slice_Figures";
const double DATA_POT = 4.54e+19;
const double MC_POT = 1.30503e+21;

 
void convertMapToArrays(const std::map<double, std::vector<double>>& inputMap, std::vector<Double_t*>& outputArrays);
void FillBinnMap(std::map<int , BinMap> &InputBinMap, std::map<double, std::vector<double>>  input_BinMap );
void scaleTHStack(THStack* stack, double scaleFactor);
void drawString_extra(std::string inputString, double text_size, int color, bool left, double add_Xpostion, double add_Ypostion );
//void DrawFakeData();
/*
 void PlotMC2D_Stack(
 std::map<Binning2D, TH1D*> InputHistmap_DATA,
 std::map<Binning2D, TH1D*> InputHistmap_DATA_BeamOff,
 std::map<Binning2D, TH1D*> InputHistmap_MC,
 std::map<std::pair<Binning2D, EventCategory>, TH1D*> StackMap_input,
 std::map<Binning2D , std::string > BinStringMap,
	char *pdf_label, char *histotitle, char *xaxislabel,
	char *legendTitle, char *zaxislabel_units,
	double Ymax, bool setMaxY, double Ymin,  bool doMultipliers,
	std::vector<double> YMultipliers,
	bool dontDo_bin_width_norm,
	double text_size,
	float POT_DATA,
float POT_scaler_MC,
float BG_Trigger_scaler);
*/
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

  void PlotFractionofEvents2D_new(
   std::map<std::pair<Binning2D, EventCategory>, TH1D*> StackMap_input,
   std::map<Binning2D, TH1D*> InputHistmap_MC,
   std::map<Binning2D , std::string > BinStringMap,
   double POT_DATA,
   double text_size,
	 char *pdf_label, 
	 char *histotitle, 
	 char *xaxislabel,
	 char *zaxislabel);





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
 void PlotFractionofEvents2D(
   std::map<std::pair<Binning2D, EventCategory>, TH1D*> StackMap_input,
   std::map<Binning2D, TH1D*> InputHistmap_MC,
   std::map<Binning2D , std::string > BinStringMap,
   double POT_DATA,
   double text_size,
	 char *pdf_label, 
	 char *histotitle, 
	 char *xaxislabel,
	 char *zaxislabel);



void MakePlots(){

///////////////////////////////////////////////////
// Getting the POT Only looking at run 1 for now need to Change this Method  
///////////////////////////////////////////////
   std::cout<<"Starting 2D MakePlots"<< std::endl;
    TROOT root("myROOT", "My ROOT Title");
   UBTH2Poly instance;
   gInterpreter->Declare("#include \"includes/UBTH2Poly.h\""); 
   TCanvas *can = new TCanvas("can");
   gStyle->SetOptStat(0); 
    char text_title_pdf1[2024];
    char text_title_pdf2[2024];
    char text_title_pdf3[2024];
    char text_title_pdf4[2024];

  sprintf(text_title_pdf1, "Make_Plots_result_2D.pdf(","" );
  can -> Print(text_title_pdf1);
  sprintf(text_title_pdf2, "Make_Plots_result_2D.pdf","" );
  sprintf(text_title_pdf3, "Make_Plots_result_2D.pdf)","" );
  sprintf(text_title_pdf4, "Make_Plots_result_2D","" );
  std::string text_title_pdf4_string(text_title_pdf4);
  std::string text_title_pdf2_string(text_title_pdf2);
   auto& fpm = FilePropertiesManager::Instance();
    const auto& EventCategory_tool = EventCategoryInterpreter::Instance();
    std::vector<NamedCategory<EventCategory>> EventCategory_NamedCategory_vector= EventCategory_tool.EventSelectionGroup_categories_;


    std::vector<NamedCategory<FidVol_Category>> Containment_NamedCategory_vector= EventCategory_tool.ContainedGroup_categories_;
    std::vector<NamedCategory<CCZeroPi_type>> Topology_NamedCategory_vector= EventCategory_tool.topology_categories_;

    std::vector<NamedCategory<Particle_type>> ReduceParticle_NamedCategory_vector= EventCategory_tool.ParticleGroup_reduced_categories_;

   
   
   std::string InputPath_Data = "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/Data/";
//std::string Data_file1 =  "DATA_Selected_newCC0pi_addthresholds_rustv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root";

std::string Data_file1 =  "DATA_Selected_newCC0pi_rustv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_NuWro_reco2_reco2.root";

std::string Data_file_bkg = "DATA_Selected_newCC0pi_rustv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";
//addthresholds_WithProtonPIDCUT
std::string InputPath_MC = "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/MC/";
std::string MC_file1 =  "MC_Selected_newCC0pi_rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
std::string MC_file2 =  "MC_Selected_rustv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";
std::string MC_file3 =  "MC_Selected_rustv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root";

std::string MC_file4 = "MC_Selected_EventSelection_5_13_2024_CC0pi_ru-stv_overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root";
                        //MC_Selected_newCC0pi_addthresholds_rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
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

std::string file_name_MC =  "/exp/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
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
float bd_scale = 10080350. / 65498807.;
std::cout << " bnb_trigs = " << bnb_trigs << std::endl;
std::cout << " ext_trigs = " << ext_trigs << std::endl;
std::cout << " bg scale = " << bd_scale << std::endl;
std::cout << " bg scale2 = " << bd_scale2 << std::endl;


char inputName[1024];
sprintf(inputName, "%s",InputMCFile_name.c_str());
TFile *TFile_MC = new TFile(inputName);
  //std::vector<TH1D*> check_vector = MakeList_TH1D(*TFile_MC_test, "name",  EventCategory_NamedCategory_vector );


//std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH1D(*TFile_MC_test, "h_p_EventCategory", EventCategory_NamedCategory_vector );

std::vector<std::string> HistNames; 

for(auto bins : MUON_2D_BIN_EDGES){

//std::cout<< "bins = "<< bins.first<< std::endl;
//for(auto xbins:bins.second ){
//std::cout<< xbins<< " , "; 

//}
std::cout<<""<<std::endl;
 

}


    std::vector<Double_t*> outputArrays;

    // Convert the map to arrays
    convertMapToArrays(MUON_2D_BIN_EDGES, outputArrays);

    // Access the arrays
    int xIndex = 0;

   
  std::map<int , BinMap> TH1Poly_binMap_test;
  std::map<int , BinMap> TH1Poly_binMap_test2;
  std::map<int , BinMap> TH1Poly_binMap_test_inclusive1; 
  std::map<int , BinMap> TH1Poly_binMap_test_inclusive2; 
   
   
   
  UBTH2Poly *h2p_UB = Make2DHist_UB(MUON_2D_BIN_EDGES, TH1Poly_binMap_test2, "h2p_UB");
  TH2Poly *h2p_test_old = Make2DHist(MUON_2D_BIN_EDGES, TH1Poly_binMap_test, "h2p_tes");

  std::cout<<"TH1Poly_binMap_test.size() = "<< TH1Poly_binMap_test.size()<< std::endl;

  
  TH2Poly *h2p_test_inclusive = Make2DHist_inclusive(MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test_inclusive1, "h2p_test_inclusive");
  UBTH2Poly  *h2p_test_inclusive_UB = Make2DHist_inclusive_UB(MUON_2D_BIN_EDGES_inclusive, TH1Poly_binMap_test_inclusive2, "h2p_test_inclusive_UB");
  
  
  
  
std::map<int , BinMap> TH1Poly_binMap;
std::map<int , BinMap> TH1Poly_binMap_inclusive;
/////////////////////

///////////////////
 std::cout<<"filling"<<std::endl;
 double count=0.0;
  Double_t countt=0.0;
  for (Int_t bin = 0; bin < h2p_UB->GetNumberOfBins(); ++bin) {
        // Get the bin's information
        Double_t x, y;
       // h2p->GetBin(bin, x, y);
        count=count+1.0;
        countt=countt+1.0;
        // Optionally, you can get more information about the bin, like the area
        //h2p->SetBinContent(bin,count);
        x = TH1Poly_binMap_test[bin].centerx;
        y = TH1Poly_binMap_test[bin].centery;

        
      //h2p->Fill(x, y, count);
      h2p_UB->Fill(x, y, countt);
      //h2p_test->Fill(x, y, countt);
      h2p_test_old->Fill(x, y, count);

  

        // You can perform additional operations on each bin here
  }
    
    
    
gStyle->SetPalette(kTemperatureMap);

//////////////////////////////////////////////////////

  h2p_test_old->GetXaxis()->SetTitle("cos#theta_{#mu}");
  h2p_test_old->GetXaxis()->CenterTitle();
  h2p_test_old->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  h2p_test_old->GetYaxis()->CenterTitle();
  h2p_test_old->SetTitle("2D Binning TH1POLY function"); 
  h2p_test_old->Draw("text");
  
  DrawBinningInfo(TH1Poly_binMap_test);

 can -> Print(text_title_pdf2);
 
/////////////////////////////////////h2p_UB/


std::cout << "h2p_UB"<< std::endl;

Int_t nBins = h2p_UB->GetNumberOfBins();
for(Int_t bin = 1; bin< nBins; bin++){
std::cout<< "bin "<< bin << " content = " << h2p_UB->GetBinContent(bin)<< std::endl;

}

  h2p_UB->GetXaxis()->SetTitle("cos#theta_{#mu}");
  h2p_UB->GetXaxis()->CenterTitle();
  h2p_UB->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  h2p_UB->GetYaxis()->CenterTitle();
  h2p_UB->SetTitle("2D Binning h2p UB using Function"); 
  h2p_UB->Draw("text");
  DrawBinningInfo(TH1Poly_binMap_test);
 can -> Print(text_title_pdf2);
 
 std::cout<<"I think calling this function "<< std::endl;
 
  UBTH2Poly* h2p_UB_copy = h2p_UB->GetCopyWithBinNumbers("h2p_UB_copy");
  std::cout<<"Not "<< std::endl;
  h2p_UB_copy->GetXaxis()->SetTitle("cos#theta_{#mu}");
  h2p_UB_copy->GetXaxis()->CenterTitle();
  h2p_UB_copy->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  h2p_UB_copy->GetYaxis()->CenterTitle();
  h2p_UB_copy->SetTitle("2D Binning h2p_UB local COPY"); 
  h2p_UB_copy->Draw("COLZ");

 can -> Print(text_title_pdf2);
 auto sepvector = h2p_UB->GetNBinsX();
 std::cout<<"~~~~~~~~~~~~~~~~~~"<< std::endl;
std::cout<<"sepvector[0] = "<<sepvector<< std::endl;
/////////////////////////////////////
count =0.0;

 for (Int_t bin = 0; bin < h2p_test_inclusive->GetNumberOfBins(); ++bin) {
        // Get the bin's information
        Double_t x, y;
       // h2p->GetBin(bin, x, y);
        count=count+1.0;
        // Optionally, you can get more information about the bin, like the area
        //h2p->SetBinContent(bin,count);
        x = TH1Poly_binMap_test_inclusive1[bin].centerx;
        y = TH1Poly_binMap_test_inclusive1[bin].centery;

        
     
      h2p_test_inclusive->Fill(x, y, count);
     h2p_test_inclusive_UB->Fill(x, y, count);

        // You can perform additional operations on each bin here
 }

//int precision = 2;  // Change this value based on your requirement
//gStyle->SetPaintTextFormat(Form(".%df", precision));

  h2p_test_inclusive->GetXaxis()->SetTitle("cos#theta_{#mu}");
  h2p_test_inclusive->GetXaxis()->CenterTitle();
  h2p_test_inclusive->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
  h2p_test_inclusive->GetYaxis()->CenterTitle();
  h2p_test_inclusive->SetTitle("2D Binning Inclusive"); 
  h2p_test_inclusive->Draw("text");

DrawBinningInfo(TH1Poly_binMap_test_inclusive1);

 
can -> Print(text_title_pdf2);


h2p_test_inclusive_UB->GetXaxis()->SetTitle("cos#theta_{#mu}");
h2p_test_inclusive_UB->GetXaxis()->CenterTitle();
h2p_test_inclusive_UB->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
h2p_test_inclusive_UB->GetYaxis()->CenterTitle();
h2p_test_inclusive_UB->SetTitle("2D Binning Inclusive"); 
h2p_test_inclusive_UB->SetMarkerSize(0);
h2p_test_inclusive_UB->Draw("text");

DrawBinningInfo(TH1Poly_binMap_test_inclusive2);
DrawBinningNum(TH1Poly_binMap_test_inclusive2);



can -> Print(text_title_pdf2);

h2p_test_inclusive_UB->Draw("COLZ");
can -> Print(text_title_pdf2);
bool includeOverFlow= false; 
//////////////////////////////////////
///// testing projection method 
////////////////////////////////////
std::vector<Binning2D> BinningVector = GetProjectBinVector();
std::map<Binning2D, TH1D*> EffbinMap = ConstructProjectionMap(*TFile_MC,BinningVector, "EfficiencyProjections" );
double max =  1.15* MaxYofMap(EffbinMap);

auto binstringmap = Projection9Bins_StringMap(MUON_2D_BIN_EDGES, "P_{#mu} [GeV/c]");
std::vector<double> YMultipliers{0};
bool dontDo_bin_width_norm = false; 
bool Do_bin_width_norm = true;

PlotMC2D(
EffbinMap,
binstringmap,
text_title_pdf2, "", "cos#theta_{#mu}",
"Efficiency", "Efficiency",
.59, false,.15, false,
 YMultipliers,
 dontDo_bin_width_norm,
.018,
1.0);
//char *pdf_label, char *histotitle, char *xaxislabel,
//char *legendTitle, char *zaxislabel_units,

std::map<Binning2D, TH1D*> PuritybinMap = ConstructProjectionMap(*TFile_MC,BinningVector, "PurityProjections" );
//double max =  1.15* MaxYofMap( PuritybinMap);




PlotMC2D(
 PuritybinMap,
binstringmap,
text_title_pdf2, "", "cos#theta_{#mu}",
"Purity", "Purity",
1.3, false, .15, false,
 YMultipliers,
 dontDo_bin_width_norm,
.018,
1.0);


PlotMC2D_purity_Eff(
EffbinMap,
PuritybinMap,
binstringmap,
	text_title_pdf2, "", "cos#theta_{#mu}",
	"Efficiency*Purity",
	1.2, true, 0.15,  false,
	YMultipliers,
 dontDo_bin_width_norm,
	.03);



  auto effIter = EffbinMap.begin();
    auto purityIter = PuritybinMap.begin();

    // Iterate through both maps simultaneously
    while (effIter != EffbinMap.end() && purityIter != PuritybinMap.end()) {
        // Access and modify elements from EffbinMap
        auto& HistPurity = purityIter->second;
        auto& HistEff = effIter->second;  // Using a reference to modify the value directly

        // Perform an operation (e.g., double the value)
       HistEff->Multiply(HistPurity);

        // Access elements from PuritybinMap


        // Advance both iterators
        ++effIter;
        ++purityIter;
    }



PlotMC2D(
EffbinMap,
binstringmap,
text_title_pdf2, "" , "c",
"Eff*Purity", "Eff*Purity",
.55, false,.15, false,
 YMultipliers,
 dontDo_bin_width_norm,
.017,
1.0);




for(auto Projection: EffbinMap){
sprintf(inputName, "%s",binstringmap[Projection.first].c_str());
Projection.second->SetTitle(inputName);
Projection.second->SetMaximum(max);
Projection.second->SetMinimum(.15);
Projection.second->Draw("Hist");
 can -> Print(text_title_pdf2);
}


std::map<Binning2D, TH1D*> Eff_inclusivebinMap = ConstructProjectionMap(*TFile_MC,BinningVector, "EfficiencyProjections_Inclusive" );
//double max =  1.15* MaxYofMap(EffbinMap);
std::map<Binning2D, TH1D*> Purity_inclusivebinMap = ConstructProjectionMap(*TFile_MC,BinningVector, "PurityProjections_Inclusive" );

auto binstringmap_inclusive = Projection9Bins_StringMap(MUON_2D_BIN_EDGES_inclusive, "cos#theta_{#mu}");
 

PlotMC2D(
 Eff_inclusivebinMap,
 binstringmap_inclusive,
	text_title_pdf2, "", "P_{#mu} [GeV/c]",
	"Efficiency", "Efficiency",
	.7, false, 0.0, false,
	 YMultipliers,
	 dontDo_bin_width_norm,
	.017,
	1.0);


PlotMC2D(
 Purity_inclusivebinMap,
binstringmap_inclusive,
	text_title_pdf2, "", "P_{#mu} [GeV/c]",
	"Purity", "Purity",
	1.1, false,.15, false,
	 YMultipliers,
	 dontDo_bin_width_norm,
	.017,
	1.0);

PlotMC2D_purity_Eff(
Eff_inclusivebinMap,
Purity_inclusivebinMap,
binstringmap_inclusive,
	text_title_pdf2, "", "P_{#mu} [GeV/c]",
	"Efficiency*Purity",
	1.2, true, 0.1,  false,
	YMultipliers,
 dontDo_bin_width_norm,
	.03);







auto eff_inclusiveIter = Eff_inclusivebinMap.begin();
    auto purity_inclusiveIter = Purity_inclusivebinMap.begin();

    // Iterate through both maps simultaneously
    while (eff_inclusiveIter != Eff_inclusivebinMap.end() && purity_inclusiveIter != Purity_inclusivebinMap.end()) {
        // Access and modify elements from EffbinMap
        auto& HistPurity = purity_inclusiveIter->second;
        auto& HistEff = eff_inclusiveIter->second;  // Using a reference to modify the value directly

        // Perform an operation (e.g., double the value)
       HistEff->Multiply(HistPurity);

        // Access elements from PuritybinMap


        // Advance both iterators
        ++eff_inclusiveIter;
        ++purity_inclusiveIter;
    }


PlotMC2D(
 Eff_inclusivebinMap,
 binstringmap_inclusive,
	text_title_pdf2, "", "P_{#mu} [GeV/c]",
	"Eff*Purity", "Eff*Purity",
	.59, false, 0.0, false,
	 YMultipliers,
	 dontDo_bin_width_norm,
	.017,
	1.0);


std::cout<<"making new tfiles "<<std::endl;

sprintf(inputName, "%s",InputDataFile_name.c_str());
//TFile *TFile_Data = new TFile(inputName);
 //std::map<Binning2D, TH1D*> InputHistmap_DATA = ConstructProjectionMap(*TFile_Data,BinningVector, "h_costheta_Pmu_UBTH2Poly" );
     sprintf(inputName, "%s",InputData_BG_File_name.c_str());
//TFile *TFile_Data_beamOff = new TFile(inputName);
  //std::map<Binning2D, TH1D*> InputHistmap_DATA_BeamOff = ConstructProjectionMap(*TFile_Data_beamOff,BinningVector, "h_costheta_Pmu_UBTH2Poly" );
   std::map<Binning2D, TH1D*> InputHistmap_MC = ConstructProjectionMap(*TFile_MC,BinningVector, "h_costheta_Pmu_UBTH2Poly" );
   
   std::map<Binning2D, TH1D*> InputHistmap_MC_Contained = ConstructProjectionMap(*TFile_MC,BinningVector, "h_costheta_Pmu_UBTH2Poly_Contained" );
   std::map<Binning2D, TH1D*> InputHistmap_MC_NonContained = ConstructProjectionMap(*TFile_MC,BinningVector, "h_costheta_Pmu_UBTH2Poly_NonContained" );
      std::map<Binning2D, TH1D*> InputHistmap_MC_Signal = ConstructProjectionMap(*TFile_MC,BinningVector, "h_costheta_Pmu_UBTH2Poly_Signal" );

   auto& eci = EventCategoryInterpreter::Instance();
   auto stackvector = eci.ReturnCategoryVector();
   
   std::cout<<"Finshed new tfiles "<<std::endl;
   
std::map<std::pair<Binning2D, EventCategory>, TH1D*> SuperStack = Construct_CategoryProjectionMap(
*TFile_MC,
 BinningVector,
 stackvector, 
 "h_costheta_Pmu_UBTH2Poly");

std::map<std::pair<Binning2D, EventCategory>, TH1D*> SuperStack_Contained = Construct_CategoryProjectionMap(
*TFile_MC,
 BinningVector,
 stackvector, 
 "h_costheta_Pmu_UBTH2Poly_Contained");


std::map<std::pair<Binning2D, EventCategory>, TH1D*> SuperStack_NonContained = Construct_CategoryProjectionMap(
*TFile_MC,
 BinningVector,
 stackvector, 
 "h_costheta_Pmu_UBTH2Poly_NonContained");
 
 std::map<std::pair<Binning2D, EventCategory>, TH1D*> SuperStack_Signal = Construct_CategoryProjectionMap(
*TFile_MC,
 BinningVector,
 stackvector, 
 "h_costheta_Pmu_UBTH2Poly_Signal");

 std::cout<<"Finshed super stack  "<<std::endl;

std::vector<double> YMultipliers_1{10,4,4,1,1,1,1,5,10};

//PlotMC2D_Stack(
// InputHistmap_DATA,
// InputHistmap_DATA_BeamOff,
//  InputHistmap_MC,
//SuperStack,
//binstringmap,
//	text_title_pdf2, "",
//	"cos#theta_{#mu}","", "Events",
//	250, false, 0,  true,
//	YMultipliers_1,
//	 !Do_bin_width_norm,
//	.017,
//	 DATA_POT, 
// POT_MC_Scale,
// bd_scale);


PlotFractionofEvents2D_new(
   SuperStack,
   InputHistmap_MC,
   binstringmap,
   DATA_POT, 
   .025,
	 text_title_pdf2, 
	 " ", 
	 "cos#theta_{#mu}",
	 "Fraction of Events");

PlotFractionofEvents2D(
   SuperStack_Contained,
   InputHistmap_MC_Contained,
   binstringmap,
   DATA_POT, 
   .018,
	 text_title_pdf2, 
	 " Muon Contained ", 
	 "cos#theta_{#mu}",
	 "Fraction of Events");

PlotFractionofEvents2D(
   SuperStack_NonContained,
   InputHistmap_MC_NonContained,
   binstringmap,
   DATA_POT, 
   .018,
	 text_title_pdf2, 
	 " Muon  NonContained ", 
	 "cos#theta_{#mu}",
	 "Fraction of Events");



PlotFractionofEvents2D(
   SuperStack_Signal,
   InputHistmap_MC_Signal,
   binstringmap,
   DATA_POT, 
   .018,
	 text_title_pdf2, 
	 " ", 
	 "cos#theta_{#mu}",
	 "Fraction of Events");


//std::map<Binning2D, TH1D*> InputHistmap_DATA_inclusive = ConstructProjectionMap(*TFile_Data,BinningVector, "h_costheta_Pmu_UBTH2Poly_inclusive" );
//std::map<Binning2D, TH1D*> InputHistmap_DATA_BeamOff_inclusive = ConstructProjectionMap(*TFile_Data_beamOff,BinningVector, "h_costheta_Pmu_UBTH2Poly_inclusive" );
std::map<Binning2D, TH1D*> InputHistmap_MC_inclusive = ConstructProjectionMap(*TFile_MC,BinningVector, "h_costheta_Pmu_UBTH2Poly_inclusive" );

std::map<std::pair<Binning2D, EventCategory>, TH1D*> SuperStack_inclusive = Construct_CategoryProjectionMap(
*TFile_MC,
 BinningVector,
 stackvector, 
 "h_costheta_Pmu_UBTH2Poly_inclusive");


std::vector<double> YMultipliers_2{4,4,2,2,1,1,1,1,1};

//PlotMC2D_Stack(
// InputHistmap_DATA_inclusive,
// InputHistmap_DATA_BeamOff_inclusive,
//  InputHistmap_MC_inclusive,
//SuperStack_inclusive,
//binstringmap_inclusive,
//	text_title_pdf2, "",
//	"P_{#mu} [GeV/c]","", "Events",
//	275, false, 0,  true,
//	YMultipliers_2,
//	 !Do_bin_width_norm,
//	.017,
//	 DATA_POT, 
// POT_MC_Scale,
// bd_scale);

PlotFractionofEvents2D_new(
   SuperStack_inclusive,
   InputHistmap_MC_inclusive,
   binstringmap_inclusive,
   DATA_POT, 
   .025,
	 text_title_pdf2, 
	 " ", 
	 "P_{#mu} [GeV/c]",
	 "Fraction of Events");


/*
 std::map<Binning2D, TH1D*> InputHistmap_DATA,
 std::map<Binning2D, TH1D*> InputHistmap_DATA_BeamOff,
 std::map<Binning2D, TH1D*> InputHistmap_MC,
 std::map<std::pair<Binning2D, EventCategory>, TH1D*> StackMap_input,
 std::map<Binning2D , std::string > BinStringMap,
	char *pdf_label, char *histotitle, char *xaxislabel,
	char *legendTitle, char *zaxislabel_units,
	double Ymax, bool setMaxY, double Ymin,  bool doMultipliers,
	std::vector<double> YMultipliers,
	bool dontDo_bin_width_norm,
	double text_size,
	float POT_DATA, 
float POT_scaler_MC,
float BG_Trigger_scaler)
*/


 can -> Print(text_title_pdf3);

  
  
}/////End Of Function 
//////////////////////





void convertMapToArrays(const std::map<double, std::vector<double>>& inputMap, std::vector<Double_t*>& outputArrays) {
    // Iterate through the map
    for (const auto& entry : inputMap) {
        // Copy the vector values into a new array
        Double_t* newArray = new Double_t[entry.second.size() + 1];
        std::copy(entry.second.begin(), entry.second.end(), newArray);
        newArray[entry.second.size()] = {};  // Null-terminate the array

        // Add the new array to the output
        outputArrays.push_back(newArray);
    }
}


//////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////


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
  
  //TH1D* h_extBNB =(TH1D*)extBNB_input->Clone("h_extBNB");
  //EventInterp.set_ext_histogram_style( h_extBNB );
  
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
    slice_pred_stack,
    "NEvents / Bin Width",
    "",
    makeNormWidth, 
    "title",
    ymax,
    c1);


   //lg1_stacked->AddEntry(slice_ext->hist_.get(),"EXT BNB" , "f" );
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
  TLegend* lg1_stacked = new TLegend( 0.38, 0.63, 0.88, 0.87 );
  lg1_stacked->SetNColumns(2);
  lg1_stacked->SetBorderSize(0);
  Double_t defaultTextSize = lg1_stacked->GetTextSize();
  lg1_stacked->SetTextSize(.03); //
  
  const auto& DrawStackMCandData_EventCategory_tool = EventCategoryInterpreter::Instance();
  std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH1D(*TFile_MC, BaseStackName, inputEventCateory_vector );
  // Set Hist to Stevens Style 
  
  
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
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
  EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB


  for(auto &HistInMap : h_p_EventCategory_StackMap){
    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
    lg1_stacked->AddEntry(HistInMap.second, HistInMap.second->GetTitle() , "f" );
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
      EventCategory_stack,
      Yaxis_title,
      Xaxis_title,
     DoBinWidthNorm, 
      Title,
      YMax,
      Canvas);
  
  
  

   lg1_stacked->AddEntry(h_Data_BeamOFF,"EXT BNB" , "f" );
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
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
  EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB



   
   for (auto it = h_p_EventCategory_StackMap.rbegin(); it != h_p_EventCategory_StackMap.rend(); ++it) {
    // 'it' is a reverse iterator pointing to the current element in reverse order
    auto& HistInMap = *it;

    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
}
   
     for(auto &HistInMap : h_p_EventCategory_StackMap){
  
    //DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    //EventCategory_stack->Add( HistInMap.second );
    lg1_stacked->AddEntry(HistInMap.second, HistInMap.second->GetTitle() , "f" );
  
  }
   
   h_MC_Total->GetXaxis()->SetTitle(Xaxis_title);
   h_Data_BeamOn->GetXaxis()->SetTitle( Xaxis_title );
   h_Data_BeamOFF->GetXaxis()->SetTitle( Xaxis_title );
  
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
  //pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
  
  
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
  
  
   lg1_stacked->AddEntry(h_Data_BeamOFF,"EXT BNB" , "f" );
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

  //h_MC_Total->Add(h_Data_BeamOFF);
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
  //EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB


  for(auto &HistInMap : h_p_EventCategory_StackMap){
    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
    lg1_stacked->AddEntry(HistInMap.second, HistInMap.second->GetTitle() , "f" );
  }
   h_MC_Total->GetXaxis()->SetTitle(Xaxis_title);
   h_Data_BeamOn->GetXaxis()->SetTitle( Xaxis_title );
   h_Data_BeamOFF->GetXaxis()->SetTitle( Xaxis_title );
  
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
  //pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
  
  
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
  
  
   lg1_stacked->AddEntry(h_Data_BeamOFF,"EXT BNB" , "f" );
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
  std::map<Particle_type, TH1D*> h_p_EventCategory_StackMap =  MakeNamedCategoryMap_TH1D(*TFile_MC, BaseStackName, inputEventCateory_vector );
  // Set Hist to Stevens Style 
  
  
  TH1D* h_Data_BeamOn_1 =  GetTH1DHist(*TFile_Data, InputHisBaseName );
  TH1D* h_Data_BeamOFF_1 =  GetTH1DHist(*TFile_Data_BeamOff, InputHisBaseName );
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

  //h_MC_Total->Add(h_Data_BeamOFF);
  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }
  
  TH1D* h_Data_BeamOn_clone =(TH1D*)h_Data_BeamOn->Clone("h_Data_BeamOn_clone");
  TH1D* h_MC_Total_clone =(TH1D*)h_MC_Total->Clone("h_MC_Total_clone");
  
  
  
  int Nbins =  h_MC_Total->GetNbinsX();
  
  THStack* EventCategory_stack = new THStack( "mc+ext", "" );
  
  DrawStackMCandData_EventCategory_tool.set_ext_histogram_style(h_Data_BeamOFF  );
  DrawStackMCandData_EventCategory_tool.set_bnb_data_histogram_style( h_Data_BeamOn );
  
  lg1_stacked->AddEntry(h_Data_BeamOn, "Data", "pe" );
  lg1_stacked->AddEntry(h_MC_Total, "Total MC", "le" );
  //EventCategory_stack->Add( h_Data_BeamOFF ); // extBNB


  for(auto &HistInMap : h_p_EventCategory_StackMap){
    DrawStackMCandData_EventCategory_tool.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
    lg1_stacked->AddEntry(HistInMap.second, HistInMap.second->GetTitle() , "f" );
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
      EventCategory_stack,
      Yaxis_title,
      Xaxis_title,
     DoBinWidthNorm, 
      Title,
      YMax,
      Canvas);
  
  
   lg1_stacked->AddEntry(h_Data_BeamOFF,"EXT BNB" , "f" );
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


int main() {

 
std::cout<<" Testing Slice Plots "<< std::endl;

  //tutorial_slice_plots();
  MakePlots();
  return 0;
}

/*
void DrawBinningInfo(std::map<int , BinMap> TH1Poly_binMap_input){
 for(auto binmap :TH1Poly_binMap_input){
  //double centerx; 
  //double centery; 
    Double_t xLowEdge =binmap.second.xlow;
    Double_t xUpEdge = binmap.second.xhigh;
    Double_t yLowEdge = binmap.second.ylow;
    Double_t yUpEdge =binmap.second.yhigh;

    // Offset values for demonstration
    Double_t offset1 = 0.025;
    Double_t offset2 = -0.034;
    Double_t offset2y = -0.01;
    
     TText* textrange1 = new TText(binmap.second.centerx,
                             binmap.second.centery + offset2 + offset2y,
                             Form("X (%1.2f,%1.2f)",xLowEdge,xUpEdge));
            textrange1->SetTextSize(0.005);
            textrange1->SetTextAlign(22);  // Centered
            textrange1->SetTextColor(kBlue);
            textrange1->Draw(); 
            
            TText* textrange = new TText(binmap.second.centerx,
                                     binmap.second.centery + offset2 ,
                                     Form("(%1.2f,%1.2f)",yLowEdge,yUpEdge));
            textrange->SetTextSize(0.005);
            textrange->SetTextAlign(22);  // Centered
            textrange->SetTextColor(kBlue );
            textrange->Draw();             
 }// End of Loop 
  
 return; 

}


void DrawBinningNum(std::map<int , BinMap> TH1Poly_binMap_input){
 
 for(auto binmap :TH1Poly_binMap_input){
  //double centerx; 
  //double centery; 
    Double_t xLowEdge =binmap.second.xlow;
    Double_t xUpEdge = binmap.second.xhigh;
    Double_t yLowEdge = binmap.second.ylow;
    Double_t yUpEdge =binmap.second.yhigh;

    // Offset values for demonstration
    Double_t offset1 = 0.025;
    Double_t offset2 = -0.034;
    Double_t offset2y = -0.01;
    
     TText* textrange1 = new TText(binmap.second.centerx,
                             binmap.second.centery,
                             Form("%i",binmap.first));
            textrange1->SetTextSize(0.02);
            textrange1->SetTextAlign(22);  // Centered
            textrange1->SetTextColor(kRed);
            textrange1->Draw(); 
                       
 }// End of Loop 
  
 return; 

}
*/
/*
void AddandFillBinMap(std::map<int , BinMap> &InputBinMap, std::map<double, std::vector<double>>  input_BinMap ){

    for (auto it = input_BinMap.begin(); it !=input_BinMap.end(); ++it) {
        double currentKey = it->first;
        auto  Xvector = it->second;
        int vector_size = Xvector.size();   
        // Check if there is a next element
        auto nextIt = std::next(it);
        if (nextIt != input_BinMap.end()) {
            double nextKey = nextIt->first;

            for (int k = 0; k <vector_size-1; k++ )
            {
            
            double avg_X = .5 * (Xvector.at(k) +Xvector.at(k+1) );
            double avg_Y = .5 * (currentKey + nextKey); 
            
            BinMap BinMapSingle{currentKey,nextKey,Xvector.at(k),Xvector.at(k+1),avg_X,avg_Y};
            int binN = h2p->AddBin(Xvector.at(k),currentKey, Xvector.at(k+1),nextKey);
            Double_t x1 = Xvector.at(k);
            Double_t y1 = currentKey;
            Double_t x2 = Xvector.at(k+1);
            Double_t y2 = nextKey;
            
            
            h2p_UB->AddBin(x1,y1,x2,y2);
            
            
            //std::cout<<"BinN = "<< binN<< std::endl;
                TH1Poly_binMap.insert(std::pair<int,BinMap >(binN,BinMapSingle)); 
            }
            //std::cout<< " "<< std::endl;
            
            
        }
        else {
            // If there is no next element, handle accordingly
           // std::cout << "Finished last bin  element with Key=" << currentKey << std::endl;
        }
    
    //std::cout<<"~~~~~~~~~~~~~~~~~"<< std::endl;
    
    
 }


}// End of function 
/////////////////////


void FillBinnMap_inclusive(std::map<int , BinMap> &InputBinMap, std::map<double, std::vector<double>>  input_BinMap ){



 for (auto it = input_BinMap.begin(); it != input_BinMap.end(); ++it) {
        double currentKey = it->first;
        auto  Xvector = it->second;
        int vector_size = Xvector.size();   
        // Check if there is a next element
        auto nextIt = std::next(it);
        if (nextIt != input_BinMap.end()) {
            double nextKey = nextIt->first;

            for (int k = 0; k <vector_size-1; k++ )
            {
            //std::cout<<Xvector.at(k)<< ", "<< Xvector.at(k+1)<< std::endl;
            double avg_X = .5 * (Xvector.at(k) +Xvector.at(k+1) );
            double avg_Y = .5 * (currentKey + nextKey); 
            
            BinMap BinMapSingle{currentKey,nextKey,Xvector.at(k),Xvector.at(k+1),avg_X,avg_Y};
            int binN = h2p_inclusive->AddBin(currentKey,Xvector.at(k),nextKey, Xvector.at(k+1));
            //std::cout<<"BinN = "<< binN<< std::endl;
                InputBinMap.insert(std::pair<int,BinMap >(binN,BinMapSingle)); 
            }
            //std::cout<< " "<< std::endl;
            
            
        } else {
            // If there is no next element, handle accordingly
            //std::cout << "Finished last bin  element with Key=" << currentKey << std::endl;
        }
    
    //
    
    
   // std::cout<<"~~~~~~~~~~~~~~~~~"<< std::endl;
    
    
}



}// End of function 
/////////////////////
*/
/*
void PlotMC2D_Stack(
 std::map<Binning2D, TH1D*> InputHistmap_DATA,
 std::map<Binning2D, TH1D*> InputHistmap_DATA_BeamOff,
 std::map<Binning2D, TH1D*> InputHistmap_MC,
 std::map<std::pair<Binning2D, EventCategory>, TH1D*> StackMap_input,
 std::map<Binning2D , std::string > BinStringMap,
 char *pdf_label, char *histotitle, char *xaxislabel,
 char *legendTitle, char *zaxislabel_units,
 double Ymax, bool setMaxY, double Ymin,  bool doMultipliers,
 std::vector<double> YMultipliers,
 bool dontDo_bin_width_norm,
 double text_size,
 float POT_DATA, 
float POT_scaler_MC,
float BG_Trigger_scaler)
{

 std::cout<<"inside:PlotMC2D_Stack "<< std::endl;
	auto BinVector = GetProjectBinVector();
  auto& CategoryInterpreter = EventCategoryInterpreter::Instance();
  std::vector<EventCategory>  Category_vector =  CategoryInterpreter.ReturnCategoryVector();
	int nbins_X = InputHistmap_DATA.size()-1;
	int nbins_Y = InputHistmap_DATA.size()-1;
	//std::cout << "The NEntries for TObjArray: " << nVars << std::endl;

	 TLegend *leg = new TLegend(0.41, 0.1, 0.96, 0.25);
    std::string Leg_title_string = get_legend_title(POT_DATA);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.02);
	leg->SetNColumns(3);
  leg->SetHeader(Leg_title_string.c_str());
	

	double min_XAxis = InputHistmap_DATA[BinVector.at(0)]->GetXaxis()->GetXmin();
	double max_XAxis = InputHistmap_DATA[BinVector.at(0)]->GetXaxis()->GetBinUpEdge(InputHistmap_DATA[BinVector.at(0)]->GetNbinsX());
 std::cout<< "max_XAxis = "<< max_XAxis<< std::endl;
 double maxmax = -1; 

	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	grid_x = sqrt(nbins_X + 1);
	grid_y = nbins_X / (grid_x - 1);
	Nbins = InputHistmap_DATA.size();

	 std::cout << nbins_X - grid_x * grid_y << std::endl;
	 std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");
	gc->SetRightMargin(0.06);
	gc->SetLeftMargin(0.06);
	gc->ResetPads();
	
	std::cout<<"starting Loop "<< std::endl;

	for (int i = 1; i <= Nbins; ++i)
	{
	
	   std::cout<<"Bin index: "<< i << " BinVector.at(i-1) = "<< BinVector.at(i-1) << std::endl; 
	
		TH1D *h_Data_BeamOn = (TH1D *)InputHistmap_DATA[BinVector.at(i-1)]->Clone(uniq());
		std::cout<<" 1"<< std::endl; 
		TH1D *h_Data_BeamOFF = (TH1D *)InputHistmap_DATA_BeamOff[BinVector.at(i-1)]->Clone(uniq());
		std::cout<<"2 "<< std::endl; 
		TH1D *h_MC_Total = (TH1D *)InputHistmap_MC[BinVector.at(i-1)]->Clone(uniq());
		
		 std::cout<<"Made clones "<< std::endl; 
		
		h_MC_Total->Scale(POT_scaler_MC); 
    h_Data_BeamOFF->Scale(BG_Trigger_scaler); 
		
		 THStack* EventCategory_stack = new THStack( "mc+ext", "" );
		
		 
		
		  CategoryInterpreter.set_ext_histogram_style(h_Data_BeamOFF);
      CategoryInterpreter.set_bnb_data_histogram_style( h_Data_BeamOn );
		
		h_MC_Total->Add(h_Data_BeamOFF);
		
		std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap;
		//std::cout<<" about to do this loop : Category_vector "<< std::endl; 
		for(auto cat: Category_vector){
		  	//std::cout<<" Getting (,) =  ("<< BinVector.at(i-1)<< " , "<< cat<<std::endl; 
		std::string leg_stack_check = CategoryInterpreter.label( cat );
     	//std::cout<<" Getting (,) =  ("<< BinVector.at(i-1)<< " , "<< cat<<std::endl; 
       TH1D* hist = (TH1D*) StackMap_input[{BinVector.at(i-1),cat}]->Clone(uniq());
       //hist->SetTitle(map.second->GetTitle());
       h_p_EventCategory_StackMap.insert(std::pair<EventCategory, TH1D*>(cat,hist)); 
     }
		
	  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }

		
		/////////////////
		/// Want Beam Off First 
		////////////////
		EventCategory_stack->Add( h_Data_BeamOFF ); 
		
	for (auto it = h_p_EventCategory_StackMap.rbegin(); it != h_p_EventCategory_StackMap.rend(); ++it) {
    // 'it' is a reverse iterator pointing to the current element in reverse order
    auto& HistInMap = *it;

    CategoryInterpreter.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
 }
		
		
		 if(i ==1){
		for(auto &HistInMap : h_p_EventCategory_StackMap){
     std::string leg_stack_label = CategoryInterpreter.label( HistInMap.first );
    leg->AddEntry(HistInMap.second, leg_stack_label.c_str() , "f" );
    }
    
  }
		
		
		if (i == 1) {
		leg->AddEntry(h_Data_BeamOn, "Data", "pe" );
    leg->AddEntry(h_MC_Total, "Total MC", "le" );
		leg->AddEntry(h_Data_BeamOFF, "bnbext", "f");
		
		}
		
		
		
	//////////////////////////////	
		gc->cd(i);
	/////////////////////////////


		//if (dontDo_bin_width_norm == true) {
		//h_Data->Scale(1., "width");
		//h_DataBeamOff->Scale(1., "width");
		//h_mc->Scale(1., "width");
		//}
		
		if (doMultipliers == true) {
		h_Data_BeamOn ->Scale(YMultipliers.at(i - 1));
		scaleTHStack(EventCategory_stack, YMultipliers.at(i - 1));
		//EventCategory_stack->Scale(YMultipliers.at(i - 1));
		
		h_MC_Total->Scale(YMultipliers.at(i - 1));
		}
		
		DrawStackMCandData(
    h_Data_BeamOn, 
    h_MC_Total,
    h_Data_BeamOFF,
    	EventCategory_stack,
    dontDo_bin_width_norm, 
    Ymax);
   
   
   
   
   drawString(BinStringMap[BinVector.at(i-1)], text_size, false );
		
		if (doMultipliers && YMultipliers.at(i - 1) != 1)
		{
			auto pad = gc->cd(i);
			TLatex *la2 = new TLatex(1 - pad->GetRightMargin() - 0.02,
				1 - pad->GetTopMargin() - 0.05,
				TString::Format("#times %.1f", YMultipliers.at(i - 1)));
			la2->SetTextAlign(33);	// top right
			la2->SetNDC();
			la2->SetTextFont(42);
			la2->SetTextSize(0.03);
			la2->Draw();
		}

 if (i == 1){DrawFakeData();}


		//mcprojyArr->Clear();
	}	/// End of Loop

    std::cout<<"Finished 2D drawing loop "<< std::endl;

	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.015);

	
	double MaxY = gc->GetPadMax();
	std::cout<<"Got Maxium Pad"<< std::endl;
	if (setMaxY == false) gc->SetYLimits(0, MaxY *1.35);
	else
	{
		gc->SetYLimits(0, Ymax);
	}
	
	
	gc->SetInterpadSpace(.005);
	//gc->SetRightMargin(0.05);
	//gc->SetLeftMargin(0.08);
  gc->SetYLimits(Ymin, Ymax);
  gc->SetXLimits(min_XAxis,max_XAxis);
  gc->ResetPads();
	gc->SetXTitle(xaxislabel);
	gc->SetYTitleSize(25);
	gc->SetXTitleSize(20);  
	gc->SetYTitle(zaxislabel_units);
	gc->SetTitleAlignmentFor6Hist();
	leg->Draw("SAME");
	gc->Modified();
	gc->Print(pdf_label);
	//leg->Clear();
	delete gc;

}
*/
////////////////////////////////////////////////////////////////////////////////
/*
void scaleTHStack(THStack* stack, double scaleFactor) {
    if (!stack) {
        std::cerr << "Error: Invalid THStack pointer." << std::endl;
        return;
    }

    // Get the list of histograms from the THStack
    TList* histograms = stack->GetHists();

    if (!histograms) {
        std::cerr << "Error: No histograms found in THStack." << std::endl;
        return;
    }

    // Iterate through the histograms and scale each one
    TIter next(histograms);
    TObject* obj;

    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom("TH1")) {
            TH1* hist = dynamic_cast<TH1*>(obj);
            hist->Scale(scaleFactor);
        }
    }
}
*/
/*
void DrawFakeData(){
    TLatex* textfake = new TLatex;
    textfake->SetNDC();
    textfake->SetTextSize(0.015);
    textfake->SetTextColor(kBlue);
    textfake->DrawLatex(0.14, 0.86, "USING FAKE DATA (NuWro)");

}
*/


  void PlotFractionofEvents2D(
   std::map<std::pair<Binning2D, EventCategory>, TH1D*> StackMap_input,
   std::map<Binning2D, TH1D*> InputHistmap_MC,
   std::map<Binning2D , std::string > BinStringMap,
   double POT_DATA,
   double text_size,
	 char *pdf_label, 
	 char *histotitle, 
	 char *xaxislabel,
	 char *zaxislabel)
{
	
	auto BinVector = GetProjectBinVector();
  auto& CategoryInterpreter = EventCategoryInterpreter::Instance();
  std::vector<EventCategory>  Category_vector =  CategoryInterpreter.ReturnCategoryVector();
  double CompleteTotal = 0.0; 
  
	for(auto hist: InputHistmap_MC){CompleteTotal = hist.second->Integral() + CompleteTotal;}
		
	int nbins_X = InputHistmap_MC.size()-1;
	int nbins_Y = InputHistmap_MC.size()-1;
	
	double min_XAxis = InputHistmap_MC[BinVector.at(0)]->GetXaxis()->GetXmin();
	double max_XAxis = InputHistmap_MC[BinVector.at(0)]->GetXaxis()->GetBinUpEdge(InputHistmap_MC[BinVector.at(0)]->GetNbinsX());

	 TLegend *leg = new TLegend(0.41, 0.1, 0.96, 0.25);
    //std::string Leg_title_string = get_legend_title(POT_DATA);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.02);
	leg->SetNColumns(3);
  leg->SetHeader("MicroBooNE Simulation");
	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	//   grid_x = sqrt(nbins_Y)+1;
	// grid_y = nbins_Y/(grid_x-1);
	// Nbins = dataHist->GetNbinsY();

	grid_x = sqrt(nbins_X) + 1;
	grid_y = nbins_X / (grid_x - 1);
	Nbins = InputHistmap_MC.size();

	//if (DeBug == true) std::cout << nbins_X - grid_x * grid_y << std::endl;
	//ZZif (DeBug == true) std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	//std::vector<int> bins;	// = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");
	//gc->SetRightMargin(0.06);
	//gc->SetLeftMargin(0.06);
	gc->ResetPads();

	for (int i = 1; i <= Nbins; ++i)
	{
    auto BIN_index = BinVector.at(i-1);
		TObjArray *mcprojyArr = new TObjArray();
		
		
	 TH1D *h_MC_Total = (TH1D *)InputHistmap_MC[BIN_index]->Clone(uniq());
   TH1D *h_MC_Total_clone = (TH1D *)h_MC_Total->Clone(uniq());

		double area = h_MC_Total->Integral();
		h_MC_Total->Scale(1.0 / area);
		h_MC_Total_clone->Scale(1.0 / CompleteTotal);

		h_MC_Total->SetMaximum(1.25);
		h_MC_Total->SetMarkerStyle(20);
		h_MC_Total->SetMarkerSize(.5);
		h_MC_Total_clone->SetMarkerStyle(21);
		h_MC_Total_clone->SetMarkerSize(.3);
		
		if (i == 1)
		{
			leg->AddEntry(h_MC_Total, "Bin Fraction", "p");
			leg->AddEntry(h_MC_Total_clone, "Total Fraction", "p");
		}

 ////////////////////////////////////////////

		
		for(auto Category: Category_vector ){
		
		TH1D* mcHist_Category = (TH1D*)StackMap_input[{BIN_index,Category}]->Clone(uniq());
		CategoryInterpreter.set_mc_histogram_style(Category, mcHist_Category );
		
			if (i == 1)
			{
			  std::string words = CategoryInterpreter.label(Category);
				char words_char[words.length() + 1];
				strcpy(words_char, words.c_str());
				leg->AddEntry(mcHist_Category, words_char, "f");
			}
		
				mcprojyArr->Add(mcHist_Category);
		
		}
 ///////////////////////////////////////////
		gc->cd(i);

		BinNormalizeTOFractionOF_Events(*mcprojyArr);
		
		THStack* EventCategory_stack = new THStack( "mc+ext", "" );
		//int i = 0; i < mcprojyArr->GetEntries(); ++i
		 for (int i = mcprojyArr->GetEntries() - 1; i >= 0; --i) {
        TObject* obj = mcprojyArr->At(i);

        if (obj && obj->IsA()->InheritsFrom("TH1D")) {
            TH1D* hist = dynamic_cast<TH1D*>(obj);
            EventCategory_stack->Add(hist);
        }
    }
		
		

   DrawStack(
    h_MC_Total, 
    EventCategory_stack,
    false, 
    "",
    1.25);
    
		h_MC_Total_clone->Draw("SAME E1 X0");
   drawString(BinStringMap[BIN_index], text_size, false );

	}	/// End of Loop
   //////////////////////////////////////////

	gc->SetInterpadSpace(.005);
	gc->SetXLimits(min_XAxis, max_XAxis);
	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.015);
	gc->SetYLimits(0, 1.25);
	gc->SetXTitle(xaxislabel);
	gc->SetYTitleSize(25);
	gc->SetXTitleSize(20);  
	gc->SetYTitle(zaxislabel);
	gc->SetTitleAlignmentFor6Hist();
	gc->SetHistTexts();
	leg->Draw("SAME");
	gc->Modified();
	gc->Print(pdf_label);

	delete gc;

}
/////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////
  void PlotFractionofEvents2D_new(
   std::map<std::pair<Binning2D, EventCategory>, TH1D*> StackMap_input,
   std::map<Binning2D, TH1D*> InputHistmap_MC,
   std::map<Binning2D , std::string > BinStringMap,
   double POT_DATA,
   double text_size,
	 char *pdf_label, 
	 char *histotitle, 
	 char *xaxislabel,
	 char *zaxislabel)
{
	
	std::map<EventCategory, double> Combined_SignalType_stats; 
	auto BinVector = GetProjectBinVector();
  auto& CategoryInterpreter = EventCategoryInterpreter::Instance();
  std::vector<EventCategory>  Category_vector =  CategoryInterpreter.ReturnCategoryVector();
  double CompleteTotal = 0.0; 
  
	for(auto hist: InputHistmap_MC){
	double amount_Panel = hist.second->Integral();
	CompleteTotal = amount_Panel + CompleteTotal;
	}
	
	for(auto cateogry:Category_vector ){Combined_SignalType_stats.insert(std::pair<EventCategory,double>(cateogry,0.0));}
	
	
	for(const auto& entry : StackMap_input){
	      const auto& binningEventPair = entry.first;
        const Binning2D& binning = binningEventPair.first;
        const EventCategory& eventCategory = binningEventPair.second;

        // Extract TH1D pointer
        TH1D* hist = entry.second;
        double area = hist->Integral();
        Combined_SignalType_stats[eventCategory] = Combined_SignalType_stats[eventCategory] + area; 
	}
	
	
		
	int nbins_X = InputHistmap_MC.size()-1;
	int nbins_Y = InputHistmap_MC.size()-1;
	
	double min_XAxis = InputHistmap_MC[BinVector.at(0)]->GetXaxis()->GetXmin();
	double max_XAxis = InputHistmap_MC[BinVector.at(0)]->GetXaxis()->GetBinUpEdge(InputHistmap_MC[BinVector.at(0)]->GetNbinsX());

	 TLegend *leg = new TLegend(0.25, 0.01, 0.96, 0.21);
    //std::string Leg_title_string = get_legend_title(POT_DATA);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetNColumns(3);
  leg->SetHeader("MicroBooNE Simulation");
	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	//   grid_x = sqrt(nbins_Y)+1;
	// grid_y = nbins_Y/(grid_x-1);
	// Nbins = dataHist->GetNbinsY();

	grid_x = sqrt(nbins_X) + 1;
	grid_y = nbins_X / (grid_x - 1);
	Nbins = InputHistmap_MC.size();

	//if (DeBug == true) std::cout << nbins_X - grid_x * grid_y << std::endl;
	//ZZif (DeBug == true) std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	//std::vector<int> bins;	// = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");
	gc->SetBottomMargin(.00);
  gc->SetTopMargin(.02);
  gc->SetRightMargin(.05);
	//gc->SetRightMargin(0.06);
	//gc->SetLeftMargin(0.06);
	gc->ResetPads();



	for (int i = 1; i <= Nbins; ++i)
	{
    auto BIN_index = BinVector.at(i-1);
		TObjArray *mcprojyArr = new TObjArray();
		
		
	 TH1D *h_MC_Total = (TH1D *)InputHistmap_MC[BIN_index]->Clone(uniq());
   TH1D *h_MC_Total_clone = (TH1D *)h_MC_Total->Clone(uniq());

		double area = h_MC_Total->Integral();
		h_MC_Total->Scale(1.0 / area);
		h_MC_Total_clone->Scale(1.0 / CompleteTotal);

		h_MC_Total->SetMaximum(1.25);
		h_MC_Total->SetMarkerStyle(20);
		h_MC_Total->SetMarkerSize(.8);
		h_MC_Total_clone->SetMarkerStyle(21);
		h_MC_Total_clone->SetMarkerSize(.5);
		
		if (i == 1)
		{
			leg->AddEntry(h_MC_Total, "Statical Panel Fraction", "p");
			//leg->AddEntry(h_MC_Total_clone, "Total Fraction", "p");
		}

 ////////////////////////////////////////////

		
		for(auto Category: Category_vector ){
		
		TH1D* mcHist_Category = (TH1D*)StackMap_input[{BIN_index,Category}]->Clone(uniq());
		CategoryInterpreter.set_mc_histogram_style(Category, mcHist_Category );
		
			if (i == 1)
			{
			
			double totalEventAmount = Combined_SignalType_stats[Category];
			double percentage = (totalEventAmount/CompleteTotal) * 100; 
			
			  std::string words = CategoryInterpreter.label(Category) + " (" + to_string_with_precision(percentage, 1) + " %) ";
				char words_char[words.length() + 1];
				strcpy(words_char, words.c_str());
				leg->AddEntry(mcHist_Category, words_char, "f");
			}
		
				mcprojyArr->Add(mcHist_Category);
		
		}
 ///////////////////////////////////////////
		gc->cd(i);

		BinNormalizeTOFractionOF_Events(*mcprojyArr);
		
		THStack* EventCategory_stack = new THStack( "mc+ext", "" );
		//int i = 0; i < mcprojyArr->GetEntries(); ++i
		 for (int i = mcprojyArr->GetEntries() - 1; i >= 0; --i) {
        TObject* obj = mcprojyArr->At(i);

        if (obj && obj->IsA()->InheritsFrom("TH1D")) {
            TH1D* hist = dynamic_cast<TH1D*>(obj);
            EventCategory_stack->Add(hist);
        }
    }
		
		

   DrawStack(
    h_MC_Total, 
    EventCategory_stack,
    false, 
    "",
    1.25);
    
	//h_MC_Total_clone->Draw("SAME E1 X0");
   drawString(BinStringMap[BIN_index], text_size, false );

   float fraction = (area / CompleteTotal)*100;
   std::string Ratio = "[" +    to_string_with_precision(fraction, 1)  + " %]";
   drawString_extra(Ratio, .025, kRed, true, .01, .02 );
   

	}	/// End of Loop
   //////////////////////////////////////////

	gc->SetInterpadSpace(.005);
	gc->SetXLimits(min_XAxis, max_XAxis);
	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.015);
	gc->SetYLimits(0, 1.25);
	gc->SetXTitle(xaxislabel);
	gc->SetYTitleSize(25);
	gc->SetXTitleSize(20);  
	gc->SetYTitle(zaxislabel);
	gc->SetTitleAlignmentFor6Hist();
	gc->SetHistTexts();
	leg->Draw("SAME");
	gc->Modified();
	gc->Print(pdf_label);

	delete gc;

}

void drawString_extra(std::string inputString, double text_size, int color, bool left, double add_Xpostion, double add_Ypostion )
{
 
  TLatex* la=0;
  TString text(inputString);
  
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+add_Xpostion,
                  1-gPad->GetTopMargin()-add_Ypostion,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-add_Xpostion,
                  1-gPad->GetTopMargin()-add_Ypostion,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(text_size);
  la->SetTextColor(color);
  la->Draw();
}
