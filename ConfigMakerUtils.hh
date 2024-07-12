#pragma once
// Provides variable definitions for possible use in multiple configuration
// file writer scripts

// Standard library includes
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

// STV analysis includes
//#include "HistUtils.hh"
#include "includes/SliceBinning.hh"
#include "includes/HistUtils_cc0pi.hh"
// Helper function used to get the indices for the variables used to define a
// SliceBinning object
int find_slice_var_index( const std::string& name,
  const std::vector< SliceVariable >& sv_vec )
{
  constexpr int BOGUS_INDEX = -1;
  auto iter = std::find_if( sv_vec.cbegin(), sv_vec.cend(),
    [ &name ]( const SliceVariable& svar )
      -> bool { return name == svar.name_; }
  );
  if ( iter == sv_vec.cend() ) return BOGUS_INDEX;
  else return std::distance( sv_vec.cbegin(), iter );
}

// Helper function for adding a new slice to a SliceBinning object. This
// is useful when defining slices simultaneously with a binning scheme.
Slice& add_slice( SliceBinning& sb, const std::vector< double >& bin_edges_1d,
  int active_var_idx, int other_var_idx = -1, double other_low = DBL_MAX,
  double other_high = DBL_MAX )
{
  sb.slices_.emplace_back();
  auto& new_slice = sb.slices_.back();

  // Create the slice histogram
  int num_bins = bin_edges_1d.size() - 1;
  TH1D* slice_hist = new TH1D( "slice_hist", ";;events", num_bins,
    bin_edges_1d.data() );
  slice_hist->SetDirectory( nullptr );
  new_slice.hist_.reset( slice_hist );

  // Also set up the slice variable definitions
  new_slice.active_var_indices_.push_back( active_var_idx );

  if ( other_var_idx >= 0 ) {
    new_slice.other_vars_.emplace_back( other_var_idx, other_low, other_high );
  }

  return new_slice;
}
////////////////////////////////////////////////////////////////////////////////////
std::vector<double> get_bin_low_edges_local( double xmin, double xmax, int Nbins )
{
  std::vector<double> bin_low_edges;
  double bin_step = ( xmax - xmin ) / Nbins;
  for ( int b = 0; b <= Nbins; ++b ) {
    double low_edge = xmin + b*bin_step;
    bin_low_edges.push_back( low_edge );
  }

  return bin_low_edges;
}//////////////////////////////// End of Function 
////////////////////////////////////////////////////////////////////////////////////

Slice& add_slice( SliceBinning& sb, int num_bins, double active_low,
  double active_high, int active_var_idx, int other_var_idx = -1,
  double other_low = DBL_MAX, double other_high = DBL_MAX )
{
  auto edges = get_bin_low_edges_local( active_low, active_high, num_bins );

  return add_slice( sb, edges, active_var_idx, other_var_idx, other_low,
    other_high );
}




// These strings provide non-default selections (in TTree::Draw format) for
// sideband control samples.

const std::string DIRT_SIDEBAND_SELECTION =
  "!sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_muon_passed_mom_cuts"
  " && !sel_has_pion_candidate";

/*
const std::string NC_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && !sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers";
  
  const std::string NC_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && !sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_has_p_candidate"
  " && sel_passed_proton_pid_cut && sel_protons_contained"
  " && sel_lead_p_passed_mom_cuts";

*/

const std::string NC_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && !sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_has_p_candidate"
  " && sel_protons_contained"
  " && sel_lead_p_passed_mom_cuts"
   "&& sel_num_proton_candidates > 0";


const std::string CCNPI_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_muon_passed_mom_cuts"
  " && sel_has_pion_candidate"
  " && sel_num_pion_candidates > 0";

/*
const std::string CCNP_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_muon_passed_mom_cuts"
  " && sel_muon_contained && sel_muon_quality_ok"
  " && sel_has_p_candidate && sel_protons_contained"
  " && sel_lead_p_passed_mom_cuts"
  " && BTD_Muon_Candidate_prediction != 9999 && !TMath::IsNaN(BTD_Muon_Candidate_prediction) && BTD_Muon_Candidate_prediction > .14 ";
*/


//const std::string BDT_Bogus_SIDEBAND_SELECTION =
//  "sel_CC0pi && !sel_BDT_NotBOGUS  && sel_BDT_predicts_1plusTrks_tobeProtons";
//
//
//const std::string BDT_Else_SIDEBAND_SELECTION =
//"sel_CC0pi && sel_BDT_NotBOGUS  && !sel_BDT_predicts_1plusTrks_tobeProtons";
//
//
//const std::string BDT_Bogus_Else_SIDEBAND_SELECTION =
//"sel_CC0pi && !sel_BDT_NotBOGUS  && !sel_BDT_predicts_1plusTrks_tobeProtons";


// Bin definitions for the 2D cross-section measurements

  // Using floating-point numbers as std::map keys is admittedly evil, but
  // it's safe in this case: all we'll do with this map is iterate over the
  // elements. Keys are muon momentum bin edges, values are muon scattering
  // cosine bin edges.
  /* 
  // No need for an underflow bin: due to the signal definition, all muons
    // with reco momentum below 0.1 GeV/c will be lost
    
    //Steven type binning 
    
    
    { 0.1, { -1, -0.55, -0.25, 0., 0.25, 0.45, 0.7, 1.00 }, },
    { 0.24, { -1, -0.55, -0.25, 0., 0.25, 0.45, 0.7, 1.00 } },
    { 0.3,  { -1, -0.4, -0.1, 0.1, 0.35, 0.5, 0.7, 0.85, 1. } },
    { 0.38, { -1, 0, 0.5, 0.65, 0.8, 0.92, 1.00 } },
    { 0.48, { -1, 0.2, 0.5, 0.65, 0.8, 0.875, 0.950, 1.00 } },
    { 0.7, { -1, 0.65, 0.8, 0.875, 0.950, 1.00 } },
    { 0.85, { -1, 0.85, 0.9, 0.950, 1.00 } },
 1.28, 2.50
    // Upper edge of the last bin. Due to the signal definition, no overflow
    // bin is needed for muons above 1.2 GeV/c
    { 1.2, {} }
2 GeV/c
    { 1.2, {} }


cc-incluive binning 

{ −1.00, {0.00, 0.18, 0.30, 0.45, 0.77, 2.50} },
{ −0.50, {0.00, 0.18, 0.30, 0.45, 0.77, 2.50} },
 { 0.00, {0.00, 0.18, 0.30, 0.45, 0.77, 2.50} },
 { 0.27, {0.00, 0.30, 0.45, 0.77, 2.50} },
 { 0.45, {0.00, 0.30, 0.45, 0.77, 2.50} },
 { 0.62, {0.00, 0.30, 0.45, 0.77, 2.50} },
 { 0.76, {0.00, 0.30, 0.45, 0.77, 1.28, 2.50} },
 { 0.86, {0.00, 0.30, 0.45, 0.77, 1.28, 2.50} },
 { 0.94, {0.00, 0.30, 0.45, 0.77, 1.28, 2.50} },
 { 1.0, {} }


// Took Steves binning added the 1.585
 
 OLD Binning
  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES = {
    { 0.1, {  -1, -0.55, -0.25,    0., 0.25,  0.45,  0.7, 1.00 } },
    { 0.24, { -1, -0.55, -0.25,    0., 0.25,  0.45,  0.7, 1.00 } },
    { 0.3,  { -1,  -0.4,  -0.1,   0.1, 0.35,  0.5,   0.7, 0.85, 1. } },
    { 0.38, { -1,     0,   0.5,  0.65, 0.8,   0.92,  1.00 } },
    { 0.48, { -1,   0.2,   0.5,  0.65, 0.8,   0.875, 0.950, 1.00 } },
    { 0.7, {  -1,  0.65,   0.8, 0.875, 0.950, 1.00 } },
    { 0.85, { -1,  0.85,   0.9, 0.950, 1.00} },
    { 1.28, { -1,  0.85,   0.9, 0.950, 1.00} },
    { 1.58, { -1,  0.85,   0.9, 0.950, 1.00} },
    { 2.0, {} }

  };



*/
  
  
  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES = {
    { 0.1, {  -1, -0.55,    0.,  0.45, 1.00 } },
    { 0.24, { -1, -0.55, -0.25,    0., 0.25,  0.45,  0.7, 1.00 } },
    { 0.3,  { -1,  -0.4,  -0.1,   0.1, 0.35,  0.5,   0.7, 0.85, 1. } },
    { 0.38, { -1,     0,   0.5,  0.65, 0.8,   0.92,  1.00 } },
    { 0.48, { -1,   0.2,   0.5,  0.65, 0.8,   0.875, 0.950, 1.00 } },
    { 0.7,  { -1,  0.65,   0.8, 0.875, 0.950, 1.00 } },
    { 0.85, { -1,  0.85,   0.9, 0.950, 1.00} },
    { 1.28, { -1,  0.85,  1.00} },
    { 1.58, { -1,  0.85,  1.00} },
    { 2.0, {} }

  };



  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES_for1DSidBand = {
    { 0.1, {-1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1. } },
    { 0.275, {-1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1. } },
    { 0.375,  {-1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1.  } },
    { 0.525, { -1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1. } },
    { 0.48, { -1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1. } },
    { 0.7, {  -1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1. } },
    { 0.9, { -1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1.} },
    { 1.1, {-1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1. } },
    { 2.0, {-1., -0.525, -0.375, -0.225, -0.1, 0.025, 0.15, 0.275, 0.45, 0.6, 0.775, 0.925, 1.} }

  };

  /*
  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES_inclusive = {
// the 2D binning of inclusive but took the max limit of 2 GeV and not 2.5 GeV which its given  


{ -1.00, {0.00, 0.18, 0.30, 0.45, 0.77, 2.0} },
{ -0.50, {0.00, 0.18, 0.30, 0.45, 0.77, 2.0} },
 { 0.00, {0.00, 0.18, 0.30, 0.45, 0.77, 2.0} },
 { 0.27, {0.00, 0.30, 0.45, 0.77, 2.0} },
 { 0.45, {0.00, 0.30, 0.45, 0.77, 2.0} },
 { 0.62, {0.00, 0.30, 0.45, 0.77, 2.0} },
 { 0.76, {0.00, 0.30, 0.45, 0.77, 1.28, 2.0} },
 { 0.86, {0.00, 0.30, 0.45, 0.77, 1.28, 2.0} },
 { 0.94, {0.00, 0.30, 0.45, 0.77, 1.28, 2.0} },
 { 1.0, {} }

  };
*/

  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES_inclusive = {
// the 2D binning of inclusive but took the max limit of 2 GeV and not 2.5 GeV which its given  
{ -1.00, 
     {0.1, 0.24, 0.38, 2.0} },
{ -0.50,
    {0.1, 0.24, 0.3, 0.38, 0.48, 2.0} },
 { 0.00, 
   {0.1, 0.24, 0.3, 0.38, 0.48, 2.0} },
 { 0.27, 
     {0.1, 0.3, 0.48, 0.85, 2.0} },
 { 0.45, 
     {0.1, 0.38, 0.48, .68, 0.85, 2.0} },
 { 0.62, 
     {0.1, 0.38, 0.48, 0.68, 0.85, 1.28, 2.0} },
 { 0.76, 
     {0.1, 0.3, 0.48, 0.68, 0.85, 2.0} },
 { 0.86, 
     {0.1, 0.48, 0.68, 0.85, 1.28,  2.0} },
 { 0.94, 
     {0.1, 0.38, 0.48, 0.68, 0.85, 1.28, 2.0} },
 { 1.0, {} }
};




/*

  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES_inclusive = {
// the 2D binning of inclusive but took the max limit of 2 GeV and not 2.5 GeV which its given  
{ -1.00, 
     {0.1, 0.24, 0.38, 2.0} },
{ -0.50,
    {0.1, 0.24, 0.3, 0.38, 0.48, 2.0} },
 { 0.00, 
   {0.1, 0.24, 0.3, 0.38, 0.48, 2.0} },
 { 0.27, 
     {0.1, 0.24, 0.3, 0.48, 0.85, 2.0} },
 { 0.45, 
     {0.1, 0.24, 0.3, 0.38, 0.48, .68, 0.85, 2.0} },
 { 0.62, 
     {0.1, 0.24, 0.3, 0.38, 0.48, 0.68, 0.85, 1.28, 2.0} },
 { 0.76, 
     {0.1, 0.3, 0.48, 0.68, 0.85, 1.28, 1.58, 2.0} },
 { 0.86, 
     {0.1, 0.3, 0.48, 0.68, 0.85, 1.28,  2.0} },
 { 0.94, 
     {0.1, 0.38, 0.48, 0.68, 0.85, 1.28, 2.0} },
 { 1.0, {} }
};
*/
 


  // Using floating-point numbers as std::map keys is admittedly evil, but it's
  // safe in this case: all we'll do with this map is iterate over the
  // elements. Keys are proton momentum bin edges, values are proton cosine bin
  // edges.
  std::map< double, std::vector<double> > PROTON_2D_BIN_EDGES = {

    // No need for an underflow bin: due to the signal definition, all leading
    // protons with reco momentum below 0.25 GeV/c will be lost
    { 0.250, { -1, 0., 1.0 } },
    { 0.325, { -1, -0.5, 0, 0.5, 0.8, 1.0 } },
    { 0.4,   { -1, -0.6, -0.2, 0.2, 0.5, 0.65, 0.85, 1.0 } },
    { 0.5, { -1, -0.2, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.6,   { -1, 0.1, 0.37, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.7,   { -1, 0.45, 0.65, 0.75, 0.82, 0.9, 1.0 } },

    // Upper edge of the last bin. We don't need an overflow bin because the
    // signal definition excludes any leading protons with momenta above
    // 1 GeV/c
    { 1., {} }

  };
  
  
  
  

///////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> MakeKeyVectorFromMap(std::map< double, std::vector<double> > inputMap){


std::vector<double> return_vector; 

 for (auto it = inputMap.begin(); it != inputMap.end(); ++it) {
     return_vector.push_back(it->first); 
    }

   return return_vector; 

}
///////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<double>> MakeVectorofVectorsFromMap(std::map< double, std::vector<double> > inputMap){
  
  
  std::vector<std::vector<double>> return_vector;
  
   for (auto it = inputMap.begin(); it != inputMap.end(); ++it) {
     return_vector.push_back(it->second); 
    }
  
  return return_vector; 
  
  
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////
  
  int NBinsFromMap(std::map< double, std::vector<double> > inputMap){
  int binN = 1;
  
     for (auto it = inputMap.begin(); it != inputMap.end(); ++it) {
     binN = (it->second.size() -1) + binN; 
    }
  return binN; 
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void make_config_mcc9_2D_muon_CC0Pi();
  void make_config_mcc9_2D_muon_CC0Pi_inclusive();
  std::map<int, double> readTextFileToMap(const char* fileName);
  std::map<int, double> readTextFileToMap_2(const char* fileName);
  