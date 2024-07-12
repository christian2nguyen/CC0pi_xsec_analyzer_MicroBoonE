#include "ConfigMakerUtils.hh"
#include "UniverseMaker.hh"
#include "includes/HistUtils_cc0pi.hh"
///////////////////////////////
// redefined for cc0pi
//////////////////////////////
std::vector<double> generateDoubleVector(int start, int end);
std::vector<double> generateDoubleVectorWithIncrement(double start, double end, double increment);
void make_config_RECOCuting();


void make_config_mcc9_2D_muon_CC0Pi() {

  std::cout<<"Inside: :make_config_mcc9_2D_muon_CC0Pi"<< std::endl;
  std::string selection = "sel_CC0pi";
  std::string signal_def = "mc_is_cc0pi_signal";
//&& sel_BDT_predicts_1plusTrks_tobeProtons
  // By construction, MC event categories 5-11 contain all beam-correlated
  // backgrounds. This list is therefore comprehensive apart from cosmic
  // overlay stuff which is directly measured in a dedicated sample.
  std::vector< std::string > background_defs = {
    "category == 0","category == 5", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };


  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;
  
  std::map<int, double> Eff_BinMap = readTextFileToMap("Binning_Eff.txt");
  std::map<int, double> FakeBinRate_BinMap = readTextFileToMap_2("Num_RECO_Entries.txt");

 std::cout<<"Eff_BinMap.size() = "<< Eff_BinMap.size()<< std::endl;
 for(auto bin : Eff_BinMap){
 std::cout<<bin.first<< " , "<<bin.second << std::endl;
 }
std::cout<<"fake Data"<< std::endl;

 for(auto bin : FakeBinRate_BinMap){
 std::cout<<bin.first<< " , "<<bin.second << std::endl;
 }

  // Also save the bin definitions to an output file that will be used to
  // make a LaTeX table
  std::ofstream tex_bin_table_file( "mybintable_mcc9_muon2D_CC0pi_v1_july10_2024.tex" );
  tex_bin_table_file << "\\documentclass{standalone}\n"
    << "\\usepackage{booktabs}\n"
    << "\\usepackage{siunitx}\n"
    << "\\DeclareSIUnit\\clight{\\text{\\ensuremath{c}}}\n"
    << "\\sisetup{per-mode=symbol}\n"
    << "\\begin{document}\n"
    << "\\begin{tabular}{cSSSScc}\n"
    << "\\toprule\n"
    << "bin number\n"
    << "& {$p_\\mu^\\mathrm{low}$ (\\si{\\GeV\\per\\clight})}"
    << "& {$p_\\mu^\\mathrm{high}$ (\\si{\\GeV\\per\\clight})}"
    << " & {$\\cos\\theta_\\mu^\\mathrm{low}$}\n"
    << "& {$\\cos\\theta_\\mu^\\mathrm{high}$} & efficiency & occupancy \\\\\n"
    << "\\midrule\n";

  // Configure kinematic limits for all of the signal bins

  // The reco bins are numbered in the order that their edges appear in the
  // map, so just keep a running counter here to keep track of which reco
  // bin we are on.
  size_t cur_reco_bin = 1u;

  // Get an iterator to the last map element. They are sorted numerically,
  // so this will be the upper edge of the last non-overflow bin.
  auto last = MUON_2D_BIN_EDGES.cend();
  --last;

  auto iter = MUON_2D_BIN_EDGES.cbegin();
  for ( auto iter = MUON_2D_BIN_EDGES.cbegin(); iter != last; ++iter ) {

    // Get an iterator to the map element after the current one. Due to
    // the automatic sorting, this is guaranteed to contain the upper edge
    // of the current muon momentum bin.
    auto next = iter;
    ++next;

    // Get the current muon momentum bin limits
    double pmu_low = iter->first;
    double pmu_high = next->first;

    // Now iterate over the scattering cosine bins associated with the
    // current momentum bin. Note that we will skip any situations in
    // which the binning is undefined (i.e., because there are less than
    // two bin edges given)
    const auto& cosine_bin_edges = iter->second;

    size_t num_cosine_edges = cosine_bin_edges.size();
    size_t num_cosine_bins = 0u;
    if ( num_cosine_edges >= 2u ) num_cosine_bins = num_cosine_edges - 1u;

    for ( size_t b = 0u; b < num_cosine_bins; ++b ) {

      double cosmu_low = cosine_bin_edges.at( b );
      double cosmu_high = cosine_bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def
        << " && mc_p3_mu.Mag() >= " << pmu_low
        << " && mc_p3_mu.Mag() < " << pmu_high
        << " && mc_p3_mu.CosTheta() >= " << cosmu_low
        << " && mc_p3_mu.CosTheta() < " << cosmu_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, 0 );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && p3_mu.Mag() >= " << pmu_low
        << " && p3_mu.Mag() < " << pmu_high
        << " && p3_mu.CosTheta() >= " << cosmu_low
        << " && p3_mu.CosTheta() < " << cosmu_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, 0 );

      tex_bin_table_file << cur_reco_bin << " & ";
      
      if ( b == 0u ) {
        tex_bin_table_file << pmu_low << " & " << pmu_high << " & ";
      }
      else tex_bin_table_file << " & & ";

      tex_bin_table_file << cosmu_low << " & " << cosmu_high <<  " & "<< Eff_BinMap[cur_reco_bin] << " & " << .24*FakeBinRate_BinMap[cur_reco_bin]
        << " \\\\"; 
      // Add extra space at the end of each momentum bin
      if ( b == num_cosine_bins - 1 ) tex_bin_table_file << "[2mm]";
      tex_bin_table_file << '\n';
      ++cur_reco_bin;
      // We don't need an overflow cosine bin because the entire angular
      // range is covered. We'll use a single bin for the overflow in pmu.

    } // loop over scattering cosine bins

  } // loop over muon momentum bins

  // No overflow bin is needed due to the muon momentum upper limit in the
  // signal definition

  // Add true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin );
  }

  // We're done with all the "ordinary" bins. Add some extra reco bins to use
  // in the sideband control samples.

  // The control samples have significantly lower statistics, so bin in 1D
  // momentum slices only. Use the muon candidate momentum where possible, but
  // for NC switch to the leading proton candidate momentum. This is done
  // because the NC sideband selection excludes events in which a muon candidate
  // has been identified.
  //
  // Keys are selections to use for sidebands, values are the branch names for
  // the reconstructed momentum to use in each case. This is a pretty hacky way
  // to organize the information, but it is simple.
  std::map< std::string, std::string > sideband_selection_to_momentum_map = {
   /* {  DIRT_SIDEBAND_SELECTION, "p3_mu"     },*/
    {    NC_SIDEBAND_SELECTION, "p3_lead_p" },
    { CCNPI_SIDEBAND_SELECTION, "p3_mu"     },
  };

  // Loop over the sideband selection definitions. Prepare new reco bin
  // definitions for each in the appropriate 1D reconstructed momentum space.
  for ( const auto& sel_mom_pair : sideband_selection_to_momentum_map ) {
    const auto& side_sel = sel_mom_pair.first;
    const auto& mom_branch = sel_mom_pair.second;
    std::map< double, std::vector<double> >* bin_edge_map = nullptr;
    if ( mom_branch == "p3_mu" ) bin_edge_map = &MUON_2D_BIN_EDGES;
    else if ( mom_branch == "p3_lead_p" ) bin_edge_map = &PROTON_2D_BIN_EDGES;
    else throw std::runtime_error( "Unimplemented sideband momentum!" );

    // Get an iterator to the last map element. They are sorted numerically,
    // so this will be the upper edge of the last non-overflow momentum bin.
    auto last = bin_edge_map->cend();
    --last;

    for ( auto iter = bin_edge_map->cbegin(); iter != last; ++iter ) {

      // Get an iterator to the map element after the current one. Due to the
      // automatic sorting, this is guaranteed to contain the upper edge of the
      // current momentum bin.
      auto next = iter;
      ++next;

      // Get the current momentum bin limits
      double p_low = iter->first;
      double p_high = next->first;

      std::stringstream reco_ss;
      reco_ss << side_sel
        << " && " << mom_branch << ".Mag() >= " << p_low
        << " && " << mom_branch << ".Mag() < " << p_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );
    } // 1D reco momentum bins

    // For sidebands that use bins of muon candidate momentum, create the
    // overflow bin. This isn't needed for the proton momentum due to the
    // upper limit imposed in the signal definition.
    if ( mom_branch != "p3_mu" ) continue;

    double pmu_overflow_min = bin_edge_map->crbegin()->first;

    std::stringstream reco_ss;
    reco_ss << side_sel << " && " << mom_branch << ".Mag() >= "
      << pmu_overflow_min;

    std::string reco_bin_def = reco_ss.str();
    reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );

  } // sideband selection definitions

  // Dump this information to the output file
  std::ofstream out_file( "myconfig_mcc9_2D_muon_v1_july10_2024.txt" );
  out_file << "Muon2D\n";
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_mcc9_2D_muon_v1_july10_2024.txt" );
  sb_file << "3\n";
  sb_file << "\"reco p_{#mu}\" \"GeV/c\" \"reco $p_{\\mu}$\" \"GeV$/c$\"\n";
  sb_file << "\"reco cos#theta_{#mu}\" \"\" \"reco $\\cos\\theta_{\\mu}$\""
    " \"\"\n";
  sb_file << "\"reco bin number\" \"\" \"reco bin number\" \"\"\n";
  // Includes a slice for the overflow bin and three extra slices. One
  // for everything in terms of reco bin number, one integrated over angles,
  // and one showing the sideband control sample results (in terms of reco
  // bin number).
  size_t num_slices = MUON_2D_BIN_EDGES.size() + 3;
  sb_file << num_slices << '\n';

  // Get an iterator to the final entry in the edge map (this is the
  // upper edge of the last bin)
  auto last_edge = MUON_2D_BIN_EDGES.cend();
  --last_edge;

  // The reco bins are numbered in the order that their edges appear in the
  // map, so just keep a running counter here to keep track of which reco
  // bin we are on.
  cur_reco_bin = 0u;

  for ( auto iter = MUON_2D_BIN_EDGES.cbegin(); iter != last_edge; ++iter ) {
    // Each 1D slice uses the same y-axis units (reco events)
    sb_file << "\"events\"\n";
    const auto& edges = iter->second;
    int num_edges = edges.size();
    int num_bins = num_edges - 1;
    // The muon cosine is the sole "active variable" in each slice
    sb_file << "1 1 " << num_edges;
    for ( const auto& edge : edges ) {
      sb_file << ' ' << edge;
    }
    // The muon momentum is the sole "other variable" in each slice
    double pmu_low = iter->first;
    auto next = iter;
    ++next;
    double pmu_high = next->first;
    sb_file << "\n1 0 " << pmu_low << ' ' << pmu_high << '\n';
    sb_file << num_bins;

    for ( int b = 0; b < num_bins; ++b ) {
      int root_bin_idx = b + 1;
      sb_file << '\n' << cur_reco_bin << " 1 " << root_bin_idx;
      ++cur_reco_bin;
    } // cosine bins
    sb_file << '\n';

  } // pmu slices

  // Handle the overflow bin separately
  sb_file << "\"events\"\n"; // y-axis label
  // Still treat the muon cosine as the single active variable in a single bin
  // spanning the entire range
  sb_file << "1 1 2 -1.0 1.0\n";
  // The muon momentum is the sole "other" variable. This is an overflow bin,
  // which we signal by making the lower and upper edges equal
  double last_pmu_edge_value = last_edge->first;
  sb_file << "1 0 " << last_pmu_edge_value
    << ' ' << last_pmu_edge_value << '\n';
  // A single UniverseMaker reco bin contributes to the sole ROOT bin
  // in this histogram
  sb_file << "1\n" << cur_reco_bin << " 1 1";

  // Make a final slice with everything expressed in terms of reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Acount for the zero-based UniverseMaker bin indices
  size_t num_reco_bins = cur_reco_bin + 1;
  // There is one more edge than the number of bins
  sb_file << num_reco_bins + 1;
  for ( size_t e = 0u; e <= num_reco_bins; ++e ) {
    sb_file << ' ' << e;
  }
  sb_file << '\n';
  // For the "overall slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker bin and assign it to the matching
  // ROOT histogram bin
  sb_file << num_reco_bins << '\n';
  for ( size_t b = 0u; b < num_reco_bins; ++b ) {
    sb_file << b << " 1 " << b + 1 << '\n';
  }

  // Make a 1D slice in which the angles have been integrated out. This is a
  // measurement of the leading muon momentum distribution. For this slice,
  // we're still working in terms of reco event counts
  sb_file << "\"events\"\n";
  int num_pmu_edges = MUON_2D_BIN_EDGES.size();
  int num_pmu_bins = num_pmu_edges - 1;
  // The muon momentum is the sole "active variable" in each slice
  sb_file << "1 0 " << num_pmu_edges;
  for ( const auto& pmu_edge_pair : MUON_2D_BIN_EDGES ) {
    const auto pmu_edge = pmu_edge_pair.first;
    sb_file << ' ' << pmu_edge;
  }
  // There is no "other" variable for this slice since we've integrated out
  // the angular information
  sb_file << "\n0\n";
  // Now we're ready to build the 1D muon momentum bins from the 2D ones. We
  // need one entry in the list per reco bin (apart from the overflow bin),
  // although multiple reco bins will contribute to each slice bin in this
  // case.
  sb_file << num_reco_bins - 1;

  // Iterate through the 2D reco bins, noting that they are numbered in the
  // order that their edges appear in the map. In this case, all angular reco
  // bins with the same reco muon momentum should contribute to a particular
  // slice p_mu bin.
  cur_reco_bin = 0u;

  // Keep track of the ROOT slice bin index (one-based) with this counter
  int cur_slice_bin_idx = 1;

  for ( auto iter = MUON_2D_BIN_EDGES.cbegin(); iter != last_edge; ++iter ) {
    const auto& angle_bin_edges = iter->second;
    int num_angle_bins = angle_bin_edges.size() - 1;
    for ( int b = 0; b < num_angle_bins; ++b ) {
      sb_file << '\n' << cur_reco_bin << " 1 " << cur_slice_bin_idx;
      ++cur_reco_bin;
    } // cosine bins

    // Move to the next muon momentum bin in the slice
    ++cur_slice_bin_idx;
  } // pmu slices

  sb_file << '\n';

  // Make a slice containing the sideband results organized by reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Count the number of sideband reco bins. Also find the index of the
  // first one.
  size_t num_sideband_reco_bins = 0u;
  size_t first_sideband_bin_idx = 0u;
  bool found_first_sideband_bin = false;
  
  for ( size_t b = 0u; b < reco_bins.size(); ++b ) {
    const auto& rbin = reco_bins.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) {
      ++num_sideband_reco_bins;
      if ( !found_first_sideband_bin ) {
        found_first_sideband_bin = true;
        first_sideband_bin_idx = b;
      }
    }
  }
  // There is one more edge than the number of sideband bins
  sb_file << num_sideband_reco_bins + 1;
  for ( size_t e = 0u; e <= num_sideband_reco_bins; ++e ) {
    sb_file << ' ' << e + first_sideband_bin_idx;
  }
  sb_file << '\n';
  // For the "sideband slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker sideband bin and assign it to the
  // matching ROOT histogram bin
  sb_file << num_sideband_reco_bins << '\n';
  for ( size_t b = 0u; b < num_sideband_reco_bins; ++b ) {
    sb_file << b + first_sideband_bin_idx << " 1 " << b + 1 << '\n';
  }






}
void make_config_mcc9_2D_muon_CC0Pi_inclusive() {

  std::string selection = "sel_CC0pi ";
  std::string signal_def = "mc_is_cc0pi_signal";

  // By construction, MC event categories 5-11 contain all beam-correlated
  // backgrounds. This list is therefore comprehensive apart from cosmic
  // overlay stuff which is directly measured in a dedicated sample.
  std::vector< std::string > background_defs = {
    "category == 0", "category == 5", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };

  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;

  std::map<int, double> Eff_BinMap = readTextFileToMap("Binning_Eff_inclusive.txt");
  std::map<int, double> FakeBinRate_BinMap = readTextFileToMap_2("Num_RECO_Entries_Inclusive.txt");

  // Also save the bin definitions to an output file that will be used to
  // make a LaTeX table
  std::ofstream tex_bin_table_file( "mybintable_mcc9_muon2D_CC0pi_inclusive_v1_july10_2024.tex" );
  tex_bin_table_file << "\\documentclass{standalone}\n"
    << "\\usepackage{booktabs}\n"
    << "\\usepackage{siunitx}\n"
    << "\\DeclareSIUnit\\clight{\\text{\\ensuremath{c}}}\n"
    << "\\sisetup{per-mode=symbol}\n"
    << "\\begin{document}\n"
    << "\\begin{tabular}{cSSSScc}\n"
    << "\\toprule\n"
    << "bin number\n"//p_\\mu^
    << "& {$\\cos\\theta_\\mu^\\mathrm{low}$ (\\si{\\rad\\per\\clight})}"
    << "& {$\\cos\\theta_\\mu^\\mathrm{high}$ (\\si{\\rad\\per\\clight})}"
    << " & {$p_\\mu^\\mathrm{low}$}\n"
    << "& {$\\p_\\mu^\\mathrm{high}$} & efficiency & occupancy \\\\\n"
    << "\\midrule\n";

  // Configure kinematic limits for all of the signal bins

  // The reco bins are numbered in the order that their edges appear in the
  // map, so just keep a running counter here to keep track of which reco
  // bin we are on.
  size_t cur_reco_bin = 0u;

  // Get an iterator to the last map element. They are sorted numerically,
  // so this will be the upper edge of the last non-overflow bin.
  auto last = MUON_2D_BIN_EDGES_inclusive.cend();
  --last;

  auto iter = MUON_2D_BIN_EDGES_inclusive.cbegin();
  for ( auto iter = MUON_2D_BIN_EDGES_inclusive.cbegin(); iter != last; ++iter ) {

    // Get an iterator to the map element after the current one. Due to
    // the automatic sorting, this is guaranteed to contain the upper edge
    // of the current muon momentum bin.
    auto next = iter;
    ++next;

    // Get the current muon momentum bin limits
    double pmu_low = iter->first;
    double pmu_high = next->first;

    // Now iterate over the scattering cosine bins associated with the
    // current momentum bin. Note that we will skip any situations in
    // which the binning is undefined (i.e., because there are less than
    // two bin edges given)
    const auto& cosine_bin_edges = iter->second;

    size_t num_cosine_edges = cosine_bin_edges.size();
    size_t num_cosine_bins = 0u;
    if ( num_cosine_edges >= 2u ) num_cosine_bins = num_cosine_edges - 1u;

    for ( size_t b = 0u; b < num_cosine_bins; ++b ) {

      double cosmu_low = cosine_bin_edges.at( b );
      double cosmu_high = cosine_bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def
        << " && mc_p3_mu.CosTheta() >= " << pmu_low
        << " && mc_p3_mu.CosTheta() < " << pmu_high
        << " && mc_p3_mu.Mag() >= " << cosmu_low
        << " && mc_p3_mu.Mag() < " << cosmu_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, 0 );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && p3_mu.CosTheta() >= " << pmu_low
        << " && p3_mu.CosTheta() < " << pmu_high
        << " && p3_mu.Mag() >= " << cosmu_low
        << " && p3_mu.Mag() < " << cosmu_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, 0 );

      tex_bin_table_file << cur_reco_bin << " & ";
      ++cur_reco_bin;
      if ( b == 0u ) {
        tex_bin_table_file << pmu_low << " & " << pmu_high << " & ";
      }
      else tex_bin_table_file << " & & ";

      tex_bin_table_file << cosmu_low << " & " << cosmu_high << " & "<< Eff_BinMap[cur_reco_bin] << " & " << .24*FakeBinRate_BinMap[cur_reco_bin]
        << " \\\\";
      // Add extra space at the end of each momentum bin
      if ( b == num_cosine_bins - 1 ) tex_bin_table_file << "[2mm]";
      tex_bin_table_file << '\n';

      // We don't need an overflow cosine bin because the entire angular
      // range is covered. We'll use a single bin for the overflow in pmu.

    } // loop over scattering cosine bins

  } // loop over muon momentum bins

  // No overflow bin is needed due to the muon momentum upper limit in the
  // signal definition

  // Add true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin );
  }

  // We're done with all the "ordinary" bins. Add some extra reco bins to use
  // in the sideband control samples.

  // The control samples have significantly lower statistics, so bin in 1D
  // momentum slices only. Use the muon candidate momentum where possible, but
  // for NC switch to the leading proton candidate momentum. This is done
  // because the NC sideband selection excludes events in which a muon candidate
  // has been identified.
  //
  // Keys are selections to use for sidebands, values are the branch names for
  // the reconstructed momentum to use in each case. This is a pretty hacky way
  // to organize the information, but it is simple.
  std::map< std::string, std::string > sideband_selection_to_momentum_map = {
    {    NC_SIDEBAND_SELECTION, "p3_lead_p" },
    { CCNPI_SIDEBAND_SELECTION, "p3_mu"    },
  };

  // Loop over the sideband selection definitions. Prepare new reco bin
  // definitions for each in the appropriate 1D reconstructed momentum space.
  for ( const auto& sel_mom_pair : sideband_selection_to_momentum_map ) {
    const auto& side_sel = sel_mom_pair.first;
    const auto& mom_branch = sel_mom_pair.second;
    std::map< double, std::vector<double> >* bin_edge_map = nullptr;
    if ( mom_branch == "p3_mu" ) bin_edge_map = &MUON_2D_BIN_EDGES_inclusive;
    else if ( mom_branch == "p3_lead_p" ) bin_edge_map = &PROTON_2D_BIN_EDGES;
    else throw std::runtime_error( "Unimplemented sideband momentum!" );

    // Get an iterator to the last map element. They are sorted numerically,
    // so this will be the upper edge of the last non-overflow momentum bin.
    auto last = bin_edge_map->cend();
    --last;

    for ( auto iter = bin_edge_map->cbegin(); iter != last; ++iter ) {

      // Get an iterator to the map element after the current one. Due to the
      // automatic sorting, this is guaranteed to contain the upper edge of the
      // current momentum bin.
      auto next = iter;
      ++next;

      // Get the current momentum bin limits
      double p_low = iter->first;
      double p_high = next->first;

      std::stringstream reco_ss;
      reco_ss << side_sel
        << " && " << mom_branch << ".CosTheta() >= " << p_low
        << " && " << mom_branch << ".CosTheta() < " << p_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );
    } // 1D reco momentum bins

    // For sidebands that use bins of muon candidate momentum, create the
    // overflow bin. This isn't needed for the proton momentum due to the
    // upper limit imposed in the signal definition.
    if ( mom_branch != "p3_mu" ) continue;

    double pmu_overflow_min = bin_edge_map->crbegin()->first;

    std::stringstream reco_ss;
    reco_ss << side_sel << " && " << mom_branch << ".CosTheta() >= "
      << pmu_overflow_min;

    std::string reco_bin_def = reco_ss.str();
    reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );

  } // sideband selection definitions

  // Dump this information to the output file
  std::ofstream out_file( "myconfig_mcc9_2D_muon_inclusive_v1_july10_2024.txt" );
  out_file << "Muon2D\n";
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_mcc9_2D_muon_inclusive_v1_july10_2024.txt" );
  sb_file << "3\n";
  sb_file << "\"reco cos#theta_{#mu}\" \"\" \"reco $\\cos\\theta_{\\mu}$\"\n";
  sb_file << "\"reco p_{#mu}\" \"GeV/c\" \"reco $p_{\\mu}$\" \"GeV$/c$\"\""
    " \"\"\n";
  sb_file << "\"reco bin number\" \"\" \"reco bin number\" \"\"\n";
  // Includes a slice for the overflow bin and three extra slices. One
  // for everything in terms of reco bin number, one integrated over angles,
  // and one showing the sideband control sample results (in terms of reco
  // bin number).
  size_t num_slices = MUON_2D_BIN_EDGES_inclusive.size() + 3;
  sb_file << num_slices << '\n';

  // Get an iterator to the final entry in the edge map (this is the
  // upper edge of the last bin)
  auto last_edge = MUON_2D_BIN_EDGES_inclusive.cend();
  --last_edge;

  // The reco bins are numbered in the order that their edges appear in the
  // map, so just keep a running counter here to keep track of which reco
  // bin we are on.
  cur_reco_bin = 0u;

  for ( auto iter = MUON_2D_BIN_EDGES_inclusive.cbegin(); iter != last_edge; ++iter ) {
    // Each 1D slice uses the same y-axis units (reco events)
    sb_file << "\"events\"\n";
    const auto& edges = iter->second;
    int num_edges = edges.size();
    int num_bins = num_edges - 1;
    // The muon cosine is the sole "active variable" in each slice
    sb_file << "1 1 " << num_edges;
    for ( const auto& edge : edges ) {
      sb_file << ' ' << edge;
    }
    // The muon momentum is the sole "other variable" in each slice
    double pmu_low = iter->first;
    auto next = iter;
    ++next;
    double pmu_high = next->first;
    sb_file << "\n1 0 " << pmu_low << ' ' << pmu_high << '\n';
    sb_file << num_bins;

    for ( int b = 0; b < num_bins; ++b ) {
      int root_bin_idx = b + 1;
      sb_file << '\n' << cur_reco_bin << " 1 " << root_bin_idx;
      ++cur_reco_bin;
    } // cosine bins
    sb_file << '\n';

  } // pmu slices

  // Handle the overflow bin separately
  sb_file << "\"events\"\n"; // y-axis label
  // Still treat the muon cosine as the single active variable in a single bin
  // spanning the entire range
  sb_file << "1 1 2 .1 2.0\n";
  // The muon momentum is the sole "other" variable. This is an overflow bin,
  // which we signal by making the lower and upper edges equal
  double last_pmu_edge_value = last_edge->first;
  sb_file << "1 0 " << last_pmu_edge_value
    << ' ' << last_pmu_edge_value << '\n';
  // A single UniverseMaker reco bin contributes to the sole ROOT bin
  // in this histogram
  sb_file << "1\n" << cur_reco_bin << " 1 1";

  // Make a final slice with everything expressed in terms of reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Acount for the zero-based UniverseMaker bin indices
  size_t num_reco_bins = cur_reco_bin + 1;
  // There is one more edge than the number of bins
  sb_file << num_reco_bins + 1;
  for ( size_t e = 0u; e <= num_reco_bins; ++e ) {
    sb_file << ' ' << e;
  }
  sb_file << '\n';
  // For the "overall slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker bin and assign it to the matching
  // ROOT histogram bin
  sb_file << num_reco_bins << '\n';
  for ( size_t b = 0u; b < num_reco_bins; ++b ) {
    sb_file << b << " 1 " << b + 1 << '\n';
  }

  // Make a 1D slice in which the angles have been integrated out. This is a
  // measurement of the leading muon momentum distribution. For this slice,
  // we're still working in terms of reco event counts
  sb_file << "\"events\"\n";
  int num_pmu_edges = MUON_2D_BIN_EDGES_inclusive.size();
  int num_pmu_bins = num_pmu_edges - 1;
  // The muon momentum is the sole "active variable" in each slice
  sb_file << "1 0 " << num_pmu_edges;
  for ( const auto& pmu_edge_pair : MUON_2D_BIN_EDGES_inclusive ) {
    const auto pmu_edge = pmu_edge_pair.first;
    sb_file << ' ' << pmu_edge;
  }
  // There is no "other" variable for this slice since we've integrated out
  // the angular information
  sb_file << "\n0\n";
  // Now we're ready to build the 1D muon momentum bins from the 2D ones. We
  // need one entry in the list per reco bin (apart from the overflow bin),
  // although multiple reco bins will contribute to each slice bin in this
  // case.
  sb_file << num_reco_bins - 1;

  // Iterate through the 2D reco bins, noting that they are numbered in the
  // order that their edges appear in the map. In this case, all angular reco
  // bins with the same reco muon momentum should contribute to a particular
  // slice p_mu bin.
  cur_reco_bin = 0u;

  // Keep track of the ROOT slice bin index (one-based) with this counter
  int cur_slice_bin_idx = 1;

  for ( auto iter = MUON_2D_BIN_EDGES_inclusive.cbegin(); iter != last_edge; ++iter ) {
    const auto& angle_bin_edges = iter->second;
    int num_angle_bins = angle_bin_edges.size() - 1;
    for ( int b = 0; b < num_angle_bins; ++b ) {
      sb_file << '\n' << cur_reco_bin << " 1 " << cur_slice_bin_idx;
      ++cur_reco_bin;
    } // cosine bins

    // Move to the next muon momentum bin in the slice
    ++cur_slice_bin_idx;
  } // pmu slices

  sb_file << '\n';

  // Make a slice containing the sideband results organized by reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Count the number of sideband reco bins. Also find the index of the
  // first one.
  size_t num_sideband_reco_bins = 0u;
  size_t first_sideband_bin_idx = 0u;
  bool found_first_sideband_bin = false;
  for ( size_t b = 0u; b < reco_bins.size(); ++b ) {
    const auto& rbin = reco_bins.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) {
      ++num_sideband_reco_bins;
      if ( !found_first_sideband_bin ) {
        found_first_sideband_bin = true;
        first_sideband_bin_idx = b;
      }
    }
  }
  // There is one more edge than the number of sideband bins
  sb_file << num_sideband_reco_bins + 1;
  for ( size_t e = 0u; e <= num_sideband_reco_bins; ++e ) {
    sb_file << ' ' << e + first_sideband_bin_idx;
  }
  sb_file << '\n';
  // For the "sideband slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker sideband bin and assign it to the
  // matching ROOT histogram bin
  sb_file << num_sideband_reco_bins << '\n';
  for ( size_t b = 0u; b < num_sideband_reco_bins; ++b ) {
    sb_file << b + first_sideband_bin_idx << " 1 " << b + 1 << '\n';
  }

}
//////////////////////////////////////////////////////////////
///
//////////////////////////////////////////////////////////////
std::vector<double> generateDoubleVector(int start, int end) {
    std::vector<double> result;
    for (int i = start; i <= end; ++i) {
        result.push_back(static_cast<double>(i));
    }
    return result;
}


//////////////////////////////////////////////////////////////
///
//////////////////////////////////////////////////////////////
std::vector<double> generateDoubleVectorWithIncrement(double start, double end, double increment) {
    std::vector<double> result;

    for (double value = start; value <= end; value += increment) {
        result.push_back(value);
    }

    return result;
}

//////////////////////////////////////////////////////////////
///
//////////////////////////////////////////////////////////////

void make_config_mcc8() {
   constexpr int DUMMY_BLOCK_INDEX = -1; 

   // Keys are the reco STV ntuple branch names of interest. Values
   // are vectors of bin edges (taken from the MCC8 CCNp0pi paper)
  
  
  std::vector <double> Costheta_binning = GetCCZeroPi_Binning(kBINNING_Costheta);
  std::vector <double> Pm_binning = GetCCZeroPi_Binning(kBINNING_Pmu);
  
  

  
  std::map< std::string, std::vector<double> > mcc8_bin_edge_map = {

    { "p3_mu.CosTheta()", {Costheta_binning} },
    { "p3_mu.Mag()", {Pm_binning} }
        
        
        
  // -1.0, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43,
      //0.59, 0.73, 0.83, 0.91, 1.0 
   // { "theta_mu_p", { 0.0, 0.8, 1.2, 1.57, 1.94, 2.34, M_PI } },

   // { "p3_lead_p.Mag()", { 0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87,
   //   0.93, 1.2 } },


    // 0.1, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5
   // { "p3_lead_p.CosTheta()", { -1.0, -0.5, 0.0, 0.27, 0.45, 0.62, 0.76, 0.86,
   //   0.94, 1.0 } },

  };
  
  



  // Used for converting from the variable names given in the bin edge
  // map to the ones used by the SliceBinning object
  std::map< std::string, std::string > var_name_map = {
    { "p3_mu.CosTheta()", "reco cos#theta_{#mu}" },
    { "p3_mu.Mag()", "reco p_{#mu}" }
  };

  // Keys are the same reco variable branch expressions as above. Values
  // are bool pairs indicating whether an underflow and overflow bin
  // should be produced.
  std::map< std::string, std::pair<bool,bool> > mcc8_under_overflow_map = {

    // Restricted by valid angular ranges to lie within the bin limits
    { "p3_mu.CosTheta()", { false, false } },

    //{ "p3_lead_p.CosTheta()", { false, false } },

    //{ "theta_mu_p", { false, false } },

    // Restricted by the signal definition to lie within the proton momentum
    // bins
    //{ "p3_lead_p.Mag()", { false, false } },

    // No underflow bin is needed due to the signal definition. No upper limit
    // is imposed, however, so we will make an overflow bin.
    { "p3_mu.Mag()", { false, true } }

  };

  // Set up an initially empty container to hold the slice definitions. We'll
  // populate it in parallel with defining the bins themselves.
  SliceBinning sb;

  // Set the variables to use when defining phase-space slices
  sb.slice_vars_ = {
    { "reco p_{#mu}", "GeV/c", "reco $p_{\\mu}$", "GeV$/c$" },
    //{ "reco cos#theta_{#mu}", "", "reco $\\cos\\theta_{\\mu}$", "" },
    //{ "reco p_{p}", "GeV/c", "reco $p_{\\mu}$", "GeV$/c$" },
    //{ "reco cos#theta_{p}", "", "reco $\\cos\\theta_{\\mu}$", "" },
    { "reco #theta_{#mup}", "", "reco $\\theta_{\\mu p}$", "" },
    { "reco bin number", "", "reco bin number", "" }
  };

  // NOTE: this script assumes that the definitions for the selection flag
  // (sel_CCNp0pi) and signal flag (mc_is_signal) have been suitably changed to
  // match the MCC8 analysis. The user is responsible for using ntuples that
  // have been post-processed consistently. Strictly speaking, this includes
  // the tiny difference between the pionless and mesonless signal definitions
  // (a ~0.06% effect).
  std::string selection = "sel_CC0pi";
  std::string signal_def = "mc_is_cc0pi_signal";

  // By construction, MC event categories 5-11 contain all beam-correlated
  // backgrounds. This list is therefore comprehensive apart from cosmic
  // overlay stuff which is directly measured in a dedicated sample.
  // NOTE: We add an extra background bin for events that the MCC9 analysis
  // considers signal and the MCC8 analysis considers background.
  std::vector< std::string > background_defs = {
   "category == 5", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };
  
  
   std::map< std::string, std::string > sideband_selection_to_momentum_map = {
    {  DIRT_SIDEBAND_SELECTION,           "p3_mu"    },
    {  NC_SIDEBAND_SELECTION,             "p3_lead_p"},
    { CCNPI_SIDEBAND_SELECTION,           "p3_mu"    }
  /*    { BDT_Bogus_SIDEBAND_SELECTION ,      "p3_mu"},
    { BDT_Else_SIDEBAND_SELECTION ,       "p3_mu"},
    { BDT_Bogus_Else_SIDEBAND_SELECTION , "p3_mu"}*/
      };
    

  std::map<std::string, std::vector<double> > sidebandBinning{

   { DIRT_SIDEBAND_SELECTION,           Pm_binning },
   { NC_SIDEBAND_SELECTION,             {0.250, 0.325, 0.4,  0.5,  0.6,  0.7, 1.0}},
   { CCNPI_SIDEBAND_SELECTION,           Pm_binning}
  /*  { BDT_Bogus_SIDEBAND_SELECTION ,      Pm_binning},
   { BDT_Else_SIDEBAND_SELECTION ,       Pm_binning},
   { BDT_Bogus_Else_SIDEBAND_SELECTION , Pm_binning}*/
  };

  std::map<std::string, size_t > sidebandNBins{

   { DIRT_SIDEBAND_SELECTION,            7 },
   { NC_SIDEBAND_SELECTION,              6},
   { CCNPI_SIDEBAND_SELECTION,           7},
  /*  { BDT_Bogus_SIDEBAND_SELECTION ,      7},
   { BDT_Else_SIDEBAND_SELECTION ,       7},
   { BDT_Bogus_Else_SIDEBAND_SELECTION , 7}*/
 };
    
    
     std::map<std::string, size_t > sideband_OverFlowBinN{

   { DIRT_SIDEBAND_SELECTION,            27},
   { NC_SIDEBAND_SELECTION,              34},
   { CCNPI_SIDEBAND_SELECTION,           48},
  /*  { BDT_Bogus_SIDEBAND_SELECTION ,      55},
   { BDT_Else_SIDEBAND_SELECTION ,       61},
   { BDT_Bogus_Else_SIDEBAND_SELECTION , 71}*/
  };
    
   std::vector<double >OverFlowBinN_vector{27,34,48,55,61,71};

  

  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;

  // Create separate blocks of bins for each kinematic variable using the
  // bin definitions from the MCC8 CCNp0pi analysis
  int block_idx = -1;
  for ( const auto& pair : mcc8_bin_edge_map ) {
    // Start a new block of related bins
    ++block_idx;

    std::string reco_branchexpr = pair.first;
    std::string true_branchexpr = "mc_" + reco_branchexpr;

    // Get the index for the "active" variable in the current block. We will
    // use it below to make a new slice while also defining the bins in the
    // block.
    const std::string& act_var_name = var_name_map.at( reco_branchexpr );
    int act_var_idx = find_slice_var_index( act_var_name, sb.slice_vars_ );

    const auto flag_pair = mcc8_under_overflow_map.at( reco_branchexpr );
    bool needs_underflow_bin = flag_pair.first;
    bool needs_overflow_bin = flag_pair.second;

    // Require at least two bin edges to be present in the input vector.
    // Any variables for which this is not true will be skipped entirely.
    const auto& bin_edges = pair.second;

    size_t num_edges = bin_edges.size();
    size_t num_bins = 0u;
    if ( num_edges >= 2u ) num_bins = num_edges - 1u;
    else continue;

    // Before defining each bin, make a new Slice object and set up the
    // corresponding ROOT histogram within it
    auto& cur_slice = add_slice( sb, bin_edges, act_var_idx );

    // If needed, then create the underflow bin in both true and reco space
    if ( needs_underflow_bin ) {
      double var_underflow_max = bin_edges.front();

      std::stringstream true_ss;
      true_ss << signal_def << " && " << true_branchexpr
        << " < " << var_underflow_max;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_branchexpr
        << " < " << var_underflow_max;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so the underflow bin is always at index zero.
      cur_slice.bin_map_[ 0 ].insert( ana_bin_idx );

      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );
    }

    // Create the ordinary signal bins using the requested edges
    for ( size_t b = 0u; b < num_bins; ++b ) {

      double var_low = bin_edges.at( b );
      double var_high = bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def
        << " && " << true_branchexpr << " >= " << var_low
        << " && " << true_branchexpr << " < "  << var_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && " << reco_branchexpr << " >= " << var_low
        << " && " << reco_branchexpr << " < "  << var_high;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so we correct for that in the line below.
      cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

      // Add the completed reco bin definition to the vector
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );

    } // loop over ordinary bins for the current variable

    // If needed, then create the overflow bin in both true and reco space
    if ( needs_overflow_bin ) {
      double var_overflow_min = bin_edges.back();

      std::stringstream true_ss;
      true_ss << signal_def << " && " << true_branchexpr
        << " >= " << var_overflow_min;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_branchexpr
        << " >= " << var_overflow_min;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so the overflow bin has an index equal to
      // the number of bin edges.
      cur_slice.bin_map_[ bin_edges.size() ].insert( ana_bin_idx );

      // Add the completed reco bin definition to the vector
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );
    }

  } // loop over kinematic variables

  // Add a single set of true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin, DUMMY_BLOCK_INDEX );
  }

  // Create a slice showing all blocks together as a function of bin number
  int num_reco_bins = reco_bins.size();
  int bin_number_var_idx = find_slice_var_index( "reco bin number",
    sb.slice_vars_ );

  auto& bin_num_slice = add_slice( sb, num_reco_bins, 0, num_reco_bins,
    bin_number_var_idx );
  for ( int ab = 0; ab < num_reco_bins; ++ab ) {
    // The ROOT histogram bins are one-based, so we correct for this here
    bin_num_slice.bin_map_[ ab + 1 ].insert( ab );
  }



  for ( const auto& sel_mom_pair : sideband_selection_to_momentum_map ) {
    const auto& side_sel = sel_mom_pair.first;
    const auto& mom_branch = sel_mom_pair.second;
    std::map< double, std::vector<double> >* bin_edge_map = nullptr;
    if ( mom_branch == "p3_mu" ) bin_edge_map = &MUON_2D_BIN_EDGES_for1DSidBand;
    else if ( mom_branch == "p3_lead_p" ) bin_edge_map = &PROTON_2D_BIN_EDGES;
    else throw std::runtime_error( "Unimplemented sideband momentum!" );

    // Get an iterator to the last map element. They are sorted numerically,
    // so this will be the upper edge of the last non-overflow momentum bin.
    auto last = bin_edge_map->cend();
    --last;

    for ( auto iter = bin_edge_map->cbegin(); iter != last; ++iter ) {

      // Get an iterator to the map element after the current one. Due to the
      // automatic sorting, this is guaranteed to contain the upper edge of the
      // current momentum bin.
      auto next = iter;
      ++next;

      // Get the current momentum bin limits
      double p_low = iter->first;
      double p_high = next->first;

      std::stringstream reco_ss;
      reco_ss << side_sel
        << " && " << mom_branch << ".Mag() >= " << p_low
        << " && " << mom_branch << ".Mag() < " << p_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );
    } // 1D reco momentum bins

    // For sidebands that use bins of muon candidate momentum, create the
    // overflow bin. This isn't needed for the proton momentum due to the
    // upper limit imposed in the signal definition.
    if ( mom_branch != "p3_mu" ) continue;

    double pmu_overflow_min = bin_edge_map->crbegin()->first;

    std::stringstream reco_ss;
    reco_ss << side_sel << " && " << mom_branch << ".Mag() >= "
      << pmu_overflow_min;

    std::string reco_bin_def = reco_ss.str();
    reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );

  } // sideband selection definitions



  // Dump this information to the output file
  std::ofstream out_file( "myconfig_mcc8_CC0pi_1D_NoBDTproton_new.txt" );
  out_file << "mcc8_all" << '\n';
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_mcc8_1D_NoBDTproton.txt" );
  sb_file << sb;
  
  // Make a slice containing the sideband results organized by reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Count the number of sideband reco bins. Also find the index of the
  // first one.
  size_t num_sideband_reco_bins = 0u;
  size_t first_sideband_bin_idx = 0u;
  bool found_first_sideband_bin = false;
  for ( size_t b = 0u; b < reco_bins.size(); ++b ) {
    const auto& rbin = reco_bins.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) {
      ++num_sideband_reco_bins;
      if ( !found_first_sideband_bin ) {
        found_first_sideband_bin = true;
        first_sideband_bin_idx = b;
      }
    }
  }
  // There is one more edge than the number of sideband bins
  sb_file << num_sideband_reco_bins + 1;
  for ( size_t e = 0u; e <= num_sideband_reco_bins; ++e ) {
    sb_file << ' ' << e + first_sideband_bin_idx;
  }
  sb_file << '\n';
  // For the "sideband slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker sideband bin and assign it to the
  // matching ROOT histogram bin
  sb_file << num_sideband_reco_bins << '\n';
  for ( size_t b = 0u; b < num_sideband_reco_bins; ++b ) {
    sb_file << b + first_sideband_bin_idx << " 1 " << b + 1 << '\n';
  }

  
  
  
  /*/// added sideband 
  sb_file << '\n';
  sb_file << "\"events \"\n"; // y-axis label
  sb_file << "1 2 ";
  // Count the number of sideband reco bins. Also find the index of the
  // first one.
  size_t num_sideband_reco_bins = 0u;
  
  size_t first_sideband_bin_idx = 0u;
  bool found_first_sideband_bin = false;
  for ( size_t b = 0u; b < reco_bins.size(); ++b ) {
    const auto& rbin = reco_bins.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) {
      ++num_sideband_reco_bins;
      if ( !found_first_sideband_bin ) {
        found_first_sideband_bin = true;
        first_sideband_bin_idx = b;
      }
    }
  }
  
  // There is one more edge than the number of sideband bins
  sb_file << num_sideband_reco_bins + 1;
  for ( size_t e = 0u; e <= num_sideband_reco_bins; ++e ) {
    sb_file << ' ' << e + first_sideband_bin_idx;
  }
  
  sb_file << '\n';
  // For the "sideband slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker sideband bin and assign it to the
  // matching ROOT histogram bin
  sb_file << num_sideband_reco_bins << '\n';
  for ( size_t b = 0u; b < num_sideband_reco_bins; ++b ) {
    sb_file << b + first_sideband_bin_idx << " 1 " << b + 1 << '\n';
  }
 //// I want to add plots for each sideband sepertely for mybins
 
 
 
 int starting_sideBinplots = first_sideband_bin_idx;
 int stop_sideBinplots = first_sideband_bin_idx + 7; // number of bins for sideband 
 
 
 for(auto sideband_iter : sidebandBinning){
  
  sb_file << "\"events \"\n"; // y-axis label
  sb_file << "1 0 ";
  
  auto myVector = sideband_iter.second;
  
  sb_file << myVector.size()-1 << ' ';
  std::copy(myVector.begin(), myVector.end(), std::ostream_iterator<double>(sb_file, " "));
   sb_file << '\n'; // add line
  sb_file << "0 \n";
  size_t Nbinlength = stop_sideBinplots -  starting_sideBinplots;
  sb_file << Nbinlength << " \n";
  
  int SliceBinIndex = 1; 
   
   for ( size_t b = starting_sideBinplots; b <= stop_sideBinplots; ++b ) {
   if(b == sideband_OverFlowBinN[sideband_iter.first]) continue; 
    sb_file << b << " 1 " << SliceBinIndex << '\n';
    SliceBinIndex++;
  }
 
 
 starting_sideBinplots = starting_sideBinplots + Nbinlength; 
 
 stop_sideBinplots = starting_sideBinplots + sidebandNBins[sideband_iter.first];
 
 } 
  */

  //sb_file << "\"events \"\n"; // y-axis label
  //sb_file << "1 2 6 0 1 2 3 4 5 6";
  //sb_file << '\n'; // add line
  //sb_file << "0 \n";
  //sb_file << OverFlowBinN_vector.size() << ' ';
  //int SliceBin_overlow_Index = 1;
  //for(auto overflowsidebinN:OverFlowBinN_vector ){
  // sb_file << overflowsidebinN << " 1 " << SliceBin_overlow_Index << '\n';
  //SliceBin_overlow_Index++;
  //
  //}
  
  
  
}//End of function
/////////////////////////////////////////////////////////////

void make_config_RECOCuting() {
   constexpr int DUMMY_BLOCK_INDEX = -1; 

   // Keys are the reco STV ntuple branch names of interest. Values
   // are vectors of bin edges (taken from the MCC8 CCNp0pi paper)  
  
  //std::vector<double> Zero_to_1_score_bins = generateBins();
  
  std::vector< double > Mult= {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
  
  //vector<double> Track_distance_vectex_bins{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  //vector<double> Track_distance_bins{0,1,2,4,5,6,7,8,9,10,11,12,13,14,15};
  //std::vector<double> TrackLength = generateDoubleVector(0, 100);
  //std::vector<double> LLR_score =generateDoubleVectorWithIncrement(-1.0, 1.0, .1);


  std::map< std::string, std::vector<double> > mcc8_bin_edge_map = {

         { "proton_multiplicity", {Mult} }//,
        // { "trk_score_v", {Zero_to_1_score_bins} },
        // { "trk_distance_v", {Track_distance_vectex_bins} },
        // { "trk_len_v", {TrackLength} },
        // { "trk_llr_pid_score_v", {LLR_score} },
        // { "trk_llr_pid_score_v", {LLR_score} },
        
  };
  
  



  // Used for converting from the variable names given in the bin edge
  // map to the ones used by the SliceBinning object
  std::map< std::string, std::string > var_name_map = {
        { "proton_multiplicity","Track Multiplicity" }
         //{ "trk_score_v", "Track Score" },
         //{ "trk_distance_v", "Distance to Vertex [cm]" },
         //{ "trk_len_v", "Tracklength" },
         //{ "trk_llr_pid_score_v", "LLR PID Track score"}
  };

  // Keys are the same reco variable branch expressions as above. Values
  // are bool pairs indicating whether an underflow and overflow bin
  // should be produced.
  std::map< std::string, std::pair<bool,bool> > mcc8_under_overflow_map = {

    // Restricted by valid angular ranges to lie within the bin limits
    // Not going worry about over/under flows
        { "proton_multiplicity",{ false, false } }
         //{ "trk_score_v", { false, false } },
         //{ "trk_distance_v", { false, false } },
         //{ "trk_len_v", { false, false } },
         //{ "trk_llr_pid_score_v", { false, false }}

  };

  // Set up an initially empty container to hold the slice definitions. We'll
  // populate it in parallel with defining the bins themselves.
  SliceBinning sb;

  // Set the variables to use when defining phase-space slices
  sb.slice_vars_ = {
    { "proton multiplicity", "", "proton multiplicity", "" },
    //{ "Track Score", "", "Track Score", "" },
    //{ "Distance to Vertex", "cm", "Distance to Vertex", "cm" },
    //{ "Tracklength", "cm", "Tracklength", "cn" },
    //{ "LLR PID Track score", "", "LLR PID Track score", "" }
 
  };

  // NOTE: this script assumes that the definitions for the selection flag
  // (sel_CCNp0pi) and signal flag (mc_is_signal) have been suitably changed to
  // match the MCC8 analysis. The user is responsible for using ntuples that
  // have been post-processed consistently. Strictly speaking, this includes
  // the tiny difference between the pionless and mesonless signal definitions
  // (a ~0.06% effect).
  std::string selection = "sel_CC0pi && sel_BDT_NotBOGUS";
  std::string signal_def = "mc_is_cc0pi_signal";

  // By construction, MC event categories 5-11 contain all beam-correlated
  // backgrounds. This list is therefore comprehensive apart from cosmic
  // overlay stuff which is directly measured in a dedicated sample.
  // NOTE: We add an extra background bin for events that the MCC9 analysis
  // considers signal and the MCC8 analysis considers background.
  std::vector< std::string > background_defs = {
   "category == 5", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };
  
  
  

  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;

  // Create separate blocks of bins for each kinematic variable using the
  // bin definitions from the MCC8 CCNp0pi analysis
  int block_idx = -1;
  for ( const auto& pair : mcc8_bin_edge_map ) {
    // Start a new block of related bins
    ++block_idx;

    std::string reco_branchexpr = pair.first;
    std::string true_branchexpr = pair.first; //"mc_" + reco_branchexpr; // I think there is not truth branch for topolgy

    // Get the index for the "active" variable in the current block. We will
    // use it below to make a new slice while also defining the bins in the
    // block.
    const std::string& act_var_name = var_name_map.at( reco_branchexpr );
    int act_var_idx = find_slice_var_index( act_var_name, sb.slice_vars_ );

    const auto flag_pair = mcc8_under_overflow_map.at( reco_branchexpr );
    bool needs_underflow_bin = flag_pair.first;
    bool needs_overflow_bin = flag_pair.second;

    // Require at least two bin edges to be present in the input vector.
    // Any variables for which this is not true will be skipped entirely.
    const auto& bin_edges = pair.second;

    size_t num_edges = bin_edges.size();
    size_t num_bins = 0u;
    if ( num_edges >= 2u ) num_bins = num_edges - 1u;
    else continue;

    // Before defining each bin, make a new Slice object and set up the
    // corresponding ROOT histogram within it
    auto& cur_slice = add_slice( sb, bin_edges, act_var_idx );

    // If needed, then create the underflow bin in both true and reco space
    if ( needs_underflow_bin ) {
      double var_underflow_max = bin_edges.front();

      std::stringstream true_ss;
      true_ss << signal_def << " && " << true_branchexpr
        << " < " << var_underflow_max;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_branchexpr
        << " < " << var_underflow_max;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so the underflow bin is always at index zero.
      cur_slice.bin_map_[ 0 ].insert( ana_bin_idx );

      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );
    }

    // Create the ordinary signal bins using the requested edges
    for ( size_t b = 0u; b < num_bins; ++b ) {

      double var_low = bin_edges.at( b );
      double var_high = bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_def
        << " && " << true_branchexpr << " >= " << var_low
        << " && " << true_branchexpr << " < "  << var_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection
        << " && " << reco_branchexpr << " >= " << var_low
        << " && " << reco_branchexpr << " < "  << var_high;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so we correct for that in the line below.
      cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

      // Add the completed reco bin definition to the vector
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );

    } // loop over ordinary bins for the current variable

    // If needed, then create the overflow bin in both true and reco space
    if ( needs_overflow_bin ) {
      double var_overflow_min = bin_edges.back();

      std::stringstream true_ss;
      true_ss << signal_def << " && " << true_branchexpr
        << " >= " << var_overflow_min;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection << " && " << reco_branchexpr
        << " >= " << var_overflow_min;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so the overflow bin has an index equal to
      // the number of bin edges.
      cur_slice.bin_map_[ bin_edges.size() ].insert( ana_bin_idx );

      // Add the completed reco bin definition to the vector
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );
    }

  } // loop over kinematic variables

  // Add a single set of true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin, DUMMY_BLOCK_INDEX );
  }

  // Create a slice showing all blocks together as a function of bin number
  int num_reco_bins = reco_bins.size();
  int bin_number_var_idx = find_slice_var_index( "reco bin number",
    sb.slice_vars_ );

  auto& bin_num_slice = add_slice( sb, num_reco_bins, 0, num_reco_bins,
    bin_number_var_idx );
  for ( int ab = 0; ab < num_reco_bins; ++ab ) {
    // The ROOT histogram bins are one-based, so we correct for this here
    bin_num_slice.bin_map_[ ab + 1 ].insert( ab );
  }

// Dont need a sideband to make this plot 
/*
  for ( const auto& sel_mom_pair : sideband_selection_to_momentum_map )
  {
    const auto& side_sel = sel_mom_pair.first;
    const auto& mom_branch = sel_mom_pair.second;
    std::map< double, std::vector<double> >* bin_edge_map = nullptr;
    if ( mom_branch == "p3_mu" ) bin_edge_map = &MUON_2D_BIN_EDGES_for1DSidBand;
    else if ( mom_branch == "p3_lead_p" ) bin_edge_map = &PROTON_2D_BIN_EDGES;
    else throw std::runtime_error( "Unimplemented sideband momentum!" );

    // Get an iterator to the last map element. They are sorted numerically,
    // so this will be the upper edge of the last non-overflow momentum bin.
    auto last = bin_edge_map->cend();
    --last;

    for ( auto iter = bin_edge_map->cbegin(); iter != last; ++iter ) {

      // Get an iterator to the map element after the current one. Due to the
      // automatic sorting, this is guaranteed to contain the upper edge of the
      // current momentum bin.
      auto next = iter;
      ++next;

      // Get the current momentum bin limits
      double p_low = iter->first;
      double p_high = next->first;

      std::stringstream reco_ss;
      reco_ss << side_sel
        << " && " << mom_branch << ".Mag() >= " << p_low
        << " && " << mom_branch << ".Mag() < " << p_high;

      std::string reco_bin_def = reco_ss.str();

      reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );
    } // 1D reco momentum bins

    // For sidebands that use bins of muon candidate momentum, create the
    // overflow bin. This isn't needed for the proton momentum due to the
    // upper limit imposed in the signal definition.
    if ( mom_branch != "p3_mu" ) continue;

    double pmu_overflow_min = bin_edge_map->crbegin()->first;

    std::stringstream reco_ss;
    reco_ss << side_sel << " && " << mom_branch << ".Mag() >= "
      << pmu_overflow_min;

    std::string reco_bin_def = reco_ss.str();
    reco_bins.emplace_back( reco_bin_def, kSidebandRecoBin );

  } // sideband selection definitions
*/


  // Dump this information to the output file
  std::ofstream out_file( "myconfig_mcc8_CC0pi_NTracks.txt" );
  out_file << "mcc8_all" << '\n';
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_mcc8_1D_NTracks.txt" );
  sb_file << sb;
  
  // Make a slice containing the sideband results organized by reco bin number
  sb_file << "\"events\"\n"; // y-axis label
  sb_file << "1 2 ";
  // Count the number of sideband reco bins. Also find the index of the
  // first one.
  size_t num_sideband_reco_bins = 0u;
  size_t first_sideband_bin_idx = 0u;
  bool found_first_sideband_bin = false;
  for ( size_t b = 0u; b < reco_bins.size(); ++b ) {
    const auto& rbin = reco_bins.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) {
      ++num_sideband_reco_bins;
      if ( !found_first_sideband_bin ) {
        found_first_sideband_bin = true;
        first_sideband_bin_idx = b;
      }
    }
  }
  
// I dont think I need this , 
/*
  // There is one more edge than the number of sideband bins
  sb_file << num_sideband_reco_bins + 1;
  for ( size_t e = 0u; e <= num_sideband_reco_bins; ++e ) {
    sb_file << ' ' << e + first_sideband_bin_idx;
  }
  sb_file << '\n';
  // For the "sideband slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker sideband bin and assign it to the
  // matching ROOT histogram bin
  sb_file << num_sideband_reco_bins << '\n';
  for ( size_t b = 0u; b < num_sideband_reco_bins; ++b ) {
    sb_file << b + first_sideband_bin_idx << " 1 " << b + 1 << '\n';
  }
*/
  
  
  
  /*/// added sideband 
  sb_file << '\n';
  sb_file << "\"events \"\n"; // y-axis label
  sb_file << "1 2 ";
  // Count the number of sideband reco bins. Also find the index of the
  // first one.
  size_t num_sideband_reco_bins = 0u;
  
  size_t first_sideband_bin_idx = 0u;
  bool found_first_sideband_bin = false;
  for ( size_t b = 0u; b < reco_bins.size(); ++b ) {
    const auto& rbin = reco_bins.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) {
      ++num_sideband_reco_bins;
      if ( !found_first_sideband_bin ) {
        found_first_sideband_bin = true;
        first_sideband_bin_idx = b;
      }
    }
  }
  
  // There is one more edge than the number of sideband bins
  sb_file << num_sideband_reco_bins + 1;
  for ( size_t e = 0u; e <= num_sideband_reco_bins; ++e ) {
    sb_file << ' ' << e + first_sideband_bin_idx;
  }
  
  sb_file << '\n';
  // For the "sideband slice," there is no other variable apart from reco bin
  // number
  sb_file << "0\n";
  // Loop over each UniverseMaker sideband bin and assign it to the
  // matching ROOT histogram bin
  sb_file << num_sideband_reco_bins << '\n';
  for ( size_t b = 0u; b < num_sideband_reco_bins; ++b ) {
    sb_file << b + first_sideband_bin_idx << " 1 " << b + 1 << '\n';
  }
 //// I want to add plots for each sideband sepertely for mybins
 
 
 
 int starting_sideBinplots = first_sideband_bin_idx;
 int stop_sideBinplots = first_sideband_bin_idx + 7; // number of bins for sideband 
 
 
 for(auto sideband_iter : sidebandBinning){
  
  sb_file << "\"events \"\n"; // y-axis label
  sb_file << "1 0 ";
  
  auto myVector = sideband_iter.second;
  
  sb_file << myVector.size()-1 << ' ';
  std::copy(myVector.begin(), myVector.end(), std::ostream_iterator<double>(sb_file, " "));
   sb_file << '\n'; // add line
  sb_file << "0 \n";
  size_t Nbinlength = stop_sideBinplots -  starting_sideBinplots;
  sb_file << Nbinlength << " \n";
  
  int SliceBinIndex = 1; 
   
   for ( size_t b = starting_sideBinplots; b <= stop_sideBinplots; ++b ) {
   if(b == sideband_OverFlowBinN[sideband_iter.first]) continue; 
    sb_file << b << " 1 " << SliceBinIndex << '\n';
    SliceBinIndex++;
  }
 
 
 starting_sideBinplots = starting_sideBinplots + Nbinlength; 
 
 stop_sideBinplots = starting_sideBinplots + sidebandNBins[sideband_iter.first];
 
 } 
  */

  //sb_file << "\"events \"\n"; // y-axis label
  //sb_file << "1 2 6 0 1 2 3 4 5 6";
  //sb_file << '\n'; // add line
  //sb_file << "0 \n";
  //sb_file << OverFlowBinN_vector.size() << ' ';
  //int SliceBin_overlow_Index = 1;
  //for(auto overflowsidebinN:OverFlowBinN_vector ){
  // sb_file << overflowsidebinN << " 1 " << SliceBin_overlow_Index << '\n';
  //SliceBin_overlow_Index++;
  //
  //}
  
  
  
}//End of function
/////////////////////////////////////////////////////////////
int main() {

 
std::cout<<" Testing Function  "<< std::endl;
make_config_mcc9_2D_muon_CC0Pi();
make_config_mcc9_2D_muon_CC0Pi_inclusive();
  //tutorial_slice_plots();
  //make_config_mcc8();
  //make_config_RECOCuting();
  std::cout<<"Finsihed "<< std::endl;
  return 0;
}

std::map<int, double> readTextFileToMap(const char* fileName) {
    std::map<int, double> dataMap;

    // Open the file for reading
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open the file " << fileName << " for reading." << std::endl;
        return dataMap;  // Return an empty map on error
    }
    // Read each line of the file

    std::string line;
    while(std::getline(inputFile, line)) {
        std::istringstream iss{line};
        std::string entry;
        std::vector<std::string> Row_entries;
        // Tokenize the line using a comma as a delimiter
        while (std::getline(iss, entry, ',')) {
          Row_entries.push_back(entry);
        }
          
          int key = std::stoi(Row_entries[0]);
        
          double value = std::stod(Row_entries[1]);
          double value_rounded = std::round(value * 1000.0) / 1000.0;
            // Try to read the key and value from the entry
            dataMap[key] = value_rounded;
        
        
    }

    // Close the file
    inputFile.close();

    return dataMap;
}
std::map<int, double> readTextFileToMap_2(const char* fileName) {
    std::map<int, double> dataMap;

    // Open the file for reading
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open the file " << fileName << " for reading." << std::endl;
        return dataMap;  // Return an empty map on error
    }
    // Read each line of the file

    std::string line;
    while(std::getline(inputFile, line)) {
        std::istringstream iss{line};
        std::string entry;
        std::vector<std::string> Row_entries;
        // Tokenize the line using a comma as a delimiter
        while (std::getline(iss, entry, ',')) {
          Row_entries.push_back(entry);
        }
          
          int key = std::stoi(Row_entries[0]);
        
          double value = std::stod(Row_entries[1]);
          double value_rounded = std::round(value * 10.0) / 10.0;
            // Try to read the key and value from the entry
            dataMap[key] = value_rounded;
        
        
    }

    // Close the file
    inputFile.close();

    return dataMap;
}