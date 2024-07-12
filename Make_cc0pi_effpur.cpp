// Efficiencies and purities for the MCC9 CCNp0pi STV analysis
//
// ATM 2022/07/15: Modified for CC0pi

// For a given analysis TTree (containing the events to use), compute the
// efficiency and purity given definitions for what events count as signal
// (based on MC truth information) and what events are selected (based on
// reconstructed information) in terms of our analysis TTree branch variables.
// Store the results in the output variables eff and pur.
/// Rewritten as a excluable 
// Christian Nguyen 
// Jan/2024

#include "Make_cc0pi_eff.hh"

#include <TFormula.h>


void compute_eff_pur( TTree& stv_tree, const std::string& signal_cuts,
  const std::string& selection_cuts, double& eff, double& pur )
{
  // For computing efficiencies and purities, we need to only use MC events.
  // Unconditionally add this requirement to the cuts defined above.
  std::string signal = signal_cuts + " && is_mc ";
  std::string selection = selection_cuts + " && is_mc ";

  // These are actually integer counts, but since we will use them below to
  // compute ratios, intrinsically cast them to double-precision values for
  // convenience.
  double num_signal = stv_tree.Draw( "", signal.c_str(), "goff" );
  double num_selected = stv_tree.Draw( "", selection.c_str(), "goff" );
  double num_selected_signal = stv_tree.Draw( "",
    (signal + " && " + selection).c_str(), "goff" );

  eff = num_selected_signal / num_signal;
  pur = num_selected_signal / num_selected;

  //std::cout << "signal = " << num_signal << '\n';
  //std::cout << "selected = " << num_selected << '\n';
  //std::cout << "selected_signal = " << num_selected_signal << '\n';

}

void cc0pi_effpur() {


 std::cout<<"Starting MakePlots"<< std::endl;
   TCanvas* c1 = new TCanvas;

   
    char text_title_pdf1[2024];
    char text_title_pdf2[2024];
    char text_title_pdf3[2024];
    char text_title_pdf4[2024];

  sprintf(text_title_pdf1, "CC0pi_CUT_eff_pur_Figures_v2.pdf(");
  c1 -> Print(text_title_pdf1);
  sprintf(text_title_pdf2, "CC0pi_CUT_eff_pur_Figures_v2.pdf" );
  sprintf(text_title_pdf3, "CC0pi_CUT_eff_pur_Figures_v2.pdf)");
  sprintf(text_title_pdf4, "CC0pi_CUT_eff_pur_Figures_v2");
  std::string text_title_pdf4_string(text_title_pdf4);
  std::string text_title_pdf2_string(text_title_pdf2);




  const std::vector< std::string > signal_defs = {
    "1",
    "mc_vertex_in_FV",
    "mc_vertex_in_FV && mc_neutrino_is_numu",
    "mc_vertex_in_FV && mc_neutrino_is_numu && mc_muon_in_mom_range",
    "mc_vertex_in_FV && mc_neutrino_is_numu && mc_muon_in_mom_range && mc_no_fs_pi0",
    "mc_vertex_in_FV && mc_neutrino_is_numu && mc_muon_in_mom_range && mc_no_fs_pi0 && mc_no_charged_pi_above_threshold",
    "mc_is_cc0pi_signal"
  };




  const std::vector< std::string > selection_defs = {
    "1",
    "sel_reco_vertex_in_FV",
    "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV",
    "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV && sel_topo_cut_passed",
    "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV && sel_topo_cut_passed && sel_no_reco_showers",
    "sel_presel && sel_has_muon_candidate",
    "sel_nu_mu_cc && sel_muon_passed_mom_cuts", 
    "sel_nu_mu_cc && sel_muon_passed_mom_cuts && !sel_has_pion_candidate;"
    /*"sel_CC0pi && sel_BDT_NotBOGUS && sel_BDT_predicts_1plusTrks_tobeProtons",
    "sel_CC0pi && sel_BDT_NotBOGUS && sel_BDT_predicts_1plusTrks_tobeProtons && sel_cosmic_ip_cut_passed",
    "sel_CC0pi && sel_BDT_NotBOGUS && sel_BDT_predicts_1plusTrks_tobeProtons && sel_cosmic_ip_cut_passed && sel_muon_distance_from_Vertex_ok;"*/
  };
  

  
    const std::vector< std::string > bin_labels = {
    "No Cuts",
    "In FV",
    "Starts Contained",
    "Topology Cut",
    "No Showers",
    "CCincl (BDT #mu Candidate)",
    "#mu Momentum Limits",
    "No #pi^{#pm} Above Threshold"
  };
  
  
      const std::vector< std::string > selection_defs_latex = {
       "1",
    "sel\\_reco\\_vertex\\_in\\_FV",
    "sel\\_reco\\_vertex\\_in\\_FV \\&\\& sel\\_pfp\\_starts\\_in\\_PCV",
    "sel\\_reco\\_vertex\\_in\\_FV \\&\\& sel\\_pfp\\_starts\\_in\\_PCV \\&\\& sel\\_topo\\_cut\\_passed",
    "sel\\_reco\\_vertex\\_in\\_FV \\&\\& sel\\_pfp\\_starts\\_in\\_PCV \\&\\& sel\\_topo\\_cut\\_passed \\&\\& sel\\_no\\_reco\\_showers",
    "sel\\_presel\\ \\&\\& sel\\_has\\_muon\\_candidate",
    "sel\\_nu\\_mu\\_cc \\&\\& sel\\_muon\\_passed\\_mom\\_cuts", 
    "sel\\_nu\\_mu\\_cc \\&\\& sel\\_muon\\_passed\\_mom\\_cuts \\&\\& !sel\\_has\\_pion\\_candidate"  };
  //&& !std::isnan(BTD_Muon_Candidate_prediction) && !std::isinf(BTD_Muon_Candidate_prediction) 
/*
  sel_CC0pi_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_muon_passed_mom_cuts_ && sel_muon_quality_ok_ //&& sel_muon_contained_
    && !sel_pions_above_threshold_;
*/
  TChain stv_ch( "stv_tree" );
  //stv_ch.Add( "/uboone/data/users/mastbaum/stv-analysis-ru3/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root" );
  //stv_ch.Add( "test.root" );
  stv_ch.Add( "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_4_29_2024_PmuCorrection/ru-stv_overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root");
             ///uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_12_24_24/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
  ///uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_12_15_24/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
  
  //TFormula formula("CheckVector_Entries_greaterthanConstant", "CheckVector_Entries_greaterthanConstant(x, [constant])");
  //formula.SetParameter(0, .28);
     std::ofstream outFile("Latex_Cut_table_newTuples.txt");

    // Check if the file is successfully opened
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open the file Latex_Cut_table.txt " << std::endl;
        return;
    }
  
     // Write the LaTeX table header
    outFile << "\\begin{table}" << std::endl;
    outFile << "\\centering" << std::endl;
    outFile << "\\begin{tabular}{||c||c|c|c||}" << std::endl;
    outFile << "\\hline" << std::endl;
    outFile << "Label & Eff & Purity & bool cut pars \\\\ \\hline" << std::endl;


  
  size_t num_points = selection_defs.size();
  TGraph* eff_graph = new TGraph( num_points );
  TGraph* pur_graph = new TGraph( num_points );

  std::string signal = signal_defs.back();
  double eff, pur;
  for ( size_t k = 0u; k < num_points; ++k  ) {

    const auto& selection = selection_defs.at( k );
    const auto& selection_latex =selection_defs_latex.at( k );
    const auto& Label = bin_labels.at( k );
    
    compute_eff_pur( stv_ch, signal, selection, eff, pur );

    eff_graph->SetPoint( k, k + 1, eff );
    pur_graph->SetPoint( k, k + 1, pur );

    std::cout << "Label = " << Label << '\n';
    std::cout << "eff = " << eff << '\n';
    std::cout << "pur = " << pur << '\n';
    std::cout << "\n\n";
    

        outFile << Label << " & " << eff << " & " << pur << " & " << selection_latex << "\\\\" << std::endl;
        outFile << "\\hline" << std::endl;
    
    
  }


    //"sel_nu_mu_cc && sel_no_reco_showers && sel_muon_passed_mom_cuts && sel_muon_quality_ok_ && !sel_pions_above_threshold",



  
  c1->SetBottomMargin(0.27);
  c1->SetGrid();
  eff_graph->SetTitle( ";;Efficiency or Purity" );
  eff_graph->SetLineColor(kBlue);
  eff_graph->SetMarkerColor(kBlue);
  eff_graph->SetLineWidth(3);
  eff_graph->SetMarkerStyle(20);
  eff_graph->GetYaxis()->SetRangeUser( 0., 1.2 );
  eff_graph->Draw( "alp" );

  for ( int b = 1; b <= bin_labels.size(); ++b ) {
    eff_graph->GetHistogram()->GetXaxis()
      ->SetBinLabel( eff_graph->GetHistogram()->FindBin(b),
      bin_labels.at(b - 1).c_str() );
  }
  eff_graph->Draw( "same" );

  pur_graph->SetLineColor(kRed);
  pur_graph->SetMarkerColor(kRed);
  pur_graph->SetLineWidth(3);
  pur_graph->SetMarkerStyle(20);
  pur_graph->Draw("same lp");

  TLegend* lg = new TLegend(0.65, 0.75, 0.85, 0.85);
  lg->AddEntry( eff_graph, "efficiency", "lp" );
  lg->AddEntry( pur_graph, "purity", "lp" );
  lg->SetTextSize(.03); //
  lg->Draw("same");
  
 c1 -> Print(text_title_pdf2);
 c1 -> Print(text_title_pdf3);


    // Loop to add rows to the table


    // Write the table footer
    outFile << "\\hline" << std::endl;
    outFile << "\\end{tabular}" << std::endl;
    outFile << "\\caption{Your table caption here.}" << std::endl;
    outFile << "\\label{your_label}" << std::endl;
    outFile << "\\end{table}" << std::endl;

    outFile.close(); 

}



int main() {


std::cout<<" Running  Make_cc0pi_effpur "<< std::endl;

  //tutorial_slice_plots();
  cc0pi_effpur();
  return 0;
}


// Custom function to be used in the TFormula
Double_t CheckVector_Entries_greaterthanConstant(const std::vector<double>& values, Double_t constant) {
    for (const auto& value : values) {
        if (value <= constant) {
            return 0.0;  // Return 0 if any value is not greater than the constant
        }
    }
    return 1.0;  // Return 1 if all values are greater than the constant
}