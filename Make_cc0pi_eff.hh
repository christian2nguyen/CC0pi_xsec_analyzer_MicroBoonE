//////////////
///// I want the Eff and Purity macrosn as excluables 
///////////
// Chistian Nguyen 
// 1/2024

#pragma once
// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"
#include <TGraphAsymmErrors.h>
#include "TEfficiency.h"


#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TPad.h"
#include "TGraph.h"
#include <algorithm>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TChain.h"
#include "TLegendEntry.h"
#include <TGraphErrors.h>
#include "TLine.h"
#include "TAttLine.h"
#include "TF1.h"


// STV analysis includes
#include "includes/EventCategory.hh"
#include "FiducialVolume.hh"
#include "includes/TreeUtils.hh"
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "includes/PlotUtils.hh"
#include "includes/HistUtils_cc0pi.hh"
#include "includes/SliceBinning.hh"
#include "SliceHistogram.hh"

void RunEfficnecy();
void cc0pi_effpur();
void compute_eff_pur( TTree& stv_tree, const std::string& signal_cuts,
  const std::string& selection_cuts, double& eff, double& pur );
void SetMarkerColorTEfficiency(TEfficiency *efficiency, Color_t color);

void DrawEfficiency(TH1D* h_TRUE_input,
TH1D* h_TRUE_RECO_input , 
TH1D* h_RECO_input,
TH1D* h_RECO_PURE_input, 
std::string Xaxis,
std::string Title, 
TCanvas* Can, std::string PDF_name );