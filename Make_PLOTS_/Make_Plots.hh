#pragma once
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
#include "TLine.h"
#include "../includes/GridCanvas.hh"
#include "TH2Poly.h"
#include "../includes/UBTH2Poly.h" 
#include "../includes/HistUtils_cc0pi.hh"
//#include "HistUtils.hh"
#include "../ConfigMakerUtils.hh"

//using NFT = NtupleFileType;
std::string MicroBooNE_LegendTitle();
