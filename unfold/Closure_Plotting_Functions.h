//Header of GENIE Closure
#pragma once
// Standard library includes
// Standard library includes
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"
#include <TStyle.h> // Include TStyle header for gStyle
#include <THStack.h>
// STV analysis includes
#include "../DAgostiniUnfolder.hh"
#include "../FiducialVolume.hh"
#include "../MatrixUtils.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../NormShapeCovMatrix.hh"
#include "../PGFPlotsDumpUtils.hh"
#include "../includes/SliceBinning.hh"
//#include "../SliceHistogram.hh"
#include "../WienerSVDUnfolder.hh"
#include "../includes/PlotUtils.hh"

#include "../ConfigMakerUtils.hh"

/////////////////////////////////////////////////////////
const std::string MicroBooNEType_string =  "MicroBooNE Tune";
constexpr double BIG_DOUBLE = 1e300;
constexpr bool PrintStatement_Debug = false;
constexpr bool PrintStatement_Uncertainty_Debug = false;

/////////////////////////////////////////////////////////////
// Function
/////////////////////////////////////////////////////////////
void multiply_1d_hist_by_matrix( TMatrixD* mat, TH1* hist ) ;
void IncreaseTitleTH1(TH1& hist, double input);
void printMatrixAsLatexTable(const TMatrixD& matrix, const std::string& fileName);
void unfoldingGENIEClosure();