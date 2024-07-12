// Starting Authors: Serenity Engel, Andy Mastbaum
// rewriiten : Christian Nguyen 
// Jan/2024

#include "Make_cc0pi_eff.hh"


   


void RunEfficnecy(){

TCanvas* c2 = new TCanvas;

 std::cout<<"Starting MakePlots"<< std::endl;
   

   
    char text_title_pdf1[2024];
    char text_title_pdf2[2024];
    char text_title_pdf3[2024];
    char text_title_pdf4[2024];

  sprintf(text_title_pdf1, "CC0pi_eff_pur_Figures_noProtonBDTcut.pdf(","" );
  c2 -> Print(text_title_pdf1);
  sprintf(text_title_pdf2, "CC0pi_eff_pur_Figures_noProtonBDTcut.pdf","" );
  sprintf(text_title_pdf3, "CC0pi_eff_pur_Figures_noProtonBDTcut.pdf)","" );
  sprintf(text_title_pdf4, "CC0pi_eff_pur_Figures_noProtonBDTcut","" );
  std::string text_title_pdf4_string(text_title_pdf4);
  std::string text_title_pdf2_string(text_title_pdf2);


///uboone/data/users/mastbaum/stv-analysis-ru3/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
//TFile* f = TFile::Open("/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_12_15_24/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root");


//TFile* f_new= TFile::Open("/uboone/data/users/cnguyen/CC0Pi_Selection/MC/MC_Selected_newCC0pi_addthresholds_rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root");
TFile* f_new= TFile::Open("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/MC/MC_Selected_EventSelection_5_13_2024_CC0pi_ru-stv_overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root");
                           
//TTree* stv_tree = (TTree*) f->Get("stv_tree");



//stv_tree->Draw("p3_mu.Mag()>>h(20,0,2)", "mc_is_cc0pi_signal", "goff");
TH1D* h = GetTH1DHist(*f_new ,"h_p_TRUE" );//(TH1F*) gDirectory->Get("h");
TH1D* h_input = (TH1D*)h->Clone(uniq());
//stv_tree->Draw("p3_mu.Mag()>>h2(20,0,2)", "mc_is_cc0pi_signal&sel_CC0pi", "goff");
//TH1F* h2 = (TH1F*) gDirectory->Get("h2");
TH1D* h2 = GetTH1DHist(*f_new ,"h_p_TRUE_RECO" );

TH1D* h_Eff_clone = (TH1D*)h2->Clone(uniq());
TH1D* h2_input = (TH1D*)h2->Clone(uniq());
TH1D* h2_clone = (TH1D*)h->Clone(uniq());

TH1D* Dom_eff_times_pur_1 = (TH1D*)h->Clone(uniq());
TH1D* Nom_eff_times_pur_1 = (TH1D*)h2->Clone(uniq());

h_Eff_clone->Divide(h_Eff_clone,h2_clone);



TEfficiency* eff = new TEfficiency(*h2, *h);
eff->SetTitle(";Muon momentum p_{#mu} [GeV/c];Efficiency");
eff->Draw();
c2 -> Print(text_title_pdf2);

// Purity = (selected signal events) / (all selected events)

//TCanvas* c3 = new TCanvas;

//stv_tree->Draw("p3_mu.Mag()>>h3(20,0,2)", "sel_CC0pi", "goff");

//TH1F* h3 = (TH1F*) gDirectory->Get("h3");
TH1D* h3 = GetTH1DHist(*f_new ,"h_p" );
h3->SetMarkerColor(kRed);
//stv_tree->Draw("p3_mu.Mag()>>h4(20,0,2)", "mc_is_cc0pi_signal&sel_CC0pi", "goff");
//TH1F* h4 = (TH1F*) gDirectory->Get("h4");
TH1D* h4 = GetTH1DHist(*f_new ,"h_p_purity" );
TH1D* h_Purity_clone = (TH1D*)h4->Clone(uniq());
TH1D* h4_clone = (TH1D*)h3->Clone(uniq());


TH1D* h3_input = GetTH1DHist(*f_new ,"h_p" );
TH1D* h4_input = GetTH1DHist(*f_new ,"h_p_purity" );

TH1D* Dom_eff_times_pur_2 = (TH1D*)h3->Clone("Dom_eff_times_pur_2");
TH1D* Nom_eff_times_pur_2 = (TH1D*)h4->Clone("Nom_eff_times_pur_2");

Dom_eff_times_pur_1->Multiply(Dom_eff_times_pur_2);

Nom_eff_times_pur_1->Multiply(Nom_eff_times_pur_2);

h_Purity_clone->Divide(h_Purity_clone,h4_clone);




TEfficiency* pur = new TEfficiency(*h4, *h3);
pur->SetTitle(";Muon momentum p_{#mu} [GeV/c]; Purity");
pur->Draw();
//c2->SaveAs("cc0pi_pur.pdf");
c2 -> Print(text_title_pdf2);
//////////////////////////////////////////////
//
/////////////////////////////////////////////
//SetMarkerColorTEfficiency(eff, kRed);
//eff->SetMarkerColor(kRed);
//pur->SetFillColor(kBlue);
//eff->SetMarkerColor(4);
//eff->SetFillColor(4);
//SetMarkerColorTEfficiency(pur, 4);
    // Access the histograms of the TEfficiency objects
   //TH1D *hist1_Dom = (TH1D*)eff->GetCopyTotalHisto();
   //TH1D *hist2_Dom = (TH1D*)pur->GetCopyTotalHisto();
//
   //TH1D *hist1_nom = (TH1D*)eff->GetCopyPassedHisto();
   //TH1D *hist2_nom = (TH1D*)pur->GetCopyPassedHisto();
//
   // // Multiply the histograms bin by bin
   // hist1_Dom->Multiply(hist1_Dom,hist2_Dom);
   //hist1_nom->Multiply(hist1_nom,hist2_nom);
   
 // Access the underlying TGraphAsymmErrors
    TGraphAsymmErrors *graph2 = pur->GetPaintedGraph();

    // Create a new TGraphAsymmErrors with the same data
    TGraphAsymmErrors *newGraph2 = new TGraphAsymmErrors(*graph2);

    // Set the marker color for the new graph
 
   newGraph2->SetLineColor(kRed);


 // Access the underlying TGraphAsymmErrors
    TGraphAsymmErrors *graph = eff->GetPaintedGraph();

    // Create a new TGraphAsymmErrors with the same data
    TGraphAsymmErrors *newGraph = new TGraphAsymmErrors(*graph);

    // Set the marker color for the new graph
    newGraph->SetLineColor(kBlue);
    std::cout<<"here"<< std::endl;
    

    
   TEfficiency* purtimesEff_Pmu = new TEfficiency(*Nom_eff_times_pur_1, *Dom_eff_times_pur_1);
   
   /*
   
     //TGraphAsymmErrors *newGraph3 = new TGraphAsymmErrors(*graph3);
    std::cout<<"here2"<< std::endl;
  TGraphAsymmErrors *graph3 = purtimesEff_Pmu->GetPaintedGraph();
      std::cout<<"here3"<< std::endl;
      
      if (graph3) {
    
 std::cout<<"graph3 okay"<<std::endl;
} else {
std::cout<<"graph3 is null"<<std::endl;
    // handle the case where graph3 is null
}
      

    std::cout<<"here4"<< std::endl;
    */
    // Draw the new graph
 // "AP" option to draw markers
    //Nom_eff_times_pur_1->Divide(Dom_eff_times_pur_1);
    newGraph2->GetYaxis()->SetTitle("Efficiency*Purity");
    newGraph2->SetMaximum(1.2);
    newGraph2->SetMinimum(0.0);
    c2->SetGrid();
    newGraph2->Draw("ALP");  
    newGraph->Draw("LP SAME");
    purtimesEff_Pmu->Draw("LP SAME");
    //Nom_eff_times_pur_1->Draw("Hist");
    //Dom_eff_times_pur_1->Draw("Same");
   //newGraph3->Draw("LP same"); 

TLegend* Legend4 = new TLegend(0.25, 0.75, 0.35, 0.88);
  Legend4->SetNColumns(1);
  Legend4->SetBorderSize(0);
  Legend4->SetTextSize(.03); //
  Legend4->AddEntry(newGraph2, "Purity", "l" );
  Legend4->AddEntry(newGraph, "Eff", "l" );
  Legend4->AddEntry(purtimesEff_Pmu, "Eff*Purity", "l" );
  Legend4->Draw("same");


//eff->Draw("Same");
c2 -> Print(text_title_pdf2);

delete pur;
delete eff;
delete purtimesEff_Pmu;
///////////////////////////////////////////////////////////
//DrawEfficiency(h_input,
//h2_input, 
//h3_input,
//h4_input, 
//"f",
//"f", 
//c2, text_title_pdf2_string );

//////////////////////////////////////////


h_Eff_clone->SetTitle(";Muon momentum p_{#mu} [GeV/c];Efficiency*Purity");
h_Eff_clone->SetMarkerStyle(21);
h_Eff_clone->Draw("pl");
h_Eff_clone->SetMaximum(1.3);
TH1F* h_Efftimespurity_clone = (TH1F*)h_Eff_clone->Clone("h_Efftimespurity_clone");
h_Purity_clone->SetMarkerStyle(4);
h_Purity_clone->Draw("P L same");

h_Efftimespurity_clone->Multiply(h_Efftimespurity_clone,h_Purity_clone);
h_Efftimespurity_clone->SetMarkerStyle(22);
h_Efftimespurity_clone->Draw("P L Same");
   gStyle->SetOptStat(0);
TLegend* Legend = new TLegend( 0.3, 0.7, 0.6, 0.85 );
  Legend->SetNColumns(1);
  Legend->SetBorderSize(0);
  Legend->SetTextSize(.03); //
  Legend->AddEntry(h_Eff_clone, "Eff", "p" );
  Legend->AddEntry(h_Purity_clone, "Purity", "p" );
  Legend->AddEntry(h_Efftimespurity_clone, "Eff*Purity", "p" );
  Legend->Draw("same");

c2 -> Print(text_title_pdf2);


 
//////////////////////////////////////////////////////////////////////////////
/// Muon BDT Prediction 
//////////////////////////////////////////////////////////////////////////////
/*
 TH1D* h_Dem_eff_BDT =  GetTH1DHist(*f_new ,"h_CCZeroP_MuonCandidate_BTDPrediction_TRUE" );
 TH1D* h_Nom_eff_BDT =  GetTH1DHist(*f_new ,"h_CCZeroP_MuonCandidate_BTDPrediction_TRUE_RECO" );

TH1D* h_Dem_eff_BDT_Clone = (TH1D*)h_Dem_eff_BDT->Clone(uniq());
TH1D* h_Nom_eff_BDT_Clone = (TH1D*)h_Nom_eff_BDT->Clone(uniq());
 
 TEfficiency* eff_BDT = new TEfficiency(*h_Nom_eff_BDT, *h_Dem_eff_BDT);

 eff_BDT->SetTitle("BDT Predicted Muon Candidate;BDT SCORE; Efficiency");
 eff_BDT->Draw();
 c2 -> Print(text_title_pdf2);
 
 
 TH1D* h_Dem_purity_BDT =  GetTH1DHist(*f_new ,"h_CCZeroP_MuonCandidate_BTDPrediction" );
 TH1D* h_Nom_purity_BDT =  GetTH1DHist(*f_new ,"h_CCZeroP_MuonCandidate_BTDPrediction_purity" );
 
 TH1D* h_Dem_purity_BDT_Clone = (TH1D*)h_Dem_purity_BDT->Clone(uniq());
 TH1D* h_Nom_purity_BDT_Clone = (TH1D*)h_Nom_purity_BDT->Clone(uniq());
 
 h_Dem_purity_BDT_Clone->Multiply(h_Dem_eff_BDT_Clone);
 h_Nom_purity_BDT_Clone->Multiply(h_Nom_eff_BDT_Clone);
 
 
 
 
 TEfficiency* purity_BDT = new TEfficiency(*h_Nom_purity_BDT, *h_Dem_purity_BDT);
 purity_BDT->SetTitle("BDT Predicted Muon Candidate; BDT SCORE; Purity");
 purity_BDT->Draw();
 c2 -> Print(text_title_pdf2);
 
h_Nom_eff_BDT->Divide(h_Nom_eff_BDT,h_Dem_eff_BDT);
h_Nom_eff_BDT->SetMarkerStyle(21);
h_Nom_eff_BDT->SetMaximum(1.2);
h_Nom_eff_BDT->SetMinimum(0.0);
h_Nom_eff_BDT->SetTitle("BDT Predicted Muon Candidate");
h_Nom_eff_BDT->GetYaxis()->SetTitle("Efficiency*Purity");
h_Nom_eff_BDT->SetLineStyle(1);
h_Nom_eff_BDT->SetLineWidth(2);
h_Nom_eff_BDT->Draw("pl HIST");
h_Nom_eff_BDT->GetXaxis()->SetTitle("BDT score");

h_Nom_purity_BDT->Divide(h_Nom_purity_BDT,h_Dem_purity_BDT);
h_Nom_purity_BDT->SetMarkerStyle(4);
h_Nom_purity_BDT->SetLineStyle(1);
h_Nom_purity_BDT->SetLineWidth(2);
h_Nom_purity_BDT->Draw("pl HIST SAME");
TH1D* h_Efftimespurity_BDT_clone = (TH1D*)h_Nom_eff_BDT->Clone("h_Efftimespurity_BDT_clone");


h_Efftimespurity_BDT_clone->Multiply(h_Efftimespurity_BDT_clone,h_Nom_purity_BDT);
h_Efftimespurity_BDT_clone->SetMarkerStyle(22);
h_Efftimespurity_BDT_clone->SetLineStyle(2);
h_Efftimespurity_BDT_clone->SetLineWidth(2);
h_Efftimespurity_BDT_clone->Draw("pl HIST SAME");

//AddCutArrow(.14, 0.0, 0.8, .06, false,4,1,8);
	char textplace[1024];
  //  TLatex* text = new TLatex;
  //    text->SetNDC();
  //    text->SetTextSize(0.03);
  //    text->SetTextColor(kRed);
  //    sprintf(textplace, "Cut Value =  %.2f",.14 );
  //    text->DrawLatex(0.12, 0.66, textplace);



TLegend* Legend2 = new TLegend( 0.2, 0.7, 0.5, 0.85 );
  Legend2->SetNColumns(1);
  Legend2->SetBorderSize(0);
  Legend2->SetTextSize(.03); //
  Legend2->AddEntry(h_Nom_eff_BDT, "Eff", "pl" );
  Legend2->AddEntry(h_Nom_purity_BDT, "Purity", "pl" );
  Legend2->AddEntry(h_Efftimespurity_BDT_clone, "Eff*Purity", "pl" );
  Legend2->Draw("same");
c2 -> Print(text_title_pdf2);
delete Legend2;
/////////////////////////////////////////////////////////////////////

    TGraphAsymmErrors *graph4 = purity_BDT->GetPaintedGraph();
    // Create a new TGraphAsymmErrors with the same data
    TGraphAsymmErrors *newGraph4 = new TGraphAsymmErrors(*graph4);
     newGraph4->SetLineColor(kRed);


    TGraphAsymmErrors *graph5 = eff_BDT->GetPaintedGraph();
    // Create a new TGraphAsymmErrors with the same data
    TGraphAsymmErrors *newGraph5 = new TGraphAsymmErrors(*graph5);
     newGraph5->SetLineColor(kBlue);

   TEfficiency* purtimesEff_Pmu_BDT_muon = new TEfficiency(*h_Nom_purity_BDT_Clone, *h_Dem_purity_BDT_Clone);



    newGraph4->GetYaxis()->SetTitle("Efficiency*Purity");
    newGraph4->SetMaximum(1.2);
    newGraph4->SetMinimum(0.0);
    c2->SetGrid();
    newGraph4->Draw("ALP");  
    newGraph5->Draw("LP SAME");
    purtimesEff_Pmu_BDT_muon->Draw("LP SAME");

//AddCutArrow(.14, 0.0, 0.8, .06, false,4,1,8);
  //    text->SetNDC();
  //    text->SetTextSize(0.03);
  //    text->SetTextColor(kRed);
  //    sprintf(textplace, "Cut Value =  %.2f",.14 );
  //    text->DrawLatex(0.12, 0.66, textplace);



TLegend* Legend5 = new TLegend(0.25, 0.75, 0.35, 0.88);
  Legend5->SetNColumns(1);
  Legend5->SetBorderSize(0);
  Legend5->SetTextSize(.03); //
  Legend5->AddEntry(newGraph4, "Purity", "l" );
  Legend5->AddEntry(newGraph5, "Eff", "l" );
  Legend5->AddEntry(purtimesEff_Pmu_BDT_muon, "Eff*Purity", "l" );
  Legend5->Draw("same");
 
  c2 -> Print(text_title_pdf2);


  delete purtimesEff_Pmu_BDT_muon;
  delete purity_BDT;
  delete eff_BDT;
  delete Legend5;

*/
/////////////////////////////////////////////////////////////////////
// Proton BDT predction 
/////////////////////////////////////////////////////////////////////
/*

 TH1D* h_Dem_eff_BDT_Proton =  GetTH1DHist(*f_new ,"h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE" );
 TH1D* h_Nom_eff_BDT_Proton =  GetTH1DHist(*f_new ,"h_CCZeroP_ALLProtonCandidate_BTDPrediction_TRUE_RECO" );
 
 TH1D* h_Dem_eff_BDT_Proton_Clone = (TH1D*) h_Dem_eff_BDT_Proton->Clone(uniq());
 TH1D* h_Nom_eff_BDT_Proton_Clone = (TH1D*) h_Nom_eff_BDT_Proton->Clone(uniq());
 
 
 
 TEfficiency* eff_BDT_Proton = new TEfficiency(*h_Nom_eff_BDT_Proton, *h_Dem_eff_BDT_Proton);
 
 eff_BDT_Proton->SetTitle("BTD Predicted Proton Candidates;BDT SCORE; Efficiency");
 eff_BDT_Proton->Draw();
 c2 -> Print(text_title_pdf2);


 TH1D* h_Dem_purity_BDT_Proton =  GetTH1DHist(*f_new ,"h_CCZeroP_ALLProtonCandidate_BTDPrediction" );
 TH1D* h_Nom_purity_BDT_Proton =  GetTH1DHist(*f_new ,"h_CCZeroP_ALLProtonCandidate_BTDPrediction_Purity" );
 
 TH1D* h_Dem_purity_BDT_Proton_Clone = (TH1D*) h_Dem_purity_BDT_Proton->Clone(uniq());
 TH1D* h_Nom_purity_BDT_Proton_Clone = (TH1D*) h_Nom_purity_BDT_Proton->Clone(uniq());
 
TEfficiency* purity_BDT_Proton = new TEfficiency(*h_Nom_purity_BDT_Proton, *h_Dem_purity_BDT_Proton);
 purity_BDT_Proton->SetTitle("BTD Predicted Proton Candidates; BDT SCORE; Purity");
 purity_BDT_Proton->Draw();
 c2 -> Print(text_title_pdf2);

h_Nom_eff_BDT_Proton->Divide(h_Nom_eff_BDT_Proton,h_Dem_eff_BDT_Proton);
h_Nom_eff_BDT_Proton->SetMarkerStyle(21);
h_Nom_eff_BDT_Proton->SetMaximum(1.2);
h_Nom_eff_BDT_Proton->SetMinimum(0.0);
h_Nom_eff_BDT_Proton->SetTitle("BTD Predicted Proton Candidates");
h_Nom_eff_BDT_Proton->GetYaxis()->SetTitle("Efficiency*Purity");
h_Nom_eff_BDT_Proton->SetLineStyle(1);
h_Nom_eff_BDT_Proton->SetLineWidth(2);
h_Nom_eff_BDT_Proton->Draw("pl HIST");
h_Nom_eff_BDT_Proton->GetXaxis()->SetTitle("BDT score");

h_Nom_purity_BDT_Proton->Divide(h_Nom_purity_BDT_Proton,h_Dem_purity_BDT_Proton);
h_Nom_purity_BDT_Proton->SetMarkerStyle(4);
h_Nom_purity_BDT_Proton->SetLineStyle(1);
h_Nom_purity_BDT_Proton->SetLineWidth(2);
h_Nom_purity_BDT_Proton->Draw("pl HIST SAME");
TH1D* h_Efftimespurity_BDT_Proton_clone = (TH1D*)h_Nom_eff_BDT_Proton->Clone("h_Efftimespurity_BDT_Proton_clone");


h_Efftimespurity_BDT_Proton_clone->Multiply(h_Efftimespurity_BDT_Proton_clone,h_Nom_purity_BDT_Proton);
h_Efftimespurity_BDT_Proton_clone->SetMarkerStyle(22);
h_Efftimespurity_BDT_Proton_clone->SetLineStyle(2);
h_Efftimespurity_BDT_Proton_clone->SetLineWidth(2);
h_Efftimespurity_BDT_Proton_clone->Draw("pl HIST SAME");

TLegend* Legend3 = new TLegend( 0.15, 0.7, 0.4, 0.85 );
  Legend3->SetNColumns(1);
  Legend3->SetBorderSize(0);
  Legend3->SetTextSize(.03); //
  Legend3->AddEntry(h_Nom_eff_BDT_Proton, "Eff", "pl" );
  Legend3->AddEntry(h_Nom_purity_BDT_Proton, "Purity", "pl" );
  Legend3->AddEntry(h_Efftimespurity_BDT_Proton_clone, "Eff*Purity", "pl" );
  Legend3->Draw("same");

//AddCutArrow(.28, 0.0, 0.8, .06, false,4,1,8);

  //  TLatex* text1 = new TLatex;
  //    text1->SetNDC();
  //    text1->SetTextSize(0.03);
  //    text1->SetTextColor(kRed);
  //    sprintf(textplace, "Cut Value =  %.2f",.28 );
  //    text1->DrawLatex(0.24, 0.66, textplace);
       c2 -> Print(text_title_pdf2);
       delete Legend3;
//////////////////////////////////////////////////////////
    TGraphAsymmErrors *graph6 = purity_BDT_Proton->GetPaintedGraph();
    // Create a new TGraphAsymmErrors with the same data
    TGraphAsymmErrors *newGraph6 = new TGraphAsymmErrors(*graph6);
     newGraph6->SetLineColor(kRed);

    TGraphAsymmErrors *graph7 = eff_BDT_Proton->GetPaintedGraph();
    // Create a new TGraphAsymmErrors with the same data
    TGraphAsymmErrors *newGraph7 = new TGraphAsymmErrors(*graph7);
     newGraph7->SetLineColor(kBlue);
     
     
     h_Dem_eff_BDT_Proton_Clone->Multiply(h_Dem_purity_BDT_Proton_Clone);
     h_Nom_eff_BDT_Proton_Clone->Multiply(h_Nom_purity_BDT_Proton_Clone);
     TEfficiency* purtimesEff_Pmu_BDT_Proton = new TEfficiency(*h_Nom_eff_BDT_Proton_Clone, *h_Dem_eff_BDT_Proton_Clone);

    newGraph6->GetYaxis()->SetTitle("Efficiency*Purity");
    newGraph6->SetMaximum(1.2);
    newGraph6->SetMinimum(0.0);
    c2->SetGrid();
    newGraph6->Draw("ALP");  
    newGraph7->Draw("LP SAME");
    purtimesEff_Pmu_BDT_Proton->Draw("LP SAME");

   //AddCutArrow(.28, 0.0, 0.8, .06, false,4,1,8);
      
      //text1->SetNDC();
      //text1->SetTextSize(0.03);
      //text1->SetTextColor(kRed);
      //sprintf(textplace, "Cut Value =  %.2f",.28 );
     // text1->DrawLatex(0.24, 0.66, textplace);



TLegend* Legend6 = new TLegend(0.2, 0.75, 0.3, 0.88);
  Legend6->SetNColumns(1);
  Legend6->SetBorderSize(0);
  Legend6->SetTextSize(.03); //
  Legend6->AddEntry(newGraph6, "Purity", "l" );
  Legend6->AddEntry(newGraph7, "Eff", "l" );
  Legend6->AddEntry(purtimesEff_Pmu_BDT_Proton, "Eff*Purity", "l" );
  Legend6->Draw("same");
c2 -> Print(text_title_pdf2);

  delete purtimesEff_Pmu_BDT_Proton;
  delete purity_BDT_Proton;
  delete eff_BDT_Proton;
  delete Legend6;
  
  */
///////////////////////////////////////////////////
// Angle
/////////////////////////////////////////////////
 TH1D* h_Dem_eff_Angle =  GetTH1DHist(*f_new ,"h_costheta_p_TRUE" );
 TH1D* h_Nom_eff_Angle =  GetTH1DHist(*f_new ,"h_costheta_p_TRUE_RECO" );
 
  TH1D* h_Dem_eff_Angle_Clone = (TH1D*) h_Dem_eff_Angle->Clone(uniq());
  TH1D* h_Nom_eff_Angle_Clone = (TH1D*) h_Nom_eff_Angle->Clone(uniq());
 
  TH1D* h_Dem_pur_Angle =  GetTH1DHist(*f_new ,"h_costheta_p" );
  TH1D* h_Nom_pur_Angle =  GetTH1DHist(*f_new ,"h_costheta_p_purity" );
 
  TH1D* h_Dem_pur_Angle_Clone = (TH1D*) h_Dem_pur_Angle->Clone(uniq());
  TH1D* h_Nom_pur_Angle_Clone = (TH1D*) h_Nom_pur_Angle->Clone(uniq());
 
 h_Dem_pur_Angle_Clone->Multiply(h_Dem_eff_Angle_Clone);
 h_Nom_pur_Angle_Clone->Multiply(h_Nom_eff_Angle_Clone);
 
 TEfficiency* eff_Angle = new TEfficiency(*h_Nom_eff_Angle, *h_Dem_eff_Angle);
 TEfficiency* pure_Angle = new TEfficiency(*h_Nom_pur_Angle, *h_Dem_pur_Angle);
TEfficiency* puretimeseff_Angle = new TEfficiency(*h_Nom_pur_Angle_Clone, *h_Dem_pur_Angle_Clone);
  pure_Angle->SetTitle("; Cos(#theta_{#mu}); Purity");
  pure_Angle->Draw();
  c2 -> Print(text_title_pdf2);
  
  eff_Angle->SetTitle("; Cos(#theta_{#mu}); Efficiency");
  eff_Angle->Draw();
  c2 -> Print(text_title_pdf2);
  //////
  TGraphAsymmErrors *graph8 = pure_Angle->GetPaintedGraph();
    TGraphAsymmErrors *newGraph8 = new TGraphAsymmErrors(*graph8);
     newGraph8->SetLineColor(kRed);
     
     
      TGraphAsymmErrors *graph9 = eff_Angle->GetPaintedGraph();
    TGraphAsymmErrors *newGraph9 = new TGraphAsymmErrors(*graph9);
     newGraph9->SetLineColor(kBlue); 
     

    newGraph8->GetYaxis()->SetTitle("Efficiency*Purity");
    newGraph8->SetMinimum(0.0);
    newGraph8->SetMaximum(1.2);
//    c2->SetGrid();
   newGraph8->Draw("ALP");  
   newGraph9->Draw("LP SAME");
  puretimeseff_Angle->Draw("LP SAME");
  
  
  TLegend* Legend7 = new TLegend(0.25, 0.75, 0.35, 0.88);
  Legend7->SetNColumns(1);
  Legend7->SetBorderSize(0);
  Legend7->SetTextSize(.03); //
  Legend7->AddEntry(newGraph8, "Purity", "l" );
  Legend7->AddEntry(newGraph9, "Eff", "l" );
  Legend7->AddEntry(puretimeseff_Angle, "Eff*Purity", "l" );
  Legend7->Draw("same");
  
  
  c2 -> Print(text_title_pdf2);
  
 delete eff_Angle;
 delete pure_Angle;
 delete puretimeseff_Angle;
 delete Legend7;
///////////////////////////////////////////////////
 TH2D* h_Dem_eff_Pmu_costheta =  GetTH2DHist(*f_new ,"h_costheta_Pmu_TRUE" );
 TH2D* h_Nom_eff_Pmu_costheta =  GetTH2DHist(*f_new ,"h_costheta_Pmu_TRUE_RECO" );


h_Nom_eff_Pmu_costheta->Divide(h_Dem_eff_Pmu_costheta); 
Draw_heatMap(
  h_Nom_eff_Pmu_costheta, 
  "Cos(#theta_{#mu})",
  "Pmu [GeV/c]",
  "2D Efficiency",
 text_title_pdf2,
  0,
  c2,
 false );


 TH2D* h_Dem_pure_Pmu_costheta =  GetTH2DHist(*f_new ,"h_costheta_Pmu" );
 TH2D* h_Nom_pure_Pmu_costheta =  GetTH2DHist(*f_new ,"h_costheta_Pmu_Purity" );


h_Nom_pure_Pmu_costheta->Divide(h_Dem_pure_Pmu_costheta); 


Draw_heatMap(
 h_Nom_pure_Pmu_costheta, 
  "Cos(#theta_{#mu})",
  "Pmu [GeV/c]",
  "2D Purity",
 text_title_pdf2,
  0,
  c2,
 false );

h_Nom_eff_Pmu_costheta->Multiply(h_Nom_pure_Pmu_costheta);

Draw_heatMap(
 h_Nom_eff_Pmu_costheta, 
  "Cos(#theta_{#mu})",
  "Pmu [GeV/c]",
  "2D Efficiency*Purity",
 text_title_pdf2,
  0,
  c2,
 false );
 ////////////////////////////////////////////////////
 // Proton Mulit
 //////////////////////////////////////////////////
 
  TH1D* h_Dem_eff_ProtonMult =  GetTH1DHist(*f_new ,"h_Proton_mult_TRUE" );
  TH1D* h_Nom_eff_ProtonMult =  GetTH1DHist(*f_new ,"h_Proton_mult_TRUE_RECO" );
 
  TH1D* h_Dem_eff_ProtonMult_Clone = (TH1D*) h_Dem_eff_ProtonMult->Clone(uniq());
  TH1D* h_Nom_eff_ProtonMult_Clone = (TH1D*) h_Nom_eff_ProtonMult->Clone(uniq());
 
  TH1D* h_Dem_pur_ProtonMult =  GetTH1DHist(*f_new ,"h_Proton_mult" );
  TH1D* h_Nom_pur_ProtonMult =  GetTH1DHist(*f_new ,"h_Proton_mult_Purity" );
 
  TH1D* h_Dem_pur_ProtonMult_Clone = (TH1D*) h_Dem_pur_ProtonMult->Clone(uniq());
  TH1D* h_Nom_pur_ProtonMult_Clone = (TH1D*) h_Nom_pur_ProtonMult->Clone(uniq());
 
 h_Dem_pur_ProtonMult_Clone->Multiply(h_Dem_eff_ProtonMult_Clone);
 h_Nom_pur_ProtonMult_Clone->Multiply(h_Nom_eff_ProtonMult_Clone);
 
 TEfficiency* eff_ProtonMult = new TEfficiency(*h_Nom_eff_ProtonMult, *h_Dem_eff_ProtonMult);
 TEfficiency* pure_ProtonMult = new TEfficiency(*h_Nom_pur_ProtonMult, *h_Dem_pur_ProtonMult);
 TEfficiency* puretimeseff_ProtonMult = new TEfficiency(*h_Nom_pur_ProtonMult_Clone, *h_Dem_pur_ProtonMult_Clone);
  pure_ProtonMult->SetTitle("; Proton Multiplicity [N]; Purity");
  pure_ProtonMult->Draw();
  c2 -> Print(text_title_pdf2);
  eff_ProtonMult->SetTitle("; Proton Multiplicity [N]; Efficiency");
  eff_ProtonMult->Draw();
  c2 -> Print(text_title_pdf2);
  
  
    TGraphAsymmErrors *graph10 = pure_ProtonMult->GetPaintedGraph();
    TGraphAsymmErrors *newGraph10 = new TGraphAsymmErrors(*graph10);
     newGraph10->SetLineColor(kRed);
     
     
      TGraphAsymmErrors *graph11 = eff_ProtonMult->GetPaintedGraph();
    TGraphAsymmErrors *newGraph11 = new TGraphAsymmErrors(*graph11);
     newGraph11->SetLineColor(kBlue); 
     

    newGraph10->GetYaxis()->SetTitle("Eff#timesPurity");
    newGraph10->GetYaxis()->SetTitleSize(0.05); 
    newGraph10->GetYaxis()->SetTitleOffset(0.9); 
    newGraph10->SetMinimum(0.0);
    newGraph10->SetMaximum(1.3);
   newGraph10->SetLineWidth(2);
   newGraph10->Draw("ALP");  
   newGraph11->Draw("LP SAME");
   newGraph11->SetLineWidth(2);
   puretimeseff_ProtonMult->Draw("LP SAME");
   puretimeseff_ProtonMult->SetLineWidth(2);
  TLegend* Legend8 = new TLegend(0.25, 0.75, 0.45, 0.88);
  Legend8->SetNColumns(1);
  Legend8->SetBorderSize(0);
  Legend8->SetTextSize(.035); //
  Legend8->AddEntry(newGraph10, "Purity", "l" );
  Legend8->AddEntry(newGraph11, "Efficiency", "l" );
  Legend8->AddEntry(puretimeseff_ProtonMult, "Eff*Purity", "l" );
  Legend8->Draw("same");
  c2 -> Print(text_title_pdf2);
  
  
delete eff_ProtonMult;
delete pure_ProtonMult; 
delete puretimeseff_ProtonMult;
 ////////////////////////////////////////////////////
 // Proton eff per track 
 //////////////////////////////////////////////////
 
  TH1D* h_Dem_eff_Cosmic =  GetTH1DHist(*f_new ,"h_Cosmic_Cut_TRUE" );
  TH1D* h_Nom_eff_Cosmic =  GetTH1DHist(*f_new ,"h_Cosmic_Cut_TRUE_RECO" );
 
  TH1D* h_Dem_eff_Cosmic_Clone = (TH1D*) h_Dem_eff_Cosmic->Clone(uniq());
  TH1D* h_Nom_eff_Cosmic_Clone = (TH1D*) h_Nom_eff_Cosmic->Clone(uniq());
 
 
  TH1D* h_Dem_pur_Cosmic =  GetTH1DHist(*f_new ,"h_Cosmic_Cut" );
  TH1D* h_Nom_pur_Cosmic =  GetTH1DHist(*f_new ,"h_Cosmic_Cut_Purity" );
 
  TH1D* h_Dem_pur_Cosmic_Clone = (TH1D*) h_Dem_pur_Cosmic->Clone(uniq());
  TH1D* h_Nom_pur_Cosmic_Clone = (TH1D*) h_Nom_pur_Cosmic->Clone(uniq());
 
 h_Dem_pur_Cosmic_Clone->Multiply(h_Dem_eff_Cosmic_Clone);
 h_Nom_eff_Cosmic_Clone->Multiply(h_Nom_eff_Cosmic_Clone);
 
 TEfficiency* eff_Cosmics = new TEfficiency(*h_Nom_eff_Cosmic, *h_Dem_eff_Cosmic);
 TEfficiency* pure_Cosmics = new TEfficiency(*h_Nom_pur_Cosmic, *h_Dem_pur_Cosmic);
 TEfficiency* puretimeseff_Cosmics = new TEfficiency(*h_Nom_pur_Cosmic_Clone, *h_Dem_pur_Cosmic_Clone);
  pure_Cosmics->SetTitle("Cosmics; Cosmics [cm]; Purity");
  pure_Cosmics->Draw();
  c2 -> Print(text_title_pdf2);
  eff_Cosmics->SetTitle("Cosmics; Cosmics [cm]; Efficiency");
  eff_Cosmics->Draw();
  c2 -> Print(text_title_pdf2);
  
  
    TGraphAsymmErrors *graph12 = pure_Cosmics->GetPaintedGraph();
    TGraphAsymmErrors *newGraph12 = new TGraphAsymmErrors(*graph12);
     newGraph12->SetLineColor(kRed);
     
     
      TGraphAsymmErrors *graph13 = eff_Cosmics->GetPaintedGraph();
    TGraphAsymmErrors *newGraph13 = new TGraphAsymmErrors(*graph13);
     newGraph13->SetLineColor(kBlue); 
     

    newGraph13->GetYaxis()->SetTitle("Efficiency*Purity");
    newGraph12->SetMinimum(0.0);
   newGraph12->SetMaximum(1.0);
   newGraph12->Draw("ALP");  
   newGraph13->Draw("LP SAME");
   puretimeseff_Cosmics->Draw("LP SAME");
   
  TLegend* Legend9 = new TLegend(0.25, 0.75, 0.35, 0.88);
  Legend9->SetNColumns(1);
  Legend9->SetBorderSize(0);
  Legend9->SetTextSize(.035); //
  Legend9->AddEntry(newGraph12, "Purity", "l" );
  Legend9->AddEntry(newGraph13, "Eff", "l" );
  Legend9->AddEntry(puretimeseff_Cosmics, "Eff*Purity", "l" );
  Legend9->Draw("same");
  c2 -> Print(text_title_pdf2);
  
  
delete eff_Cosmics;
delete pure_Cosmics; 
delete puretimeseff_Cosmics;

 
 ////////////////////////////////////////////////////
 // Proton eff per track 
 //////////////////////////////////////////////////
 
  TH1D* h_Dem_eff_trk_len_v =  GetTH1DHist(*f_new ,"h_trk_len_v_TRUE" );
  TH1D* h_Nom_eff_trk_len_v =  GetTH1DHist(*f_new ,"h_trk_len_v_TRUE_RECO" );
 
  TH1D* h_Dem_eff_trk_len_v_Clone = (TH1D*) h_Dem_eff_trk_len_v->Clone(uniq());
  TH1D* h_Nom_eff_trk_len_v_Clone = (TH1D*) h_Nom_eff_trk_len_v->Clone(uniq());
 
 
  TH1D* h_Dem_pur_trk_len_v =  GetTH1DHist(*f_new ,"h_trk_len_v" );
  TH1D* h_Nom_pur_trk_len_v =  GetTH1DHist(*f_new ,"h_trk_len_v_Purity" );
 
  TH1D* h_Dem_pur_trk_len_v_Clone = (TH1D*) h_Dem_pur_trk_len_v->Clone(uniq());
  TH1D* h_Nom_pur_trk_len_v_Clone = (TH1D*) h_Nom_pur_trk_len_v->Clone(uniq());
 
 h_Dem_pur_trk_len_v_Clone->Multiply(h_Dem_eff_trk_len_v_Clone);
 h_Nom_eff_trk_len_v_Clone->Multiply(h_Nom_eff_trk_len_v_Clone);
 
 TEfficiency* eff_trk_len_v = new TEfficiency(*h_Nom_eff_trk_len_v, *h_Dem_eff_trk_len_v);
 TEfficiency* pure_trk_len_v = new TEfficiency(*h_Nom_pur_trk_len_v, *h_Dem_pur_trk_len_v);
 TEfficiency* puretimeseff_trk_len_v = new TEfficiency(*h_Nom_pur_trk_len_v_Clone, *h_Dem_pur_trk_len_v_Clone);
  pure_trk_len_v->SetTitle("TrkLength; TrkLength [cm]; Purity");
  pure_trk_len_v->Draw();
  c2 -> Print(text_title_pdf2);
  eff_trk_len_v->SetTitle("TrkLength; TrkLength [cm]; Efficiency");
  eff_trk_len_v->Draw();
  c2 -> Print(text_title_pdf2);
  
  
    TGraphAsymmErrors *graph14 = pure_trk_len_v->GetPaintedGraph();
    TGraphAsymmErrors *newgraph14 = new TGraphAsymmErrors(*graph14);
     newgraph14->SetLineColor(kRed);
     
     
      TGraphAsymmErrors *graph15 = eff_trk_len_v->GetPaintedGraph();
    TGraphAsymmErrors *newGraph15 = new TGraphAsymmErrors(*graph15);
     newGraph15->SetLineColor(kBlue); 
     

    newGraph15->GetYaxis()->SetTitle("Efficiency*Purity");
    newgraph14->SetMinimum(0.0);
   newgraph14->SetMaximum(1.0);
   newgraph14->Draw("ALP");  
   newGraph15->Draw("LP SAME");
   puretimeseff_trk_len_v->Draw("LP SAME");
   
  TLegend* Legend10 = new TLegend(0.25, 0.75, 0.35, 0.88);
  Legend10->SetNColumns(1);
  Legend10->SetBorderSize(0);
  Legend10->SetTextSize(.035); //
  Legend10->AddEntry(newgraph14, "Purity", "l" );
  Legend10->AddEntry(newGraph15, "Eff", "l" );
  Legend10->AddEntry(puretimeseff_trk_len_v, "Eff*Purity", "l" );
  Legend10->Draw("same");
  c2 -> Print(text_title_pdf2);
  
  
delete eff_trk_len_v;
delete pure_trk_len_v; 
delete puretimeseff_trk_len_v;
 
 
 
 
 
///////////////////////////////////////////////////
// Final Page
/////////////////////////////////////////////////
c2 -> Print(text_title_pdf3);


}//////////////end of Function 
/////////////////////////////////////////////

int main() {


std::cout<<" Running  Make_cc0pi_effpur_pmu "<< std::endl;

  //tutorial_slice_plots();
 RunEfficnecy();
  return 0;
}




void SetMarkerColorTEfficiency(TEfficiency *efficiency, Color_t color) {
    if (!efficiency) {
        // Check if the TEfficiency object is valid
        return;
    }

    // Access the underlying TGraphAsymmErrors
    TGraphAsymmErrors *graph = efficiency->GetPaintedGraph();

    // Set the marker color
    graph->SetMarkerColor(color);
}




void DrawEfficiency(TH1D* h_TRUE_input,
TH1D* h_TRUE_RECO_input , 
TH1D* h_RECO_input,
TH1D* h_RECO_PURE_input, 
std::string Xaxis,
std::string Title, 
TCanvas* Can, std::string PDF_name ){

TH1D* h_TRUE = (TH1D*)h_TRUE_input->Clone(uniq()); 
TH1D* h_TRUE_RECO = (TH1D*)h_TRUE_RECO_input->Clone(uniq()); 

TH1D* h_RECO = (TH1D*)h_RECO_input->Clone(uniq()); 
TH1D* h_RECO_PURE = (TH1D*)h_RECO_PURE_input->Clone(uniq());


TH1D* Dom_eff_times_pur_1 = (TH1D*)h_TRUE_input->Clone(uniq());
TH1D* Nom_eff_times_pur_1 = (TH1D*)h_TRUE_RECO_input->Clone(uniq());


TH1D* Dom_eff_times_pur_2 = (TH1D*)h_RECO_input->Clone(uniq());
TH1D* Nom_eff_times_pur_2 = (TH1D*)h_RECO_PURE_input->Clone(uniq());
std::cout<<"here 1"<< std::endl;

Dom_eff_times_pur_1->Multiply(Dom_eff_times_pur_2);
Nom_eff_times_pur_1->Multiply(Nom_eff_times_pur_2);

TEfficiency* purtimesEff_Pmu = new TEfficiency(*Nom_eff_times_pur_1, *Dom_eff_times_pur_1);

TEfficiency* T_Purity = new TEfficiency(*h_RECO_PURE, *h_RECO);
std::cout<<"here 2"<< std::endl;
//TGraphAsymmErrors *graph2 = T_Purity->GetPaintedGraph();
std::cout<<"here 2.1"<< std::endl;

//if (graph2 && graph2->GetN() > 0) {
//    std::cout <<"seems okay"<< std::endl;
//    // Rest of the code
//} else {
//    std::cout << "Error: graph2 is either nullptr or empty." << std::endl;
//    // Handle the case where graph2 is nullptr or empty
//}


//TGraphAsymmErrors *newGraph2 = new TGraphAsymmErrors(*graph2);
std::cout<<"here 2.3"<< std::endl;
//newGraph2->SetLineColor(kBlue);
std::cout<<"here 2.4"<< std::endl;
TEfficiency* T_Effiency = new TEfficiency(*h_TRUE_RECO, *h_TRUE);
std::cout<<"here 2.5"<< std::endl;
TGraphAsymmErrors *graph1 = T_Effiency->GetPaintedGraph();
TGraphAsymmErrors *newGraph1 = new TGraphAsymmErrors(*graph1);
std::cout<<"here 3"<< std::endl;
newGraph1->SetLineColor(kRed);

//newGraph2->Draw("ALP");  
newGraph1->Draw("ALP"); //LP SAME
purtimesEff_Pmu->Draw("LP SAME");

TLegend* Legend = new TLegend( 0.3, 0.7, 0.6, 0.85 );
  Legend->SetNColumns(1);
  Legend->SetBorderSize(0);
  Legend->SetTextSize(.03); //
  //Legend->AddEntry(newGraph2, "Eff", "l" );
  Legend->AddEntry(newGraph1, "Eff", "l" );
  Legend->AddEntry(purtimesEff_Pmu, "Eff*Purity", "l" );
  Legend->Draw("same");

  Can -> Print(PDF_name.c_str());

}/////End OF Function
////////////////////////