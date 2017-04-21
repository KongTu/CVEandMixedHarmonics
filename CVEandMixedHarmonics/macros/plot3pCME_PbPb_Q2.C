#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {0,35,60,90,120,150,185,220,260};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
int ntrkBinCenter[] = {17.5, 47.5, 75, 105, 135, 167.5, 202.5, 240};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};
double PbPb_5TeV_ntrkBinCenter[] = {103.921, 134.061, 166.539, 201.615, 239.096, 279.13, 324.044, 374.081, 448.17};

double PbPb_ntrkCentralityBinCenter[] = {625.275, 811.118, 1025.79, 1257.64};

double PbPb_centralityBinCenter[] = {57.5, 52.5, 47.5, 42.5, 37.5, 32.5};

const int Nmults = 11;

double total_systematics_pPb = 0.0000;
double total_systematics_PbPb = 0.0000;

void plot3pCME_PbPb_Q2(){

	gStyle->SetErrorX(0);

	TFile* file1[7];
	file1[0] = new TFile("../dataPoints/PbPb_5TeV_cme_q2_cent_3.root");

	TGraphErrors* gr1[4];
	for(int i = 0; i < 4; i++){

		gr1[i] = (TGraphErrors*) file1[0]->Get(Form("Graph;%d", i+1));
	}
    // gr1[2] = (TGraphErrors*) file1[0]->Get(Form("Graph;%d", 5));
    // gr1[3] = (TGraphErrors*) file1[0]->Get(Form("Graph;%d", 6));



//start plotting

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "v_{2,tracker}", "#gamma_{112}", 1000,0,1.0,kBlack);
	base1->GetYaxis()->SetRangeUser(-0.0008, 0.0005);
	base1->GetXaxis()->SetRangeUser(0.05, 0.13);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.6);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.6);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.6);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.6);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	TCanvas* c1 = new TCanvas("c1","c1",1,1,700,700);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	base1->Draw();

	gr1[0]->SetMarkerStyle(20);
	gr1[0]->SetMarkerSize(1.4);
	gr1[0]->SetMarkerColor(kRed);
	gr1[0]->SetLineColor(kRed);
	gr1[0]->Draw("Psame");

	gr1[1]->SetMarkerStyle(21);
	gr1[1]->SetMarkerSize(1.4);
	gr1[1]->SetMarkerColor(kBlue);
	gr1[1]->SetLineColor(kBlue);
	gr1[1]->Draw("Psame");

	TLegend *w1 = new TLegend(0.2,0.33,0.5,0.40);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(18);
    w1->SetTextFont(43);
    w1->SetNColumns(2);

    w1->AddEntry(gr1[0], "  ", "P");
    w1->AddEntry(gr1[1], "  PbPb 5 TeV", "P");
    // w1->AddEntry(gr1[2], "  ", "P");
    // w1->AddEntry(gr1[3], "  #phi_{c}(p-going)", "P");


    w1->Draw("same");

    TLatex* r4 = new TLatex(0.18, 0.78, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);

    TLatex* lmult = new TLatex(0.18, 0.72, "185 #leq N^{offline}_{trk} < 250");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);

	TLatex* latex1 = new TLatex(0.2, 0.42, "SS");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");
    TLatex* latex2 = new TLatex(0.26, 0.42, "OS");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    latex2->Draw("same");


   	TLatex* r11 = new TLatex(0.18,0.84, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);

    TLatex* r22 = new TLatex(0.27,0.84, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
    
    r4->Draw("same");
    lmult->Draw("same");
    r11->Draw("same");
    //r22->Draw("same");

	TH1D* base2 = makeHist("base2", "Pb-going", "v_{2,tracker}", "#gamma_{112} (os-ss)", 1000,0,1.0,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.0004, 0.0015);
	base2->GetXaxis()->SetRangeUser(0.0, 0.16);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.6);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.6);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetNdivisions(4,6,0);
	base2->GetYaxis()->SetNdivisions(4,6,0);


    TCanvas* c3 = new TCanvas("c3","c3",700,700);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

    TGraphErrors* gr_new1 = new TGraphErrors(Nmults);
    TGraphErrors* gr_new2 = new TGraphErrors(Nmults);


	for(int mult = 0; mult < Nmults; mult++){

    	double ye = total_systematics_pPb;

    	double x1;
    	double value1;
    	double value1_error;
    	gr1[0]->GetPoint(mult, x1, value1);
    	value1_error = gr1[0]->GetErrorY(mult);

    	double x2;
    	double value2;
    	double value2_error;
    	gr1[1]->GetPoint(mult, x2, value2);
		value2_error = gr1[1]->GetErrorY(mult);

    	double x = x1;
    	double y = value2 - value1;
    	double ey = sqrt(value1_error*value1_error + value2_error*value2_error);

    	gr_new1->SetPoint(mult, x, y);
    	gr_new1->SetPointError(mult, 0, ey);

    }

	gStyle->SetErrorX(0);

    base2->Draw();
    gr_new1->SetMarkerStyle(34);
    gr_new1->SetMarkerSize(1.5);
    gr_new1->SetLineColor(kRed);
    gr_new1->SetMarkerColor(kRed);

    gr_new1->Fit("pol1");
    TF1 * myFunc1 = gr_new1->GetFunction("pol1");
    myFunc1->SetLineStyle(2);
    double intersect_1 = myFunc1->GetParameter(0);
    double intersect_1_error = myFunc1->GetParError(0);
    double slope_1 = myFunc1->GetParameter(1);
    double slope_1_error = myFunc1->GetParError(1);
    myFunc1->SetRange(0,1);

    TH1D *hint = new TH1D("hint","Fitted gaussian with .95 conf.band", 10000, 0.00, 2);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint, 0.68);
     //Now the "hint" histogram has the fitted function values as the
     //bin contents and the confidence intervals as bin errors
    hint->SetStats(kFALSE);
    hint->SetFillColor(kBlack);
    hint->SetFillStyle(3001);
    hint->Draw("e3 same");

    TLatex* latex3 = new TLatex(0.18, 0.66, Form("slope: %.6f +/- %.6f",slope_1, slope_1_error ));
    latex3->SetNDC();
    latex3->SetTextSize(20);
    latex3->SetTextFont(43);
    latex3->SetTextColor(kRed);
    latex3->Draw("same");

    TLatex* latex4 = new TLatex(0.18, 0.63, Form("intersect: %.6f +/- %.6f",intersect_1, intersect_1_error ));
    latex4->SetNDC();
    latex4->SetTextSize(20);
    latex4->SetTextFont(43);
    latex4->SetTextColor(kRed);
    latex4->Draw("same");

    gr_new1->Draw("Psame");

    r4->Draw("same");
    lmult->Draw("same");
    r11->Draw("same");
    //r22->Draw("same");

    TLegend *w2 = new TLegend(0.6,0.14,0.8,0.24);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(22);
    w2->SetTextFont(43);
    w2->AddEntry(gr_new1, "PbPb", "P");
    //w2->Draw("same");

    TH1D* base3 = makeHist("base3", "Pb-going", "v_{2,tracker}", "#delta (os-ss)", 1000,0,1.0,kBlack);
    base3->GetYaxis()->SetRangeUser(-0.004, 0.015);
    base3->GetXaxis()->SetRangeUser(0.0, 0.16);
    base3->GetXaxis()->SetTitleColor(kBlack);
    
    fixedFontHist1D(base3,1.1,1.25);

    base3->GetYaxis()->SetTitleOffset(1.3);
    base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.6);
    base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.6);
    base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.6);
    base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.6);
    base3->GetXaxis()->SetNdivisions(4,6,0);
    base3->GetYaxis()->SetNdivisions(4,6,0);

    TCanvas* c4 = new TCanvas("c4","c4",700,700);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    gStyle->SetPadBorderMode(0.1);
    gStyle->SetOptTitle(0);

    TGraphErrors* gr_delta = new TGraphErrors(Nmults);

    for(int mult = 0; mult < Nmults; mult++){

        double x1;
        double value1;
        double value1_error;
        gr1[2]->GetPoint(mult, x1, value1);
        value1_error = gr1[2]->GetErrorY(mult);

        double x2;
        double value2;
        double value2_error;
        gr1[3]->GetPoint(mult, x2, value2);
        value2_error = gr1[3]->GetErrorY(mult);

        double x = x1;
        double y = value2 - value1;
        double ey = sqrt(value1_error*value1_error + value2_error*value2_error);

        gr_delta->SetPoint(mult, x, y);
        gr_delta->SetPointError(mult, 0, ey);

    }

    gStyle->SetErrorX(0);
    base3->Draw();
    gr_delta->SetMarkerStyle(34);
    gr_delta->SetMarkerSize(1.5);
    gr_delta->SetLineColor(kRed);
    gr_delta->SetMarkerColor(kRed);
    gr_delta->Draw("Psame");
    r4->Draw("same");

    // TFile f1("../dataPoints/PbPb_ESE_points_5.root", "RECREATE");
    // gr_new1->Write();
    // myFunc1->Write();

    // hint->Write();

    // c1->Print("../figures/CME_figure7.pdf");
    // c3->Print("../figures/CME_figure8.pdf");


}