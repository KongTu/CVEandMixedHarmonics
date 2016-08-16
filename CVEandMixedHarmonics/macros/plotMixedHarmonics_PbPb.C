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

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};

double PbPb_5TeV_ntrkBinCenter[] = {103.921, 134.061, 166.539, 201.615, 239.096, 279.13, 324.044, 374.081, 448.17};
//double PbPb_5TeV_ntrkBinCenter[] = {73.69, 103.921, 134.061, 166.539, 201.615, 239.096, 279.13, 324.044, 374.081, 448.17};

double PbPb_ntrkCentralityBinCenter[] = {625.275, 811.118, 1025.79, 1257.64};

//double PbPb_centralityBinCenter[] = {75, 65, 57.5, 52.5, 47.5, 42.5, 37.5, 32.5};
//double PbPb_centralityBinCenter[] = {75,67.5,62.5,57.5,52.5,47.5,42.5,37.5,32.5};
double PbPb_centralityBinCenter_tracker[] = {65, 55, 45, 35};
double PbPb_centralityBinCenter[] = {35, 45, 55, 65, 75};

const int Nmults = 5;

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;

void plotMixedHarmonics_PbPb(){


	gStyle->SetErrorX(0);

	TFile* file[10];
	file[0] = new TFile("../dataPoints/PbPb_data_centrality.root");

	TGraphErrors* gr1[4];
	for(int i = 0; i < 2; i++){

		gr1[i] = (TGraphErrors*) file[0]->Get(Form("Graph;%d", i+1));
	}


	//start plotting

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "N^{offline}_{trk}", "#LTcos(n_{1}#phi_{#alpha}+n_{2}#phi_{#beta}+n_{3}#phi_{c})#GT/v_{n_{3},c}", 5000,0.1,10000,kBlack);
	TH1D* base2 = makeHist("base2", "p-going", "N^{offline}_{trk}", "#LTcos(n_{1}#phi_{#alpha}+n_{2}#phi_{#beta}+n_{3}#phi_{c})#GT/v_{n_{3},c}", 5000,0.1,10000,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0009, 0.0006);
	base1->GetXaxis()->SetRangeUser(81, 1700);

	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.0009, 0.013);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.23);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetNdivisions(8,18,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	base2->GetYaxis()->SetTitleOffset(1.23);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetTitle("#LTcos(n_{1}#phi_{#alpha}+n_{2}#phi_{#beta}+n_{3}#phi_{c})#GT/v_{n_{3},c} (OS - SS)");
	base3->GetYaxis()->SetRangeUser(-0.0007,0.0016);
	base3->GetYaxis()->SetTitleOffset(1.23);
	base3->GetXaxis()->SetTitleOffset(1.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.2);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.2);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetNdivisions(5,6,0);
	
	TH1D* base4 = (TH1D*) base2->Clone("base4");
	base4->GetYaxis()->SetRangeUser(-0.0006,0.0016);
	base4->GetYaxis()->SetTitleOffset(1.9);
	base4->GetXaxis()->SetTitleOffset(3.1);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.0);
	base4->GetYaxis()->SetNdivisions(6);
	

	TCanvas* c1 = new TCanvas("c1","c1",1,1,650,650);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base10 = makeHist("base10", "", "Centrality (%)", "#LTcos(n_{1}#phi_{#alpha}+n_{2}#phi_{#beta}+n_{3}#phi_{c})#GT/v_{n_{3},c}", 100,0,100,kBlack);

	// Remove the current axis
	
	base10->GetXaxis()->SetLabelOffset(999);
	base10->GetXaxis()->SetTickLength(0);

	// Redraw the new axis
	gPad->Update();
	TGaxis *newaxis1 = new TGaxis(80,
	                            -0.004,
	                            0,
	                            -0.004,
	                            0,
	                            80,
	                            510,"-");
	newaxis1->SetLabelOffset(-0.03);
	newaxis1->SetLabelFont(42);

	base10->GetYaxis()->SetRangeUser(-0.004, 0.0006);
	base10->GetXaxis()->SetRangeUser(0, 80);
	base10->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base10,1.1,1.25);

	base10->GetYaxis()->SetTitleOffset(1.23);
	base10->GetYaxis()->SetTitleSize(base10->GetYaxis()->GetTitleSize()*1.4);
	base10->GetXaxis()->SetTitleSize(base10->GetXaxis()->GetTitleSize()*1.4);
	base10->GetYaxis()->SetLabelSize(base10->GetYaxis()->GetLabelSize()*1.5);
	base10->GetXaxis()->SetLabelSize(base10->GetXaxis()->GetLabelSize()*1.5);
	base10->GetXaxis()->SetNdivisions(8,18,0);
	base10->GetYaxis()->SetNdivisions(4,6,0);

	base10->Draw("");
	newaxis1->Draw("Asame");

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

	TLatex* r41 = new TLatex(0.15, 0.84, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r41->SetNDC();
    r41->SetTextSize(25);
    r41->SetTextFont(43);
    r41->SetTextColor(kBlack);
    r41->Draw("same");

    TLatex* r43 = new TLatex(0.65,0.92, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.05);
    r43->Draw("same");
    
    TLatex* r44 = new TLatex(0.77,0.92, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(24);
    r44->SetTextFont(53);
    r44->Draw("same");
    

    TLegend *w4 = new TLegend(0.5,0.2,0.80,0.35);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(23);
    w4->SetTextFont(45);
    w4->AddEntry(gr1[0], "  SS", "P");
    w4->AddEntry(gr1[1], "  OS", "P");
    w4->Draw("same");

}

