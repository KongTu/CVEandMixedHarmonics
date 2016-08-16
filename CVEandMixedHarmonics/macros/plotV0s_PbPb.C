#include "RiceStyle.h"

using namespace std;


void plotV0s_PbPb(){

	gStyle->SetErrorX(0);

	TFile* file[10];
	file[0] = new TFile("../dataPoints/V0s_PbPb.root");

	TGraphErrors* gr1[5];
	for(int i = 0; i < 5; i++){

		gr1[i] = (TGraphErrors*) file[0]->Get(Form("Graph;%d", i+1));
	}

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "Centrality", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 100,0,100,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.09, 0.06);
	base1->GetXaxis()->SetRangeUser(0,100);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.23);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetNdivisions(8,8,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);
	
	TCanvas* c1 = new TCanvas("c1","c1",1,1,650,650);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
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

	gr1[2]->SetMarkerStyle(24);
	gr1[2]->SetMarkerSize(1.4);
	gr1[2]->SetMarkerColor(kRed);
	gr1[2]->SetLineColor(kRed);
	gr1[2]->Draw("Psame");

	gr1[3]->SetMarkerStyle(25);
	gr1[3]->SetMarkerSize(1.4);
	gr1[3]->SetMarkerColor(kBlue);
	gr1[3]->SetLineColor(kBlue);
	gr1[3]->Draw("Psame");

	gr1[4]->SetMarkerStyle(34);
	gr1[4]->SetMarkerSize(1.4);
	gr1[4]->SetMarkerColor(kBlack);
	gr1[4]->SetLineColor(kBlack);
	gr1[4]->Draw("Psame");


}