#include "RiceStyle.h"

using namespace std;


void plotV0s_PbPb(){

	gStyle->SetErrorX(0);

	TFile* file[10];
	file[0] = new TFile("../dataPoints/V0s_PbPb_test.root");

	TGraphErrors* gr1_Pb[3];
	TGraphErrors* gr1_p[3];
	for(int i = 0; i < 3; i++){

		gr1_Pb[i] = (TGraphErrors*) file[0]->Get(Form("la_la_Pb_%d", i));
		gr1_p[i] = (TGraphErrors*) file[0]->Get(Form("la_la_p_%d", i));

	}

	TGraphErrors* gr2_Pb;
	TGraphErrors* gr2_p;

	gr2_Pb = (TGraphErrors*) file[0]->Get("ks_ks_Pb");
	gr2_p = (TGraphErrors*) file[0]->Get("ks_ks_p");

	TGraphErrors* gr3_Pb[2];
	TGraphErrors* gr3_p[2];

	TGraphErrors* gr4_Pb[2];
	TGraphErrors* gr4_p[2];
	for(int i = 0; i < 2; i++){

		gr3_Pb[i] = (TGraphErrors*) file[0]->Get(Form("la_ks_Pb_%d", i));
		gr3_p[i] = (TGraphErrors*) file[0]->Get(Form("la_ks_p_%d", i));

		gr4_Pb[i] = (TGraphErrors*) file[0]->Get(Form("la_h_Pb_%d", i));
		gr4_p[i] = (TGraphErrors*) file[0]->Get(Form("la_h_p_%d", i));

	}

	TGraphErrors* gr5_Pb;
	TGraphErrors* gr5_p;

	gr5_Pb = (TGraphErrors*) file[0]->Get("ks_h_Pb");
	gr5_p = (TGraphErrors*) file[0]->Get("ks_h_p");

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "Centrality", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 100,0,100,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.9, 10.0);
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


	TGraphErrors* gr1_new = new TGraphErrors(6);
	TGraphErrors* gr2_new = new TGraphErrors(6);
	TGraphErrors* gr3_new = new TGraphErrors(6);
	TGraphErrors* gr4_new = new TGraphErrors(6);
	TGraphErrors* gr5_new = new TGraphErrors(6);

	for(int i = 0; i < 6; i++){

		double x1, value1 = 0.0;
		gr1_Pb[0]->GetPoint(i, x1, value1);
		double value1_err = gr1_Pb[0]->GetErrorY(i);
		
		double x2, value2 = 0.0;
		gr1_Pb[1]->GetPoint(i, x2, value2);
		double value2_err = gr1_Pb[1]->GetErrorY(i);

		double x3, value3 = 0.0;
		gr1_p[0]->GetPoint(i, x3, value3);
		double value3_err = gr1_p[0]->GetErrorY(i);

		double x4, value4 = 0.0;
		gr1_p[1]->GetPoint(i, x4, value4);
		double value4_err = gr1_p[1]->GetErrorY(i);

		double total_value = (value1 + value2 + value3 + value4)/4.0;
		double total_error = sqrt(value1_err*value1_err + value2_err*value2_err + value3_err*value3_err + value4_err*value4_err)/4.0;

		gr1_new->SetPoint(i, x1, total_value);
		gr1_new->SetPointError(i, 0.0, total_error);

		double x1, value1 = 0.0;
		gr2_Pb->GetPoint(i, x1, value1);
		double value1_err = gr2_Pb->GetErrorY(i);
		
		double x2, value2 = 0.0;
		gr2_p->GetPoint(i, x2, value2);
		double value2_err = gr2_p->GetErrorY(i);

		double total_value = (value1 + value2)/2.0;
		double total_error = sqrt(value1_err*value1_err + value2_err*value2_err)/2.0;

		gr2_new->SetPoint(i, x1, total_value);
		gr2_new->SetPointError(i, 0.0, total_error);

		double x1, value1 = 0.0;
		gr3_Pb[0]->GetPoint(i, x1, value1);
		double value1_err = gr3_Pb[0]->GetErrorY(i);
		
		double x2, value2 = 0.0;
		gr3_Pb[1]->GetPoint(i, x2, value2);
		double value2_err = gr3_Pb[1]->GetErrorY(i);

		double x3, value3 = 0.0;
		gr3_p[0]->GetPoint(i, x3, value3);
		double value3_err = gr3_p[0]->GetErrorY(i);

		double x4, value4 = 0.0;
		gr3_p[1]->GetPoint(i, x4, value4);
		double value4_err = gr3_p[1]->GetErrorY(i);

		double total_value = (value1 + value2 + value3 + value4)/4.0;
		double total_error = sqrt(value1_err*value1_err + value2_err*value2_err + value3_err*value3_err + value4_err*value4_err)/4.0;

		gr3_new->SetPoint(i, x1, total_value);
		gr3_new->SetPointError(i, 0.0, total_error);

		double x1, value1 = 0.0;
		gr4_Pb[0]->GetPoint(i, x1, value1);
		double value1_err = gr4_Pb[0]->GetErrorY(i);
		
		double x2, value2 = 0.0;
		gr4_Pb[1]->GetPoint(i, x2, value2);
		double value2_err = gr4_Pb[1]->GetErrorY(i);

		double x3, value3 = 0.0;
		gr4_p[0]->GetPoint(i, x3, value3);
		double value3_err = gr4_p[0]->GetErrorY(i);

		double x4, value4 = 0.0;
		gr4_p[1]->GetPoint(i, x4, value4);
		double value4_err = gr4_p[1]->GetErrorY(i);

		double total_value = (value1 + value2 + value3 + value4)/4.0;
		double total_error = sqrt(value1_err*value1_err + value2_err*value2_err + value3_err*value3_err + value4_err*value4_err)/4.0;

		gr4_new->SetPoint(i, x1, total_value);
		gr4_new->SetPointError(i, 0.0, total_error);

		double x1, value1 = 0.0;
		gr5_Pb->GetPoint(i, x1, value1);
		double value1_err = gr5_Pb->GetErrorY(i);
		
		double x2, value2 = 0.0;
		gr5_p->GetPoint(i, x2, value2);
		double value2_err = gr5_p->GetErrorY(i);

		double total_value = (value1 + value2)/2.0;
		double total_error = sqrt(value1_err*value1_err + value2_err*value2_err)/2.0;

		gr5_new->SetPoint(i, x1, total_value);
		gr5_new->SetPointError(i, 0.0, total_error);

	}

	//lamba-lambda same sign
	gr1_new->SetMarkerStyle(20);
	gr1_new->SetMarkerSize(1.4);
	gr1_new->SetMarkerColor(kRed);
	gr1_new->SetLineColor(kRed);
	gr1_new->Draw("Psame");

	//lambda-antilambda opposite sign
	gr1_Pb[2]->SetMarkerStyle(21);
	gr1_Pb[2]->SetMarkerSize(1.4);
	gr1_Pb[2]->SetMarkerColor(kBlue);
	gr1_Pb[2]->SetLineColor(kBlue);
	gr1_Pb[2]->Draw("Psame");

	//K0s-K0s
	gr2_new->SetMarkerStyle(24);
	gr2_new->SetMarkerSize(1.4);
	gr2_new->SetMarkerColor(kRed);
	gr2_new->SetLineColor(kRed);
	gr2_new->Draw("Psame");

	//Lambda(antiLambda)-K0s
	gr3_new->SetMarkerStyle(25);
	gr3_new->SetMarkerSize(1.4);
	gr3_new->SetMarkerColor(kBlack);
	gr3_new->SetLineColor(kBlack);
	gr3_new->Draw("Psame");

	//Lambda(antiLambda)-H
	gr4_new->SetMarkerStyle(27);
	gr4_new->SetMarkerSize(1.4);
	gr4_new->SetMarkerColor(kGreen-1);
	gr4_new->SetLineColor(kGreen-1);
	gr4_new->Draw("Psame");

	//K0s-H
	gr5_new->SetMarkerStyle(28);
	gr5_new->SetMarkerSize(1.4);
	gr5_new->SetMarkerColor(kPink-1);
	gr5_new->SetLineColor(kPink-1);
	gr5_new->Draw("Psame");


	TLegend *w1 = new TLegend(0.47,0.72,0.82,0.89);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(17);
    w1->SetTextFont(43);
    w1->AddEntry(gr1_new, "#Lambda-#Lambda", "P");
    w1->AddEntry(gr1_Pb[2], "#Lambda-#bar{#Lambda}", "P");
    w1->AddEntry(gr2_new, "K^{0}_{s}-K^{0}_{s}", "P");
    w1->AddEntry(gr3_new, "#Lambda-K^{0}_{s}", "P");
    w1->AddEntry(gr4_new, "#Lambda-H", "P");
    w1->AddEntry(gr5_new, "K^{0}_{s}-H", "P");
    w1->Draw("same");


}