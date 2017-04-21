#include "RiceStyle.h"

using namespace std;


void plotV0s_pPb(){

	gStyle->SetErrorX(0);

	TFile* file[10];
	file[0] = new TFile("../dataPoints/V0s_pPb_test.root");

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

	TH1D* base1 = makeHist("base1", "Pb-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 500,0,500,kBlack);

	base1->GetYaxis()->SetRangeUser(-1.0, 1.0);
	base1->GetXaxis()->SetRangeUser(0,500);
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

	TGraphErrors* gr1_new_Pb = new TGraphErrors(6);
	TGraphErrors* gr1_new_p = new TGraphErrors(6);

	TGraphErrors* gr2_new_Pb = new TGraphErrors(6);
	TGraphErrors* gr2_new_p = new TGraphErrors(6);

	TGraphErrors* gr3_new_Pb = new TGraphErrors(6);
	TGraphErrors* gr3_new_p = new TGraphErrors(6);

	TGraphErrors* gr4_new_Pb = new TGraphErrors(6);
	TGraphErrors* gr4_new_p = new TGraphErrors(6);

	TGraphErrors* gr5_new_Pb = new TGraphErrors(6);
	TGraphErrors* gr5_new_p = new TGraphErrors(6);

	for(int i = 0; i < 10; i++){

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

		double total_value_Pb = (value1 + value2 )/2.0;
		double total_error_Pb = sqrt(value1_err*value1_err + value2_err*value2_err )/2.0;

		gr1_new_Pb->SetPoint(i, x1, total_value_Pb);
		gr1_new_Pb->SetPointError(i, 0.0, total_error_Pb);

		double total_value_p = (value3 + value4)/2.0;
		double total_error_p = sqrt(value3_err*value3_err + value4_err*value4_err)/2.0;

		gr1_new_p->SetPoint(i, x1, total_value_p);
		gr1_new_p->SetPointError(i, 0.0, total_error_p);


		double x1, value1 = 0.0;
		gr2_Pb->GetPoint(i, x1, value1);
		double value1_err = gr2_Pb->GetErrorY(i);
		
		double x2, value2 = 0.0;
		gr2_p->GetPoint(i, x2, value2);
		double value2_err = gr2_p->GetErrorY(i);
		
		gr2_new_Pb->SetPoint(i, x1, value1);
		gr2_new_Pb->SetPointError(i, 0.0, value1_err);

		gr2_new_p->SetPoint(i, x1, value2);
		gr2_new_p->SetPointError(i, 0.0, value2_err);


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

		double total_value_Pb = (value1 + value2 )/2.0;
		double total_error_Pb = sqrt(value1_err*value1_err + value2_err*value2_err)/2.0;

		gr3_new_Pb->SetPoint(i, x1, total_value_Pb);
		gr3_new_Pb->SetPointError(i, 0.0, total_error_Pb);

		double total_value_p = (value3 + value4)/2.0;
		double total_error_p = sqrt(value3_err*value3_err + value4_err*value4_err)/2.0;

		gr3_new_p->SetPoint(i, x1, total_value_p);
		gr3_new_p->SetPointError(i, 0.0, total_error_p);

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

		double total_value_Pb = (value1 + value2)/2.0;
		double total_error_Pb = sqrt(value1_err*value1_err + value2_err*value2_err)/2.0;

		gr4_new_Pb->SetPoint(i, x1, total_value_Pb);
		gr4_new_Pb->SetPointError(i, 0.0, total_error_Pb);

		double total_value_p = (value3 + value4)/2.0;
		double total_error_p = sqrt(value3_err*value3_err + value4_err*value4_err)/2.0;

		gr4_new_p->SetPoint(i, x1, total_value_p);
		gr4_new_p->SetPointError(i, 0.0, total_error_p);
		

		double x1, value1 = 0.0;
		gr5_Pb->GetPoint(i, x1, value1);
		double value1_err = gr5_Pb->GetErrorY(i);
		
		double x2, value2 = 0.0;
		gr5_p->GetPoint(i, x2, value2);
		double value2_err = gr5_p->GetErrorY(i);

		gr5_new_Pb->SetPoint(i, x1, value1);
		gr5_new_Pb->SetPointError(i, 0.0, value1_err);

		gr5_new_p->SetPoint(i, x1, value2);
		gr5_new_p->SetPointError(i, 0.0, value2_err);

	}

	//lamba-lambda same sign
	gr1_new_Pb->SetMarkerStyle(20);
	gr1_new_Pb->SetMarkerSize(1.4);
	gr1_new_Pb->SetMarkerColor(kRed);
	gr1_new_Pb->SetLineColor(kRed);
	gr1_new_Pb->Draw("Psame");

	//lambda-antilambda opposite sign
	gr1_Pb[2]->SetMarkerStyle(21);
	gr1_Pb[2]->SetMarkerSize(1.4);
	gr1_Pb[2]->SetMarkerColor(kBlue);
	gr1_Pb[2]->SetLineColor(kBlue);
	gr1_Pb[2]->Draw("Psame");


	gr1_new_p->SetMarkerStyle(24);
	gr1_new_p->SetMarkerSize(1.4);
	gr1_new_p->SetMarkerColor(kRed);
	gr1_new_p->SetLineColor(kRed);
	gr1_new_p->Draw("Psame");

	//lambda-antilambda opposite sign
	gr1_p[2]->SetMarkerStyle(25);
	gr1_p[2]->SetMarkerSize(1.4);
	gr1_p[2]->SetMarkerColor(kBlue);
	gr1_p[2]->SetLineColor(kBlue);
	gr1_p[2]->Draw("Psame");

	TLegend *w1 = new TLegend(0.47,0.72,0.82,0.89);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(17);
    w1->SetTextFont(43);
    w1->AddEntry(gr1_new_Pb, "#Lambda-#Lambda Pb-going", "P");
    w1->AddEntry(gr1_Pb[2], "#Lambda-#bar{#Lambda} Pb-going", "P");
    w1->AddEntry(gr1_new_p, "#Lambda-#Lambda p-going", "P");
    w1->AddEntry(gr1_p[2], "#Lambda-#bar{#Lambda} p-going", "P");

    // w1->AddEntry(gr2_new_Pb, "K^{0}_{s}-K^{0}_{s}", "P");
    // w1->AddEntry(gr3_new_Pb, "#Lambda-K^{0}_{s}", "P");
    // w1->AddEntry(gr4_new_Pb, "#Lambda-H", "P");
    // w1->AddEntry(gr5_new_Pb, "K^{0}_{s}-H", "P");
    w1->Draw("same");


	TH1D* base2 = makeHist("base2", "Pb-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 500,0,500,kBlack);

	base2->GetYaxis()->SetRangeUser(-1.0, 1.0);
	base2->GetXaxis()->SetRangeUser(0,500);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.23);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetNdivisions(8,8,0);
	base2->GetYaxis()->SetNdivisions(4,6,0);
	
	TCanvas* c2 = new TCanvas("c2","c2",1,1,650,650);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	base2->Draw();

	TGraphErrors* temp = new TGraphErrors(10);
	TGraphErrors* temp1 = new TGraphErrors(10);

	for(int i = 0; i < 10; i++){
		double xx = 0.0;
		double yy = 0.0;
		gr1_new_Pb->GetPoint(i, xx, yy);
		double ey = gr1_new_Pb->GetErrorY(i);
		
		double x = 0.0; 
		double y = 0.0;
		gr1_Pb[2]->GetPoint(i, x, y);
		double ey1 = gr1_Pb[2]->GetErrorY(i);

		temp->SetPoint(i, x, y-yy);
		temp->SetPointError(i, 0, sqrt(ey1*ey1 + ey*ey));

		double xx = 0.0;
		double yy = 0.0;
		gr1_new_p->GetPoint(i, xx, yy);
		double ey = gr1_new_p->GetErrorY(i);
		
		double x = 0.0; 
		double y = 0.0;
		gr1_p[2]->GetPoint(i, x, y);
		double ey1 = gr1_p[2]->GetErrorY(i);

		temp1->SetPoint(i, x, y-yy);
		temp1->SetPointError(i, 0, sqrt(ey1*ey1 + ey*ey));

	}

	temp->SetMarkerStyle(20);
	temp->SetMarkerSize(1.4);
	temp->SetMarkerColor(kRed);
	temp->SetLineColor(kRed);
	temp->Draw("Psame");


	temp1->SetMarkerStyle(20);
	temp1->SetMarkerSize(1.4);
	temp1->SetMarkerColor(kBlue);
	temp1->SetLineColor(kBlue);
	temp1->Draw("Psame");


	// TFile f1("../dataPoints/pPb_Ntrk.root","RECREATE");
	// temp->Write();
	// temp1->Write();


	TH1D* base3 = makeHist("base3", "Pb-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 500,0,500,kBlack);

	base3->GetYaxis()->SetRangeUser(-1.0, 1.0);
	base3->GetXaxis()->SetRangeUser(0,500);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.23);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.4);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.4);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.5);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.4);
	base3->GetXaxis()->SetNdivisions(8,8,0);
	base3->GetYaxis()->SetNdivisions(4,6,0);
	
	TCanvas* c3 = new TCanvas("c3","c3",1,1,650,650);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	base3->Draw();

	//K0s-K0s
	gr2_new_Pb->SetMarkerStyle(24);
	gr2_new_Pb->SetMarkerSize(1.4);
	gr2_new_Pb->SetMarkerColor(kRed);
	gr2_new_Pb->SetLineColor(kRed);
	gr2_new_Pb->Draw("Psame");

	//Lambda(antiLambda)-K0s
	gr3_new_Pb->SetMarkerStyle(25);
	gr3_new_Pb->SetMarkerSize(1.4);
	gr3_new_Pb->SetMarkerColor(kBlack);
	gr3_new_Pb->SetLineColor(kBlack);
	gr3_new_Pb->Draw("Psame");

	//Lambda(antiLambda)-H
	gr4_new_Pb->SetMarkerStyle(27);
	gr4_new_Pb->SetMarkerSize(1.4);
	gr4_new_Pb->SetMarkerColor(kGreen-1);
	gr4_new_Pb->SetLineColor(kGreen-1);
	gr4_new_Pb->Draw("Psame");

	//K0s-H
	gr5_new_Pb->SetMarkerStyle(28);
	gr5_new_Pb->SetMarkerSize(1.4);
	gr5_new_Pb->SetMarkerColor(kPink-1);
	gr5_new_Pb->SetLineColor(kPink-1);
	gr5_new_Pb->Draw("Psame");


	TLegend *w2 = new TLegend(0.47,0.72,0.82,0.89);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(17);
    w2->SetTextFont(43);
    w2->AddEntry(gr2_new_Pb, "K^{0}_{s}-K^{0}_{s}", "P");
    w2->AddEntry(gr3_new_Pb, "#Lambda-K^{0}_{s}", "P");
    w2->AddEntry(gr4_new_Pb, "#Lambda-H", "P");
    w2->AddEntry(gr5_new_Pb, "K^{0}_{s}-H", "P");
    w2->Draw("same");

	TFile f1("../dataPoints/pPb_Ntrk_K0s.root","RECREATE");
	gr2_new_Pb->Write();
	gr3_new_Pb->Write();
	gr4_new_Pb->Write();
	gr5_new_Pb->Write();

	TH1D* base4 = makeHist("base4", "p-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 500,0,500,kBlack);

	base4->GetYaxis()->SetRangeUser(-1.0, 1.0);
	base4->GetXaxis()->SetRangeUser(0,500);
	base4->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base4,1.1,1.25);

	base4->GetYaxis()->SetTitleOffset(1.23);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.4);
	base4->GetXaxis()->SetTitleSize(base4->GetXaxis()->GetTitleSize()*1.4);
	base4->GetYaxis()->SetLabelSize(base4->GetYaxis()->GetLabelSize()*1.5);
	base4->GetXaxis()->SetLabelSize(base4->GetXaxis()->GetLabelSize()*1.4);
	base4->GetXaxis()->SetNdivisions(8,8,0);
	base4->GetYaxis()->SetNdivisions(4,6,0);
	
	TCanvas* c4 = new TCanvas("c4","c4",1,1,650,650);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	base4->Draw();

	//K0s-K0s
	gr2_new_p->SetMarkerStyle(24);
	gr2_new_p->SetMarkerSize(1.4);
	gr2_new_p->SetMarkerColor(kRed);
	gr2_new_p->SetLineColor(kRed);
	gr2_new_p->Draw("Psame");

	//Lambda(antiLambda)-K0s
	gr3_new_p->SetMarkerStyle(25);
	gr3_new_p->SetMarkerSize(1.4);
	gr3_new_p->SetMarkerColor(kBlack);
	gr3_new_p->SetLineColor(kBlack);
	gr3_new_p->Draw("Psame");

	//Lambda(antiLambda)-H
	gr4_new_p->SetMarkerStyle(27);
	gr4_new_p->SetMarkerSize(1.4);
	gr4_new_p->SetMarkerColor(kGreen-1);
	gr4_new_p->SetLineColor(kGreen-1);
	gr4_new_p->Draw("Psame");

	//K0s-H
	gr5_new_p->SetMarkerStyle(28);
	gr5_new_p->SetMarkerSize(1.4);
	gr5_new_p->SetMarkerColor(kPink-1);
	gr5_new_p->SetLineColor(kPink-1);
	gr5_new_p->Draw("Psame");


	TLegend *w2 = new TLegend(0.47,0.72,0.82,0.89);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(17);
    w2->SetTextFont(43);
    w2->AddEntry(gr2_new_p, "K^{0}_{s}-K^{0}_{s}", "P");
    w2->AddEntry(gr3_new_p, "#Lambda-K^{0}_{s}", "P");
    w2->AddEntry(gr4_new_p, "#Lambda-H", "P");
    w2->AddEntry(gr5_new_p, "K^{0}_{s}-H", "P");
    w2->Draw("same");


}