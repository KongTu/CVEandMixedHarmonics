#include "RiceStyle.h"

using namespace std;

double PbPb_centralityBinCenter[] = {75,67.5,62.5,57.5,52.5,47.5,42.5,37.5,32.5};
double PbPb_centralityBinCenter_tracker[] = {65, 55, 45, 35};

double total_systematics_pPb = 0.000041;
double total_systematics_PbPb = 0.000028;


void plot3pCMEIntegratedResults_diff(){

	gStyle->SetErrorX(0);

	TFile* file = new TFile("../dataPoints/CME_Ntrk_112_pPb_8TeV.root");
	TFile* file1 = new TFile("../dataPoints/CME_Ntrk_123_pPb_8TeV.root");
	TFile* file_PbPb_1 = new TFile("../dataPoints/CME_Ntrk_112_PbPb_5TeV.root");
	TFile* file_PbPb_2 = new TFile("../dataPoints/CME_Ntrk_123_PbPb_5TeV.root");
	
	TGraphErrors* gr112[10];
	TGraphErrors* gr123[10];

	TGraphErrors* gr112_PbPb[10];
	TGraphErrors* gr123_PbPb[10];

//all files are as function of Ntrk
/*8 TGraph in each file
	- graph1, 3p Pb-going SS
	- graph2, 3p Pb-going OS
	- graph3, 3p p-going SS
	- graph4, 3p p-going OS
	- graph5, vn (n = 2 for 112, n=3 for 123)
	- graph6, dummy, same as graph5
	- graph7, 2p SS
	- graph8, 2p OS

	in PbPb, graph3-4 is not in use, only 1 and 2 are enough. 
*/
	for(int i = 0 ; i < 8; i++){

		gr112[i] = (TGraphErrors*) file->Get( Form("Graph;%d",i+1) );
		gr123[i] = (TGraphErrors*) file1->Get( Form("Graph;%d",i+1) );

		gr112_PbPb[i] = (TGraphErrors*) file_PbPb_1->Get( Form("Graph;%d",i+1) );
		gr123_PbPb[i] = (TGraphErrors*) file_PbPb_2->Get( Form("Graph;%d",i+1) );

	}

	TFile* file2[10];
	file2[0] = new TFile("~/Dropbox/CMEplottingMacros/dataPoints/pPb_data.root");
	file2[1] = new TFile("~/Dropbox/CMEplottingMacros/dataPoints/PbPb5TeV_data.root");
	file2[2] = new TFile("~/Dropbox/CMEplottingMacros/dataPoints/PbPb_5TeV_centrality_data.root");


	TGraphErrors* gr1[4];
	for(int i = 0; i < 4; i++){

		gr1[i] = (TGraphErrors*) file2[0]->Get(Form("Graph;%d", i+1));
	}

	TGraphErrors* gr2[3];
	for(int i = 0; i < 3; i++){

		gr2[i] = (TGraphErrors*) file2[1]->Get(Form("Graph;%d", i+1));
	}
	TGraphErrors* gr3[2];
	for(int i = 0; i < 2; i++){

		gr3[i] = (TGraphErrors*) file2[2]->Get(Form("Graph;%d", i+1));
	}

//fine tune the binning here
//pPb Ntrk
	for(int i = 0; i < 4; i++){
		for(int mult = 0; mult < 3; mult++){
			double x, y;
			gr1[i]->GetPoint(mult, x, y);
			gr1[i]->SetPoint(mult, x, 100);
		}

		gr1[i]->SetPoint(9,x,100);
	}

/*
plotting difference Fig.3
*/

		TGaxis::SetMaxDigits(3);

	
	TH1D* base5 = makeHist("base5", "Pb-going", "N^{offline}_{trk}", "#gamma_{112}(OS#minusSS)", 5000,0.1,10000,kBlack);

	base5->GetYaxis()->SetRangeUser(-0.00006, 0.0014);
	base5->GetXaxis()->SetRangeUser(72, 1700);

	base5->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base5,0.9,1.25);

	base5->GetYaxis()->SetTitleOffset(1.3);
	base5->GetYaxis()->SetTitleSize(base5->GetYaxis()->GetTitleSize()*1.4);
	base5->GetXaxis()->SetTitleSize(base5->GetXaxis()->GetTitleSize()*1.4);
	base5->GetYaxis()->SetLabelSize(base5->GetYaxis()->GetLabelSize()*1.4);
	base5->GetXaxis()->SetLabelSize(base5->GetXaxis()->GetLabelSize()*1.4);
	base5->GetXaxis()->SetNdivisions(8,18,0);
	base5->GetYaxis()->SetNdivisions(4,6,0);
	
	TBox *box7[50];
    TBox *box8[50];
    TBox *box9[50];
    TBox *box10[50];

    TCanvas* c3 = new TCanvas("c3","c3",1,1,1200,450);
	c3->Divide(3,1,0.01,0.01);
	c3->cd(1);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.2);
	gPad->SetTopMargin(0.07);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.0);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	base5->GetXaxis()->SetLabelOffset(999);
	base5->GetXaxis()->SetTickLength(0);

	TGaxis *newaxis2 = new TGaxis(81,
	                            -0.00006,
	                            1700,
	                            -0.00006,
	                            81,
	                            1700,
	                            510,"G");
	newaxis2->SetLabelOffset(0.01);
	newaxis2->SetLabelFont(42);
	//newaxis2->SetLabelSize(0.045);

	base5->Draw();
	newaxis2->Draw("Asame");

	TGraphErrors* new_gr1_pPb_8TeV = new TGraphErrors(6);
	TGraphErrors* new_gr2_pPb_8TeV = new TGraphErrors(6);
	
	TGraphErrors* new_gr1_pPb = new TGraphErrors(9);
	TGraphErrors* new_gr2_pPb = new TGraphErrors(9);
	
	TGraphErrors* new_gr1_PbPb = new TGraphErrors(9);
	TGraphErrors* new_gr2_PbPb = new TGraphErrors(9);

	double xe1[6];
    double xe2[9]; 
    double xe3[4];

	for(int mult = 0; mult < 6; mult++){

		xe1[mult] = 7*log(1.1*(mult+1));
    	if(mult == 0) xe1[mult] = 6;
    	double ye = total_systematics_pPb;

//pPb Pb-going
		double x1, y1, y1_error;
		gr1[0]->GetPoint(mult+3, x1, y1);
		y1_error = gr1[0]->GetErrorY(mult+3);

		double x2, y2, y2_error;
		gr1[1]->GetPoint(mult+3, x2, y2);
		y2_error = gr1[1]->GetErrorY(mult+3);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr1_pPb->SetPoint(mult+3, x1, y2-y1);
		new_gr1_pPb->SetPointError(mult+3, 0, total_error);

		box7[mult] = new TBox(x1-xe1[mult],y2-y1-ye,x1+xe1[mult],y2-y1+ye);
		box7[mult]->SetFillColor(kRed);
	 	box7[mult]->SetFillColorAlpha(kGray+2,0.4);
        box7[mult]->SetFillStyle(1001);
    	box7[mult]->SetLineWidth(0);
    	box7[mult]->SetLineColor(kRed);
        box7[mult]->Draw("SAME");

//pPb p-going
		double x1, y1, y1_error;
		gr1[2]->GetPoint(mult+3, x1, y1);
		y1_error = gr1[2]->GetErrorY(mult+3);

		double x2, y2, y2_error;
		gr1[3]->GetPoint(mult+3, x2, y2);
		y2_error = gr1[3]->GetErrorY(mult+3);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		if(mult == 5) y2 = 100;
		new_gr2_pPb->SetPoint(mult+3, x1, y2-y1);
		new_gr2_pPb->SetPointError(mult+3, 0, total_error);

		box8[mult] = new TBox(x1-xe1[mult],y2-y1-ye,x1+xe1[mult],y2-y1+ye);
		box8[mult]->SetFillColor(kRed);
	 	box8[mult]->SetFillColorAlpha(kGray+2,0.4);
        box8[mult]->SetFillStyle(1001);
    	box8[mult]->SetLineWidth(0);
    	box8[mult]->SetLineColor(kRed);
        if(mult != 5) box8[mult]->Draw("SAME");


//pPb 8 TeV Pb-going
		double x1, y1, y1_error;
		double v2;
		gr112[4]->GetPoint(mult, x1, v2);

		gr112[0]->GetPoint(mult, x1, y1);
		y1_error = gr112[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr112[1]->GetPoint(mult, x2, y2);
		y2_error = gr112[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);
		
		new_gr1_pPb_8TeV->SetPoint(mult, x1, (y2-y1));
		new_gr1_pPb_8TeV->SetPointError(mult, 0, total_error);

//pPb 8 TeV p-going
		double x1, y1, y1_error;
		gr112[2]->GetPoint(mult, x1, y1);
		y1_error = gr112[2]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr112[3]->GetPoint(mult, x2, y2);
		y2_error = gr112[3]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		if(mult == 5) {
			y2 = 100;}
		new_gr2_pPb_8TeV->SetPoint(mult, x1, (y2-y1));
		new_gr2_pPb_8TeV->SetPointError(mult, 0, total_error);

	}


	for(int mult = 0; mult < 9; mult++){

		xe2[mult] = 7*log(1.1*(mult+1));
    	if(mult == 0) xe2[mult] = 6;

    	double ye = total_systematics_PbPb;

//PbPb Ntrk 
		double x1, y1, y1_error;
		gr2[0]->GetPoint(mult, x1, y1);
		y1_error = gr2[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr2[1]->GetPoint(mult, x2, y2);
		y2_error = gr2[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr1_PbPb->SetPoint(mult, x1, y2-y1);
		new_gr1_PbPb->SetPointError(mult, 0, total_error);

		box9[mult] = new TBox(x1-xe2[mult],y2-y1-ye,x1+xe2[mult],y2-y1+ye);
		box9[mult]->SetFillColor(kRed);
	 	box9[mult]->SetFillColorAlpha(kGray+2,0.4);
        box9[mult]->SetFillStyle(1001);
    	box9[mult]->SetLineWidth(0);
    	box9[mult]->SetLineColor(kRed);
        box9[mult]->Draw("SAME");

	}

	for(int mult = 0; mult < 4; mult++){

		xe3[mult] = 30*log(1.9*(mult+1));
    	if(mult == 0) xe3[mult] = 35;
    	double ye = total_systematics_PbPb;
	
//PbPb Centrality		
		double x1, y1, y1_error;
		gr3[0]->GetPoint(mult, x1, y1);
		y1_error = gr3[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr3[1]->GetPoint(mult, x2, y2);
		y2_error = gr3[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr2_PbPb->SetPoint(mult, x1, y2-y1);
		new_gr2_PbPb->SetPointError(mult, 0, total_error);

		box10[mult] = new TBox(x1-xe3[mult],y2-y1-ye,x1+xe3[mult],y2-y1+ye);
		box10[mult]->SetFillColor(kRed);
	 	box10[mult]->SetFillColorAlpha(kGray+2,0.4);
        box10[mult]->SetFillStyle(1001);
    	box10[mult]->SetLineWidth(0);
    	box10[mult]->SetLineColor(kRed);
        box10[mult]->Draw("SAME");

	}

	new_gr1_pPb_8TeV->SetMarkerStyle(29);
	new_gr1_pPb_8TeV->SetMarkerSize(1.8);
	new_gr1_pPb_8TeV->SetMarkerColor(kBlack);
	new_gr1_pPb_8TeV->SetLineColor(kBlack);
	new_gr1_pPb_8TeV->Draw("Psame");

	new_gr2_pPb_8TeV->SetMarkerStyle(23);
	new_gr2_pPb_8TeV->SetMarkerSize(1.8);
	new_gr2_pPb_8TeV->SetMarkerColor(kGreen-2);
	new_gr2_pPb_8TeV->SetLineColor(kGreen-2);
	new_gr2_pPb_8TeV->Draw("Psame");

	new_gr1_pPb->SetMarkerStyle(20);
	new_gr1_pPb->SetMarkerSize(1.8);
	new_gr1_pPb->SetMarkerColor(kRed);
	new_gr1_pPb->SetLineColor(kRed);
	//new_gr1_pPb->Draw("Psame");

	new_gr2_pPb->SetMarkerStyle(21);
	new_gr2_pPb->SetMarkerSize(1.8);
	new_gr2_pPb->SetMarkerColor(kBlue);
	new_gr2_pPb->SetLineColor(kBlue);
	//new_gr2_pPb->Draw("Psame");

	new_gr1_PbPb->SetMarkerStyle(28);
	new_gr1_PbPb->SetMarkerSize(2.0);
	new_gr1_PbPb->SetMarkerColor(kBlack);
	new_gr1_PbPb->SetLineColor(kBlack);
	new_gr1_PbPb->Draw("Psame");

	new_gr2_PbPb->SetMarkerStyle(28);
	new_gr2_PbPb->SetMarkerSize(2.0);
	new_gr2_PbPb->SetMarkerColor(kBlack);
	new_gr2_PbPb->SetLineColor(kBlack);
	new_gr2_PbPb->Draw("Psame");

 	TLegend *ww3 = new TLegend(0.43,0.65,0.6,0.9);
    ww3->SetLineColor(kWhite);
    ww3->SetFillColor(0);
    ww3->SetTextSize(17);
    ww3->SetTextFont(43);
    ww3->AddEntry(new_gr1_pPb, "pPb 5.02 TeV, #phi_{c}(Pb-going)", "P" );
    ww3->AddEntry(new_gr2_pPb, "pPb 5.02 TeV, #phi_{c}(p-going)", "P" );
    ww3->AddEntry(new_gr1_pPb_8TeV, "pPb 8.16 TeV, #phi_{c}(Pb-going)", "P" );
    ww3->AddEntry(new_gr2_pPb_8TeV, "pPb 8.16 TeV, #phi_{c}(p-going)", "P" );
    ww3->AddEntry(new_gr1_PbPb, "PbPb", "P" );
    ww3->Draw("same");

   	TLatex* r11 = new TLatex(0.63,0.94, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.05);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.73,0.94, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(20);
    r22->SetTextFont(53);
    r22->Draw("same");

    c3->cd(2);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.2);
	gPad->SetTopMargin(0.07);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.0);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	TH1D* base6 = (TH1D*) base5->Clone("base6");
	base6->GetYaxis()->SetTitle("#gamma_{123}(OS#minusSS)");

	TGaxis *newaxis3 = new TGaxis(81,
	                            -0.00006,
	                            1700,
	                            -0.00006,
	                            81,
	                            1700,
	                            510,"G");
	newaxis3->SetLabelOffset(0.01);
	newaxis3->SetLabelFont(42);
	//newaxis3->SetLabelSize(0.045);

	base6->Draw();
	newaxis3->Draw("Asame");

	TGraphErrors* new_gr1_123_pPb_8TeV = new TGraphErrors(6);
	TGraphErrors* new_gr2_123_pPb_8TeV = new TGraphErrors(6);

	TGraphErrors* new_gr1_123_PbPb = new TGraphErrors(6);


	for(int mult = 0; mult < 6; mult++){
		
		//pPb 8 TeV Pb-going
		double x1, y1, y1_error;
		double v3;
		gr123[4]->GetPoint(mult, x1, v3);
		gr123[0]->GetPoint(mult, x1, y1);
		y1_error = gr123[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr123[1]->GetPoint(mult, x2, y2);
		y2_error = gr123[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);
		
		new_gr1_123_pPb_8TeV->SetPoint(mult, x1, (y2-y1));
		new_gr1_123_pPb_8TeV->SetPointError(mult, 0, total_error);

		//pPb 8 TeV p-going
		double x1, y1, y1_error;
		gr123[2]->GetPoint(mult, x1, y1);
		y1_error = gr123[2]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr123[3]->GetPoint(mult, x2, y2);
		y2_error = gr123[3]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		// if(mult == 5) {
		// 	y2 = 100;}
		new_gr2_123_pPb_8TeV->SetPoint(mult, x1, (y2-y1));
		new_gr2_123_pPb_8TeV->SetPointError(mult, 0, total_error);



		//PbPb Ntrk 
		double x1, y1, y1_error;
		double v3;
		gr123_PbPb[4]->GetPoint(mult, x1, v3);

		gr123_PbPb[0]->GetPoint(mult, x1, y1);
		y1_error = gr123_PbPb[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr123_PbPb[1]->GetPoint(mult, x2, y2);
		y2_error = gr123_PbPb[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr1_123_PbPb->SetPoint(mult, x1, (y2-y1));
		new_gr1_123_PbPb->SetPointError(mult, 0, total_error);


		}

		new_gr1_123_pPb_8TeV->SetMarkerStyle(20);
		new_gr1_123_pPb_8TeV->SetMarkerSize(1.8);
		new_gr1_123_pPb_8TeV->SetMarkerColor(kRed);
		new_gr1_123_pPb_8TeV->SetLineColor(kRed);
		new_gr1_123_pPb_8TeV->Draw("Psame");

		new_gr2_123_pPb_8TeV->SetMarkerStyle(21);
		new_gr2_123_pPb_8TeV->SetMarkerSize(1.8);
		new_gr2_123_pPb_8TeV->SetMarkerColor(kBlue);
		new_gr2_123_pPb_8TeV->SetLineColor(kBlue);
		new_gr2_123_pPb_8TeV->Draw("Psame");

		new_gr1_123_PbPb->SetMarkerStyle(28);
		new_gr1_123_PbPb->SetMarkerSize(1.8);
		new_gr1_123_PbPb->SetMarkerColor(kBlack);
		new_gr1_123_PbPb->SetLineColor(kBlack);
		new_gr1_123_PbPb->Draw("Psame");

	TLegend *ww4 = new TLegend(0.43,0.65,0.6,0.9);
    ww4->SetLineColor(kWhite);
    ww4->SetFillColor(0);
    ww4->SetTextSize(17);
    ww4->SetTextFont(43);
    ww4->AddEntry(new_gr1_123_pPb_8TeV, "pPb 8.16 TeV, #phi_{c}(Pb-going)", "P" );
    ww4->AddEntry(new_gr2_123_pPb_8TeV, "pPb 8.16 TeV, #phi_{c}(p-going)", "P" );
    ww4->AddEntry(new_gr1_123_PbPb, "PbPb", "P" );
    ww4->Draw("same");

	c3->cd(3);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.2);
	gPad->SetTopMargin(0.07);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.0);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	TH1D* base7 = (TH1D*) base5->Clone("base7");
	base7->GetYaxis()->SetRangeUser(0,0.01);
	base7->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}#minus#phi_{#beta})#GT (OS#minusSS)");

	TGaxis *newaxis4 = new TGaxis(81,
	                            0,
	                            1700,
	                            0,
	                            81,
	                            1700,
	                            510,"G");
	newaxis4->SetLabelOffset(0.01);
	newaxis4->SetLabelFont(42);
	//newaxis4->SetLabelSize(0.045);

	base7->Draw();
	newaxis4->Draw("Asame");

	TGraphErrors* new_delta_pPb_8TeV = new TGraphErrors(6);
	TGraphErrors* new_delta_PbPb = new TGraphErrors(6);


	for(int mult = 0; mult < 6; mult++){
		
		//pPb 8 TeV
		double x1, y1, y1_error;
		gr112[6]->GetPoint(mult, x1, y1);
		y1_error = gr112[6]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr112[7]->GetPoint(mult, x2, y2);
		y2_error = gr112[7]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);
		
		new_delta_pPb_8TeV->SetPoint(mult, x1, y2-y1);
		new_delta_pPb_8TeV->SetPointError(mult, 0, total_error);


		//PbPb Ntrk 
		double x1, y1, y1_error;
		gr112_PbPb[6]->GetPoint(mult, x1, y1);
		y1_error = gr112_PbPb[6]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr112_PbPb[7]->GetPoint(mult, x2, y2);
		y2_error = gr112_PbPb[7]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_delta_PbPb->SetPoint(mult, x1, y2-y1);
		new_delta_PbPb->SetPointError(mult, 0, total_error);


		}

		new_delta_pPb_8TeV->SetMarkerStyle(20);
		new_delta_pPb_8TeV->SetMarkerSize(1.8);
		new_delta_pPb_8TeV->SetMarkerColor(kRed);
		new_delta_pPb_8TeV->SetLineColor(kRed);
		new_delta_pPb_8TeV->Draw("Psame");

		new_delta_PbPb->SetMarkerStyle(28);
		new_delta_PbPb->SetMarkerSize(1.8);
		new_delta_PbPb->SetMarkerColor(kBlack);
		new_delta_PbPb->SetLineColor(kBlack);
		new_delta_PbPb->Draw("Psame");

	TLegend *ww5 = new TLegend(0.43,0.65,0.6,0.9);
    ww5->SetLineColor(kWhite);
    ww5->SetFillColor(0);
    ww5->SetTextSize(17);
    ww5->SetTextFont(43);
    ww5->AddEntry(new_delta_pPb_8TeV, "pPb 8.16 TeV", "P" );
    ww5->AddEntry(new_delta_PbPb, "PbPb", "P" );
    ww5->Draw("same");


}
