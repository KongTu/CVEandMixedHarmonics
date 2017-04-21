#include "RiceStyle.h"

using namespace std;

double PbPb_centralityBinCenter[] = {75,67.5,62.5,57.5,52.5,47.5,42.5,37.5,32.5};
double PbPb_centralityBinCenter_tracker[] = {65, 55, 45, 35};

double total_systematics_pPb = 0.000041;
double total_systematics_PbPb = 0.000028;


void plot3pCMEIntegratedResults(){

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


	TGaxis::SetMaxDigits(3);

    TLatex* r41 = new TLatex(0.48, 0.87, "pPb #sqrt{s_{NN}} = 8.16 TeV");
    r41->SetNDC();
    r41->SetTextSize(23);
    r41->SetTextFont(43);
    r41->SetTextColor(kBlack);

    TLatex* r42 = new TLatex(0.43, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.61,0.95, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.06);
    
    TLatex* r45 = new TLatex(0.5, 0.80, "185 #leq N^{offline}_{trk} < 250");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

    TLatex* r46 = new TLatex(0.5, 0.78, "185 #leq N^{offline}_{trk} < 250");
    r46->SetNDC();
    r46->SetTextSize(23);
    r46->SetTextFont(43);
    r46->SetTextColor(kBlack);

	TLatex* latex1 = new TLatex(0.51, 0.31, "  SS");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    TLatex* latex2 = new TLatex(0.62, 0.31, "OS");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);

/*
plotting Ntrk dependence, Fig.2
*/

	TH1D* base2 = makeHist("base2", "Pb-going", "N^{offline}_{trk}", "#gamma_{112}", 5000,0.1,10000,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.0009, 0.0006);
	base2->GetXaxis()->SetRangeUser(81, 1700);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,0.9,1.25);

	base2->GetYaxis()->SetTitleOffset(1.0);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.8);
	base2->GetXaxis()->SetNdivisions(8,18,0);
	base2->GetYaxis()->SetNdivisions(4,6,0);

	
	TCanvas* c2 = new TCanvas("c2","c2",1,1,1200,450);
	c2->Divide(3,1,0.01,0.01);
	c2->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.16);
	gPad->SetBottomMargin(0.16);
	gPad->SetRightMargin(0.0);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	base2->GetXaxis()->SetLabelOffset(999);
	base2->GetXaxis()->SetTickLength(0);
	
	TGaxis *newaxis2 = new TGaxis(81,
	                            -0.0009,
	                            1700,
	                            -0.0009,
	                            81,
	                            1700,
	                            510,"G");
	newaxis2->SetLabelOffset(0.01);
	newaxis2->SetLabelFont(42);
	newaxis2->SetLabelSize(newaxis2->GetLabelSize()*1.5);
	
	base2->Draw();
	newaxis2->Draw("SS");
	gPad->Update();
	gPad->SetLogx(1);


//systematic box:

    TBox *box1[50];
    TBox *box2[50];
    TBox *box3[50];
    TBox *box4[50];
    TBox *box5[50];
    TBox *box6[50];

    double xe[9];

    for(int mult = 0; mult < 6; mult++){

    	xe[mult] = 7*log(1.1*(mult+1));
    	if(mult == 0) xe[mult] = 6;
    	double ye = total_systematics_pPb;

    	double x1;
    	double value1;
    	gr1[0]->GetPoint(mult+3, x1, value1);

    	double x2;
    	double value2;
    	gr1[1]->GetPoint(mult+3, x2, value2);

    	box1[mult] = new TBox(x1-xe[mult],value1-ye,x1+xe[mult],value1+ye);
		box1[mult]->SetFillColor(kRed);
	 	box1[mult]->SetFillColorAlpha(kGray+2,0.4);
        box1[mult]->SetFillStyle(1001);
    	box1[mult]->SetLineWidth(0);
    	box1[mult]->SetLineColor(kRed);
        box1[mult]->Draw("SAME");

		box2[mult] = new TBox(x2-xe[mult],value2-ye,x2+xe[mult],value2+ye);
		box2[mult]->SetFillColor(kBlue);
	 	box2[mult]->SetFillColorAlpha(kGray+2,0.4);
        box2[mult]->SetFillStyle(1001);
    	box2[mult]->SetLineWidth(0);
    	box2[mult]->SetLineColor(kBlue);
        box2[mult]->Draw("SAME");

    }

    double xe2[11]; 

    for(int mult = 0; mult < 9; mult++){

    	xe2[mult] = 7*log(1.1*(mult+1));
    	if(mult == 0) xe2[mult] = 6;

    	double ye = total_systematics_PbPb;

    	double x1;
    	double value1;
    	gr2[0]->GetPoint(mult, x1, value1);

    	double x2;
    	double value2;
    	gr2[1]->GetPoint(mult, x2, value2);

    	box3[mult] = new TBox(x1-xe2[mult],value1[mult]-ye,x1+xe2[mult],value1[mult]+ye);
		box3[mult]->SetFillColor(kRed);
        box3[mult]->SetFillColorAlpha(kGray+2,0.4);
        box3[mult]->SetFillStyle(1001);
    	box3[mult]->SetLineWidth(0);
    	box3[mult]->SetLineColor(kRed);
        box3[mult]->Draw("SAME");

		box4[mult] = new TBox(x2-xe2[mult],value2[mult]-ye,x2+xe2[mult],value2[mult]+ye);
		box4[mult]->SetFillColor(kBlue);
        box4[mult]->SetFillColorAlpha(kGray+2,0.4);
        box4[mult]->SetFillStyle(1001);
    	box4[mult]->SetLineWidth(0);
    	box4[mult]->SetLineColor(kBlue);
        box4[mult]->Draw("SAME");
    }

    double xe3[4];
    for(int mult = 0; mult < 4; mult++){

    	xe3[mult] = 30*log(1.9*(mult+1));
    	if(mult == 0) xe3[mult] = 35;
    	double ye = total_systematics_PbPb;

    	double x1;
    	double value1;
    	gr3[0]->GetPoint(mult, x1, value1);

    	double x2;
    	double value2;
    	gr3[1]->GetPoint(mult, x2, value2);


    	box5[mult] = new TBox(x1-xe3[mult],value1[mult]-ye,x1+xe3[mult],value1[mult]+ye);
		box5[mult]->SetFillColor(kRed);
        box5[mult]->SetFillColorAlpha(kGray+2,0.4);
        box5[mult]->SetFillStyle(1001);
    	box5[mult]->SetLineWidth(0);
    	box5[mult]->SetLineColor(kRed);
        box5[mult]->Draw("SAME");

		box6[mult] = new TBox(x1-xe3[mult],value2[mult]-ye,x1+xe3[mult],value2[mult]+ye);
		box6[mult]->SetFillColor(kBlue);
        box6[mult]->SetFillColorAlpha(kGray+2,0.4);
        box6[mult]->SetFillStyle(1001);
    	box6[mult]->SetLineWidth(0);
    	box6[mult]->SetLineColor(kBlue);
        box6[mult]->Draw("SAME");
    }

//end of systematic box

	gr112[0]->SetMarkerStyle(20);
	gr112[0]->SetMarkerSize(1.6);
	gr112[0]->SetMarkerColor(kRed);
	gr112[0]->SetLineColor(kRed);
	gr112[0]->Draw("Psame");

	gr112[1]->SetMarkerStyle(21);
	gr112[1]->SetMarkerSize(1.6);
	gr112[1]->SetMarkerColor(kBlue);
	gr112[1]->SetLineColor(kBlue);
	gr112[1]->Draw("Psame");

	gr1[0]->SetMarkerStyle(24);
	gr1[0]->SetMarkerSize(1.4);
	gr1[0]->SetMarkerColor(kBlack);
	gr1[0]->SetLineColor(kBlack);
	gr1[0]->Draw("Psame");

	gr1[1]->SetMarkerStyle(25);
	gr1[1]->SetMarkerSize(1.4);
	gr1[1]->SetMarkerColor(kGreen-2);
	gr1[1]->SetLineColor(kGreen-2);
	gr1[1]->Draw("Psame");

	// gr1[2]->SetMarkerStyle(20);
	// gr1[2]->SetMarkerColor(kRed);
	// gr1[2]->SetLineColor(kRed);
	// gr1[2]->Draw("Psame");

	// gr1[3]->SetMarkerStyle(21);
	// gr1[3]->SetMarkerColor(kBlue);
	// gr1[3]->SetLineColor(kBlue);
	// gr1[3]->Draw("Psame");

	gr2[0]->SetMarkerStyle(24);
	gr2[0]->SetMarkerSize(1.4);
	gr2[0]->SetMarkerColor(kRed);
	gr2[0]->SetLineColor(kRed);
	gr2[0]->Draw("Psame");

	gr2[1]->SetMarkerStyle(25);
	gr2[1]->SetMarkerSize(1.4);
	gr2[1]->SetMarkerColor(kBlue);
	gr2[1]->SetLineColor(kBlue);
	gr2[1]->Draw("Psame");	

	gr3[0]->SetMarkerStyle(24);
	gr3[0]->SetMarkerSize(1.4);
	gr3[0]->SetMarkerColor(kRed);
	gr3[0]->SetLineColor(kRed);
	gr3[0]->Draw("Psame");

	gr3[1]->SetMarkerStyle(25);
	gr3[1]->SetMarkerSize(1.4);
	gr3[1]->SetMarkerColor(kBlue);
	gr3[1]->SetLineColor(kBlue);
	gr3[1]->Draw("Psame");

	TLegend *w2 = new TLegend(0.48,0.2,0.88,0.3);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(gr1[1], "SS, pPb 5.02 TeV", "P");
    w2->AddEntry(gr1[0], "OS, pPb 5.02 TeV", "P");

    w2->Draw("same");

    TLatex* r3 = new TLatex(0.63, 0.96, "PbPb centrality(%)");
    r3->SetNDC();
    r3->SetTextSize(18);
    r3->SetTextFont(43);
    r3->SetTextColor(kBlack);
    r3->Draw("same");

    TLatex* cent1[7];
    cent1[0] = new TLatex(0.58, 0.91, "55");
    cent1[1] = new TLatex(0.73, 0.91, "45");
    cent1[2] = new TLatex(0.86, 0.91, "35");
    cent1[3] = new TLatex(0.38, 0.91, "65");
    cent1[4] = new TLatex(0.75, 0.91, "15");
    cent1[5] = new TLatex(0.80, 0.91, "7.5");
    cent1[6] = new TLatex(0.84, 0.91, "2.5");

    for(int i = 0; i < 4; i++){
    	cent1[i]->SetNDC();
    	cent1[i]->SetTextSize(20);
	    cent1[i]->SetTextFont(43);
	    cent1[i]->SetTextColor(kBlack);
	    cent1[i]->Draw("same");
    }

    TLine* l1[7];
    l1[0] = new TLine(404.1,0.00057, 404.1, 0.0006);
    l1[0]->SetLineWidth(2);
    l1[0]->Draw("Lsame");

    l1[1] = new TLine(717.6,0.00057, 717.6, 0.0006);
    l1[1]->SetLineWidth(2);
    l1[1]->Draw("Lsame");

    l1[2] = new TLine(1141,0.00057, 1141, 0.0006);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    l1[3] = new TLine(81.412,0.00057, 81.412, 0.0006);
    l1[3]->SetLineWidth(2);
    //l1[3]->Draw("Lsame");

    l1[4] = new TLine(197,0.00057, 197, 0.0006);
    l1[4]->SetLineWidth(2);
    l1[4]->Draw("Lsame");

    l1[5] = new TLine(3577.6,0.00057, 3577.6, 0.0006);
    l1[5]->SetLineWidth(2);
    //l1[5]->Draw("Lsame");

    l1[6] = new TLine(4474,0.00057, 4474, 0.0006);
    l1[6]->SetLineWidth(2);
    //l1[6]->Draw("Lsame");

   	TLatex* r11 = new TLatex(0.16,0.96, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.05);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.28,0.96, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(18);
    r22->SetTextFont(53);
	r22->Draw("same");

	c2->cd(2);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.16);
	gPad->SetBottomMargin(0.16);
	gPad->SetRightMargin(0.0);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	TH1D* base3 = (TH1D*) base2->Clone("base3");
	base3->GetYaxis()->SetTitle("#gamma_{123}");
	base3->GetYaxis()->SetRangeUser(-0.002, 0.0007);
	base3->GetXaxis()->SetLabelOffset(999);
	base3->GetXaxis()->SetTickLength(0);
	
	TGaxis *newaxis3 = new TGaxis(81,
	                            -0.002,
	                            1700,
	                            -0.002,
	                            81,
	                            1700,
	                            510,"G");
	newaxis3->SetLabelOffset(0.01);
	newaxis3->SetLabelFont(42);
	newaxis3->SetLabelSize(newaxis3->GetLabelSize()*1.5);

	base3->Draw();
	newaxis3->Draw("SS");
	gPad->Update();
	gPad->SetLogx(1);

	gr123[0]->SetMarkerStyle(20);
	gr123[0]->SetMarkerSize(1.6);
	gr123[0]->SetMarkerColor(kRed);
	gr123[0]->SetLineColor(kRed);
	gr123[0]->Draw("Psame");

	gr123[1]->SetMarkerStyle(21);
	gr123[1]->SetMarkerSize(1.6);
	gr123[1]->SetMarkerColor(kBlue);
	gr123[1]->SetLineColor(kBlue);
	gr123[1]->Draw("Psame");

	gr123_PbPb[0]->SetMarkerStyle(24);
	gr123_PbPb[0]->SetMarkerSize(1.6);
	gr123_PbPb[0]->SetMarkerColor(kRed);
	gr123_PbPb[0]->SetLineColor(kRed);
	gr123_PbPb[0]->Draw("Psame");

	gr123_PbPb[1]->SetMarkerStyle(25);
	gr123_PbPb[1]->SetMarkerSize(1.6);
	gr123_PbPb[1]->SetMarkerColor(kBlue);
	gr123_PbPb[1]->SetLineColor(kBlue);
	gr123_PbPb[1]->Draw("Psame");

	r3->Draw("same");
	for(int i = 0; i < 4; i++){
    	cent1[i]->SetNDC();
    	cent1[i]->SetTextSize(20);
	    cent1[i]->SetTextFont(43);
	    cent1[i]->SetTextColor(kBlack);
	    cent1[i]->Draw("same");
    }

    l1[0] = new TLine(404.1,0.00067, 404.1, 0.0007);
    l1[0]->SetLineWidth(2);
    l1[0]->Draw("Lsame");

    l1[1] = new TLine(717.6,0.00067, 717.6, 0.0007);
    l1[1]->SetLineWidth(2);
    l1[1]->Draw("Lsame");

    l1[2] = new TLine(1141,0.00067, 1141, 0.0007);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    l1[3] = new TLine(81.412,0.00067, 81.412, 0.0007);
    l1[3]->SetLineWidth(2);
    //l1[3]->Draw("Lsame");

    l1[4] = new TLine(197,0.00067, 197, 0.0007);
    l1[4]->SetLineWidth(2);
    l1[4]->Draw("Lsame");
	
	TLegend *w1 = new TLegend(0.3,0.71,0.75,0.83);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(20);
    w1->SetTextFont(43);
    w1->SetNColumns(2);
    w1->AddEntry(gr112[0], "  ", "P");
    w1->AddEntry(gr112[1], "  pPb 8 TeV, #phi_{c}(Pb-going)", "P");
    w1->AddEntry(gr2[0], "  ", "P");
    w1->AddEntry(gr2[1], "  PbPb 5 TeV", "P");
    w1->Draw("same");

    TLatex* r4 = new TLatex(0.19, 0.84, "#sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(26);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
    //r4->Draw("same");

	TLatex* latex1 = new TLatex(0.29, 0.84, "SS");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");
    TLatex* latex2 = new TLatex(0.38, 0.84, "OS");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    latex2->Draw("same");


    c2->cd(3);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.16);
	gPad->SetBottomMargin(0.16);
	gPad->SetRightMargin(0.0);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	TH1D* base4 = (TH1D*) base2->Clone("base4");
	base4->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}#minus#phi_{#beta})#GT");
	base4->GetYaxis()->SetRangeUser(-0.0011, 0.01);
	base4->GetXaxis()->SetLabelOffset(999);
	base4->GetXaxis()->SetTickLength(0);
	
	TGaxis *newaxis4 = new TGaxis(81,
	                            -0.0011,
	                            1700,
	                            -0.0011,
	                            81,
	                            1700,
	                            510,"G");
	newaxis4->SetLabelOffset(0.01);
	newaxis4->SetLabelFont(42);
	newaxis4->SetLabelSize(newaxis4->GetLabelSize()*1.5);

	base4->Draw();
	newaxis4->Draw("SS");
	gPad->Update();
	gPad->SetLogx(1);

	gr112[6]->SetMarkerStyle(20);
	gr112[6]->SetMarkerSize(1.6);
	gr112[6]->SetMarkerColor(kRed);
	gr112[6]->SetLineColor(kRed);
	gr112[6]->Draw("Psame");

	gr112[7]->SetMarkerStyle(21);
	gr112[7]->SetMarkerSize(1.6);
	gr112[7]->SetMarkerColor(kBlue);
	gr112[7]->SetLineColor(kBlue);
	gr112[7]->Draw("Psame");

	gr112_PbPb[6]->SetMarkerStyle(24);
	gr112_PbPb[6]->SetMarkerSize(1.6);
	gr112_PbPb[6]->SetMarkerColor(kRed);
	gr112_PbPb[6]->SetLineColor(kRed);
	gr112_PbPb[6]->Draw("Psame");

	gr112_PbPb[7]->SetMarkerStyle(25);
	gr112_PbPb[7]->SetMarkerSize(1.6);
	gr112_PbPb[7]->SetMarkerColor(kBlue);
	gr112_PbPb[7]->SetLineColor(kBlue);
	gr112_PbPb[7]->Draw("Psame");

	r3->Draw("same");
	for(int i = 0; i < 4; i++){
    	cent1[i]->SetNDC();
    	cent1[i]->SetTextSize(20);
	    cent1[i]->SetTextFont(43);
	    cent1[i]->SetTextColor(kBlack);
	    cent1[i]->Draw("same");
    }

	l1[0] = new TLine(404.1,0.0098, 404.1, 0.01);
    l1[0]->SetLineWidth(2);
    l1[0]->Draw("Lsame");

    l1[1] = new TLine(717.6,0.0098, 717.6, 0.01);
    l1[1]->SetLineWidth(2);
    l1[1]->Draw("Lsame");

    l1[2] = new TLine(1141,0.0098, 1141, 0.01);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    l1[3] = new TLine(81.412,0.0098, 81.412, 0.01);
    l1[3]->SetLineWidth(2);
    //l1[3]->Draw("Lsame");

    l1[4] = new TLine(197,0.0098, 197, 0.01);
    l1[4]->SetLineWidth(2);
    l1[4]->Draw("Lsame");





}