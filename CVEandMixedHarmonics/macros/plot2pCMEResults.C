#include "RiceStyle.h"

using namespace std;

double PbPb_centralityBinCenter[] = {75,67.5,62.5,57.5,52.5,47.5,42.5,37.5,32.5};
double PbPb_centralityBinCenter_tracker[] = {65, 55, 45, 35};

double total_systematics_pPb = 0.000041;
double total_systematics_PbPb = 0.000028;


void plot2pCMEResults(){

		gStyle->SetErrorX(0);

	TFile* file = new TFile("../dataPoints/CME_deltaEta_pPb_8TeV.root");
	TFile* file1[10];
	file1[0] = new TFile("~/cernbox/2015RUN2work/2015Analysis/CMEandCorrelation/ThreePointCorrelator/dataPoints/pPb_data_2p.root");
	file1[1] = new TFile("~/cernbox/2015RUN2work/2015Analysis/CMEandCorrelation/ThreePointCorrelator/dataPoints/PbPb5TeV_data_2p.root");

	TGraphErrors* gr1[8];
	for(int i = 0; i < 8; i++){

		gr1[i] = (TGraphErrors*) file1[0]->Get(Form("Graph;%d", i+1));
	}

	TGraphErrors* gr2[7];
	for(int i = 0; i < 7; i++){

		gr2[i] = (TGraphErrors*) file1[1]->Get(Form("Graph;%d", i+1));
	}

	TH1D* temp[4];
	TGraphErrors* gr[10];

	temp[0] = (TH1D*) file->Get("temp1");
	temp[1] = (TH1D*) file->Get("temp2");
	temp[2] = (TH1D*) file->Get("temp3");
	temp[3] = (TH1D*) file->Get("temp4");

	for(int i = 0 ; i < 8; i++){

		gr[i] = (TGraphErrors*) file->Get( Form("Graph;%d",i+1) );
	}


//start plotting

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}-#phi_{#beta})#GT", 5000,0.1,10000,kBlack);
	TH1D* base2 = makeHist("base2", "p-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 5000,0.1,10000,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0009, 0.011);
	base1->GetXaxis()->SetRangeUser(81, 800);

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
	base3->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c} (OS - SS)");
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
	
	TCanvas* c3 = new TCanvas("c3", "c3", 1,1,650,650);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	base1->GetXaxis()->SetLabelOffset(999);
	base1->GetXaxis()->SetTickLength(0);
	
	TGaxis *newaxis2 = new TGaxis(81,
	                            -0.0009,
	                            800,
	                            -0.0009,
	                            81,
	                            800,
	                            510,"G");
	newaxis2->SetLabelOffset(0.01);
	newaxis2->SetLabelFont(42);
	
	base1->Draw();
	newaxis2->Draw("SS");
	gPad->Update();
	gPad->SetLogx(1);

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
    	gr1[6]->GetPoint(mult+3, x1, value1);

    	double x2;
    	double value2;
    	gr1[7]->GetPoint(mult+3, x2, value2);

    	box1[mult] = new TBox(x1-xe[mult],value1-ye,x1+xe[mult],value1+ye);
		box1[mult]->SetFillColor(kRed);
	 	box1[mult]->SetFillColorAlpha(kGray+1,0.3);
        box1[mult]->SetFillStyle(1001);
    	box1[mult]->SetLineWidth(0);
    	box1[mult]->SetLineColor(kRed);
        //box1[mult]->Draw("SAME");

		box2[mult] = new TBox(x2-xe[mult],value2-ye,x2+xe[mult],value2+ye);
		box2[mult]->SetFillColor(kBlue);
	 	box2[mult]->SetFillColorAlpha(kGray+1,0.3);
        box2[mult]->SetFillStyle(1001);
    	box2[mult]->SetLineWidth(0);
    	box2[mult]->SetLineColor(kBlue);
        //box2[mult]->Draw("SAME");

    }

	gr[6]->SetMarkerStyle(20);
	gr[6]->SetMarkerSize(1.4);
	gr[6]->SetMarkerColor(kBlack);
	gr[6]->SetLineColor(kBlack);
	gr[6]->Draw("Lsame");

	gr[7]->SetMarkerStyle(21);
	gr[7]->SetMarkerSize(1.4);
	gr[7]->SetMarkerColor(kGreen-2);
	gr[7]->SetLineColor(kGreen-2);
	gr[7]->Draw("Lsame");

	gr1[6]->SetMarkerStyle(20);
	gr1[6]->SetMarkerSize(1.4);
	gr1[6]->SetMarkerColor(kRed);
	gr1[6]->SetLineColor(kRed);
	gr1[6]->Draw("Psame");

	gr1[7]->SetMarkerStyle(21);
	gr1[7]->SetMarkerSize(1.4);
	gr1[7]->SetMarkerColor(kBlue);
	gr1[7]->SetLineColor(kBlue);
	gr1[7]->Draw("Psame");

	gr2[5]->SetMarkerStyle(24);
	gr2[5]->SetMarkerSize(1.4);
	gr2[5]->SetMarkerColor(kRed);
	gr2[5]->SetLineColor(kRed);
	gr2[5]->Draw("Psame");

	gr2[6]->SetMarkerStyle(25);
	gr2[6]->SetMarkerSize(1.4);
	gr2[6]->SetMarkerColor(kBlue);
	gr2[6]->SetLineColor(kBlue);
	gr2[6]->Draw("Psame");


	TLegend *w1 = new TLegend(0.48,0.70,0.88,0.82);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(18);
    w1->SetTextFont(43);
    //w1->SetHeader("CMS");

    w1->SetNColumns(2);

    w1->AddEntry(gr1[6], "  ", "P");
    w1->AddEntry(gr1[7], "  pPb 5.02 TeV", "P");

    w1->AddEntry(gr2[5], "  ", "P");
    w1->AddEntry(gr2[6], "  PbPb 5.02 TeV", "P");
  
    w1->AddEntry(gr[6], "  ", "L");
    w1->AddEntry(gr[7], "  pPb 8.16 TeV", "L");


    w1->Draw("same");

    TLatex* latex3 = new TLatex(0.49, 0.84, "SS");
    latex3->SetNDC();
    latex3->SetTextSize(20);
    latex3->SetTextFont(43);
    latex3->SetTextColor(kBlack);
    latex3->Draw("same");
    TLatex* latex4 = new TLatex(0.56, 0.84, "OS");
    latex4->SetNDC();
    latex4->SetTextSize(20);
    latex4->SetTextFont(43);
    latex4->SetTextColor(kBlack);
    latex4->Draw("same");


    //c3->Print("../figures/CME_figure4.pdf");
}

