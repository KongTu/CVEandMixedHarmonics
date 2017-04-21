#include "RiceStyle.h"

using namespace std;

double PbPb_centralityBinCenter[] = {75,67.5,62.5,57.5,52.5,47.5,42.5,37.5,32.5};
double PbPb_centralityBinCenter_tracker[] = {65, 55, 45, 35};

double total_systematics_pPb = 0.000041;
double total_systematics_PbPb = 0.000028;


void plotCMEResults_differential_diff(){

	gStyle->SetErrorX(0);

/*
18 histograms with 6 in deta, 6 in dpt, and 6 inptAve
*/
	TFile* file_differential_112 = new TFile("../dataPoints/CME_differentials_112_pPb_8TeV.root");
	TFile* file_differential_123 = new TFile("../dataPoints/CME_differentials_123_pPb_8TeV.root");
	TFile* file_differential_PbPb_112 = new TFile("../dataPoints/CME_differentials_112_PbPb_5TeV.root");
	TFile* file_differential_PbPb_123 = new TFile("../dataPoints/CME_differentials_123_PbPb_5TeV.root");

//for correlator with 112 coefficients
/*6 in deta:
	- temp1: 3p Pb-going SS
	- temp2: 3p Pb-going OS
	- temp3: 3p p-going SS
	- temp4: 3p p-going OS
	- temp5: 2p correlator SS
	- temp6: 2p correlator OS
*/

 	TH1D* temp1 = (TH1D*) file_differential_112->Get("temp1");
 	TH1D* temp2 = (TH1D*) file_differential_112->Get("temp2");
 	TH1D* temp3 = (TH1D*) file_differential_112->Get("temp3");
 	TH1D* temp4 = (TH1D*) file_differential_112->Get("temp4");
 	TH1D* temp5 = (TH1D*) file_differential_112->Get("temp5");
 	TH1D* temp6 = (TH1D*) file_differential_112->Get("temp6");

//6 in dpT, similar structure as in deta
 	TH1D* temp11 = (TH1D*) file_differential_112->Get("temp11");
 	TH1D* temp12 = (TH1D*) file_differential_112->Get("temp12");
 	TH1D* temp13 = (TH1D*) file_differential_112->Get("temp13");
 	TH1D* temp14 = (TH1D*) file_differential_112->Get("temp14");
 	TH1D* temp15 = (TH1D*) file_differential_112->Get("temp15");
 	TH1D* temp16 = (TH1D*) file_differential_112->Get("temp16");

//6 in pTave, similar structure as in deta
 	TH1D* temp21 = (TH1D*) file_differential_112->Get("temp21");
 	TH1D* temp22 = (TH1D*) file_differential_112->Get("temp22");
 	TH1D* temp23 = (TH1D*) file_differential_112->Get("temp23");
 	TH1D* temp24 = (TH1D*) file_differential_112->Get("temp24");
 	TH1D* temp25 = (TH1D*) file_differential_112->Get("temp25");
 	TH1D* temp26 = (TH1D*) file_differential_112->Get("temp26");


//for correlator with 112 coefficients
//same structure as above

	TH1D* temp123_1 = (TH1D*) file_differential_123->Get("temp1");
 	TH1D* temp123_2 = (TH1D*) file_differential_123->Get("temp2");
 	TH1D* temp123_3 = (TH1D*) file_differential_123->Get("temp3");
 	TH1D* temp123_4 = (TH1D*) file_differential_123->Get("temp4");
 	TH1D* temp123_5 = (TH1D*) file_differential_123->Get("temp5");
 	TH1D* temp123_6 = (TH1D*) file_differential_123->Get("temp6");

 	TH1D* temp123_11 = (TH1D*) file_differential_123->Get("temp11");
 	TH1D* temp123_12 = (TH1D*) file_differential_123->Get("temp12");
 	TH1D* temp123_13 = (TH1D*) file_differential_123->Get("temp13");
 	TH1D* temp123_14 = (TH1D*) file_differential_123->Get("temp14");
 	TH1D* temp123_15 = (TH1D*) file_differential_123->Get("temp15");
 	TH1D* temp123_16 = (TH1D*) file_differential_123->Get("temp16");

 	TH1D* temp123_21 = (TH1D*) file_differential_123->Get("temp21");
 	TH1D* temp123_22 = (TH1D*) file_differential_123->Get("temp22");
 	TH1D* temp123_23 = (TH1D*) file_differential_123->Get("temp23");
 	TH1D* temp123_24 = (TH1D*) file_differential_123->Get("temp24");
 	TH1D* temp123_25 = (TH1D*) file_differential_123->Get("temp25");
 	TH1D* temp123_26 = (TH1D*) file_differential_123->Get("temp26");


/*
PbPb
*/

 	TH1D* tempPbPb1 = (TH1D*) file_differential_PbPb_112->Get("temp1");
 	TH1D* tempPbPb2 = (TH1D*) file_differential_PbPb_112->Get("temp2");
 
 	TH1D* tempPbPb5 = (TH1D*) file_differential_PbPb_112->Get("temp5");
 	TH1D* tempPbPb6 = (TH1D*) file_differential_PbPb_112->Get("temp6");

 	TH1D* tempPbPb11 = (TH1D*) file_differential_PbPb_112->Get("temp11");
 	TH1D* tempPbPb12 = (TH1D*) file_differential_PbPb_112->Get("temp12");

 	TH1D* tempPbPb15 = (TH1D*) file_differential_PbPb_112->Get("temp15");
 	TH1D* tempPbPb16 = (TH1D*) file_differential_PbPb_112->Get("temp16");

 	TH1D* tempPbPb21 = (TH1D*) file_differential_PbPb_112->Get("temp21");
 	TH1D* tempPbPb22 = (TH1D*) file_differential_PbPb_112->Get("temp22");

 	TH1D* tempPbPb25 = (TH1D*) file_differential_PbPb_112->Get("temp25");
 	TH1D* tempPbPb26 = (TH1D*) file_differential_PbPb_112->Get("temp26");

//123 coeffients. 
 	TH1D* tempPbPb123_1 = (TH1D*) file_differential_PbPb_123->Get("temp1");
 	TH1D* tempPbPb123_2 = (TH1D*) file_differential_PbPb_123->Get("temp2");
 
 	TH1D* tempPbPb123_5 = (TH1D*) file_differential_PbPb_123->Get("temp5");
 	TH1D* tempPbPb123_6 = (TH1D*) file_differential_PbPb_123->Get("temp6");

 	TH1D* tempPbPb123_11 = (TH1D*) file_differential_PbPb_123->Get("temp11");
 	TH1D* tempPbPb123_12 = (TH1D*) file_differential_PbPb_123->Get("temp12");

 	TH1D* tempPbPb123_15 = (TH1D*) file_differential_PbPb_123->Get("temp15");
 	TH1D* tempPbPb123_16 = (TH1D*) file_differential_PbPb_123->Get("temp16");

 	TH1D* tempPbPb123_21 = (TH1D*) file_differential_PbPb_123->Get("temp21");
 	TH1D* tempPbPb123_22 = (TH1D*) file_differential_PbPb_123->Get("temp22");

 	TH1D* tempPbPb123_25 = (TH1D*) file_differential_PbPb_123->Get("temp25");
 	TH1D* tempPbPb123_26 = (TH1D*) file_differential_PbPb_123->Get("temp26");


//HIN-16-009 data points. 
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
plotting deltaEta dependence, Fig.1
*/

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "", "|#Delta#eta|", "#gamma_{112}(os-ss)", 48,0,4.8,kBlack);
	base1->GetYaxis()->SetRangeUser(-0.0002,0.003);

	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.0,1.25);

	base1->GetYaxis()->SetTitleOffset(1.4);
	base1->GetXaxis()->SetTitleOffset(1.2);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetNdivisions(5,6,0);

	TH1D* base1_1 = base1->Clone("base1_1");
	base1_1->GetXaxis()->SetRangeUser(-0.2,3.0);
	base1_1->GetXaxis()->SetTitle("|#Deltap_{T}|");
	base1_1->GetYaxis()->SetRangeUser(-0.0002,0.003);

	TH1D* base1_2 = base1->Clone("base1_2");
	base1_2->GetXaxis()->SetRangeUser(-0.2,3.0);
	base1_2->GetYaxis()->SetRangeUser(-0.0002,0.003);

	base1_2->GetXaxis()->SetTitle("(p_{T,#alpha}+p_{T,#beta})/2");

    TLatex* r41 = new TLatex(0.41, 0.87, "pPb #sqrt{s_{NN}} = 8.16 TeV");
    r41->SetNDC();
    r41->SetTextSize(23);
    r41->SetTextFont(43);
    r41->SetTextColor(kBlack);

    TLatex* r42 = new TLatex(0.43, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.47,0.95, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.06);
    
    TLatex* r44 = new TLatex(0.61,0.95, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(20);
    r44->SetTextFont(53);
    
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


//begin figure1

    TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,460);
    c1->Divide(3,1,0,0);
	c1->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base1->Draw();

	temp2->SetMarkerStyle(20);
	temp2->SetLineColor(kRed);
	temp2->SetMarkerColor(kRed);
	temp2->Add(temp1, -1);
	temp2->Draw("Psame");
	temp4->SetMarkerStyle(21);
	temp4->SetMarkerColor(kBlue);
	temp4->Add(temp3, -1);
	temp4->Draw("Psame");

	tempPbPb2->SetMarkerStyle(28);
	tempPbPb2->SetMarkerColor(kBlack);
	tempPbPb2->SetLineColor(kBlack);
	tempPbPb2->Add(tempPbPb1, -1);
	tempPbPb2->Draw("Psame");

	r41->Draw("same");
    r46->Draw("same");
    

	c1->cd(2);
	gPad->SetLeftMargin(0);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base1_1->Draw();

	temp12->SetMarkerStyle(20);
	temp12->SetMarkerColor(kRed);
	temp12->SetLineColor(kRed);
	temp12->Add(temp11, -1);
	temp12->Draw("Psame");
	
	temp14->SetMarkerStyle(21);
	temp14->SetMarkerColor(kBlue);
	temp14->Add(temp13, -1);
	temp14->Draw("Psame");

	tempPbPb12->SetMarkerStyle(28);
	tempPbPb12->SetMarkerColor(kBlack);
	tempPbPb12->SetLineColor(kBlack);
	tempPbPb12->Add(tempPbPb11, -1);
	tempPbPb12->Draw("Psame");

	TLegend *w4 = new TLegend(0.1,0.7,0.6,0.9);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(21);
    w4->SetTextFont(45);
    w4->AddEntry(temp2, "pPb 8TeV, #phi_{c}(Pb-going)","P");
    w4->AddEntry(temp4, "pPb 8TeV, #phi_{c}(p-going)","P");
    w4->AddEntry(tempPbPb12, "PbPb 5TeV","P");

	w4->Draw("same");

	c1->cd(3);
	gPad->SetLeftMargin(0);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.1);
	gPad->SetTicks();
	base1_2->Draw();
	temp22->SetMarkerStyle(20);
	temp22->SetMarkerColor(kRed);
	temp22->SetLineColor(kRed);
	temp22->Add(temp21, -1);
	temp22->Draw("Psame");
	temp24->SetMarkerStyle(21);
	temp24->SetMarkerColor(kBlue);
	temp24->Add(temp23, -1);
	temp24->Draw("Psame");

	tempPbPb22->SetMarkerStyle(28);
	tempPbPb22->SetMarkerColor(kBlack);
	tempPbPb22->SetLineColor(kBlack);
	tempPbPb22->Add(tempPbPb21, -1);
	tempPbPb22->Draw("Psame");

	r43->Draw("same");
    r44->Draw("same");

//end of figure1


    TH1D* base2 = makeHist("base2", "", "|#Delta#eta|", "#delta(os-ss)", 48,0,4.8,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.0005,0.012);
	base2->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.0,1.25);

	base2->GetYaxis()->SetTitleOffset(1.4);
	base2->GetXaxis()->SetTitleOffset(1.2);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetNdivisions(5,6,0);

	TH1D* base2_1 = base2->Clone("base2_1");
	base2_1->GetXaxis()->SetRangeUser(0,3.0);
	base2_1->GetYaxis()->SetRangeUser(-0.0005,0.012);
	base2_1->GetXaxis()->SetTitle("|#Deltap_{T}|");

	TH1D* base2_2 = base2->Clone("base2_2");
	base2_2->GetXaxis()->SetRangeUser(0,3.0);
	base2_2->GetYaxis()->SetRangeUser(-0.0005,0.012);
	base2_2->GetXaxis()->SetTitle("(p_{T,#alpha}+p_{T,#beta})/2");


//delta scale by v2.

	TCanvas* c2 = new TCanvas("c2","c2",1,1,1000,460);
    c2->Divide(3,1,0,0);
	c2->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base2->Draw();
	
	temp6->SetMarkerStyle(20);
	temp6->SetMarkerColor(kRed);
	temp6->SetLineColor(kRed);
	temp6->Add(temp5, -1);
	temp6->Draw("Psame");

	tempPbPb6->SetMarkerStyle(28);
	tempPbPb6->SetMarkerColor(kBlack);
	tempPbPb6->SetLineColor(kBlack);
	tempPbPb6->Add(tempPbPb5, -1);
	tempPbPb6->Draw("Psame");
	
	r41->Draw("same");    
    r46->Draw("same");
    

	c2->cd(2);
	gPad->SetLeftMargin(0.0);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base2_1->Draw();

	temp16->SetMarkerStyle(20);
	temp16->SetMarkerColor(kRed);
	temp16->SetLineColor(kRed);
	temp16->Add(temp15, -1);
	temp16->Draw("Psame");

	tempPbPb16->SetMarkerStyle(28);
	tempPbPb16->SetMarkerColor(kBlack);
	tempPbPb16->SetLineColor(kBlack);
	tempPbPb16->Add(tempPbPb15, -1);
	tempPbPb16->Draw("Psame");

	TLegend *w5 = new TLegend(0.1,0.75,0.6,0.9);
    w5->SetLineColor(kWhite);
    w5->SetFillColor(0);
    w5->SetTextSize(21);
    w5->SetTextFont(45);
    w5->AddEntry(temp16, " pPb 8TeV","P");
    w5->AddEntry(tempPbPb16, " PbPb 8TeV","P");
	w5->Draw("same");

	c2->cd(3);
	gPad->SetLeftMargin(0.0);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.1);
	gPad->SetTicks();
	base2_2->Draw();

	temp26->SetMarkerStyle(20);
	temp26->SetMarkerColor(kRed);
	temp26->SetLineColor(kRed);
	temp26->Add(temp25, -1);
	temp26->Draw("Psame");

	tempPbPb26->SetMarkerStyle(28);
	tempPbPb26->SetMarkerColor(kBlack);
	tempPbPb26->SetLineColor(kBlack);
	tempPbPb26->Add(tempPbPb25, -1);
	tempPbPb26->Draw("Psame");

	r43->Draw("same");
    r44->Draw("same");

	TH1D* base3 = makeHist("base3", "", "|#Delta#eta|", "#gamma_{123}(os-ss)", 48,0,4.8,kBlack);
	base3->GetYaxis()->SetRangeUser(-0.0005,0.0012);
	base3->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.0,1.25);

	base3->GetYaxis()->SetTitleOffset(1.4);
	base3->GetXaxis()->SetTitleOffset(1.2);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.3);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.4);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.4);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.4);
	base3->GetXaxis()->SetNdivisions(5,6,0);

	TH1D* base3_1 = base3->Clone("base3_1");
	base3_1->GetXaxis()->SetRangeUser(0,3.0);
	base3_1->GetYaxis()->SetRangeUser(-0.0005,0.0012);
	base3_1->GetXaxis()->SetTitle("|#Deltap_{T}|");

	TH1D* base3_2 = base3->Clone("base3_2");
	base3_2->GetXaxis()->SetRangeUser(0,3.0);
	base3_2->GetYaxis()->SetRangeUser(-0.0005,0.0012);
	base3_2->GetXaxis()->SetTitle("(p_{T,#alpha}+p_{T,#beta})/2");


	TCanvas* c3 = new TCanvas("c3","c3",1,1,1000,460);
    c3->Divide(3,1,0,0);
	c3->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base3->Draw();
	
	temp123_2->SetMarkerStyle(20);
	temp123_2->SetMarkerColor(kRed);
	temp123_2->SetLineColor(kRed);
	temp123_2->Add(temp123_1, -1);	
	temp123_2->Draw("Psame");
	temp123_4->SetMarkerStyle(21);
	temp123_4->SetMarkerColor(kBlue);
	temp123_4->Add(temp123_3, -1);
	temp123_4->Draw("Psame");

	tempPbPb123_2->SetMarkerStyle(28);
	tempPbPb123_2->SetMarkerColor(kBlack);
	tempPbPb123_2->SetLineColor(kBlack);
	tempPbPb123_2->Add(tempPbPb123_1, -1);
	tempPbPb123_2->Draw("Psame");
	
	r41->Draw("same");
    r46->Draw("same");
   

	c3->cd(2);
	gPad->SetLeftMargin(0.0);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base3_1->Draw();

	temp123_12->SetMarkerStyle(20);
	temp123_12->SetMarkerColor(kRed);
	temp123_12->SetLineColor(kRed);
	temp123_12->Add(temp123_11, -1);
	temp123_12->Draw("Psame");
	temp123_14->SetMarkerStyle(21);
	temp123_14->SetMarkerColor(kBlue);
	temp123_14->Add(temp123_13, -1);
	temp123_14->Draw("Psame");

	tempPbPb123_12->SetMarkerStyle(28);
	tempPbPb123_12->SetMarkerColor(kBlack);
	tempPbPb123_12->SetLineColor(kBlack);
	tempPbPb123_12->Add(tempPbPb123_11, -1);
	tempPbPb123_12->Draw("Psame");

	TLegend *w6 = new TLegend(0.15,0.7,0.60,0.9);
    w6->SetLineColor(kWhite);
    w6->SetFillColor(0);
    w6->SetTextSize(21);
    w6->SetTextFont(45);
    w6->AddEntry(temp123_12, "pPb 8TeV, #phi_{c}(Pb-going)","P");
    w6->AddEntry(temp123_14, "pPb 8TeV, #phi_{c}(p-going)","P");
	w6->AddEntry(tempPbPb123_12, "PbPb 5TeV","P");
	w6->Draw("same");

	c3->cd(3);
	gPad->SetLeftMargin(0.0);
	gPad->SetBottomMargin(0.17);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.1);
	gPad->SetTicks();
	base3_2->Draw();

	temp123_22->SetMarkerStyle(20);
	temp123_22->SetMarkerColor(kRed);
	temp123_22->SetLineColor(kRed);
	temp123_22->Add(temp123_21, -1);
	temp123_22->Draw("Psame");
	temp123_24->SetMarkerStyle(21);
	temp123_24->SetMarkerColor(kBlue);
	temp123_24->Add(temp123_23, -1);
	temp123_24->Draw("Psame");
 	
	tempPbPb123_22->SetMarkerStyle(28);
	tempPbPb123_22->SetMarkerColor(kBlack);
	tempPbPb123_22->SetLineColor(kBlack);
	tempPbPb123_22->Add(tempPbPb123_21, -1);
	tempPbPb123_22->Draw("Psame");

 	r43->Draw("same");
    r44->Draw("same");
	TGaxis::SetMaxDigits(3);
	// c1->Print("../figures/CME_figure_2_a.pdf");
	// c2->Print("../figures/CME_figure_2_b.pdf");
	// c3->Print("../figures/CME_figure_2_c.pdf");


}