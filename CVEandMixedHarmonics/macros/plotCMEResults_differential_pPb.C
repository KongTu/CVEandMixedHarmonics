#include "RiceStyle.h"

using namespace std;

double PbPb_centralityBinCenter[] = {75,67.5,62.5,57.5,52.5,47.5,42.5,37.5,32.5};
double PbPb_centralityBinCenter_tracker[] = {65, 55, 45, 35};

double total_systematics_pPb = 0.000041;
double total_systematics_PbPb = 0.000028;


void plotCMEResults_differential_pPb(){

	gStyle->SetErrorX(0);

/*
18 histograms with 6 in deta, 6 in dpt, and 6 inptAve
*/
	TFile* file_differential_112 = new TFile("../dataPoints/CME_differentials_112_pPb_8TeV.root");
	TFile* file_differential_123 = new TFile("../dataPoints/CME_differentials_123_pPb_8TeV.root");

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

	TH1D* base1 = makeHist("base1", "", "|#Delta#eta|", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 48,0,4.8,kBlack);
	base1->GetYaxis()->SetRangeUser(-0.0013,0.0012);

	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetNdivisions(5,6,0);

	TH1D* base1_1 = base1->Clone("base1_1");
	base1_1->GetXaxis()->SetRangeUser(0,3.0);
	base1_1->GetXaxis()->SetTitle("|#Deltap_{T}|");
	base1_1->GetYaxis()->SetRangeUser(-0.0005,0.001);

	TH1D* base1_2 = base1->Clone("base1_2");
	base1_2->GetXaxis()->SetRangeUser(0,3.0);
	base1_2->GetYaxis()->SetRangeUser(-0.003,0.002);

	base1_2->GetXaxis()->SetTitle("(p_{T,#alpha}+p_{T,#beta})/2");

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
    r43->SetTextSize(0.05);
    
    TLatex* r44 = new TLatex(0.72,0.95, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(21);
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

	TCanvas* c1 = makeMultiCanvas("c1","c1", 3,1);
	c1->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base1->Draw();
	temp1->Draw("Psame");
	temp2->Draw("Psame");
	temp3->Draw("Psame");
	temp4->Draw("Psame");

	r41->Draw("same");
    
    r46->Draw("same");
    r43->Draw("same");
    r44->Draw("same");

	c1->cd(2);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base1_1->Draw();
	temp11->Draw("Psame");
	temp12->Draw("Psame");
	temp13->Draw("Psame");
	temp14->Draw("Psame");

	c1->cd(3);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base1_2->Draw();
	temp21->Draw("Psame");
	temp22->Draw("Psame");
	temp23->Draw("Psame");
	temp24->Draw("Psame");

	TLegend *w4 = new TLegend(0.25,0.17,0.70,0.3);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(23);
    w4->SetTextFont(45);
    w4->SetNColumns(2);
    w4->AddEntry(temp1, "  ", "P");
    w4->AddEntry(temp2, "  #phi_{c}(Pb-going)","P");
    
    w4->AddEntry(temp3, "  ", "P");
    w4->AddEntry(temp4, "  #phi_{c}(p-going)","P");

	TLatex* latex1 = new TLatex(0.23, 0.31, "  SS");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    TLatex* latex2 = new TLatex(0.33, 0.31, "OS");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);

	w4->Draw("same");
    latex1->Draw("same");
    latex2->Draw("same");
//end of figure1


    TH1D* base2 = makeHist("base2", "", "|#Delta#eta|", "#LTcos(#phi_{#alpha}-#phi_{#beta})#GT", 48,0,4.8,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.007,0.02);
	base2->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetXaxis()->SetTitleOffset(0.95);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetNdivisions(5,6,0);

	TH1D* base2_1 = base2->Clone("base2_1");
	base2_1->GetXaxis()->SetRangeUser(0,3.0);
	base2_1->GetYaxis()->SetRangeUser(-0.007,0.02);
	base2_1->GetXaxis()->SetTitle("|#Deltap_{T}|");

	TH1D* base2_2 = base2->Clone("base2_2");
	base2_2->GetXaxis()->SetRangeUser(0,3.0);
	base2_2->GetYaxis()->SetRangeUser(-0.01,0.01);
	base2_2->GetXaxis()->SetTitle("(p_{T,#alpha}+p_{T,#beta})/2");


    TCanvas* c2 = makeMultiCanvas("c2","c2", 3,1);
	c2->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base2->Draw();
	temp5->Draw("Psame");
	temp6->Draw("Psame");
	
	r41->Draw("same");    
    r46->Draw("same");
    r43->Draw("same");
    r44->Draw("same");

	c2->cd(2);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base2_1->Draw();

	temp15->Draw("Psame");
	temp16->Draw("Psame");

	c2->cd(3);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base2_2->Draw();

	temp25->Draw("Psame");
	temp26->Draw("Psame");

	TLegend *w5 = new TLegend(0.27,0.17,0.45,0.3);
    w5->SetLineColor(kWhite);
    w5->SetFillColor(0);
    w5->SetTextSize(23);
    w5->SetTextFont(45);
    w5->SetNColumns(2);
    w5->AddEntry(temp5, "  ", "P");
    // w5->AddEntry(temp5, "  #phi_{c}(Pb-going)","P");
    
    w5->AddEntry(temp6, "  ", "P");
    // w5->AddEntry(temp6, "  #phi_{c}(p-going)","P");

	TLatex* latex1 = new TLatex(0.22, 0.31, "  SS");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    TLatex* latex2 = new TLatex(0.33, 0.31, "OS");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);

	w5->Draw("same");
    latex1->Draw("same");
    latex2->Draw("same");

	TH1D* base3 = makeHist("base3", "", "|#Delta#eta|", "#LTcos(#phi_{#alpha}+2#phi_{#beta}-3#phi_{c})#GT/v_{3,c}", 48,0,4.8,kBlack);
	base3->GetYaxis()->SetRangeUser(-0.002,0.0006);
	base3->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.3);
	base3->GetXaxis()->SetTitleOffset(0.95);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.3);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.4);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.4);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.4);
	base3->GetXaxis()->SetNdivisions(5,6,0);

	TH1D* base3_1 = base3->Clone("base3_1");
	base3_1->GetXaxis()->SetRangeUser(0,3.0);
	base3_1->GetYaxis()->SetRangeUser(-0.003,0.0003);
	base3_1->GetXaxis()->SetTitle("|#Deltap_{T}|");

	TH1D* base3_2 = base3->Clone("base3_2");
	base3_2->GetXaxis()->SetRangeUser(0,3.0);
	base3_2->GetYaxis()->SetRangeUser(-0.01,0.002);
	base3_2->GetXaxis()->SetTitle("(p_{T,#alpha}+p_{T,#beta})/2");


    TCanvas* c3 = makeMultiCanvas("c3","c3", 3,1);
	c3->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base3->Draw();
	temp123_1->Draw("Psame");
	temp123_2->Draw("Psame");
	temp123_3->Draw("Psame");
	temp123_4->Draw("Psame");
	
	r41->Draw("same");
    r46->Draw("same");
    r43->Draw("same");
    r44->Draw("same");

	c3->cd(2);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base3_1->Draw();

	temp123_11->Draw("Psame");
	temp123_12->Draw("Psame");
	temp123_13->Draw("Psame");
	temp123_14->Draw("Psame");

	c3->cd(3);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.0);
	gPad->SetTicks();
	base3_2->Draw();

	temp123_21->Draw("Psame");
	temp123_22->Draw("Psame");
	temp123_23->Draw("Psame");
	temp123_24->Draw("Psame");

	TLegend *w6 = new TLegend(0.25,0.17,0.70,0.3);
    w6->SetLineColor(kWhite);
    w6->SetFillColor(0);
    w6->SetTextSize(23);
    w6->SetTextFont(45);
    w6->SetNColumns(2);
    w6->AddEntry(temp1, "  ", "P");
    w6->AddEntry(temp2, "  #phi_{c}(Pb-going)","P");
    
    w6->AddEntry(temp3, "  ", "P");
    w6->AddEntry(temp4, "  #phi_{c}(p-going)","P");

	TLatex* latex1 = new TLatex(0.23, 0.31, "  SS");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    TLatex* latex2 = new TLatex(0.33, 0.31, "OS");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);

	w6->Draw("same");
    latex1->Draw("same");
    latex2->Draw("same");

	TGaxis::SetMaxDigits(3);
	c1->Print("../figures/CME_figure_1_a.pdf");
	c2->Print("../figures/CME_figure_1_b.pdf");
	c3->Print("../figures/CME_figure_1_c.pdf");




}