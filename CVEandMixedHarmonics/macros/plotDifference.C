#include "RiceStyle.h"

using namespace std;

void plotDifference(){
	

	TFile* file1 = new TFile("../dataPoints/pPb_Ntrk.root");
	TFile* file2 = new TFile("../dataPoints/PbPb_Cent.root");
	TFile* file3 = new TFile("../dataPoints/CME_data.root");
	TFile* file4 = new TFile("../dataPoints/pPb_Ntrk_K0s.root");


	TGraphErrors* CME_pPb_Pb = (TGraphErrors*) file3->Get("Graph;1");
	TGraphErrors* CME_pPb_p = (TGraphErrors*) file3->Get("Graph;2");
	TGraphErrors* CME_PbPb = (TGraphErrors*) file3->Get("Graph;3");

	TGraphErrors* gr1_Pb = (TGraphErrors*) file1->Get("Graph;1");
	TGraphErrors* gr1_p = (TGraphErrors*) file1->Get("Graph;2");
	TGraphErrors* gr2 = (TGraphErrors*) file2->Get("Graph");

	TGraphErrors* gr2_K0s = (TGraphErrors*) file4->Get("Graph;1");
	TGraphErrors* gr2_Lam = (TGraphErrors*) file4->Get("Graph;2");
	TGraphErrors* gr2_LamH = (TGraphErrors*) file4->Get("Graph;3");
	TGraphErrors* gr2_K0sH = (TGraphErrors*) file4->Get("Graph;4");

	TH1D* base2 = makeHist("base2", "", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c} (#Lambda#bar{#Lambda} - #Lambda#Lambda/#bar{#Lambda}#bar{#Lambda})", 1500,0,1500,kBlack);

	base2->GetYaxis()->SetRangeUser(-0.05, 0.05);
	base2->GetXaxis()->SetRangeUser(55,1500);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.23);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetNdivisions(505);
	base2->GetYaxis()->SetNdivisions(505);

	TCanvas* c2 = new TCanvas("c2","c2",1,1,650,650);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx();

	base2->Draw();
	//gr1_p->Draw("Psame");
	gr1_Pb->Draw("Psame");
	gr2->Draw("Psame");

	CME_pPb_Pb->SetMarkerStyle(25);
	CME_pPb_Pb->SetMarkerColor(kBlue);
	CME_pPb_Pb->SetLineColor(kBlue);

	CME_pPb_Pb->Draw("Psame");
	CME_pPb_p->SetMarkerStyle(25);
	//CME_pPb_p->Draw("Psame");

	gr2_K0s->SetMarkerStyle(28);
	// gr2_K0s->Draw("Psame");
	// gr2_Lam->Draw("Psame");
	// gr2_LamH->Draw("Psame");
	// gr2_K0sH->Draw("Psame");

	TLegend *w1 = new TLegend(0.47,0.2,0.82,0.33);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(19);
    w1->SetTextFont(43);
    w1->AddEntry(gr1_Pb, "CVE pPb 8 TeV, Pb-going", "P");
    //w1->AddEntry(gr1_p, "pPb 8 TeV, p-going", "P");
    w1->AddEntry(gr2, "CVE PbPb 5 TeV", "P");
    w1->AddEntry(CME_pPb_Pb, "CME pPb 5 TeV, Pb-going", "P");
    //w1->AddEntry(CME_pPb_p, "PbPb 5 TeV", "P");
    w1->Draw("same");

	TLatex* r33 = new TLatex(0.16,0.91, "CMS");
	r33->SetNDC();
	r33->SetTextSize(0.045);
	r33->Draw("same");

	TLatex* r44 = new TLatex(0.26,0.91, "Preliminary");
	r44->SetNDC();
	r44->SetTextSize(24);
	r44->SetTextFont(53);
	r44->Draw("same");

	TLatex* r1 = new TLatex(0.18,0.83, "K^{0}_{s} (0.3,5.0) GeV/c");
	r1->SetNDC();
	r1->SetTextSize(24);
	r1->SetTextFont(44);
	//r1->Draw("same");

	TLatex* r2 = new TLatex(0.18,0.83, "#Lambda (0.6,5.0) GeV/c");
	r2->SetNDC();
	r2->SetTextSize(24);
	r2->SetTextFont(44);
	r2->Draw("same");








}