#include "RiceStyle.h"

using namespace std;


void plotCME_ESE(){

	TFile* file[2];

	file[0] = new TFile("../dataPoints/pPb_ESE_points.root");
	file[1] = new TFile("../dataPoints/PbPb_ESE_points.root");

	TGraphErrors* gr1_pPb = (TGraphErrors*) file[0]->Get("Graph;1");
	TGraphErrors* gr2_pPb = (TGraphErrors*) file[0]->Get("Graph;2");

	TF1 * f1_pPb = (TF1*) file[0]->Get("pol1;1");
	TF1 * f2_pPb = (TF1*) file[0]->Get("pol1;2");

	TGraphErrors* gr1_PbPb = (TGraphErrors*) file[1]->Get("Graph;1");
	TF1 * f1_PbPb = (TF1*) file[1]->Get("pol1;1");

	TGaxis::SetMaxDigits(3);

	TCanvas* c3 = new TCanvas("c3","c3",700,700);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "Pb-going", "v_{2,tracker}", "#gamma_{112} (os-ss)", 1000,0,1.0,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.0004, 0.0015);
	base2->GetXaxis()->SetRangeUser(0.0, 0.13);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.6);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.6);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetNdivisions(4,6,0);
	base2->GetYaxis()->SetNdivisions(4,6,0);

	base2->Draw();

	gr1_PbPb->SetMarkerColor(kBlue);
	gr1_PbPb->SetLineColor(kBlue);
	f1_PbPb->SetLineColor(kBlue);

    gr1_pPb->SetMarkerSize(1.6);
	gr1_pPb->Draw("Psame");
	f1_pPb->Draw("same");

    gr1_PbPb->SetMarkerSize(1.6);
	gr1_PbPb->Draw("Psame");
	f1_PbPb->Draw("same");

 	f1_pPb->SetLineStyle(2);
    double intersect_1 = f1_pPb->GetParameter(0);
    double intersect_1_error = f1_pPb->GetParError(0);
    double slope_1 = f1_pPb->GetParameter(1);
    double slope_1_error = f1_pPb->GetParError(1);
    f1_pPb->SetRange(0,1);

    TLatex* latex3 = new TLatex(0.18, 0.66, Form("slope: %.5f +/- %.5f",slope_1, slope_1_error ));
    latex3->SetNDC();
    latex3->SetTextSize(20);
    latex3->SetTextFont(43);
    latex3->SetTextColor(kRed);
    latex3->Draw("same");

    TLatex* latex4 = new TLatex(0.18, 0.63, Form("intersect: %.5f +/- %.5f",intersect_1, intersect_1_error ));
    latex4->SetNDC();
    latex4->SetTextSize(20);
    latex4->SetTextFont(43);
    latex4->SetTextColor(kRed);
    latex4->Draw("same");

	f1_PbPb->SetLineColor(kBlue);
    f1_PbPb->SetLineStyle(2);
    double intersect_2 = f1_PbPb->GetParameter(0);
    double intersect_2_error = f1_PbPb->GetParError(0);
    double slope_2 = f1_PbPb->GetParameter(1);
    double slope_2_error = f1_PbPb->GetParError(1);
    f1_PbPb->SetRange(0,1);

    TLatex* latex5 = new TLatex(0.18, 0.59, Form("slope: %.5f +/- %.5f",slope_2, slope_2_error ));
    latex5->SetNDC();
    latex5->SetTextSize(20);
    latex5->SetTextFont(43);
    latex5->SetTextColor(kBlue);
    latex5->Draw("same");

    TLatex* latex6 = new TLatex(0.18, 0.56, Form("intersect: %.5f +/- %.5f",intersect_2, intersect_2_error ));
    latex6->SetNDC();
    latex6->SetTextSize(20);
    latex6->SetTextFont(43);
    latex6->SetTextColor(kBlue);
    latex6->Draw("same");

    TLatex* r11 = new TLatex(0.18,0.84, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);

    TLatex* r22 = new TLatex(0.27,0.84, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);

    TLatex* lmult = new TLatex(0.18, 0.78, "185 #leq N^{offline}_{trk} < 250");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);

    TLatex* r33 = new TLatex(0.18, 0.72, "|#eta_{#alpha} - #eta_{#beta}| < 1.6");
    r33->SetNDC();
    r33->SetTextSize(23);
    r33->SetTextFont(43);
    r33->SetTextColor(kBlack);
    r33->Draw("same");

    r11->Draw("same");
    //r22->Draw("same");
    lmult->Draw("same");

    TLegend *w1 = new TLegend(0.5,0.16,0.7,0.26);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(21);
    w1->SetTextFont(43);
    w1->AddEntry(gr1_pPb, "pPb 8 TeV", "P");
    w1->AddEntry(gr1_PbPb, "PbPb 5 TeV", "P");
    w1->Draw("same");

    //c3->Print("../figures/CME_figure9.pdf");
}