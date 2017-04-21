#include "RiceStyle.h"

using namespace std;

void findQ2(){


	TFile* file = new TFile("../rootfiles/CMEandMixedHarmonics_PbPb_Centrality_full_v5_5.root");
	TH1D* q2_mag = (TH1D*) file->Get("ana/q2_mag");
	
	q2_mag->Draw();

	double total = q2_mag->Integral();

	cout << "total: " << total << endl;
	
	double ratio_value[] = {1.0, 0.95, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.00095};
	for(int i = 0; i < q2_mag->GetNbinsX(); i++){

		double num = q2_mag->Integral(i+1, q2_mag->GetNbinsX());
		double ratio = num/total;

		for(int j = 0; j < 11; j++){
			if( ratio < ratio_value[j] && ratio > ratio_value[j+1] ){

				cout << "ESE top class " << j+1 << endl;
				cout << "Q2 value " << q2_mag->GetBinCenter( i+1 ) << endl;
				cout << "bin number: " << i+1 << endl;
				cout << "ratio value : " << q2_mag->Integral(i+1, 1000)/total << endl;

	
			}
		}

	}

	TCanvas* c1 = new TCanvas("c1","c1",1,1,700,700);
	q2_mag->SetMarkerColor(kBlack);
	q2_mag->SetLineColor(kBlack);
	q2_mag->GetXaxis()->SetRangeUser(0,0.65);
	q2_mag->SetStats(kFALSE);
	q2_mag->SetTitle(" ");
	q2_mag->GetXaxis()->SetTitle("Q_{2}");
	q2_mag->GetXaxis()->CenterTitle();
	q2_mag->GetXaxis()->SetTitleSize(0.04);
	q2_mag->GetYaxis()->SetRangeUser(0.1, 1000000000);
	q2_mag->Draw();

	gPad->SetTicks();
	gPad->SetLogy(1);

	TLine* l1[10];
	double temp = q2_mag->GetBinContent(1029);
    l1[0] = new TLine(0.0285,0, 0.0285, temp);

    double temp = q2_mag->GetBinContent(1060);
    l1[1] = new TLine(0.0595,0, 0.0595, temp);

    double temp = q2_mag->GetBinContent(1091);
    l1[2] = new TLine(0.0905,0, 0.0905, temp);

    double temp = q2_mag->GetBinContent(1106);
    l1[3] = new TLine(0.1055,0, 0.1055, temp);

    double temp = q2_mag->GetBinContent(1122);
    l1[4] = new TLine(0.1215,0, 0.1215, temp);

    double temp = q2_mag->GetBinContent(1139);
    l1[5] = new TLine(0.1385,0, 0.1385, temp);

    double temp = q2_mag->GetBinContent(1161);
    l1[6] = new TLine(0.1605,0, 0.1605, temp);

    double temp = q2_mag->GetBinContent(1193);
    l1[7] = new TLine(0.1925,0, 0.1925, temp);

    double temp = q2_mag->GetBinContent(1220);
    l1[8] = new TLine(0.2195,0, 0.2195, temp);

    double temp = q2_mag->GetBinContent(1274);
    l1[9] = new TLine(0.2735,0, 0.2735, temp);

   	TLatex* r1 = new TLatex(0.55,0.84, "Q_{2} ESE classes:");
    r1->SetNDC();
    r1->SetTextSize(0.04);
    r1->Draw("same");

	TLatex* s1[11];
	s1[0] = new TLatex(0.11,0.2, "1");
	s1[1] = new TLatex(0.15,0.2, "2");
	s1[2] = new TLatex(0.18,0.2, "3");
	s1[3] = new TLatex(0.21,0.2, "4");
	s1[4] = new TLatex(0.23,0.2, "5");
	s1[5] = new TLatex(0.25,0.2, "6");
	s1[6] = new TLatex(0.28,0.2, "7");
	s1[7] = new TLatex(0.31,0.2, "8");
	s1[8] = new TLatex(0.35,0.2, "9");
	s1[9] = new TLatex(0.38,0.2, "10");
	s1[10] = new TLatex(0.50,0.2, "11");
    s1[10]->SetNDC();
	s1[10]->SetTextSize(0.04);
	s1[10]->Draw("same");

	TLatex* s2[11];
	s2[0] = new TLatex(0.70,0.78, "1, 95-100%");
	s2[1] = new TLatex(0.70,0.75, "2, 80-95%");
	s2[2] = new TLatex(0.70,0.72, "3, 60-80%");
	s2[3] = new TLatex(0.70,0.69, "4, 50-60%");
	s2[4] = new TLatex(0.70,0.66, "5, 40-50%");
	s2[5] = new TLatex(0.70,0.63, "6, 30-40%");
	s2[6] = new TLatex(0.70,0.60, "7, 20-30%");
	s2[7] = new TLatex(0.70,0.57, "8, 10-20%");
	s2[8] = new TLatex(0.70,0.54, "9, 5-10%");
	s2[9] = new TLatex(0.70,0.51, "10, 1-5%");
	s2[10] = new TLatex(0.70,0.48, "11, 0-1%");
	s2[10]->SetNDC();
	s2[10]->SetTextSize(0.03);
	s2[10]->Draw("same");

    for(int i = 0; i < 10; i++){
    	l1[i]->SetLineWidth(1.3);
    	l1[i]->SetLineStyle(2);
    	l1[i]->SetLineColor(kRed);
    	l1[i]->Draw("Lsame");

		s1[i]->SetNDC();
		s1[i]->SetTextSize(0.04);
		s1[i]->Draw("same");

		s2[i]->SetNDC();
		s2[i]->SetTextSize(0.03);
		s2[i]->Draw("same");
    }

    TLatex* r11 = new TLatex(0.80,0.92, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

	TLatex* r4 = new TLatex(0.13, 0.84, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
    r4->Draw("same");

    TLatex* r5 = new TLatex(0.13, 0.79, "3.0 < |#Delta#eta| < 5.0");
    r5->SetNDC();
    r5->SetTextSize(23);
    r5->SetTextFont(43);
    r5->SetTextColor(kBlack);
    r5->Draw("same");

    TLatex* lmult = new TLatex(0.13, 0.74, "185 #leq N^{offline}_{trk} < 250");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");

    return;


	TCanvas*c2 = new TCanvas("c2","c2",1,1,700,700);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	//pPb
	//double v2_tracker[] = {0.0618466,0.0629346,0.0650344,0.0670603,0.0686977,0.0705805,0.0728225,0.0759878,0.0794297,0.0830854,0.0868916};
	//double q2_HF[] = {0.0189014,0.0453839,0.0753163,0.0979606,0.113383,0.129792,0.149015,0.17511,0.204744,0.24038,0.303132};
	//PbPb
	double v2_tracker[] = {0.0831979,0.0852045,0.0889806,0.0922329,0.0947941,0.0976101,0.101027,0.105152,0.109602,0.114184,0.118274};
	double q2_HF[] = {0.0195701,0.0475269,0.078339,0.100965,0.116397,0.133288,0.153026,0.17871,0.207737,0.242607,0.302706};
	
	double xwidth[] = {0,0,0,0,0,0,0,0,0,0,0};
	TGraphErrors* gr1 = new TGraphErrors(11, q2_HF, v2_tracker, xwidth, xwidth);

	TH1D* base1 = makeHist("base1","","Q_{2}", "v_{2,tracker}", 1000,0.0,0.4,kBlack);
	base1->GetYaxis()->SetRangeUser(0.07, 0.13);
	base1->GetXaxis()->SetRangeUser(0.0, 0.34);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.6);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.6);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.6);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.6);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);
	base1->Draw();
	

	gr1->SetMarkerStyle(20);
	gr1->SetMarkerSize(1.4);
	gr1->SetMarkerColor(kRed);
	gr1->SetLineColor(kRed);
	gr1->Draw("Psame");


    TLatex* r11 = new TLatex(0.80,0.92, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

	TLatex* r4 = new TLatex(0.17, 0.84, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
    r4->Draw("same");

    TLatex* r5 = new TLatex(0.17, 0.79, "3.0 < |#Delta#eta| < 5.0");
    r5->SetNDC();
    r5->SetTextSize(23);
    r5->SetTextFont(43);
    r5->SetTextColor(kBlack);
    r5->Draw("same");

    TLatex* lmult = new TLatex(0.17, 0.74, "185 #leq N^{offline}_{trk} < 250");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");


}