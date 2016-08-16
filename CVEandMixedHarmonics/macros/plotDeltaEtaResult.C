#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
//rebin option1:
double dEtaReBins[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.4,4.2,4.8};
const int NdEtaReBins = sizeof(dEtaReBins) / sizeof(dEtaReBins[0]) - 1;

double dEtaReBinCenter[] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.2,2.6,3.1,3.8,4.5};

//rebin option2:
double dEtaReBins2[] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.2,2.8,3.8,4.8};
const int NdEtaReBins2 = sizeof(dEtaReBins2) / sizeof(dEtaReBins2[0]) - 1;

double dEtaReBinCenter2[] = {0.15,0.45,0.75,1.05,1.35,1.65,2.0,2.5,3.3,4.3};

double ntrkBins[] = {0,35,60,90,120,150,185,220,260};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
const int Nmults = 2;

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;

double weightedAverage(double a1, double a2, double a3, double eta1, double eta2, double eta3){

	double temp1 = a1*eta1 + a2*eta2 + a3*eta3;
	double temp2 = (a1+a2+a3);

	return temp1/temp2;
}

double weightedAverageError(double a1, double a2, double a3, double etaError1, double etaError2, double etaError3){

	double temp1 = (a1/(a1+a2+a3))*(a1/(a1+a2+a3));
	double temp2 = etaError1*etaError1;
	double temp3 = (a2/(a1+a2+a3))*(a2/(a1+a2+a3));
	double temp4 = etaError2*etaError2;
	double temp5 = (a3/(a1+a2+a3))*(a3/(a1+a2+a3));
	double temp6 = etaError3*etaError3;

	double total = temp1*temp2 + temp3*temp4 + temp5*temp6;

	return sqrt(total);

}

void plotDeltaEtaResult(){

	gStyle->SetErrorX(0);

	TFile* file[16];

	file[0] = new TFile("../rootfiles/CVEandMixedHarmonics_PbPb_30_100_v2_3.root");
	file[1] = new TFile("../rootfiles/CVEandMixedHarmonics_PbPb_30_100_v2_4.root");	
	file[2] = new TFile("/Users/kongkong/2015RUN2work/2015Analysis/CMEandCorrelation/ThreePointCorrelator/rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v17_1.root");

 	string multrange = "185 #leq N^{offline}_{trk} < 220";

	TH1D* QvsdEta[16][48][3][2];

	TH1D* delEta3p[16][3][2];

	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				delEta3p[mult][sign][HF] = (TH1D*) file[2]->Get(Form("ana/delEta3p_%d_%d",sign,HF));
			}
		}
	}

	TH1D* QaQb[16]; TH1D* QaQc[16]; TH1D* QcQb[16];
	TH1D* aveQ3[16][2][2];

	for(int mult = 0; mult < Nmults; mult++){

		QaQb[mult] = (TH1D*)file[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file[mult]->Get("ana/c2_cb");
	}


	for(int mult = 0; mult < Nmults; mult++){
		for(int deta = 0; deta < NdEtaBins; deta++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){
			  
				  QvsdEta[mult][deta][sign][HF] = (TH1D*) file[mult]->Get( Form("ana/c3_real_%d_%d_%d",deta,sign,HF) );
				  
				}
			}
		}
	}

	double v2[16][3];//get corrected v2_3

	for(int mult = 0; mult < Nmults; mult++){
		
		double meanQaQb = QaQb[mult]->GetMean();
		double meanQaQc = QaQc[mult]->GetMean();
		double meanQcQb = QcQb[mult]->GetMean();

		double c2_a = meanQaQb*meanQaQc/meanQcQb;
		double c2_b = meanQaQb*meanQcQb/meanQaQc;
		double c2_ab = meanQaQb;

		v2[mult][0] = sqrt(c2_b );
		v2[mult][1] = sqrt(c2_a );
		v2[mult][2] = sqrt(c2_ab );

	}

	TH1D* hist1[8][3][2];
	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){
				hist1[mult][sign][HF] = new TH1D(Form("hist1_%d_%d_%d",mult,sign,HF),"test", NdEtaReBins2, dEtaReBins2);
			}
		}
	}

	for(int mult = 0; mult < Nmults; mult++){
		for(int deta = 0; deta < NdEtaReBins2; deta++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){

					if(deta < 9){

						double Q_total_real_dEta1 = QvsdEta[mult][3*deta][sign][HF]->GetMean();
						double Q_total_real_dEta_error1 = QvsdEta[mult][3*deta][sign][HF]->GetMeanError();

						double Q_total_real_dEta2 = QvsdEta[mult][3*deta+1][sign][HF]->GetMean();
						double Q_total_real_dEta_error2 = QvsdEta[mult][3*deta+1][sign][HF]->GetMeanError();

						double Q_total_real_dEta3 = QvsdEta[mult][3*deta+2][sign][HF]->GetMean();
						double Q_total_real_dEta_error3 = QvsdEta[mult][3*deta+2][sign][HF]->GetMeanError();

						double weight1 = delEta3p[mult][sign][HF]->GetBinContent( 3*deta+1 );
						double weight2 = delEta3p[mult][sign][HF]->GetBinContent( 3*deta+2 );
						double weight3 = delEta3p[mult][sign][HF]->GetBinContent( 3*deta+3 );
						
						double value = weightedAverage(weight1, weight2, weight3, Q_total_real_dEta1, Q_total_real_dEta2, Q_total_real_dEta3);
						double error = weightedAverageError(weight1, weight2, weight3, Q_total_real_dEta_error1, Q_total_real_dEta_error2, Q_total_real_dEta_error3 );
						
						hist1[mult][sign][HF]->SetBinContent(deta+1, value );
						hist1[mult][sign][HF]->SetBinError(deta+1, error );


					}
					else{

						double Q_total_real_dEta1 = QvsdEta[mult][27][sign][HF]->GetMean();
						double Q_total_real_dEta_error1 = QvsdEta[mult][27][sign][HF]->GetMeanError();
						
						double Q_total_real_dEta2 = QvsdEta[mult][28][sign][HF]->GetMean();
						double Q_total_real_dEta_error2 = QvsdEta[mult][28][sign][HF]->GetMeanError();

						double weight1 = delEta3p[mult][sign][HF]->GetBinContent( 28 );
						double weight2 = delEta3p[mult][sign][HF]->GetBinContent( 29 );

						double value = weightedAverage(weight1, weight2, 0, Q_total_real_dEta1, Q_total_real_dEta2, 0);
						double error = weightedAverageError(weight1, weight2, 0, Q_total_real_dEta_error1, Q_total_real_dEta_error2, 0 );
						
						hist1[mult][sign][HF]->SetBinContent(deta+1, value );
						hist1[mult][sign][HF]->SetBinError(deta+1,  error);
					}


				}
			}
		}
	}

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "", "#Delta#eta", "#LTcos(n_{1}#phi_{#alpha}+n_{2}#phi_{#beta}+n_{3}#phi_{c})#GT/v_{n_{3},c}", 48,0,4.8,kBlack);
	TH1D* base2 = makeHist("base2", "", "#Delta#eta", "#LTcos(n_{1}#phi_{#alpha}+n_{2}#phi_{#beta}+n_{3}#phi_{c})#GT/v_{n_{3},c}", 48,0,4.8,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0012,0.00);
	//base1->GetYaxis()->SetRangeUser(-0.01, 0.03);
	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.0012,0.00);
	//base2->GetYaxis()->SetRangeUser(-0.01, 0.03);

	base2->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetXaxis()->SetTitleOffset(0.95);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetTitle("#LTcos(n_{1}#phi_{#alpha}+n_{2}#phi_{#beta}+n_{3}#phi_{c})#GT/v_{n_{3},c} (oppo - same)");
	
	base3->GetYaxis()->SetRangeUser(-0.0002,0.001);
	
	base3->GetYaxis()->SetTitleOffset(1.2);
	base3->GetXaxis()->SetTitleOffset(1.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.2);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.2);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetNdivisions(5,6,0);
	

	TH1D* temp1_plot = (TH1D*)hist1[0][0][0]->Clone("temp1_plot");
	temp1_plot->Add(hist1[0][1][0], +1);
	temp1_plot->Add(hist1[0][0][1], +1);
	temp1_plot->Add(hist1[0][1][1], +1);
	temp1_plot->Scale(0.25);
	temp1_plot->Scale(1.0/v2[0][2]);
	temp1_plot->SetMarkerStyle(20);
	temp1_plot->SetMarkerSize(1.4);
	temp1_plot->SetMarkerColor(kRed);
	temp1_plot->SetLineColor(kRed);

	TH1D* temp2 = (TH1D*) hist1[0][2][0]->Clone("temp2");
	temp2->Add(hist1[0][2][1], +1);
	temp2->Scale(0.5);
	temp2->SetMarkerStyle(21);
	temp2->Scale(1.0/v2[0][2]);
	temp2->SetMarkerColor(kBlue);
	temp2->SetMarkerSize(1.4);
	temp2->SetLineColor(kBlue);

	TH1D* temp5_plot = (TH1D*)hist1[1][0][0]->Clone("temp5_plot");
	temp5_plot->Add(hist1[1][0][1], +1);
	temp5_plot->Add(hist1[1][1][0], +1);
	temp5_plot->Add(hist1[1][1][1], +1);

	temp5_plot->Scale(0.25);
	temp5_plot->Scale(1.0/v2[1][2]);
	temp5_plot->SetMarkerStyle(20);
	temp5_plot->SetMarkerSize(1.4);
	temp5_plot->SetMarkerColor(kRed);
	temp5_plot->SetLineColor(kRed);

	TH1D* temp9 = (TH1D*) hist1[1][2][0]->Clone("temp9");
	TH1D* temp10 = (TH1D*) hist1[1][2][1]->Clone("temp10");
	temp9->Add(temp10, +1);
	temp9->Scale(0.5);
	temp9->Scale(1.0/v2[1][2]);
	temp9->SetMarkerStyle(21);
	temp9->SetMarkerColor(kBlue);
	temp9->SetMarkerSize(1.4);
	temp9->SetLineColor(kBlue);

//plotting difference between unlike and like sign:

    TCanvas* c3 = new TCanvas("c3","c3", 600, 600);
 	gPad->SetLeftMargin(0.16);
 	gPad->SetBottomMargin(0.13);
 	gPad->SetTopMargin(0.1);
	gPad->SetTicks();
	base3->Draw();

	TH1D* diff1 = (TH1D*) temp2->Clone("diff1");
	diff1->SetMarkerStyle(20);
	diff1->SetMarkerSize(1.3);
	diff1->SetMarkerColor(kRed);
	diff1->SetLineColor(kRed);
	diff1->Add(temp1_plot, -1);
	
    TLatex* r33 = new TLatex(0.25, 0.82, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r33->SetNDC();
    r33->SetTextSize(23);
    r33->SetTextFont(43);
    r33->SetTextColor(kBlack);
    //r33->Draw("same");

    TLatex* lmult = new TLatex(0.20, 0.82, multrange.c_str());
    lmult->SetNDC();
    lmult->SetTextSize(26);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    //lmult->Draw("same");

	TH1D* diff3 = (TH1D*) temp9->Clone("diff3");
	diff3->Add(temp5_plot, -1);
	diff3->SetMarkerStyle(28);
	diff3->SetMarkerSize(1.4);
	diff3->SetMarkerColor(kBlack);
	diff3->SetLineColor(kBlack);

	diff1->Draw("Psame");
	diff3->Draw("Psame");

   	TLatex* r11 = new TLatex(0.62,0.91, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.71,0.91, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
    r22->Draw("same");

	TLegend *w3 = new TLegend(0.40,0.5,0.6,0.7);
    w3->SetLineColor(kWhite);
    w3->SetFillColor(0);
    w3->SetTextSize(20);
    w3->SetTextFont(43);
    w3->AddEntry(diff1, "PbPb 5.02 TeV, 50-60%");
    w3->AddEntry(diff3, "PbPb 5.02 TeV, 60-70%");
    w3->Draw("same");

	TH1D* base5 = (TH1D*) base1->Clone("base3");
	base5->GetYaxis()->SetTitleOffset(1.3);
	base5->GetXaxis()->SetTitleOffset(0.95);
	base5->GetYaxis()->SetTitleSize(base5->GetYaxis()->GetTitleSize()*1.3);
	base5->GetXaxis()->SetTitleSize(base5->GetXaxis()->GetTitleSize()*1.3);
	base5->GetYaxis()->SetLabelSize(base5->GetYaxis()->GetLabelSize()*1.3);
	base5->GetXaxis()->SetLabelSize(base5->GetXaxis()->GetLabelSize()*1.3);
	base5->GetXaxis()->SetNdivisions(5,6,0);
	base5->GetYaxis()->SetRangeUser(-0.003,0.0016);

	TH1D* base6 = (TH1D*) base5->Clone("base6");
	base6->GetYaxis()->SetRangeUser(-0.01, 0.04);
	base6->GetXaxis()->SetRangeUser(-0.2, 4.8);

    TLatex* r41 = new TLatex(0.24, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r41->SetNDC();
    r41->SetTextSize(23);
    r41->SetTextFont(43);
    r41->SetTextColor(kBlack);

    TLatex* r42 = new TLatex(0.05, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.55,0.95, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.052);
    
    TLatex* r44 = new TLatex(0.69,0.95, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(24);
    r44->SetTextFont(53);
    
    TLatex* r45 = new TLatex(0.05, 0.80, "60-70%");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

    TLatex* r46 = new TLatex(0.24, 0.80, "50-60%");
    r46->SetNDC();
    r46->SetTextSize(23);
    r46->SetTextFont(43);
    r46->SetTextColor(kBlack);

    TLegend *w4 = new TLegend(0.5,0.2,0.95,0.35);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(23);
    w4->SetTextFont(45);
    w4->SetNColumns(2);
    w4->AddEntry(temp1_plot, "  ");
    w4->AddEntry(temp2, " PbPb");
    

	TLatex* latex1 = new TLatex(0.48, 0.37, "same");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    TLatex* latex2 = new TLatex(0.59, 0.37, "oppo");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    
    TCanvas* c4 = new TCanvas("c4","c4",1000,600);
	c4->Divide(2,1,0,0);
	c4->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base5->Draw();

    temp1_plot->Draw("Psame");
	temp2->Draw("Psame");
	r41->Draw("same");
    //w4->Draw("same");
    //latex1->Draw("same");
    //latex2->Draw("same");
    r46->Draw("same");


    c4->cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.06);
	gPad->SetTicks();
	base5->Draw();
	r45->Draw("same");

    TLegend *w40 = new TLegend(0.5,0.22,0.73,0.4);
    w40->SetLineColor(kWhite);
    w40->SetFillColor(0);
    w40->SetTextSize(23);
    w40->SetTextFont(45);
    w40->SetNColumns(2);
    w40->SetTextAlign();
    w40->AddEntry(temp5_plot, " ");
    w40->AddEntry(temp9, " ");

    TLatex* r434 = new TLatex(0.70,0.29, "PbPb");
    r434->SetNDC();
    r434->SetTextSize(21);
    r434->SetTextFont(43);

    temp5_plot->Draw("Psame");
	temp9->Draw("Psame");
	r42->Draw("same");
	r43->Draw("same");
	r44->Draw("same");
	w40->Draw("same");
	latex1->Draw("same");
    latex2->Draw("same");
    r434->Draw("same");



	// TFile t1("../dataPoints/diff.root","RECREATE");
	// diff1->Write();
	// diff2->Write();
	// diff3->Write();
	// temp1->Write();
	// temp2->Write();
	// temp3->Write();
	// temp4->Write();
	// temp5->Write();
	// temp9->Write();

	//c2->Print("../results/deltaEtaResults.pdf");

	return;

}