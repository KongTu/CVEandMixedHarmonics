#include "RiceStyle.h"

using namespace std;

double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {120,150,185,250,300,350,400};

//rebin option2:
double dEtaReBins2[] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.2,2.8,3.8,4.8};
const int NdEtaReBins2 = sizeof(dEtaReBins2) / sizeof(dEtaReBins2[0]) - 1;
double dEtaReBinCenter2[] = {0.15,0.45,0.75,1.05,1.35,1.65,2.0,2.5,3.3,4.3};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0};

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

//both Ntrk and deta are saved from this output

void makeTGraph_cme_pPb_Ntrk(){

	double pPb_ntrkBinCenter[6];

	
	TFile* file[10];
	TH1D* QaQb[16]; TH1D* QaQc[16]; TH1D* QcQb[16];
	TH1D* Ntrk[10];
	TH1D* delEta3p[10][3];
	TH1D* delEta2p[10][3];
	TH1D* c2_tracker[10];

	for(int mult = 0; mult < 6; mult++){

		file[mult] = new TFile(Form("../rootfiles/CMEandMixedHarmonics_pPb_full_HM_v2_%d.root", mult+1) ); 
		Ntrk[mult] = (TH1D*) file[mult]->Get("ana/Ntrk");
		pPb_ntrkBinCenter[mult] = Ntrk[mult]->GetMean();

		QaQb[mult] = (TH1D*)file[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file[mult]->Get("ana/c2_cb");
		c2_tracker[mult] = (TH1D*) file[mult]->Get("ana/cn_tracker");

		for(int sign = 0; sign < 3; sign++){
			delEta2p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delEta2p_%d",sign));
			delEta3p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delEta3p_%d",sign));
		}
	}

	TH1D* QvsdEta[30][48][3][2];
	TH1D* PvsdEta[30][48][3];

	for(int mult = 0; mult < 6; mult++){
		for(int deta = 0; deta < NdEtaBins; deta++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){
			  
				  QvsdEta[mult][deta][sign][HF] = (TH1D*) file[mult]->Get( Form("ana/c3_real_%d_%d_%d",deta,sign,HF) );
				  
				}
			}
		}
	}

	for(int mult = 0; mult < 6; mult++){
		for(int deta = 0; deta < NdEtaBins; deta++){
			for(int sign = 0; sign < 3; sign++){
			  
			  PvsdEta[mult][deta][sign] = (TH1D*) file[mult]->Get( Form("ana/c2_real_%d_%d",deta,sign) );
				  
			}
		}
	}

/*
v2 (event plane resolution)
*/

	double v2[10][3];
	double v2_tracker[10];

	for(int mult = 0; mult < 6; mult++){
		
		double meanQaQb = QaQb[mult]->GetMean();
		double meanQaQc = QaQc[mult]->GetMean();
		double meanQcQb = QcQb[mult]->GetMean();

		double c2_a = meanQaQb*meanQaQc/meanQcQb;
		double c2_b = meanQaQb*meanQcQb/meanQaQc;
		double c2_ab = meanQaQb;

		v2[mult][0] = sqrt(c2_b );
		v2[mult][1] = sqrt(c2_a  );
		v2[mult][2] = sqrt(c2_ab );

		double c2 = c2_tracker[mult]->GetMean();
		v2_tracker[mult] =  sqrt( c2 );
		cout << "v2_tracker: " << v2_tracker[mult] << endl;
	}

/*
integrated results |deta| < 1.6
*/

	TH1D* total1[3][2];
	TH1D* total2[3];

	for(int sign = 0; sign < 3; sign++){
		total2[sign] = new TH1D(Form("total2_%d",sign), "", 6, ntrkBins);

		for(int HF = 0; HF < 2; HF++){
			total1[sign][HF] = new TH1D(Form("total1_%d_%d",sign,HF), "", 6, ntrkBins);
		}
	}

	double threeParticleNtrk[10][3][2];
	double threeParticleNtrkError[10][3][2];
	double totalWeight[10][3][2];

	for(int mult = 0; mult < 6; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				for(int deta = 0; deta < 16; deta++){

					double Q_total_real_dEta = QvsdEta[mult][deta][sign][HF]->GetMean();
					double Q_total_real_dEta_error = QvsdEta[mult][deta][sign][HF]->GetMeanError();
					double deltaEtaWeight = delEta3p[mult][sign]->GetBinContent( deta+1 );

					threeParticleNtrk[mult][sign][HF] += Q_total_real_dEta*deltaEtaWeight;
					double temp = (Q_total_real_dEta_error*Q_total_real_dEta_error)*(deltaEtaWeight*deltaEtaWeight);
					threeParticleNtrkError[mult][sign][HF] += temp;
					totalWeight[mult][sign][HF] += deltaEtaWeight;

				}
			}
		}
	}	

	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){
			for(int mult = 0; mult < 6; mult++){

				//pPb(0,7)
				double value = threeParticleNtrk[mult][sign][HF]/totalWeight[mult][sign][HF];
				value = value/v2[mult][HF];
				total1[sign][HF]->SetBinContent( mult+1, value);
				
				double error = threeParticleNtrkError[mult][sign][HF]/(totalWeight[mult][sign][HF]*totalWeight[mult][sign][HF]);
				error = sqrt(error)/v2[mult][HF];
				total1[sign][HF]->SetBinError( mult+1, error);

			}
		}
	}

	double twoParticleNtrk[10][3];
	double twoParticleNtrkError[10][3];
	double total2pWeight[10][3];

	for(int mult = 0; mult < 6; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int deta = 0; deta < 16; deta++){

				double P_total_real_dEta = PvsdEta[mult][deta][sign]->GetMean();
				double P_total_real_dEta_error = PvsdEta[mult][deta][sign]->GetMeanError();
				double deltaEtaWeight = delEta2p[mult][sign]->GetBinContent( deta+1 );

				twoParticleNtrk[mult][sign] += P_total_real_dEta*deltaEtaWeight;
				double temp = (P_total_real_dEta_error*P_total_real_dEta_error)*(deltaEtaWeight*deltaEtaWeight);
				twoParticleNtrkError[mult][sign] += temp;
				total2pWeight[mult][sign] += deltaEtaWeight;

			}
		}
	}	

	for(int sign = 0; sign < 3; sign++){
		for(int mult = 0; mult < 6; mult++){

			//pPb(0,7)
			double value = twoParticleNtrk[mult][sign]/total2pWeight[mult][sign];
			total2[sign]->SetBinContent( mult+1, value);
			
			double error = twoParticleNtrkError[mult][sign]/(total2pWeight[mult][sign]*total2pWeight[mult][sign]);
			error = sqrt(error);
			total2[sign]->SetBinError( mult+1, error);

		}
	}

	TH1D* temp11 = (TH1D*)total1[0][0]->Clone("temp11");
	temp11->Add(total1[1][0], +1);
	temp11->Scale(0.5);
	temp11->SetMarkerStyle(24);
	temp11->SetMarkerColor(kRed);
	temp11->SetLineColor(kRed);

	TH1D* temp22 = (TH1D*) total1[2][0]->Clone("temp22");
	temp22->SetMarkerStyle(25);
	temp22->SetMarkerColor(kBlue);
	temp22->SetLineColor(kBlue);

	TH1D* temp33 = (TH1D*)total1[0][1]->Clone("temp33");
	temp33->Add(total1[1][1], +1);
	temp33->Scale(0.5);
	temp33->SetMarkerStyle(20);
	temp33->SetMarkerColor(kRed);
	temp33->SetLineColor(kRed);

	TH1D* temp44 = (TH1D*) total1[2][1]->Clone("temp44");
	temp44->SetMarkerStyle(21);
	temp44->SetMarkerColor(kBlue);
	temp44->SetLineColor(kBlue);

	TH1D* temp55 = (TH1D*)total2[0]->Clone("temp55");
	temp55->Add(total2[1], +1);
	temp55->Scale(0.5);
	temp55->SetMarkerStyle(20);
	temp55->SetMarkerColor(kRed);
	temp55->SetLineColor(kRed);

	TH1D* temp66 = (TH1D*) total2[2]->Clone("temp66");
	temp66->SetMarkerStyle(21);
	temp66->SetMarkerColor(kBlue);
	temp66->SetLineColor(kBlue);

    double value1[6];
    double value1_error[6];
    double value2[6];
    double value2_error[6];
    double value3[6];
    double value3_error[6];
    double value4[6];
    double value4_error[6];

    double value5[6];
    double value5_error[6];
    double value6[6];
    double value6_error[6];

    double value7[6];
    double value7_error[6];
    double value8[6];
    double value8_error[6];

    for(int mult = 0; mult < 6; mult++){

    	value1[mult] = temp11->GetBinContent(mult+1);
    	value1_error[mult] = temp11->GetBinError(mult+1);

    	value2[mult] = temp22->GetBinContent(mult+1);
    	value2_error[mult] = temp22->GetBinError(mult+1);

    	value3[mult] = temp33->GetBinContent(mult+1);
    	value3_error[mult] = temp33->GetBinError(mult+1);

    	value4[mult] = temp44->GetBinContent(mult+1);
    	value4_error[mult] = temp44->GetBinError(mult+1);

    	value5[mult] = v2_tracker[mult];//0
    	value5_error[mult] = QcQb[mult]->GetMeanError();

    	value6[mult] = v2_tracker[mult];//1
    	value6_error[mult] = QaQc[mult]->GetMeanError();

    	value7[mult] = temp55->GetBinContent(mult+1);
    	value7_error[mult] = temp55->GetBinError(mult+1);

    	value8[mult] = temp66->GetBinContent(mult+1);
    	value8_error[mult] = temp66->GetBinError(mult+1);
    }


    TGraphErrors* gr1 = new TGraphErrors(6, pPb_ntrkBinCenter, value1, xbinwidth, value1_error);
    TGraphErrors* gr2 = new TGraphErrors(6, pPb_ntrkBinCenter, value2, xbinwidth, value2_error);
    TGraphErrors* gr3 = new TGraphErrors(6, pPb_ntrkBinCenter, value3, xbinwidth, value3_error);
    TGraphErrors* gr4 = new TGraphErrors(6, pPb_ntrkBinCenter, value4, xbinwidth, value4_error);
    TGraphErrors* gr5 = new TGraphErrors(6, pPb_ntrkBinCenter, value5, xbinwidth, value5_error);
    TGraphErrors* gr6 = new TGraphErrors(6, pPb_ntrkBinCenter, value6, xbinwidth, value6_error);
	TGraphErrors* gr7 = new TGraphErrors(6, pPb_ntrkBinCenter, value7, xbinwidth, value7_error);
    TGraphErrors* gr8 = new TGraphErrors(6, pPb_ntrkBinCenter, value8, xbinwidth, value8_error);


	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	gPad->SetLogx();

	TH1D* base3 = makeHist("base3", "", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 10000,0,10000,kBlack);

	base3->GetYaxis()->SetRangeUser(-0.0009, 0.0006);
	base3->GetXaxis()->SetRangeUser(80, 1700);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.3);
	base3->GetXaxis()->SetTitleOffset(0.95);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.3);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.4);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.4);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.4);

	base3->Draw();

	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(kRed);
	gr1->Draw("Psame");

	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(kBlue);
	gr2->Draw("Psame");

	TCanvas* c4 = new TCanvas("c4","c4",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	gPad->SetLogx();

	TH1D* base4 = makeHist("base4", "", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}-#phi_{#beta})#GT", 10000,0,10000,kBlack);

	base4->GetYaxis()->SetRangeUser(-0.009, 0.02);
	base4->GetXaxis()->SetRangeUser(80, 1700);
	base4->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base4,1.1,1.25);

	base4->GetYaxis()->SetTitleOffset(1.3);
	base4->GetXaxis()->SetTitleOffset(0.95);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.3);
	base4->GetXaxis()->SetTitleSize(base4->GetXaxis()->GetTitleSize()*1.4);
	base4->GetYaxis()->SetLabelSize(base4->GetYaxis()->GetLabelSize()*1.4);
	base4->GetXaxis()->SetLabelSize(base4->GetXaxis()->GetLabelSize()*1.4);

	base4->Draw();

	gr7->SetMarkerStyle(20);
	gr7->SetMarkerColor(kRed);
	gr7->Draw("Psame");

	gr8->SetMarkerStyle(21);
	gr8->SetMarkerColor(kBlue);
	gr8->Draw("Psame");

	// TFile f1("../dataPoints/CME_Ntrk_112_pPb_8TeV.root", "RECREATE");

	// gr1->Write();
	// gr2->Write();
	// gr3->Write();
	// gr4->Write();
	// gr5->Write();
	// gr6->Write();
	// gr7->Write();
	// gr8->Write();



}