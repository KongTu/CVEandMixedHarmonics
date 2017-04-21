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

void makeTGraph_cme_pPb_v2Eta(){

	double pPb_ntrkBinCenter[6];
	
	TFile* file[10];
	TH1D* QaQb[16]; TH1D* QaQc[16]; TH1D* QcQb[16];
	TH1D* Ntrk[10];
	TH1D* delEta3p[10][3];
	TH1D* delEta2p[10][3];
	TH1D* c2_tracker[10];
	TH1D* cn_eta[48][2];

	for(int mult = 0; mult < 1; mult++){

		file[mult] = new TFile(Form("../rootfiles/CMEandMixedHarmonics_PbPb_Centrality_full_v8.root", mult+1) ); 
		Ntrk[mult] = (TH1D*) file[mult]->Get("ana/Ntrk");
		pPb_ntrkBinCenter[mult] = Ntrk[mult]->GetMean();

		QaQb[mult] = (TH1D*)file[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file[mult]->Get("ana/c2_cb");
		c2_tracker[mult] = (TH1D*) file[mult]->Get("ana/cn_tracker");


		for(int eta = 0; eta < 48; eta++){
			for(int HF = 0; HF < 2; HF++){

				cn_eta[eta][HF] = (TH1D*) file[mult]->Get(Form("ana/cn_eta_%d_%d", eta, HF));
			
			}
		}
	}

	cout << "v2: " << cn_eta[24][0]->GetMean() << endl;


/*
v2 (event plane resolution)
*/

	double v2[10][3];
	double v2_tracker[10];

	for(int mult = 0; mult < 1; mult++){
		
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

		cout << "resolution: " << v2[mult][0] << endl;
	}

	TH1D* v2_hist_abs_Pb = new TH1D("v2_hist_abs_Pb","v2_hist_abs_Pb", Nbins, etabins);
	TH1D* v2_hist_abs_p = new TH1D("v2_hist_abs_p","v2_hist_abs_p", Nbins, etabins);
	
	for(int eta = 0; eta < Nbins; eta++){

		double v2_temp = cn_eta[eta][0]->GetMean();
		double v2_temp_err = cn_eta[eta][0]->GetMeanError();

		double value = v2_temp/v2[0][0];
		double error = v2_temp_err/v2[0][0];

		v2_hist_abs_Pb->SetBinContent(eta+1, value);
		v2_hist_abs_Pb->SetBinError(eta+1, error);

		double v2_temp = cn_eta[eta][1]->GetMean();
		double v2_temp_err = cn_eta[eta][1]->GetMeanError();

		double value = v2_temp/v2[0][1];
		double error = v2_temp_err/v2[0][1];

		v2_hist_abs_p->SetBinContent(eta+1, value);
		v2_hist_abs_p->SetBinError(eta+1, error);

	}

	
	TCanvas* c3 = makeCanvas("c3","c3");
	TH1D* hist2 = makeHist("hist2","hist2","#eta", "v_{2}", 1000,-2.4,2.4,kBlack);
	hist2->GetYaxis()->SetRangeUser(0,0.15);
	hist2->Draw();

	v2_hist_abs_Pb->SetMarkerStyle(20);
	v2_hist_abs_Pb->SetMarkerColor(kRed);
	v2_hist_abs_Pb->Draw("Psame");

	v2_hist_abs_p->SetMarkerStyle(21);
	v2_hist_abs_p->SetMarkerColor(kBlue);
	v2_hist_abs_p->Draw("Psame");

	TLegend *w4 = new TLegend(0.2,0.18,0.8,0.35);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(23);
    w4->SetTextFont(45);
    w4->AddEntry(v2_hist_abs_Pb, "pPb 8 TeV, #phi_{c}(Pb-going)","P");
    w4->AddEntry(v2_hist_abs_p, "pPb 8 TeV, #phi_{c}(p-going)","P");
 	w4->Draw("same");

	TH1D* v2_hist_Pb = new TH1D("v2_hist_Pb","v2_hist_Pb", NdEtaBins, dEtaBins);
	TH1D* v2_hist_p = new TH1D("v2_hist_p","v2_hist_p", NdEtaBins, dEtaBins);
	
	for(int ieta = 0; ieta < Nbins; ieta++){
		for(int jeta = 0; jeta < Nbins; jeta++){

			double deltaEta = fabs(etabins[ieta] - etabins[jeta]);

			for(int deta = 0; deta < NdEtaBins; deta++){
				if( deltaEta > dEtaBins[deta] && deltaEta < dEtaBins[deta+1] ){

					TH1D* temp_sum_Pb = (TH1D*) cn_eta[ieta][0];
					TH1D* temp_sum1_Pb = (TH1D*) cn_eta[jeta][0];

					temp_sum_Pb->Add(temp_sum1_Pb, +1);
					double v2_temp = temp_sum_Pb->GetMean();
					double v2_temp_err = temp_sum_Pb->GetMeanError();

					double value = v2_temp/v2[0][0];
					double error = v2_temp_err/v2[0][0];

					v2_hist_Pb->SetBinContent(deta+1, value);
					v2_hist_Pb->SetBinError(deta+1, error);

					TH1D* temp_sum_p = (TH1D*) cn_eta[ieta][1];
					TH1D* temp_sum1_p = (TH1D*) cn_eta[jeta][1];

					temp_sum_p->Add(temp_sum1_p, +1);
					double v2_temp = temp_sum_p->GetMean();
					double v2_temp_err = temp_sum_p->GetMeanError();

					double value = v2_temp/v2[0][1];
					double error = v2_temp_err/v2[0][1];

					v2_hist_p->SetBinContent(deta+1, value);
					v2_hist_p->SetBinError(deta+1, error);

				}

			}
		}
	}

	TCanvas* c1 = makeCanvas("c1","c1");
	TH1D* hist1 = makeHist("hist1","hist1","|#Delta#eta|", "v_{2}", 1000,0.0,4.8,kBlack);
	hist1->GetYaxis()->SetRangeUser(0,0.15);
	hist1->Draw();

	v2_hist_Pb->SetMarkerStyle(20);
	v2_hist_Pb->SetMarkerColor(kRed);
	v2_hist_Pb->Draw("Psame");

	v2_hist_p->SetMarkerStyle(21);
	v2_hist_p->SetMarkerColor(kBlue);
	v2_hist_p->Draw("Psame");

	TLegend *w4 = new TLegend(0.2,0.18,0.8,0.35);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(23);
    w4->SetTextFont(45);
    w4->AddEntry(v2_hist_Pb, "pPb 8 TeV, #phi_{c}(Pb-going)","P");
    w4->AddEntry(v2_hist_p, "pPb 8 TeV, #phi_{c}(p-going)","P");
 	w4->Draw("same");



}