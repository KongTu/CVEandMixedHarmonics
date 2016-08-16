#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {0,35,60,90,120,150,185,220,260,300};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
int ntrkBinCenter[] = {17.5, 47.5, 75, 105, 135, 167.5, 202.5, 240};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5, 270.4, 310.9};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};

double PbPb_centralityBinCenter[] = {35, 45, 55, 65, 75};
const int Nmults = 5;

double total_systematics_pPb = 0.00015;
double total_systematics_PbPb = 0.00014;

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

void makeTGraph_PbPb_MixedHarmonics(){
	
	TFile* file1[5];
	
	for(int i = 0; i < 5; i++){

		file1[i] = new TFile(Form("../rootfiles/CVEandMixedHarmonics_PbPb_30_100_v2_%d.root", i+1));
	}

	TFile* file = new TFile("~/2015RUN2work/2015Analysis/CMEandCorrelation/ThreePointCorrelator/rootfiles/CME_QvsdEta_pPb_HM_v32_3.root");

	TH1D* QvsdEta[5][48][3][2];

	TH1D* delEta3p[5][3][2];

	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				delEta3p[mult][sign][HF] = (TH1D*) file->Get(Form("ana/delEta3p_%d_%d",sign,HF));
			}
		}
	}

	TH1D* QaQb[5]; TH1D* QaQc[5]; TH1D* QcQb[5];

	for(int mult = 0; mult < Nmults; mult++){

		QaQb[mult] = (TH1D*)file1[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file1[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file1[mult]->Get("ana/c2_cb");
	}

	for(int mult = 0; mult < Nmults; mult++){
		for(int deta = 0; deta < NdEtaBins; deta++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){
			  
				  QvsdEta[mult][deta][sign][HF] = (TH1D*) file1[mult]->Get( Form("ana/c3_real_%d_%d_%d",deta,sign,HF) );
				  
				}
			}
		}
	}

	double v2[23][3];//get corrected v2_3

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

	TH1D* hist1[5][2];
	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF), "", NntrkBins, ntrkBins);
		}
	}

	double threeParticleNtrk[10][3][2];
	double threeParticleNtrkError[10][3][2];
	double totalWeight[10][3][2];

	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				for(int deta = 0; deta < 16; deta++){

					double Q_total_real_dEta = QvsdEta[mult][deta][sign][HF]->GetMean();
					double Q_total_real_dEta_error = QvsdEta[mult][deta][sign][HF]->GetMeanError();
					double deltaEtaWeight = delEta3p[mult][sign][HF]->GetBinContent( deta+1 );

					threeParticleNtrk[mult][sign][HF] += Q_total_real_dEta*deltaEtaWeight;
					double temp = (Q_total_real_dEta_error*Q_total_real_dEta_error)*(deltaEtaWeight*deltaEtaWeight);
					threeParticleNtrkError[mult][sign][HF] += temp;
					totalWeight[mult][sign][HF] += deltaEtaWeight;

				}
			}
		}
	}	

	//pPb:
	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){
			for(int mult = 0; mult < Nmults; mult++){

				//pPb(0,7)
				double value = threeParticleNtrk[mult][sign][HF]/totalWeight[mult][sign][HF];
				value = value/v2[mult][HF];
				hist1[sign][HF]->SetBinContent( mult+1, value);
				
				double error = threeParticleNtrkError[mult][sign][HF]/(totalWeight[mult][sign][HF]*totalWeight[mult][sign][HF]);
				error = sqrt(error)/v2[mult][HF];
				hist1[sign][HF]->SetBinError( mult+1, error);

			}
		}
	}

	TH1D* temp1 = (TH1D*)hist1[0][0]->Clone("temp1");
	TH1D* temp2 = (TH1D*)hist1[0][1]->Clone("temp1");
	TH1D* temp3 = (TH1D*)hist1[1][0]->Clone("temp1");
	TH1D* temp4 = (TH1D*)hist1[1][1]->Clone("temp1");
	
	temp1->Add(temp2, +1);
	temp1->Add(temp3, +1);
	temp1->Add(temp4, +1);

	temp1->Scale(0.25);
	temp1->SetMarkerStyle(24);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp5 = (TH1D*) hist1[2][0]->Clone("temp5");
	TH1D* temp6 = (TH1D*) hist1[2][1]->Clone("temp6");
	
	temp5->Add(temp6, +1);
	
	temp5->Scale(0.5);
	temp5->SetMarkerStyle(25);
	temp5->SetMarkerColor(kBlue);
	temp5->SetLineColor(kBlue);


    double value1[10];
    double value1_error[10];
    double value2[10];
    double value2_error[10];
	
	double PbPb_centralityBinCenter_fill[10];

    for(int mult = 0; mult < Nmults; mult++){

    	value1[mult] = temp1->GetBinContent(mult+1);
    	value1_error[mult] = temp1->GetBinError(mult+1);

    	value2[mult] = temp5->GetBinContent(mult+1);
    	value2_error[mult] = temp5->GetBinError(mult+1);

    	PbPb_centralityBinCenter_fill[mult] = 80 - PbPb_centralityBinCenter[mult];


    }

    TGraphErrors* gr1 = new TGraphErrors(Nmults, PbPb_centralityBinCenter_fill, value1, xbinwidth, value1_error);
    TGraphErrors* gr2 = new TGraphErrors(Nmults, PbPb_centralityBinCenter_fill, value2, xbinwidth, value2_error);

    TFile t1("../dataPoints/PbPb_data_centrality.root","RECREATE");
    gr1->Write();
    gr2->Write();




}