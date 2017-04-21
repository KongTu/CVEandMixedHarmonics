#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {0,35,60,90,120,150,185,220,260,300,350,400};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
int ntrkBinCenter[] = {17.5, 47.5, 75, 105, 135, 167.5, 202.5, 240};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

const int Nmults = 5;

double total_systematics_pPb = 0.00015;
double total_systematics_PbPb = 0.00014;


void makeTGraph_cme_pPb_Q2(){

	double pPb_v2_bincenter[11];
	

	TFile* file[11];

	for(int mult = 0; mult < 11; mult++){

		file[mult] = new TFile(Form("../rootfiles/CMEandMixedHarmonics_pPb_HM300_400_q2_v2_%d.root", mult+1));
	}

	TH1D* delEta3p[11][3];
	TH1D* delEta2p[11][3];
	for(int mult = 0; mult < 11; mult++){
		for(int sign = 0; sign < 3; sign++){
			delEta2p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delEta2p_%d",sign));
			delEta3p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delEta3p_%d",sign));
			
		}
	}

	TH1D* QvsdEta[30][48][3][2];
	TH1D* PvsdEta[30][48][3];

	TH1D* QaQb[30]; TH1D* QaQc[30]; TH1D* QcQb[30];
	TH1D* vn_tracker[30];
	TH1D* q2_mag[30];

	for(int mult = 0; mult < 11; mult++){

		QaQb[mult] = (TH1D*)file[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file[mult]->Get("ana/c2_cb");
		vn_tracker[mult] = (TH1D*)file[mult]->Get("ana/cn_tracker");
		q2_mag[mult] = (TH1D*) file[mult]->Get("ana/q2_mag");

	}

	for(int mult = 0; mult < 11; mult++){
		for(int deta = 0; deta < NdEtaBins; deta++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){
			  
				  QvsdEta[mult][deta][sign][HF] = (TH1D*) file[mult]->Get( Form("ana/c3_real_%d_%d_%d",deta,sign,HF) );
				  
				}
			}
		}
	}

	for(int mult = 0; mult < 11; mult++){
		for(int deta = 0; deta < NdEtaBins; deta++){
			for(int sign = 0; sign < 3; sign++){
			  
			  PvsdEta[mult][deta][sign] = (TH1D*) file[mult]->Get( Form("ana/c2_real_%d_%d",deta,sign) );
				  
			}
		}
	}

	double v2[23][3];//get corrected v2_3
	double v2_tracker[23];


	cout << "v2: {";
	for(int mult = 0; mult < 11; mult++){
		
		double meanQaQb = QaQb[mult]->GetMean();
		double meanQaQc = QaQc[mult]->GetMean();
		double meanQcQb = QcQb[mult]->GetMean();

		double c2_a = meanQaQb*meanQaQc/meanQcQb;
		double c2_b = meanQaQb*meanQcQb/meanQaQc;
		double c2_ab = meanQaQb;
		double q2 = q2_mag[mult]->GetMean();

		v2[mult][0] = sqrt(c2_b );
		v2[mult][1] = sqrt(c2_a );
		v2[mult][2] = sqrt(c2_ab );

		double c2 = vn_tracker[mult]->GetMean();
		v2_tracker[mult] =  sqrt( c2 );
		pPb_v2_bincenter[mult] = v2_tracker[mult];
		cout << pPb_v2_bincenter[mult] << ",";
	}
	cout << "}" << endl;

	TH1D* hist1[3][2];
	TH1D* total2[3];

	for(int sign = 0; sign < 3; sign++){
		total2[sign] = new TH1D(Form("total2_%d",sign), "", NntrkBins, ntrkBins);

		for(int HF = 0; HF < 2; HF++){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF), "", NntrkBins, ntrkBins);
		}
	}

	double threeParticleNtrk[11][3][2];
	double threeParticleNtrkError[11][3][2];
	double totalWeight[11][3][2];

	for(int mult = 0; mult < 11; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				threeParticleNtrk[mult][sign][HF] = 0.0;
				threeParticleNtrkError[mult][sign][HF] = 0.0;
				totalWeight[mult][sign][HF] = 0.0;

			}
		}
	}

	for(int mult = 0; mult < 11; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				for(int deta = 0; deta < 16; deta++){

					double Q_total_real_dEta = QvsdEta[mult][deta][sign][HF]->GetMean();
					double Q_total_real_dEta_error = QvsdEta[mult][deta][sign][HF]->GetMeanError();
					double deltaEtaWeight = delEta3p[mult][sign]->GetBinContent( deta+1 );

					threeParticleNtrk[mult][sign][HF] += Q_total_real_dEta*deltaEtaWeight;
					threeParticleNtrkError[mult][sign][HF] += (Q_total_real_dEta_error*Q_total_real_dEta_error)*(deltaEtaWeight*deltaEtaWeight);
					totalWeight[mult][sign][HF] += deltaEtaWeight;

				}
			}
		}
	}
	
	//pPb:
	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){
			for(int mult = 0; mult < 11; mult++){

				//pPb(0,7)
				double value = threeParticleNtrk[mult][sign][HF]/totalWeight[mult][sign][HF];
				value = value/v2[mult][HF];
				hist1[sign][HF]->SetBinContent( mult+1, value);
				double error = threeParticleNtrkError[mult][sign][HF]/(totalWeight[mult][sign][HF]*totalWeight[mult][sign][HF]);
				error = sqrt(error)/v2[mult][HF];

				//error = sqrt(error);
				hist1[sign][HF]->SetBinError( mult+1, error);

			}
		}
	}

	double twoParticleNtrk[11][3];
	double twoParticleNtrkError[11][3];
	double total2pWeight[11][3];

	for(int mult = 0; mult < 11; mult++){
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
		for(int mult = 0; mult < 11; mult++){

			//pPb(0,7)
			double value = twoParticleNtrk[mult][sign]/total2pWeight[mult][sign];
			total2[sign]->SetBinContent( mult+1, value);
			
			double error = twoParticleNtrkError[mult][sign]/(total2pWeight[mult][sign]*total2pWeight[mult][sign]);
			error = sqrt(error);
			total2[sign]->SetBinError( mult+1, error);

		}
	}

	TH1D* temp1 = (TH1D*)hist1[0][0]->Clone("temp1");
	temp1->Add(hist1[1][0], +1);
	temp1->Scale(0.5);
	temp1->SetMarkerStyle(24);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp2 = (TH1D*) hist1[2][0]->Clone("temp2");
	temp2->SetMarkerStyle(25);
	temp2->SetMarkerColor(kBlue);
	temp2->SetLineColor(kBlue);

	TH1D* temp3 = (TH1D*)hist1[0][1]->Clone("temp3");
	temp3->Add(hist1[1][1], +1);
	temp3->Scale(0.5);
	temp3->SetMarkerStyle(20);
	temp3->SetMarkerColor(kRed);
	temp3->SetLineColor(kRed);

	TH1D* temp4 = (TH1D*) hist1[2][1]->Clone("temp4");
	temp4->SetMarkerStyle(21);
	temp4->SetMarkerColor(kBlue);
	temp4->SetLineColor(kBlue);

	TH1D* temp5 = (TH1D*)total2[0]->Clone("temp5");
	temp5->Add(total2[1], +1);
	temp5->Scale(0.5);
	temp5->SetMarkerStyle(20);
	temp5->SetMarkerColor(kRed);
	temp5->SetLineColor(kRed);

	TH1D* temp6 = (TH1D*) total2[2]->Clone("temp6");
	temp6->SetMarkerStyle(21);
	temp6->SetMarkerColor(kBlue);
	temp6->SetLineColor(kBlue);

    double value1[11];
    double value1_error[11];
    double value2[11];
    double value2_error[11];
    double value3[11];
    double value3_error[11];
    double value4[11];
    double value4_error[11];
    double value5[11];
    double value5_error[11];
    double value6[11];
    double value6_error[11];

    for(int mult = 0; mult < 11; mult++){

    	value1[mult] = temp1->GetBinContent(mult+1);
    	value1_error[mult] = temp1->GetBinError(mult+1);

    	value2[mult] = temp2->GetBinContent(mult+1);
    	value2_error[mult] = temp2->GetBinError(mult+1);

    	value3[mult] = temp3->GetBinContent(mult+1);
    	value3_error[mult] = temp3->GetBinError(mult+1);

    	value4[mult] = temp4->GetBinContent(mult+1);
    	value4_error[mult] = temp4->GetBinError(mult+1);

    	value5[mult] = temp5->GetBinContent(mult+1);
    	value5_error[mult] = temp5->GetBinError(mult+1);

    	value6[mult] = temp6->GetBinContent(mult+1);
    	value6_error[mult] = temp6->GetBinError(mult+1);
    }

    TGraphErrors* gr1 = new TGraphErrors(11, pPb_v2_bincenter, value1, xbinwidth, value1_error);
    TGraphErrors* gr2 = new TGraphErrors(11, pPb_v2_bincenter, value2, xbinwidth, value2_error);
    TGraphErrors* gr3 = new TGraphErrors(11, pPb_v2_bincenter, value3, xbinwidth, value3_error);
    TGraphErrors* gr4 = new TGraphErrors(11, pPb_v2_bincenter, value4, xbinwidth, value4_error);
   	TGraphErrors* gr5 = new TGraphErrors(11, pPb_v2_bincenter, value5, xbinwidth, value5_error);
    TGraphErrors* gr6 = new TGraphErrors(11, pPb_v2_bincenter, value6, xbinwidth, value6_error);

    gr5->Draw();

    gr6->SetLineColor(kRed);
    gr6->Draw("same");

    TFile t1("../dataPoints/pPb_8TeV_cme_q2_5.root","RECREATE");
    gr1->Write();
    gr2->Write();
    gr3->Write();
    gr4->Write();
 	gr5->Write();
    gr6->Write();


}
