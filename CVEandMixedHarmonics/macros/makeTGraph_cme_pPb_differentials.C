#include "RiceStyle.h"

using namespace std;

double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;


double dPtBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
                    1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,
                    2.4,2.6,2.8,3.0,3.4,4.0,5.0,6.0};
const int NdPtBins = sizeof(dPtBins) / sizeof(dPtBins[0]) - 1;


double dPtBinsRebins1[] = {0.0,0.2,0.4,0.6,0.8,1.0,
                    1.2,1.4,1.6,1.8,2.0,
                    2.4,2.8};
const int NdPtBinsRebins1 = sizeof(dPtBinsRebins1) / sizeof(dPtBinsRebins1[0]) - 1;

double dPtBinsRebins2[] = {0.3,0.5,0.7,0.9,
                    1.1,1.3,1.5,1.7,1.9,2.2,
                    2.6,3.0};
const int NdPtBinsRebins2 = sizeof(dPtBinsRebins2) / sizeof(dPtBinsRebins2[0]) - 1;


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

void makeTGraph_cme_pPb_differentials(){

	TFile* file[106];
	TH1D* QaQb[16]; TH1D* QaQc[16]; TH1D* QcQb[16];
	TH1D* Ntrk[10];
	TH1D* delEta3p[16][3];
	TH1D* delEta2p[16][3];
	TH1D* delPt3p[16][3];
	TH1D* delPt2p[16][3];
	TH1D* ptAve3p[16][3];
	TH1D* ptAve2p[16][3];

	double pPb_ntrkBinCenter[6];


	for(int mult = 0; mult < 1; mult++){

		file[mult] = new TFile(Form("../rootfiles/CMEandMixedHarmonics_AMPT_pPb_v1.root", mult+1) ); 
		Ntrk[mult] = (TH1D*) file[mult]->Get("ana/Ntrk");
		pPb_ntrkBinCenter[mult] = Ntrk[mult]->GetMean();

		QaQb[mult] = (TH1D*)file[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file[mult]->Get("ana/c2_cb");

		for(int sign = 0; sign < 3; sign++){
			
			delEta2p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delEta2p_%d",sign));
			delEta3p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delEta3p_%d",sign));
			
			delPt2p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delPt2p_%d",sign));
			delPt3p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delPt3p_%d",sign));
			
			ptAve2p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/ptAve2p_%d",sign));
			ptAve3p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/ptAve3p_%d",sign));
			
		}
	}

	TH1D* QvsdEta[30][48][3][2];
	TH1D* QvsdPt[30][48][3][2];
	TH1D* QvsPtAve[30][48][3][2];

	TH1D* PvsdEta[30][48][3];
	TH1D* PvsdPt[30][48][3];
	TH1D* PvsPtAve[30][48][3];

	for(int mult = 0; mult < 1; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){
				for(int deta = 0; deta < NdEtaBins; deta++){
			  
					QvsdEta[mult][deta][sign][HF] = (TH1D*) file[mult]->Get( Form("ana/c3_real_%d_%d_%d",deta,sign,HF) ); 
				}
				for(int dpt = 0; dpt < NdPtBins; dpt++){

					//if( dpt >= 23 ) continue;//where the nonzero entries are
					QvsdPt[mult][dpt][sign][HF] = (TH1D*) file[mult]->Get( Form("ana/c3_dpT_real_%d_%d_%d",dpt,sign,HF) ); 
				}
				for(int dpt = 0; dpt < NdPtBins; dpt++){

					//if (dpt < 3 || dpt > 24 ) continue;
					QvsPtAve[mult][dpt][sign][HF] = (TH1D*) file[mult]->Get( Form("ana/c3_pTave_real_%d_%d_%d",dpt,sign,HF) ); 
				}

			}
		}
	}

	for(int mult = 0; mult < 1; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int deta = 0; deta < NdEtaBins; deta++){
			  
				PvsdEta[mult][deta][sign] = (TH1D*) file[mult]->Get( Form("ana/c2_real_%d_%d",deta,sign) );
				  
			}
			for(int dpt = 0; dpt < NdPtBins; dpt++){

				//if( dpt >= 23 ) continue;//where the nonzero entries are
				PvsdPt[mult][dpt][sign] = (TH1D*) file[mult]->Get( Form("ana/c2_dpT_real_%d_%d",dpt,sign) ); 
			}
			for(int dpt = 0; dpt < NdPtBins; dpt++){

				//if (dpt < 3 || dpt > 24 ) continue;
				PvsPtAve[mult][dpt][sign] = (TH1D*) file[mult]->Get( Form("ana/c2_pTave_real_%d_%d",dpt,sign) ); 
			}
		}
	}

/*
v2 (event plane resolution)
*/

	double v2[10][3];

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
	}

	TH1D* hist1[8][3][2];
	TH1D* hist2[8][3];

	TH1D* hist1_dpT[8][3][2];
	TH1D* hist2_dpT[8][3];

	TH1D* hist1_pTave[8][3][2];
	TH1D* hist2_pTave[8][3];

	for(int mult = 0; mult < 1; mult++){
		for(int sign = 0; sign < 3; sign++){

			hist2[mult][sign] = new TH1D(Form("hist2_%d_%d",mult,sign),"test1", NdEtaReBins2, dEtaReBins2);
			hist2_dpT[mult][sign] = new TH1D(Form("hist2_dpT_%d_%d",mult,sign),"test1", NdPtBinsRebins1, dPtBinsRebins1);
			hist2_pTave[mult][sign] = new TH1D(Form("hist2_pTave_%d_%d",mult,sign),"test1", NdPtBinsRebins2, dPtBinsRebins2);

			for(int HF = 0; HF < 2; HF++){
				hist1[mult][sign][HF] = new TH1D(Form("hist1_%d_%d_%d",mult,sign,HF),"test", NdEtaReBins2, dEtaReBins2);
				hist1_dpT[mult][sign][HF] = new TH1D(Form("hist1_dpT_%d_%d_%d",mult,sign,HF),"test", NdPtBinsRebins1, dPtBinsRebins1);
				hist1_pTave[mult][sign][HF] = new TH1D(Form("hist1_pTave_%d_%d_%d",mult,sign,HF),"test", NdPtBinsRebins2, dPtBinsRebins2);

			}
		}
	}

	for(int mult = 0; mult < 1; mult++){
		for(int deta = 0; deta < NdEtaReBins2; deta++){
			for(int sign = 0; sign < 3; sign++){

				if(deta < 9){

					double P_total_real_dEta1 = PvsdEta[mult][3*deta][sign]->GetMean();
					double P_total_real_dEta_error1 = PvsdEta[mult][3*deta][sign]->GetMeanError();

					double P_total_real_dEta2 = PvsdEta[mult][3*deta+1][sign]->GetMean();
					double P_total_real_dEta_error2 = PvsdEta[mult][3*deta+1][sign]->GetMeanError();

					double P_total_real_dEta3 = PvsdEta[mult][3*deta+2][sign]->GetMean();
					double P_total_real_dEta_error3 = PvsdEta[mult][3*deta+2][sign]->GetMeanError();

					double weight1 = delEta3p[mult][sign]->GetBinContent( 3*deta+1 );
					double weight2 = delEta3p[mult][sign]->GetBinContent( 3*deta+2 );
					double weight3 = delEta3p[mult][sign]->GetBinContent( 3*deta+3 );
					
					double value = weightedAverage(weight1, weight2, weight3, P_total_real_dEta1, P_total_real_dEta2, P_total_real_dEta3);
					double error = weightedAverageError(weight1, weight2, weight3, P_total_real_dEta_error1, P_total_real_dEta_error2, P_total_real_dEta_error3 );
					
					hist2[mult][sign]->SetBinContent(deta+1, value );
					hist2[mult][sign]->SetBinError(deta+1, error );

				}
				else{

					double P_total_real_dEta1 = PvsdEta[mult][27][sign]->GetMean();
					double P_total_real_dEta_error1 = PvsdEta[mult][27][sign]->GetMeanError();
					
					double P_total_real_dEta2 = PvsdEta[mult][28][sign]->GetMean();
					double P_total_real_dEta_error2 = PvsdEta[mult][28][sign]->GetMeanError();

					double weight1 = delEta3p[mult][sign]->GetBinContent( 28 );
					double weight2 = delEta3p[mult][sign]->GetBinContent( 29 );

					double value = weightedAverage(weight1, weight2, 0, P_total_real_dEta1, P_total_real_dEta2, 0);
					double error = weightedAverageError(weight1, weight2, 0, P_total_real_dEta_error1, P_total_real_dEta_error2, 0 );
					
					hist2[mult][sign]->SetBinContent(deta+1, value );
					hist2[mult][sign]->SetBinError(deta+1,  error);
				}

				for(int HF = 0; HF < 2; HF++){

					if(deta < 9){

						double Q_total_real_dEta1 = QvsdEta[mult][3*deta][sign][HF]->GetMean();
						double Q_total_real_dEta_error1 = QvsdEta[mult][3*deta][sign][HF]->GetMeanError();

						double Q_total_real_dEta2 = QvsdEta[mult][3*deta+1][sign][HF]->GetMean();
						double Q_total_real_dEta_error2 = QvsdEta[mult][3*deta+1][sign][HF]->GetMeanError();

						double Q_total_real_dEta3 = QvsdEta[mult][3*deta+2][sign][HF]->GetMean();
						double Q_total_real_dEta_error3 = QvsdEta[mult][3*deta+2][sign][HF]->GetMeanError();

						double weight1 = delEta3p[mult][sign]->GetBinContent( 3*deta+1 );
						double weight2 = delEta3p[mult][sign]->GetBinContent( 3*deta+2 );
						double weight3 = delEta3p[mult][sign]->GetBinContent( 3*deta+3 );
						
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

						double weight1 = delEta3p[mult][sign]->GetBinContent( 28 );
						double weight2 = delEta3p[mult][sign]->GetBinContent( 29 );

						double value = weightedAverage(weight1, weight2, 0, Q_total_real_dEta1, Q_total_real_dEta2, 0);
						double error = weightedAverageError(weight1, weight2, 0, Q_total_real_dEta_error1, Q_total_real_dEta_error2, 0 );
						
						hist1[mult][sign][HF]->SetBinContent(deta+1, value );
						hist1[mult][sign][HF]->SetBinError(deta+1,  error);
					}


				}
			}
		}


		for(int dpt = 0; dpt < NdPtBinsRebins1; dpt++){

			for(int sign = 0; sign < 3; sign++){

				double P_dpT_value_1 = PvsdPt[mult][2*dpt][sign]->GetMean();
				double P_dpT_error_1 = PvsdPt[mult][2*dpt][sign]->GetMeanError();
				
				double P_dpT_value_2 = PvsdPt[mult][2*dpt+1][sign]->GetMean();
				double P_dpT_error_2 = PvsdPt[mult][2*dpt+1][sign]->GetMeanError();

				double weight1 = delPt2p[mult][sign]->GetBinContent( 2*dpt );
				double weight2 = delPt2p[mult][sign]->GetBinContent( 2*dpt+1 );

				double value = weightedAverage(weight1, weight2, 0, P_dpT_value_1, P_dpT_value_2, 0);
				double error = weightedAverageError(weight1, weight2, 0, P_dpT_error_1, P_dpT_error_2, 0 );

				hist2_dpT[mult][sign]->SetBinContent(dpt+1, value);
				hist2_dpT[mult][sign]->SetBinError(dpt+1, error);

				for(int HF = 0; HF < 2; HF++){

					double Q_dpT_value_1 = QvsdPt[mult][2*dpt][sign][HF]->GetMean();
					double Q_dpT_error_1 = QvsdPt[mult][2*dpt][sign][HF]->GetMeanError();

					double Q_dpT_value_2 = QvsdPt[mult][2*dpt+1][sign][HF]->GetMean();
					double Q_dpT_error_2 = QvsdPt[mult][2*dpt+1][sign][HF]->GetMeanError();

					double weight1 = delPt3p[mult][sign]->GetBinContent( 2*dpt );
					double weight2 = delPt3p[mult][sign]->GetBinContent( 2*dpt+1 );

					double value = weightedAverage(weight1, weight2, 0, Q_dpT_value_1, Q_dpT_value_2, 0);
					double error = weightedAverageError(weight1, weight2, 0, Q_dpT_error_1, Q_dpT_error_2, 0 );
						
					hist1_dpT[mult][sign][HF]->SetBinContent(dpt+1, value);
					hist1_dpT[mult][sign][HF]->SetBinError(dpt+1, error);

				}

			}
		}
	
		for(int dpt = 0; dpt < NdPtBinsRebins2; dpt++){

			for(int sign = 0; sign < 3; sign++){

				double P_pTave_value_1 = PvsPtAve[mult][2*dpt+3][sign]->GetMean();
				double P_pTave_error_1 = PvsPtAve[mult][2*dpt+3][sign]->GetMeanError();

				double P_pTave_value_2 = PvsPtAve[mult][2*dpt+1+3][sign]->GetMean();
				double P_pTave_error_2 = PvsPtAve[mult][2*dpt+1+3][sign]->GetMeanError();

				double weight1 = ptAve2p[mult][sign]->GetBinContent( 2*dpt+3 );
				double weight2 = ptAve2p[mult][sign]->GetBinContent( 2*dpt+1+3 );

				double value = weightedAverage(weight1, weight2, 0, P_pTave_value_1, P_pTave_value_2, 0);
				double error = weightedAverageError(weight1, weight2, 0, P_pTave_error_1, P_pTave_error_2, 0 );

				hist2_pTave[mult][sign]->SetBinContent(dpt+1, value);
				hist2_pTave[mult][sign]->SetBinError(dpt+1, error);

				for(int HF = 0; HF < 2; HF++){

					double Q_pTave_value_1 = QvsPtAve[mult][2*dpt+3][sign][HF]->GetMean();
					double Q_pTave_error_1 = QvsPtAve[mult][2*dpt+3][sign][HF]->GetMeanError();

					double Q_pTave_value_2 = QvsPtAve[mult][2*dpt+1+3][sign][HF]->GetMean();
					double Q_pTave_error_2 = QvsPtAve[mult][2*dpt+1+3][sign][HF]->GetMeanError();

					double weight1 = ptAve3p[mult][sign]->GetBinContent( 2*dpt+3 );
					double weight2 = ptAve3p[mult][sign]->GetBinContent( 2*dpt+1+3 );

					double value = weightedAverage(weight1, weight2, 0, Q_pTave_value_1, Q_pTave_value_2, 0);
					double error = weightedAverageError(weight1, weight2, 0, Q_pTave_error_1, Q_pTave_error_2, 0 );

					hist1_pTave[mult][sign][HF]->SetBinContent(dpt+1, value);
					hist1_pTave[mult][sign][HF]->SetBinError(dpt+1, error);
				}

			}
		}
	}

//begin of deta

	TH1D* base1 = makeHist("base1", "", "|#Delta#eta|", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 48,0,4.8,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0012,0.0014);
	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);


	TH1D* temp1 = (TH1D*)hist1[0][0][0]->Clone("temp1");
	temp1->Add(hist1[0][1][0], +1);
	temp1->Scale(0.5);
	temp1->Scale(1.0/v2[0][0]);
	temp1->SetMarkerStyle(20);
	temp1->SetMarkerSize(1.4);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp2 = (TH1D*) hist1[0][2][0]->Clone("temp2");
	temp2->SetMarkerStyle(21);
	temp2->Scale(1.0/v2[0][0]);
	temp2->SetMarkerColor(kBlue);
	temp2->SetMarkerSize(1.4);
	temp2->SetLineColor(kBlue);

	TH1D* temp3 = (TH1D*)hist1[0][0][1]->Clone("temp3");
	temp3->Add(hist1[0][1][1], +1);
	temp3->Scale(0.5);
	temp3->Scale(1.0/v2[0][1]);
	temp3->SetMarkerStyle(24);
	temp3->SetMarkerSize(1.4);
	temp3->SetMarkerColor(kRed);
	temp3->SetLineColor(kRed);

	TH1D* temp4 = (TH1D*) hist1[0][2][1]->Clone("temp4");
	temp4->SetMarkerStyle(25);
	temp4->Scale(1.0/v2[0][1]);	
	temp4->SetMarkerSize(1.3);
	temp4->SetMarkerColor(kBlue);
	temp4->SetLineColor(kBlue);


	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	base1->Draw();
	temp1->Draw("Psame");
	temp2->Draw("Psame");
	temp3->Draw("Psame");
	temp4->Draw("Psame");

	TH1D* base2 = makeHist("base2", "", "|#Delta#eta|", "#LTcos(#phi_{#alpha}-#phi_{#beta})#GT", 48,0,4.8,kBlack);

	base2->GetYaxis()->SetRangeUser(-0.012,0.014);
	base2->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetXaxis()->SetTitleOffset(0.95);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);


	TH1D* temp5 = (TH1D*)hist2[0][0]->Clone("temp5");
	temp5->Add(hist2[0][1], +1);
	temp5->Scale(0.5);
	temp5->SetMarkerStyle(20);
	temp5->SetMarkerSize(1.4);
	temp5->SetMarkerColor(kRed);
	temp5->SetLineColor(kRed);

	TH1D* temp6 = (TH1D*) hist2[0][2]->Clone("temp6");
	temp6->SetMarkerStyle(21);
	temp6->SetMarkerColor(kBlue);
	temp6->SetMarkerSize(1.4);
	temp6->SetLineColor(kBlue);


	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	base2->Draw();
	temp5->Draw("Psame");
	temp6->Draw("Psame");

//end of deta

//begin of dpT

	TH1D* base3 = makeHist("base3", "", "|#Deltap_{T}|", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 48,0,3.0,kBlack);

	base3->GetYaxis()->SetRangeUser(-0.0012,0.0014);
	base3->GetXaxis()->SetRangeUser(-0.2, 3.0);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.3);
	base3->GetXaxis()->SetTitleOffset(0.95);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.3);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.4);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.4);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.4);

	TH1D* temp11 = (TH1D*)hist1_dpT[0][0][0]->Clone("temp11");
	temp11->Add(hist1_dpT[0][1][0], +1);
	temp11->Scale(0.5);
	temp11->Scale(1.0/v2[0][0]);
	temp11->SetMarkerStyle(20);
	temp11->SetMarkerSize(1.4);
	temp11->SetMarkerColor(kRed);
	temp11->SetLineColor(kRed);

	TH1D* temp12 = (TH1D*) hist1_dpT[0][2][0]->Clone("temp12");
	temp12->SetMarkerStyle(21);
	temp12->Scale(1.0/v2[0][0]);
	temp12->SetMarkerColor(kBlue);
	temp12->SetMarkerSize(1.4);
	temp12->SetLineColor(kBlue);

	TH1D* temp13 = (TH1D*)hist1_dpT[0][0][1]->Clone("temp13");
	temp13->Add(hist1_dpT[0][1][1], +1);
	temp13->Scale(0.5);
	temp13->Scale(1.0/v2[0][1]);
	temp13->SetMarkerStyle(24);
	temp13->SetMarkerSize(1.4);
	temp13->SetMarkerColor(kRed);
	temp13->SetLineColor(kRed);

	TH1D* temp14 = (TH1D*) hist1_dpT[0][2][1]->Clone("temp14");
	temp14->SetMarkerStyle(25);
	temp14->Scale(1.0/v2[0][1]);	
	temp14->SetMarkerSize(1.3);
	temp14->SetMarkerColor(kBlue);
	temp14->SetLineColor(kBlue);


	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	base3->Draw();
	temp11->Draw("Psame");
	temp12->Draw("Psame");
	temp13->Draw("Psame");
	temp14->Draw("Psame");

	TH1D* base4 = makeHist("base4", "", "|#Deltap_{T}|", "#LTcos(#phi_{#alpha}-#phi_{#beta})#GT", 48,0,3.0,kBlack);

	base4->GetYaxis()->SetRangeUser(-0.012,0.014);
	base4->GetXaxis()->SetRangeUser(-0.2, 3.0);
	base4->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base4,1.1,1.25);

	base4->GetYaxis()->SetTitleOffset(1.3);
	base4->GetXaxis()->SetTitleOffset(0.95);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.3);
	base4->GetXaxis()->SetTitleSize(base4->GetXaxis()->GetTitleSize()*1.4);
	base4->GetYaxis()->SetLabelSize(base4->GetYaxis()->GetLabelSize()*1.4);
	base4->GetXaxis()->SetLabelSize(base4->GetXaxis()->GetLabelSize()*1.4);


	TH1D* temp15 = (TH1D*)hist2_dpT[0][0]->Clone("temp15");
	temp15->Add(hist2_dpT[0][1], +1);
	temp15->Scale(0.5);
	temp15->SetMarkerStyle(20);
	temp15->SetMarkerSize(1.4);
	temp15->SetMarkerColor(kRed);
	temp15->SetLineColor(kRed);

	TH1D* temp16 = (TH1D*) hist2_dpT[0][2]->Clone("temp16");
	temp16->SetMarkerStyle(21);
	temp16->SetMarkerColor(kBlue);
	temp16->SetMarkerSize(1.4);
	temp16->SetLineColor(kBlue);


	TCanvas* c4 = new TCanvas("c4","c4",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	base4->Draw();
	temp15->Draw("Psame");
	temp16->Draw("Psame");

//end of dpT

//begin of pTave

	TH1D* base5 = makeHist("base5", "", "(p_{T,#alpha} #plus p_{T,#beta})/2", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 48,0,3.0,kBlack);

	base5->GetYaxis()->SetRangeUser(-0.0022,0.0044);
	base5->GetXaxis()->SetRangeUser(-0.2, 3.0);
	base5->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base5,1.1,1.25);

	base5->GetYaxis()->SetTitleOffset(1.3);
	base5->GetXaxis()->SetTitleOffset(0.95);
	base5->GetYaxis()->SetTitleSize(base5->GetYaxis()->GetTitleSize()*1.3);
	base5->GetXaxis()->SetTitleSize(base5->GetXaxis()->GetTitleSize()*1.4);
	base5->GetYaxis()->SetLabelSize(base5->GetYaxis()->GetLabelSize()*1.4);
	base5->GetXaxis()->SetLabelSize(base5->GetXaxis()->GetLabelSize()*1.4);

	TH1D* temp21 = (TH1D*)hist1_pTave[0][0][0]->Clone("temp21");
	temp21->Add(hist1_pTave[0][1][0], +1);
	temp21->Scale(0.5);
	temp21->Scale(1.0/v2[0][0]);
	temp21->SetMarkerStyle(20);
	temp21->SetMarkerSize(1.4);
	temp21->SetMarkerColor(kRed);
	temp21->SetLineColor(kRed);

	TH1D* temp22 = (TH1D*) hist1_pTave[0][2][0]->Clone("temp22");
	temp22->SetMarkerStyle(21);
	temp22->Scale(1.0/v2[0][0]);
	temp22->SetMarkerColor(kBlue);
	temp22->SetMarkerSize(1.4);
	temp22->SetLineColor(kBlue);

	TH1D* temp23 = (TH1D*)hist1_pTave[0][0][1]->Clone("temp23");
	temp23->Add(hist1_pTave[0][1][1], +1);
	temp23->Scale(0.5);
	temp23->Scale(1.0/v2[0][1]);
	temp23->SetMarkerStyle(24);
	temp23->SetMarkerSize(1.4);
	temp23->SetMarkerColor(kRed);
	temp23->SetLineColor(kRed);

	TH1D* temp24 = (TH1D*) hist1_pTave[0][2][1]->Clone("temp24");
	temp24->SetMarkerStyle(25);
	temp24->Scale(1.0/v2[0][1]);	
	temp24->SetMarkerSize(1.3);
	temp24->SetMarkerColor(kBlue);
	temp24->SetLineColor(kBlue);


	TCanvas* c5 = new TCanvas("c5","c5",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	base5->Draw();
	temp21->Draw("Psame");
	temp22->Draw("Psame");
	temp23->Draw("Psame");
	temp24->Draw("Psame");

	TH1D* base6 = makeHist("base6", "", "(p_{T,#alpha} #plus p_{T,#beta})/2", "#LTcos(#phi_{#alpha}-#phi_{#beta})#GT", 48,0,3.0,kBlack);

	base6->GetYaxis()->SetRangeUser(-0.022,0.024);
	base6->GetXaxis()->SetRangeUser(-0.2, 3.0);
	base6->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base6,1.1,1.25);

	base6->GetYaxis()->SetTitleOffset(1.3);
	base6->GetXaxis()->SetTitleOffset(0.95);
	base6->GetYaxis()->SetTitleSize(base6->GetYaxis()->GetTitleSize()*1.3);
	base6->GetXaxis()->SetTitleSize(base6->GetXaxis()->GetTitleSize()*1.4);
	base6->GetYaxis()->SetLabelSize(base6->GetYaxis()->GetLabelSize()*1.4);
	base6->GetXaxis()->SetLabelSize(base6->GetXaxis()->GetLabelSize()*1.4);


	TH1D* temp25 = (TH1D*)hist2_pTave[0][0]->Clone("temp25");
	temp25->Add(hist2_pTave[0][1], +1);
	temp25->Scale(0.5);
	temp25->SetMarkerStyle(20);
	temp25->SetMarkerSize(1.4);
	temp25->SetMarkerColor(kRed);
	temp25->SetLineColor(kRed);

	TH1D* temp26 = (TH1D*) hist2_pTave[0][2]->Clone("temp26");
	temp26->SetMarkerStyle(21);
	temp26->SetMarkerColor(kBlue);
	temp26->SetMarkerSize(1.4);
	temp26->SetLineColor(kBlue);


	TCanvas* c6 = new TCanvas("c6","c6",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	base6->Draw();
	temp25->Draw("Psame");
	temp26->Draw("Psame");

//end of pTave

/*
Save in a root file
*/

 // 	TFile f1("../dataPoints/CME_differentials_123_pPb_8TeV.root", "RECREATE");
 
	// temp1->Write();
	// temp2->Write();
	// temp3->Write();
	// temp4->Write();
	// temp5->Write();
	// temp6->Write();

	// temp11->Write();
	// temp12->Write();
	// temp13->Write();
	// temp14->Write();
	// temp15->Write();
	// temp16->Write();

	// temp21->Write();
	// temp22->Write();
	// temp23->Write();
	// temp24->Write();
	// temp25->Write();
	// temp26->Write();


}