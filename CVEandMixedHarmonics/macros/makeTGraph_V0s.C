#include "RiceStyle.h"
#include "../V0s/RooFit_La.h"
#include "../V0s/RooFit_Poly_Ks.h"
#include "../V0s/fitting.h"

using namespace std;

double PbPb_centralityBinCenter[] = {35, 45, 55, 65, 75};
int mult_ = 1;

double V0sV0s_correlator( double gamma_obs, double gamma_sig_bkg, double gamma_bkg_sig, double gamma_bkg_bkg, double f_sig_1, double f_sig_2){

	double f_bkg_1 = 1.0-f_sig_1;
	double f_bkg_2 = 1.0-f_sig_2;

	double numerator = gamma_obs - f_sig_1*f_bkg_2*gamma_sig_bkg + f_sig_2*f_bkg_1*gamma_bkg_sig - f_bkg_1*f_bkg_2*gamma_bkg_bkg;
	double denominator = f_sig_1*f_sig_2;

	return numerator/denominator;

}

double V0sH_correlator(double gamma_obs, double gamma_bkg, double f_sig){

	double f_bkg = 1.0 - f_sig;
	double numerator = gamma_obs - f_bkg*gamma_bkg;
	return numerator/f_sig;
}

void makeTGraph_V0s(){

	TFile* file1[6];
/*
Fit V0s under-fly 
*/
	
	TH3D* ks_mass[6];
	TH3D* la_mass_1[6];
	TH3D* la_mass_2[6];

	TH1D* ks_mass_1D[6];
	TH1D* la_mass_1_1D[6];
	TH1D* la_mass_2_1D[6];
	
	for(int i = 0; i < mult_; i++){

		file1[i] = new TFile(Form("../rootfiles/CVEandMixedHarmonics_PbPb_30_100_v3_%d.root", i+1));
		ks_mass[i] = (TH3D*) file1[i]->Get("ana/ks_mass");
		la_mass_1[i] = (TH3D*) file1[i]->Get("ana/la_mass_1");
		la_mass_2[i] = (TH3D*) file1[i]->Get("ana/la_mass_2");

		ks_mass_1D[i] = (TH1D*) ks_mass[i]->ProjectionZ(Form("ks_mass_1D_%d", i), 11,59,1,50);//|eta| < 2.4, 0.0<pT<5.0
		la_mass_1D_1[i] = (TH1D*) la_mass_1[i]->ProjectionZ(Form("la_mass_1D_1_%d", i), 11,59,4,50);//|eta| < 2.4, 0.4<pT<5.0
		la_mass_1D_2[i] = (TH1D*) la_mass_2[i]->ProjectionZ(Form("la_mass_1D_2_%d", i), 11,59,4,50);//|eta| < 2.4, 0.4<pT<5.0
	}

	double f_sig_K0s[6];
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(3,2);
	for(int i = 0; i < mult_; i++){

		c1->cd(i+1);
		vector<double> fitParameters_K0s;
		InitialFit(c1, ks_mass_1D[i], 1, 1, 1, fitParameters_K0s);
		f_sig_K0s[i] = fitParameters_K0s[6];
	}
	double f_sig_Lam[6];
	TCanvas* c2 = new TCanvas("c2","c2",800,800);
	for(int i = 0; i < mult_; i++){

		vector<double> fitParameters_Lam;
		InitialFit_la(c2, la_mass_1D_1[i], 1, 1, 1, fitParameters_Lam);
		f_sig_Lam[i] = fitParameters_Lam[5];
	}

	TH1D* QaQb[5]; TH1D* QaQc[5]; TH1D* QcQb[5];

	for(int mult = 0; mult < mult_; mult++){

		QaQb[mult] = (TH1D*)file1[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file1[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file1[mult]->Get("ana/c2_cb");
	}

	double v2[23][3];//get corrected v2_3

	for(int mult = 0; mult < mult_; mult++){
		
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

	TH1D* c3_LL_real[6][3][4][2];
	double c3_LL_value[6][3][4][2];
	double c3_LL_error[6][3][4][2];

	TH1D* c3_KK_real[6][4][2];
	double c3_KK_value[6][4][2];
	double c3_KK_error[6][4][2];

	TH1D* c3_LK_real[6][3][4][2];
	double c3_LK_value[6][3][4][2];
	double c3_LK_error[6][3][4][2];

	TH1D* c3_LH_real[6][3][4][2];
	double c3_LH_value[6][3][4][2];
	double c3_LH_error[6][3][4][2];

	TH1D* c3_KH_real[6][4][2];
	double c3_KH_value[6][4][2];
	double c3_KH_error[6][4][2];

	for(int mult = 0; mult < mult_; mult++){
		for(int sign = 0; sign < 3; sign++){//lambda baryon number
			for(int comb = 0; comb < 4; comb++){//sig bkg combination
				for(int HF = 0; HF < 2; HF++){//HF sides

					c3_LL_real[mult][sign][comb][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_LL_real_%d_%d_%d", sign,comb,HF));
					c3_LL_value[mult][sign][comb][HF] = c3_LL_real[mult][sign][comb][HF]->GetMean();
					c3_LL_error[mult][sign][comb][HF] = c3_LL_real[mult][sign][comb][HF]->GetMeanError();
				
					c3_LK_real[mult][sign][comb][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_LK_real_%d_%d_%d", sign,comb,HF));
					c3_LK_value[mult][sign][comb][HF] = c3_LK_real[mult][sign][comb][HF]->GetMean();
					c3_LK_error[mult][sign][comb][HF] = c3_LK_real[mult][sign][comb][HF]->GetMeanError();

				}
			}

			for(int comb = 0; comb < 2; comb++){//sig bkg combination
				for(int HF = 0; HF < 2; HF++){//HF sides

					c3_LH_real[mult][sign][comb][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_LH_real_%d_%d_%d", sign,comb,HF));
					c3_LH_value[mult][sign][comb][HF] = c3_LH_real[mult][sign][comb][HF]->GetMean();
					c3_LH_error[mult][sign][comb][HF] = c3_LH_real[mult][sign][comb][HF]->GetMeanError();
				}
			}

		}
		
		for(int comb = 0; comb < 4; comb++){//sig bkg combination
			for(int HF = 0; HF < 2; HF++){//HF sides

				c3_KK_real[mult][comb][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_KK_real_%d_%d",comb,HF));
				c3_KK_value[mult][comb][HF] = c3_KK_real[mult][comb][HF]->GetMean();
				c3_KK_error[mult][comb][HF] = c3_KK_real[mult][comb][HF]->GetMeanError();
		
			}
		}
		
		for(int comb = 0; comb < 2; comb++){//sig bkg combination
			for(int HF = 0; HF < 2; HF++){//HF sides

				c3_KH_real[mult][comb][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_KH_real_%d_%d",comb,HF));
				c3_KH_value[mult][comb][HF] = c3_KH_real[mult][comb][HF]->GetMean();
				c3_KH_error[mult][comb][HF] = c3_KH_real[mult][comb][HF]->GetMeanError();
		
			}
		}
	}

/*
Setting up the TGraphs
*/


	TGraphErrors* gr1_Pb[3];
	TGraphErrors* gr1_p[3];

	TGraphErrors* gr3_Pb[3];
	TGraphErrors* gr3_p[3];

	TGraphErrors* gr4_Pb[2];
	TGraphErrors* gr4_p[2];

	for(int i = 0; i < 3; i++){
		
		gr1_Pb[i] = new TGraphErrors(6);
		gr1_Pb[i]->SetName(Form("la_la_Pb_%d", i));

		gr1_p[i] = new TGraphErrors(6);
		gr1_p[i]->SetName(Form("la_la_p_%d", i));

		gr3_Pb[i] = new TGraphErrors(6);
		gr3_Pb[i]->SetName(Form("la_ks_Pb_%d", i));

		gr3_p[i] = new TGraphErrors(6);
		gr3_p[i]->SetName(Form("la_ks_p_%d", i));

	}

	for(int i = 0; i < 2; i++){

		gr4_Pb[i] = new TGraphErrors(6);
		gr4_Pb[i]->SetName(Form("la_h_Pb_%d", i));

		gr4_p[i] = new TGraphErrors(6);
		gr4_p[i]->SetName(Form("la_h_p_%d", i));		
	}
	

	TGraphErrors* gr2_Pb = new TGraphErrors(6);
	gr2_Pb->SetName("ks_ks_Pb");
	TGraphErrors* gr2_p = new TGraphErrors(6);
	gr2_p->SetName("ks_ks_p");

	TGraphErrors* gr5_Pb = new TGraphErrors(6);
	gr5_Pb->SetName("ks_h_Pb");
	TGraphErrors* gr5_p = new TGraphErrors(6);
	gr5_p->SetName("ks_h_p");
	
/*
Filling the TGraphs with place-holder errors
*/

	for(int mult = 0; mult < mult_; mult++){

		for(int i = 0; i < 3; i++){
			//lambda-lambda correlation:
			//Pb-going:
			double Pb_value = V0sV0s_correlator(c3_LL_value[mult][i][0][0], c3_LL_value[mult][i][1][0], c3_LL_value[mult][i][2][0], c3_LL_value[mult][i][3][0], f_sig_Lam[mult], f_sig_Lam[mult]);
			double Pb_error = c3_LL_error[mult][i][0][0];

			gr1_Pb[i]->SetPoint(mult, PbPb_centralityBinCenter[mult], Pb_value/v2[mult][0]);
			gr1_Pb[i]->SetPointError(mult, PbPb_centralityBinCenter[mult], Pb_error/v2[mult][0]);

			//p-going:
			double p_value = V0sV0s_correlator(c3_LL_value[mult][i][0][1], c3_LL_value[mult][i][1][1], c3_LL_value[mult][i][2][1], c3_LL_value[mult][i][3][1], f_sig_Lam[mult], f_sig_Lam[mult]);
			double p_error = c3_LL_error[mult][i][0][1];

			gr1_p[i]->SetPoint(mult, PbPb_centralityBinCenter[mult], p_value/v2[mult][1]);
			gr1_p[i]->SetPointError(mult, PbPb_centralityBinCenter[mult], p_error/v2[mult][1]);

		}

		for(int i = 0; i < 2; i++){

			//lambda-ks correlation:
			//Pb-going:
			double Pb_value = V0sV0s_correlator(c3_LK_value[mult][i][0][0], c3_LK_value[mult][i][1][0], c3_LK_value[mult][i][2][0], c3_LK_value[mult][i][3][0], f_sig_Lam[mult], f_sig_K0s[mult]);
			double Pb_error = c3_LK_error[mult][i][0][0];

			gr3_Pb[i]->SetPoint(mult, PbPb_centralityBinCenter[mult], Pb_value/v2[mult][0]);
			gr3_Pb[i]->SetPointError(mult, PbPb_centralityBinCenter[mult], Pb_error/v2[mult][0]);

			//p-going:
			double p_value = V0sV0s_correlator(c3_LK_value[mult][i][0][1], c3_LK_value[mult][i][1][1], c3_LK_value[mult][i][2][1], c3_LK_value[mult][i][3][1], f_sig_Lam[mult], f_sig_K0s[mult]);
			double p_error = c3_LK_error[mult][i][0][1];

			gr3_p[i]->SetPoint(mult, PbPb_centralityBinCenter[mult], p_value/v2[mult][1]);
			gr3_p[i]->SetPointError(mult, PbPb_centralityBinCenter[mult], p_error/v2[mult][1]);

			//lambda-H correlation:
			//Pb-going:
			double Pb_value = V0sH_correlator(c3_LH_value[mult][i][0][0], c3_LH_value[mult][i][1][0], f_sig_Lam[mult]);
			double Pb_error = c3_LH_error[mult][i][0][0];

			gr4_Pb[i]->SetPoint(mult, PbPb_centralityBinCenter[mult], Pb_value/v2[mult][0]);
			gr4_Pb[i]->SetPointError(mult, PbPb_centralityBinCenter[mult], Pb_error/v2[mult][0]);

			//p-going:
			double p_value = V0sH_correlator(c3_LH_value[mult][i][0][1], c3_LH_value[mult][i][1][1], f_sig_K0s[mult]);
			double p_error = c3_LH_error[mult][i][0][1];

			gr4_p[i]->SetPoint(mult, PbPb_centralityBinCenter[mult], p_value/v2[mult][1]);
			gr4_p[i]->SetPointError(mult, PbPb_centralityBinCenter[mult], p_error/v2[mult][1]);
		}
	
		//K0s-K0s correlation:
		//Pb-going:
		double Pb_value = V0sV0s_correlator(c3_KK_value[mult][0][0], c3_KK_value[mult][1][0], c3_KK_value[mult][2][0], c3_KK_value[mult][3][0], f_sig_K0s[mult], f_sig_K0s[mult]);
		double Pb_error = c3_KK_error[mult][0][0];

		gr2_Pb->SetPoint(mult, PbPb_centralityBinCenter[mult], Pb_value/v2[mult][0]);
		gr2_Pb->SetPointError(mult, PbPb_centralityBinCenter[mult], Pb_error/v2[mult][0]);

		//p-going:
		double p_value = V0sV0s_correlator(c3_KK_value[mult][0][1], c3_KK_value[mult][1][1], c3_KK_value[mult][2][1], c3_KK_value[mult][3][1], f_sig_K0s[mult], f_sig_K0s[mult]);
		double p_error = c3_KK_error[mult][0][1];

		gr2_p->SetPoint(mult, PbPb_centralityBinCenter[mult], p_value/v2[mult][1]);
		gr2_p->SetPointError(mult, PbPb_centralityBinCenter[mult], p_error/v2[mult][1]);


		//K0s-H correlation:
		//Pb-going:
		double Pb_value = V0sH_correlator(c3_KH_value[mult][0][0], c3_KH_value[mult][1][0], f_sig_K0s[mult]);
		double Pb_error = c3_KH_error[mult][0][0];

		gr5_Pb->SetPoint(mult, PbPb_centralityBinCenter[mult], Pb_value/v2[mult][0]);
		gr5_Pb->SetPointError(mult, PbPb_centralityBinCenter[mult], Pb_error/v2[mult][0]);

		//p-going:
		double p_value = V0sH_correlator(c3_KH_value[mult][0][1], c3_KH_value[mult][1][1], f_sig_K0s[mult]);
		double p_error = c3_KH_error[mult][0][1];

		gr5_p->SetPoint(mult, PbPb_centralityBinCenter[mult], p_value/v2[mult][1]);
		gr5_p->SetPointError(mult, PbPb_centralityBinCenter[mult], p_error/v2[mult][1]);

		
	}

	TFile f1("../dataPoints/V0s_PbPb_test.root","RECREATE");
	for(int i = 0; i < 3; i++){
		gr1_Pb[i]->Write();
		gr1_p[i]->Write();
	}
	for(int i = 0; i < 2; i++){

		gr3_Pb[i]->Write();
		gr3_p[i]->Write();
		
		gr4_Pb[i]->Write();
		gr4_p[i]->Write();
	}

	gr2_Pb->Write();
	gr2_p->Write();

	gr5_Pb->Write();
	gr5_p->Write();


}