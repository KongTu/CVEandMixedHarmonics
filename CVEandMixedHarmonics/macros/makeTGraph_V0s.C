#include "RiceStyle.h"

using namespace std;

double PbPb_centralityBinCenter[] = {35, 45, 55, 65, 75};

double f_sig = 0.7;

void makeTGraph_V0s(){

	TFile* file1[5];
	
	for(int i = 0; i < 5; i++){

		file1[i] = new TFile(Form("../rootfiles/CVEandMixedHarmonics_PbPb_30_100_v2_%d.root", i+1));
	}

	TH1D* QaQb[5]; TH1D* QaQc[5]; TH1D* QcQb[5];

	for(int mult = 0; mult < 5; mult++){

		QaQb[mult] = (TH1D*)file1[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file1[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file1[mult]->Get("ana/c2_cb");
	}

	double v2[23][3];//get corrected v2_3

	for(int mult = 0; mult < 5; mult++){
		
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

	TH1D* c3_LL_real[5][2][2];
	TH1D* c3_KK_real[5][2][2];
	TH1D* c3_LK_real[5][2][2];
	TH1D* c3_LH_real[5][2][2];
	TH1D* c3_KH_real[5][2][2];

	for(int mult = 0; mult < 5; mult++){
		for(int sig = 0; sig < 2; sig++){
			for(int HF = 0; HF < 2; HF++){

				c3_LL_real[mult][sig][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_LL_real_%d_%d", sig,HF));
				c3_KK_real[mult][sig][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_KK_real_%d_%d", sig,HF));
				c3_LK_real[mult][sig][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_LK_real_%d_%d", sig,HF));
				c3_LH_real[mult][sig][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_LH_real_%d_%d", sig,HF));
				c3_KH_real[mult][sig][HF] = (TH1D*) file1[mult]->Get(Form("ana/c3_KH_real_%d_%d", sig,HF));
			
			}
		}
	}

	TGraphErrors* gr1 = new TGraphErrors(5);
	TGraphErrors* gr2 = new TGraphErrors(5);
	TGraphErrors* gr3 = new TGraphErrors(5);
	TGraphErrors* gr4 = new TGraphErrors(5);
	TGraphErrors* gr5 = new TGraphErrors(5);

	for(int mult = 0; mult < 5; mult++){

		//LL
		TH1D* temp_signal = (TH1D*) c3_LL_real[mult][0][0]->Clone(Form("test1_%d", mult));
		temp_signal->Add(c3_LL_real[mult][0][1], +1);

		double Q_sig = temp_signal->GetMean();
		double error_sig  = temp_signal->GetMeanError();

		TH1D* temp_bkg = (TH1D*) c3_LL_real[mult][1][0]->Clone(Form("test11_%d", mult));
		temp_bkg->Add(c3_LL_real[mult][1][1], +1);

		double bkg = temp_bkg->GetMean();
		double error_bkg  = temp_bkg->GetMeanError();

		double total = (Q_sig + (1-f_sig)*bkg)/f_sig;
		total = total/v2[mult][2];
		error_sig = error_sig/v2[mult][2];

		gr1->SetPoint(mult, PbPb_centralityBinCenter[mult], total);
		gr1->SetPointError(mult, 0.0, error_sig);

		//KK
		TH1D* temp_signal = (TH1D*) c3_KK_real[mult][0][0]->Clone(Form("test2_%d", mult));
		temp_signal->Add(c3_KK_real[mult][0][1], +1);

		double Q_sig = temp_signal->GetMean();
		double error_sig  = temp_signal->GetMeanError();

		TH1D* temp_bkg = (TH1D*) c3_KK_real[mult][1][0]->Clone(Form("test22_%d", mult));
		temp_bkg->Add(c3_KK_real[mult][1][1], +1);

		double bkg = temp_bkg->GetMean();
		double error_bkg  = temp_bkg->GetMeanError();

		double total = (Q_sig + (1-f_sig)*bkg)/f_sig;
		total = total/v2[mult][2];
		error_sig = error_sig/v2[mult][2];

		gr2->SetPoint(mult, PbPb_centralityBinCenter[mult], total);
		gr2->SetPointError(mult, 0.0, error_sig);

		//LK
		TH1D* temp_signal = (TH1D*) c3_LK_real[mult][0][0]->Clone(Form("test3_%d", mult));
		temp_signal->Add(c3_LK_real[mult][0][1], +1);

		double Q_sig = temp_signal->GetMean();
		double error_sig  = temp_signal->GetMeanError();

		TH1D* temp_bkg = (TH1D*) c3_LK_real[mult][1][0]->Clone(Form("test33_%d", mult));
		temp_bkg->Add(c3_LK_real[mult][1][1], +1);

		double bkg = temp_bkg->GetMean();
		double error_bkg  = temp_bkg->GetMeanError();

		double total = (Q_sig + (1-f_sig)*bkg)/f_sig;
		total = total/v2[mult][2];
		error_sig = error_sig/v2[mult][2];

		gr3->SetPoint(mult, PbPb_centralityBinCenter[mult], total);
		gr3->SetPointError(mult, 0.0, error_sig);

		//LH
		TH1D* temp_signal = (TH1D*) c3_LH_real[mult][0][0]->Clone(Form("test4_%d", mult));
		temp_signal->Add(c3_LH_real[mult][0][1], +1);

		double Q_sig = temp_signal->GetMean();
		double error_sig  = temp_signal->GetMeanError();

		TH1D* temp_bkg = (TH1D*) c3_LH_real[mult][1][0]->Clone(Form("test44_%d", mult));
		temp_bkg->Add(c3_LH_real[mult][1][1], +1);

		double bkg = temp_bkg->GetMean();
		double error_bkg  = temp_bkg->GetMeanError();

		double total = (Q_sig + (1-f_sig)*bkg)/f_sig;
		total = total/v2[mult][2];
		error_sig = error_sig/v2[mult][2];

		gr4->SetPoint(mult, PbPb_centralityBinCenter[mult], total);
		gr4->SetPointError(mult, 0.0, error_sig);

		//KH
		TH1D* temp_signal = (TH1D*) c3_KH_real[mult][0][0]->Clone(Form("test5_%d", mult));
		temp_signal->Add(c3_KH_real[mult][0][1], +1);

		double Q_sig = temp_signal->GetMean();
		double error_sig  = temp_signal->GetMeanError();

		TH1D* temp_bkg = (TH1D*) c3_KH_real[mult][1][0]->Clone(Form("test55_%d", mult));
		temp_bkg->Add(c3_KH_real[mult][1][1], +1);

		double bkg = temp_bkg->GetMean();
		double error_bkg  = temp_bkg->GetMeanError();

		double total = (Q_sig + (1-f_sig)*bkg)/f_sig;
		total = total/v2[mult][2];
		error_sig = error_sig/v2[mult][2];

		gr5->SetPoint(mult, PbPb_centralityBinCenter[mult], total);
		gr5->SetPointError(mult, 0.0, error_sig);


	}

	TFile f1("../dataPoints/V0s_PbPb.root","RECREATE");
	gr1->Write();
	gr2->Write();
	gr3->Write();
	gr4->Write();
	gr5->Write();






}