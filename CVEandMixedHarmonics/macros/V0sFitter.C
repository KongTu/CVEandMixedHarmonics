#include "RiceStyle.h"
//#include "../V0s/RooFit_La.h"
#include "../V0s/RooFit_La_original.h"
#include "../V0s/RooFit_Poly_Ks.h"
#include "../V0s/fitting.h"

using namespace std;


void V0sFitter(){

	TFile* file1 = new TFile("../dataPoints/FullyScanK0sCuts_v1.root");
	TFile* file2 = new TFile("../dataPoints/FullyScanLambdaCuts_v1.root");
	
	TH1D* ks_hist[6][6][6][6];
	TH1D* la_hist[6][6][6][6];
	for(int mult = 4; mult < 6; mult++){
		for(int i = 0; i < 6; i++){
			for(int j = 0; j < 6; j++){
				for(int k = 0; k < 6; k++){
					ks_hist[mult][i][j][k] = (TH1D*) file1->Get(Form("ks_hist_%d_%d_%d_%d", mult,i,j,k));
					la_hist[mult][i][j][k] = (TH1D*) file2->Get(Form("la_hist_%d_%d_%d_%d", mult,i,j,k));

				}
			}
		}
	}

	double maxSignif = 0.0;
	int i_max = 0, j_max = 0, k_max = 0;

	TH1D* signif[2];
	signif[0] = new TH1D("signif1","",120,0,1200);
	signif[1] = new TH1D("signif2","",120,0,1200);
	TCanvas* c2[6][6];
	
	for(int mult = 4; mult < 6; mult++){
		for(int i = 0; i < 6; i++){

			c2[mult][i] = new TCanvas(Form("c2_%d_%d",mult,i),Form("c2_%d_%d",mult,i),900,600);
			c2[mult][i]->Divide(6,6);

			int l = 1;
			for(int j = 0; j < 6; j++){
				for(int k = 0; k < 6; k++){

					c2[mult][i]->cd(l);l++;
			
					vector<double> fitParameters;
					InitialFit(c2[0][i], ks_hist[mult][i][j][k], i, j, k, 1, fitParameters);
					signif[0]->Fill( fitParameters[7] );

					if( fitParameters[7] > maxSignif ){
						maxSignif = fitParameters[7];
						i_max = i;
						j_max = j;
						k_max = k;
					}

				}
			}

			if(mult == 4) c2[mult][i]->Print(Form("../fits/K0sFits_40_50_%d.pdf", i));
			if(mult == 5) c2[mult][i]->Print(Form("../fits/K0sFits_30_40_%d.pdf", i));

		}
	}


	// for(int mult = 4; mult < 6; mult++){
	// 	for(int i = 0; i < 6; i++){

	// 		c2[mult][i] = new TCanvas(Form("c2_%d_%d",mult,i),Form("c2_%d_%d",mult,i),900,600);
	// 		c2[mult][i]->Divide(6,6);

	// 		int l = 1;
	// 		for(int j = 0; j < 6; j++){
	// 			for(int k = 0; k < 6; k++){

	// 				c2[mult][i]->cd(l);l++;
	// 				vector<double> fitParameters;
	// 				fitParameters = la_YieldCal(la_hist[mult][i][j][k], i, j, k);
	// 				signif[0]->Fill( fitParameters[2] );

	// 				if( fitParameters[2] > maxSignif ){
	// 					maxSignif = fitParameters[2];
	// 					i_max = i;
	// 					j_max = j;
	// 					k_max = k;
	// 				}

	// 			}
	// 		}

	// 		if(mult == 4) c2[mult][i]->Print(Form("../fits/LambdaFits_40_50_%d.pdf", i));
	// 		if(mult == 5) c2[mult][i]->Print(Form("../fits/LambdaFits_30_40_%d.pdf", i));

	// 	}

	// }

	cout << "max signal significance: " << endl;
	cout << maxSignif << " with i = " << i_max << ", j = " << j_max << ", k = " << k_max << endl;

}