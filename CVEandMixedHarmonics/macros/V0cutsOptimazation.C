#include "RiceStyle.h"
//#include "../V0s/RooFit_La.h"
#include "../V0s/RooFit_La_original.h"
#include "../V0s/RooFit_Poly_Ks.h"
#include "../V0s/fitting.h"

using namespace std;

void V0cutsOptimazation(){

	TFile* file = new TFile("../rootfiles/V0Ntuple_small.root");
	TNtuple* ks_tree = (TNtuple*) file->Get("ana/V0AnalyzerSimpleNtuple_ks");
	TNtuple* la_tree = (TNtuple*) file->Get("ana/V0AnalyzerSimpleNtuple_la");

	stringstream common_cut[6][6][6][6];
	stringstream la_cut_DecayLength[10];
	stringstream la_cut_PAngle[10];
	stringstream la_cut_DCA[10];
	stringstream mult_cut[6];

	double dlos_cut[] = {5.0,6.0,7.0,8.0,9.0,10.0,11.0};
	double pangle_cut[] = {0.999, 0.9992, 0.9994, 0.9996, 0.9998, 0.9999};
	double dca_cut[] = {1.0,1.2,1.5,1.8,2.0,2.2};

	for(int mult = 4; mult < 6; mult++){
		for(int i = 0; i < 6; i++){
			for(int j = 0; j < 6; j++){
				for(int k = 0; k < 6; k++){
					
					common_cut[mult][i][j][k].str("");
					common_cut[mult][i][j][k] << "pt > 0.3 && pt < 5.0 && TMath::Abs(eta) < 2.4 && trkNHits1 > 3 &&  trkNHits2 > 3";
				}
			}
		}
	}

	for(int i = 0; i < 6; i++){

		la_cut_DecayLength[i].str("");
		la_cut_PAngle[i].str("");
		la_cut_DCA[i].str("");
		mult_cut[i].str("");

		la_cut_DecayLength[i] << " && L > "; la_cut_DecayLength[i] << dlos_cut[i];
		la_cut_PAngle[i] << " && PAngle > "; la_cut_PAngle[i] << pangle_cut[i];
		la_cut_DCA[i] << " && TMath::Abs(trkDCA1z) > " << dca_cut[i] << " && TMath::Abs(trkDCA2z) > " << dca_cut[i] << " && TMath::Abs(trkDCA1xy) > "<< dca_cut[i] << " && TMath::Abs(trkDCA2xy) > " << dca_cut[i];
	}

	
	mult_cut[0] << "&& ntrk > 0 && ntrk < 50";
	mult_cut[1] << "&& ntrk > 20 && ntrk < 160";
	mult_cut[2] << "&& ntrk > 100 && ntrk < 300";
	mult_cut[3] << "&& ntrk > 200 && ntrk < 600";
	mult_cut[4] << "&& ntrk > 500 && ntrk < 900";
	mult_cut[5] << "&& ntrk > 800 && ntrk < 1500";

	TH1D* ks_hist[6][6][6][6];
	TH1D* la_hist[6][6][6][6];
	for(int mult = 4; mult < 6; mult++){
		for(int i = 0; i < 6; i++){
			for(int j = 0; j < 6; j++){
				for(int k = 0; k < 6; k++){
					
					ks_hist[mult][i][j][k] = new TH1D(Form("ks_hist_%d_%d_%d_%d",mult,i,j,k), Form("ks_hist_%d_%d_%d_%d",mult,i,j,k), 360, 0.44,0.56);
					//la_hist[mult][i][j][k] = new TH1D(Form("la_hist_%d_%d_%d_%d",mult,i,j,k), Form("la_hist_%d_%d_%d_%d",mult,i,j,k), 360,1.08,1.16);

					common_cut[mult][i][j][k] << la_cut_DecayLength[i].str().c_str() << la_cut_PAngle[j].str().c_str() << la_cut_DCA[k].str().c_str() << mult_cut[mult].str().c_str();
					ks_tree->Draw(Form("mass>>ks_hist_%d_%d_%d_%d",mult,i,j,k), (common_cut[mult][i][j][k]).str().c_str(), "goff" );

				}
			}
		}
	}

	TFile f1("../dataPoints/FullyScanK0sCuts_v1.root", "RECREATE");
	for(int mult = 4; mult < 6; mult++){
		for(int i = 0; i < 6; i++){
			for(int j = 0; j < 6; j++){
				for(int k = 0; k < 6; k++){

					ks_hist[mult][i][j][k]->Write();
				}
			}
		}
	}


	



}