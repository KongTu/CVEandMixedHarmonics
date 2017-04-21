#include "RiceStyle.h"

using namespace std;


void compareV0Producer(){

	TFile* file[10];

	file[0] = new TFile("../rootfiles/CVEandMixedHarmonics_pPb_HM_v1_185_250.root");
	file[1] = new TFile("../rootfiles/CVEandMixedHarmonics_pPb_HM_v2_3.root");

	TH3D* la_mass_old[10];
	TH3D* la_mass_new[10];

	la_mass_old[0] = (TH3D*) file[0]->Get("ana/la_mass_1");
	la_mass_new[0] = (TH3D*) file[1]->Get("ana/la_mass_1");


	double range[] = {0.,4.,6.,8.,10,12,16,20,40,60,100};
	
	TH1D* la_mass_old_1D[20];
	TH1D* la_mass_new_1D[20];


	TCanvas* c1 = makeMultiCanvas("c1","c1", 3,3);

	for(int pt = 0; pt < 9; pt++){

		c1->cd(pt+1);
		gPad->SetLeftMargin(0.16);
	 	gPad->SetBottomMargin(0.13);
	 	gPad->SetTopMargin(0.1);
		gPad->SetTicks();

		la_mass_old_1D[pt] = (TH1D*) la_mass_old[0]->ProjectionZ(Form("la_mass_old_1D_%d", pt), 1, 70, range[pt], range[pt+1]);
		la_mass_new_1D[pt] = (TH1D*) la_mass_new[0]->ProjectionZ(Form("la_mass_new_1D_%d", pt), 1, 70, range[pt], range[pt+1]);

		la_mass_old_1D[pt]->SetMarkerStyle(20);
		la_mass_old_1D[pt]->SetMarkerColor(kBlack);

		la_mass_new_1D[pt]->SetMarkerStyle(20);
		la_mass_new_1D[pt]->SetMarkerColor(kRed);
		la_mass_new_1D[pt]->SetTitle(Form("%.1f < p_{T} < %.1f", range[pt]*0.1., range[pt+1]*0.1));
		la_mass_new_1D[pt]->Draw("p");
		la_mass_old_1D[pt]->Draw("psame");

	}

	c1->cd(1);
	TLegend *w3 = new TLegend(0.40,0.35,0.6,0.7);
    w3->SetLineColor(kWhite);
    w3->SetFillColor(0);
    w3->SetTextSize(20);
    w3->SetTextFont(43);
    w3->SetHeader("pPb 8TeV HM 185-250");
    w3->AddEntry(la_mass_old_1D[0], "generalV0sCollection");
    w3->AddEntry(la_mass_new_1D[0], "generalV0sCollectionNew");
    w3->Draw("same");





}