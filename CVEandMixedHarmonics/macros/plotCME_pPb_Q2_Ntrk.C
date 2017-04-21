#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {0,35,60,90,120,150,185,220,260};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
int ntrkBinCenter[] = {17.5, 47.5, 75, 105, 135, 167.5, 202.5, 240};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {132.1, 163, 202.4, 264, 314.9};

const int Nmults = 11;

double total_systematics_pPb = 0.0000;
double total_systematics_PbPb = 0.0000;

void plotCME_pPb_Q2_Ntrk(){

	gStyle->SetErrorX(0);

	TFile* file1[7];
    TGraphErrors* gr1[10][4];

    for(int mult = 0; mult < 5; mult++){
        
        file1[mult] = new TFile(Form("../dataPoints/pPb_8TeV_cme_q2_%d.root", mult+1));   
    
        for(int i = 0; i < 4; i++){

            gr1[mult][i] = (TGraphErrors*) file1[mult]->Get(Form("Graph;%d", i+1));
        }
    }
	
//start plotting

    TGraphErrors* gr_pPb_Pb[10];
    TGraphErrors* gr_pPb_p[10];

    for(int mult = 0; mult < 5; mult++){

        gr_pPb_Pb[mult] = new TGraphErrors(Nmults);
        gr_pPb_p[mult] = new TGraphErrors(Nmults);

        for(int q2 = 0; q2 < 11; q2++){

            double ye = total_systematics_pPb;

            double x1;
            double value1;
            double value1_error;
            gr1[mult][0]->GetPoint(q2, x1, value1);
            value1_error = gr1[mult][0]->GetErrorY(q2);

            double x2;
            double value2;
            double value2_error;
            gr1[mult][1]->GetPoint(q2, x2, value2);
            value2_error = gr1[mult][1]->GetErrorY(q2);

            double x = x1;
            double y = value2 - value1;
            double ey = sqrt(value1_error*value1_error + value2_error*value2_error);

            gr_pPb_Pb[mult]->SetPoint(q2, x, y);
            gr_pPb_Pb[mult]->SetPointError(q2, 0, ey);


            double x1;
            double value1;
            double value1_error;
            gr1[mult][2]->GetPoint(q2, x1, value1);
            value1_error = gr1[mult][2]->GetErrorY(q2);

            double x2;
            double value2;
            double value2_error;
            gr1[mult][3]->GetPoint(q2, x2, value2);
            value2_error = gr1[mult][3]->GetErrorY(q2);

            double x = x1;
            double y = value2 - value1;
            double ey = sqrt(value1_error*value1_error + value2_error*value2_error);

            gr_pPb_p[mult]->SetPoint(q2, x, y);
            gr_pPb_p[mult]->SetPointError(q2, 0, ey);

        }

    }



	gStyle->SetErrorX(0);


    TGraphErrors* f_cme_Pb = new TGraphErrors(5);
    TGraphErrors* f_cme_p = new TGraphErrors(5);

    for(int mult = 0; mult < 5; mult++){

        gr_pPb_Pb[mult]->Fit("pol1");
        
        TF1 * myFunc1 = gr_pPb_Pb[mult]->GetFunction("pol1");
        myFunc1->SetLineStyle(2);
        double intersect_1 = myFunc1->GetParameter(0);
        double intersect_1_error = myFunc1->GetParError(0);
        double slope_1 = myFunc1->GetParameter(1);
        double slope_1_error = myFunc1->GetParError(1);
        myFunc1->SetRange(0,1);

        f_cme_Pb->SetPoint(mult, pPb_ntrkBinCenter[mult], intersect_1);
        f_cme_Pb->SetPointError(mult, 0, intersect_1_error);


        gr_pPb_p[mult]->Fit("pol1");

        TF1 * myFunc2 = gr_pPb_p[mult]->GetFunction("pol1");
        myFunc2->SetLineColor(kBlue);
        myFunc2->SetLineStyle(2);
        double intersect_2 = myFunc2->GetParameter(0);
        double intersect_2_error = myFunc2->GetParError(0);
        double slope_2 = myFunc2->GetParameter(1);
        double slope_2_error = myFunc2->GetParError(1);
        myFunc2->SetRange(0,1);
    
        f_cme_p->SetPoint(mult, pPb_ntrkBinCenter[mult], intersect_2);
        f_cme_p->SetPointError(mult, 0, intersect_2_error);

    }   

    TCanvas* c1 = new TCanvas("c1","c1",700,700);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    gStyle->SetPadBorderMode(0.1);
    gStyle->SetOptTitle(0);

    TH1D* base2 = makeHist("base2", "Pb-going", "N^{offline}_{trk}", "f^{ cme}", 5000,0,5000,kBlack);
    base2->GetYaxis()->SetRangeUser(-0.0009, 0.0009);
    base2->GetXaxis()->SetRangeUser(81, 1700);
    base2->GetXaxis()->SetTitleColor(kBlack);
    
    fixedFontHist1D(base2,1.1,1.25);

    base2->GetYaxis()->SetTitleOffset(1.3);
    base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.6);
    base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.6);
    base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.6);
    base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.6);
    base2->GetXaxis()->SetNdivisions(4,6,0);
    base2->GetYaxis()->SetNdivisions(4,6,0);

    base2->GetXaxis()->SetLabelOffset(999);
    base2->GetXaxis()->SetTickLength(0);
    
    TGaxis *newaxis2 = new TGaxis(81,
                                -0.0009,
                                1700,
                                -0.0009,
                                81,
                                1700,
                                510,"G");
    newaxis2->SetLabelOffset(0.01);
    newaxis2->SetLabelFont(42);
    newaxis2->SetLabelSize(newaxis2->GetLabelSize()*1.0);
    
    base2->Draw();
    newaxis2->Draw("SS");
    gPad->Update();
    gPad->SetLogx(1);

    f_cme_Pb->SetMarkerStyle(20);
    f_cme_Pb->SetMarkerSize(1.5);
    f_cme_Pb->SetLineColor(kRed);
    f_cme_Pb->SetMarkerColor(kRed);
    f_cme_Pb->Draw("Psame");

    f_cme_p->SetMarkerStyle(21);
    f_cme_p->SetMarkerSize(1.5);
    f_cme_p->SetLineColor(kBlue);
    f_cme_p->SetMarkerColor(kBlue);
    //f_cme_p->Draw("Psame");
    

    TLegend *w2 = new TLegend(0.3,0.17,0.5,0.27);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(22);
    w2->SetTextFont(43);
    w2->AddEntry(f_cme_Pb, "pPb 8.16 TeV, #phi_{c}(Pb-going)", "P");
    w2->AddEntry(f_cme_p, "pPb 8.16 TeV, #phi_{c}(p-going)", "P");
    w2->Draw("same");

    TLatex* r11 = new TLatex(0.79,0.84, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

    TLatex* r3 = new TLatex(0.69, 0.94, "PbPb centrality(%)");
    r3->SetNDC();
    r3->SetTextSize(20);
    r3->SetTextFont(43);
    r3->SetTextColor(kBlack);
    r3->Draw("same");

    TLatex* cent1[7];
    cent1[0] = new TLatex(0.53, 0.91, "55");
    cent1[1] = new TLatex(0.67, 0.91, "45");
    cent1[2] = new TLatex(0.78, 0.91, "35");
    cent1[3] = new TLatex(0.35, 0.91, "65");
    cent1[4] = new TLatex(0.75, 0.91, "15");
    cent1[5] = new TLatex(0.80, 0.91, "7.5");
    cent1[6] = new TLatex(0.84, 0.91, "2.5");

    for(int i = 0; i < 4; i++){
        cent1[i]->SetNDC();
        cent1[i]->SetTextSize(20);
        cent1[i]->SetTextFont(43);
        cent1[i]->SetTextColor(kBlack);
        cent1[i]->Draw("same");
    }

    TLine* l1[7];
    l1[0] = new TLine(404.1,0.0009-0.00003, 404.1, 0.0009);
    l1[0]->SetLineWidth(2);
    l1[0]->Draw("Lsame");

    l1[1] = new TLine(717.6,0.0009-0.00003, 717.6, 0.0009);
    l1[1]->SetLineWidth(2);
    l1[1]->Draw("Lsame");

    l1[2] = new TLine(1141,0.0009-0.00003, 1141, 0.0009);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    l1[3] = new TLine(81.412,0.0009-0.00003, 81.412, 0.0009);
    l1[3]->SetLineWidth(2);
    //l1[3]->Draw("Lsame");

    l1[4] = new TLine(197,0.0009-0.00003, 197, 0.0009);
    l1[4]->SetLineWidth(2);
    l1[4]->Draw("Lsame");

    l1[5] = new TLine(3577.6,0.0009-0.00003, 3577.6, 0.0009);
    l1[5]->SetLineWidth(2);
    //l1[5]->Draw("Lsame");

    l1[6] = new TLine(4474,0.0009-0.00003, 4474, 0.0009);
    l1[6]->SetLineWidth(2);
    //l1[6]->Draw("Lsame");

    TLatex* r33 = new TLatex(0.18, 0.78, "|#eta_{#alpha} - #eta_{#beta}| < 1.6");
    r33->SetNDC();
    r33->SetTextSize(23);
    r33->SetTextFont(43);
    r33->SetTextColor(kBlack);
    r33->Draw("same");

    TLatex* r48 = new TLatex(0.18, 0.84, "0.3 < p_{T} < 3.0 GeV/c");
    r48->SetNDC();
    r48->SetTextSize(23);
    r48->SetTextFont(43);
    r48->SetTextColor(kBlack);
    r48->Draw("same");

   //  c1->Print("../figures/CME_figure5.pdf");
   //  c3->Print("../figures/CME_figure6.pdf");


}