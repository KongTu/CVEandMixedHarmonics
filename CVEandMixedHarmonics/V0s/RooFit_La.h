#include <iostream>
#include <fstream>

using namespace RooFit;

const bool Do_IntegrateFunc = true;

double ksptbin[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,5.0};//28 bins
double laptbin[] = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,5.0};

std::vector<double> N_sig_la[8][5];// 4 mult bin
std::vector<double> N_err_la[8][5];// error from fitting
std::vector<double> N_StatE_la[8][5];//error sqrt(sig+2*bkg)
std::vector<double> N_gen_la[8][5];//number of gen 
double Fit_CovQual;

void InitialFit_la(TCanvas* c,TH1* h, int ipt ,int j,int imult, vector<double>& FitParameter)
{
  cout<<"##########Fitting Bin j #"<<j<<" Pt Bin #"<< ipt<<endl;

  RooMsgService::instance().Print();
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  
  RooRealVar x("x","x",1.0,1.2);
  RooDataHist data("data","dataset with x",x,h);

  //-----double gaussian signal-----
  RooRealVar mean("mean","mean",1.1156,1.114,1.117);

  RooRealVar sigma1("sigma1","sigma1",0.002,0.001,0.002);//0.001, 0.0095
  RooRealVar sigma2("sigma2","sigma2",0.002,0.001,0.007);

  RooGaussian sig1("sig1","sig1",x,mean,sigma1); 
  RooGaussian sig2("sig2","sig2",x,mean,sigma2);

  string fitFunc1 = "A*pow(x-1.077,B)";
  string fitFunc2 = "A*pow(x-1.115,0.5)+B*pow(x-1.115,1.5)";

  RooClassFactory::makePdf("MyPdfV3","x,A,B","",fitFunc2.c_str(),kFALSE,kFALSE,"");

  #ifdef __CINT__
      gROOT->ProcessLineSync(".x MyPdfV3.cxx+") ;
  #endif

  //guess initial value of Nsig and Nbkg
  int peakbin1 = h->FindBin(1.09);
  int peakbin2 = h->FindBin(1.14);
  double Ntot_guess = h->Integral(peakbin1,peakbin2);
  int sidebin1 = h->FindBin(1.08); 
  int sidebin2 = h->FindBin(1.15);

  double Nside1_guess = h->Integral(sidebin1,peakbin1);
  double Nside2_guess = h->Integral(peakbin2,sidebin2);
  double Nbkgpeak_guess = 0.6*(Nside1_guess+Nside2_guess);
  double Nsig_guess = Ntot_guess-Nbkgpeak_guess;

  int startbin = h->FindBin(1.0);
  int endbin = h->FindBin(1.2);
  double Nfull = h->Integral(startbin,endbin);
  double Nbkg_guess = Nfull-Nsig_guess;


  //-----Add signal components-----
  RooRealVar sig1frac("sig1frac","sig1frac",0.3,0.,1);
  RooAddPdf sig("sig","sig",RooArgList(sig1,sig2),sig1frac);

  //-----number of signal/background counts-----
  RooRealVar Nsig("Nsig","Nsig",Nsig_guess,1000,1E9);
  RooRealVar Nbkg("Nbkg","Nbkg",Nbkg_guess*2,100,1E8);

  // Creat instance of MyPdfV3 class
  RooRealVar a0("a0","a0",0.1,-10,10) ;
  RooRealVar a1("a1","a1",0.1,-10,10) ;

  MyPdfV3 bkg("bkg","bkg",x,a0,a1) ;

  //-----Add signal and background-----
  RooAddPdf model("model","g1+g2+p",RooArgList(bkg,sig),RooArgList(Nbkg,Nsig));


  x.setRange("cut",1.097,1.13);

  model.fitTo(data,Minos(kTRUE),PrintLevel(-1),Range("cut"));

  RooFitResult *r = model.fitTo(data,Save(),Minos(kTRUE),PrintLevel(-1));

  int CovQuality = r->covQual();
  FitParameter.push_back(a0.getVal());
  FitParameter.push_back(a1.getVal());
  FitParameter.push_back(Nbkg.getVal());
  FitParameter.push_back(Nsig.getVal());
  FitParameter.push_back(mean.getVal());
  FitParameter.push_back(sig1frac.getVal());
  FitParameter.push_back(sigma1.getVal());
  FitParameter.push_back(sigma2.getVal());

  Fit_CovQual = CovQuality; 
  if(CovQuality < 3) return; 
  r->Print();

  RooPlot* frame = x.frame();
  data.plotOn(frame,RooFit::Name("data"),DrawOption("EX0"));
  model.plotOn(frame,RooFit::Name("model"),LineColor(kRed));
  model.plotOn(frame,Components(bkg),LineStyle(kDashed),LineColor(kBlack));// can add LineWidth(2.0) option

  frame->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c",laptbin[ipt],laptbin[ipt+1]));
  frame->SetXTitle("M(P#pi)[GeV]");
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->SetYTitle("Events(per 0.002GeV)");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->Draw();
  h->Draw("same");

  //get chi2
  double chi2 = frame->chiSquare("model","data",8);// the number of your free parameter

  double sigma1_val = sigma1.getVal();       
  double sigma2_val = sigma2.getVal();            
  double sig1frac_val = sig1frac.getVal();
  double ave_sigma = TMath::Sqrt(sig1frac_val*sigma1_val*sigma1_val+(1-sig1frac_val)*sigma2_val*sigma2_val);

  double mean_val = mean.getVal();
  
  double LowBand = mean_val-5*ave_sigma;
  double UpBand = mean_val+5*ave_sigma;

  double LowBand_InBin = h->FindBin(LowBand);// LowBand locates in which bin
  double UpBand_InBin = h->FindBin(UpBand);

  double bin_width = h->GetBinWidth(UpBand_InBin);

  double LowBand_lowedge = h->GetBinLowEdge(LowBand_InBin);
  double UpBand_lowedge = h->GetBinLowEdge(UpBand_InBin);
  double LowBand_upedge = LowBand_lowedge+bin_width;

  double LowBand_percent = (LowBand_upedge-LowBand)/bin_width;
  double UpBand_percent = (UpBand-UpBand_lowedge)/bin_width;
  
  double Counts_InLowBin = h->GetBinContent(LowBand_InBin)*LowBand_percent;
  double Counts_InUpBin = h->GetBinContent(UpBand_InBin)*UpBand_percent;

  int LowBin = LowBand_InBin + 1; // Bins for Integrate
  int UpBin = UpBand_InBin -1;

  //-----integrate to get background-----
  x.setRange("signal",LowBand,UpBand);
  RooAbsReal* fracsig = sig.createIntegral(x,x,"signal");
  double nsig_frac = Nsig.getVal()*fracsig->getVal();
 
  RooAbsReal* fracbkg = bkg.createIntegral(x,x,"signal");
  double nbkg_frac = Nbkg.getVal()*fracbkg->getVal();
  
  //-----integrate histogram over 5 sigma to get total yield-----
  double FullCounts = h->Integral(LowBin,UpBin)+Counts_InLowBin+Counts_InUpBin;

  if(Do_IntegrateFunc) double SigCounts = nsig_frac;// integrate fitting function
  else double SigCounts = FullCounts - nbkg_frac; // integrate histogram  

  double Yield_Error = Nsig.getError();//error from fitting 
  double Stat_Error = sqrt(SigCounts+2*nbkg_frac);// sqrt(sig+2bkg)

  TLatex latex;
  latex.SetTextAlign(12);
  latex.SetTextSize(0.06);
  latex.SetNDC();
  // latex.DrawLatex(0.15,0.80,Form("mean = %.4f GeV/c",mean_val));
  // latex.DrawLatex(0.15,0.70,Form("ave_sigma = %.4f GeV/c",ave_sigma));
  // latex.DrawLatex(0.15,0.60,Form("chi^{2}/ndof = %.3f",chi2));
  // latex.DrawLatex(0.15,0.50,Form("CovQuality = %d",CovQuality));

  N_sig_la[imult][j].push_back(SigCounts);
  N_err_la[imult][j].push_back(Yield_Error);
  N_StatE_la[imult][j].push_back(Stat_Error);

}

void RecursiveFit_la(TCanvas* c,TH1* h, int ipt ,int j,int imult, vector<double>& FitParameter)
{
  cout<<"##########Fitting Bin j #"<<j<<" Pt Bin #"<< ipt<<endl;

  RooMsgService::instance().Print();
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  RooRealVar x("x","x",1.0,1.2);
  RooDataHist data("data","dataset with x",x,h);

  double a0_G = 0.1;
  double a1_G = 0.1;
  double Nbkg_G = FitParameter[2];
  double Nsig_G = FitParameter[3];
  double mean_G = FitParameter[4];
  double sig1_frac_G = 0.4; 
  double sigma1_G = FitParameter[6];
  double sigma2_G = FitParameter[7];

  //-----double gaussian signal-----
  RooRealVar mean("mean","mean",mean_G,1.114,1.117);

  RooRealVar sigma1("sigma1","sigma1",sigma1_G,0.001,0.005);//0.001,0.0095
  RooRealVar sigma2("sigma2","sigma2",sigma2_G,0.001,0.005);

  RooGaussian sig1("sig1","sig1",x,mean,sigma1); 
  RooGaussian sig2("sig2","sig2",x,mean,sigma2);   

  //-----Add signal components-----
  RooRealVar sig1frac("sig1frac","sig1frac",sig1_frac_G,0.,1);
  RooAddPdf sig("sig","sig",RooArgList(sig1,sig2),sig1frac);

  //-----number of signal/background counts-----
  RooRealVar Nsig("Nsig","Nsig",Nsig_G,0,1E7);
  RooRealVar Nbkg("Nbkg","Nbkg",Nbkg_G,0,1E7);

  // Creat instance of MyPdfV3 class
  RooRealVar a0("a0","a0",0.1,-1,1) ;
  RooRealVar a1("a1","a1",0.1,-1,1) ;

  MyPdfV3 bkg("bkg","bkg",x,a0,a1) ;

  //-----Add signal and background-----
  RooAddPdf model("model","g1+g2+p",RooArgList(bkg,sig),RooArgList(Nbkg,Nsig));
	  
  model.fitTo(data,Minos(kTRUE),PrintLevel(-1));

  RooFitResult *r = model.fitTo(data,Save(),Minos(kTRUE),PrintLevel(-1));

  FitParameter[0] = a0.getVal();
  FitParameter[1] = a1.getVal();
  FitParameter[2] = Nbkg.getVal();
  FitParameter[3] = Nsig.getVal();
  FitParameter[4] = mean.getVal();
  FitParameter[5] = sig1frac.getVal();
  FitParameter[6] = sigma1.getVal();
  FitParameter[7] = sigma2.getVal();

  int CovQuality = r->covQual();
  Fit_CovQual = CovQuality; 
  if(CovQuality < 3) return; 
  r->Print();

  RooPlot* frame = x.frame();
  data.plotOn(frame,RooFit::Name("data"),DrawOption("EX0"));
  model.plotOn(frame,RooFit::Name("model"),LineColor(kRed));
  model.plotOn(frame,Components(bkg),LineStyle(kDashed),LineColor(kBlack));// can add LineWidth(2.0) option

  frame->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c",laptbin[ipt],laptbin[ipt+1]));
  frame->SetXTitle("M(P#pi)[GeV]");
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->SetYTitle("Events(per 0.002GeV)");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->Draw();
  h->Draw("same");

  //get chi2
  double chi2 = frame->chiSquare("model","data",8);// the number of your free parameter

  double sigma1_val = sigma1.getVal();       
  double sigma2_val = sigma2.getVal();            
  double sig1frac_val = sig1frac.getVal();
  double ave_sigma = TMath::Sqrt(sig1frac_val*sigma1_val*sigma1_val+(1-sig1frac_val)*sigma2_val*sigma2_val);

  double mean_val = mean.getVal();
  
  double LowBand = mean_val-5*ave_sigma;
  double UpBand = mean_val+5*ave_sigma;

  double LowBand_InBin = h->FindBin(LowBand);// LowBand locates in which bin
  double UpBand_InBin = h->FindBin(UpBand);

  double bin_width = h->GetBinWidth(UpBand_InBin);

  double LowBand_lowedge = h->GetBinLowEdge(LowBand_InBin);
  double UpBand_lowedge = h->GetBinLowEdge(UpBand_InBin);
  double LowBand_upedge = LowBand_lowedge+bin_width;

  double LowBand_percent = (LowBand_upedge-LowBand)/bin_width;
  double UpBand_percent = (UpBand-UpBand_lowedge)/bin_width;
  
  double Counts_InLowBin = h->GetBinContent(LowBand_InBin)*LowBand_percent;
  double Counts_InUpBin = h->GetBinContent(UpBand_InBin)*UpBand_percent;

  int LowBin = LowBand_InBin + 1; // Bins for Integrate
  int UpBin = UpBand_InBin -1;

  //-----integrate to get background-----
  x.setRange("signal",LowBand,UpBand);
  RooAbsReal* fracsig = sig.createIntegral(x,x,"signal");
  double nsig_frac = Nsig.getVal()*fracsig->getVal();
 
  RooAbsReal* fracbkg = bkg.createIntegral(x,x,"signal");
  double nbkg_frac = Nbkg.getVal()*fracbkg->getVal();
  
  //-----integrate histogram over 5 sigma to get total yield-----
  double FullCounts = h->Integral(LowBin,UpBin)+Counts_InLowBin+Counts_InUpBin;

  if(Do_IntegrateFunc) double SigCounts = nsig_frac;// integrate fitting function
  else double SigCounts = FullCounts - nbkg_frac; // integrate histogram  

  double Yield_Error = Nsig.getError();//error from fitting 
  double Stat_Error = sqrt(SigCounts+2*nbkg_frac);// sqrt(sig+2bkg)

  TLatex latex;
  latex.SetTextAlign(12);
  latex.SetTextSize(0.06);
  latex.SetNDC();
  latex.DrawLatex(0.15,0.80,Form("mean = %.4f GeV/c",mean_val));
  latex.DrawLatex(0.15,0.70,Form("ave_sigma = %.4f GeV/c",ave_sigma));
  latex.DrawLatex(0.15,0.60,Form("chi^{2}/ndof = %.3f",chi2));
  latex.DrawLatex(0.15,0.50,Form("CovQuality = %d",CovQuality));

  N_sig_la[imult][j].push_back(SigCounts);
  N_err_la[imult][j].push_back(Yield_Error);
  N_StatE_la[imult][j].push_back(Stat_Error);
  
}

