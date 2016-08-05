// -*- C++ -*-
//
// Package:    CVEandMixedHarmonics/CVEandMixedHarmonics
// Class:      CVEandMixedHarmonics
// 
/**\class CVEandMixedHarmonics CVEandMixedHarmonics.cc CVEandMixedHarmonics/CVEandMixedHarmonics/plugins/CVEandMixedHarmonics.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Mon, 01 Aug 2016 09:01:02 GMT
//
//


#include "CVEandMixedHarmonics/CVEandMixedHarmonics/interface/CVEandMixedHarmonicsBase.h"


CVEandMixedHarmonics::CVEandMixedHarmonics(const edm::ParameterSet& iConfig)
{

  trackName_  =  iConfig.getParameter<edm::InputTag>("trackName");
  vertexName_ =  iConfig.getParameter<edm::InputTag>("vertexName");
  towerName_ =  iConfig.getParameter<edm::InputTag>("towerName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);

  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");

  n1_ = iConfig.getUntrackedParameter<int>("n1");
  n2_ = iConfig.getUntrackedParameter<int>("n2");
  n3_ = iConfig.getUntrackedParameter<int>("n3");
  n4_ = iConfig.getUntrackedParameter<int>("n4");
  
  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality");
  reverseBeam_ = iConfig.getUntrackedParameter<bool>("reverseBeam");
  doEffCorrection_ = iConfig.getUntrackedParameter<bool>("doEffCorrection");
  useEtaGap_ = iConfig.getUntrackedParameter<bool>("useEtaGap");

  eff_ = iConfig.getUntrackedParameter<int>("eff");

  etaTracker_ = iConfig.getUntrackedParameter<double>("etaTracker");

  gapValue_ = iConfig.getUntrackedParameter<double>("gapValue");
  
  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");
  
  vzLow_ = iConfig.getUntrackedParameter<double>("vzLow");
  vzHigh_ = iConfig.getUntrackedParameter<double>("vzHigh");
  
  ptLow_ = iConfig.getUntrackedParameter<double>("ptLow");
  ptHigh_ = iConfig.getUntrackedParameter<double>("ptHigh");

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);
  offlineChi2_ = iConfig.getUntrackedParameter<double>("offlineChi2", 0.0);
  offlinenhits_ = iConfig.getUntrackedParameter<double>("offlinenhits", 0.0);
  
  etaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
  dEtaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("dEtaBins");
  ptBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ptBins");
  centBins_ = iConfig.getUntrackedParameter<std::vector<double>>("centBins");

}


CVEandMixedHarmonics::~CVEandMixedHarmonics()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CVEandMixedHarmonics::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); 
  bestvx = vtx.x(); 
  bestvy = vtx.y();
  bestvzError = vtx.zError(); 
  bestvxError = vtx.xError(); 
  bestvyError = vtx.yError();

  //first selection; vertices
  if( fabs(bestvz) < vzLow_ || fabs(bestvz) > vzHigh_ ) return;

  vtxZ->Fill( bestvz );

  Handle<CaloTowerCollection> towers;
  iEvent.getByToken(towerSrc_, towers);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackSrc_, tracks);

  int nTracks = 0;
  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];
  
     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError); 

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > 0.1 ) continue;
        if(fabs(dzvtx/dzerror) > 3.0) continue;
        if(fabs(dxyvtx/dxyerror) > 3.0) continue;
        if(trk.pt() < 0.4 || fabs(trk.eta()) > 2.4) continue;
        nTracks++;//count multiplicity

  }

  if( !useCentrality_ ) if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;

  double etHFtowerSumPlus = 0.0;
  double etHFtowerSumMinus = 0.0;
  double etHFtowerSum = 0.0;
  
  if( useCentrality_ ){

    for( unsigned i = 0; i<towers->size(); ++ i){
       const CaloTower & tower = (*towers)[ i ];
       double eta = tower.eta();
       bool isHF = tower.ietaAbs() > 29;
          if(isHF && eta > 0){
            etHFtowerSumPlus += tower.pt();
          }
          if(isHF && eta < 0){
            etHFtowerSumMinus += tower.pt();
          }
    }
    etHFtowerSum=etHFtowerSumPlus + etHFtowerSumMinus;

    int bin = -1;
    for(int j=0; j<200; j++){
      if( etHFtowerSum >= centBins_[j] ){
         bin = j; break;
      }
    }

    int hiBin = bin;
    if( hiBin < Nmin_ || hiBin >= Nmax_ ) return;
    cbinHist->Fill( hiBin );

  }

  Ntrk->Fill( nTracks );

  const int NetaBins = etaBins_.size() - 1 ;
  const int NdEtaBins = dEtaBins_.size() - 1;
  double dEtaBinsArray[100];

  for(unsigned i = 0; i < dEtaBins_.size(); i++){

    dEtaBinsArray[i] = dEtaBins_[i]-0.0001;
  }

/*
define all the ingredients for CME correlator <cos(n1*phi_1 + n2*phi_2 + n3*phi_3)> ,
where Q_coefficient_power is used in the following names 
 */

//2-particle correlator

  TComplex Q_n1_1[NetaBins][2], Q_n2_1[NetaBins][2];
  
  TComplex Q_n1n2_2[NetaBins][2];
  
  TComplex Q_0_1[NetaBins][2], Q_0_2[NetaBins][2];

//------------------------------------------------------------------

//Start filling Q-vectors;

  //track loop to fill charged particles Q-vectors
  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];
  
     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError); 
        double nhits = trk.numberOfValidHits();
        double chi2n = trk.normalizedChi2();
        double nlayers = trk.hitPattern().trackerLayersWithMeasurement();
        chi2n = chi2n/nlayers;
        double nPixelLayers = trk.hitPattern().pixelLayersWithMeasurement();//only pixel layers
        double phi = trk.phi();
        double trkEta = trk.eta();

        double weight = 1.0;
        if( doEffCorrection_ ) { weight = 1.0/effTable[eff_]->GetBinContent( effTable[eff_]->FindBin(trk.eta(), trk.pt()) );}

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
        if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
        if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
        if(chi2n > offlineChi2_ ) continue;
        if(nhits < offlinenhits_ ) continue;
        if( nPixelLayers <= 0 ) continue;
        if(trk.pt() < ptLow_ || trk.pt() > ptHigh_ ) continue;
        if(fabs(trkEta) > etaTracker_ ) continue;

        trkPhi->Fill(phi, weight);
        trkPt->Fill(trk.pt(), weight);
        trk_eta->Fill(trkEta, weight);

        for(int eta = 0; eta < NetaBins; eta++){
          if( trkEta > etaBins_[eta] && trkEta < etaBins_[eta+1] ){

            if( trk.charge() == +1 ){

              Q_n1_1[eta][0] += q_vector(n1_, 1, weight, phi);
              Q_n2_1[eta][0] += q_vector(n2_, 1, weight, phi);

              Q_n1n2_2[eta][0] += q_vector(n1_+n2_, 2, weight, phi);

              Q_0_1[eta][0] += q_vector(0, 1, weight, phi);
              Q_0_2[eta][0] += q_vector(0, 2, weight, phi);

            }
            if( trk.charge() == -1 ){

              Q_n1_1[eta][1] += q_vector(n1_, 1, weight, phi);
              Q_n2_1[eta][1] += q_vector(n2_, 1, weight, phi);

              Q_n1n2_2[eta][1] += q_vector(n1_+n2_, 2, weight, phi);

              Q_0_1[eta][1] += q_vector(0, 1, weight, phi);
              Q_0_2[eta][1] += q_vector(0, 2, weight, phi);

            }
          }
        }
  }

  TComplex  Q_n3_1_HFplus, Q_n3_1_HFminus, Q_0_1_HFplus, Q_0_1_HFminus;
  //HF towers loop to fill the towers' Q-vectors:
  for(unsigned i = 0; i < towers->size(); ++i){

          const CaloTower & hit= (*towers)[i];

          double caloEta = hit.eta();
          double caloPhi = hit.phi();
          double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );

          hfPhi->Fill(caloPhi, w);
          
          if( reverseBeam_ ) caloEta = -hit.eta();          
          if( caloEta < etaHighHF_ && caloEta > etaLowHF_ ){
            
              Q_n3_1_HFplus += q_vector(n3_, 1, w, caloPhi);
              Q_0_1_HFplus += q_vector(0, 1, w, caloPhi);

          }
          else if( caloEta < -etaLowHF_ && caloEta > -etaHighHF_ ){

              Q_n3_1_HFminus += q_vector(n3_, 1, w, caloPhi);
              Q_0_1_HFminus += q_vector(0, 1, w, caloPhi); 

          }
          else{continue;}
  }
    
/*
calculate the 3-particles correlator with the charged-particles
 */

  for(int ieta = 0; ieta < NetaBins; ieta++){
    for(int jeta = 0; jeta < NetaBins; jeta++){

      double deltaEta = fabs( etaBins_[ieta] - etaBins_[jeta] );

      for(int deta = 0; deta < NdEtaBins; deta++){
        if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta+1] ){
          
          TComplex N_2;
          TComplex D_2;

          TComplex N_3_HFplus, N_3_HFminus;
          TComplex D_3_HFplus, D_3_HFminus;

          //same sign correlator:
          for(int sign = 0; sign < 2; sign++){
            if( ieta == jeta ){

              N_2 = Q_n1_1[ieta][sign]*Q_n2_1[ieta][sign] - Q_n1n2_2[ieta][sign];
              D_2 = Q_0_1[ieta][sign]*Q_0_1[ieta][sign] - Q_0_2[ieta][sign];

            }
            else{
              
              N_2 = Q_n1_1[ieta][sign]*Q_n2_1[jeta][sign];
              D_2 = Q_0_1[ieta][sign]*Q_0_1[jeta][sign];

            }        

            N_3_HFplus = N_2*Q_n3_1_HFplus;
            D_3_HFplus = D_2*Q_0_1_HFplus;

            c3_real[deta][sign][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re() );
            c3_imag[deta][sign][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re() );

            N_3_HFminus = N_2*Q_n3_1_HFminus;
            D_3_HFminus = D_2*Q_0_1_HFminus;

            c3_real[deta][sign][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re() );
            c3_imag[deta][sign][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re() );
  
          }
          //opposite sign correlator:
          N_2 = Q_n1_1[ieta][0]*Q_n2_1[jeta][1];
          D_2 = Q_0_1[ieta][0]*Q_0_1[jeta][1];

          N_3_HFplus = N_2*Q_n3_1_HFplus;
          D_3_HFplus = D_2*Q_0_1_HFplus;

          c3_real[deta][2][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re() );
          c3_imag[deta][2][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re() );

          N_3_HFminus = N_2*Q_n3_1_HFminus;
          D_3_HFminus = D_2*Q_0_1_HFminus;

          c3_real[deta][2][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re() );
          c3_imag[deta][2][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re() );

        }
      }
    }
  }

}
// ------------ method called once each job just before starting event loop  ------------
void 
CVEandMixedHarmonics::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  const int NdEtaBins = dEtaBins_.size() - 1;
  const int NetaBins = etaBins_.size() - 1;
  double etaBinsArray[100];
  for(unsigned i = 0; i < etaBins_.size(); i++){
    etaBinsArray[i] = etaBins_[i];
  }
  const int Nptbins = ptBins_.size() - 1;
  double ptBinsArray[100];
  for(unsigned i = 0; i < ptBins_.size(); i++){
    ptBinsArray[i] = ptBins_[i];
  }

  edm::FileInPath fip1("CVEandMixedHarmonics/CVEandMixedHarmonics/data/Hydjet_eff_mult_v1.root");
  TFile f1(fip1.fullPath().c_str(),"READ");
  for(int i = 0; i < 5; i++){
     effTable[i] = (TH2D*)f1.Get(Form("rTotalEff3D_%d",i));
  }

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  vtxZ = fs->make<TH1D>("vtxZ",";vz", 400,-20,20);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);
  trkPt = fs->make<TH1D>("trkPt", ";p_{T}(GeV)", Nptbins,ptBinsArray);
  trk_eta = fs->make<TH1D>("trk_eta", ";#eta", NetaBins, etaBinsArray);

  for(int deta = 0; deta < NdEtaBins; deta++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < 2; HF++){

        c3_real[deta][sign][HF] = fs->make<TH1D>(Form("c3_real_%d_%d_%d", deta, sign, HF),";c3", 20000,-1,1);
        c3_imag[deta][sign][HF] = fs->make<TH1D>(Form("c3_imag_%d_%d_%d", deta, sign, HF),";c3", 20000,-1,1);

      }
    }    
  }
}

TComplex 
CVEandMixedHarmonics::q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CVEandMixedHarmonics::endJob() 
{
}
void 
CVEandMixedHarmonics::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CVEandMixedHarmonics::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CVEandMixedHarmonics::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CVEandMixedHarmonics::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CVEandMixedHarmonics::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CVEandMixedHarmonics);
