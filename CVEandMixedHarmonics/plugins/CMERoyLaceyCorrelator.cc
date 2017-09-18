// -*- C++ -*-
//
// Package:    CMERoyLaceyCorrelator/CMERoyLaceyCorrelator
// Class:      CMERoyLaceyCorrelator
// 
/**\class CMERoyLaceyCorrelator CMERoyLaceyCorrelator.cc CMERoyLaceyCorrelator/CMERoyLaceyCorrelator/plugins/CMERoyLaceyCorrelator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Mon, 01 Aug 2016 09:01:02 GMT
//
//


#include "CVEandMixedHarmonics/CVEandMixedHarmonics/interface/CMERoyLaceyCorrelatorBase.h"


CMERoyLaceyCorrelator::CMERoyLaceyCorrelator(const edm::ParameterSet& iConfig)
{

  trackName_  =  iConfig.getParameter<edm::InputTag>("trackName");
  vertexName_ =  iConfig.getParameter<edm::InputTag>("vertexName");
  towerName_ =  iConfig.getParameter<edm::InputTag>("towerName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);

  NsubSamples_ = iConfig.getUntrackedParameter<int>("NsubSamples");
  Nembedded_ = iConfig.getUntrackedParameter<int>("Nembedded");

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
  dopPb_ = iConfig.getUntrackedParameter<bool>("dopPb");
  doBothSide_ = iConfig.getUntrackedParameter<bool>("doBothSide");
  doLightWeight_ = iConfig.getUntrackedParameter<bool>("doLightWeight");

  eff_ = iConfig.getUntrackedParameter<int>("eff");

  etaTracker_ = iConfig.getUntrackedParameter<double>("etaTracker");
  gapValue_ = iConfig.getUntrackedParameter<double>("gapValue");
  
  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");
  
  etaLowQ2_ = iConfig.getUntrackedParameter<double>("etaLowQ2");
  etaHighQ2_ = iConfig.getUntrackedParameter<double>("etaHighQ2");
   
  vzLow_ = iConfig.getUntrackedParameter<double>("vzLow");
  vzHigh_ = iConfig.getUntrackedParameter<double>("vzHigh");
  
  ptLow_ = iConfig.getUntrackedParameter<double>("ptLow");
  ptHigh_ = iConfig.getUntrackedParameter<double>("ptHigh");

  q2max_ = iConfig.getUntrackedParameter<double>("q2max");
  q2min_ = iConfig.getUntrackedParameter<double>("q2min");

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);
  offlineChi2_ = iConfig.getUntrackedParameter<double>("offlineChi2", 0.0);
  offlinenhits_ = iConfig.getUntrackedParameter<double>("offlinenhits", 0.0);

  etaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
  dEtaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("dEtaBins");
  ptBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ptBins");
  dPtBins_ = iConfig.getUntrackedParameter<std::vector<double>>("dPtBins");  
  centBins_ = iConfig.getUntrackedParameter<std::vector<double>>("centBins");

}


CMERoyLaceyCorrelator::~CMERoyLaceyCorrelator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CMERoyLaceyCorrelator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

/*
q2 calculation at HF and selections:
*/

  double qHFcos = 0.;
  double qHFsin = 0.;
  double qHF_count = 0.;
  for(unsigned i = 0; i < towers->size(); ++i){

          const CaloTower & hit= (*towers)[i];

          double caloEta = hit.eta();
          double caloPhi = hit.phi();
          double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );
          if( reverseBeam_ ) caloEta = -hit.eta();

    if( dopPb_ && !doBothSide_ ) {if( caloEta < etaLowQ2_ || caloEta > etaHighQ2_ ) continue;}
    else{ if( caloEta < etaLowQ2_ || caloEta > etaHighQ2_ ) continue; }

    qHFcos += w*cos(-n3_*caloPhi);
    qHFsin += w*sin(-n3_*caloPhi);
    qHF_count += w;

  }

  double q2HF_real = qHFcos/qHF_count;
  double q2HF_imag = qHFsin/qHF_count;
  double magnitude_HF = sqrt(q2HF_imag*q2HF_imag + q2HF_real*q2HF_real);
  double Psi_RP_HF = TMath::ATan(q2HF_imag/q2HF_real);

  if( magnitude_HF > q2max_ || magnitude_HF < q2min_ ) return;//q2 selections. 
  
  q2_mag->Fill( magnitude_HF );
  Ntrk_q2->Fill(nTracks);


  double s_hp = 0.0;
  double s_hn = 0.0;
  double s_mp = 0.0;
  double s_mn = 0.0;

  double Np = 0.0;
  double Nn = 0.0;
  double Nmp = 0.0;
  double Nmn = 0.0;

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
    if( doEffCorrection_ ) { 

      if( dopPb_ ){
        weight = 1.0/effTable_pPb[eff_]->GetBinContent( effTable_pPb[eff_]->FindBin(trk.eta(), trk.pt()) );
      }
      else{
        weight = 1.0/effTable[eff_]->GetBinContent( effTable[eff_]->FindBin(trk.eta(), trk.pt()) );
      }

    }

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

    if( trk.charge() == 1 ){
      s_hp += weight*sin( phi - Psi_RP_HF );
      Np += weight;
    }
    if( trk.charge() == -1 ){
      s_hn += weight*sin( phi - Psi_RP_HF );
      Nn += weight;
    }

    unsigned int random_charge;
    random_charge = gRandom->Integer(2);

    if( random_charge == 0 ){
      s_mp += weight*sin( phi - Psi_RP_HF );
      Nmp++;
    }
    if( random_charge == 1 ){
      s_mn += weight*sin( phi - Psi_RP_HF );
      Nmn++;
    }

  }


  double numerator = s_hp/Np - s_hn/Nn;
  double denominator = s_mp/Nmp - s_mn/Nmn;

  CcS->Fill(numerator/denominator, Np+Nn);


}
// ------------ method called once each job just before starting event loop  ------------
void 
CMERoyLaceyCorrelator::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  double ptBinsArray[100];
  const int Nptbins = ptBins_.size() - 1;
  for(unsigned i = 0; i < ptBins_.size(); i++){
    ptBinsArray[i] = ptBins_[i];
  }

  edm::FileInPath fip1("CVEandMixedHarmonics/CVEandMixedHarmonics/data/Hydjet_eff_mult_v1.root");
  TFile f1(fip1.fullPath().c_str(),"READ");
  for(int i = 0; i < 5; i++){
     effTable[i] = (TH2D*)f1.Get(Form("rTotalEff3D_%d",i));
  }

  edm::FileInPath fip2("CVEandMixedHarmonics/CVEandMixedHarmonics/data/Hijing_8TeV_CME_dataBS.root");
  TFile f2(fip2.fullPath().c_str(),"READ");
  for(int i = 0; i < 1; i++){
     effTable_pPb[i] = (TH2D*)f2.Get(Form("rTotalEff3D_%d",i));
  }
  
  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  vtxZ = fs->make<TH1D>("vtxZ",";vz", 400,-20,20);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);
  trkPt = fs->make<TH1D>("trkPt", ";p_{T}(GeV)", Nptbins,ptBinsArray);
  trk_eta = fs->make<TH1D>("trk_eta", ";#eta", 50,-2.5,2.5);
  q2_mag = fs->make<TH1D>("q2_mag", "q2", 2000,-1,1);
  Ntrk_q2 = fs->make<TH1D>("Ntrk_q2",";Ntrk",5000,0,5000);

  CcS = fs->make<TH1D>("CcS", ";#DeltaS", 2000,-2,2);
}
vector<double> 
CMERoyLaceyCorrelator::get4Momentum(double pt, double eta, double phi, double mass)
{
  double polar_angle = 2*TMath::ATan( TMath::Exp(-eta) );
  double pz = pt/TMath::Tan( polar_angle );
  double px = sqrt(pt*pt/( 1+TMath::Tan(phi)*TMath::Tan(phi) ) );
  double py = sqrt(pt*pt - px*px);
  double E = sqrt(px*px+py*py+pz*pz + mass*mass);

  vector<double> temp;

  if( phi > 0 && phi < PI){
    py = py; 
  }
  else if( phi < 0 && phi > -PI){
    py = -py;
  }

  if( phi > -PI/2.0 && phi < PI/2.0 ){
    px = px;
  }
  else if( phi < -PI/2.0 || phi > PI/2.0 ){
    px = -px;
  }
  
  temp.push_back( E );
  temp.push_back( px );
  temp.push_back( py );
  temp.push_back( pz );
  temp.push_back( polar_angle ); 

  return temp;

}
vector<double> 
CMERoyLaceyCorrelator::getLightConeVar(double px, double py, double pz){

  double pt = sqrt(px*px + py*py);
  double phi = TMath::ATan(py/px);
  double three_momentum = sqrt(px*px+py*py+pz*pz);
  double eta = TMath::ATanH( pz/three_momentum );

  vector<double> temp;

  if( px > 0 && py > 0){
    phi = phi;
  }
  else if( px > 0 && py < 0){
    phi = phi;
  }
  else if( px < 0 && py > 0){
    phi = PI + phi;
  }
  else if( px < 0 && py < 0){
    phi = -PI + phi;
  }

  temp.push_back( pt );
  temp.push_back( eta );
  temp.push_back( phi );

  return temp; 
}
TComplex 
CMERoyLaceyCorrelator::q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}
// ------------ method called once each job just after ending the event loop  ------------
void 
CMERoyLaceyCorrelator::endJob() 
{
}
void 
CMERoyLaceyCorrelator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CMERoyLaceyCorrelator::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CMERoyLaceyCorrelator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CMERoyLaceyCorrelator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CMERoyLaceyCorrelator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CMERoyLaceyCorrelator);
