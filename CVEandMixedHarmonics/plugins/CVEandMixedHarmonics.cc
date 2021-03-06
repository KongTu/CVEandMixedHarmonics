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
  generalV0_ksName_  =  iConfig.getParameter<edm::InputTag>("generalV0_ksName");
  generalV0_laName_  =  iConfig.getParameter<edm::InputTag>("generalV0_laName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);
  generalV0_ks_ = consumes<reco::VertexCompositeCandidateCollection>(generalV0_ksName_);
  generalV0_la_ = consumes<reco::VertexCompositeCandidateCollection>(generalV0_laName_);

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

  v0EtaTracker_ = iConfig.getUntrackedParameter<double>("v0EtaTracker",2.4);
  v0sNhitsCut_ = iConfig.getUntrackedParameter<double>("v0sNhitsCut",3.0);
  
  K0s_dcaCut_ = iConfig.getUntrackedParameter<double>("K0s_dcaCut",1.0);
  K0s_decayLengthCut_ = iConfig.getUntrackedParameter<double>("K0s_decayLengthCut",5.0);
  K0s_pointingAngleCut_ = iConfig.getUntrackedParameter<double>("K0s_pointingAngleCut",0.999);

  Lam_dcaCut_ = iConfig.getUntrackedParameter<double>("Lam_dcaCut",1.0);
  Lam_decayLengthCut_ = iConfig.getUntrackedParameter<double>("Lam_decayLengthCut",5.0);
  Lam_pointingAngleCut_ = iConfig.getUntrackedParameter<double>("Lam_pointingAngleCut",0.999);

  lambdaMassWindow_ = iConfig.getUntrackedParameter<double>("lambdaMassWindow",0.01);
  ksMassWindow_ = iConfig.getUntrackedParameter<double>("ksMassWindow",0.02);
  
  K0sPtLow_ = iConfig.getUntrackedParameter<double>("K0sPtLow");
  K0sPtHigh_ = iConfig.getUntrackedParameter<double>("K0sPtHigh");
  LamPtLow_ = iConfig.getUntrackedParameter<double>("LamPtLow");
  LamPtHigh_ = iConfig.getUntrackedParameter<double>("LamPtHigh");

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

  Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
  iEvent.getByToken(generalV0_ks_,v0candidates_ks);
    
  Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
  iEvent.getByToken(generalV0_la_,v0candidates_la);

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

//2-particle correlator (charge dependent for 0 and 1, and charge independent for 2 in the second argument)

  TComplex Q_n1_1[NetaBins][2], Q_n2_1[NetaBins][2];
  
  TComplex Q_n1n2_2[NetaBins][2];
  
  TComplex Q_0_1[NetaBins][2], Q_0_2[NetaBins][2];


//Scalar product, charge independent, |eta|<2.4

  TComplex Q_n3_trk, Q_0_trk;

//3 particle to pair with V0s, charge independent

  TComplex Q_forV0s_trk, Q_0_forV0s_trk;

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

    Q_n3_trk += q_vector(-n3_, 1, weight, phi);//for scalar product in tracker
    Q_0_trk += q_vector(0, 1, weight, phi);

    Q_forV0s_trk += q_vector(n2_, 1, weight, phi);//for charged particles to be compared with V0s
    Q_0_forV0s_trk += q_vector(0, 1, weight, phi);

    for(int eta = 0; eta < NetaBins; eta++){
      if( trkEta > etaBins_[eta] && trkEta < etaBins_[eta+1] ){

        if( trk.charge() == +1 ){//positive charge

          Q_n1_1[eta][0] += q_vector(n1_, 1, weight, phi);
          Q_n2_1[eta][0] += q_vector(n2_, 1, weight, phi);

          Q_n1n2_2[eta][0] += q_vector(n1_+n2_, 2, weight, phi);

          Q_0_1[eta][0] += q_vector(0, 1, weight, phi);
          Q_0_2[eta][0] += q_vector(0, 2, weight, phi);

        }
        if( trk.charge() == -1 ){//negative charge

          Q_n1_1[eta][1] += q_vector(n1_, 1, weight, phi);
          Q_n2_1[eta][1] += q_vector(n2_, 1, weight, phi);

          Q_n1n2_2[eta][1] += q_vector(n1_+n2_, 2, weight, phi);

          Q_0_1[eta][1] += q_vector(0, 1, weight, phi);
          Q_0_2[eta][1] += q_vector(0, 2, weight, phi);

        }
      }
    }
  }

  //Declaring TComplex varaibles for Q-vectors of particle c (HF)

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

  //Declaring TComplex varaibles for Q-vectors of V0s
  
  TComplex Q_n1_1_K0s_sig, Q_n2_1_K0s_sig, Q_n1n2_2_K0s_sig, Q_0_1_K0s_sig, Q_0_2_K0s_sig;
  TComplex Q_n1_1_K0s_bkg, Q_n2_1_K0s_bkg, Q_n1n2_2_K0s_bkg, Q_0_1_K0s_bkg, Q_0_2_K0s_bkg;

  //lambda has particle [0] and anti-particle [1]
  TComplex Q_n1_1_Lambda_sig[2], Q_n2_1_Lambda_sig[2], Q_n1n2_2_Lambda_sig[2], Q_0_1_Lambda_sig[2], Q_0_2_Lambda_sig[2];
  TComplex Q_n1_1_Lambda_bkg[2], Q_n2_1_Lambda_bkg[2], Q_n1n2_2_Lambda_bkg[2], Q_0_1_Lambda_bkg[2], Q_0_2_Lambda_bkg[2];

  //filling Q-vectors of K0s
  for(unsigned it=0; it<v0candidates_ks->size(); ++it){ 

    const reco::VertexCompositeCandidate & V0s = (*v0candidates_ks)[it];
    bool isK0s = true;
    if( !passV0sCut(V0s, vtx, isK0s) ) continue;

    double weight = 1.0;
    double mass = V0s.mass();
    double phi = V0s.phi();
    double pt = V0s.pt();
    double eta = V0s.eta();

    ks_mass->Fill(eta, pt, mass, weight);

    if( mass < K0sMass+ksMassWindow_ && mass > K0sMass-ksMassWindow_ ){

      Q_n1_1_K0s_sig += q_vector(n1_, 1, weight, phi);
      Q_n2_1_K0s_sig += q_vector(n2_, 1, weight, phi);

      Q_n1n2_2_K0s_sig += q_vector(n1_+n2_, 2, weight, phi);

      Q_0_1_K0s_sig += q_vector(0, 1, weight, phi);
      Q_0_2_K0s_sig += q_vector(0, 2, weight, phi);

    }
    if( mass > K0sMass+ksMassWindow_ || mass < K0sMass-ksMassWindow_ ){

      Q_n1_1_K0s_bkg += q_vector(n1_, 1, weight, phi);
      Q_n2_1_K0s_bkg += q_vector(n2_, 1, weight, phi);

      Q_n1n2_2_K0s_bkg += q_vector(n1_+n2_, 2, weight, phi);

      Q_0_1_K0s_bkg += q_vector(0, 1, weight, phi);
      Q_0_2_K0s_bkg += q_vector(0, 2, weight, phi);

    }
  }

  //filling Q-vectors of Lambda
  for(unsigned it=0; it<v0candidates_la->size(); ++it){ 

    const reco::VertexCompositeCandidate & V0s = (*v0candidates_la)[it];
    bool isK0s = false;//isLambda
    if( !passV0sCut(V0s, vtx, isK0s) ) continue;

    double weight = 1.0;
    double mass = V0s.mass();
    double phi = V0s.phi();
    double pt = V0s.pt();
    double eta = V0s.eta();
    int pdgId = V0s.pdgId();

    if( pdgId == 3122 ){

      la_mass_1->Fill(eta, pt, mass, weight);

      if( mass < LambdaMass+lambdaMassWindow_ && mass > LambdaMass-lambdaMassWindow_ ){

      Q_n1_1_Lambda_sig[0] += q_vector(n1_, 1, weight, phi);
      Q_n2_1_Lambda_sig[0] += q_vector(n2_, 1, weight, phi);

      Q_n1n2_2_Lambda_sig[0] += q_vector(n1_+n2_, 2, weight, phi);

      Q_0_1_Lambda_sig[0] += q_vector(0, 1, weight, phi);
      Q_0_2_Lambda_sig[0] += q_vector(0, 2, weight, phi);

      }
      if( mass > LambdaMass+lambdaMassWindow_ || mass < LambdaMass-lambdaMassWindow_ ){

        Q_n1_1_Lambda_bkg[0] += q_vector(n1_, 1, weight, phi);
        Q_n2_1_Lambda_bkg[0] += q_vector(n2_, 1, weight, phi);

        Q_n1n2_2_Lambda_bkg[0] += q_vector(n1_+n2_, 2, weight, phi);

        Q_0_1_Lambda_bkg[0] += q_vector(0, 1, weight, phi);
        Q_0_2_Lambda_bkg[0] += q_vector(0, 2, weight, phi);

      }
    }  
    if( pdgId == -3122 ){

      la_mass_2->Fill(eta, pt, mass, weight);

      if( mass < LambdaMass+lambdaMassWindow_ && mass > LambdaMass-lambdaMassWindow_ ){

      Q_n1_1_Lambda_sig[1] += q_vector(n1_, 1, weight, phi);
      Q_n2_1_Lambda_sig[1] += q_vector(n2_, 1, weight, phi);

      Q_n1n2_2_Lambda_sig[1] += q_vector(n1_+n2_, 2, weight, phi);

      Q_0_1_Lambda_sig[1] += q_vector(0, 1, weight, phi);
      Q_0_2_Lambda_sig[1] += q_vector(0, 2, weight, phi);

      }
      if( mass > LambdaMass+lambdaMassWindow_ || mass < LambdaMass-lambdaMassWindow_ ){

        Q_n1_1_Lambda_bkg[1] += q_vector(n1_, 1, weight, phi);
        Q_n2_1_Lambda_bkg[1] += q_vector(n2_, 1, weight, phi);

        Q_n1n2_2_Lambda_bkg[1] += q_vector(n1_+n2_, 2, weight, phi);

        Q_0_1_Lambda_bkg[1] += q_vector(0, 1, weight, phi);
        Q_0_2_Lambda_bkg[1] += q_vector(0, 2, weight, phi);

      }
    } 
  }

/*
calculate the Scalar product denominator, v_{2,c}
*/

  TComplex N_2_trk, D_2_trk;

  N_2_trk = Q_n3_trk*Q_n3_1_HFplus;
  D_2_trk = Q_0_trk*Q_0_1_HFplus;

  c2_cb->Fill(N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

  N_2_trk = Q_n3_trk*Q_n3_1_HFminus;
  D_2_trk = Q_0_trk*Q_0_1_HFminus;

  c2_ac->Fill(N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

  N_2_trk = Q_n3_1_HFplus*TComplex::Conjugate(Q_n3_1_HFminus);
  D_2_trk = Q_0_1_HFplus*Q_0_1_HFminus;

  c2_ab->Fill(N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());
    
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

/*
calculate the 3-particle correlator with V0s, generally 5 cases. Baryon number separation correlator
 */
  
//temp complex number for each fill
TComplex N_3_HFplus, N_3_HFminus;
TComplex D_3_HFplus, D_3_HFminus;

//Lambda-Lambda
  //2-particles with V0s

  for(int sign = 0; sign < 2; sign++){

    TComplex N_2_sig_LL, D_2_sig_LL; //sig-sig correlator
    TComplex N_2_bkg_LL, D_2_bkg_LL; //bkg-bkg correlatro
    TComplex N_2_sig_bkg_LL, D_2_sig_bkg_LL;//sig-bkg correlator (note that the n1 and n2 would make a different from the one below)
    TComplex N_2_bkg_sig_LL, D_2_bkg_sig_LL;//bkg-sig correlator

    //--------------------------------------------
    N_2_sig_LL = Q_n1_1_Lambda_sig[sign]*Q_n2_1_Lambda_sig[sign] - Q_n1n2_2_Lambda_sig[sign];
    D_2_sig_LL = Q_0_1_Lambda_sig[sign]*Q_0_1_Lambda_sig[sign] - Q_0_2_Lambda_sig[sign];

    N_2_sig_bkg_LL = Q_n1_1_Lambda_sig[sign]*Q_n2_1_Lambda_bkg[sign];//sig-bkg doesn't overlap, no need to subtract the 3rd term
    D_2_sig_bkg_LL = Q_0_1_Lambda_sig[sign]*Q_0_1_Lambda_bkg[sign];

    N_2_bkg_sig_LL = Q_n1_1_Lambda_bkg[sign]*Q_n2_1_Lambda_sig[sign];//bkg-sig doesn't overlap, no need to subtract the 3rd term
    D_2_bkg_sig_LL = Q_0_1_Lambda_bkg[sign]*Q_0_1_Lambda_sig[sign];

    N_2_bkg_LL = Q_n1_1_Lambda_bkg[sign]*Q_n2_1_Lambda_bkg[sign] - Q_n1n2_2_Lambda_bkg[sign];
    D_2_bkg_LL = Q_0_1_Lambda_bkg[sign]*Q_0_1_Lambda_bkg[sign] - Q_0_2_Lambda_bkg[sign];
    //--------------------------------------------

    //mutiplying particle c Q-vectors

    for( int j = 0; j < 4; j++){
      if( j == 0 ){
        N_3_HFplus = N_2_sig_LL*Q_n3_1_HFplus;
        D_3_HFplus = D_2_sig_LL*Q_0_1_HFplus;

        N_3_HFminus = N_2_sig_LL*Q_n3_1_HFminus;
        D_3_HFminus = D_2_sig_LL*Q_0_1_HFminus;
      }
      if( j == 1 ){
        N_3_HFplus = N_2_sig_bkg_LL*Q_n3_1_HFplus;
        D_3_HFplus = D_2_sig_bkg_LL*Q_0_1_HFplus;

        N_3_HFminus = N_2_sig_bkg_LL*Q_n3_1_HFminus;
        D_3_HFminus = D_2_sig_bkg_LL*Q_0_1_HFminus;
      }
      if( j == 2 ){
        N_3_HFplus = N_2_bkg_sig_LL*Q_n3_1_HFplus;
        D_3_HFplus = D_2_bkg_sig_LL*Q_0_1_HFplus;

        N_3_HFminus = N_2_bkg_sig_LL*Q_n3_1_HFminus;
        D_3_HFminus = D_2_bkg_sig_LL*Q_0_1_HFminus;
      }
      if( j == 3 ){
        N_3_HFplus = N_2_bkg_LL*Q_n3_1_HFplus;
        D_3_HFplus = D_2_bkg_LL*Q_0_1_HFplus;

        N_3_HFminus = N_2_bkg_LL*Q_n3_1_HFminus;
        D_3_HFminus = D_2_bkg_LL*Q_0_1_HFminus;
      }
      
      c3_LL_real[sign][j][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
      c3_LL_real[sign][j][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
      c3_LL_imag[sign][j][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
      c3_LL_imag[sign][j][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());
    }
  }
  //lambda-antilambda correlator, no overlaps is needed:

    TComplex N_2_sig_LL, D_2_sig_LL; //sig-sig correlator
    TComplex N_2_bkg_LL, D_2_bkg_LL; //bkg-bkg correlatro
    TComplex N_2_sig_bkg_LL, D_2_sig_bkg_LL;//sig-bkg correlator (note that the n1 and n2 would make a different from the one below)
    TComplex N_2_bkg_sig_LL, D_2_bkg_sig_LL;//bkg-sig correlator

    //--------------------------------------------
    N_2_sig_LL = Q_n1_1_Lambda_sig[0]*Q_n2_1_Lambda_sig[1];
    D_2_sig_LL = Q_0_1_Lambda_sig[0]*Q_0_1_Lambda_sig[1];

    N_2_sig_bkg_LL = Q_n1_1_Lambda_sig[0]*Q_n2_1_Lambda_bkg[1];//sig-bkg doesn't overlap, no need to subtract the 3rd term
    D_2_sig_bkg_LL = Q_0_1_Lambda_sig[0]*Q_0_1_Lambda_bkg[1];

    N_2_bkg_sig_LL = Q_n1_1_Lambda_bkg[0]*Q_n2_1_Lambda_sig[1];//bkg-sig doesn't overlap, no need to subtract the 3rd term
    D_2_bkg_sig_LL = Q_0_1_Lambda_bkg[0]*Q_0_1_Lambda_sig[1];

    N_2_bkg_LL = Q_n1_1_Lambda_bkg[0]*Q_n2_1_Lambda_bkg[1];
    D_2_bkg_LL = Q_0_1_Lambda_bkg[0]*Q_0_1_Lambda_bkg[1];
    //--------------------------------------------

    //mutiplying particle c Q-vectors

    for( int j = 0; j < 4; j++){
      if( j == 0 ){
        N_3_HFplus = N_2_sig_LL*Q_n3_1_HFplus;
        D_3_HFplus = D_2_sig_LL*Q_0_1_HFplus;

        N_3_HFminus = N_2_sig_LL*Q_n3_1_HFminus;
        D_3_HFminus = D_2_sig_LL*Q_0_1_HFminus;
      }
      if( j == 1 ){
        N_3_HFplus = N_2_sig_bkg_LL*Q_n3_1_HFplus;
        D_3_HFplus = D_2_sig_bkg_LL*Q_0_1_HFplus;

        N_3_HFminus = N_2_sig_bkg_LL*Q_n3_1_HFminus;
        D_3_HFminus = D_2_sig_bkg_LL*Q_0_1_HFminus;
      }
      if( j == 2 ){
        N_3_HFplus = N_2_bkg_sig_LL*Q_n3_1_HFplus;
        D_3_HFplus = D_2_bkg_sig_LL*Q_0_1_HFplus;

        N_3_HFminus = N_2_bkg_sig_LL*Q_n3_1_HFminus;
        D_3_HFminus = D_2_bkg_sig_LL*Q_0_1_HFminus;
      }
      if( j == 3 ){
        N_3_HFplus = N_2_bkg_LL*Q_n3_1_HFplus;
        D_3_HFplus = D_2_bkg_LL*Q_0_1_HFplus;

        N_3_HFminus = N_2_bkg_LL*Q_n3_1_HFminus;
        D_3_HFminus = D_2_bkg_LL*Q_0_1_HFminus;
      }
      
      c3_LL_real[2][j][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
      c3_LL_real[2][j][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
      c3_LL_imag[2][j][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
      c3_LL_imag[2][j][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());
    }
//end of lambda-lambda 


//K0s-K0s
  //2-particles with V0s
  TComplex N_2_sig_KK, D_2_sig_KK, N_2_bkg_KK, D_2_bkg_KK;
  TComplex N_2_sig_bkg_KK, D_2_sig_bkg_KK;//sig-bkg correlator (note that the n1 and n2 would make a different from the one below)
  TComplex N_2_bkg_sig_KK, D_2_bkg_sig_KK;//bkg-sig correlator

//--------------------------------------------
  N_2_sig_KK = Q_n1_1_K0s_sig*Q_n2_1_K0s_sig - Q_n1n2_2_K0s_sig;
  D_2_sig_KK = Q_0_1_K0s_sig*Q_0_1_K0s_sig - Q_0_2_K0s_sig;

  N_2_sig_bkg_KK = Q_n1_1_K0s_sig*Q_n2_1_K0s_bkg;//sig-bkg doesn't overlap, no need to subtract the 3rd term
  D_2_sig_bkg_KK = Q_0_1_K0s_sig*Q_0_1_K0s_bkg;

  N_2_bkg_sig_KK = Q_n1_1_K0s_bkg*Q_n2_1_K0s_sig;//bkg-sig doesn't overlap, no need to subtract the 3rd term
  D_2_bkg_sig_KK = Q_0_1_K0s_bkg*Q_0_1_K0s_sig;

  N_2_bkg_KK = Q_n1_1_K0s_bkg*Q_n2_1_K0s_bkg - Q_n1n2_2_K0s_bkg;
  D_2_bkg_KK = Q_0_1_K0s_bkg*Q_0_1_K0s_bkg - Q_0_2_K0s_bkg;
//--------------------------------------------

  //mutiplying particle c Q-vectors
  for( int j = 0; j < 4; j++){
    if( j == 0 ){
      N_3_HFplus = N_2_sig_KK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_sig_KK*Q_0_1_HFplus;

      N_3_HFminus = N_2_sig_KK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_sig_KK*Q_0_1_HFminus;
    }
    if( j == 1 ){
      N_3_HFplus = N_2_sig_bkg_KK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_sig_bkg_KK*Q_0_1_HFplus;

      N_3_HFminus = N_2_sig_bkg_KK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_sig_bkg_KK*Q_0_1_HFminus;
    }
    if( j == 2 ){
      N_3_HFplus = N_2_bkg_sig_KK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_bkg_sig_KK*Q_0_1_HFplus;

      N_3_HFminus = N_2_bkg_sig_KK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_bkg_sig_KK*Q_0_1_HFminus;
    }
    if( j == 3 ){
      N_3_HFplus = N_2_bkg_KK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_bkg_KK*Q_0_1_HFplus;

      N_3_HFminus = N_2_bkg_KK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_bkg_KK*Q_0_1_HFminus;
    }
    
    c3_KK_real[j][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
    c3_KK_real[j][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
    c3_KK_imag[j][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
    c3_KK_imag[j][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());
  }
//end of K0s-K0s

//Lambda-K0s:
  TComplex N_2_sig_LK, D_2_sig_LK, N_2_bkg_LK, D_2_bkg_LK;
  TComplex N_2_sig_bkg_LK, D_2_sig_bkg_LK;//sig-bkg correlator (note that the n1 and n2 would make a different from the one below)
  TComplex N_2_bkg_sig_LK, D_2_bkg_sig_LK;//bkg-sig correlator

//--------------------------------------------
  N_2_sig_LK = Q_n1_1_Lambda_sig[0]*Q_n2_1_K0s_sig;
  D_2_sig_LK = Q_0_1_Lambda_sig[0]*Q_0_1_K0s_sig;

  N_2_sig_bkg_LK = Q_n1_1_Lambda_sig[0]*Q_n2_1_K0s_bkg;//sig-bkg doesn't overlap, no need to subtract the 3rd term
  D_2_sig_bkg_LK = Q_0_1_Lambda_sig[0]*Q_0_1_K0s_bkg;

  N_2_bkg_sig_LK = Q_n1_1_K0s_bkg*Q_n2_1_Lambda_sig[0];//bkg-sig doesn't overlap, no need to subtract the 3rd term
  D_2_bkg_sig_LK = Q_0_1_K0s_bkg*Q_0_1_Lambda_sig[0];

  N_2_bkg_LK = Q_n1_1_Lambda_bkg[0]*Q_n2_1_K0s_bkg;
  D_2_bkg_LK = Q_0_1_Lambda_bkg[0]*Q_0_1_K0s_bkg;
//--------------------------------------------

//mutiplying particle c Q-vectors
  for( int j = 0; j < 4; j++){
    if( j == 0 ){
      N_3_HFplus = N_2_sig_LK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_sig_LK*Q_0_1_HFplus;

      N_3_HFminus = N_2_sig_LK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_sig_LK*Q_0_1_HFminus;
    }
    if( j == 1 ){
      N_3_HFplus = N_2_sig_bkg_LK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_sig_bkg_LK*Q_0_1_HFplus;

      N_3_HFminus = N_2_sig_bkg_LK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_sig_bkg_LK*Q_0_1_HFminus;
    }
    if( j == 2 ){
      N_3_HFplus = N_2_bkg_sig_LK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_bkg_sig_LK*Q_0_1_HFplus;

      N_3_HFminus = N_2_bkg_sig_LK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_bkg_sig_LK*Q_0_1_HFminus;
    }
    if( j == 3 ){
      N_3_HFplus = N_2_bkg_LK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_bkg_LK*Q_0_1_HFplus;

      N_3_HFminus = N_2_bkg_LK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_bkg_LK*Q_0_1_HFminus;
    }
    
    c3_LK_real[0][j][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
    c3_LK_real[0][j][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
    c3_LK_imag[0][j][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
    c3_LK_imag[0][j][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());
  }

//--------------------------------------------
  N_2_sig_LK = Q_n1_1_Lambda_sig[1]*Q_n2_1_K0s_sig;
  D_2_sig_LK = Q_0_1_Lambda_sig[1]*Q_0_1_K0s_sig;

  N_2_sig_bkg_LK = Q_n1_1_Lambda_sig[1]*Q_n2_1_K0s_bkg;//sig-bkg doesn't overlap, no need to subtract the 3rd term
  D_2_sig_bkg_LK = Q_0_1_Lambda_sig[1]*Q_0_1_K0s_bkg;

  N_2_bkg_sig_LK = Q_n1_1_K0s_bkg*Q_n2_1_Lambda_sig[1];//bkg-sig doesn't overlap, no need to subtract the 3rd term
  D_2_bkg_sig_LK = Q_0_1_K0s_bkg*Q_0_1_Lambda_sig[1];

  N_2_bkg_LK = Q_n1_1_Lambda_bkg[1]*Q_n2_1_K0s_bkg;
  D_2_bkg_LK = Q_0_1_Lambda_bkg[1]*Q_0_1_K0s_bkg;
//--------------------------------------------

//mutiplying particle c Q-vectors
  for( int j = 0; j < 4; j++){
    if( j == 0 ){
      N_3_HFplus = N_2_sig_LK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_sig_LK*Q_0_1_HFplus;

      N_3_HFminus = N_2_sig_LK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_sig_LK*Q_0_1_HFminus;
    }
    if( j == 1 ){
      N_3_HFplus = N_2_sig_bkg_LK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_sig_bkg_LK*Q_0_1_HFplus;

      N_3_HFminus = N_2_sig_bkg_LK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_sig_bkg_LK*Q_0_1_HFminus;
    }
    if( j == 2 ){
      N_3_HFplus = N_2_bkg_sig_LK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_bkg_sig_LK*Q_0_1_HFplus;

      N_3_HFminus = N_2_bkg_sig_LK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_bkg_sig_LK*Q_0_1_HFminus;
    }
    if( j == 3 ){
      N_3_HFplus = N_2_bkg_LK*Q_n3_1_HFplus;
      D_3_HFplus = D_2_bkg_LK*Q_0_1_HFplus;

      N_3_HFminus = N_2_bkg_LK*Q_n3_1_HFminus;
      D_3_HFminus = D_2_bkg_LK*Q_0_1_HFminus;
    }
    
    c3_LK_real[1][j][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
    c3_LK_real[1][j][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
    c3_LK_imag[1][j][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
    c3_LK_imag[1][j][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());
  }
//end of Lambda-K0s


//Lambda-h
  TComplex N_2_sig_LH, D_2_sig_LH, N_2_bkg_LH, D_2_bkg_LH;

//--------------------------------------------
  N_2_sig_LH = Q_n1_1_Lambda_sig[0]*Q_forV0s_trk;
  D_2_sig_LH = Q_0_1_Lambda_sig[0]*Q_0_forV0s_trk;

  N_2_bkg_LH = Q_n1_1_Lambda_bkg[0]*Q_forV0s_trk;
  D_2_bkg_LH = Q_0_1_Lambda_bkg[0]*Q_0_forV0s_trk;
//--------------------------------------------

  //mutiplying particle c Q-vectors
  N_3_HFplus = N_2_sig_LH*Q_n3_1_HFplus;
  D_3_HFplus = D_2_sig_LH*Q_0_1_HFplus;

  N_3_HFminus = N_2_sig_LH*Q_n3_1_HFminus;
  D_3_HFminus = D_2_sig_LH*Q_0_1_HFminus;

  c3_LH_real[0][0][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
  c3_LH_real[0][0][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
  c3_LH_imag[0][0][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
  c3_LH_imag[0][0][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());

  N_3_HFplus = N_2_bkg_LH*Q_n3_1_HFplus;
  D_3_HFplus = D_2_bkg_LH*Q_0_1_HFplus;

  N_3_HFminus = N_2_bkg_LH*Q_n3_1_HFminus;
  D_3_HFminus = D_2_bkg_LH*Q_0_1_HFminus;

  c3_LH_real[0][1][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[bkg][HF]
  c3_LH_real[0][1][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
  c3_LH_imag[0][1][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[bkg][HF]
  c3_LH_imag[0][1][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());


  //anti-lambda correlates with hadron. 

  //--------------------------------------------
  N_2_sig_LH = Q_n1_1_Lambda_sig[1]*Q_forV0s_trk;
  D_2_sig_LH = Q_0_1_Lambda_sig[1]*Q_0_forV0s_trk;

  N_2_bkg_LH = Q_n1_1_Lambda_bkg[1]*Q_forV0s_trk;
  D_2_bkg_LH = Q_0_1_Lambda_bkg[1]*Q_0_forV0s_trk;
  //--------------------------------------------

  //mutiplying particle c Q-vectors
  N_3_HFplus = N_2_sig_LH*Q_n3_1_HFplus;
  D_3_HFplus = D_2_sig_LH*Q_0_1_HFplus;

  N_3_HFminus = N_2_sig_LH*Q_n3_1_HFminus;
  D_3_HFminus = D_2_sig_LH*Q_0_1_HFminus;

  c3_LH_real[1][0][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
  c3_LH_real[1][0][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
  c3_LH_imag[1][0][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
  c3_LH_imag[1][0][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());

  N_3_HFplus = N_2_bkg_LH*Q_n3_1_HFplus;
  D_3_HFplus = D_2_bkg_LH*Q_0_1_HFplus;

  N_3_HFminus = N_2_bkg_LH*Q_n3_1_HFminus;
  D_3_HFminus = D_2_bkg_LH*Q_0_1_HFminus;

  c3_LH_real[1][1][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[bkg][HF]
  c3_LH_real[1][1][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
  c3_LH_imag[1][1][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[bkg][HF]
  c3_LH_imag[1][1][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());
//end of Lambda-h

//K0s-h
  TComplex N_2_sig_KH, D_2_sig_KH, N_2_bkg_KH, D_2_bkg_KH;

//--------------------------------------------
  N_2_sig_KH = Q_n1_1_K0s_sig*Q_forV0s_trk;
  D_2_sig_KH = Q_0_1_K0s_sig*Q_0_forV0s_trk;

  N_2_bkg_KH = Q_n1_1_K0s_bkg*Q_forV0s_trk;
  D_2_bkg_KH = Q_0_1_K0s_bkg*Q_0_forV0s_trk;
//--------------------------------------------

  //mutiplying particle c Q-vectors
  N_3_HFplus = N_2_sig_KH*Q_n3_1_HFplus;
  D_3_HFplus = D_2_sig_KH*Q_0_1_HFplus;

  N_3_HFminus = N_2_sig_KH*Q_n3_1_HFminus;
  D_3_HFminus = D_2_sig_KH*Q_0_1_HFminus;

  c3_KH_real[0][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
  c3_KH_real[0][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
  c3_KH_imag[0][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[signal][HF]
  c3_KH_imag[0][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());

  N_3_HFplus = N_2_bkg_KH*Q_n3_1_HFplus;
  D_3_HFplus = D_2_bkg_KH*Q_0_1_HFplus;

  N_3_HFminus = N_2_bkg_KH*Q_n3_1_HFminus;
  D_3_HFminus = D_2_bkg_KH*Q_0_1_HFminus;

  c3_KH_real[1][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re());//[bkg][HF]
  c3_KH_real[1][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re());
  c3_KH_imag[1][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re());//[bkg][HF]
  c3_KH_imag[1][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re());

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

  edm::FileInPath fip2("CVEandMixedHarmonics/CVEandMixedHarmonics/data/Hijing_8TeV_dataBS.root");
  TFile f2(fip2.fullPath().c_str(),"READ");
  for(int i = 0; i < 1; i++){
     effTable_pPb[i] = (TH2D*)f2.Get(Form("rTotalEff3D_%d",i));
  }

  ks_mass = fs->make<TH3D>("ks_mass",";#eta;pT(GeV/c);mass(GeV/c^{2})",70,-3.5,3.5,120,0,12,360,0.44,0.56);
  la_mass_1 = fs->make<TH3D>("la_mass_1",";#eta;pT(GeV/c);mass(GeV/c^{2})",70,-3.5,3.5,120,0,12,360,1.08,1.16);
  la_mass_2 = fs->make<TH3D>("la_mass_2",";#eta;pT(GeV/c);mass(GeV/c^{2})",70,-3.5,3.5,120,0,12,360,1.08,1.16);

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  vtxZ = fs->make<TH1D>("vtxZ",";vz", 400,-20,20);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);
  trkPt = fs->make<TH1D>("trkPt", ";p_{T}(GeV)", Nptbins,ptBinsArray);
  trk_eta = fs->make<TH1D>("trk_eta", ";#eta", NetaBins, etaBinsArray);

  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 20000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 20000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 20000,-1,1);
  
  for(int deta = 0; deta < NdEtaBins; deta++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < 2; HF++){

        c3_real[deta][sign][HF] = fs->make<TH1D>(Form("c3_real_%d_%d_%d", deta, sign, HF),";c3", 20000,-1,1);
        c3_imag[deta][sign][HF] = fs->make<TH1D>(Form("c3_imag_%d_%d_%d", deta, sign, HF),";c3", 20000,-1,1);

      }
    }    
  }

  for(int sig = 0; sig < 4; sig++){
    for(int HF = 0; HF < 2; HF++){

      c3_KK_real[sig][HF] = fs->make<TH1D>(Form("c3_KK_real_%d_%d", sig, HF), ";c3", 20000,-1,1);
      c3_KK_imag[sig][HF] = fs->make<TH1D>(Form("c3_KK_imag_%d_%d", sig, HF), ";c3", 20000,-1,1);

      c3_KH_real[sig][HF] = fs->make<TH1D>(Form("c3_KH_real_%d_%d", sig, HF), ";c3", 20000,-1,1);
      c3_KH_imag[sig][HF] = fs->make<TH1D>(Form("c3_KH_imag_%d_%d", sig, HF), ";c3", 20000,-1,1);
      
      for(int sign = 0; sign < 3; sign++){

        c3_LL_real[sign][sig][HF] = fs->make<TH1D>(Form("c3_LL_real_%d_%d_%d", sign, sig, HF), ";c3", 20000,-1,1);
        c3_LL_imag[sign][sig][HF] = fs->make<TH1D>(Form("c3_LL_imag_%d_%d_%d", sign, sig, HF), ";c3", 20000,-1,1);

        c3_LK_real[sign][sig][HF] = fs->make<TH1D>(Form("c3_LK_real_%d_%d_%d", sign, sig, HF), ";c3", 20000,-1,1);
        c3_LK_imag[sign][sig][HF] = fs->make<TH1D>(Form("c3_LK_imag_%d_%d_%d", sign, sig, HF), ";c3", 20000,-1,1);

        c3_LH_real[sign][sig][HF] = fs->make<TH1D>(Form("c3_LH_real_%d_%d_%d", sign, sig, HF), ";c3", 20000,-1,1);
        c3_LH_imag[sign][sig][HF] = fs->make<TH1D>(Form("c3_LH_imag_%d_%d_%d", sign, sig, HF), ";c3", 20000,-1,1);
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
//V0s selections
bool
CVEandMixedHarmonics::passV0sCut(const reco::VertexCompositeCandidate & trk, const reco::Vertex & vtx, bool isK0s){

  double bestvz = vtx.z(); 
  double bestvx = vtx.x(); 
  double bestvy = vtx.y();
  double bestvzError = vtx.zError(); 
  double bestvxError = vtx.xError(); 
  double bestvyError = vtx.yError();

  const reco:: Candidate * d1 = trk.daughter(0);
  const reco:: Candidate * d2 = trk.daughter(1);

  auto dau1 = d1->get<reco::TrackRef>();
  auto dau2 = d2->get<reco::TrackRef>();

  double px_dau1 = d1->px();
  double py_dau1 = d1->py();
  double pz_dau1 = d1->pz();
 
  double px_dau2 = d2->px();
  double py_dau2 = d2->py();
  double pz_dau2 = d2->pz();    

  double v0s_pt = trk.pt();
  double v0s_px = trk.px();
  double v0s_py = trk.py();
  double v0s_pz = trk.pz();
  double v0s_eta = trk.eta();

  //PAngle
  double secvz = trk.vz();
  double secvx = trk.vx();
  double secvy = trk.vy();
  TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
  TVector3 secvec(v0s_px,v0s_py,v0s_pz);

  double agl = cos(secvec.Angle(ptosvec));

 //Decay length
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
  SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);

  double dl = ROOT::Math::Mag(distanceVector);
  double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
  double dlos = dl/dlerror;
  
  //NumberofValidHits for two daughters"
  double dau1_Nhits = dau1->numberOfValidHits();
  double dau2_Nhits = dau2->numberOfValidHits();

  //DCA
  math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
  
  double dzbest1 = dau1->dz(bestvtx);
  double dxybest1 = dau1->dxy(bestvtx);
  double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
  double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);

  double dzos1 = dzbest1/dzerror1;
  double dxyos1 = dxybest1/dxyerror1;
  
  double dzbest2 = dau2->dz(bestvtx);
  double dxybest2 = dau2->dxy(bestvtx);
  double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
  double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
  
  double dzos2 = dzbest2/dzerror2;
  double dxyos2 = dxybest2/dxyerror2;

  if( dau1_Nhits <= v0sNhitsCut_ || dau2_Nhits <= v0sNhitsCut_ ) return false;
  if( fabs(v0s_eta) > v0EtaTracker_ ) return false;

  if( isK0s == true ){

    if( v0s_pt < K0sPtLow_ || v0s_pt > K0sPtHigh_ ) return false;
    if( dlos < K0s_decayLengthCut_ ) return false;
    if( agl < K0s_pointingAngleCut_ ) return false;
    if( fabs(dzos1) < K0s_dcaCut_ || fabs(dzos2) < K0s_dcaCut_ || fabs(dxyos1) < K0s_dcaCut_ || fabs(dxyos2) < K0s_dcaCut_ ) return false;

    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
    double temp_reverse = Mass_ks(px_dau2,py_dau2,pz_dau2,px_dau1,py_dau1,pz_dau1);

    if ( (temp < LambdaMass+lambdaMassWindow_ && temp > LambdaMass-lambdaMassWindow_) ) return false;
    if ((temp_reverse < LambdaMass+lambdaMassWindow_ && temp_reverse > LambdaMass-lambdaMassWindow_)) return false;
    if ( temp_e < 0.015) return false;
  }
  else{

    if( v0s_pt < LamPtLow_ || v0s_pt > LamPtHigh_ ) return false;
    if( dlos < Lam_decayLengthCut_ ) return false;
    if( agl < Lam_pointingAngleCut_ ) return false;
    if( fabs(dzos1) < Lam_dcaCut_ || fabs(dzos2) < Lam_dcaCut_ || fabs(dxyos1) < Lam_dcaCut_ || fabs(dxyos2) < Lam_dcaCut_ ) return false;

    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
    
    if ( temp < K0sMass+ksMassWindow_ && temp > K0sMass-ksMassWindow_ ) return false;
    if ( temp_e < 0.015) return false;
  }

  return true;  

}
double 
CVEandMixedHarmonics::Mass_ks(double px_1,double py_1,double pz_1,double px_2,double py_2,double pz_2)
{
  double temp = 0.0;
  double E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.93827203*0.93827203));
  double E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.13957018*0.13957018));
  double E_tot = E1+E2;
  temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
  return sqrt(temp);
}

double 
CVEandMixedHarmonics::Mass_la(double px_1,double py_1,double pz_1,double px_2,double py_2,double pz_2)
{    
  double temp = 0.0;
  double E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.13957018*0.13957018));
  double E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.13957018*0.13957018));
  double E_tot = E1+E2;
  temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
  return sqrt(temp);
}

double 
CVEandMixedHarmonics::Mass_e(double px_1,double py_1,double pz_1,double px_2,double py_2,double pz_2)
{
  double temp = 0.0;
  double E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.000511*0.000511));
  double E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.000511*0.000511));
  double E_tot = E1+E2;
  temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
  return sqrt(temp);
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
