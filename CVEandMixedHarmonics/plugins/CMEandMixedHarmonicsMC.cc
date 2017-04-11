// -*- C++ -*-
//
// Package:    CMEandMixedHarmonicsMC/CMEandMixedHarmonicsMC
// Class:      CMEandMixedHarmonicsMC
// 
/**\class CMEandMixedHarmonicsMC CMEandMixedHarmonicsMC.cc CMEandMixedHarmonicsMC/CMEandMixedHarmonicsMC/plugins/CMEandMixedHarmonicsMC.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Mon, 01 Aug 2016 09:01:02 GMT
//
//


#include "CVEandMixedHarmonics/CVEandMixedHarmonics/interface/CMEandMixedHarmonicsMCBase.h"


CMEandMixedHarmonicsMC::CMEandMixedHarmonicsMC(const edm::ParameterSet& iConfig)
{

  trackName_  =  iConfig.getParameter<edm::InputTag>("trackName");
  vertexName_ =  iConfig.getParameter<edm::InputTag>("vertexName");
  towerName_ =  iConfig.getParameter<edm::InputTag>("towerName");
  genName_ =  iConfig.getParameter<edm::InputTag>("genName");
  hepName_ = iConfig.getParameter<edm::InputTag>("hepName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);
  genSrc_ = consumes<reco::GenParticleCollection>(genName_);
  hepSrc_ = consumes<HepMCProduct>(hepName_);

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
  doGenOnly_ = iConfig.getUntrackedParameter<bool>("doGenOnly");

  eff_ = iConfig.getUntrackedParameter<int>("eff");

  etaTracker_ = iConfig.getUntrackedParameter<double>("etaTracker");
  gapValue_ = iConfig.getUntrackedParameter<double>("gapValue");
  
  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");
  
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


CMEandMixedHarmonicsMC::~CMEandMixedHarmonicsMC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CMEandMixedHarmonicsMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace HepMC;

  edm::Handle<HepMCProduct> mc;
  iEvent.getByToken(hepSrc_,mc);

  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  iEvent.getByToken(genSrc_, genParticleCollection);
 
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexSrc_,vertices);

  Handle<CaloTowerCollection> towers;
  iEvent.getByToken(towerSrc_, towers);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackSrc_, tracks);

  if( doGenOnly_ ){

    int nTracks_gen = 0;
    const GenEvent* evt = mc->GetEvent();

    HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
    HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
    for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){

      if((*it)->status() != 1) continue;

       int pdg_id = (*it)->pdg_id();
       const ParticleData * part = pdt->particle(pdg_id);
       int charge = static_cast<int>(part->charge());
       if(charge == 0) continue;

       if((*it)->momentum().perp()<0.01) continue;
       if(fabs((*it)->momentum().eta())>20) continue;

       double eta = (*it)->momentum().eta();
       double pt = (*it)->momentum().perp();

       if( pt < 0.4 || fabs(eta) > 2.4 ) continue;
       nTracks_gen++;
    }

    if( nTracks_gen < Nmin_ || nTracks_gen >= Nmax_ ) return;
  }
  else{

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

    if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;
    Ntrk->Fill( nTracks );

  }

  const int NetaBins = etaBins_.size() - 1 ;
  const int NdEtaBins = dEtaBins_.size() - 1;

  const int NptBins = ptBins_.size() - 1;
  const int NdPtBins = dPtBins_.size() - 1;

  double dEtaBinsArray[100];
  double dPtBinsArray[100];

  for(unsigned i = 0; i < dEtaBins_.size(); i++){

    dEtaBinsArray[i] = dEtaBins_[i]-0.0001;
  }
  for(unsigned i = 0; i < dPtBins_.size(); i++){

    dPtBinsArray[i] = dPtBins_[i]-0.0001;
  }

/*
define all the ingredients for CME correlator (gamma) <cos(n1*phi_1 + n2*phi_2 + n3*phi_3)>
and delta correlator <cos(phi_1 - phi_2)> ,
where Q_coefficient_power (P_coefficient_power) is used in the following names 
 */

/*
1. eta dimension:
*/

//3-particle correlator q-vectors (charge dependent for 0 and 1)

  TComplex Q_n1_1[NetaBins][2], Q_n2_1[NetaBins][2];
  
  TComplex Q_n1n2_2[NetaBins][2];
  
  TComplex Q_0_1[NetaBins][2], Q_0_2[NetaBins][2];

//2-particle correlator q-vectors

  TComplex P_n1_1[NetaBins][2], P_n2_1[NetaBins][2];

  TComplex P_n1n2_2[NetaBins][2];

  TComplex P_0_1[NetaBins][2], P_0_2[NetaBins][2];

/*
2. pT dimension:
*/

//3-particle correlator q-vectors (charge dependent for 0 and 1)

  TComplex Q_pT_n1_1[NptBins][2], Q_pT_n2_1[NptBins][2];
  
  TComplex Q_pT_n1n2_2[NptBins][2];
  
  TComplex Q_pT_0_1[NptBins][2], Q_pT_0_2[NptBins][2];

//2-particle correlator q-vectors

  TComplex P_pT_n1_1[NptBins][2], P_pT_n2_1[NptBins][2];

  TComplex P_pT_n1n2_2[NptBins][2];

  TComplex P_pT_0_1[NptBins][2], P_pT_0_2[NptBins][2];

/*
Share Q_n3 for both dimensions:
*/

//Scalar product, charge independent, |eta|<2.4

  TComplex Q_n3_trk, Q_0_trk;

//Tracker v2

  TComplex Q_nC_trk[NetaBins], Q_0_nC_trk[NetaBins];


//------------------------------------------------------------------

  const GenEvent* evt = mc->GetEvent();

  HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
  HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
  for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){

    if((*it)->status() != 1) continue;

    int pdg_id = (*it)->pdg_id();
    const ParticleData * part = pdt->particle(pdg_id);
    int gencharge = static_cast<int>(part->charge());
    if(gencharge == 0) continue;

    if((*it)->momentum().perp()<0.01) continue;
    if(fabs((*it)->momentum().eta())>20) continue;

    double geneta = (*it)->momentum().eta();
    double genphi = (*it)->momentum().phi();
    double genpt = (*it)->momentum().perp();
    double weight = 1.0;

    if( genpt < ptLow_ || genpt > ptHigh_) continue;
    if( fabs(geneta) > etaTracker_ ) continue;

    trkPhi->Fill(genphi, weight);
    trkPt->Fill(genpt, weight);
    trk_eta->Fill(geneta, weight);

    Q_n3_trk += q_vector(-n3_, 1, weight, genphi);//for scalar product in tracker
    Q_0_trk += q_vector(0, 1, weight, genphi);

    for(int eta = 0; eta < NetaBins; eta++){
      if( geneta > etaBins_[eta] && geneta < etaBins_[eta+1] ){

        Q_nC_trk[eta] += q_vector(-n3_, 1, weight, genphi);
        Q_0_nC_trk[eta] += q_vector(0, 1, weight, genphi);

        if( gencharge == +1 ){//positive charge

          //3p:
          Q_n1_1[eta][0] += q_vector(n1_, 1, weight, genphi);
          Q_n2_1[eta][0] += q_vector(n2_, 1, weight, genphi);

          Q_n1n2_2[eta][0] += q_vector(n1_+n2_, 2, weight, genphi);

          Q_0_1[eta][0] += q_vector(0, 1, weight, genphi);
          Q_0_2[eta][0] += q_vector(0, 2, weight, genphi);

          //2p: (similar way but be careful of the order of harmonics)

          P_n1_1[eta][0] += q_vector(n1_, 1, weight, genphi);
          P_n2_1[eta][0] += q_vector(-n2_, 1, weight, genphi);//it is a minus n2_ because n2_ = 1

          P_n1n2_2[eta][0] += q_vector(n1_-n2_, 2, weight, genphi);

          P_0_1[eta][0] += q_vector(0, 1, weight, genphi);
          P_0_2[eta][0] += q_vector(0, 2, weight, genphi);


        }
        if( gencharge == -1 ){//negative charge

          Q_n1_1[eta][1] += q_vector(n1_, 1, weight, genphi);
          Q_n2_1[eta][1] += q_vector(n2_, 1, weight, genphi);

          Q_n1n2_2[eta][1] += q_vector(n1_+n2_, 2, weight, genphi);

          Q_0_1[eta][1] += q_vector(0, 1, weight, genphi);
          Q_0_2[eta][1] += q_vector(0, 2, weight, genphi);

          //2p: (similar way but be careful of the order of harmonics)

          P_n1_1[eta][1] += q_vector(n1_, 1, weight, genphi);
          P_n2_1[eta][1] += q_vector(-n2_, 1, weight, genphi);//it is a minus n2_ because n2_ = 1

          P_n1n2_2[eta][1] += q_vector(n1_-n2_, 2, weight, genphi);

          P_0_1[eta][1] += q_vector(0, 1, weight, genphi);
          P_0_2[eta][1] += q_vector(0, 2, weight, genphi);

        }
      }
    }//end of eta dimension

    //begin of pT dimension
    for(int pt = 0; pt < NptBins; pt++){
      if( genpt > ptBins_[pt] && genpt < ptBins_[pt+1] ){

        if( gencharge == +1 ){//positive charge

          //3p:
          Q_pT_n1_1[pt][0] += q_vector(n1_, 1, weight, genphi);
          Q_pT_n2_1[pt][0] += q_vector(n2_, 1, weight, genphi);

          Q_pT_n1n2_2[pt][0] += q_vector(n1_+n2_, 2, weight, genphi);

          Q_pT_0_1[pt][0] += q_vector(0, 1, weight, genphi);
          Q_pT_0_2[pt][0] += q_vector(0, 2, weight, genphi);

          //2p: (similar way but be careful of the order of harmonics)

          P_pT_n1_1[pt][0] += q_vector(n1_, 1, weight, genphi);
          P_pT_n2_1[pt][0] += q_vector(-n2_, 1, weight, genphi);//it is a minus n2_ because n2_ = 1

          P_pT_n1n2_2[pt][0] += q_vector(n1_-n2_, 2, weight, genphi);

          P_pT_0_1[pt][0] += q_vector(0, 1, weight, genphi);
          P_pT_0_2[pt][0] += q_vector(0, 2, weight, genphi);


        }
        if( gencharge == -1 ){//negative charge

          Q_pT_n1_1[pt][1] += q_vector(n1_, 1, weight, genphi);
          Q_pT_n2_1[pt][1] += q_vector(n2_, 1, weight, genphi);

          Q_pT_n1n2_2[pt][1] += q_vector(n1_+n2_, 2, weight, genphi);

          Q_pT_0_1[pt][1] += q_vector(0, 1, weight, genphi);
          Q_pT_0_2[pt][1] += q_vector(0, 2, weight, genphi);

          //2p: (similar way but be careful of the order of harmonics)

          P_pT_n1_1[pt][1] += q_vector(n1_, 1, weight, genphi);
          P_pT_n2_1[pt][1] += q_vector(-n2_, 1, weight, genphi);//it is a minus n2_ because n2_ = 1

          P_pT_n1n2_2[pt][1] += q_vector(n1_-n2_, 2, weight, genphi);

          P_pT_0_1[pt][1] += q_vector(0, 1, weight, genphi);
          P_pT_0_2[pt][1] += q_vector(0, 2, weight, genphi);

        }
      }
    }//end of pT dimension
  }//end of gen loop
  
  //Declaring TComplex varaibles for Q-vectors of particle c (HF)

  TComplex  Q_n3_1_HFplus, Q_n3_1_HFminus, Q_0_1_HFplus, Q_0_1_HFminus;
  
  for(unsigned it=0; it<genParticleCollection->size(); ++it) {

    const reco::GenParticle & genCand = (*genParticleCollection)[it];
    int status = genCand.status();
    double geneta = genCand.eta();
    double genphi = genCand.phi();
    int gencharge = genCand.charge();
    double w = 1.0;

    if( status != 1 || gencharge == 0 ) continue;

    if( reverseBeam_ ) geneta = -genCand.eta();
    if( geneta < etaHighHF_ && geneta > etaLowHF_ ){

        Q_n3_1_HFplus += q_vector(n3_, 1, w, genphi);
        Q_0_1_HFplus += q_vector(0, 1, w, genphi);

    }
    else if( geneta < -etaLowHF_ && geneta > -etaHighHF_ ){

        Q_n3_1_HFminus += q_vector(n3_, 1, w, genphi);
        Q_0_1_HFminus += q_vector(0, 1, w, genphi);
    }
    else{continue;}

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
1. eta dimension:
*/

/*
calculate the 3-particles correlator with the charged-particles
 */

  for(int ieta = 0; ieta < NetaBins; ieta++){
    for(int jeta = 0; jeta < NetaBins; jeta++){

      double deltaEta = fabs( etaBins_[ieta] - etaBins_[jeta] );

      if( deltaEta > 2.0 ){//calculate the tracker v2 with a gap of 2

          TComplex N_2;
          TComplex D_2;

          N_2 = Q_nC_trk[ieta]*TComplex::Conjugate(Q_nC_trk[jeta]);
          D_2 = Q_0_nC_trk[ieta]*Q_0_nC_trk[jeta];

          cn_tracker->Fill(N_2.Re()/D_2.Re(), D_2.Re());
      }

      for(int deta = 0; deta < NdEtaBins; deta++){
        if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta+1] ){
          
          TComplex N_2;
          TComplex D_2;

          TComplex N_3_HFplus, N_3_HFminus;
          TComplex D_3_HFplus, D_3_HFminus;

          //same sign correlator:
          for(int sign = 0; sign < 2; sign++){
            if( ieta == jeta ){

              delEta3p[sign]->Fill( deltaEta );

              N_2 = Q_n1_1[ieta][sign]*Q_n2_1[ieta][sign] - Q_n1n2_2[ieta][sign];
              D_2 = Q_0_1[ieta][sign]*Q_0_1[ieta][sign] - Q_0_2[ieta][sign];

            }
            else{

              delEta3p[sign]->Fill( deltaEta );

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

          delEta3p[2]->Fill( deltaEta );

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
calculate the 2-particles correlator with the charged-particles
*/

  for(int ieta = 0; ieta < NetaBins; ieta++){
    for(int jeta = 0; jeta < NetaBins; jeta++){

      double deltaEta = fabs( etaBins_[ieta] - etaBins_[jeta] );

      for(int deta = 0; deta < NdEtaBins; deta++){
        if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta+1] ){
          
          TComplex N_2;
          TComplex D_2;

          //same sign correlator:
          for(int sign = 0; sign < 2; sign++){
            if( ieta == jeta ){

              delEta2p[sign]->Fill( deltaEta );


              N_2 = P_n1_1[ieta][sign]*P_n2_1[ieta][sign] - P_n1n2_2[ieta][sign];
              D_2 = P_0_1[ieta][sign]*P_0_1[ieta][sign] - P_0_2[ieta][sign];

            }
            else{

              delEta2p[sign]->Fill( deltaEta );

              N_2 = P_n1_1[ieta][sign]*P_n2_1[jeta][sign];
              D_2 = P_0_1[ieta][sign]*P_0_1[jeta][sign];

            }        

            c2_real[deta][sign]->Fill(N_2.Re()/D_2.Re(), D_2.Re() );
            c2_imag[deta][sign]->Fill(N_2.Im()/D_2.Re(), D_2.Re() );

          }

          delEta2p[2]->Fill( deltaEta );

          //opposite sign correlator:
          N_2 = P_n1_1[ieta][0]*P_n2_1[jeta][1];
          D_2 = P_0_1[ieta][0]*P_0_1[jeta][1];

          c2_real[deta][2]->Fill(N_2.Re()/D_2.Re(), D_2.Re() );
          c2_imag[deta][2]->Fill(N_2.Im()/D_2.Re(), D_2.Re() );

        }
      }
    }
  }

/*
2. pt dimension:
*/

/*
calculate the 3-particles correlator with the charged-particles
 */

  for(int ipt = 0; ipt < NptBins; ipt++){
    for(int jpt = 0; jpt < NptBins; jpt++){

      double deltaPt = fabs( ptBins_[ipt] - ptBins_[jpt] );
      double pTave = (ptBins_[ipt] + ptBins_[jpt])/2.0;

      for(int dpt = 0; dpt < NdPtBins; dpt++){
        //begin of delta pT
        if( deltaPt > dPtBinsArray[dpt] && deltaPt < dPtBinsArray[dpt+1] ){
          
          TComplex N_2;
          TComplex D_2;

          TComplex N_3_HFplus, N_3_HFminus;
          TComplex D_3_HFplus, D_3_HFminus;

          //same sign correlator:
          for(int sign = 0; sign < 2; sign++){
            if( ipt == jpt ){

              delPt3p[sign]->Fill( deltaPt );

              N_2 = Q_pT_n1_1[ipt][sign]*Q_pT_n2_1[ipt][sign] - Q_pT_n1n2_2[ipt][sign];
              D_2 = Q_pT_0_1[ipt][sign]*Q_pT_0_1[ipt][sign] - Q_pT_0_2[ipt][sign];

            }
            else{

              delPt3p[sign]->Fill( deltaPt );

              N_2 = Q_pT_n1_1[ipt][sign]*Q_pT_n2_1[jpt][sign];
              D_2 = Q_pT_0_1[ipt][sign]*Q_pT_0_1[jpt][sign];

            }        

            N_3_HFplus = N_2*Q_n3_1_HFplus;
            D_3_HFplus = D_2*Q_0_1_HFplus;

            c3_dpT_real[dpt][sign][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re() );
            c3_dpT_imag[dpt][sign][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re() );

            N_3_HFminus = N_2*Q_n3_1_HFminus;
            D_3_HFminus = D_2*Q_0_1_HFminus;

            c3_dpT_real[dpt][sign][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re() );
            c3_dpT_imag[dpt][sign][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re() );
  
          }
          //opposite sign correlator:

          delPt3p[2]->Fill( deltaPt );

          N_2 = Q_pT_n1_1[ipt][0]*Q_pT_n2_1[jpt][1];
          D_2 = Q_pT_0_1[ipt][0]*Q_pT_0_1[jpt][1];

          N_3_HFplus = N_2*Q_n3_1_HFplus;
          D_3_HFplus = D_2*Q_0_1_HFplus;

          c3_dpT_real[dpt][2][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re() );
          c3_dpT_imag[dpt][2][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re() );

          N_3_HFminus = N_2*Q_n3_1_HFminus;
          D_3_HFminus = D_2*Q_0_1_HFminus;

          c3_dpT_real[dpt][2][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re() );
          c3_dpT_imag[dpt][2][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re() );

        }//end of delta pT

        //begin of pTave
        if( pTave > dPtBinsArray[dpt] && pTave < dPtBinsArray[dpt+1] ){
          
          TComplex N_2;
          TComplex D_2;

          TComplex N_3_HFplus, N_3_HFminus;
          TComplex D_3_HFplus, D_3_HFminus;

          //same sign correlator:
          for(int sign = 0; sign < 2; sign++){
            if( ipt == jpt ){

              ptAve3p[sign]->Fill( pTave );

              N_2 = Q_pT_n1_1[ipt][sign]*Q_pT_n2_1[ipt][sign] - Q_pT_n1n2_2[ipt][sign];
              D_2 = Q_pT_0_1[ipt][sign]*Q_pT_0_1[ipt][sign] - Q_pT_0_2[ipt][sign];

            }
            else{

              ptAve3p[sign]->Fill( pTave );

              N_2 = Q_pT_n1_1[ipt][sign]*Q_pT_n2_1[jpt][sign];
              D_2 = Q_pT_0_1[ipt][sign]*Q_pT_0_1[jpt][sign];

            }        

            N_3_HFplus = N_2*Q_n3_1_HFplus;
            D_3_HFplus = D_2*Q_0_1_HFplus;

            c3_pTave_real[dpt][sign][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re() );
            c3_pTave_imag[dpt][sign][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re() );

            N_3_HFminus = N_2*Q_n3_1_HFminus;
            D_3_HFminus = D_2*Q_0_1_HFminus;

            c3_pTave_real[dpt][sign][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re() );
            c3_pTave_imag[dpt][sign][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re() );
  
          }
          //opposite sign correlator:

          ptAve3p[2]->Fill( pTave );

          N_2 = Q_pT_n1_1[ipt][0]*Q_pT_n2_1[jpt][1];
          D_2 = Q_pT_0_1[ipt][0]*Q_pT_0_1[jpt][1];

          N_3_HFplus = N_2*Q_n3_1_HFplus;
          D_3_HFplus = D_2*Q_0_1_HFplus;

          c3_pTave_real[dpt][2][0]->Fill(N_3_HFplus.Re()/D_3_HFplus.Re(), D_3_HFplus.Re() );
          c3_pTave_imag[dpt][2][0]->Fill(N_3_HFplus.Im()/D_3_HFplus.Re(), D_3_HFplus.Re() );

          N_3_HFminus = N_2*Q_n3_1_HFminus;
          D_3_HFminus = D_2*Q_0_1_HFminus;

          c3_pTave_real[dpt][2][1]->Fill(N_3_HFminus.Re()/D_3_HFminus.Re(), D_3_HFminus.Re() );
          c3_pTave_imag[dpt][2][1]->Fill(N_3_HFminus.Im()/D_3_HFminus.Re(), D_3_HFminus.Re() );

        }//end of pTave
      }
    }
  }

/*
calculate the 2-particles correlator with the charged-particles
*/

  for(int ipt = 0; ipt < NptBins; ipt++){
    for(int jpt = 0; jpt < NptBins; jpt++){

      double deltaPt = fabs( ptBins_[ipt] - ptBins_[jpt] );
      double pTave = (ptBins_[ipt] + ptBins_[jpt])/2.0;
      
      for(int dpt = 0; dpt < NdPtBins; dpt++){
        //begin of delta pT
        if( deltaPt > dPtBinsArray[dpt] && deltaPt < dPtBinsArray[dpt+1] ){
          
          TComplex N_2;
          TComplex D_2;

          //same sign correlator:
          for(int sign = 0; sign < 2; sign++){
            if( ipt == jpt ){

              delPt2p[sign]->Fill( deltaPt );


              N_2 = P_pT_n1_1[ipt][sign]*P_pT_n2_1[ipt][sign] - P_pT_n1n2_2[ipt][sign];
              D_2 = P_pT_0_1[ipt][sign]*P_pT_0_1[ipt][sign] - P_pT_0_2[ipt][sign];

            }
            else{

              delPt2p[sign]->Fill( deltaPt );

              N_2 = P_pT_n1_1[ipt][sign]*P_pT_n2_1[jpt][sign];
              D_2 = P_pT_0_1[ipt][sign]*P_pT_0_1[jpt][sign];

            }        

            c2_dpT_real[dpt][sign]->Fill(N_2.Re()/D_2.Re(), D_2.Re() );
            c2_dpT_imag[dpt][sign]->Fill(N_2.Im()/D_2.Re(), D_2.Re() );

          }

          delPt2p[2]->Fill( deltaPt );

          //opposite sign correlator:
          N_2 = P_pT_n1_1[ipt][0]*P_pT_n2_1[jpt][1];
          D_2 = P_pT_0_1[ipt][0]*P_pT_0_1[jpt][1];

          c2_dpT_real[dpt][2]->Fill(N_2.Re()/D_2.Re(), D_2.Re() );
          c2_dpT_imag[dpt][2]->Fill(N_2.Im()/D_2.Re(), D_2.Re() );

        }//end of delta pT
        
        //begin of pTave
        if( pTave > dPtBinsArray[dpt] && pTave < dPtBinsArray[dpt+1] ){
          
          TComplex N_2;
          TComplex D_2;

          //same sign correlator:
          for(int sign = 0; sign < 2; sign++){
            if( ipt == jpt ){

              ptAve2p[sign]->Fill( pTave );

              N_2 = P_pT_n1_1[ipt][sign]*P_pT_n2_1[ipt][sign] - P_pT_n1n2_2[ipt][sign];
              D_2 = P_pT_0_1[ipt][sign]*P_pT_0_1[ipt][sign] - P_pT_0_2[ipt][sign];

            }
            else{

              ptAve2p[sign]->Fill( pTave );

              N_2 = P_pT_n1_1[ipt][sign]*P_pT_n2_1[jpt][sign];
              D_2 = P_pT_0_1[ipt][sign]*P_pT_0_1[jpt][sign];

            }        

            c2_pTave_real[dpt][sign]->Fill(N_2.Re()/D_2.Re(), D_2.Re() );
            c2_pTave_imag[dpt][sign]->Fill(N_2.Im()/D_2.Re(), D_2.Re() );

          }

          ptAve2p[2]->Fill( pTave );

          //opposite sign correlator:
          N_2 = P_pT_n1_1[ipt][0]*P_pT_n2_1[jpt][1];
          D_2 = P_pT_0_1[ipt][0]*P_pT_0_1[jpt][1];

          c2_pTave_real[dpt][2]->Fill(N_2.Re()/D_2.Re(), D_2.Re() );
          c2_pTave_imag[dpt][2]->Fill(N_2.Im()/D_2.Re(), D_2.Re() );

        }//end of pTave
      }
    }
  }


}
// ------------ method called once each job just before starting event loop  ------------
void 
CMEandMixedHarmonicsMC::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  const int NdEtaBins = dEtaBins_.size() - 1;
  const int NetaBins = etaBins_.size() - 1;

  const int NdPtBins = dPtBins_.size() - 1;
  
  double etaBinsArray[100];
  for(unsigned i = 0; i < etaBins_.size(); i++){
    etaBinsArray[i] = etaBins_[i];
  }

  double dEtaBinsArray[100];
  for(unsigned i = 0; i < dEtaBins_.size(); i++){
    dEtaBinsArray[i] = dEtaBins_[i]-0.0001;
  }

  double ptBinsArray[100];
  const int Nptbins = ptBins_.size() - 1;
  for(unsigned i = 0; i < ptBins_.size(); i++){
    ptBinsArray[i] = ptBins_[i];
  }

  double dPtBinsArray[100];
  for(unsigned i = 0; i < dPtBins_.size(); i++){

    dPtBinsArray[i] = dPtBins_[i]-0.0001;
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
  trk_eta = fs->make<TH1D>("trk_eta", ";#eta", NetaBins, etaBinsArray);
  q2_mag = fs->make<TH1D>("q2_mag", "q2", 2000,-1,1);
  Ntrk_q2 = fs->make<TH1D>("Ntrk_q2",";Ntrk",5000,0,5000);

  for(int sign = 0; sign < 3; sign++){

    delEta2p[sign] = fs->make<TH1D>(Form("delEta2p_%d",sign),";#Delta#eta", NdEtaBins, dEtaBinsArray);
    delEta3p[sign] = fs->make<TH1D>(Form("delEta3p_%d",sign),";#Delta#eta", NdEtaBins, dEtaBinsArray);
    
    delPt2p[sign] = fs->make<TH1D>(Form("delPt2p_%d",sign),";#Deltap_{T}", NdPtBins, dPtBinsArray);
    delPt3p[sign] = fs->make<TH1D>(Form("delPt3p_%d",sign),";#Deltap_{T}", NdPtBins, dPtBinsArray);
    
    ptAve2p[sign] = fs->make<TH1D>(Form("ptAve2p_%d",sign),";bar{p_{T}}", NdPtBins, dPtBinsArray);
    ptAve3p[sign] = fs->make<TH1D>(Form("ptAve3p_%d",sign),";bar{p_{T}}", NdPtBins, dPtBinsArray);
    
  }

  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 20000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 20000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 20000,-1,1);

  cn_tracker = fs->make<TH1D>("cn_tracker",";cn_tracker", 20000,-1,1);

  
  for(int deta = 0; deta < NdEtaBins; deta++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < 2; HF++){

        c3_real[deta][sign][HF] = fs->make<TH1D>(Form("c3_real_%d_%d_%d", deta, sign, HF),";c3", 20000,-1,1);
        c3_imag[deta][sign][HF] = fs->make<TH1D>(Form("c3_imag_%d_%d_%d", deta, sign, HF),";c3", 20000,-1,1);
      }
    }    
  }

  for(int dpt = 0; dpt < NdPtBins; dpt++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < 2; HF++){

        c3_dpT_real[dpt][sign][HF] = fs->make<TH1D>(Form("c3_dpT_real_%d_%d_%d", dpt, sign, HF),";c3", 20000,-1,1);
        c3_dpT_imag[dpt][sign][HF] = fs->make<TH1D>(Form("c3_dpT_imag_%d_%d_%d", dpt, sign, HF),";c3", 20000,-1,1);

        c3_pTave_real[dpt][sign][HF] = fs->make<TH1D>(Form("c3_pTave_real_%d_%d_%d", dpt, sign, HF),";c3", 20000,-1,1);
        c3_pTave_imag[dpt][sign][HF] = fs->make<TH1D>(Form("c3_pTave_imag_%d_%d_%d", dpt, sign, HF),";c3", 20000,-1,1);

      }
    }    
  }

  for(int deta = 0; deta < NdEtaBins; deta++){
    for(int sign = 0; sign < 3; sign++){

      c2_real[deta][sign] = fs->make<TH1D>(Form("c2_real_%d_%d", deta, sign),";c2", 20000,-1,1);
      c2_imag[deta][sign] = fs->make<TH1D>(Form("c2_imag_%d_%d", deta, sign),";c2", 20000,-1,1);
    }    
  }

  for(int dpt = 0; dpt < NdPtBins; dpt++){
    for(int sign = 0; sign < 3; sign++){

      c2_dpT_real[dpt][sign] = fs->make<TH1D>(Form("c2_dpT_real_%d_%d", dpt, sign),";c2", 20000,-1,1);
      c2_dpT_imag[dpt][sign] = fs->make<TH1D>(Form("c2_dpT_imag_%d_%d", dpt, sign),";c2", 20000,-1,1);
   
      c2_pTave_real[dpt][sign] = fs->make<TH1D>(Form("c2_pTave_real_%d_%d", dpt, sign),";c2", 20000,-1,1);
      c2_pTave_imag[dpt][sign] = fs->make<TH1D>(Form("c2_pTave_imag_%d_%d", dpt, sign),";c2", 20000,-1,1);
   
    }    
  }

}
TComplex 
CMEandMixedHarmonicsMC::q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}
// ------------ method called once each job just after ending the event loop  ------------
void 
CMEandMixedHarmonicsMC::endJob() 
{
}
void 
CMEandMixedHarmonicsMC::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CMEandMixedHarmonicsMC::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CMEandMixedHarmonicsMC::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CMEandMixedHarmonicsMC::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CMEandMixedHarmonicsMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CMEandMixedHarmonicsMC);
