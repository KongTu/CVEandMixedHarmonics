// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>


#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//////////////////////////////////////////////
// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerLayerIdAccessor.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"

#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Heavyion
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations     
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//


//
// class decleration
//

#define PI 3.1416
using namespace std;

class V0AnalyzerSimpleNtuple : public edm::EDAnalyzer {
public:
  explicit V0AnalyzerSimpleNtuple(const edm::ParameterSet&);
  ~V0AnalyzerSimpleNtuple();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
    
    TNtuple* V0AnalyzerSimpleNtuple_ks;
    TNtuple* V0AnalyzerSimpleNtuple_la;
    TNtuple* V0AnalyzerSimpleNtuple_genks;
    TNtuple* V0AnalyzerSimpleNtuple_genla;
    
    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double multMax_;
    double multMin_;
    int bkgFactor_;
    bool doGenParticle_;

    edm::EDGetTokenT<edm::View<reco::Track> > trackSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleSrc_;

    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> generalV0_ks_;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> generalV0_la_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

V0AnalyzerSimpleNtuple::V0AnalyzerSimpleNtuple(const edm::ParameterSet& iConfig)
{

  //now do what ever initialization is needed

    trackSrc_ = consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("trackSrc"));
    vertexSrc_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"));
    genParticleSrc_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"));

    generalV0_ks_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("generalV0_ks"));
    generalV0_la_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("generalV0_la"));

    etaMin_trg_ = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_ = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_ = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_ = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 220);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 185);
    bkgFactor_ = iConfig.getUntrackedParameter<int>("bkgFactor", 10);

    doGenParticle_ = iConfig.getUntrackedParameter<bool>("doGenParticle", false);

    
}


V0AnalyzerSimpleNtuple::~V0AnalyzerSimpleNtuple()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
V0AnalyzerSimpleNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& 
iSetup)
{
    using std::vector;
    using namespace edm;
    
    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertexSrc_,vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    if(bestvz < -15.0 || bestvz>15.0) return;

    edm::Handle<edm::View<reco::Track>> tracks;
    iEvent.getByToken(trackSrc_, tracks);

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
    iEvent.getByToken(generalV0_ks_,v0candidates_ks);
    if(!v0candidates_ks.isValid()) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
    iEvent.getByToken(generalV0_la_,v0candidates_la);
    if(!v0candidates_la.isValid()) return;

    if( doGenParticle_ ){

      edm::Handle<reco::GenParticleCollection> genParticleCollection;
      iEvent.getByToken(genParticleSrc_, genParticleCollection);

      for(unsigned it=0; it<genParticleCollection->size(); ++it) {

        const reco::GenParticle & genCand = (*genParticleCollection)[it];
        int id = genCand.pdgId();
        int status = genCand.status();
        double genpt = genCand.pt();
        double geneta = genCand.eta();

      if ( geneta < -2.4 || geneta > 2.4 ) continue;
    
        if ( status == 1 ){

          if( id == 310 ){

            V0AnalyzerSimpleNtuple_genks->Fill(genpt,geneta,genCand.mass());

          }

    //Finding mother:
        int mid = 0;
          if( TMath::Abs(id) == 3122 ){

            if(genCand.numberOfMothers()==1){
              const reco::Candidate * mom = genCand.mother();
              mid = mom->pdgId();
              if(mom->numberOfMothers()==1){
                const reco::Candidate * mom1 = mom->mother();
                mid = mom1->pdgId();
              }
            }

            if (TMath::Abs(mid) != 3322 && TMath::Abs(mid) != 3312 && TMath::Abs(mid) != 3324 && TMath::Abs(mid) != 3314 && TMath::Abs(mid) != 3334){

             V0AnalyzerSimpleNtuple_genla->Fill(genpt, geneta, genCand.mass());

            }
          }
     
        }
      }

  } 
    
    for(unsigned it=0; it<v0candidates_ks->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
                    
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        
        auto dau1 = d1->get<reco::TrackRef>();
        auto dau2 = d2->get<reco::TrackRef>();
        
        //pt,mass
        double eta = trk.eta();
        double pt = trk.pt();
        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        double mass = trk.mass();
        
        if(eta > 2.4 || eta < -2.4 ) continue;
        
        //PAngle
        secvz = trk.vz();
        secvx = trk.vx();
        secvy = trk.vy();
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);

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
        
        //Fill
        V0AnalyzerSimpleNtuple_ks->Fill(pt,eta,mass,dzos1,dzos2,dxyos1,dxyos2,dau1_Nhits,dau2_Nhits,dlos,agl);
            
    }
    
    
    for(unsigned it=0; it<v0candidates_la->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        
        auto dau1 = d1->get<reco::TrackRef>();
        auto dau2 = d2->get<reco::TrackRef>();
        
        //pt,mass
        double eta = trk.eta();
        double pt = trk.pt();
        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        double mass = trk.mass();
        
        if(eta > 2.4 || eta < -2.4 ) continue;
       
        //PAngle
        secvz = trk.vz();
        secvx = trk.vx();
        secvy = trk.vy();
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);

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
        
        //Fill
        V0AnalyzerSimpleNtuple_la->Fill(pt,eta,mass,dzos1,dzos2,dxyos1,dxyos2,dau1_Nhits,dau2_Nhits,dlos,agl);
        
    }

}


// ------------ method called once each job just before starting event
//loop  ------------
void 
V0AnalyzerSimpleNtuple::beginJob()
{
    edm::Service<TFileService> fs;
        
    TH1D::SetDefaultSumw2();
    
    V0AnalyzerSimpleNtuple_ks = fs->make< TNtuple>("V0AnalyzerSimpleNtuple_ks","V0AnalyzerSimpleNtuple_ks","pt:eta:mass:trkDCA1z:trkDCA2z:trkDCA1xy:trkDCA2xy:trkNHits1:trkNHits2:L:PAngle");
    V0AnalyzerSimpleNtuple_la = fs->make< TNtuple>("V0AnalyzerSimpleNtuple_la","V0AnalyzerSimpleNtuple_la","pt:eta:mass:trkDCA1z:trkDCA2z:trkDCA1xy:trkDCA2xy:trkNHits1:trkNHits2:L:PAngle");

    if( doGenParticle_ ){

        V0AnalyzerSimpleNtuple_genks = fs->make< TNtuple>("V0AnalyzerSimpleNtuple_genks","V0AnalyzerSimpleNtuple_genks","ks_pt:ks_eta:ks_mass");
        V0AnalyzerSimpleNtuple_genla = fs->make< TNtuple>("V0AnalyzerSimpleNtuple_genla","V0AnalyzerSimpleNtuple_genla","la_pt:la_eta:la_mass");
    }
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
V0AnalyzerSimpleNtuple::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0AnalyzerSimpleNtuple);


