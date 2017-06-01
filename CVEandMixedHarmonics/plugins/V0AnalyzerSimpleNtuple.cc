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
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "FWCore/Utilities/interface/GCC11Compatibility.h"
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
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
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

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
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
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
//
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

//
// class decleration
//

#define PI 3.1416
using namespace std;
using namespace reco;
using namespace edm;


float e_mass = 0.000511;

float pion_mass = 0.13957018;
float proton_mass = 0.938272013;
float kaon_mass = 0.493677;
float pion_sigma = pion_mass*1.e-6;
float proton_sigma = proton_mass*1.e-6;
float kaon_sigma = kaon_mass*1.e-6;
float lambda_mass = 1.115683;
float lambda_sigma = 0.000006;

float ks_mass = 0.497614;
float xi_mass = 1.32171;

const int PDGID[4] = {310,3122,3312,3334};
const double PDGMASS[4] = {0.497614,1.115683,1.32171,1.67245};
const string PNAMES[4] = {"ks","lambda","xi","omega"};



class V0AnalyzerSimpleNtuple : public edm::EDAnalyzer {
public:
  explicit V0AnalyzerSimpleNtuple(const edm::ParameterSet&);
  ~V0AnalyzerSimpleNtuple();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool MatchV0(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::VertexCompositeCandidate & trk, int pdgid, string genPLabel, double & genpt_matched);
   

  // ----------member data ---------------------------
    
    TNtuple* V0AnalyzerSimpleNtuple_ks;
    TNtuple* V0AnalyzerSimpleNtuple_la;
    TNtuple* V0AnalyzerSimpleNtuple_genks;
    TNtuple* V0AnalyzerSimpleNtuple_genla;

    bool doGenParticle_;
    bool doSimParticle_;

    edm::EDGetTokenT<edm::View<reco::Track> > trackSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleSrc_;
    edm::EDGetTokenT<TrackingParticleCollection> tpHitsSrc_;

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

V0AnalyzerSimpleNtuple::V0AnalyzerSimpleNtuple(const edm::ParameterSet& iConfig):
    tpHitsSrc_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpHitsSrc")))
{

  //now do what ever initialization is needed

    trackSrc_ = consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("trackSrc"));
    vertexSrc_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"));
    genParticleSrc_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"));

    generalV0_ks_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("generalV0_ks"));
    generalV0_la_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("generalV0_la"));

    doGenParticle_ = iConfig.getUntrackedParameter<bool>("doGenParticle", false);
    doSimParticle_ = iConfig.getUntrackedParameter<bool>("doSimParticle", false);

    
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
    using namespace reco;
    
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

    double nMult_ass_good = 0.0;
    for(unsigned it = 0; it < tracks->size(); it++){

        const reco::Track & trk = (*tracks)[it];

        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > 0.1 ) continue;
        if(fabs(dzvtx/dzerror) > 3.0 ) continue;
        if(fabs(dxyvtx/dxyerror) > 3.0 ) continue;
        if( fabs(trk.eta()) < 2.4 && trk.pt() > 0.4 ){nMult_ass_good++;}// NtrkOffline        

        const reco::HitPattern & p = trk.hitPattern();

        HitCategory hitCat = 0;
        //loop over the hits of the track
        for (int i=0; i<p.numberOfHits(hitCat); i++) {
        uint32_t hit = p.getHitPattern(hitCat,i);
        cout << "hit: " << hit << endl;
        }



    }

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
    iEvent.getByToken(generalV0_ks_,v0candidates_ks);
    if(!v0candidates_ks.isValid()) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
    iEvent.getByToken(generalV0_la_,v0candidates_la);
    if(!v0candidates_la.isValid()) return;

    edm::Handle<TrackingParticleCollection> tpCollection;
    iEvent.getByToken(tpHitsSrc_, tpCollection);


    if( doSimParticle_ ){

        for(TrackingParticleCollection::size_type i=0; i<tpCollection->size(); i++){
            
            TrackingParticleRef tpr(tpCollection, i);
            TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

            //TrackingParticle selections:
            if(tp->status() < 0 || tp->charge()==0 || tp->pt()<1 || TMath::Abs(tp->eta())>2.4) continue;

            cout << "test" << endl;
        }

    }

    if( doGenParticle_ ){

      edm::Handle<reco::GenParticleCollection> genParticleCollection;
      iEvent.getByToken(genParticleSrc_, genParticleCollection);

      for(unsigned it=0; it<genParticleCollection->size(); ++it) {

        const reco::GenParticle & genCand = (*genParticleCollection)[it];
        int id = genCand.pdgId();
        int status = genCand.status();
        double genpt = genCand.pt();
        double geneta = genCand.eta();

        if(abs(genCand.pdgId()) != 310) continue;

        if(genCand.numberOfDaughters() != 2) continue;
        if(genCand.daughter(0)->vertex() != genCand.daughter(1)->vertex()) continue;

        int posGenDauNdx = 0;
        int negGenDauNdx = 1;
        if(genCand.daughter(0)->charge() < 0){
         
          posGenDauNdx = 1;
          negGenDauNdx = 0;
        }

        const Candidate* genPosDau = genCand.daughter(posGenDauNdx);
        const Candidate* genNegDau = genCand.daughter(negGenDauNdx);


    //   if ( geneta < -2.4 || geneta > 2.4 ) continue;
    
    //     if ( status == 1 ){

    //       if( id == 310 ){

    //         V0AnalyzerSimpleNtuple_genks->Fill(genpt,geneta,genCand.mass());

    //       }

    // //Finding mother:
    //     int mid = 0;
    //       if( TMath::Abs(id) == 3122 ){

    //         if(genCand.numberOfMothers()==1){
    //           const reco::Candidate * mom = genCand.mother();
    //           mid = mom->pdgId();
    //           if(mom->numberOfMothers()==1){
    //             const reco::Candidate * mom1 = mom->mother();
    //             mid = mom1->pdgId();
    //           }
    //         }

    //         if (TMath::Abs(mid) != 3322 && TMath::Abs(mid) != 3312 && TMath::Abs(mid) != 3324 && TMath::Abs(mid) != 3314 && TMath::Abs(mid) != 3334){

    //          V0AnalyzerSimpleNtuple_genla->Fill(genpt, geneta, genCand.mass());

    //         }
    //       }
     
    //     }
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

        if( fabs(dzos1) < 1.0 || fabs(dzos2) < 1.0 || fabs(dxyos1) < 1.0 || fabs(dxyos2) < 1.0 ) continue;
        if( dlos < 5.0 ) continue;
        if( agl < 0.999 ) continue;

        double genPt = -0.99;
        double matchedK0s =  MatchV0(iEvent, iSetup, trk, 310, "genParticlesPlusSim", genPt);

        // hit pattern of the track
        // const reco::HitPattern& p = dau1->hitPattern();

        // // loop over the hits of the track
        // for (int i=0; i<p.numberOfHits(); i++) {
        // uint32_t hit = p.getHitPattern(i);

        // // if the hit is valid and in pixel barrel, print out the layer
        // if (p.validHitFilter(hit) && p.pixelBarrelHitFilter(hit))
        // std::cout << "valid hit found in pixel barrel layer "
        //       << p.getLayer(hit) << std::endl;
        // }

        // // count the number of valid tracker *** hits ***
        // std::cout << "number of of valid tracker hits is "
        //   << p.numberOfValidTrackerHits() << std::endl;

        // // count the number of tracker *** layers *** with measurement
        // std::cout << "number of tracker layers with measurement is "
        //   << p.trackerLayersWithMeasurement() << std::endl;

        //algo
        int dau1_algo = dau1->algo();
        int dau2_algo = dau2->algo();

        int dau1_originAlgo = dau1->originalAlgo();
        int dau2_originAlgo = dau2->originalAlgo();

        //inner most hit Det id
        unsigned int id_1;
        id_1 = dau1->innerDetId();
        unsigned int id_2;
        id_2 = dau2->innerDetId();

        DetId detId(id_1);
        unsigned int subid1 = detId.subdetId();
        DetId detId2(id_2);
        unsigned int subid2 = detId2.subdetId();

        int dau1_layer = 0;
        int dau2_layer = 0;

        if( subid1 == 1){
            PXBDetId pxbid(id_1);
            dau1_layer = (int)pxbid.layer();
        }
        if( subid1 == 2){
            PXFDetId pxfid(id_1);
            dau1_layer = (int)pxfid.disk();
        }
        if( subid1 == 3){
            TIBDetId tibid(id_1);
            dau1_layer = tibid.layer();
        }
        if( subid1 == 4){
            TIDDetId tidid(id_1);
            dau1_layer = (int)tidid.wheel();
        }
        if( subid1 == 5){
            TOBDetId tobid(id_1);
            dau1_layer = tobid.layer();
        }
        if( subid1 == 6){
            TECDetId tecid(id_1);
            dau1_layer = (int)tecid.wheel();
        }


        if( subid2 == 1){
            PXBDetId pxbid(id_2);
            dau2_layer = (int)pxbid.layer();
        }
        if( subid2 == 2){
            PXFDetId pxfid(id_2);
            dau2_layer = (int)pxfid.disk();
        }
        if( subid2 == 3){
            TIBDetId tibid(id_2);
            dau2_layer = tibid.layer();
        }
        if( subid2 == 4){
            TIDDetId tidid(id_2);
            dau2_layer = (int)tidid.wheel();
        }
        if( subid2 == 5){
            TOBDetId tobid(id_2);
            dau2_layer = tobid.layer();
        }
        if( subid2 == 6){
            TECDetId tecid(id_2);
            dau2_layer = (int)tecid.wheel();
        }

        //Fill
        V0AnalyzerSimpleNtuple_ks->Fill(pt,eta,mass,nMult_ass_good,dau1_Nhits,dau2_Nhits,dau1_algo,dau2_algo,dau1_originAlgo,dau2_originAlgo,subid1,dau1_layer,subid2,dau2_layer, matchedK0s);
            
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
        
        if( fabs(dzos1) < 1.0 || fabs(dzos2) < 1.0 || fabs(dxyos1) < 1.0 || fabs(dxyos2) < 1.0 ) continue;
        if( dlos < 5.0 ) continue;
        if( agl < 0.999 ) continue;

        double genPt = -0.99;
        double matchedLam =  MatchV0(iEvent, iSetup, trk, 3122, "genParticlesPlusSim", genPt);

        //algo
        int dau1_algo = dau1->algo();
        int dau2_algo = dau2->algo();
       
        int dau1_originAlgo = dau1->originalAlgo();
        int dau2_originAlgo = dau2->originalAlgo();

        //inner most hit Det id
        unsigned int id_1;
        id_1 = dau1->innerDetId();
        unsigned int id_2;
        id_2 = dau2->innerDetId();

        DetId detId(id_1);
        unsigned int subid1 = detId.subdetId();
        DetId detId2(id_2);
        unsigned int subid2 = detId2.subdetId();

        int dau1_layer = 0;
        int dau2_layer = 0;

        if( subid1 == 1){
            PXBDetId pxbid(id_1);
            dau1_layer = (int)pxbid.layer();
        }
        if( subid1 == 2){
            PXFDetId pxfid(id_1);
            dau1_layer = (int)pxfid.disk();
        }
        if( subid1 == 3){
            TIBDetId tibid(id_1);
            dau1_layer = tibid.layer();
        }
        if( subid1 == 4){
            TIDDetId tidid(id_1);
            dau1_layer = (int)tidid.wheel();
        }
        if( subid1 == 5){
            TOBDetId tobid(id_1);
            dau1_layer = tobid.layer();
        }
        if( subid1 == 6){
            TECDetId tecid(id_1);
            dau1_layer = (int)tecid.wheel();
        }


        if( subid2 == 1){
            PXBDetId pxbid(id_2);
            dau2_layer = (int)pxbid.layer();
        }
        if( subid2 == 2){
            PXFDetId pxfid(id_2);
            dau2_layer = (int)pxfid.disk();
        }
        if( subid2 == 3){
            TIBDetId tibid(id_2);
            dau2_layer = tibid.layer();
        }
        if( subid2 == 4){
            TIDDetId tidid(id_2);
            dau2_layer = (int)tidid.wheel();
        }
        if( subid2 == 5){
            TOBDetId tobid(id_2);
            dau2_layer = tobid.layer();
        }
        if( subid2 == 6){
            TECDetId tecid(id_2);
            dau2_layer = (int)tecid.wheel();
        }



        //Fill
        V0AnalyzerSimpleNtuple_la->Fill(pt,eta,mass,nMult_ass_good,dau1_Nhits,dau2_Nhits,dau1_algo,dau2_algo,dau1_originAlgo,dau2_originAlgo,subid1,dau1_layer,subid2,dau2_layer,matchedLam);
        
    }

}


// ------------ method called once each job just before starting event
//loop  ------------
void 
V0AnalyzerSimpleNtuple::beginJob()
{
    edm::Service<TFileService> fs;
        
    TH1D::SetDefaultSumw2();
    
    V0AnalyzerSimpleNtuple_ks = fs->make< TNtuple>("V0AnalyzerSimpleNtuple_ks","V0AnalyzerSimpleNtuple_ks","pt:eta:mass:ntrk:trkNHits1:trkNHits2:trkAlgo1:trkAlgo2:trkOriginalAlgo1:trkOriginalAlgo2:trkHitDet1:trkHitDet2:trkHitLayer1:trkHitLayer2:matched");
    V0AnalyzerSimpleNtuple_la = fs->make< TNtuple>("V0AnalyzerSimpleNtuple_la","V0AnalyzerSimpleNtuple_la","pt:eta:mass:ntrk:trkNHits1:trkNHits2:trkAlgo1:trkAlgo2:trkOriginalAlgo1:trkOriginalAlgo2:trkHitDet1:trkHitDet2:trkHitLayer1:trkHitLayer2:matched");

    if( doGenParticle_ ){

        V0AnalyzerSimpleNtuple_genks = fs->make< TNtuple>("V0AnalyzerSimpleNtuple_genks","V0AnalyzerSimpleNtuple_genks","ks_pt:ks_eta:ks_mass");
        V0AnalyzerSimpleNtuple_genla = fs->make< TNtuple>("V0AnalyzerSimpleNtuple_genla","V0AnalyzerSimpleNtuple_genla","la_pt:la_eta:la_mass");
    }
}
bool
V0AnalyzerSimpleNtuple::MatchV0(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::VertexCompositeCandidate & trk, int pdgid, string genPLabel, double & genpt_matched)
{
        using namespace edm;

        edm::Handle<reco::GenParticleCollection> genCol;
        iEvent.getByLabel(genPLabel, genCol);

        bool isMatched = false;
        int posRecoDauNdx = 0;
        int negRecoDauNdx = 1;

        if(trk.daughter(0)->charge() < 0)
        {      
               posRecoDauNdx = 1;
               negRecoDauNdx = 0;
        }

        const Candidate* recoPosDau = trk.daughter(posRecoDauNdx);
        const Candidate* recoNegDau = trk.daughter(negRecoDauNdx);

        double recoPosDauPhi = recoPosDau->momentum().phi();
        double recoNegDauPhi = recoNegDau->momentum().phi();

        double recoPosDauEta = recoPosDau->momentum().eta();
        double recoNegDauEta = recoNegDau->momentum().eta();

        GlobalPoint recoVtx(trk.vx(), trk.vy(), trk.vz());

        double deltaL = -1;
        double posDeltaR = -1;
        double negDeltaR = -1;
        for(unsigned igen = 0; igen < genCol->size(); igen++)
        {
               const reco::GenParticle & genCand = (*genCol)[igen];
               if(abs(genCand.pdgId()) != pdgid) continue;

               if(genCand.numberOfDaughters() != 2) continue;
               if(genCand.daughter(0)->vertex() != genCand.daughter(1)->vertex()) continue;

               int posGenDauNdx = 0;
               int negGenDauNdx = 1;
               if(genCand.daughter(0)->charge() < 0)
               {
                      posGenDauNdx = 1;
                      negGenDauNdx = 0;
               }

               const Candidate* genPosDau = genCand.daughter(posGenDauNdx);
               const Candidate* genNegDau = genCand.daughter(negGenDauNdx);

               double genPosDauPhi = genPosDau->momentum().phi();
               double genNegDauPhi = genNegDau->momentum().phi();

               double genPosDauEta = genPosDau->momentum().eta();
               double genNegDauEta = genNegDau->momentum().eta();

               double posDeltaPhi = reco::deltaPhi(genPosDauPhi, recoPosDauPhi);
               double negDeltaPhi = reco::deltaPhi(genNegDauPhi, recoNegDauPhi);
               double posDeltaEta = genPosDauEta - recoPosDauEta;
               double negDeltaEta = genNegDauEta - recoNegDauEta;

               posDeltaR = sqrt(posDeltaPhi*posDeltaPhi + posDeltaEta*posDeltaEta);
               negDeltaR = sqrt(negDeltaPhi*negDeltaPhi + negDeltaEta*negDeltaEta);
               GlobalPoint genVtx(genPosDau->vx(), genPosDau->vy(), genPosDau->vz());
               deltaL = (genVtx - recoVtx).mag();
                
               // if(pdgid == 310) 
               // {
               //        map_posDeltaR["ks"]->Fill(posDeltaR, genCand.pt());
               //        map_negDeltaR["ks"]->Fill(negDeltaR, genCand.pt());
               //        map_batDeltaR["ks"]->Fill(0., genCand.pt());
               //        map_deltaL["ks"]->Fill(deltaL, genCand.pt()); 
               // }
              
               // if(pdgid == 3122)
               // {
               //        map_posDeltaR["lambda"]->Fill(posDeltaR, genCand.pt());
               //        map_negDeltaR["lambda"]->Fill(negDeltaR, genCand.pt());
               //        map_batDeltaR["lambda"]->Fill(0., genCand.pt());
               //        map_deltaL["lambda"]->Fill(deltaL, genCand.pt());                            }

               if(posDeltaR < 0.1 && negDeltaR < 0.1 && deltaL < 10.)
               {
                      isMatched = true;
                      genpt_matched = genCand.pt();
               }
        }
        return isMatched;
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
V0AnalyzerSimpleNtuple::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0AnalyzerSimpleNtuple);


