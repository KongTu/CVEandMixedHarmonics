#ifndef CMEandMixedHarmonicsMCBase_
#define CMEandMixedHarmonicsMCBase_


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
#include <TGraph.h>
#include <TComplex.h>

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

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"


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

#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#define PI 2
#define K0sMass 0.497614
#define LambdaMass 1.115683

using namespace std;
using namespace reco;
using namespace edm;


class CMEandMixedHarmonicsMC : public edm::EDAnalyzer {
   public:
      explicit CMEandMixedHarmonicsMC(const edm::ParameterSet&);
      ~CMEandMixedHarmonicsMC();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual TComplex q_vector(double n, double p, double w, double phi);

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
      edm::EDGetTokenT<reco::TrackCollection> trackSrc_;
      edm::EDGetTokenT<CaloTowerCollection> towerSrc_;
      edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;
      edm::EDGetTokenT<HepMCProduct> hepSrc_;

      edm::InputTag vertexName_;
      edm::InputTag trackName_;
      edm::InputTag pfCandName_;
      edm::InputTag towerName_;
      edm::InputTag genName_;
      edm::InputTag hepName_;

      edm::ESHandle<ParticleDataTable> pdt;

      //correction table
      TH2D* effTable[5];
      TH2D* effTable_pPb[5];

      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

      TH1D* cn_tracker;

      TH1D* cn_eta[48][2];

      TH1D* c2_real[48][3];
      TH1D* c2_imag[48][3];

      TH1D* c3_real[48][3][2];
      TH1D* c3_imag[48][3][2];

      TH1D* c2_dpT_real[48][3];
      TH1D* c2_dpT_imag[48][3];

      TH1D* c3_dpT_real[48][3][2];
      TH1D* c3_dpT_imag[48][3][2];

      TH1D* c2_pTave_real[48][3];
      TH1D* c2_pTave_imag[48][3];

      TH1D* c3_pTave_real[48][3][2];
      TH1D* c3_pTave_imag[48][3][2];

      TH1D* delEta3p[3];
      TH1D* delEta2p[3];

      TH1D* delPt3p[3];
      TH1D* delPt2p[3];

      TH1D* ptAve3p[3];
      TH1D* ptAve2p[3];

      TH1D* Ntrk;
      TH1D* Ntrk_MB;
      TH1D* Ntrk_HF;
      TH1D* vtxZ;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* trkPt;
      TH1D* trk_eta;
      TH1D* cbinHist;
      TH1D* q2_mag;
      TH1D* Ntrk_q2;

      TH2D* mother_Spectra;

      int Nmin_;
      int Nmax_;

      int eff_;

      int n1_;
      int n2_;
      int n3_;
      int n4_;

      double etaTracker_;
      double gapValue_;
      double etaLowHF_;
      double etaHighHF_;
      
      double vzLow_;
      double vzHigh_;
      
      double ptLow_;
      double ptHigh_;

      double q2max_;
      double q2min_;
      
      double offlineptErr_;
      double offlineDCA_;
      double offlineChi2_;
      double offlinenhits_;

      bool useCentrality_;
      bool reverseBeam_;
      bool doEffCorrection_;
      bool useEtaGap_;
      bool dopPb_;
      bool doGenOnly_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;
      std::vector<double> ptBins_;
      std::vector<double> dPtBins_;
      std::vector<double> centBins_;

};

#endif