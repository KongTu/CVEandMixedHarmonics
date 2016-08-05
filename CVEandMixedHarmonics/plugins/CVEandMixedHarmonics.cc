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
//         Created:  Thu, 04 Aug 2016 13:09:20 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class CVEandMixedHarmonics : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit CVEandMixedHarmonics(const edm::ParameterSet&);
      ~CVEandMixedHarmonics();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
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
CVEandMixedHarmonics::CVEandMixedHarmonics(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");

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



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
CVEandMixedHarmonics::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CVEandMixedHarmonics::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CVEandMixedHarmonics::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CVEandMixedHarmonics);
