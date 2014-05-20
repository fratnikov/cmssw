// -*- C++ -*-
//
// Package:    ShashlikAnalyzer
// Class:      ShashlikAnalyzer
// 
/**\class ShashlikAnalyzer ShashlikAnalyzer.cc Validation/ShashlikAnalyzer/plugins/ShashlikAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shervin Nourbakhsh
//         Created:  Mon, 12 May 2014 07:26:43 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EcalDigi/interface/EKDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "CalibCalorimetry/EcalTrivialCondModules/interface/EcalTrivialConditionRetriever.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"



#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/ShashlikNumberingRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalCommonData/interface/ShashlikDDDConstants.h"
#include "Geometry/CaloTopology/interface/ShashlikTopology.h"
#include "Geometry/CaloTopology/interface/ShashlikGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

//
// class declaration
//

class ShashlikAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ShashlikAnalyzer(const edm::ParameterSet&);
      ~ShashlikAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  TFile *tree_file;
  TTree *tree;

  Int_t     	runNumber;   ///< 
  Long64_t      eventNumber; ///<
  Int_t         lumiBlock;   ///< lumi section
  UInt_t 	runTime;     ///< unix time

  Int_t genIndex;
  Float_t genEnergy[10];
  Float_t genEta[10];
  Float_t genPhi[10];

  Int_t         sampleADC[10], sampleGain[10];
  Float_t         adcToGeV;

  Int_t         ix, iy,iz;
  UInt_t        rawId;

  Float_t hitEnergy, hitTime, hitEta, hitPhi;

  edm::InputTag EKdigiCollection_;
  edm::InputTag EKrecHitCollection_;

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
ShashlikAnalyzer::ShashlikAnalyzer(const edm::ParameterSet& iConfig):
  EKdigiCollection_(iConfig.getParameter<edm::InputTag>("EKdigiCollection")),
  EKrecHitCollection_(iConfig.getParameter<edm::InputTag>("EKrecHitCollection"))
{
}


ShashlikAnalyzer::~ShashlikAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ShashlikAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle< std::vector<reco::GenParticle> > genParticlesHandle;
   iEvent.getByLabel("genParticles", genParticlesHandle);

   Handle<EKDigiCollection> EcalDigiEK;
   iEvent.getByLabel(EKdigiCollection_ , EcalDigiEK );
   
   Handle<EKRecHitCollection> EcalRecHitsEK;
  iEvent.getByLabel(EKrecHitCollection_ , EcalRecHitsEK );
  
  // ADC -> GeV Scale
  edm::ESHandle<EcalADCToGeVConstant> pAgc;
  iSetup.get<EcalADCToGeVConstantRcd>().get(pAgc);
  const EcalADCToGeVConstant* adctogev = pAgc.product();

//   edm::ESHandle<CaloSubdetectorGeometry> shgeo;
//   iSetup.get<ShashlikGeometryRecord>().get( shgeo );
//   if(! shgeo.isValid())
//     std::cout << "Cannot get a valid ShashlikGeometry Object\n";

  
  edm::ESHandle<CaloGeometry>               hGeometry   ;
  iSetup.get<CaloGeometryRecord>().get( hGeometry ) ;

  //  assert(hGeometry->getSubdetectorGeometry( DetId::Ecal, EcalShashlik  ) == shgeo.product());

  
  //const CaloGeometry* pGeometry = &*hGeometry;
  //  const CaloSubdetectorGeometry* geometry = shgeo.product();
  const CaloSubdetectorGeometry* geometry = hGeometry->getSubdetectorGeometry( DetId::Ecal, EcalShashlik  );

  //  const CaloSubdetectorGeometry *subGeo = 
  //  assert(subGeo==geometry);
  //const std::vector<DetId>& ids = geometry->getValidDetIds();

//   edm::ESHandle<ShashlikTopology> topo;
//   iSetup.get<ShashlikNumberingRecord>().get(topo);
//   const ShashlikTopology* topology = topo.product();

  // Return if no Shashlik data available
  if( !EcalDigiEK.isValid() ) return;

  edm::Timestamp runTime_  = iEvent.eventAuxiliary().time();
  runTime = runTime_.unixTime();
  runNumber = iEvent.id().run();
  eventNumber = iEvent.id().event();
  if( iEvent.isRealData() ) {
    lumiBlock = iEvent.luminosityBlock();
  } else {
    lumiBlock = -1;
  }

  assert(EcalDigiEK->size()!=0);
  assert(EcalDigiEK->size()==EcalRecHitsEK->size());

  genIndex=0;
  assert(genParticlesHandle->size()>0);
  for(std::vector<reco::GenParticle>::const_iterator gen_itr = genParticlesHandle->begin();
      gen_itr != genParticlesHandle->end();
      gen_itr++){
    if((abs(gen_itr->pdgId())!=11 && abs(gen_itr->pdgId())!=13)) continue; // || (gen_itr->status()!=3))    continue;
    assert(genIndex<10);
    genEnergy[genIndex] = gen_itr->energy();
    genEta[genIndex]    = gen_itr->eta();
    genPhi[genIndex]    = gen_itr->phi();
    genIndex++;
  }

  for (EKDigiCollection::const_iterator digi_itr = EcalDigiEK->begin();
       digi_itr!=EcalDigiEK->end();
       digi_itr++){
    
    EKDataFrame eedf=*digi_itr;
    assert(eedf.size()==10);
    
    EKDetId eeid = eedf.id();
    ix=eeid.ix();
    iy=eeid.iy();
    iz=eeid.zside();
    rawId = eeid.rawId();

    adcToGeV = adctogev->getEKValue();
        
    for (int sample = 0 ; sample < 10; ++sample) {
      
      EcalMGPASample mySample = eedf[sample];
      
      sampleADC[sample] = mySample.adc();
      //sampleGain[sample]  = mySample.gainId() ;
      if(mySample.gainId()==1) sampleGain[sample]=1;
      else if( mySample.gainId()==2) sampleGain[sample]=2;
      else sampleGain[sample]=12;
    }
    EKRecHitCollection::const_iterator recHit_itr = EcalRecHitsEK->find(eeid);
    hitEnergy = recHit_itr->energy();
    hitTime   = recHit_itr->time();

  
    
    const CaloCellGeometry * cell = geometry->getGeometry(eeid);
    hitEta = cell->etaPos();
    hitPhi = cell->phiPos();
  
    if(iz>0)
      tree->Fill(); // fill for every hit/digi

  }


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
ShashlikAnalyzer::beginJob()
{
  tree_file = new TFile("shashlik.root", "recreate");
  if(tree_file->IsZombie()){
    throw cms::Exception("OutputError") <<  "Output tree not created (Zombie): ";
    return;
  }
  tree_file->cd();
  
  //now do what ever initialization is needed
  tree = new TTree("shashlik","");
  //tree = fs->make<TTree>("selected","selected"); //no otherwise you have the extraCalib in the same file
  
  tree->SetDirectory(tree_file);
  // controllo su tree==NULL
   //now do what ever initialization is needed

  //tree->Branch("runNumber",     &runNumber,     "runNumber/I");
  tree->Branch("eventNumber",   &eventNumber, "eventNumber/l");
  //tree->Branch("lumiBlock",     &lumiBlock,     "lumiBlock/I");
  //tree->Branch("runTime",       &runTime,         "runTime/i");

  tree->Branch("genIndex", &genIndex, "genIndex/I");
  tree->Branch("genEnergy", genEnergy, "genEnergy[genIndex]/F");
  tree->Branch("genEta", genEta, "genEta[genIndex]/F");
  tree->Branch("genPhi", genPhi, "genPhi[genIndex]/F");

  tree->Branch("sampleADC", sampleADC, "sampleADC[10]/I");
  tree->Branch("sampleGain", sampleGain, "sampleGain[10]/I");
  tree->Branch("adcToGeV", &adcToGeV, "adcToGeV/F");

  tree->Branch("ix", &ix, "ix/I");
  tree->Branch("iy", &iy, "iy/I");
  tree->Branch("iz", &iz, "iz/I");
  tree->Branch("rawId", &rawId, "rawId/i");

  tree->Branch("hitEnergy", &hitEnergy, "hitEnergy/F");
  tree->Branch("hitTime",   &hitTime,   "hitTime/F");
  tree->Branch("hitEta", &hitEta, "hitEta/F");
  tree->Branch("hitPhi", &hitPhi, "hitPhi/F");
  


}

// ------------ method called once each job just after ending the event loop  ------------
void 
ShashlikAnalyzer::endJob() 
{
  // save the tree into the file
  tree_file->cd();
  tree->Write();
  tree_file->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
ShashlikAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ShashlikAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ShashlikAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ShashlikAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ShashlikAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ShashlikAnalyzer);
