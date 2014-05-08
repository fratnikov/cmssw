
#include "SimCalorimetry/EcalZeroSuppressionProducers/interface/EcalZeroSuppressionProducer.h"

EcalZeroSuppressionProducer::EcalZeroSuppressionProducer(const edm::ParameterSet& params) 
{
  digiProducer_   = params.getParameter<std::string>("digiProducer");
  EBdigiCollection_ = params.getParameter<std::string>("EBdigiCollection");
  EEdigiCollection_ = params.getParameter<std::string>("EEdigiCollection");
  EKdigiCollection_ = params.getParameter<std::string>("EKdigiCollection");
  EBZSdigiCollection_ = params.getParameter<std::string>("EBZSdigiCollection");
  EEZSdigiCollection_ = params.getParameter<std::string>("EEZSdigiCollection");
  EKZSdigiCollection_ = params.getParameter<std::string>("EKZSdigiCollection");

  // initialize the default values for the thresholds in number of noise sigmas

  glbBarrelThreshold_ = params.getUntrackedParameter<double>("glbBarrelThreshold",0.2);
  glbEndcapThreshold_ = params.getUntrackedParameter<double>("glbEndcapThreshold",0.4);
  glbShashlikThreshold_ = params.getUntrackedParameter<double>("glbShashlikThreshold",0.4);

  produces<EBDigiCollection>(EBZSdigiCollection_);
  produces<EEDigiCollection>(EEZSdigiCollection_);
  produces<EKDigiCollection>(EKZSdigiCollection_);

}


EcalZeroSuppressionProducer::~EcalZeroSuppressionProducer() 
{ }


void EcalZeroSuppressionProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup) 
{
  
  // Get Inputs

  initCalibrations(eventSetup);
  
  edm::Handle< EBDigiCollection > pEBDigis;
  edm::Handle< EEDigiCollection > pEEDigis;
  edm::Handle< EKDigiCollection > pEKDigis;
  
  const EBDigiCollection* fullBarrelDigis =0;
  const EEDigiCollection* fullEndcapDigis =0;
  const EKDigiCollection* fullShashlikDigis =0;
  
  event.getByLabel( digiProducer_, pEBDigis);
  if (pEBDigis.isValid()){ 
    fullBarrelDigis = pEBDigis.product(); // get a ptr to the produc
    edm::LogInfo("ZeroSuppressionInfo") << "total # fullBarrelDigis: " << fullBarrelDigis->size() ;
  } else {
    edm::LogError("ZeroSuppressionError") << "Error! can't get the product " << EBdigiCollection_.c_str() ;
  }
  
  event.getByLabel( digiProducer_, pEEDigis);
  if (pEEDigis.isValid()){ 
    fullEndcapDigis = pEEDigis.product(); // get a ptr to the product
    edm::LogInfo("ZeroSuppressionInfo") << "total # fullEndcapDigis: " << fullEndcapDigis->size() ;
  } else {
    edm::LogError("ZeroSuppressionError") << "Error! can't get the product " << EEdigiCollection_.c_str() ;
  }

  event.getByLabel( digiProducer_, pEKDigis);
  if (pEKDigis.isValid()){ 
    fullShashlikDigis = pEKDigis.product(); // get a ptr to the product
    edm::LogInfo("ZeroSuppressionInfo") << "total # fullShashlikDigis: " << fullShashlikDigis->size() ;
  } else {
    edm::LogError("ZeroSuppressionError") << "Error! can't get the product " << EKdigiCollection_.c_str() ;
  }
  
  // collection of zero suppressed digis to put in the event
  
  std::auto_ptr< EBDigiCollection > gzsBarrelDigis(new EBDigiCollection());
  std::auto_ptr< EEDigiCollection > gzsEndcapDigis(new EEDigiCollection());
  std::auto_ptr< EKDigiCollection > gzsShashlikDigis(new EKDigiCollection());

  CaloDigiCollectionSorter sorter(5);

  // Barrel zero suppression

  if (fullBarrelDigis) {

    for(EBDigiCollection::const_iterator digiItr = (*fullBarrelDigis).begin();
        digiItr != (*fullBarrelDigis).end(); ++digiItr)
      {
        
        bool isAccepted = theBarrelZeroSuppressor_.accept(*digiItr, glbBarrelThreshold_);
        if (isAccepted) {
          (*gzsBarrelDigis).push_back(digiItr->id(), digiItr->begin());
        }
        
      }
    edm::LogInfo("ZeroSuppressionInfo") << "EB Digis: " << gzsBarrelDigis->size();


    //std::vector<EBDataFrame> sortedDigisEB = sorter.sortedVector(*gzsBarrelDigis);
    //LogDebug("ZeroSuppressionDump") << "Top 10 EB digis";
    //for(int i = 0; i < std::min(10,(int) sortedDigisEB.size()); ++i) 
    //  {
    //    LogDebug("ZeroSuppressionDump") << sortedDigisEB[i];
    //  }
  }
  
  // Endcap zero suppression

  if (fullEndcapDigis) {

    for(EEDigiCollection::const_iterator digiItr = (*fullEndcapDigis).begin();
        digiItr != (*fullEndcapDigis).end(); ++digiItr)
      {
        
        bool isAccepted = theEndcapZeroSuppressor_.accept(*digiItr, glbEndcapThreshold_);
        if (isAccepted) {
          (*gzsEndcapDigis).push_back(digiItr->id(), digiItr->begin());
        }
        
      }
    edm::LogInfo("ZeroSuppressionInfo") << "EE Digis: " << gzsBarrelDigis->size();
  }

  if (fullShashlikDigis) {
    for(EKDigiCollection::const_iterator digiItr = (*fullShashlikDigis).begin();
        digiItr != (*fullShashlikDigis).end(); ++digiItr)
      {
        bool isAccepted = theShashlikZeroSuppressor_.accept(*digiItr, glbShashlikThreshold_, true);
        if (isAccepted || glbShashlikThreshold_==0) {
          (*gzsShashlikDigis).push_back(digiItr->id(), digiItr->begin());
        }
      }
    edm::LogInfo("ZeroSuppressionInfo") << "EK Digis: " << gzsShashlikDigis->size();
  }

  // Step D: Put outputs into event
  event.put(gzsBarrelDigis, EBZSdigiCollection_);
  event.put(gzsEndcapDigis, EEZSdigiCollection_);
  event.put(gzsShashlikDigis, EKZSdigiCollection_);

}


void EcalZeroSuppressionProducer::initCalibrations(const edm::EventSetup & eventSetup) {

  // Pedestals from event setup
                                                                                                                                                             
  edm::ESHandle<EcalPedestals> dbPed;
  eventSetup.get<EcalPedestalsRcd>().get( dbPed );
  const EcalPedestals * thePedestals=dbPed.product();

  theBarrelZeroSuppressor_.setPedestals( thePedestals );
  theEndcapZeroSuppressor_.setPedestals( thePedestals );
  theShashlikZeroSuppressor_.setPedestals( thePedestals );
                                                                                                                                                             
}
