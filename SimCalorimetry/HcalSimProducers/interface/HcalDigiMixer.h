#ifndef SimCalorimetry_HcalSimProducers_HcalDigiMixer_h
#define SimCalorimetry_HcalSimProducers_HcalDigiMixer_h

#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixMod.h"
#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalCoderFactory.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include <vector>

namespace edm {
  class ConsumesCollector;
  namespace one {
    class EDProducerBase;
  }
  class ParameterSet;
  class StreamID;
}

template <typename DIGIS>
class HcalDigiMixer : public DigiAccumulatorMixMod {
 public:
  typedef typename DIGIS::value_type DIGI;
  typedef typename DIGIS::key_type ID;
  typedef typename DIGIS::iterator iterator;
  
 HcalDigiMixer(edm::ParameterSet const& pset, edm::one::EDProducerBase& mixMod, edm::ConsumesCollector& iC) 
   : digiTag (pset.getParameter< edm::InputTag > ("src"))
    {
      std::cout << "--------------> HcalDigiMixer 1 <---------------" << std::endl;
      mixMod.produces<DIGIS> ();
    }
  
 HcalDigiMixer(edm::ParameterSet const& pset, edm::ConsumesCollector& iC)
   : digiTag (pset.getParameter< edm::InputTag > ("src"))
    {
      std::cout << "--------------> HcalDigiMixer 2 <---------------" << std::endl;
    }

  virtual void initializeEvent(edm::Event const& ev, edm::EventSetup const& es) override {
    digiCache = std::auto_ptr<DIGIS> (new DIGIS());
  }

  virtual void finalizeEvent(edm::Event& ev, edm::EventSetup const& es) override {
    std::cout << "HcalDigiMixer::finalizeEvent-> " << digiCache->size() << std::endl;
    ev.put (digiCache);
  }

  virtual void accumulate(edm::Event const& ev, edm::EventSetup const& es) override {
    std::cout << "--------------> HcalDigiMixer Accumulate Signal <---------------" << ev.eventAuxiliary() << std::endl;
    accumulateT (ev, es);
    /* int* a = 0; *a = 0;  */
  }

  virtual void accumulate(PileUpEventPrincipal const& ev, edm::EventSetup const& es, edm::StreamID const&) override {
    std::cout << "--------------> HcalDigiMixer Accumulate PU <---------------" << ev.principal().id()  << std::endl;
    accumulateT (ev, es);
  }

 private:
  std::auto_ptr<DIGIS> digiCache;
  edm::InputTag digiTag;
  void addDigis (DIGIS* result, const DIGIS& addon);
  int getCapId (const typename DIGIS::value_type& frame, int sample) const;
  template <typename EVENT>
    void accumulateT(EVENT const&, edm::EventSetup const&);
};

typedef HcalDigiMixer<HBHEDigiCollection> HcalDigiMixerHBHE;
typedef HcalDigiMixer<HODigiCollection> HcalDigiMixerHO;
typedef HcalDigiMixer<HFDigiCollection> HcalDigiMixerHF;
typedef HcalDigiMixer<ZDCDigiCollection> HcalDigiMixerZDC;
//typedef HcalDigiMixer<HBHEUpgradeDigiCollection> HcalDigiMixerHBHEUpgrade;
//typedef HcalDigiMixer<HFUpgradeDigiCollection> HcalDigiMixerHFUpgrade;


template <typename DIGIS>
void HcalDigiMixer<DIGIS>::addDigis (DIGIS* result, const DIGIS& addon) {
  if (result->empty()) {
    *result = addon;
  }
  else {
    auto coder (HcalCoderFactory (HcalCoderFactory::NOMINAL).coder(DetId(0)));
    for (auto& newDigi : addon) {
      ID id = newDigi.id();
      iterator sumDigi = result->find (id);
      if (sumDigi != result->end()) {
	CaloSamples sumCS, newCS;
	coder->adc2fC (*sumDigi, sumCS);
	coder->adc2fC (newDigi, newCS);
	int capIdOffset = getCapId (*sumDigi, 0);
	sumCS += newCS;
	DIGI newFrame;
	coder->fC2adc (sumCS, newFrame, capIdOffset);
	for (int frame = 0; frame < sumDigi->size(); ++frame) {
	  sumDigi->setSample (frame, newFrame[frame]);
	}
      }
      else {
	result->push_back (newDigi);
      }
    }
  }
}

 template <typename DIGIS>
 int HcalDigiMixer<DIGIS>::getCapId (const typename DIGIS::value_type& frame, int i) const {
   return frame.sample(i).capid();
 }

 template <>
   int HcalDigiMixer<HBHEUpgradeDigiCollection>::getCapId (const HcalDigiMixer<HBHEUpgradeDigiCollection>::DIGI& frame, int i) const {
   return frame.capId(i);
 }

 /* template <> */ 
 /*   int HcalDigiMixer<HFUpgradeDigiCollection>::getCapId (const HFUpgradeDigiCollection::value_type& frame, int i) const { */
 /*   return frame.capId(i); */
 /* } */


template <typename DIGIS>
template <typename EVENT>
  void HcalDigiMixer<DIGIS>::accumulateT(EVENT const& ev, edm::EventSetup const& es) {
  edm::Handle<DIGIS> newDigis;
  ev.getByLabel (digiTag, newDigis);
  addDigis (&*digiCache, *newDigis);
}

#endif
