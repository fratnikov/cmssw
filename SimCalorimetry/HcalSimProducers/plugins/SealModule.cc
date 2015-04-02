#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixModFactory.h"
#include "SimCalorimetry/HcalSimProducers/src/HcalHitAnalyzer.h"
#include "SimCalorimetry/HcalSimProducers/src/HcalDigiAnalyzer.h"
#include "SimCalorimetry/HcalSimProducers/interface/HcalDigiProducer.h"
#include "SimCalorimetry/HcalSimProducers/interface/HcalDigiMixer.h"

DEFINE_FWK_MODULE(HcalHitAnalyzer);
DEFINE_FWK_MODULE(HcalDigiAnalyzer);
DEFINE_DIGI_ACCUMULATOR(HcalDigiProducer);
DEFINE_DIGI_ACCUMULATOR(HcalDigiMixerHBHE);
DEFINE_DIGI_ACCUMULATOR(HcalDigiMixerHO);
DEFINE_DIGI_ACCUMULATOR(HcalDigiMixerHF);
DEFINE_DIGI_ACCUMULATOR(HcalDigiMixerZDC);
//DEFINE_DIGI_ACCUMULATOR(HcalDigiMixerHBHEUpgrade);
//DEFINE_DIGI_ACCUMULATOR(HcalDigiMixerHFUpgrade);
