#include "Adjuster.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"

namespace edm {
namespace detail {
void doTheOffset(int bunchSpace, int bcr, std::vector<SimTrack>& simtracks, unsigned int evtNr, int vertexOffset) { 

  EncodedEventId id(bcr,evtNr);
  for (auto& item : simtracks) {
    item.setEventId(id);
    if (!item.noVertex()) {
      item.setVertexIndex(item.vertIndex() + vertexOffset);
    }
  }
}

void doTheOffset(int bunchSpace, int bcr, std::vector<SimVertex>& simvertices, unsigned int evtNr, int vertexOffset) { 

  int timeOffset = bcr * bunchSpace;
  EncodedEventId id(bcr,evtNr);
  for (auto& item : simvertices) {
    item.setEventId(id);
    item.setTof(item.position().t() + timeOffset);
  }
}

void doTheOffset(int bunchSpace, int bcr, std::vector<PSimHit>& simhits, unsigned int evtNr, int vertexOffset) { 

  int timeOffset = bcr * bunchSpace;
  EncodedEventId id(bcr,evtNr);
  for (auto& item : simhits) {
    item.setEventId(id);
    item.setTof(item.timeOfFlight() + timeOffset);
  }
}

void doTheOffset(int bunchSpace, int bcr, std::vector<PCaloHit>& calohits, unsigned int evtNr, int vertexOffset) { 

  int timeOffset = bcr * bunchSpace;
  EncodedEventId id(bcr,evtNr);
  for (auto& item : calohits) {
    item.setEventId(id);
    item.setTime(item.time() + timeOffset);
  }
}

  void doTheOffset(int bunchSpace, int bcr, HBHEDigiCollection& hcaldigi, unsigned int evtNr, int vertexOffset) { 
    
    EncodedEventId id(bcr,evtNr);
    for (auto& item : hcaldigi) {
      int size = item.size();
      if (bcr > 0) {
	for (int i = size; i > 0;) {
	  --i;
	  int newBX = i + bcr;
	  if (newBX < size) {
	    item.setSample (newBX, item.sample (i));
	  }
	} 
	for (int i = 0; i < bcr; ++i) {
	  item.setSample (i, 0);
	}
      }
      else if (bcr < 0) {
	for (int i = 0; i < size; ++i) {
	  int newBX = i + bcr;
	  if (newBX >= 0) {
	    item.setSample (newBX, item.sample (i));
	  }
	} 
	for (int i = size + bcr; i < size; ++i) {
	  item.setSample (i, 0);
	}
      }
    }
  }
}
}
