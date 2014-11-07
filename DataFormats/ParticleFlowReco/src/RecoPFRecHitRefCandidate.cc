#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidate.h"

using namespace reco;
using namespace std;

void RecoPFRecHitRefCandidate::init () {
  // get square of the cell
  const std::vector< math::XYZPoint >& corners = ref_->getCornersXYZ();
  std::vector< math::XYZVector > etaPhicorners;
  for (auto corner : corners) {
    etaPhicorners.push_back(math::XYZVector (corner.eta(), corner.phi(), 0));
  }
  if (etaPhicorners.size() >= 4) {
    etaPhiSquare_ = 0.5 * ((etaPhicorners[0]-etaPhicorners[1]).Cross (etaPhicorners[0]-etaPhicorners[3]).r() +
				  (etaPhicorners[2]-etaPhicorners[1]).Cross (etaPhicorners[2]-etaPhicorners[3]).r());
  }
}
