#include "SimG4CMS/Muon/interface/MuonMe0FrameRotation.h"
#include "Geometry/MuonNumbering/interface/MuonDDDConstants.h"
#include "Geometry/MuonNumbering/interface/MuonBaseNumber.h"

#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"

//#define LOCAL_DEBUG

MuonMe0FrameRotation::MuonMe0FrameRotation(const DDCompactView& cpv) : MuonFrameRotation::MuonFrameRotation(cpv) {
  g4numbering     = new MuonG4Numbering(cpv);
  MuonDDDConstants muonConstants(cpv);
  int theLevelPart= muonConstants.getValue("level");
  theSectorLevel  = muonConstants.getValue("mg_sector")/theLevelPart;
#ifdef LOCAL_DEBUG
  std::cout << "MuonMe0FrameRotation: theSectorLevel " << theSectorLevel 
	    << std::endl;
#endif
}

MuonMe0FrameRotation::~MuonMe0FrameRotation() {
  delete g4numbering;
}

Local3DPoint MuonMe0FrameRotation::transformPoint(const Local3DPoint & point,const G4Step * aStep=0) const {
  if (!aStep) return Local3DPoint(0.,0.,0.);  

  //check if it is rotated
#ifdef LOCAL_DEBUG
  std::cout << "Position " << aStep->GetPreStepPoint()->GetPosition() << std::endl;
#endif
  MuonBaseNumber num = g4numbering->PhysicalVolumeToBaseNumber(aStep);
  bool rotated       = (num.getBaseNo(theSectorLevel)>=50);
#ifdef LOCAL_DEBUG
  std::cout << "MuonMe0FrameRotation num " << num.getBaseNo(theSectorLevel)
	    << " Rotation " << rotated << std::endl;
#endif
  if (rotated) {
    //    return Local3DPoint(-point.x(),point.z(),point.y());
    return Local3DPoint(point.x(),point.z(),-point.y());
  } else {
    return Local3DPoint(point.x(),point.z(),-point.y());
  }
}
