// -*- C++ -*-
//
// Package:    HGCGeometryTester
// Class:      HGCGeometryTester
// 
/**\class HGCGeometryTester HGCGeometryTester.cc test/HGCGeometryTester.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Sunanda Banerjee
//         Created:  Mon 2014/02/07
// $Id: HGCGeometryTester.cc,v 1.0 2014/02/07 14:06:07 sunanda Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDExpandedView.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "SimG4CMS/Calo/interface/HGCNumberingScheme.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "CoralBase/Exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/RotationInterfaces.h"

typedef CLHEP::Hep3Vector G4ThreeVector;

#include "TRandom2.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

TRandom2 rndm; 

const double k_ScaleFromDDDtoGeant = 0.1;

class GeometryTesterHistograms {
public:
  struct LayerHists {
    LayerHists (const std::string& fSuffix) 
      : suffix (fSuffix) 
    {
      std::string histName ("faultPoints"); histName += fSuffix;
      faultPoints = TH2D (histName.c_str(), histName.c_str(), 40, -1, 1, 40, -1,1);
      (histName = ("allPoints")) += fSuffix;
      allPoints = TH2D (histName.c_str(), histName.c_str(), 40, -1, 1, 40, -1,1);
      (histName = ("faultPointsWeighted")) += fSuffix;
      faultPointsWeighted = TH2D (histName.c_str(), histName.c_str(), 40, -1, 1, 40, -1,1);
      (histName = ("faultDistance")) += fSuffix;
      faultDistance = TH1D (histName.c_str(), histName.c_str(), 2000, -20, 20);
    }
    static void Delete (LayerHists* hists) {delete hists;}

    void fillFault (double dr, const double* point) {
      faultPoints.Fill (point[0], point[1]);
      faultPointsWeighted.Fill (point[0], point[1], dr);
      faultDistance.Fill (dr);
    }

    void fillAll (const double* point) {
      allPoints.Fill (point[0], point[1]);
    }

    void write () {
      faultPoints.Write();
      allPoints.Write();
      faultPointsWeighted.Write();
      faultDistance.Write();
    }
    
    TH2D faultPoints;
    TH2D allPoints;
    TH2D faultPointsWeighted;
    TH1D faultDistance;
    std::string suffix;
  };
  struct SubdetHists {
    SubdetHists (const std::string& fSuffix)
      : suffix (fSuffix) 
    {
      std::string histName ("faultPoints"); histName += fSuffix;
      faultPoints = TH2D (histName.c_str(), histName.c_str(), 40, -1, 1, 40, -1,1);
      (histName = ("allPoints")) += fSuffix;
      allPoints = TH2D (histName.c_str(), histName.c_str(), 40, -1, 1, 40, -1,1);
      (histName = ("faultPointsWeighted")) += fSuffix;
      faultPointsWeighted = TH2D (histName.c_str(), histName.c_str(), 40, -1, 1, 40, -1,1);
      (histName = ("faultDistance")) += fSuffix;
      faultDistance = TH1D (histName.c_str(), histName.c_str(), 2000, -20, 20);
    }
    ~SubdetHists () {
      std::for_each (layerHists.begin(), layerHists.end(), LayerHists::Delete);
    }
    static void Delete (SubdetHists* hists) {delete hists;}

    void fillAll (const double* point, int layer) {
      if (layer >= int(layerHists.size())) {
	layerHists.resize (layer+1, 0);
      }
      if (layerHists[layer] == 0) {
	char layerSuffix [1024];
	sprintf (layerSuffix, "%s_%d", suffix.c_str(), layer);
	layerHists[layer] = new LayerHists (layerSuffix);
      }
      layerHists[layer]->fillAll (point);
      allPoints.Fill (point[0], point[1]);
    }

    void fillFault (double dr, const double* point, int layer) {
      if (layer >= int(layerHists.size())) {
	layerHists.resize (layer+1, 0);
      }
      if (layerHists[layer] == 0) {
	char layerSuffix [1024];
	sprintf (layerSuffix, "%s_%d", suffix.c_str(), layer);
	layerHists[layer] = new LayerHists (layerSuffix);
      }
      layerHists[layer]->fillFault (dr, point);
      faultPoints.Fill (point[0], point[1]);
      faultPointsWeighted.Fill (point[0], point[1], dr);
      faultDistance.Fill (dr);
    }

  void write () {
    for (size_t i = 0; i < layerHists.size(); ++i) {
      if (layerHists[i]) layerHists[i]->write();
    }
    faultPoints.Write();
    allPoints.Write();
    faultPointsWeighted.Write();
    faultDistance.Write();
  }
    

    TH2D faultPoints;
    TH2D allPoints;
    TH2D faultPointsWeighted;
    TH1D faultDistance;
    std::string suffix;
    std::vector<LayerHists*> layerHists;
  };

  GeometryTesterHistograms () {};

  ~GeometryTesterHistograms () {
    std::for_each (subdetHists.begin(), subdetHists.end(), SubdetHists::Delete);
  }

  void fillFault (double dr, const double* point, ForwardSubdetector subdet, int layer) {
    int index = (subdet == HGCEE) ? 0 :
      (subdet == HGCHEF) ? 1 :
      (subdet == HGCHEB) ? 2 : -1;
    if (index < 0) return;
    if (int(subdetHists.size()) <= index) subdetHists.resize (index+1, 0);
    if (subdetHists[index] == 0) {
      const char* subDetSuffix [] = {"_EE", "_HEF", "_HEB"};
      subdetHists[index] = new SubdetHists (std::string (subDetSuffix[index]));
    }
    subdetHists[index]->fillFault (dr, point, layer);
  }

  void fillAll (const double* point, ForwardSubdetector subdet, int layer) {
    int index = (subdet == HGCEE) ? 0 :
      (subdet == HGCHEF) ? 1 :
      (subdet == HGCHEB) ? 2 : -1;
    if (index < 0) return;
    if (int(subdetHists.size()) <= index) subdetHists.resize (index+1, 0);
    if (subdetHists[index] == 0) {
      const char* subDetSuffix [] = {"_EE", "_HEF", "_HEB"};
      subdetHists[index] = new SubdetHists (std::string (subDetSuffix[index]));
    }
    subdetHists[index]->fillAll (point, layer);
  }

  void write () {
    for (size_t i = 0; i < subdetHists.size(); ++i) {
      if (subdetHists[i]) subdetHists[i]->write();
    }
  }

private:
  std::vector<SubdetHists*> subdetHists;
};

//
// class decleration
//

class HGCGeometryTester : public edm::EDAnalyzer {
private:
  GeometryTesterHistograms histograms;

public:
  explicit HGCGeometryTester( const edm::ParameterSet& );
  ~HGCGeometryTester();

  
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  virtual void endJob();
private:
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
HGCGeometryTester::HGCGeometryTester(const edm::ParameterSet& ) {}


HGCGeometryTester::~HGCGeometryTester() {}


//
// member functions
//

// ------------ method called to produce the data  ------------
void HGCGeometryTester::analyze( const edm::Event& iEvent, 
				 const edm::EventSetup& iSetup ) {

  edm::ESTransientHandle<DDCompactView> pDD;
  iSetup.get<IdealGeometryRecord>().get( pDD );

  // get numbering, geometry, etc
  
  std::vector <HGCNumberingScheme*> numberingScheme (size_t(ForwardSubdetector::FastTime)+1, 0);
  std::vector <edm::ESHandle<HGCalGeometry> > geometry (size_t(ForwardSubdetector::FastTime)+1);
  std::string nameX = "HGCalEESensitive";
  numberingScheme[size_t(ForwardSubdetector::HGCEE)] = 
    new HGCNumberingScheme (*pDD, nameX, true);
  iSetup.get<IdealGeometryRecord>().get(nameX, geometry[size_t(ForwardSubdetector::HGCEE)]);
  nameX = "HGCalHESiliconSensitive";
  numberingScheme[size_t(ForwardSubdetector::HGCHEF)] = 
    new HGCNumberingScheme (*pDD, nameX, true);
  iSetup.get<IdealGeometryRecord>().get(nameX, geometry[size_t(ForwardSubdetector::HGCHEF)]);
  nameX = "HGCalHEScintillatorSensitive";
  numberingScheme[size_t(ForwardSubdetector::HGCHEB)] = 
    new HGCNumberingScheme (*pDD, nameX, true);
  iSetup.get<IdealGeometryRecord>().get(nameX, geometry[size_t(ForwardSubdetector::HGCHEB)]);

  //parse the DD for sensitive volumes
  DDExpandedView eview(*pDD);
  int counterAll = 0;
  int counterFault = 0;
  do {
    const DDLogicalPart &logPart=eview.logicalPart();
    std::string name=logPart.name();

    if ((name.find("HGCal") != std::string::npos) &&
	(name.find("Sensitive") != std::string::npos)) {
    
    ForwardSubdetector subdet = ForwardSubdetector::ForwardEmpty;
      if (name.find("HGCalEESensitive")!=std::string::npos) {
        subdet = ForwardSubdetector::HGCEE;
      } else if (name.find("HGCalHESiliconSensitive")!=std::string::npos) {
        subdet = ForwardSubdetector::HGCHEF;
      } else if (name.find("HGCalHEScintillatorSensitive")!=std::string::npos) {
        subdet = ForwardSubdetector::HGCHEB;
      }


      if (subdet != ForwardEmpty) {
	std::vector<double> solidPar=eview.logicalPart().solid().parameters();
	int layer = eview.copyNumbers()[eview.copyNumbers().size()-1];
	int section = eview.copyNumbers()[eview.copyNumbers().size()-2];
// 	std::cout << name << " Section:Layer " << section <<':' << layer<<"->"
// 		  << std::setprecision(3)<<"dz/y/xl/xt/alpha"<<'\t'
// 		  <<solidPar[0]<<'\t'
// 		  <<solidPar[3]<<'\t'<<solidPar[4]<<'\t'<<solidPar[5]<<'\t'
// 		  <<solidPar[6]
// 		  << std::setprecision(6) << std::endl;
// 	std::cout << name << " Layer " << layer << " " << solidPar[3] 
// 		  << " " << solidPar[4] << " " << solidPar[5] << std::endl;
	
	// get global positions
	DD3Vector x, y, z;
	eview.rotation().GetComponents( x, y, z ) ;
	const CLHEP::HepRep3x3 rotation ( x.X(), y.X(), z.X(),
					  x.Y(), y.Y(), z.Y(),
					  x.Z(), y.Z(), z.Z() );
	const CLHEP::HepRotation hr ( rotation );
	const CLHEP::Hep3Vector h3v ( eview.translation().X(),
				      eview.translation().Y(),
				      eview.translation().Z()  ) ;
	const HepGeom::Transform3D ht3d ( hr,          // only scale translation
					  k_ScaleFromDDDtoGeant*h3v ) ;    
	
	CLHEP::Hep3Vector center = ht3d * HepGeom::Point3D<double> (0,0,0);
// 	std::cout << "Center location: " << center.x()<<':'<<center.y()<<':'<<center.z()
// 		  << " R:phi " << center.perp()<<':'<< int(floor(0.5+(center.phi()/3.1415)*180));
	
	double dz = eview.logicalPart().solid().parameters()[0]*k_ScaleFromDDDtoGeant;
	double dy = eview.logicalPart().solid().parameters()[3]*k_ScaleFromDDDtoGeant;
	double x1 = eview.logicalPart().solid().parameters()[4]*k_ScaleFromDDDtoGeant;
	double x2 = eview.logicalPart().solid().parameters()[5]*k_ScaleFromDDDtoGeant;
	double tanAlpha = tan(eview.logicalPart().solid().parameters()[6]);
	
// 	std::cout << " Corners:";
	
// 	HepGeom::Point3D<double> corner = 
// 	  ht3d * HepGeom::Point3D<double>(-x2+dy*tanAlpha, dy, 0);
// 	std::cout << ' ' << corner.x() << ':' << corner.y(); 
// 	corner = 
// 	  ht3d * HepGeom::Point3D<double>(x2+dy*tanAlpha, dy, 0);
// 	std::cout << ' ' << corner.x() << ':' << corner.y(); 
// 	corner = 
// 	  ht3d * HepGeom::Point3D<double>(-x1-dy*tanAlpha, -dy, 0);
// 	std::cout << ' ' << corner.x() << ':' << corner.y(); 
// 	corner = 
// 	  ht3d * HepGeom::Point3D<double>(x1-dy*tanAlpha, -dy, 0);
// 	std::cout << ' ' << corner.x() << ':' << corner.y();
// 	std::cout << std::endl;
	
	// get random point
	double rr [3];
	HepGeom::Point3D<double> localPos;
	for (int iRndm = 0; iRndm < 100; ++iRndm) {
	while (1) {
	  rr[2] = rndm.Uniform (-1, 1);
	  double rz = rr[2]*dz;
	  rr[1] = rndm.Uniform (-1, 1);
	  double ry = rr[1]*dy;
	  rr[0] = rndm.Uniform (-1, 1);
	  //	  double dx = 0.5*(x1+x2)+0.5*(x2-x1)*ry/dy;
	  double rx = rr[0]*std::max(x1,x2);
	  double drx = (0.5*(x1+x2)+0.5*(x2-x1)*rr[1]);
	  localPos = HepGeom::Point3D<double> (rx,ry,rz);
	  if (fabs(rx) < drx) {
	    if (tanAlpha != 0) {
	      rr[0] = 0.5*rr[0] + 0.5*(tanAlpha>0?1:-1)*drx/std::max(x1,x2); // compress rndm for HESci
	      rx += ry * tanAlpha; // offset due to tilt
	    }
	    break;
	  }
	}

	HepGeom::Point3D<double> rPoint = ht3d * localPos;
	GlobalPoint rPointG (rPoint.x(),rPoint.y(),rPoint.z());
	counterAll++;
	// convert to ID
 	int iz = rPoint.z() > 0 ? 1 : -1;
	G4ThreeVector localPosG (localPos.x()/k_ScaleFromDDDtoGeant, 
				 localPos.y()/k_ScaleFromDDDtoGeant, 
				 localPos.z()/k_ScaleFromDDDtoGeant);
 	HGCalDetId simId (numberingScheme[size_t(subdet)]->
			  getUnitID(subdet, layer, section, iz, 
				    localPosG));  // must be const G4ThreeVector&
	// convert SimId to RecoId
	const HGCalTopology &topo=geometry[size_t(subdet)]->topology();
	const HGCalDDDConstants &dddConst=topo.dddConstants();
	std::pair<int,int> recoLayerCell=dddConst.simToReco(simId.cell(), simId.layer(),topo.detectorType());
	int recoCell  = recoLayerCell.first;
	int recoLayer = recoLayerCell.second;
	if(recoLayer<0) {
	  //	  std::cout << "====> fail dddConst.simToReco cell:" << simId << " for detector " << int(subdet) << ':' << topo.subDetector() << std::endl;
	  continue;
	}
	histograms.fillAll (rr, subdet, recoLayer);

	if (recoCell != simId.cell()) {
	  std::cout << "cell change: " << simId.cell() << "->" << recoCell << std::endl;
	}
	
	//assign the RECO DetId
	DetId id( subdet == HGCEE ?
		  (uint32_t)HGCEEDetId(subdet,simId.zside(),recoLayer,simId.sector(),simId.subsector(),recoCell) :
		  (uint32_t)HGCHEDetId(subdet,simId.zside(),recoLayer,simId.sector(),simId.subsector(),recoCell) );
	// Reco Geometry
	GlobalPoint pos = geometry[size_t(subdet)]->getPosition(id);
	HGCalGeometry::CornersVec corners = geometry[size_t(subdet)]->getCorners(id);
	double cellSize = (corners[0]-corners[1]).mag();
	DetId idClose = geometry[size_t(subdet)]->getClosestCell(rPointG);
	//	std::cout << "Closest cell to "<<rPoint.x()<<':'<<rPoint.y()<<':'<<rPoint.z()<<" is "<<std::hex<<idClose.rawId()<<std::dec<<' ';
// 	if (subdet == HGCEE) {
// 	  std::cout << HGCEEDetId (idClose);
// 	}
// 	else {
// 	  std::cout << HGCHEDetId (idClose);
// 	}
	GlobalPoint posClose = geometry[size_t(subdet)]->getPosition(idClose);
	double offsetReco = 0.5 * (rPointG-pos).mag() / cellSize;
	double offsetClose = 0.5 * (rPointG-posClose).mag() / cellSize;
	if (idClose != id) {
	  counterFault++;
	  //	  if (subdet != HGCHEB) {
	  if (subdet == HGCHEB && (offsetReco - offsetClose) > 8) {
	  std::cout << "------------------------------------" << std::endl;
	  std::cout << " Detector " << name << " corners:";
	  HepGeom::Point3D<double> corner = 
	    ht3d * HepGeom::Point3D<double>(-x2+dy*tanAlpha, dy, 0);
	  std::cout << ' ' << corner.x() << ':' << corner.y(); 
	  corner = 
	    ht3d * HepGeom::Point3D<double>(x2+dy*tanAlpha, dy, 0);
	  std::cout << ' ' << corner.x() << ':' << corner.y(); 
	  corner = 
	    ht3d * HepGeom::Point3D<double>(-x1-dy*tanAlpha, -dy, 0);
	  std::cout << ' ' << corner.x() << ':' << corner.y(); 
	  corner = 
	    ht3d * HepGeom::Point3D<double>(x1-dy*tanAlpha, -dy, 0);
	  std::cout << ' ' << corner.x() << ':' << corner.y();
	  std::cout << std::endl;
	  std::cout << "Random point relative: " << rr[0]<<':'<<rr[1]<<':'<<rr[2]<<std::endl; 
	  std::cout << "Random point local: " 
		    <<  localPos.x()<<':'<<localPos.y()<<':'<<localPos.z()
		    << " global: " 
		    << rPoint.x()<<':'<<rPoint.y()<<':'<<rPoint.z()
		    << " HGCalDetId: " << simId
		    << std::endl;
	  
	  if (subdet == HGCEE) {
	    std::cout << "RecoId: " << HGCEEDetId (id) << std::endl
		      << "closest:" << HGCEEDetId (idClose) << std::endl; 
	  }
	  else {
	    std::cout << "RecoId: " << HGCHEDetId (id) << std::endl 
		      << "closest:" << HGCHEDetId (idClose) << std::endl; 
	  }
	  
	  std::cout << "Reco cell geometry center: "<<pos.x()<<':'<<pos.y()<<':'<<pos.z();
	  for (unsigned i = 0; i < 4; ++i) {
	    std::cout << " "<<i<<">"<<corners[i].x()<<':'<<corners[i].y()<<':'<<corners[i].z();
	  }
	  std::cout << std::endl;
	  
	  std::cout << " Closest cell center: "<<posClose.x()<<':'<<posClose.y()<<':'<<posClose.z()<<std::endl;
	  std::cout << "distance point<->recoCell: "<<(rPointG-pos).mag()<<"("<<int(offsetReco*100)<<"%) point-closestCell: "<<(rPointG-posClose).mag()<<"("<<int(offsetClose*100)<<"%)"<<std::endl;;
	  
	  std::cout << "===> Reco cell is not a closest cell! deltaR="<<(rPointG-pos).mag()-(rPointG-posClose).mag()<<"("<<int((offsetReco-offsetClose)*100) << "%)"
		    << " fauls/totals: " << counterFault<<'/'<<counterAll<<'='<<double(counterFault)/counterAll
		    <<std::endl;
	  }
	  histograms.fillFault (offsetReco-offsetClose, rr, subdet, recoLayer);
	}
	}
      }
    }
  } while(eview.next() );
}

void HGCGeometryTester::endJob() {
  TFile outHists ("HGCGeometryTester.root", "RECREATE");
  histograms.write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCGeometryTester);
