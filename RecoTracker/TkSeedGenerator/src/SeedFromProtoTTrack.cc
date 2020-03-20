#include "RecoTracker/TkSeedGenerator/interface/SeedFromProtoTTrack.h"

#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Units/PhysicalConstants.h"

SeedFromProtoTTrack::SeedFromProtoTTrack(const TTTrack< Ref_Phase2TrackerDigi_ >& proto, const edm::EventSetup& es) : theValid(true) {
  std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = proto.getStubRefs();

  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);

  for ( const auto& stubRef : stubRefs) {

    DetId detIdStub = tracker->idToDet( stubRef->clusterRef(0)->getDetId() )->geographicalId();
    
    auto coords = stubRef->clusterRef(0)->findAverageLocalCoordinatesCentered();  
    const GeomDet* theGeomDet = tracker->idToDet(detIdStub);
    LocalPoint lp = theGeomDet->topology().localPosition(coords);
    LocalError le;
    const auto genericDet = tracker->idToDetUnit( stubRef->clusterRef(0)->getDetId() );

    L1TrackingRecHit hit(lp, le, *genericDet, stubRef->clusterRef(0));
    theHits.push_back(hit);
  }

  init(proto, es);
}

// missing info from TTracks: charge? covariance/error?, detID of lasthit/stub?
void SeedFromProtoTTrack::init(const TTTrack< Ref_Phase2TrackerDigi_ >& proto, const edm::EventSetup& es) {
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);

  edm::ESHandle<Propagator> propagatorHandle;
  es.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", propagatorHandle);
  const Propagator* propagator = &(*propagatorHandle);

  edm::ESHandle<MagneticField> field;
  es.get<IdealMagneticFieldRecord>().get(field);  //fixme

  const MagneticField* theMagneticField = field.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();

  std::cout << "TTrack rinv " << proto.rInv() << std::endl;
  //int TTrackCharge = 1;
  //if(proto.rInv() < 0) TTrackCharge = -1;
  int dummy=0;
  float speedOfLightConverted = CLHEP::c_light/1.0E5; // B*c/2E11 - converts q/pt to track angle at some radius from beamline
  float trk_signedPt = speedOfLightConverted * mMagneticFieldStrength / proto.rInv(); // transverse curvature
  GlobalTrajectoryParameters gtp(proto.POCA(), proto.momentum(), trk_signedPt, dummy, &(*field)); 

  //Parametrization of the error matrix in the curvilinear frame.
  // *  This frame is tangent to the track at the point of definition,
  // *  with Z_T parallel to the track. X_T is in the global xy plane 
  // *  and points to the left when looking into the direction of the track,
  // *  and Y_T forms a right-handed frame with X_T and Z_T.
  // * 
  // *  The error along Z_T is therefore zero.
  // *  The parameters are <BR>
  // *    sigma^2( charge / abs_momentum) <BR>
  // *    sigma^2( lambda) <BR>
  // *    sigma^2( phi) <BR>
  // *    sigma^2( x_transverse)) <BR>
  // *    sigma^2( y_transverse)) <BR> <BR>
  // *  Please note that lambda and phi are defined in the global frame. Lambda is the helix
  // *  dip angle (pi/2 minus theta (polar angle)), while phi is the angle of 
  // *  inclination with the global x-axis in the transverse (global xy) plane.

  AlgebraicSymMatrix55 mat;
  double minC00 = 1.0;
  mat[0][0] = std::max((0.05 * trk_signedPt) * (0.05 * trk_signedPt),minC00);  // sigma^2(charge/abs_momentum)  
  mat[1][1] = std::abs(0.05 * proto.tanL());  // sigma^2(lambda)      
  mat[2][2] = 0.0002;    // sigma^2(phi)
  mat[3][3] = 20. * 20.;    // sigma^2(x_transverse))
  mat[4][4] = 20. * 20.;    // sigma^2(y_transverse))        

  CurvilinearTrajectoryError error(mat);

  FreeTrajectoryState fts(gtp, error);

  // lets try to figure out ID from stubs
  std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = proto.getStubRefs();
  // inner cluster of last stub?
  DetId detIdStub = (stubRefs.at(stubRefs.size()-1)->clusterRef(0))->getDetId();

  //DetId detIdStubOut = (stubRefs.at(stubRefs.size()-1)->clusterRef(1))->getDetId();

  TrajectoryStateOnSurface outerState =
    propagator->propagate(fts, tracker->idToDetUnit(detIdStub)->surface());

  if (!outerState.isValid()) {
    const Surface& surface = tracker->idToDetUnit(detIdStub)->surface();
    edm::LogError("SeedFromProtoTTrack") << "SeedFromProtoTTrack was trying to create a seed from:\n"
					 << fts ;
    //					 << fts << "\n propagating to stubdetid: " << std::hex << detIdStub.rawId()
    //<< std::dec << " surface position " << surface.position();
    //const Surface& surfaceOut = tracker->idToDetUnit(detIdStubOut)->surface();
    //edm::LogError("SeedFromProtoTTrack") << "SeedFromProtoTTrack surface outer position " << surfaceOut.position();
    theValid = false;
    return;
  }
  std::cout << " SeedFromProtoTTrack VALID " << std::endl;
  theValid = true;

  thePTraj = trajectoryStateTransform::persistentState(outerState, detIdStub.rawId());
}

TrajectorySeed SeedFromProtoTTrack::trajectorySeed() const {
  return TrajectorySeed(trajectoryState(), hits(), direction());
}
