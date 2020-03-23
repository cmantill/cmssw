#include "RecoTracker/TkSeedGenerator/interface/SeedFromProtoTTrack.h"

#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Units/PhysicalConstants.h"

namespace {

template <class T>
inline T sqr(T t) {
  return t * t;
}

}

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
    //std::cout << "L1_TRH: " << hit.localPosition().x() << "," << hit.localPosition().y() << " and error : "
    //<< hit.localPositionError().xx() << "," << hit.localPositionError().yy() << std::endl;
    theHits.push_back(hit);
  }

  init(proto, es);
}

CurvilinearTrajectoryError SeedFromProtoTTrack::initialError(float sin2th, const TTTrack< Ref_Phase2TrackerDigi_ >& proto) const {
  // Set initial uncertainty on track parameters, 
  // Parameters associated to the 5D curvilinear covariance matrix: <BR>
  //  * <B> (qoverp, lambda, phi, dxy, dsz) </B><BR>
  //  * defined as:  <BR>
  //  *   <DT> qoverp = q / abs(p) = signed inverse of momentum [1/GeV] </DT> 
  //  *   <DT> lambda = pi/2 - polar angle at the given point </DT>
  //  *   <DT> phi = azimuth angle at the given point </DT>
  //  *   <DT> dxy = -vx*sin(phi) + vy*cos(phi) [cm] </DT>
  //  *   <DT> dsz = vz*cos(lambda) - (vx*cos(phi)+vy*sin(phi))*sin(lambda) [cm] </DT>

  CurvilinearTrajectoryError newError;  // zeroed
  auto& C = newError.matrix();

  // FIXME: minC00. Prevent apriori uncertainty in 1/P from being too small,
  // to avoid instabilities.
  // N.B. This parameter needs optimising ...
  //auto sin2th = sin2Theta;
  //auto minC00 = MinOneOverPtError*MinOneOverPtError;
  //auto sin2Theta = gtp.momentum().perp2() / gtp.momentum().mag2();

  //  AlgebraicSymMatrix51 mat;
  //double minC00 = 1.0;
  //mat[0] = (0.005 * trk_signedPt) * (0.05 * trk_signedPt);  // sigma^2(charge/abs_momentum)                                                                                    
  //mat[1] = std::abs(0.05 * proto.tanL());  // sigma^2(lambda)                                                                                                                  
  //mat[2] = 0.0002;    // sigma^2(phi)                                                                                                                                        
  //mat[3] = 20. * 20.;    // sigma^2(x_transverse))                                                                                                                           
  //mat[4] = 20. * 20.;    // sigma^2(y_transverse))                                                                                                                          
  //CurvilinearTrajectoryError error(mat);
  //AlgebraicSymMatrix C(5,1); CurvilinearTrajectoryError error(C);

  float maxC00 = 1.0;
  C[0][0] = std::max(sin2th, maxC00);
  C[1][1] = std::abs(atan(proto.tanL()));
  C[2][2] = std::abs(proto.phi()); 
  C[3][3] = std::abs(proto.d0());
  C[4][4] = std::abs(proto.z0());

  return newError;
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
  
  int dummy=0;
  float speedOfLightConverted = CLHEP::c_light/1.0E5; // B*c/2E11 - converts q/pt to track angle at some radius from beamline
  float trk_signedPt = speedOfLightConverted * mMagneticFieldStrength / proto.rInv(); // transverse curvature
  GlobalTrajectoryParameters gtp(proto.POCA(), proto.momentum(), trk_signedPt, dummy, &(*field)); 

  auto sin2Theta = gtp.momentum().perp2() / gtp.momentum().mag2();
  CurvilinearTrajectoryError error = initialError(sin2Theta,proto);
  FreeTrajectoryState fts(gtp, error);

  // lets try to figure out ID from stubs
  std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = proto.getStubRefs();
  DetId detIdStub = (stubRefs.at(stubRefs.size()-1)->clusterRef(0))->getDetId();
  TrajectoryStateOnSurface outerState = propagator->propagate(fts, tracker->idToDetUnit(detIdStub)->surface());

  const TrackingRecHit& lastHit = theHits.back();
  const Surface& surface = tracker->idToDet(lastHit.geographicalId())->surface();

  // these are the same
  //std::cout << "stublast detid " << tracker->idToDetUnit(detIdStub)->surface().position() << " lasthit " << surface.position() << std::endl;
  if (!outerState.isValid()) {
    //const Surface& surface = tracker->idToDetUnit(detIdStub)->surface();
    //edm::LogError("SeedFromProtoTTrack") << "SeedFromProtoTTrack was trying to create a seed from:\n"
    //<< fts ;
    //					 << fts << "\n propagating to stubdetid: " << std::hex << detIdStub.rawId()
    //<< std::dec << " surface position " << surface.position();
    //const Surface& surfaceOut = tracker->idToDetUnit(detIdStubOut)->surface();
    //edm::LogError("SeedFromProtoTTrack") << "SeedFromProtoTTrack surface outer position " << surfaceOut.position();
    std::cout << " SeedFromProtoTTrack NOT VALID " << std::endl;
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
