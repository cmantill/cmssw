#ifndef RecoTracker_TkSeedGenerator_SeedFromProtoTTrack_H
#define RecoTracker_TkSeedGenerator_SeedFromProtoTTrack_H

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/L1TrackingRecHit.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

namespace edm {
  class EventSetup;
}

class SeedFromProtoTTrack {
public:
  SeedFromProtoTTrack(const TTTrack< Ref_Phase2TrackerDigi_ >& proto, const edm::EventSetup&);
  ~SeedFromProtoTTrack() {}

  TrajectorySeed trajectorySeed() const;

  bool isValid() const { return theValid; }

private:
  void init(const TTTrack< Ref_Phase2TrackerDigi_ >& proto, const edm::EventSetup& es);

  PropagationDirection direction() const { return alongMomentum; }

  PTrajectoryStateOnDet const& trajectoryState() const { return thePTraj; }

  typedef edm::OwnVector<TrackingRecHit> RecHitContainer;
  const RecHitContainer& hits() const { return theHits; }

private:
  bool theValid;
  RecHitContainer theHits;
  PTrajectoryStateOnDet thePTraj;
};
#endif
