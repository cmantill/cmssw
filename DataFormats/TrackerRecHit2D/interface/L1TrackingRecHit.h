#ifndef DataFormats_TrackerRecHit2D_L1TrackingRecHit_H
#define DataFormats_TrackerRecHit2D_L1TrackingRecHit_H

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/OwnVector.h"

#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "TkCloner.h"

class L1TrackingRecHit : public TrackerSingleRecHit {
 public:

  L1TrackingRecHit() : TrackerSingleRecHit() {}

  L1TrackingRecHit(const LocalPoint& p, const LocalError& e,
                    //const GeomDet& idet, const TTClusterRef& objref) :
		    const GeomDet& idet, const ClusterL1Ref& objref):
//edm::Ref< edmNew::DetSetVector< TTCluster< T > >, TTCluster< T > >
  TrackerSingleRecHit(p, e, idet, trackerHitRTTI::L1Tracking, objref)
    {}
  //TrackerSingleRecHit(pos,err, idet, clus) {}
				
  L1TrackingRecHit * clone() const override {return new L1TrackingRecHit( * this); }

  // things to specialize from BaseTrackerRecHit     
  bool isPhase2() const final { return true; }  
  void getKfComponents( KfComponentsHolder & holder ) const final;

  int dimension() const override {return 2;}
};

typedef edmNew::DetSetVector<L1TrackingRecHit> L1TrackingDetSetVector;
typedef edm::OwnVector<L1TrackingRecHit> L1TrackingOwnVector;

#endif

