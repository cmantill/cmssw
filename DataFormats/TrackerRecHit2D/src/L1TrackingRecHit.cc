#include "DataFormats/TrackerRecHit2D/interface/L1TrackingRecHit.h"

void L1TrackingRecHit::getKfComponents(KfComponentsHolder& holder) const {
  getKfComponents2D(holder);
}
