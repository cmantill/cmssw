#ifndef RecoTracker_TkSeedGenerator_SeedGeneratorFromProtoTTracksEDProducer_H
#define RecoTracker_TkSeedGenerator_SeedGeneratorFromProtoTTracksEDProducer_H

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

class dso_hidden SeedGeneratorFromProtoTTracksEDProducer : public edm::stream::EDProducer<> {
public:
  SeedGeneratorFromProtoTTracksEDProducer(const edm::ParameterSet& cfg);
  ~SeedGeneratorFromProtoTTracksEDProducer() override {}
  void produce(edm::Event& ev, const edm::EventSetup& es) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const edm::ParameterSet theConfig;
  const double originHalfLength;
  const double originRadius;
  const bool useProtoTrackKinematics;
  const bool useEventsWithNoVertex;
  const std::string builderName;
  const bool usePV_;
  const bool includeFourthHit_;
  const edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > theInputCollectionTag;
  const edm::EDGetTokenT<reco::VertexCollection> theInputVertexCollectionTag;
};
#endif
