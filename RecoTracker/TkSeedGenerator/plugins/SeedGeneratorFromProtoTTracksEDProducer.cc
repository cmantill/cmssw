#include "SeedGeneratorFromProtoTTracksEDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "RecoTracker/TkSeedGenerator/interface/SeedFromProtoTTrack.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"
#include "SeedFromConsecutiveHitsCreator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoTracker/TkTrackingRegions/interface/GlobalTrackingRegion.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <vector>

using namespace edm;
using namespace reco;

template <class T>
T sqr(T t) {
  return t * t;
}
typedef SeedingHitSet::ConstRecHitPointer Hit;

void SeedGeneratorFromProtoTTracksEDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<InputTag>("InputCollection", InputTag("pixelTracks"));
  desc.add<InputTag>("InputVertexCollection", InputTag(""));
  desc.add<double>("originHalfLength", 1E9);
  desc.add<double>("originRadius", 1E9);
  desc.add<bool>("useProtoTrackKinematics", true);
  desc.add<bool>("useEventsWithNoVertex", true);
  desc.add<std::string>("TTRHBuilder", "TTRHBuilderWithoutAngle4PixelTriplets");
  desc.add<bool>("usePV", false);
  desc.add<bool>("includeFourthHit", false);

  edm::ParameterSetDescription psd0;
  psd0.add<std::string>("ComponentName", std::string("SeedFromConsecutiveHitsCreator"));
  psd0.add<std::string>("propagator", std::string("PropagatorWithMaterial"));
  psd0.add<double>("SeedMomentumForBOFF", 5.0);
  psd0.add<double>("OriginTransverseErrorMultiplier", 1.0);
  psd0.add<double>("MinOneOverPtError", 1.0);
  psd0.add<std::string>("magneticField", std::string(""));
  psd0.add<std::string>("TTRHBuilder", std::string("WithTrackAngle"));
  psd0.add<bool>("forceKinematicWithRegionDirection", false);
  desc.add<edm::ParameterSetDescription>("SeedCreatorPSet", psd0);

  descriptions.add("SeedGeneratorFromProtoTTracksEDProducer", desc);
}

SeedGeneratorFromProtoTTracksEDProducer::SeedGeneratorFromProtoTTracksEDProducer(const ParameterSet& cfg)
    : theConfig(cfg),
      originHalfLength(cfg.getParameter<double>("originHalfLength")),
      originRadius(cfg.getParameter<double>("originRadius")),
      useProtoTrackKinematics(cfg.getParameter<bool>("useProtoTrackKinematics")),
      useEventsWithNoVertex(cfg.getParameter<bool>("useEventsWithNoVertex")),
      builderName(cfg.getParameter<std::string>("TTRHBuilder")),
      usePV_(cfg.getParameter<bool>("usePV")),
      includeFourthHit_(cfg.getParameter<bool>("includeFourthHit")),
      theInputCollectionTag(consumes<std::vector< TTTrack< Ref_Phase2TrackerDigi_ >>>(cfg.getParameter<InputTag>("InputCollection"))),
      theInputVertexCollectionTag(
          consumes<reco::VertexCollection>(cfg.getParameter<InputTag>("InputVertexCollection"))) {
  produces<TrajectorySeedCollection>();
}

void SeedGeneratorFromProtoTTracksEDProducer::produce(edm::Event& ev, const edm::EventSetup& es) {
  auto result = std::make_unique<TrajectorySeedCollection>();
  Handle<std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > trks;
  ev.getByToken(theInputCollectionTag, trks);

  //edm::Handle<reco::VertexCollection> vertices;
  //bool foundVertices = ev.getByToken(theInputVertexCollectionTag, vertices);
  //const reco::VertexCollection & vertices = *(h_vertices.product());

  ///
  /// need optimization: all es stuff should go out of the loop
  ///
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator it;
  for ( it = trks->begin(); it != trks->end(); it++ ) {
    const TTTrack< Ref_Phase2TrackerDigi_ >& proto = (*it);

    // we don't check the compatibility with a primary vertex
    
    // use ProtoTTrack kinematics
    SeedFromProtoTTrack seedFromProtoTTrack(proto, es);
    if (seedFromProtoTTrack.isValid()){
      std::cout << "seed from ttrack is valid " << std::endl;
      (*result).push_back(seedFromProtoTTrack.trajectorySeed());
    }
  }
  
  ev.put(std::move(result));
}
