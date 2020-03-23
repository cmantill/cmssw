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

struct HitLessByRadius {
  bool operator()(const Hit& h1, const Hit& h2) { return h1->globalPosition().perp2() < h2->globalPosition().perp2(); }
};

void SeedGeneratorFromProtoTTracksEDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<InputTag>("InputCollection", InputTag("pixelTracks"));

  edm::ParameterSetDescription psd0;
  psd0.add<std::string>("magneticField", std::string(""));

  descriptions.add("SeedGeneratorFromProtoTTracksEDProducer", desc);
}

SeedGeneratorFromProtoTTracksEDProducer::SeedGeneratorFromProtoTTracksEDProducer(const ParameterSet& cfg)
    : theConfig(cfg),
      theInputCollectionTag(consumes<std::vector< TTTrack< Ref_Phase2TrackerDigi_ >>>(cfg.getParameter<InputTag>("InputCollection"))),
  produces<TrajectorySeedCollection>();
}

void SeedGeneratorFromProtoTTracksEDProducer::produce(edm::Event& ev, const edm::EventSetup& es) {
  auto result = std::make_unique<TrajectorySeedCollection>();
  Handle<std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > trks;
  ev.getByToken(theInputCollectionTag, trks);

  edm::ESHandle<Chi2MeasurementEstimatorBase> estimatorH;
  iSetup.get<TrackingComponentsRecord>().get(theEstimatorName, estimatorH);

  edm::ESHandle<Propagator> propagatorAlongH;
  es.get<TrackingComponentsRecord>().get(thePropagatorName, propagatorAlongH);

  std::unique_ptr<Propagator> propagatorAlong = SetPropagationDirection(*propagatorAlongH, alongMomentum);

  //Loop over the L1's and make seeds for all of them:
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator it;
  for ( it = trks->begin(); it != trks->end(); it++ ) {
    const TTTrack< Ref_Phase2TrackerDigi_ >& l1 = (*it);

    std::unique_ptr<std::vector<TrajectorySeed> > out(new std::vector<TrajectorySeed>());
    std::cout << "SeedGeneratorFromTTracks::produce: L1 muon pT, eta, phi --> " << l1.momentum().mag() << " , " << l1.momentum().eta() << " , "
	      << l1.momentum().phi() << endl;

    FreeTrajectoryState fts = trajectoryStateTransform::initialFreeStateTTrack(*l1, magfieldH.product());
    dummyPlane->move(fts.position() - dummyPlane->position());
    TrajectoryStateOnSurface tsosAtIP = TrajectoryStateOnSurface(fts, *dummyPlane);
    std::cout << "SeedGeneratorFromTTracks::produce: Created TSOSatIP: " << tsosAtIP << std::endl;

    std::cout << "SeedGeneratorFromTTracks::findSeedsOnLayer: Start hitless" << endl;
    std::vector<GeometricSearchDet::DetWithState> dets;
    layer.compatibleDetsV(tsosAtIP, propagatorAlong, *estimatorH, dets);
    if (!dets.empty()) {
      auto const& detOnLayer = dets.front().first;
      auto const& tsosOnLayer = dets.front().second;
      LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: tsosOnLayer " << tsosOnLayer << endl;
      if (!tsosOnLayer.isValid()) {
	edm::LogInfo(theCategory) << "ERROR!: Hitless TSOS is not valid!";
      } else {
        // calculate SF from L2 (only once -- if needed)
        if (!analysedL2 && adjustErrorsDynamicallyForHitless_) {
          errorSFHitless = calculateSFFromL2(l2);
          analysedL2 = true;
        }

        dets.front().second.rescaleError(errorSFHitless);
        PTrajectoryStateOnDet const& ptsod =
	  trajectoryStateTransform::persistentState(tsosOnLayer, detOnLayer->geographicalId().rawId());
	TrajectorySeed::recHitContainer rHC;
        out->push_back(TrajectorySeed(ptsod, rHC, oppositeToMomentum));
        LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: TSOD (Hitless) done " << endl;
        numSeedsMade++;
      }
    }

    // if (useProtoTrackKinematics) {
    //   // use ProtoTTrack kinematics
    //   SeedFromProtoTTrack seedFromProtoTTrack(proto, es);
    //   if (seedFromProtoTTrack.isValid()){
    // 	std::cout << "seed from ttrack is valid " << std::endl;
    // 	(*result).push_back(seedFromProtoTTrack.trajectorySeed());
    //   }

    //   edm::ESHandle<TrackerGeometry> tracker;
    //   es.get<TrackerDigiGeometryRecord>().get(tracker);

    //   edm::ESHandle<Propagator> propagatorHandle;
    //   es.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", propagatorHandle);
    //   const Propagator* propagator = &(*propagatorHandle);

    //   edm::ESHandle<MagneticField> field;
    //   es.get<IdealMagneticFieldRecord>().get(field);  //fixme                                                                                                                            
    //   const MagneticField* theMagneticField = field.product();
    //   double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();

    //   int dummy=0;
    //   float speedOfLightConverted = CLHEP::c_light/1.0E5; // B*c/2E11 - converts q/pt to track angle at some radius from beamline                                                        
    //   float trk_signedPt = speedOfLightConverted * mMagneticFieldStrength / proto.rInv(); // transverse curvature                                                                        
    //   GlobalTrajectoryParameters gtp(proto.POCA(), proto.momentum(), trk_signedPt, dummy, &(*field));

    //   auto sin2Theta = gtp.momentum().perp2() / gtp.momentum().mag2();
    //   CurvilinearTrajectoryError error = initialError(sin2Theta,proto);
    //   FreeTrajectoryState fts(gtp, error);

    //   // lets try to figure out ID from stubs                                                                                                                                            
    //   std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = proto.getStubRefs();
    //   DetId detIdStub = (stubRefs.at(stubRefs.size()-1)->clusterRef(0))->getDetId();
    //   TrajectoryStateOnSurface outerState = propagator->propagate(fts, tracker->idToDetUnit(detIdStub)->surface());

    //   const TrackingRecHit& lastHit = theHits.back();
    //   const Surface& surface = tracker->idToDet(lastHit.geographicalId())->surface();


    // }
    // else {
    //   std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = proto.getStubRefs();

    //   edm::ESHandle<TrackerGeometry> tracker;
    //   es.get<TrackerDigiGeometryRecord>().get(tracker);

    //   //edm::OwnVector<TrackingRecHit> hits;
    //   std::vector<SeedingHitSet::ConstRecHitPointer> hits;
    //   //ConstRecHitContainer hits;

    //   for ( const auto& stubRef : stubRefs) {

    // 	DetId detIdStub = tracker->idToDet( stubRef->clusterRef(0)->getDetId() )->geographicalId();

    // 	auto coords = stubRef->clusterRef(0)->findAverageLocalCoordinatesCentered();
    // 	const GeomDet* theGeomDet = tracker->idToDet(detIdStub);
    // 	LocalPoint lp = theGeomDet->topology().localPosition(coords);
    // 	LocalError le;
    // 	const auto genericDet = tracker->idToDetUnit( stubRef->clusterRef(0)->getDetId() );

    //     L1TrackingRecHit hit(lp, le, *genericDet, stubRef->clusterRef(0));
    // 	if(hit.isValid()) 
    // 	  hits.push_back(hit.clone());
    //   }
      
    //   //std::sort(hits.begin(), hits.end(), HitLessByRadius());

    //   if (hits.size() > 1) {
    // 	double mom_perp =
    // 	  sqrt(proto.momentum().x() * proto.momentum().x() + proto.momentum().y() * proto.momentum().y());
    // 	GlobalPoint vtx(proto.POCA().x(), proto.POCA().y(), proto.POCA().z());
    // 	GlobalTrackingRegion region(mom_perp, vtx, 0.2, 0.2);
	
    //     edm::ParameterSet seedCreatorPSet = theConfig.getParameter<edm::ParameterSet>("SeedCreatorPSet");
    // 	SeedFromConsecutiveHitsCreator seedCreator(seedCreatorPSet);
    // 	seedCreator.init(region, es, nullptr);
    // 	SeedingHitSet seedHits(hits[0],hits[1],hits[2]);
    // 	seedCreator.makeSeed(
    // 			     *result,
    // 			     seedHits
    // 			     //SeedingHitSet(hits[0],
    // 			     //hits[1]
    // 			     //hits.size() > 2 ? hits[2] : SeedingHitSet::nullPtr(),
    // 			     //(includeFourthHit_ && hits.size() > 3) ? hits[3] : SeedingHitSet::nullPtr()));
    // 			     //)
    // 			     );
    //   }
    // }
  } // end loop over L1Tracks
  
  ev.put(std::move(result));
}
