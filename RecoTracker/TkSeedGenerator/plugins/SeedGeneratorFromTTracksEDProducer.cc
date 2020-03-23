#include "SeedGeneratorFromTTracksEDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoTracker/TkTrackingRegions/interface/GlobalTrackingRegion.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <vector>

using namespace edm;
using namespace reco;

void SeedGeneratorFromTTracksEDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<InputTag>("InputCollection", InputTag("pixelTracks"));
  desc.add<std::string>("estimator", "hltESPChi2MeasurementEstimator100");
  desc.add<std::string>("propagatorName", "PropagatorWithMaterial");
  desc.add<edm::InputTag>("MeasurementTrackerEvent", edm::InputTag("hltSiStripClusters"));
  desc.add<double>("maxEtaForTOB", 1.2);
  desc.add<double>("minEtaForTEC", 0.8);
  descriptions.add("SeedGeneratorFromTTracksEDProducer", desc);
}

SeedGeneratorFromTTracksEDProducer::SeedGeneratorFromTTracksEDProducer(const ParameterSet& cfg)
    : theConfig(cfg),
      theInputCollectionTag(consumes<std::vector< TTTrack< Ref_Phase2TrackerDigi_ >>>(cfg.getParameter<InputTag>("InputCollection"))),
      theEstimatorName(cfg.getParameter<std::string>("estimator")),
      thePropagatorName(cfg.getParameter<std::string>("propagatorName")),
      theMeasurementTrackerTag(consumes<MeasurementTrackerEvent>(cfg.getParameter<edm::InputTag>("MeasurementTrackerEvent"))),
      theMinEtaForTEC(cfg.getParameter<double>("minEtaForTEC")),
      theMaxEtaForTOB(cfg.getParameter<double>("maxEtaForTOB"))
{
  //produces<std::vector<TrajectorySeed> >();
  produces<TrajectorySeedCollection>();
}

void SeedGeneratorFromTTracksEDProducer::findSeedsOnLayer(const GeometricSearchDet& layer,
							  const TrajectoryStateOnSurface& tsosAtIP,
							  const Propagator& propagatorAlong,
							  const TTTrack< Ref_Phase2TrackerDigi_ >& l1,
							  edm::ESHandle<Chi2MeasurementEstimatorBase>& estimatorH,
							  unsigned int& numSeedsMade,
							  std::unique_ptr<std::vector<TrajectorySeed> >& out) const {
  
  std::cout << "SeedGeneratorFromTTracks::findSeedsOnLayer: Start hitless" << std::endl;
  std::vector<GeometricSearchDet::DetWithState> dets;
  layer.compatibleDetsV(tsosAtIP, propagatorAlong, *estimatorH, dets);

  if (!dets.empty()) {
    auto const& detOnLayer = dets.front().first;
    auto const& tsosOnLayer = dets.front().second;
    std::cout << "SeedGeneratorFromTTracks::findSeedsOnLayer: tsosOnLayer " << tsosOnLayer << std::endl;
    if (!tsosOnLayer.isValid()) {
      std::cout << "ERROR!: Hitless TSOS is not valid! \n";
    } else {
      //dets.front().second.rescaleError(errorSFHitless);
      PTrajectoryStateOnDet const& ptsod = trajectoryStateTransform::persistentState(tsosOnLayer, detOnLayer->geographicalId().rawId());
      TrajectorySeed::recHitContainer rHC;
      out->push_back(TrajectorySeed(ptsod, rHC, oppositeToMomentum));
      std::cout << "SeedGeneratorFromTTracks::findSeedsOnLayer: TSOD (Hitless) done " << std::endl;
      numSeedsMade++;
    }
  }
  else{
    std::cout << "SeedGeneratorFromTTracks::findSeedsOnLayer: TSOD (Hitless)  no dets " << std::endl;
  }

}

void SeedGeneratorFromTTracksEDProducer::produce(edm::Event& ev, const edm::EventSetup& es) {
  std::cout << "SeedGeneratorFromTTracks::produce start"  << std::endl;
  std::unique_ptr<std::vector<TrajectorySeed> > result(new std::vector<TrajectorySeed>());
  
  // TTrack Collection
  Handle<std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > trks;
  ev.getByToken(theInputCollectionTag, trks);

  // Trk Geometry
  edm::ESHandle<TrackerGeometry> tmpTkGeometryH;
  es.get<TrackerDigiGeometryRecord>().get(tmpTkGeometryH);

  // Mag field
  edm::ESHandle<MagneticField> magfieldH;
  es.get<IdealMagneticFieldRecord>().get(magfieldH);

  // Estimator
  edm::ESHandle<Chi2MeasurementEstimatorBase> estimatorH;
  es.get<TrackingComponentsRecord>().get(theEstimatorName, estimatorH);

  // Get Propagators
  edm::ESHandle<Propagator> propagatorAlongH;
  es.get<TrackingComponentsRecord>().get(thePropagatorName, propagatorAlongH);
  std::unique_ptr<Propagator> propagatorAlong = SetPropagationDirection(*propagatorAlongH, alongMomentum);

  // Get vector of Detector layers 
  edm::Handle<MeasurementTrackerEvent> measurementTrackerH;
  ev.getByToken(theMeasurementTrackerTag, measurementTrackerH);
  std::vector<BarrelDetLayer const*> const& tob = measurementTrackerH->geometricSearchTracker()->tobLayers();
  std::vector<ForwardDetLayer const*> const& tecPositive =
    tmpTkGeometryH->isThere(GeomDetEnumerators::P2OTEC)
    ? measurementTrackerH->geometricSearchTracker()->posTidLayers()
    : measurementTrackerH->geometricSearchTracker()->posTecLayers();
  std::vector<ForwardDetLayer const*> const& tecNegative =
    tmpTkGeometryH->isThere(GeomDetEnumerators::P2OTEC)
    ? measurementTrackerH->geometricSearchTracker()->negTidLayers()
    : measurementTrackerH->geometricSearchTracker()->negTecLayers();

  /// Surface used to make a TSOS at the PCA to the beamline
  Plane::PlanePointer dummyPlane = Plane::build(Plane::PositionType(), Plane::RotationType());

  // Loop over the L1's and make seeds for all of them:
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator it;
  for ( it = trks->begin(); it != trks->end(); it++ ) {
    const TTTrack< Ref_Phase2TrackerDigi_ >& l1 = (*it);

    std::unique_ptr<std::vector<TrajectorySeed> > out(new std::vector<TrajectorySeed>());
    std::cout << "SeedGeneratorFromTTracks::produce: L1 track pT, eta, phi --> " << l1.momentum().mag() << " , " << l1.momentum().eta() << " , "
	      << l1.momentum().phi() << std::endl;

    FreeTrajectoryState fts = trajectoryStateTransform::initialFreeStateTTrack(l1, magfieldH.product(), false);
    dummyPlane->move(fts.position() - dummyPlane->position());
    TrajectoryStateOnSurface tsosAtIP = TrajectoryStateOnSurface(fts, *dummyPlane);
    std::cout << "SeedGeneratorFromTTracks::produce: Created TSOSatIP: " << tsosAtIP << std::endl;

    unsigned int numSeedsMade = 0;
    //BARREL
    if (std::abs(l1.momentum().eta()) < theMaxEtaForTOB) {
      for (auto it = tob.rbegin(); it != tob.rend(); ++it) {  //This goes from outermost to innermost layer
	std::cout << "SeedGeneratorFromTTracks::produce: looping in TOB layer " << std::endl;
        findSeedsOnLayer(**it,
                         tsosAtIP,
                         *(propagatorAlong.get()),
                         l1,
                         estimatorH,
			 numSeedsMade,			
                         out);
      }
    }
    // //ENDCAP+
    // if (l1.momentum().eta() > theMinEtaForTEC) {
    //   for (auto it = tecPositive.rbegin(); it != tecPositive.rend(); ++it) {
    // 	std::cout << "SeedGeneratorFromTTracks::produce: looping in TEC+ layer " << std::endl;
    //     findSeedsOnLayer(**it,
    //                      tsosAtIP,
    //                      *(propagatorAlong.get()),
    //                      l1,
    //                      estimatorH,
    //                      out);
    //   }
    // }
    // //ENDCAP-
    // if (l1.momentum().eta() < -theMinEtaForTEC) {
    //   for (auto it = tecNegative.rbegin(); it != tecNegative.rend(); ++it) {
    // 	std::cout << "SeedGeneratorFromTTracks::produce: looping in TEC- layer "  << std::endl;
    //     findSeedsOnLayer(**it,
    //                      tsosAtIP,
    //                      *(propagatorAlong.get()),
    //                      l1,
    //                      estimatorH,
    //                      out);
    //   }
    // }
    for (std::vector<TrajectorySeed>::iterator it = out->begin(); it != out->end(); ++it) {
      result->push_back(*it);
    }
  } // end loop over L1Tracks
  
  std::cout << "SeedGeneratorFromTTracks::produce: number of seeds made: " << result->size() << std::endl;

  ev.put(std::move(result));
  std::cout << "SeedGeneratorFromTTracks::end " << std::endl;
}
