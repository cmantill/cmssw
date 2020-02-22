#include "L1PixelTracksProducer.h"

#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "RecoPixelVertexing/PixelTrackFitting/interface/CircleFromThreePoints.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/L1TrackingRecHit.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

using namespace std;
using namespace edm;

L1PixelTracksProducer::L1PixelTracksProducer(const edm::ParameterSet& iConfig)
{
  L1Tk_nPar          = iConfig.getParameter< int >("L1Tk_nPar");

  ttTrackToken_      = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter<edm::InputTag>("L1TrackInputTag"));
  ttPixelTrackToken_ = consumes< std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("PixelTracksInputTag"));

  produces< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >( "Level1TTTracks" ).setBranchAlias("Level1TTTracks");
  produces<L1TrackingDetSetVector>("Level1TTRecHits").setBranchAlias("Level1TTRecHits");;
}

L1PixelTracksProducer::~L1PixelTracksProducer() {}

void L1PixelTracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("L1TrackInputTag")->setComment("L1 tracks");
  desc.add<edm::InputTag>("PixelTracksInputTag")->setComment("pixel tracks");
  descriptions.add("L1PixelTracksProducer",desc);
}

void L1PixelTracksProducer::produce(edm::Event& iEvent, 
				    const edm::EventSetup& iSetup) {
  // L1 tracks
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  // Pixel tracks
  edm::Handle< std::vector<reco::Track> > PixelTrackHandle;
  iEvent.getByToken(ttPixelTrackToken_, PixelTrackHandle);

  // Geometry
  //edm::ESHandle<TrackerTopology> tTopoHandle;
  //iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  edm::ESHandle<TrackerGeometry> tGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  //const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  // tmp arrays for matching
  std::vector < double> pixel_pts,pixel_etas,pixel_phis;
  std::vector < double > curvaturesPixel, curvatureL1;
  std::vector < double > xS_center_pixel, yS_center_pixel;
  std::vector < double > xS_center_L1, yS_center_L1;

  // output
  auto outputhits = std::make_unique<L1TrackingDetSetVector>();
  auto& theoutputhits = *outputhits;
  std::unique_ptr< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > outputtracks( new std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > );

  // loop over pixel tracks
  int this_pixel = 0;
  std::vector< reco::Track >::const_iterator iterPixel;
  for ( iterPixel = PixelTrackHandle->begin(); iterPixel != PixelTrackHandle->end(); iterPixel++)
    {
      edm::Ptr< reco::Track > pixel_ptr(PixelTrackHandle, this_pixel);
      this_pixel++;
      float tmp_trk_pt  = iterPixel->pt();

      if(pixel_ptr->recHitsSize()<3)
	continue;

      if(tmp_trk_pt <= 1.9)
	continue;

      std::vector < double> X,Y,Z;

      // loop over pixel tracks rechits
      int pixCounter = 0;
      for(trackingRecHit_iterator it = iterPixel->recHitsBegin(); it!=iterPixel->recHitsEnd(); ++it)
	{
	  if(pixCounter>=3) break;
	  const TrackingRecHit* hit = &(**it);
	  //const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(hit);

	  X.push_back(hit->globalPosition().x());
	  Y.push_back(hit->globalPosition().y());
	  Z.push_back(hit->globalPosition().z());
	  pixCounter++;
	} // end loop over pixel tracks rechits  

      CircleFromThreePoints circlePixel(GlobalPoint(X[0],Y[0],0.0),GlobalPoint(X[1],Y[1],0.0),GlobalPoint(X[2],Y[2],0.0));
      curvaturesPixel.push_back(circlePixel.curvature());
      xS_center_pixel.push_back(circlePixel.center().x());
      yS_center_pixel.push_back(circlePixel.center().y());

      pixel_pts.push_back(tmp_trk_pt);
      pixel_etas.push_back(iterPixel->eta());
      pixel_phis.push_back(iterPixel->phi());
    } // end loop over pixel tracks 

  // loop over L1 tracks  (with at least 3 stubs)
  int this_l1track = 0;
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
  for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = iterL1Track->getStubRefs();
    int tmp_trk_nstub  = (int) stubRefs.size();
    if(tmp_trk_nstub < 3)
      continue;

    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
    this_l1track++;

    // L1 trk properties
    float tmp_trk_pt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
    float tmp_trk_eta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
    float tmp_trk_phi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
    //float tmp_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm                                                                                                             
    
    //float tmp_trk_d0 = -999;
    //if (L1Tk_nPar == 5) {
    //  float tmp_trk_x0   = iterL1Track->getPOCA(L1Tk_nPar).x();
    //  float tmp_trk_y0   = iterL1Track->getPOCA(L1Tk_nPar).y();
    //  tmp_trk_d0 = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
    //}

    //float tmp_trk_chi2 = iterL1Track->getChi2(L1Tk_nPar);
    //float tmp_trk_bendchi2 = iterL1Track->getStubPtConsistency(L1Tk_nPar);

    std::vector<double> xS, yS, zS;
    // loop over stubs
    for (int is=0; is<3; is++) {
      DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();
      
      MeasurementPoint coords = stubRefs.at(is)->getClusterRef(0)->findAverageLocalCoordinatesCentered();
      const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
      Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );
      
      xS.push_back(posStub.x());
      yS.push_back(posStub.y());
      zS.push_back(posStub.z());
      
    }//end loop over stubs                                                                                                                                                    

    CircleFromThreePoints circleL1(GlobalPoint(xS[0],yS[0],0.0),GlobalPoint(xS[1],yS[1],0.0),GlobalPoint(xS[2],yS[2],0.0));
    curvatureL1.push_back(circleL1.curvature());
    xS_center_L1.push_back(circleL1.center().x());
    yS_center_L1.push_back(circleL1.center().y());
    double thisCurv = circleL1.curvature();

    // match with pixel tracks (pt,eta,phi,curvature,dist)
    int counter = 0;
    for(unsigned int i = 0; i<pixel_pts.size();i++)
      {
	bool test_pt = 2 * (pixel_pts[i] - tmp_trk_pt) / (pixel_pts[i] + tmp_trk_pt) < 0.05;
	if(!test_pt) continue;
	bool test_eta = 2 * (pixel_etas[i] - tmp_trk_eta) / (pixel_etas[i] + tmp_trk_eta) < 0.05;
	if(!test_eta) continue;
	bool test_phi = 2 * (pixel_phis[i] - tmp_trk_phi) / (pixel_phis[i] + tmp_trk_phi) < 0.05;
	if(!test_phi) continue;

	bool curvature_test = 2 * (curvaturesPixel[i] - thisCurv) / (curvaturesPixel[i] + thisCurv) < 0.05;
	double dist = (circleL1.center().x()*circleL1.center().x() + circleL1.center().y()*circleL1.center().y());
	dist = dist + (xS_center_pixel[i]*xS_center_pixel[i] + yS_center_pixel[i]*yS_center_pixel[i]);
	dist = sqrt(dist);
	if(test_pt && test_eta && test_phi && curvature_test && dist < 200.0) {
	  counter++;
	}
      }

    // for matched tracks?
    // for now RecHits rely on Stubs position - should we do a CPE?
    if(counter > 0){
      for (int is=0; is<tmp_trk_nstub; is++) {
	DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();

	MeasurementPoint coords = stubRefs.at(is)->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);

	LocalPoint lp = theGeomDet->topology().localPosition(coords);
	LocalError le;

	const auto genericDet = theTrackerGeom->idToDetUnit( (stubRefs.at(is)->getClusterRef(0))->getDetId() );
	edm::Ref<edmNew::DetSetVector< TTCluster < Ref_Phase2TrackerDigi_ > >, TTCluster < Ref_Phase2TrackerDigi_ > > cluster = stubRefs.at(is)->getClusterRef(0);
	L1TrackingRecHit hit(lp, le, *genericDet, cluster);

	std::cout << "L1_TRH: " << hit.localPosition().x() << "," << hit.localPosition().y() << " : "
		  << hit.localPositionError().xx() << "," << hit.localPositionError().yy() << std::endl;

	L1TrackingDetSetVector::FastFiller recHitsOnDet(theoutputhits,(stubRefs.at(is)->getClusterRef(0))->getDetId());
	recHitsOnDet.push_back(hit);
      } // end loop over stubs?
      outputtracks->push_back(*l1track_ptr);
    }

    std::cout << "pixel tracks and L1 tracks matched " << counter << std::endl;
  } // end loop over L1 tracks 

  iEvent.put( std::move(outputhits),"Level1TTRecHit"); // for now its saving all rechits?
  iEvent.put( std::move(outputtracks),"Level1TTTracks");
}
