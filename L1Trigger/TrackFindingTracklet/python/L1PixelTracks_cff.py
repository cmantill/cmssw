import FWCore.ParameterSet.Config as cms

from RecoVertex.BeamSpotProducer.BeamSpot_cfi import *

process.load("L1Trigger.TrackFindingTracklet.Tracklet_cfi")


TTPixelTracksFromTracklet = cms.EDProducer("L1Pixel1TrackProducer",
                                           L1Tk_nPar = cms.int32(4),
                                           L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
                                           PixelTracksInputTag = cms.InputTag("pixelTracks"),
    )

TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")
TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag(cms.InputTag("TTTracksFromTracklet", "Level1TTTracks") )

L1TrackletTracks = cms.Sequence(offlineBeamSpot*TTTracksFromTracklet)
L1TrackletTracksWithAssociators = cms.Sequence(offlineBeamSpot*TTTracksFromTracklet*TrackTriggerAssociatorTracks)

