// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod {
namespace helperSpace
{
// Global bool
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, bool); //! is this ambiguos?
}
DECLARE_SOA_TABLE(HelperTable, "AOD", "HELPERTABLE", helperSpace::IsAmbiguous);
}

using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

struct myExampleTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsEta{"nBinsEta", 300, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nMax{"nMax", 200, "Maximum number of tracks per collision"};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{nBinsEta, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0., 10., "p_{T}"};
    const AxisSpec axisCounter{1, 0, 1, ""};
    const AxisSpec axisMult{1000, 0, 1.*nMax, "n_{Tracks} per event"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("trackMultiplicity", "trackMultiplicity", kTH1F, {axisMult});
    histos.add("etaHistogramAmbi", "etaHistogramAmbi", kTH1F, {axisEta});
    histos.add("ptHistogramAmbi", "ptHistogramAmbi", kTH1F, {axisPt});
  }

  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  SliceCache cache;

  Partition<myFilteredTracks> leftTracks = aod::track::eta < 0.0f;
  Partition<myFilteredTracks> rightTracks = aod::track::eta >= 0.0f;
 
  //partitions are not grouped by default!
  //auto leftTracksThisCollision = leftTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
  //auto rightTracksThisCollision = rightTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
  // For manual sliceBy
  Preslice<aod::Tracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<aod::AmbiguousTracks> ambiTrackPreslice = o2::aod::ambiguous::trackId;

  void process(aod::Collisions const& collisions, 
	       myFilteredTracks const& tracks,
	       o2::aod::AmbiguousTracks const& ambiTracks)
  {
    for (auto& collision : collisions) {
       histos.fill(HIST("eventCounter"), 0.5);

       int nTracks = 0;
       //group tracks manually
       const uint64_t collIdx = collision.globalIndex();
       auto trackTable_thisCollision = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);
       histos.fill(HIST("eventCounter"), 0.5);
       for (auto& track : trackTable_thisCollision) {
         if (track.tpcNClsCrossedRows() < 70) continue;
         histos.fill(HIST("etaHistogram"), track.eta());
         histos.fill(HIST("ptHistogram"),  track.pt());
	 nTracks ++;
	 const uint64_t trackIdx = track.globalIndex();
	 auto ambiTable = ambiTracks.sliceBy(ambiTrackPreslice, trackIdx);
	 if (ambiTable.size() > 0) {
	    histos.fill(HIST("etaHistogramAmbi"), track.eta());
	    histos.fill(HIST("ptHistogramAmbi"), track.pt());
	 }
       }
       histos.fill(HIST("trackMultiplicity"), nTracks);
    } 
  }
};

struct myExampleTask2 {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f;
  Filter trackIsAmbi = aod::helperSpace::isAmbiguous == true;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    // create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogramAmbi", "etaHistogramAmbi", kTH1F, {axisEta});
    histos.add("ptHistogramAmbi", "ptHistogramAmbi", kTH1F, {axisPt});
  }

  void process(myFilteredTracks const& tracks, o2::aod::AmbiguousTracks const& ambiTracks)
  {
    for (const auto& track : tracks) {
      if( track.tpcNClsCrossedRows() < 70 ) continue; //badly tracked
      histos.fill(HIST("etaHistogramAmbi"), track.eta());
      histos.fill(HIST("ptHistogramAmbi"), track.pt());
    }
  }
};

struct myExampleProducer {
  Produces<aod::HelperTable> myHelperTable;
  void init(InitContext const&) {}
  void process(o2::aod::Tracks const& tracks, o2::aod::AmbiguousTracks const& ambiTracks)
  {
     std::vector<bool> ambiguousOrNot;
     ambiguousOrNot.reserve(tracks.size());
     for (int ii = 0; ii < tracks.size(); ii++)
        ambiguousOrNot[ii] = false;
     for (auto const& ambiTrack : ambiTracks) {
        ambiguousOrNot[ambiTrack.trackId()] = true;
     }
     for (int ii = 0; ii < tracks.size(); ii++) {
	myHelperTable(ambiguousOrNot[ii]);
     }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc),
    adaptAnalysisTask<myExampleTask2>(cfgc),
    adaptAnalysisTask<myExampleProducer>(cfgc),
    };
}
