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
//
/// \brief
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <TParameter.h>

#include "Framework/ASoAHelpers.h"
#include "Framework/StaticFor.h"
//#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SingleTracks {

  // define global variables
  size_t count = 0;

  // loop over each single track
  void process(aod::Track const& track)
  {
    // count the tracks contained in the input file
   // LOGF(info, "Track %d: Momentum: %f", count, track.p());
    count++;
  }
};

struct FilterSingleTracks {

  HistogramRegistry registry{
    "registry",
    {
      {"pt", "pt", {HistType::kTH1F, {{400, 0., 40.}}}, true},     //
      {"eta", "#eta", {HistType::kTH1F, {{102, -2.01, 2.01}}}},     //
    }                                                               //
  };

  Filter etaFilter = aod::track::eta > 0.f;

  void process(soa::Filtered<aod::TracksIU> const& tracks)
  {
     for (auto& track : tracks) {
       registry.get<TH1>(HIST("pt"))->Fill(track.pt());
       registry.get<TH1>(HIST("eta"))->Fill(track.eta());
     }
  }
};

struct PartitionSingleTracks {

  HistogramRegistry registry{
    "registry",
    {
      {"ptleft", "pt left", {HistType::kTH1F, {{400, 0., 40.}}}, true},     //
      {"ptright", "pt right", {HistType::kTH1F, {{400, 0., 40.}}}, true},     //
    }                                                               //
  };

  Partition<o2::aod::TracksIU> leftTracks = aod::track::eta < 0.f;
  Partition<o2::aod::TracksIU> rightTracks = aod::track::eta >= 0.f;

  void process(aod::TracksIU const& tracks)
  {
     for (auto& track : leftTracks) {
       registry.get<TH1>(HIST("ptleft"))->Fill(track.pt());
     }
     for (auto& track : rightTracks) {
       registry.get<TH1>(HIST("ptright"))->Fill(track.pt());
     }
  }
};

struct AllTracks {

  // define global variables
  size_t numberDataFrames = 0;
  size_t count = 0;
  size_t totalCount = 0;

  // loop over data frames
  void process(aod::TracksIU const& tracks)
  {
    numberDataFrames++;

    // count the tracks contained in each data frame
    count = 0;
    for (auto& track : tracks) {
      LOGF(debug, "Track with momentum %f", track.pt());
      count++;
    }
    totalCount += count;

    LOGF(info, "DataFrame %d: Number of tracks: %d Accumulated number of tracks: %d", numberDataFrames, count, totalCount);
  }
};

struct EtaPhiPtHistograms {
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{
    "registry",
    {
      {"pt", "pt", {HistType::kTH1F, {{400, 0., 40.}}}, true},     //
      {"eta", "#eta", {HistType::kTH1F, {{102, -2.01, 2.01}}}},     //
      {"eta2", "#eta", {HistType::kTH1F, {{102, -2.01, 2.01}}}},     //
      {"phi", "#varphi", {HistType::kTH1F, {{100, 0., 2. * M_PI}}}, true}, //
      {"ptToPt", "#ptToPt", {HistType::kTH2F, {{100, -0.01, 10.01}, {100, -0.01, 10.01}}}} //
    }                                                               //
  };

  void process(aod::TracksIU const& tracks)
  {
    registry.fill<aod::track::Eta>(HIST("eta2"), tracks, aod::track::pt > 1.0f);
    registry.fill<aod::track::Pt, aod::track::Pt>(HIST("ptToPt"), tracks, aod::track::pt < 5.0f);
    for (auto& track : tracks) {
      registry.get<TH1>(HIST("pt"))->Fill(track.pt());
      registry.get<TH1>(HIST("eta"))->Fill(track.eta());
      registry.get<TH1>(HIST("phi"))->Fill(track.phi());
    }
  }
};


struct TracksPerCollision {
  HistogramRegistry registry{
    "registry",
    {
      {"ntracks", "ntracks", {HistType::kTH1F, {{2000, 0., 20000.}}}},     //
      {"collZ", "collisionZ", {HistType::kTH1F, {{400, -20., 20.}}}},     //
    }                                                               //
  };

//  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::TracksIU const& tracks)
  void process(aod::Collision const& collision, aod::TracksIU const& tracks)
  {
    // `tracks` contains tracks belonging to `collision`
   // LOGF(info, "Collision index : %d", collision.index());
   // LOGF(info, "Number of tracks: %d", tracks.size());

//    if (!collision.sel8()) {
//      return;
//    }

    registry.get<TH1>(HIST("ntracks"))->Fill(tracks.size());
    registry.get<TH1>(HIST("collZ"))->Fill(collision.posZ());

    // process the tracks of a given collision
//    for (auto& track : tracks) {
//      LOGF(info, "  track pT = %f GeV/c", track.pt());
//    }
  }
};

struct TracksPerDataframe {

  HistogramRegistry registry{
    "registry",
    {
      {"ntracks", "ntracks", {HistType::kTH1F, {{2000, 0., 200000.}}}},     //
      {"ncoll", "ncoll", {HistType::kTH1F, {{200, 0., 200.}}}},     //
    }                                                               //
  };

  void process(aod::Collisions const& collisions, aod::TracksIU const& tracks)
  {
    // `tracks` contains all tracks of a data frame
  //  LOGF(info, "Number of collisions: %d", collisions.size());
  //  LOGF(info, "Number of tracks    : %d", tracks.size());

    registry.get<TH1>(HIST("ncoll"))->Fill(collisions.size());
    registry.get<TH1>(HIST("ntracks"))->Fill(tracks.size());
  }
};

struct GroupByCollision {

  HistogramRegistry registry{
    "registry",
    {
      {"ntracks", "ntracks", {HistType::kTH1F, {{2000, 0., 20000.}}}},     //
      {"nv0s", "v0s", {HistType::kTH1F, {{2000, 0., 20000.}}}},     //
    }
  };

  void process(aod::Collision const& collision, aod::TracksIU const& tracks, aod::V0s const& v0s)
  {

    // `tracks` contains tracks belonging to `collision`
    // `v0s`    contains V0s    belonging to `collision`
    LOGF(info, "Collision index : %d", collision.index());
    LOGF(info, "Number of tracks: %d", tracks.size());
    LOGF(info, "Number of v0s   : %d", v0s.size());
    registry.get<TH1>(HIST("ntracks"))->Fill(tracks.size());
    registry.get<TH1>(HIST("nv0s"))->Fill(v0s.size());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SingleTracks>(cfgc),
    adaptAnalysisTask<FilterSingleTracks>(cfgc),
    adaptAnalysisTask<PartitionSingleTracks>(cfgc),
    adaptAnalysisTask<AllTracks>(cfgc),
    adaptAnalysisTask<EtaPhiPtHistograms>(cfgc),
    adaptAnalysisTask<TracksPerCollision>(cfgc),
    adaptAnalysisTask<TracksPerDataframe>(cfgc),
    adaptAnalysisTask<GroupByCollision>(cfgc),
  };
}
