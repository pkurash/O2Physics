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

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include "Framework/HistogramRegistry.h"
#include <TParameter.h>

#include "Framework/StaticFor.h"

#include <cstddef>

using namespace o2;
using namespace o2::framework;

struct SingleTracks {

  // define global variables
  size_t count = 0;

  // loop over each single track
  void process(aod::TrackIU const& track)
  {
    // count the tracks contained in the input file
    LOGF(info, "Track %d: Momentum: %f", count, track.p());
    count++;
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

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SingleTracks>(cfgc),
    adaptAnalysisTask<AllTracks>(cfgc),
    adaptAnalysisTask<EtaPhiPtHistograms>(cfgc),
  };
}
