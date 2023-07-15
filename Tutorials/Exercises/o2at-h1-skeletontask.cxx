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
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
 #include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;

//STEP 0
//This is an empty analysis skeleton: the starting point! 
struct taskskeleton {
  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "etaHistogram", 200, -1., +1)};

  void process(aod::TracksIU const& tracks)
  {
    for (auto& track : tracks) {
      etaHistogram->Fill(track.eta());
    }
  }
};

struct momentumexample {
  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "etaHistogram", 200, -1., +1)};
  OutputObj<TH1F> ptHistogram{TH1F("ptHistogram", "ptHistogram", 100, 0., 10.)};

  void process(aod::TracksIU const& tracks)
  {
    for (auto& track : tracks) {
      etaHistogram->Fill(track.eta());
      ptHistogram->Fill(track.pt());
    }
  }
};

struct configexample{
  Configurable<int> nBinsEta{"nBinsEta",200, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt histo"};

  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> etaHistogram{"etaHistogram"};
  OutputObj<TH1F> ptHistogram{"ptHistogram"};

  void init(InitContext const&)
  {
    etaHistogram.setObject(new TH1F("etaHistogram", "etaHistogram", nBinsEta, -1, +1));
    ptHistogram.setObject(new TH1F("ptHistogram", "ptHistogram", nBinsPt, 0, 10.));
  }

  void process(aod::TracksIU const& tracks)
  {
    for (auto& track : tracks) {
      etaHistogram->Fill(track.eta());
      ptHistogram->Fill(track.pt());
    }
  }
};

struct histogramregistry{
  Configurable<int> nBinsEta{"nBinsEta",200, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt histo"};

  HistogramRegistry histogramregistry{
    "registry",
    {
      {"etaHistogram", "etaHistogram", {HistType::kTH1F, {{nBinsEta, -1., 1.}}}},
      {"ptHistogram", "ptHistogram", {HistType::kTH1F, {{nBinsPt, 0., 10.}}}},
    }
  };
  void process(aod::TracksIU const& tracks)
  {
    for (auto& track : tracks) {
      histogramregistry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      histogramregistry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
    }
  }
};

struct newsubscriptions{
  Configurable<int> nBinsEta{"nBinsEta",200, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt histo"};

  HistogramRegistry histogramregistry{
    "registry",
    {
      {"etaHistogram", "#eta", {HistType::kTH1F, {{nBinsEta, -1., 1.}}}},
      {"ptHistogram", "p_{T}", {HistType::kTH1F, {{nBinsPt, 0., 10.}}}},
      {"hZVertex", "vertex z coordinate", {HistType::kTH1F, {{400, -20, 20}}}},
    }
  };
  void process(aod::Collision const& collision, 
	       soa::Join<aod::TracksIU, aod::TracksExtra> const& tracks)
  {

    histogramregistry.get<TH1>(HIST("hZVertex"))->Fill(collision.posZ());

    for (auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < 70) continue;
      histogramregistry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      histogramregistry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
    }
  }
};

struct adddca{
  Configurable<int> nBinsEta{"nBinsEta",200, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt histo"};

  HistogramRegistry histogramregistry{
    "registry",
    {
      {"etaHistogram", "#eta", {HistType::kTH1F, {{nBinsEta, -1., 1.}}}},
      {"ptHistogram", "p_{T}", {HistType::kTH1F, {{nBinsPt, 0., 10.}}}},
      {"hZVertex", "vertex z coordinate", {HistType::kTH1F, {{400, -20, 20}}}},
    }
  };
  void process(aod::Collision const& collision, 
	       soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA> const& tracks)
  {

    histogramregistry.get<TH1>(HIST("hZVertex"))->Fill(collision.posZ());

    for (auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < 70) continue;
      if (fabs(track.dcaXY()) > 2.) continue;
      histogramregistry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      histogramregistry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
    }
  }
};

struct computeresolution{
  Configurable<int> nBinsEta{"nBinsEta",200, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt histo"};

  HistogramRegistry histogramregistry{
    "registry",
    {
      {"etaHistogram", "#eta", {HistType::kTH1F, {{nBinsEta, -1., 1.}}}},
      {"ptHistogram", "p_{T}", {HistType::kTH1F, {{nBinsPt, 0., 10.}}}},
      {"hZVertex", "vertex z coordinate", {HistType::kTH1F, {{400, -20, 20}}}},
      {"resoHistogram", "p_{T} resolution", {HistType::kTH2F, {{nBinsPt, 1., 10.}, {200, -1., 1.}} }},
    }
  };
  void process(aod::Collision const& collision, 
	       soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& tracks,
	       aod::McParticles const&)
  {

    histogramregistry.get<TH1>(HIST("hZVertex"))->Fill(collision.posZ());

    for (auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < 70) continue;
      if (fabs(track.dcaXY()) > 2.) continue;
      histogramregistry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      histogramregistry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());

      auto mcParticle = track.mcParticle_as<aod::McParticles>();
      float delta = track.pt() - mcParticle.pt();
      histogramregistry.get<TH2>(HIST("resoHistogram"))->Fill(track.pt(), delta);
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskskeleton>(cfgc),
    adaptAnalysisTask<momentumexample>(cfgc),
    adaptAnalysisTask<configexample>(cfgc),
    adaptAnalysisTask<histogramregistry>(cfgc),
    adaptAnalysisTask<newsubscriptions>(cfgc),
    adaptAnalysisTask<adddca>(cfgc),
    adaptAnalysisTask<computeresolution>(cfgc),
  };
}
