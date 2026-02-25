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
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CaloClusters.h"

#include <DataFormatsPHOS/BadChannelsMap.h>

#include "CCDB/BasicCCDBManager.h"

#include <cstdio>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SingleTracks {

  OutputObj<TH1F> ptTrack{TH1F("ptTrack", "pt of tracks", 200, 0., 20.)};
  // define global variables

  // loop over each single track
  void process(aod::TrackIU const& track)
  {
    ptTrack->Fill(track.pt());
  }
};


struct SingleTracks2 {

  OutputObj<TH1F> ptTrack{TH1F("ptTrack", "pt of tracks", 200, 0., 20.)};
  OutputObj<TH1F> etaTrack{TH1F("etaTrack", "eta of tracks", 500, -5.0f, 5.0f)};

  // loop over each single track
  void process(aod::TrackIU const& track)
  {
    ptTrack->Fill(track.pt());
    etaTrack->Fill(track.eta());
  }
};

struct TracksPerCollision {

//  Filter<aod::collision> vertexzFilter = nabs(collision::posZ) < 10.0f;

  OutputObj<TH1F> ntracks{TH1F("ntracks", "number of tracks", 2000, 0., 20000.)};
  OutputObj<TH1F> collz{TH1F("collz", "z coordinate of the vertex", 400, -20., 20.)};

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, 
	       aod::TracksIU const& tracks)
  {
    if (!collision.sel8()) return; 
    
    ntracks->Fill(tracks.size());
    collz->Fill(collision.posZ());
  }
};

struct PhotonsTest {

  // configure filters
  Configurable<int> ncellMin{"ncellMin", 3, "minimum number of cells in a cluster"};
  Configurable<float> enMin{"enMin", 0.3f, "minimum cluster energy in GeV"};
  Configurable<float> tofCut{"tofCut", 12.5*pow(10.f,-9), "cluster time in ns"};

  Filter filter1 = o2::aod::calocluster::e >= enMin;
  Filter filter2 = o2::aod::calocluster::ncell >= as<uint8_t>(ncellMin);
  Filter filter3 = nabs(o2::aod::calocluster::time) <= tofCut;

  using filteredCaloClusters = soa::Filtered<o2::aod::CaloClusters>;
  
  OutputObj<TH1F> ptCluster{TH1F("ptCluster", "pt of clusters", 200, 0.f, 20.f)};
  OutputObj<TH1F> timeCluster{TH1F("timeCluster", "time of clusters", 400, -20.f, 20.f)};
  OutputObj<TH1I> mod{TH1I("mod", "module number", 5, 0.5f, 5.5f)};
  OutputObj<TH1I> ncell{TH1I("ncell", "number of cells in a cluster", 200, 0.f, 200.f)};
  OutputObj<TH1I> nclust{TH1I("nclust", "number of clusters per event", 100, 0.f, 100.f)};
  OutputObj<TH2F> hmass{TH2F("hmass", "invariant mass;m_{inv};E_{1}+E_{2}", 750, 0.f, 0.75f, 200, 0.f, 20.f)};
  OutputObj<TH2F> hmassMix{TH2F("hmassMix", "mixed invariant mass;m_{inv};E_{1}+E_{2}", 750, 0.f, 0.75f, 200, 0.f, 20.f)};
  OutputObj<TH1F> hcent{TH1F("hcent", "collision centrality", 100, 0.f, 100.f)};
  OutputObj<TH1F> hdistTrack{(TH1F("hdistTrack", "distance to the closes track", 500, 0.f, 50.f))};

  HistogramRegistry moduleHists{"Module histograms"};


  void init(o2::framework::InitContext const&)
  {
    for (int i = 0; i < int(5); i ++) {
       moduleHists.add(Form("mHistXZ_mod%d", i+1), "cluster occupancy", 
		       HistType::kTH2F, 
		       {{64, -72., 72., "x", "x (cm)"}, 
		       {56, -63., 63., "z", "z (cm)"}});
    }
  }
  
  std::deque<aod::CaloClusters>  eventsCache = {};
  
  // loop over clusters
//  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>::iterator const& collision,
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>::iterator const& collision,
                filteredCaloClusters const& clusters)
  {
    if (!collision.sel8()) return;
    
//    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
//    LOG(info) << "timestamp: " << bc.timestamp();

    hcent->Fill(collision.centFT0M());

    if (collision.centFT0M() > 90.f) return;
    if (abs(collision.posZ() > 10.f)) return;

    nclust->Fill(clusters.size());

    for (const auto& clu : clusters) {
      float pt = std::hypot(clu.px(), clu.py());
      ptCluster->Fill(pt);
      timeCluster->Fill(clu.time()/1.e-9);
      mod->Fill(clu.mod());
      ncell->Fill(clu.ncell());
      hdistTrack->Fill(clu.trackdist());
      int imod = clu.mod();
      switch (imod) {
	 case 1:      
           moduleHists.fill(HIST("mHistXZ_mod1"), clu.x(), clu.z());
	 case 2:      
           moduleHists.fill(HIST("mHistXZ_mod2"), clu.x(), clu.z());
	 case 3:      
           moduleHists.fill(HIST("mHistXZ_mod3"), clu.x(), clu.z());
	 case 4:      
           moduleHists.fill(HIST("mHistXZ_mod4"), clu.x(), clu.z());
	 default:
	   break;
      }
    }
   
    for (auto& [clu1, clu2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(clusters, clusters))) {
       float p = std::hypot(clu1.px()+clu2.px(), clu1.py()+clu2.py(), clu1.pz()+clu2.pz());
       float e = clu1.e()+clu2.e();
       float m = std::sqrt(std::pow(e, 2)-std::pow(p, 2));
       if (abs((clu1.e()-clu2.e())/(clu1.e()+clu2.e())) > 0.1) continue; 
       hmass->Fill(m, clu1.e() + clu2.e());
    }

    for (const auto& event : eventsCache) {
      for (auto& [clu1, clu2] : combinations(soa::CombinationsFullIndexPolicy(clusters, event))) {
         float p = std::hypot(clu1.px()+clu2.px(), clu1.py()+clu2.py(), clu1.pz()+clu2.pz());
         float e = clu1.e()+clu2.e();
         float m = std::sqrt(std::pow(e, 2)-std::pow(p, 2));
         if (abs((clu1.e()-clu2.e())/(clu1.e()+clu2.e())) > 0.1) continue; 
         hmassMix->Fill(m, clu1.e() + clu2.e());
      }
    }

    eventsCache.emplace_back(clusters);

    while (eventsCache.size() > 25) {
       eventsCache.pop_front();
    }

  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SingleTracks>(cfgc),
    adaptAnalysisTask<SingleTracks2>(cfgc),
    adaptAnalysisTask<TracksPerCollision>(cfgc),
    adaptAnalysisTask<PhotonsTest>(cfgc),
  };
}
