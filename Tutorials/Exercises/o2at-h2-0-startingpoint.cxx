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
/// \brief This task contains the individual steps that are to be taken
///        in the second part of the tutorial. These are 5 steps, and at the end,
///        the participant is expected to have a two-particle correlation spectrum.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


//STEP 0
//This is an example of a conveient declaration of "using"
//WARNING: THIS IS NEW
using MyCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

Double_t ComputeDeltaPhi( Double_t phi1, Double_t phi2);

//Starting point
struct startingpoint {

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"etaHistogram", "etaHistogram", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogram", "ptHistogram", {HistType::kTH1F, {{nBins, 0., 10.0}}}}
    }
  };

  void process(aod::Collision const& collision, MyCompleteTracks const& tracks) //<- this is the main change
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
    //This will take place once per event!
    for (auto& track : tracks) {
      //"Good old" imperative cuts
      if( fabs(track.eta()) > 0.5 ) continue;
      if( track.tpcNClsCrossedRows() < 70 ) continue;
      if( fabs(track.dcaXY()) > .2 ) continue;
      registry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      registry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
    }
  }
};


struct filterexample {

  Filter etaFilter = nabs(aod::track::eta) < 0.5f;
  Filter dcaFilter = nabs(aod::track::dcaXY) < 0.2f;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"etaHistogram", "etaHistogram", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogram", "ptHistogram", {HistType::kTH1F, {{nBins, 0., 10.0}}}}
    }
  };

  void process(aod::Collision const& collision, 
	       soa::Filtered<MyCompleteTracks> const& tracks) //<- this is the main change
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
    //This will take place once per event!
    for (auto& track : tracks) {
      //"Good old" imperative cuts
      if( track.tpcNClsCrossedRows() < 70 ) continue;
      registry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      registry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
    }
  }
};

struct partitionexample {

  SliceCache cache;

  //Filter etaFilter = nabs(aod::track::eta) < 0.5f;
  //Filter dcaFilter = nabs(aod::track::dcaXY) < 0.2f;

  Partition<MyCompleteTracks> leftTracks  = aod::track::eta < 0.0f;
  Partition<MyCompleteTracks> rightTracks = aod::track::eta > 0.0f;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"etaHistogramleft", "etaHistogramleft", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramleft", "ptHistogramleft", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"etaHistogramright", "etaHistogramright", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramright", "ptHistogramright", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
    }
  };

  void process(aod::Collision const& collision, 
	       /*soa::Filtered<*/MyCompleteTracks/*>*/ const& tracks) //<- this is the main change
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());

    auto leftTracksGrouped = leftTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto rightTracksGrouped = rightTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (auto& track : leftTracksGrouped) {
       registry.get<TH1>(HIST("etaHistogramleft"))->Fill(track.eta());
       registry.get<TH1>(HIST("ptHistogramleft"))->Fill(track.pt());
    }

    for (auto& track : rightTracksGrouped) {
       registry.get<TH1>(HIST("etaHistogramright"))->Fill(track.eta());
       registry.get<TH1>(HIST("ptHistogramright"))->Fill(track.pt());
    }

  }
};

struct partitionfilterexample {

  SliceCache cache;

  Filter etaFilter = nabs(aod::track::eta) < 0.5f;
  Filter dcaFilter = nabs(aod::track::dcaXY) < 0.2f;
  using MyFilteredTracks =  soa::Filtered<MyCompleteTracks>;

  Partition<MyFilteredTracks> leftTracks  = aod::track::eta < 0.0f;
  Partition<MyFilteredTracks> rightTracks = aod::track::eta > 0.0f;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"etaHistogramleft", "etaHistogramleft", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramleft", "ptHistogramleft", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"etaHistogramright", "etaHistogramright", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramright", "ptHistogramright", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
    }
  };

  void process(aod::Collision const& collision, 
	       MyFilteredTracks const& tracks) //<- this is the main change
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());

    auto leftTracksGrouped = leftTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto rightTracksGrouped = rightTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (auto& track : leftTracksGrouped) {
       if (track.tpcNClsCrossedRows() < 70) continue;
       registry.get<TH1>(HIST("etaHistogramleft"))->Fill(track.eta());
       registry.get<TH1>(HIST("ptHistogramleft"))->Fill(track.pt());
    }

    for (auto& track : rightTracksGrouped) {
       if (track.tpcNClsCrossedRows() < 70) continue;
       registry.get<TH1>(HIST("etaHistogramright"))->Fill(track.eta());
       registry.get<TH1>(HIST("ptHistogramright"))->Fill(track.pt());
    }

  }
};

struct computecorrelation {

  SliceCache cache;

  Filter dcaFilter = nabs(aod::track::dcaXY) < 0.2f;
  Filter trackFilter = nabs(aod::track::eta) < 0.8f && aod::track::pt > 2.0f;
  using MyFilteredTracks =  soa::Filtered<MyCompleteTracks>;

  Partition<MyFilteredTracks> triggerTracks  = aod::track::pt > 4.0f;
  Partition<MyFilteredTracks> assocTracks = aod::track::pt < 4.0f;
  
  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{3*nBins, -15., 15.}}}},
      {"etaHistogramtrigger", "etaHistogramtrigger", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramtrigger", "ptHistogramtrigger", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"etaHistogramassoc", "etaHistogramassoc", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramassoc", "ptHistogramassoc", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"correlationFunction", "correlationFunction", {HistType::kTH1F, {{nBins, 0, 2*TMath::Pi()}}}},
      {"correlationFunction2D", "correlationFunction2D", {HistType::kTH2F, {{nBins, 0, 2*TMath::Pi()}, {16, -0.8, 0.8}}}},
    }
  };

  void process(aod::Collision const& collision, 
	       MyFilteredTracks const& tracks) //<- this is the main change
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());

    auto triggerTracksGrouped = triggerTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto assocTracksGrouped =  assocTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);


    for (auto& [trackTrigger, trackAssoc] : combinations(soa::CombinationsFullIndexPolicy(triggerTracksGrouped, assocTracksGrouped))){
       if (trackTrigger.tpcNClsCrossedRows() < 70) continue;
       if (trackAssoc.tpcNClsCrossedRows() < 70) continue;
       registry.get<TH1>(HIST("etaHistogramtrigger"))->Fill(trackTrigger.eta());
       registry.get<TH1>(HIST("ptHistogramtrigger"))->Fill(trackTrigger.pt());
       registry.get<TH1>(HIST("etaHistogramassoc"))->Fill(trackAssoc.eta());
       registry.get<TH1>(HIST("ptHistogramassoc"))->Fill(trackAssoc.pt());

       double deltaPhi = ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi());
       double deltaEta = trackTrigger.eta() - trackAssoc.eta();
       registry.get<TH1>(HIST("correlationFunction"))->Fill(deltaPhi);
       registry.get<TH2>(HIST("correlationFunction2D"))->Fill(deltaPhi, deltaEta);
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<startingpoint>(cfgc),
    adaptAnalysisTask<filterexample>(cfgc),
    adaptAnalysisTask<partitionexample>(cfgc),
    adaptAnalysisTask<partitionfilterexample>(cfgc),
    adaptAnalysisTask<computecorrelation>(cfgc),
  };
}
Double_t ComputeDeltaPhi( Double_t phi1, Double_t phi2) {
   //To be completely sure, use inner products
   Double_t x1, y1, x2, y2;
   x1 = TMath::Cos( phi1 );
   y1 = TMath::Sin( phi1 );
   x2 = TMath::Cos( phi2 );
   y2 = TMath::Sin( phi2 );
   Double_t lInnerProd = x1*x2 + y1*y2;
   Double_t lVectorProd = x1*y2 - x2*y1;
   Double_t lReturnVal = 0;
   if( lVectorProd > 1e-8 ) {
      lReturnVal = TMath::ACos(lInnerProd);
   }
   if( lVectorProd < -1e-8 ) {
      lReturnVal = -TMath::ACos(lInnerProd);
   }
   if( lReturnVal < -TMath::Pi()/2. ) {
      lReturnVal += 2.*TMath::Pi();
   }
   return lReturnVal;
}
