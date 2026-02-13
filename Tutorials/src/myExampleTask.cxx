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
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct myExampleTask {
  // configurables
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 70, "Minimum number of crossed rows"};
  Configurable<float> maxDCAxy{"dcaXY", 0.2f, "Maximum distance of closest approach"};
  Configurable<float> maxPVz{"maxPVz", 100.0f, "max value of PVz (cm)"};

  Filter trackDCA = nabs(aod::track::dcaXY) < maxDCAxy;

  using myCompleteTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;
  using myFilteredTracksMC = soa::Filtered<myCompleteTracksMC>;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
	
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {},
  OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
     // define axes you want to use
     const AxisSpec axisEta{300, -1.5, +1.5, "#eta"};
     const AxisSpec axisPt{nBinsPt, 0., 10., "p_{T}"};
     const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};

     // create histograms
     histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
     histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt}, true);
     histos.add("hEventCounter", "event counter", kTH1F, {{1, 0.f, 1.f}});
     histos.add("ptResolution", "ptResolution", kTH2F, {axisPt, axisDeltaPt});
     histos.add("ptHistogramPion", "ptHistogramPion", kTH1F, {axisPt});
     histos.add("ptHistogramKaon", "ptHistogramKaon", kTH1F, {axisPt});
     histos.add("ptHistogramProton", "ptHistogramProton", kTH1F, {axisPt});
     histos.add("ptGeneratedPion", "ptGeneratedPion", kTH1F, {axisPt});
     histos.add("ptGeneratedKaon", "ptGeneratedKaon", kTH1F, {axisPt});
     histos.add("ptGeneratedProton", "ptGeneratedProton", kTH1F, {axisPt});
     histos.add("numberOfRecoCollisions", "numbers of collisions", kTH1F, {{11, -0.5f, 10.5f}});
  }

  template <typename TCollision, typename TTracks>
  void processStuff(TCollision const& collision, TTracks const& tracks) {
     histos.fill(HIST("hEventCounter"), 0.5f);
     for (const auto& track : tracks) {
         if ( track.tpcNClsCrossedRows() < minTPCCrossedRows ) continue; //badly tracked
         histos.fill(HIST("etaHistogram"), track.eta());
         histos.fill(HIST("ptHistogram"), track.pt());

         // this part only done if dealing with MC
         if constexpr (requires { track.has_mcParticle(); }) { // does the getter exist?
           if (track.has_mcParticle()) { // is the return 'true'? -> N.B. different question!
              auto mcParticle = track.mcParticle();
              histos.fill(HIST("ptResolution"), track.pt(), track.pt() - mcParticle.pt());
              if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.y())<0.5) {
                 if (abs(mcParticle.pdgCode()) == 221) histos.fill(HIST("ptHistogramPion"), mcParticle.pt());
                 if (abs(mcParticle.pdgCode()) == 321) histos.fill(HIST("ptHistogramKaon"), mcParticle.pt());
                 if (abs(mcParticle.pdgCode()) == 2212) histos.fill(HIST("ptHistogramProton"), mcParticle.pt());
	      }  
	   }  
	 }   
     }     
  }

 // void process(aod::TracksIU const& tracks)
//  void process(aod::Collision const& collision, myTracks const& tracks)
  void processReco(aod::Collision const& collision,
		   myFilteredTracks const& tracks)
  {
    processStuff(collision, tracks);
   // if (fabs(collision.posZ()) > maxPVz) {
   //   return;
   // }

   // histos.fill(HIST("hEventCounter"), 0.5);

   // for (auto& track : tracks) {
   //   if( track.tpcNClsCrossedRows() < minTPCCrossedRows ) continue; //badly tracked
   //   if( fabs(track.dcaXY()) > maxDCAxy) continue; //doesn’t point to primary vertex
   //    histos.fill(HIST("etaHistogram"), track.eta());
   //    histos.fill(HIST("ptHistogram"), track.pt());
   //    if (track.has_mcParticle()) {
   //      auto mcParticle = track.mcParticle();
   //      histos.fill(HIST("ptResolution"), track.pt(), track.pt() - mcParticle.pt());
   //      if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.y()) < 0.5) {
   //         if(abs(mcParticle.pdgCode())==211) histos.fill(HIST("ptHistogramPion"), mcParticle.pt());
   //         if(abs(mcParticle.pdgCode())==321) histos.fill(HIST("ptHistogramKaon"), mcParticle.pt());
   //         if(abs(mcParticle.pdgCode())==2212) histos.fill(HIST("ptHistogramProton"), mcParticle.pt());
   //      }
   //    }
   // }
  }

  PROCESS_SWITCH(myExampleTask, processReco, "process reconstructed information", true);

//  void processSim(aod::McParticles const& mcParticles)
  void processSim(aod::McCollision const& collision,
		   myFilteredTracksMC const& tracks, 
		  aod::McParticles const& mcParticles)
  {
    processStuff(collision, tracks);
//     histos.fill(HIST("numberOfRecoCollisions"), collisions.size());
//     for (const auto& mcParticle : mcParticles) {
//       if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.y()) < 0.5) { // watch out for context!!!
//         if (abs(mcParticle.pdgCode())==211) histos.fill(HIST("ptGeneratedPion"), mcParticle.pt());
//         if (abs(mcParticle.pdgCode())==321) histos.fill(HIST("ptGeneratedKaon"), mcParticle.pt());
//         if (abs(mcParticle.pdgCode())==2212) histos.fill(HIST("ptGeneratedProton"), mcParticle.pt());
//       }
//     }
  }

  PROCESS_SWITCH(myExampleTask, processSim, "process pure simulation information", true);
};

struct myExampleTask2 {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, 1, "events"};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T} (GeV/c)"};
    // create histograms
    histos.add("hEventCounter", "hEventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
  }

  void process(aod::Collision const& collision, 
		  soa::Join<aod::Tracks, aod::TracksExtra,
		  aod::TracksDCA> const& tracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      if( track.tpcNClsCrossedRows() < 70 ) continue;
      if( fabs(track.dcaXY()) > 0.2) continue;
      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("ptHistogram"), track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc),
    adaptAnalysisTask<myExampleTask2>(cfgc),
  };
}
