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
#ifndef DNDETA_H
#define DNDETA_H
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/MC.h"

#include "TDatabasePDG.h"

namespace o2::pwgmm::multiplicity
{
namespace
{
template <typename C, uint8_t TRACKTYPE>
inline auto collisionSelector(C const& collision)
{
  if constexpr (TRACKTYPE == aod::track::TrackTypeEnum::Run2Tracklet) {
    return collision.sel7();
  } else if constexpr (TRACKTYPE == o2::dataformats::GlobalTrackID::ITS) {
    return collision.sel8();
  } else {
    throw framework::runtime_error("Unsupported track selection!");
  }
}
} // namespace

using namespace o2::framework;

template <uint8_t TRACKTYPE>
struct PseudorapidityDensity {
  Service<TDatabasePDG> pdg;

  Configurable<float> etaMax{"etaMax", 2.0, "max eta value"};
  Configurable<float> etaMin{"etaMin", -2.0, "min eta value"};
  Configurable<float> vtxZMax{"vtxZMax", 15, "max z vertex"};
  Configurable<float> vtxZMin{"vtxZMin", -15, "min z vertex"};

  Configurable<bool> useDCA{"useDCA", false, "use DCA cuts"};
  Configurable<float> maxDCAXY{"maxDCAXY", 2.4, "max allowed transverse DCA"};
  Configurable<float> maxDCAZ{"maxDCAZ", 3.2, "max allowed longitudal DCA"};

  ConfigurableAxis percentileBinning{"pBins",
                                     {VARIABLE_WIDTH, 0., 0.01, 0.1, 0.5, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100},
                                     "Centrality/multiplicity percentile binning"};

  HistogramRegistry registry{
    "registry",
    {
      {"EventsNtrkZvtx", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}},    //
      {"TracksEtaZvtx", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}},           //
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},            //
      {"EventSelection", ";status;events", {HistType::kTH1F, {{3, 0.5, 3.5}}}},                                         //
      {"EventsNtrkZvtxGen", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}}, //
      {"TracksEtaZvtxGen", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}},        //
      {"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},         //
      {"EventEfficiency", "; status; events", {HistType::kTH1F, {{3, 0.5, 3.5}}}}                                       //
    }                                                                                                                   //
  };

  void init(InitContext&)
  {
    auto hstat = registry.get<TH1>(HIST("EventSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Selected");
    x->SetBinLabel(3, "Rejected");

    auto heff = registry.get<TH1>(HIST("EventEfficiency"));
    x = heff->GetXaxis();
    x->SetBinLabel(1, "Generated");
    x->SetBinLabel(2, "Reconstructed");
    x->SetBinLabel(3, "Selected");

    if (doprocessBinned) {
      registry.add({"EventsNtrkZvtxBin", "; N_{trk}; Z_{vtx}; Percentile", {HistType::kTH3F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}, percentileBinning}}});
      registry.add({"TracksEtaZvtxBin", "; #eta; Z_{vtx}; Percentile", {HistType::kTH3F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}, percentileBinning}}});
      registry.add({"TracksPhiEtaBin", "; #varphi; #eta; Percentile", {HistType::kTH3F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}, percentileBinning}}});
      registry.add({"EventSelectionBin", "; status; Percentile; events", {HistType::kTH2F, {{3, 0.5, 3.5}, percentileBinning}}});

      hstat = registry.get<TH2>(HIST("EventSelectionBin"));
      x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");
      x->SetBinLabel(3, "Rejected");
    }
  }

  template <typename C>
  bool select(C const& collision)
  {
    return collisionSelector<C, TRACKTYPE>(collision);
  }

  expressions::Filter etaFilter = (aod::track::eta < etaMax) && (aod::track::eta > etaMin);
  expressions::Filter trackTypeFilter = (aod::track::trackType == TRACKTYPE);
  expressions::Filter DCAFilter = ifnode(useDCA.node(), nabs(aod::track::dcaXY) <= maxDCAXY && nabs(aod::track::dcaZ) <= maxDCAZ, framework::expressions::LiteralNode{true});
  expressions::Filter posZFilter = (aod::collision::posZ < vtxZMax) && (aod::collision::posZ > vtxZMin);
  expressions::Filter posZFilterMC = (aod::mccollision::posZ < vtxZMax) && (aod::mccollision::posZ > vtxZMin);

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtended>> const& tracks)
  {
    registry.fill(HIST("EventSelection"), 1.);
    if (select(collision)) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();
      registry.fill(HIST("EventsNtrkZvtx"), tracks.size(), z);
      for (auto& track : tracks) {
        registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
        registry.fill(HIST("TracksPhiEta"), track.phi(), track.eta());
      }
    } else {
      registry.fill(HIST("EventSelection"), 3.);
    }
  }

  void processBinned(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtended>> const& tracks)
  {
    auto p = collision.centV0M();
    registry.fill(HIST("EventSelectionBin"), 1., p);
    if (select(collision)) {
      registry.fill(HIST("EventSelectionBin"), 2., p);
      auto z = collision.posZ();
      registry.fill(HIST("EventsNtrkZvtxBin"), tracks.size(), z, p);
      for (auto& track : tracks) {
        registry.fill(HIST("TracksEtaZvtxBin"), track.eta(), z, p);
        registry.fill(HIST("TracksPhiEtaBin"), track.phi(), track.eta(), p);
      }
    } else {
      registry.fill(HIST("EventSelectionBin"), 3., p);
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity<TRACKTYPE>, processBinned, "Process centrality/mult. percentile binned", false);

  using Particles = aod::McParticles;

  void processGen(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions, Particles const& particles, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtended>> const& tracks)
  {
    registry.fill(HIST("EventEfficiency"), 1.);
    for (auto& collision : collisions) {
      registry.fill(HIST("EventEfficiency"), 2.);
      if (select(collision)) {
        registry.fill(HIST("EventEfficiency"), 3.);
        auto stracks = tracks.sliceBy(aod::track::collisionId, collision.globalIndex());
        registry.fill(HIST("EventsNtrkZvtxGen"), stracks.size(), mcCollision.posZ());
      }
    }
    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      int charge = 0;
      if (p == nullptr) {
        // unknown particles will be skipped
        if (particle.pdgCode() > 1000000000) {
          //          auto x = (std::trunc(particle.pdgCode() / 10000) - 100000);
          //          charge = x - std::trunc(x / 1000) * 1000;
          LOGF(DEBUG, "[%d] Nucleus with PDG code %d", particle.globalIndex(), particle.pdgCode() /*, charge*/); // (charge %d)
        } else {
          LOGF(DEBUG, "[%d] Unknown particle with PDG code %d", particle.globalIndex(), particle.pdgCode());
        }
      } else {
        charge = p->Charge();
      }
      if (charge != 0 && MC::isPhysicalPrimary(particle) && (particle.eta() < etaMax) && (particle.eta() > etaMin)) {
        // FIXME: temporary before Run 3 MC is fixed
        if constexpr (TRACKTYPE == o2::dataformats::GlobalTrackID::ITS) {
          auto dcaxy = std::sqrt((particle.vx() - mcCollision.posX()) * (particle.vx() - mcCollision.posX()) +
                                 (particle.vy() - mcCollision.posY()) * (particle.vy() - mcCollision.posY()));
          auto dcaz = std::abs(particle.vz() - mcCollision.posZ());
          if (!(dcaxy <= maxDCAXY && dcaz <= maxDCAZ)) {
            continue;
          }
        }
        registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(), mcCollision.posZ());
        registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity<TRACKTYPE>, processGen, "Process generator-level info", false);
};
} // namespace o2::pwgmm::multiplicity

#endif // DNDETA_H
