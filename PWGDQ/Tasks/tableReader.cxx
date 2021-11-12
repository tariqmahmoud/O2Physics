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
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace o2::aod
{

namespace reducedevent
{
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
} // namespace reducedevent

namespace reducedtrack
{
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, uint8_t);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, uint8_t);
} // namespace reducedtrack

DECLARE_SOA_TABLE(EventCuts, "AOD", "EVENTCUTS", reducedevent::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "MIXINGHASHES", reducedevent::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "BARRELTRACKCUTS", reducedtrack::IsBarrelSelected);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "MUONTRACKCUTS", reducedtrack::IsMuonSelected);
} // namespace o2::aod

using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::EventCuts>;

using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID, aod::BarrelTrackCuts>;

using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonTrackCuts>;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

void SetUpMixing(MixingHandler* mixHand, TString mixVar);

// HACK: In order to be able to deduce which kind of aod object is transmitted to the templated VarManager::Fill functions
//         a constexpr static bit map must be defined and sent as template argument
//        The user has to include in this bit map all the tables needed in analysis, as defined in VarManager::ObjTypes
//        Additionally, one should make sure that the requested tables are actually provided in the process() function,
//       otherwise a compile time error will be thrown.
//        This is a temporary fix until the arrow/ROOT issues are solved, at which point it will be possible
//           to automatically detect the object types transmitted to the VarManager

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov;

constexpr int pairTypeEE = VarManager::kJpsiToEE;
constexpr int pairTypeMuMu = VarManager::kJpsiToMuMu;
constexpr int pairTypeEMu = VarManager::kElectronMuon;

struct DQEventSelection {
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  MixingHandler* fMixHandler;
  AnalysisCompositeCut* fEventCut;
  float* fValues;

  // TODO: make mixing binning configurable
  //  std::vector<float> fCentLimsHashing{0.0f, 10.0f, 20.0f, 30.0f, 50.0f, 70.0f, 90.0f};
  Configurable<std::string> fConfigMixingVariables{"cfgMixingVar", "Centrality3", "Mixing configuriation used, variables separated by a coma "};

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};

  //  int getHash(float centrality)
  //  {
  //    if ((centrality < *fCentLimsHashing.begin()) || (centrality > *(fCentLimsHashing.end() - 1))) {
  //      return -1;
  //    }
  //    auto cent = std::lower_bound(fCentLimsHashing.begin(), fCentLimsHashing.end(), centrality);
  //    return (cent - fCentLimsHashing.begin());
  //  }

  void init(o2::framework::InitContext&)
  {
    fValues = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    fMixHandler = new MixingHandler("mixingHandler", "mixingHandler");
    fMixHandler->Init();
    TString mixingVarStr = fConfigMixingVariables.value;
    SetUpMixing(fMixHandler, mixingVarStr.Data());

    DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    DefineCuts();
  }

  void DefineCuts()
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void process(MyEvents::iterator const& event)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables, fValues);

    VarManager::FillEvent<gkEventFillMap>(event, fValues);
    fHistMan->FillHistClass("Event_BeforeCuts", fValues); // automatically fill all the histograms in the class Event
    if (fEventCut->IsSelected(fValues)) {
      fHistMan->FillHistClass("Event_AfterCuts", fValues);
      eventSel(1);
    } else {
      eventSel(0);
    }
    //    int hh = getHash(fValues[VarManager::kCentVZERO]);
    int hh = fMixHandler->FindEventCategory(fValues);
    hash(hh);
  }
};

struct DQBarrelTrackSelection {
  Produces<aod::BarrelTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  std::vector<AnalysisCompositeCut> fTrackCuts;

  float* fValues; // array to be used by the VarManager

  Configurable<std::string> fConfigElectronCuts{"cfgElectronCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};

  // Configurable<std::string> fConfigElectronCuts{"cfgElectronCuts", "jpsiPID_tm", "Comma separated list of barrel track cuts"};
  // Configurable<std::string> fConfigHadronCuts{"cfgHadronCuts", "PionSelection", "Comma separated list of barrel track cuts"};

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    fValues = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    TString cutNames = "TrackBarrel_BeforeCuts;";
    for (auto& cut : fTrackCuts) {
      cutNames += Form("TrackBarrel_%s;", cut.GetName());
    }

    DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void DefineCuts()
  {
    // Defines track cuts for both the dielectron candidates and also for the dilepton - hadron task
    TString cutNamesStr = fConfigElectronCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    AnalysisCompositeCut correlationsHadronCut("hadronCut", "hadronCut", true);
    correlationsHadronCut.AddCut(dqcuts::GetAnalysisCut("electronStandardQuality"));
    correlationsHadronCut.AddCut(dqcuts::GetAnalysisCut("standardPrimaryTrack"));
    AnalysisCut hadronKine;
    hadronKine.AddCut(VarManager::kPt, 4.0, 100.0);
    correlationsHadronCut.AddCut(&hadronKine);
    fTrackCuts.push_back(correlationsHadronCut);

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void process(MyEvents::iterator const& event, MyBarrelTracks const& tracks)
  {
    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables, fValues);
    // fill event information which might be needed in histograms that combine track and event properties
    VarManager::FillEvent<gkEventFillMap>(event, fValues);

    uint8_t filterMap = uint8_t(0);

    trackSel.reserve(tracks.size());

    for (auto& track : tracks) {
      filterMap = uint8_t(0);
      VarManager::FillTrack<gkTrackFillMap>(track, fValues);
      fHistMan->FillHistClass("TrackBarrel_BeforeCuts", fValues);

      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(fValues)) {
          filterMap |= (uint8_t(1) << i);
          fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), fValues);
        }
      }
      trackSel(filterMap);
    }
  }
};

struct DQMuonTrackSelection {
  Produces<aod::MuonTrackCuts> trackSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  // NOTE: One single cut is implemented for muons, but multiple one can be computed
  AnalysisCompositeCut* fTrackCut;

  float* fValues;

  Configurable<float> fConfigMuonPtLow{"cfgMuonLowPt", 1.0f, "Low pt cut for muons"};

  void init(o2::framework::InitContext&)
  {
    fValues = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    DefineHistograms(fHistMan, "TrackMuon_BeforeCuts;TrackMuon_AfterCuts;"); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());                         // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    DefineCuts();
  }

  void DefineCuts()
  {
    fTrackCut = new AnalysisCompositeCut(true);
    AnalysisCut kineMuonCut;
    kineMuonCut.AddCut(VarManager::kPt, fConfigMuonPtLow, 100.0);
    fTrackCut->AddCut(&kineMuonCut);

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void process(MyEvents::iterator const& event, MyMuonTracks const& muons)
  {
    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables, fValues);
    VarManager::FillEvent<gkEventFillMap>(event, fValues);

    for (auto& muon : muons) {
      VarManager::FillTrack<gkMuonFillMap>(muon, fValues);
      fHistMan->FillHistClass("TrackMuon_BeforeCuts", fValues);

      if (fTrackCut->IsSelected(fValues)) {
        trackSel(uint8_t(1));
        fHistMan->FillHistClass("TrackMuon_AfterCuts", fValues);
      } else {
        trackSel(uint8_t(0));
      }
    }
  }
};

struct DQEventMixing {
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  float* fValues;
  // NOTE: The track filter produced by the barrel track selection contain a number of electron cut decisions and one last cut for hadrons used in the
  //           dilepton - hadron task downstream. So the bit mask is required to select pairs just based on the electron cuts
  uint8_t fTwoTrackFilterMask = 0;
  std::vector<TString> fCutNames;

  Configurable<std::string> fConfigElectronCuts{"cfgElectronCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};

  Filter filterEventSelected = aod::reducedevent::isEventSelected == 1;
  Filter filterTrackSelected = aod::reducedtrack::isBarrelSelected > uint8_t(0);
  Filter filterMuonTrackSelected = aod::reducedtrack::isMuonSelected > uint8_t(0);

  void init(o2::framework::InitContext&)
  {
    fValues = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    TString histNames = "";
    TString configCutNamesStr = fConfigElectronCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fCutNames.push_back(objArray->At(icut)->GetName());
        histNames += Form("PairsBarrelMEPM_%s;PairsBarrelMEPP_%s;PairsBarrelMEMM_%s;", fCutNames[icut].Data(), fCutNames[icut].Data(), fCutNames[icut].Data());
        histNames += Form("PairsEleMuMEPM_%s;PairsEleMuMEPP_%s;PairsEleMuMEMM_%s;", fCutNames[icut].Data(), fCutNames[icut].Data(), fCutNames[icut].Data());
        fTwoTrackFilterMask |= (uint8_t(1) << icut);
      }
    }
    histNames += "PairsMuonMEPM;PairsMuonMEPP;PairsMuonMEMM;";

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void process(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    uint8_t twoTrackFilter = 0;

    events.bindExternalIndices(&tracks);
    events.bindExternalIndices(&muons);
    auto tracksTuple = std::make_tuple(tracks);
    auto muonsTuple = std::make_tuple(muons);
    AnalysisDataProcessorBuilder::GroupSlicer slicerTracks(events, tracksTuple);
    AnalysisDataProcessorBuilder::GroupSlicer slicerMuons(events, muonsTuple);

    // Strictly upper categorised collisions, for 100 combinations per bin, skipping those in entry -1
    for (auto& [event1, event2] : selfCombinations("fMixingHash", 100, -1, events, events)) {

      // event informaiton is required to fill histograms where both event and pair information is required (e.g. inv.mass vs centrality)
      VarManager::ResetValues(0, VarManager::kNVars, fValues);
      VarManager::FillEvent<gkEventFillMap>(event1, fValues);

      auto it1 = slicerTracks.begin();
      auto it2 = slicerTracks.begin();
      for (auto& slice : slicerTracks) {
        if (slice.groupingElement().index() == event1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : slicerTracks) {
        if (slice.groupingElement().index() == event2.index()) {
          it2 = slice;
          break;
        }
      }

      auto tracks1 = std::get<soa::Filtered<MyBarrelTracksSelected>>(it1.associatedTables());
      tracks1.bindExternalIndices(&events);
      auto tracks2 = std::get<soa::Filtered<MyBarrelTracksSelected>>(it2.associatedTables());
      tracks2.bindExternalIndices(&events);

      for (auto& track1 : tracks1) {
        for (auto& track2 : tracks2) {
          twoTrackFilter = track1.isBarrelSelected() & track2.isBarrelSelected() & fTwoTrackFilterMask;

          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }

          VarManager::FillPair<pairTypeEE>(track1, track2, fValues);

          int i = 0;
          for (auto cutName = fCutNames.begin(); cutName != fCutNames.end(); cutName++, i++) {
            if (twoTrackFilter & (uint8_t(1) << i)) {
              if (track1.sign() * track2.sign() < 0) {
                fHistMan->FillHistClass(Form("PairsBarrelMEPM_%s", (*cutName).Data()), fValues);
              } else {
                if (track1.sign() > 0) {
                  fHistMan->FillHistClass(Form("PairsBarrelMEPP_%s", (*cutName).Data()), fValues);
                } else {
                  fHistMan->FillHistClass(Form("PairsBarrelMEMM_%s", (*cutName).Data()), fValues);
                }
              }
            } // end if (filter bits)
          }   // end for (cuts)
        }     // end for (track2)
      }       // end for (track1)

      auto im1 = slicerMuons.begin();
      auto im2 = slicerMuons.begin();
      for (auto& slice : slicerMuons) {
        if (slice.groupingElement().index() == event1.index()) {
          im1 = slice;
          break;
        }
      }
      for (auto& slice : slicerMuons) {
        if (slice.groupingElement().index() == event2.index()) {
          im2 = slice;
          break;
        }
      }

      auto muons1 = std::get<soa::Filtered<MyMuonTracksSelected>>(im1.associatedTables());
      muons1.bindExternalIndices(&events);
      auto muons2 = std::get<soa::Filtered<MyMuonTracksSelected>>(im2.associatedTables());
      muons2.bindExternalIndices(&events);

      for (auto& muon1 : muons1) {
        for (auto& muon2 : muons2) {
          twoTrackFilter = muon1.isMuonSelected() & muon2.isMuonSelected();
          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }
          VarManager::FillPair<pairTypeMuMu>(muon1, muon2, fValues);
          if (muon1.sign() * muon2.sign() < 0) {
            fHistMan->FillHistClass("PairsMuonMEPM", fValues);
          } else {
            if (muon1.sign() > 0) {
              fHistMan->FillHistClass("PairsMuonMEPP", fValues);
            } else {
              fHistMan->FillHistClass("PairsMuonMEMM", fValues);
            }
          }
        } // end for (muon2)
      }   // end for (muon1)

      for (auto& track1 : tracks1) {
        for (auto& muon2 : muons2) {
          twoTrackFilter = (track1.isBarrelSelected() & fTwoTrackFilterMask) && muon2.isMuonSelected();

          if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
            continue;
          }

          VarManager::FillPair<pairTypeEMu>(track1, muon2, fValues);

          int i = 0;
          for (auto cutName = fCutNames.begin(); cutName != fCutNames.end(); cutName++, i++) {
            if (twoTrackFilter & (uint8_t(1) << i)) {
              if (track1.sign() * muon2.sign() < 0) {
                fHistMan->FillHistClass(Form("PairsEleMuMEPM_%s", (*cutName).Data()), fValues);
              } else {
                if (track1.sign() > 0) {
                  fHistMan->FillHistClass(Form("PairsEleMuMEPP_%s", (*cutName).Data()), fValues);
                } else {
                  fHistMan->FillHistClass(Form("PairsEleMuMEMM_%s", (*cutName).Data()), fValues);
                }
              }
            } // end if (filter bits)
          }   // end for (cuts)
        }     // end for (muon2)
      }       // end for (track1)
    }         // end for (event combinations)
  }           // end process()
};

struct DQTableReader {
  Produces<aod::Dileptons> dileptonList;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  //NOTE: one could define also a dilepton cut, but for now basic selections can be supported using Partition

  float* fValues;
  // NOTE: The track filter produced by the barrel track selection contain a number of electron cut decisions and one last cut for hadrons used in the
  //           dilepton - hadron task downstream. So the bit mask is required to select pairs just based on the electron cuts
  uint8_t fTwoTrackFilterMask = 0;
  std::vector<TString> fCutNames;

  Filter filterEventSelected = aod::reducedevent::isEventSelected == 1;
  // NOTE: the barrel filter map contains decisions for both electrons and hadrons used in the correlation task
  Filter filterBarrelTrackSelected = aod::reducedtrack::isBarrelSelected > uint8_t(0);
  Filter filterMuonTrackSelected = aod::reducedtrack::isMuonSelected > uint8_t(0);
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> fConfigElectronCuts{"cfgElectronCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};

  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  void init(o2::framework::InitContext&)
  {
    fValues = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    TString histNames = "";
    TString configCutNamesStr = fConfigElectronCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fCutNames.push_back(objArray->At(icut)->GetName());
        histNames += Form("PairsBarrelSEPM_%s;PairsBarrelSEPP_%s;PairsBarrelSEMM_%s;", fCutNames[icut].Data(), fCutNames[icut].Data(), fCutNames[icut].Data());
        histNames += Form("PairsEleMuSEPM_%s;PairsEleMuSEPP_%s;PairsEleMuSEMM_%s;", fCutNames[icut].Data(), fCutNames[icut].Data(), fCutNames[icut].Data());
        fTwoTrackFilterMask |= (uint8_t(1) << icut);
      }
    }
    histNames += "PairsMuonSEPM;PairsMuonSEPP;PairsMuonSEMM;";
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(nolaterthan.value);
    auto histCCDB = ccdb->get<TH1F>(ccdbPath.value);
    if (!histCCDB) {
      LOGF(fatal, "CCDB histogram not found");
    }

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true); // TODO: get these parameters from Configurables
    VarManager::SetupTwoProngFwdDCAFitter(5.0f, true, 200.0f, 1.0e-3f, 0.9f, true);
  }

  void process(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    if (!event.isEventSelected()) {
      return;
    }
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNVars, fValues);
    VarManager::FillEvent<gkEventFillMap>(event, fValues);

    // Run the same event pairing for barrel tracks
    uint8_t twoTrackFilter = 0;
    uint32_t dileptonFilterMap = 0;
    for (auto& [t1, t2] : combinations(tracks, tracks)) {
      twoTrackFilter = t1.isBarrelSelected() & t2.isBarrelSelected() & fTwoTrackFilterMask;
      if (!twoTrackFilter) { // the tracks must have at least one filter bit in common to continue
        continue;
      }
      dileptonFilterMap = uint32_t(twoTrackFilter);
      VarManager::FillPair<pairTypeEE>(t1, t2, fValues);
      VarManager::FillPairVertexing<pairTypeEE>(event, t1, t2, fValues);
      dileptonList(event, fValues[VarManager::kMass], fValues[VarManager::kPt], fValues[VarManager::kEta], fValues[VarManager::kPhi], t1.sign() + t2.sign(), dileptonFilterMap);
      int i = 0;
      for (auto cutName = fCutNames.begin(); cutName != fCutNames.end(); cutName++, i++) {
        if (twoTrackFilter & (uint8_t(1) << i)) {
          if (t1.sign() * t2.sign() < 0) {
            fHistMan->FillHistClass(Form("PairsBarrelSEPM_%s", (*cutName).Data()), fValues);
          } else {
            if (t1.sign() > 0) {
              fHistMan->FillHistClass(Form("PairsBarrelSEPP_%s", (*cutName).Data()), fValues);
            } else {
              fHistMan->FillHistClass(Form("PairsBarrelSEMM_%s", (*cutName).Data()), fValues);
            }
          }
        }
      }
    } // end loop over barrel track pairs

    // same event pairing for muons
    for (auto& [muon1, muon2] : combinations(muons, muons)) {
      twoTrackFilter = muon1.isMuonSelected() & muon2.isMuonSelected();
      if (!twoTrackFilter) { // the muons must have at least one filter bit in common to continue
        continue;
      }
      // NOTE: The dimuons and electron-muon pairs in this task are pushed in the same table as the dielectrons.
      //       In order to discriminate them, the dileptonFilterMap uses the first 8 bits for dielectrons, the next 8 for dimuons and the rest for electron-muon
      // TBD:  Other implementations may be possible, for example add a column to the dilepton table to specify the pair type (dielectron, dimuon, electron-muon, etc.)
      dileptonFilterMap = uint32_t(twoTrackFilter) << 8;
      VarManager::FillPair<pairTypeMuMu>(muon1, muon2, fValues);
      VarManager::FillPairVertexing<pairTypeMuMu>(event, muon1, muon2, fValues);
      dileptonList(event, fValues[VarManager::kMass], fValues[VarManager::kPt], fValues[VarManager::kEta], fValues[VarManager::kPhi], muon1.sign() + muon2.sign(), dileptonFilterMap);
      if (muon1.sign() * muon2.sign() < 0) {
        fHistMan->FillHistClass("PairsMuonSEPM", fValues);
      } else {
        if (muon1.sign() > 0) {
          fHistMan->FillHistClass("PairsMuonSEPP", fValues);
        } else {
          fHistMan->FillHistClass("PairsMuonSEMM", fValues);
        }
      }
    } // end loop over muon track pairs

    for (auto& [trackBarrel, trackMuon] : combinations(tracks, muons)) {
      twoTrackFilter = (trackBarrel.isBarrelSelected() & fTwoTrackFilterMask) && trackMuon.isMuonSelected();
      if (!twoTrackFilter) { // the muon and barrel track must have at least one filter bit in common to continue
        continue;
      }
      // NOTE: The dimuons and electron-muon pairs in this task are pushed in the same table as the dielectrons.
      //       In order to discriminate them, the dileptonFilterMap uses the first 8 bits for dielectrons, the next 8 for dimuons and the rest for electron-muon
      dileptonFilterMap = uint32_t(twoTrackFilter) << 16;
      VarManager::FillPair<pairTypeEMu>(trackBarrel, trackMuon, fValues);
      dileptonList(event, fValues[VarManager::kMass], fValues[VarManager::kPt], fValues[VarManager::kEta], fValues[VarManager::kPhi], trackBarrel.sign() + trackMuon.sign(), dileptonFilterMap);
      int i = 0;
      for (auto cutName = fCutNames.begin(); cutName != fCutNames.end(); cutName++, i++) {
        if (twoTrackFilter & (uint8_t(1) << i)) {
          if (trackBarrel.sign() * trackMuon.sign() < 0) {
            fHistMan->FillHistClass(Form("PairsEleMuSEPM_%s", (*cutName).Data()), fValues);
          } else {
            if (trackBarrel.sign() > 0) {
              fHistMan->FillHistClass(Form("PairsEleMuSEPP_%s", (*cutName).Data()), fValues);
            } else {
              fHistMan->FillHistClass(Form("PairsEleMuSEMM_%s", (*cutName).Data()), fValues);
            }
          }
        } // end if (filter bits)
      }   // end for (cuts)
    }     // end loop over electron-muon pairs
  }
};

struct DQDileptonHadronAnalysis {
  //
  // This task combines dilepton candidates with a track and could be used for example
  //  in analyses with the dilepton as one of the decay products of a higher mass resonance (e.g. B0 -> Jpsi + K)
  //    or in dilepton + hadron correlations, etc.
  // It requires the TableReader task to be in the workflow and produce the dilepton table
  //
  //  The barrel and muon track filtering tasks can produce multiple parallel decisions, which are used to produce
  //   dileptons which inherit the logical intersection of the track level decisions (see the TableReader task). This can be used
  //   also in the dilepton-hadron correlation analysis. However, in this model of the task, we use all the dileptons produced in the
  //     TableReader task to combine them with the hadrons selected by the barrel track selection.
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesHadron;

  Filter eventFilter = aod::reducedevent::isEventSelected == 1;
  Filter dileptonFilter = aod::reducedpair::mass > 2.92f && aod::reducedpair::mass<3.16f && aod::reducedpair::pt> 0.0f && aod::reducedpair::sign == 0;
  // NOTE: the barrel track filter is shared between the filters for dilepton electron candidates (first n-bits)
  //       and the associated hadrons (n+1 bit) --> see the barrel track selection task
  //      The current condition should be replaced when bitwise operators will become available in Filter expresions
  // NOTE: the name of this configurable has to be the same as the one used in the barrel track selection task
  Configurable<std::string> fConfigElectronCuts{"cfgElectronCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  int fNHadronCutBit;

  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair;

  void init(o2::framework::InitContext&)
  {
    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    DefineHistograms(fHistMan, "DileptonsSelected;DileptonHadronInvMass;DileptonHadronCorrelation"); // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    TString configCutNamesStr = fConfigElectronCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      fNHadronCutBit = objArray->GetEntries();
    } else {
      fNHadronCutBit = 0;
    }
  }

  void process(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelected const& hadrons, soa::Filtered<aod::Dileptons> const& dileptons)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    // fill event information which might be needed in histograms that combine track/pair and event properties
    VarManager::FillEvent<gkEventFillMap>(event, fValuesHadron);
    VarManager::FillEvent<gkEventFillMap>(event, fValuesDilepton); // TODO: check if needed (just for dilepton QA which might be depending on event wise variables)

    // loop once over dileptons for QA purposes
    for (auto dilepton : dileptons) {
      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton);
      fHistMan->FillHistClass("DileptonsSelected", fValuesDilepton);
    }

    // loop over hadrons
    for (auto& hadron : hadrons) {
      if (!(hadron.isBarrelSelected() & (uint8_t(1) << fNHadronCutBit))) {
        continue;
      }
      for (auto dilepton : dileptons) {
        // TODO: At the moment there is no check on whether this hadron is one of the dilepton daughters!
        VarManager::FillDileptonHadron(dilepton, hadron, fValuesHadron);
        fHistMan->FillHistClass("DileptonHadronInvMass", fValuesHadron);
        fHistMan->FillHistClass("DileptonHadronCorrelation", fValuesHadron);
      }
    }
  }
};

//////////////////////////////////
///////////// tariq exotics
struct DQDileptonDiHadronAnalysis {

  int fAllEvents=0;  //int fEventsPassingH1=0; int fEventsPassingH2=0; 
  int fEventsPassingHH=0; int fEventsPassingLL=0; int fEventsPassingLLHH=0;int fSelEvents=0;

  int fNLL=0;  int fNLLSel=0;  int fNHH=0;  int fNHHSel=0; int fNH1=0;  int fNH1Sel=0;  int fNH2=0;  int fNH2Sel=0; 
  
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  
  Partition<MyBarrelTracksSelected> PosHadrons = aod::reducedtrack::sign > 0;// && aod::reducedtrack::isBarrelSelected > uint8_t(0);
  Partition<MyBarrelTracksSelected> NegHadrons = aod::reducedtrack::sign < 0;// && aod::reducedtrack::isBarrelSelected > uint8_t(0);
  Partition<MyBarrelTracksSelected> NeuHadrons = aod::reducedtrack::sign == 0;// && aod::reducedtrack::isBarrelSelected > uint8_t(0);
   
  // use two values array to avoid mixing up the quantities
  float* fValuesStat;
  float* fValuesDilepton;
  float* fValuesHadron;//no need
  float* fValuesDilepDihad;
  
  Filter eventFilter = aod::reducedevent::isEventSelected == 1;
  Filter dileptonFilter = aod::reducedpair::mass > 0.0f && aod::reducedpair::mass<316.0f && aod::reducedpair::pt> 0.0f && aod::reducedpair::sign == 0;
  Filter hadronFilter = aod::reducedtrack::pt> 0.1 && aod::reducedtrack::sign != 0;
  
  Configurable<std::string> fConfigElectronCuts{"cfgElectronCuts", "jpsiPID_tm", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigPionCuts{"cfgPionCuts", "PionSelection", "Comma separated list of barrel track cuts"};//move to trackselection
  
  Configurable<double> fConfigAlphaJPsiMotherUpperLimit{"cfgOpAngUL", 0.2, "upper limit of openingangle llhh-jpsi "};
  Configurable<double> fConfigLLMassLowerLimit{"cfgLLMassLowerLimit", 2.92, "lower limit of ll mass "};
  Configurable<double> fConfigLLMassUpperLimit{"cfgLLMassUpperLimit", 3.16, "upper limit of ll mass "};
  
  /////////////////////////////////////////
  bool Debug=kFALSE;

  int fNHadronCutBit;
  
  constexpr static uint32_t fgDileptonFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::Pair;
  constexpr static uint32_t fgHadronFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;
  
  
  void init(o2::framework::InitContext&)
  {
    fValuesStat = new float[VarManager::kNVars];
    fValuesDilepton = new float[VarManager::kNVars];
    fValuesHadron = new float[VarManager::kNVars];
    fValuesDilepDihad = new float[VarManager::kNVars];
    
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
    
    
    ///// FIXME: for the moment define hists here. Not with HistMan.
    DefineHistograms(fHistMan, "Statistics;Hadrons;Hadron1;Hadron2;HadronsSelected;Hadron1Selected;Hadron2Selected;HadronsII;HadronsSelectedII;Dilepton;DileptonSelected;Dihadron;DihadronSelected;DileptonDihadron"); // define all histograms
    
    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());
    
    fHistMan->AddHistogram("Statistics" , "h_nAllEvents", "# og all events", false, 8, 0.0, 8.0, VarManager::kNAllEvents); 
    fHistMan->AddHistogram("Statistics" , "h_nEventsPassingLL", "# og all events passing ll", false, 8, 0.0, 8.0, VarManager::kNEventsPassingLL); 

    fHistMan->AddHistogram("Statistics" , "h_nLL", "# ll", false, 10, 0.0, 10.0, VarManager::kNLL); 
    fHistMan->AddHistogram("Statistics" , "h_nLLSel", "# selected ll ", false, 10, 0.0, 10.0, VarManager::kNLLSel); 

    fHistMan->AddHistogram("Statistics" , "h_nH1", "# pos H", false, 50, 0.0, 50.0, VarManager::kNH1); 
    fHistMan->AddHistogram("Statistics" , "h_nH1Sel", "# Sel pos H", false, 50, 0.0, 50.0, VarManager::kNH1Sel); 
    fHistMan->AddHistogram("Statistics" , "h_nH2", "# pos H", false, 50, 0.0, 50.0, VarManager::kNH2); 
    fHistMan->AddHistogram("Statistics" , "h_nH2Sel", "# Sel pos H", false, 50, 0.0, 50.0, VarManager::kNH2Sel); 

    fHistMan->AddHistogram("Statistics" , "h_nHH", "# hh", false, 100, 0.0, 100.0, VarManager::kNHH); 
    fHistMan->AddHistogram("Statistics" , "h_nHHSel", "# selected hh ", false, 100, 0.0, 100.0, VarManager::kNHHSel); 


    /////////////////////////////
    fHistMan->AddHistogram("Hadrons" , "h_pt", "pt of h", false, 200, 0.0, 20.0, VarManager::kPt); 
    fHistMan->AddHistogram("Hadrons" , "h_eta", "eta of h", false, 100, -1.0, 1.0, VarManager::kEta); 
    fHistMan->AddHistogram("Hadron1" , "h1_pt", "pt of h1", false, 200, 0.0, 20.0, VarManager::kPt); 
    fHistMan->AddHistogram("Hadron1" , "h1_eta", "eta of h1", false, 100, -1.0, 1.0, VarManager::kEta); 
    fHistMan->AddHistogram("Hadron2" , "h2_pt", "pt of h2", false, 200, 0.0, 20.0, VarManager::kPt); 
    fHistMan->AddHistogram("Hadron2" , "h2_eta", "eta of h2", false, 100, -1.0, 1.0, VarManager::kEta); 


    fHistMan->AddHistogram("HadronsSelected" , "h_pt", "pt of h", false, 200, 0.0, 20.0, VarManager::kPt); 
    fHistMan->AddHistogram("HadronsSelected" , "h_eta", "eta of sel h", false, 100, -1.0, 1.0, VarManager::kEta); 
    fHistMan->AddHistogram("Hadron1Selected" , "h1_pt", "pt of h1", false, 200, 0.0, 20.0, VarManager::kPt); 
    fHistMan->AddHistogram("Hadron1Selected" , "h1_eta", "eta of sel h1", false, 100, -1.0, 1.0, VarManager::kEta); 
    fHistMan->AddHistogram("Hadron2Selected" , "h2_pt", "pt of h2", false, 200, 0.0, 20.0, VarManager::kPt); 
    fHistMan->AddHistogram("Hadron2Selected" , "h2_eta", "eta of sel h2", false, 100, -1.0, 1.0, VarManager::kEta); 
    //what is the pt of the ll pair: kPt or kPairPt?
    fHistMan->AddHistogram("Dilepton" , "ll_pt", "pt of ll", false, 200, 0.0, 20.0, VarManager::kPt);     
    fHistMan->AddHistogram("Dilepton" , "ll_mass", "mass of ll", false, 500, 0.0, 5.0, VarManager::kMass);     

    fHistMan->AddHistogram("DileptonSelected" , "ll_sel_pt", "pt of sel ll", false, 200, 0.0, 20.0, VarManager::kPt);     
    fHistMan->AddHistogram("DileptonSelected" , "ll_sel_mass", "mass of sel ll", false, 500, 0.0, 5.0, VarManager::kMass);     
    //DileptonSelected
    //FIXME: for the moment hh and hhsel are the same
    fHistMan->AddHistogram("Dihadron" , "diHadron_pt", "pt of dihadron", false, 200, 0.0, 20.0, VarManager::kHHPt);     
    fHistMan->AddHistogram("Dihadron" , "diHadron_mass", "inv. mass of dihadron", false, 500, 0.0, 10.0, VarManager::kHHMass);         
    fHistMan->AddHistogram("DihadronSelected" , "diHadronSel_pt", "pt of sel. dihadron", false, 200, 0.0, 20.0, VarManager::kHHPt);     
    fHistMan->AddHistogram("DihadronSelected" , "diHadronSel_mass", "inv. mass of sel. dihadron", false, 50, 0.0, 5.0, VarManager::kHHMass);     
    fHistMan->AddHistogram("DileptonDihadron" , "llhh_pt", "pt of llhh", false, 200, 0.0, 20.0, VarManager::kLLHHPt);     
    fHistMan->AddHistogram("DileptonDihadron" , "llhh_eta", "eta of llhh", false, 600, -1.0, 5.0, VarManager::kLLHHEta);     
    fHistMan->AddHistogram("DileptonDihadron" , "llhh_mass", "inv. mass of llhh", false, 1000, 0.0, 20.0, VarManager::kLLHHMass);     
    
    fHistMan->AddHistogram("DileptonDihadron" , "angle_ll", "angle between llhh and ll", false, 630, 0.0, 6.3, VarManager::kAlphaLL);     
    fHistMan->AddHistogram("DileptonDihadron" , "angle_h1", "angle between llhh and h1", false, 630, 0.0, 6.3, VarManager::kAlphaH1);     
    fHistMan->AddHistogram("DileptonDihadron" , "angle_h2", "angle between llhh and h2", false, 630, 0.0, 6.3, VarManager::kAlphaH2);     

    fHistMan->AddHistogram("DileptonDihadron" , "angle_llh1", "angle between ll and h1", false, 630, 0.0, 6.3, VarManager::kAlphaLLH1);     
    fHistMan->AddHistogram("DileptonDihadron" , "angle_llh2", "angle between ll and h2", false, 630, 0.0, 6.3, VarManager::kAlphaLLH2);     
    /////////////////////
    /////////////// test histograms plot hadron pt and lepton pt only within the jpsi-mass region
    fHistMan->AddHistogram("HadronsII" , "h_pt", "pt of h", false, 200, 0.0, 20.0, VarManager::kPt); 
    fHistMan->AddHistogram("DileptonDihadron" , "h_pt pos pi", "pt of pos pi", false, 200, 0.0, 20.0, VarManager::kAssociatedPt); 

    fHistMan->AddHistogram("HadronsSelectedII" , "h_pt", "pt of h", false, 200, 0.0, 20.0, VarManager::kPt); 
    fHistMan->AddHistogram("DileptonDihadron" , "h_pt neg pi", "pt of neg pi", false, 200, 0.0, 20.0, VarManager::kAssociated2Pt); 

    fHistMan->AddHistogram("DileptonDihadron" , "h_TrigMass", "Mass of trigger ll", false, 500, 0.0, 5.0, VarManager::kMassTrigger); //29    
    fHistMan->AddHistogram("DileptonDihadron" , "h_TrigMass2", "Mass of trigger ll", false, 1000, 2.5, 3.5, VarManager::kMassTrigger); //29    

 
    TString configCutNamesStr = fConfigElectronCuts.value;
    if (!configCutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(configCutNamesStr.Tokenize(","));
      fNHadronCutBit = objArray->GetEntries();
    } else {
      fNHadronCutBit = 0;
    }


  }
  
  //___________________________________________________
  bool IsDiLeptonSelected(double llMass){  
    //if(llMass < 2.92f || llMass > 3.16f){return false;} && aod::reducedpair::pt> 0.0f && aod::reducedpair::sign == 0;
    if(llMass < fConfigLLMassLowerLimit || llMass > fConfigLLMassUpperLimit){return false;} 
    return true;
  }
  //___________________________________________________
  void FillStatistics(int fAllEvents, int fEventsPassingLL, int fEventsPassingLLHH,int fSelEvents,  int fNLL,  int fNLLSel,  int fNHH,  int fNHHSel, int fNH1,  int fNH1Sel,  int fNH2,  int fNH2Sel, float* values){
    values[VarManager::kNAllEvents]=float(fAllEvents);
    values[VarManager::kNEventsPassingLL]=float(fEventsPassingLL);
    values[VarManager::kNLL]=float(fNLL);
    values[VarManager::kNLLSel]=float(fNLLSel);
    values[VarManager::kNH1]=float(fNH1);
    values[VarManager::kNH1Sel]=float(fNH1Sel);
    values[VarManager::kNH2]=float(fNH2);
    values[VarManager::kNH2Sel]=float(fNH2Sel);
    values[VarManager::kNHH]=float(fNHH);
    values[VarManager::kNHHSel]=float(fNHHSel); 
  }  
  
  //___________________________________________________
  void process(soa::Filtered<MyEventsVtxCovSelected>::iterator const& event, MyBarrelTracksSelected const& hadrons, soa::Filtered<aod::Dileptons> const& dileptons)
  {

    fAllEvents=0;  fEventsPassingLL=0;  fEventsPassingHH=0;  fEventsPassingLLHH=0; fSelEvents=0;
    //FIXME: does it make sence to do this here?
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepDihad);
    
    // fill event information which might be needed in histograms that combine track/pair and event properties
    VarManager::FillEvent<gkEventFillMap>(event, fValuesDilepDihad);
    VarManager::FillEvent<gkEventFillMap>(event, fValuesDilepton); // TODO: check if needed (just for dilepton QA which might be depending on event wise variables)
    VarManager::FillEvent<gkEventFillMap>(event, fValuesHadron); // 
    
    //////////////////////////////////////
    fAllEvents++;
    /////////////////////////////////////////////////////
    // loop once over dileptons for QA purposes
    
    fNLL=0; fNLLSel=0;
    bool HadronsFilled=0;
    for (auto& dilepton : dileptons) {
      VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton); 
      fHistMan->FillHistClass("Dilepton", fValuesDilepton);
      fNLL++;
       if(!IsDiLeptonSelected(dilepton.mass() )){continue;}
      fHistMan->FillHistClass("DileptonSelected", fValuesDilepton);
 
      fNLLSel++;
      if(!HadronsFilled){
	for (auto& hadron : hadrons) {
	  VarManager::FillTrack<gkTrackFillMap>(hadron, fValuesHadron); 
	  fHistMan->FillHistClass("HadronsII", fValuesHadron);
	  if (!(hadron.isBarrelSelected() & (uint8_t(1) << fNHadronCutBit))) {continue;}
	  fHistMan->FillHistClass("HadronsSelectedII", fValuesHadron);
	  HadronsFilled=1;
	}
      }
    }
    
    
    fNH1=0; fNH1Sel=0;  
    for (auto& hadron : PosHadrons) {
       VarManager::FillTrack<gkTrackFillMap>(hadron, fValuesHadron); 
      
      if(0){cout<<" nsignal = "<<fValuesHadron[VarManager::kTPCsignal]<<" nsigPi ="<<fValuesHadron[VarManager::kTPCnSigmaPi] << endl;}
      
      fHistMan->FillHistClass("Hadrons", fValuesHadron);
      fHistMan->FillHistClass("Hadron1", fValuesHadron);
      fNH1++;  
      
      if (!(hadron.isBarrelSelected() & (uint8_t(1) << fNHadronCutBit))) {continue;}
      fHistMan->FillHistClass("HadronsSelected", fValuesHadron);
      fHistMan->FillHistClass("Hadron1Selected", fValuesHadron);
      fNH1Sel++;
    }
    
    fNH2=0;  fNH2Sel=0; 
    for (auto& hadron : NegHadrons) {
      VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
      VarManager::FillTrack<gkTrackFillMap>(hadron, fValuesHadron); 
      fHistMan->FillHistClass("Hadrons", fValuesHadron);
      fHistMan->FillHistClass("Hadron2", fValuesHadron);
      fNH2++;  
      
      if (!(hadron.isBarrelSelected() & (uint8_t(1) << fNHadronCutBit))) {continue;}
      fHistMan->FillHistClass("HadronsSelected", fValuesHadron);
      fHistMan->FillHistClass("Hadron2Selected", fValuesHadron);
      fNH2Sel++;
    }
    
    if(fNH1Sel > 0 && fNH2Sel >0 ){fEventsPassingHH++;}//FIXME: not exactly they must pass HHSel, see below in the loops
 
    for (auto& dilepton : dileptons) {
      VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
      VarManager::FillTrack<fgDileptonFillMap>(dilepton, fValuesDilepton); 
      if(!IsDiLeptonSelected(dilepton.mass() )){continue;}
      
      for (auto& hadron1 : PosHadrons) {
	if (!(hadron1.isBarrelSelected() & (uint8_t(1) << fNHadronCutBit))) {  continue;}
	//////////loop over 2. had
	for (auto& hadron2 : NegHadrons) {
	  if (!(hadron2.isBarrelSelected() & (uint8_t(1) << fNHadronCutBit))) {continue;}
	  // FIXME: At the moment there is no check on whether this hadron is one of the dilepton daughters!

	  VarManager::FillDileptonDihadron(dilepton, hadron1, hadron2, fValuesDilepDihad);
	  fHistMan->FillHistClass("Dihadron", fValuesDilepDihad);
	  fHistMan->FillHistClass("DihadronSelected", fValuesDilepDihad);//same as Dihadron
	  fHistMan->FillHistClass("DileptonDihadron", fValuesDilepDihad);
	  fNHHSel++; //same as fEventsPassingLLHH //FIXME:: double counting with n-dill >1 -
	}//neg hadrons loop
      } //pos hadrons loop
    }//dilep loop
    /////////////////////////////////////////
    fEventsPassingLLHH++;
    
    VarManager::ResetValues(0, VarManager::kNVars, fValuesStat);
    FillStatistics( 
		   fAllEvents,  fEventsPassingLL,  fEventsPassingLLHH, fSelEvents,   fNLL,  fNLLSel, fNHH, fNHHSel, fNH1, fNH1Sel,  fNH2, fNH2Sel,
		   fValuesStat);
    fHistMan->FillHistClass("Statistics", fValuesStat);
    
  }
};

///////////// end tariq exotics
//////////////////////////////////

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQEventSelection>(cfgc),
      adaptAnalysisTask<DQBarrelTrackSelection>(cfgc),
      adaptAnalysisTask<DQEventMixing>(cfgc),
      adaptAnalysisTask<DQMuonTrackSelection>(cfgc),
      adaptAnalysisTask<DQTableReader>(cfgc),
      adaptAnalysisTask<DQDileptonHadronAnalysis>(cfgc),
      adaptAnalysisTask<DQDileptonDiHadronAnalysis>(cfgc, TaskName{"d-q-xcand-creator"} )};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  The histogram classes and their components histograms are defined below depending on the name of the histogram class
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,cent,muon");
    }

    if (classStr.Contains("Track")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca,tofpid");
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
      }
    }

    if (classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_barrel", "vertexing-barrel");
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_dimuon", "vertexing-forward");
      }
      if (classStr.Contains("EleMu")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair_electronmuon");
      }
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair");
    }

    if (classStr.Contains("HadronsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "kine");
    }

    if (classStr.Contains("DileptonHadronInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }

    if (classStr.Contains("DileptonHadronCorrelation")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-correlation");
    }
  } // end loop over histogram classes
}

void SetUpMixing(MixingHandler* mixHand, TString mixVar)
{
  std::unique_ptr<TObjArray> objArray(mixVar.Tokenize(","));
  for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
    dqmixing::SetUpMixing(mixHand, objArray->At(iVar)->GetName());
  } // end loop over histogram classes
}
