//
// Created by mikhail on 11/30/20.
//

#include "yield.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include <AnalysisTree/DataHeader.hpp>

TASK_IMPL(Yield)

boost::program_options::options_description Yield::GetBoostOptions() {
  using namespace boost::program_options;

  options_description desc(GetName() + " options");
  desc.add_options()
      ("tracks-branch", value(&tracks_branch_)->default_value("mdc_vtx_tracks"), "Name of branch with tracks")
      ("pdg-code", value(&pdg_code_)->default_value(2212), "PDG-Code");
  return desc;
}

void Yield::PreInit() {
  SetInputBranchNames({tracks_branch_, "event_header"});
}

void Yield::Init(std::map<std::string, void *> &Map) {
  tracks_ = static_cast<AnalysisTree::Particles*>(Map.at(tracks_branch_));
  tracks_config_  = config_->GetBranchConfig(tracks_branch_);
  event_header_ = static_cast<AnalysisTree::EventHeader*>(Map.at("event_header"));
  event_header_config_  = config_->GetBranchConfig("event_header");
  float y_axis[16];
  for(int j=0; j<16; ++j){ y_axis[j]=-0.75f+0.1f* (float) j; }
  float pt_axis[]={0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
  int p=0;
  while(p<40){
    std::string histo_name{ "yield_"+std::to_string(p)+"-"+std::to_string(p+5) };
    yields_.push_back( new TH2F( histo_name.c_str(), ";y_{cm};pT [GeV/c]", 15, y_axis, 10, pt_axis ) );
    p+=5;
  }
  centrality_classes_ = new TH1F( "centrality", ";centrality", 20, 0.0, 100.0 );
  out_file_->cd();
  std::cout << "Initialized" << std::endl;
}

void Yield::Exec() {
  auto centrality = event_header_->GetField<float>( event_header_config_.GetFieldId("selected_tof_rpc_hits_centrality") );
  centrality_classes_->Fill( centrality );
  int c_class = (int) ( (centrality-2.5)/5.0 );
  auto n_tracks = tracks_->GetNumberOfChannels();
  float y_beam_2{0.74};
  auto chi2_id = tracks_config_.GetFieldId("chi2");
  auto dca_xy_id = tracks_config_.GetFieldId("dca_xy");
  auto dca_z_id = tracks_config_.GetFieldId("dca_z");
  TH2F* histo{nullptr};
  try {
    histo = yields_.at(c_class);
  } catch (std::out_of_range&) { return; }
  for (int i_track = 0; i_track < tracks_->GetNumberOfChannels(); ++i_track) {
    auto track = tracks_->GetChannel(i_track);
    auto chi2 = track.GetField<float>(chi2_id);
    if( chi2 > 100.0 )
      continue;
    auto dca_xy = track.GetField<float>(dca_xy_id);
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    auto dca_z = track.GetField<float>(dca_z_id);
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    auto pid = track.GetPid();
    if( pid != pdg_code_ )
      continue;
    auto pT = track.GetPt();
    auto ycm = track.GetRapidity() - y_beam_2;

    histo->Fill(ycm, pT);
  }
}
void Yield::Finish() {
  centrality_classes_->Write();
  for( auto histo : yields_ ) {
    histo->Sumw2();
    histo->Write();
  }
  std::cout << "Finished" << std::endl;
}
