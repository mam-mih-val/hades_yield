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
  return desc;
}

void Yield::PreInit() {
  SetInputBranchNames({"mdc_vtx_tracks", "event_header", "sim_tracks"});
}

void Yield::UserInit(std::map<std::string, void *> &Map) {
  event_header_ = GetInBranch("event_header");
  centrality_ = GetVar("event_header/selected_tof_rpc_hits_centrality");

  tracks_ = GetInBranch("mdc_vtx_tracks");
  pdg_code_ = GetVar("mdc_vtx_tracks/pid");
  dca_xy_ = GetVar("mdc_vtx_tracks/dca_xy");
  dca_z_ = GetVar("mdc_vtx_tracks/dca_z");
  chi2_ = GetVar("mdc_vtx_tracks/chi2");

  sim_particles_ = GetInBranch("sim_tracks");
  sim_pdg_code_ = GetVar("sim_tracks/pid");
  is_primary_ = GetVar("sim_tracks/is_primary");

  float y_axis[16];
  for(int j=0; j<16; ++j){ y_axis[j]=-0.75f+0.1f* (float) j; }
  float pt_axis[]={0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
  int p=0;
  while(p<40){
    std::string histo_name{ "yield_"+std::to_string(p)+"-"+std::to_string(p+5) };
    yields_.push_back( new TH2F( histo_name.c_str(), ";y_{cm};pT [GeV/c]", 15, y_axis, 10, pt_axis ) );
    p+=5;
  }
  pt_distribution_reco_ = new TH1F( "pT_distribution_reco", ";p_{T} [GeV/c]; entries", 2000, 0, 2.0 );
  pt_distribution_sim_ = new TH1F( "pT_distribution_sim", ";p_{T} [GeV/c]; entries", 2000, 0, 2.0 );
  rapidity_true_mass_ = new TH2F( "y_tof_vs_y_pdg", ";PDG-mass y_{cm};TOF-mass y_{cm};", 200, -1.0, 1.0, 200, -1.0, 1.0 );
  centrality_classes_ = new TH1F( "centrality", ";centrality", 20, 0.0, 100.0 );
  out_file_->cd();

  std::cout << "Initialized" << std::endl;
}

void Yield::UserExec() {
  using AnalysisTree::Particle;

  auto centrality = (*event_header_)[centrality_].GetVal();
  centrality_classes_->Fill( centrality );
  int c_class = (int) ( (centrality-2.5)/5.0 );
  float y_beam_2 = data_header_->GetBeamRapidity();
  TH2F* histo{nullptr};
  if( centrality > 25 )
    return;
  if( centrality < 20 )
    return;
  try {
    histo = yields_.at(c_class);
  } catch (std::out_of_range&) { return; }
  for (auto track : tracks_->Loop()) {
    auto chi2 = track[chi2_].GetVal();
    auto dca_xy = track[dca_xy_].GetVal();
    auto dca_z = track[dca_z_].GetVal();
    if( chi2 > 100.0 )
      continue;
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    auto pid = track[pdg_code_];
    if( pid != 2212 )
      continue;
    const float p_mass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    auto mom4 = track.DataT<Particle>()->Get4MomentumByMass(p_mass);
    auto y_tof = track.DataT<Particle>()->GetRapidity() - y_beam_2;
    auto pT = mom4.Pt();
    auto y_pdg = mom4.Rapidity() - y_beam_2;
    histo->Fill(y_pdg, pT);
    rapidity_true_mass_->Fill( y_pdg, y_tof );
    if( y_pdg < 0.15 )
      continue;
    if( y_pdg > 0.25 )
      continue;
    pt_distribution_reco_->Fill(pT);
  }
  for (auto particle : sim_particles_->Loop()) {
    if( particle[sim_pdg_code_].GetInt() != 2212 )
      continue;
    if( !particle[is_primary_].GetBool() )
      continue;
    auto mass = particle.DataT<Particle>()->GetMass();
    auto mom4 = particle.DataT<Particle>()->Get4MomentumByMass(mass);
    auto y_cm = mom4.Rapidity() - y_beam_2;
    if( y_cm < 0.15 )
      continue;
    if( y_cm > 0.25 )
      continue;
    pt_distribution_sim_->Fill( mom4.Pt() );
  }
}
void Yield::UserFinish() {
  centrality_classes_->Write();
  pt_distribution_reco_->Write();
  pt_distribution_sim_->Write();
  rapidity_true_mass_->Write();
  for( auto histo : yields_ ) {
    histo->Sumw2();
    histo->Write();
  }
  std::cout << "Finished" << std::endl;
}
