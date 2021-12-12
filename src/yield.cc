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
    ("pdg-code", value(&reference_pdg_code_)->default_value(2212), "PDG-code of particle");
  return desc;
}

void Yield::PreInit() {
//  SetInputBranchNames({"mdc_vtx_tracks", "event_header", "sim_tracks"});
}

void Yield::UserInit(std::map<std::string, void *> &Map) {
  // extracting event variables
  event_header_ = GetInBranch("event_header");
  sim_header_ = GetInBranch("sim_header");
  tracks_ = GetInBranch("mdc_vtx_tracks");
  sim_particles_ = GetInBranch("sim_tracks");
  rec_sim_matching_ = static_cast<AnalysisTree::Matching*>(Map.at("mdc_vtx_tracks2sim_tracks"));

  h1_centrality_ = new TH1F( "centrality", ";TOF+RPC hits centrality (%)", 20, 0.0, 100.0 );

  int nbins[] = {12, 10, 10, 12, 10, 12};
  double xmin[] = {0., 0.0, -0.5, 0.3, -0.5, 0.3};
  double xmax[] = {60., 2.0, +0.5, 1.5, 0.5, 1.5};
  hn_tru_npairs_centrality_phi_theta_ = new THnSparseF( "hn_tru_npairs_centrality_phi_theta_",
                                                 ";centrality (%);p_{T} (GeV/c);#phi_{1} (rad);#theta_{1} (rad);#phi_{2} (rad);#theta_{2} (rad);",
                                                 6, nbins, xmin, xmax);
  hn_rec_npairs_centrality_phi_theta_ = new THnSparseF( "hn_rec_npairs_centrality_phi_theta_",
                                                 ";centrality (%);p_{T} (GeV/c);#phi_{1} (rad);#theta_{1} (rad);#phi_{2} (rad);#theta_{2} (rad);",
                                                 6, nbins, xmin, xmax);

  h3_rec_all_npart_centrality_phi_theta_ = new TH3F( "h3_rec_all_npart_centrality_phi_theta_",
                                            ";centrality (%);p_{T} (GeV/c);#theta (rad)",
                                            12, 0, 60,
                                            100, 0, 2,
                                            120, 0.3, 1.5);
  h3_tru_all_npart_centrality_phi_theta_ = new TH3F( "h3_tru_all_npart_centrality_phi_theta_",
                                            ";centrality (%);p_{T} (GeV/c);#theta (rad)",
                                            12, 0, 60,
                                            100, 0, 2,
                                            120, 0.3, 1.5);

  h3_rec_pid_npart_centrality_phi_theta_ = new TH3F( "h3_rec_pid_npart_centrality_phi_theta_",
                                            ";centrality (%);p_{T} (GeV/c);#;#theta (rad)",
                                            12, 0, 60,
                                            100, 0, 2,
                                            120, 0.3, 1.5);
  h3_tru_pid_npart_centrality_phi_theta_ = new TH3F( "h3_tru_pid_npart_centrality_phi_theta_",
                                            ";centrality (%);p_{T} (GeV/c);#;#theta (rad)",
                                            12, 0, 60,
                                            100, 0, 2,
                                            120, 0.3, 1.5);

  auto y_cm = data_header_->GetBeamRapidity();
  beta_cm_ = tanh(y_cm);
  out_file_->cd();
  std::cout << "Initialized" << std::endl;
}

void Yield::UserExec() {
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  h1_centrality_->Fill( centrality );

  this->LoopRecTracks();
  this->LoopTruParticles();
}

void Yield::LoopRecTracks() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  auto psi_rp = (*sim_header_)[GetVar( "sim_header/reaction_plane" )].GetVal();

  auto rec_chi2_var = GetVar("mdc_vtx_tracks/chi2");
  auto rec_dca_xy_var = GetVar("mdc_vtx_tracks/dca_xy");
  auto rec_dca_z_var = GetVar("mdc_vtx_tracks/dca_z");

  for (size_t idx1=0; idx1 < tracks_->size(); idx1++) {
    auto track1 = (*tracks_)[idx1];
    auto pid1 = track1.DataT<Particle>()->GetPid();
    auto mass1 = track1.DataT<Particle>()->GetMass();
    if(pid1 != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid1) )
        mass1 = TDatabasePDG::Instance()->GetParticle(pid1)->Mass();
    }
    auto mom1 = track1.DataT<Particle>()->Get4MomentumByMass(mass1);
    auto chi_sq1 = track1[rec_chi2_var].GetVal();
    auto dca_xy = track1[rec_dca_xy_var].GetVal();
    auto dca_z = track1[rec_dca_z_var].GetVal();
    if(chi_sq1 > 100.0 )
      continue;
    double phi = mom1.Phi();
    double sector1 = floor(phi / (M_PI/3.0));
    double sector_center = sector1 *(M_PI/3.0) + M_PI/6.0;
    double sector_phi1 = AngleDifference(mom1.Phi(),sector_center);
    h3_rec_all_npart_centrality_phi_theta_->Fill( centrality, mom1.Pt(),mom1.Theta() );
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    if(pid1 != reference_pdg_code_ )
      continue;
    h3_rec_pid_npart_centrality_phi_theta_->Fill( centrality, mom1.Pt(),mom1.Theta() );

    // Second Loop through particles to fill histogram for pairs
    for( size_t idx2=0; idx2 < tracks_->size(); idx2++ ){
      if( idx1 == idx2 )
        continue;
      auto track2 = (*tracks_)[idx2];
      auto pid2 = track2.DataT<Particle>()->GetPid();
      auto mass2 = track2.DataT<Particle>()->GetMass();
      if(pid2 != 0 ) {
        if( TDatabasePDG::Instance()->GetParticle(pid2) )
          mass2 = TDatabasePDG::Instance()->GetParticle(pid2)->Mass();
      }
      auto mom2 = track2.DataT<Particle>()->Get4MomentumByMass(mass2);
      auto chi_sq2 = track2[rec_chi2_var].GetVal();
      if(chi_sq2 > 100.0 )
        continue;
      double sector2 = floor(mom2.Phi() / (M_PI/3.0));
      if( fabs(sector1-sector2) > 0.001 )
        continue;
      auto sector_phi2 = AngleDifference(mom2.Phi(),sector_center);
      double args[] = {centrality, mom1.Pt(), sector_phi1, mom1.Theta(), sector_phi2, mom2.Theta()};
      hn_rec_npairs_centrality_phi_theta_->Fill(args);
    }
  }
}

void Yield::LoopTruParticles() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  auto psi_rp = (*sim_header_)[GetVar( "sim_header/reaction_plane" )].GetVal();

  auto var_is_primary = GetVar("sim_tracks/is_primary");

  for( size_t idx1=0; idx1< sim_particles_->size(); idx1++ ){
    auto particle1 = (*sim_particles_)[idx1];
    auto mass1 = particle1.DataT<Particle>()->GetMass();
    auto pid1 = particle1.DataT<Particle>()->GetPid();
    auto is_prim = particle1[var_is_primary].GetBool();
    auto mom1 = particle1.DataT<Particle>()->Get4MomentumByMass(mass1);
    double phi = mom1.Phi();
    double sector1 = floor(phi / (M_PI/3.0));
    double sector_center = sector1 *(M_PI/3.0) + M_PI/6.0;
    double sector_phi1 = AngleDifference(mom1.Phi(), sector_center);

    double charge=0.0;
    if( TDatabasePDG::Instance()->GetParticle(pid1) ){
      charge= TDatabasePDG::Instance()->GetParticle(pid1)->Charge() / 3.0;
    }
    if( fabs(charge) < 0.01 )
      continue;
    h3_tru_all_npart_centrality_phi_theta_->Fill( centrality, mom1.Pt(), mom1.Theta() );
    if( !is_prim )
      continue;
    if(pid1 !=reference_pdg_code_ )
      continue;
    h3_tru_pid_npart_centrality_phi_theta_->Fill( centrality, mom1.Pt(), mom1.Theta() );
    for( int idx2= 0; idx2 < sim_particles_->size(); idx2++ ){
      if( idx2 == idx1 )
        continue;
      auto particle2 = (*sim_particles_)[idx2];
      auto mass2 = particle2.DataT<Particle>()->GetMass();
      auto mom2 = particle2.DataT<Particle>()->Get4MomentumByMass(mass2);
      double sector2 = floor(mom2.Phi() / ( M_PI/3.0 ));
      if( fabs( sector1-sector2 ) > 0.0001 )
        continue;
      double sector_phi2 = AngleDifference(mom2.Phi(), sector_center);
      double charge2=0.0;
      auto pid2 = particle2.DataT<Particle>()->GetPid();
      if( TDatabasePDG::Instance()->GetParticle( pid2 ) ){
        charge2= TDatabasePDG::Instance()->GetParticle( pid2 )->Charge() / 3.0;
      }
      if( fabs(charge2) < 0.01 )
        continue;
      double args[] = {centrality, mom1.Pt(), sector_phi1, mom1.Theta(), sector_phi2, mom2.Theta()};
      hn_tru_npairs_centrality_phi_theta_->Fill(args);
    }
  }
}


void Yield::UserFinish() {
  out_file_->cd();
  h1_centrality_->Write();
  out_file_->cd();
  h3_rec_all_npart_centrality_phi_theta_->Write();
  h3_tru_all_npart_centrality_phi_theta_->Write();
  h3_rec_pid_npart_centrality_phi_theta_->Write();
  h3_tru_pid_npart_centrality_phi_theta_->Write();
  hn_rec_npairs_centrality_phi_theta_->Write();
  hn_tru_npairs_centrality_phi_theta_->Write();
  std::cout << "Finished" << std::endl;
}
