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

  h3_tru_delta_phi_theta_centrality_all_ = new TH3F("h3_tru_delta_phi_theta_centrality_all",
                                                    ";#Delta#phi (rad);#theta;centrality (%)",
                                                    350, -3.5, 3.5,
                                                    140, 0.2, 1.6,
                                                    12, 0.0, 60.);
  h3_tru_delta_phi_theta_centrality_pid_ = new TH3F("h3_tru_delta_phi_theta_centrality_pid",
                                                    ";#Delta#phi (rad);#theta;centrality (%)",
                                                    350, -3.5, 3.5,
                                                    140, 0.2, 1.6,
                                                    12, 0.0, 60.);
  h3_rec_delta_phi_theta_centrality_all_ = new TH3F("h3_rec_delta_phi_theta_centrality_all",
                                                    ";#Delta#phi (rad);#theta;centrality (%)",
                                                    350, -3.5, 3.5,
                                                    140, 0.2, 1.6,
                                                    12, 0.0, 60.);
  h3_rec_delta_phi_theta_centrality_pid_ = new TH3F("h3_rec_delta_phi_theta_centrality_pid",
                                                    ";#Delta#phi (rad);#theta;centrality (%)",
                                                    350, -3.5, 3.5,
                                                    140, 0.2, 1.6,
                                                    12, 0.0, 60.);

  h2_rec_theta_centrality_all_ = new TH2F( "h2_rec_theta_centrality_all",
                                          "#theta;centrality (%)",
                                          140, 0.2, 1.6,
                                          12, 0.0, 60.);
  h2_tru_theta_centrality_all_ = new TH2F( "h2_tru_theta_centrality_all",
                                          "#theta;centrality (%)",
                                          140, 0.2, 1.6,
                                          12, 0.0, 60.);

  p2_tru_v1_pid_ = new TProfile2D( "p2_tru_v1_pid", ";theta;centrality", 140, 0.2, 1.6, 12, 0.0, 60 );
  p2_rec_v1_pid_ = new TProfile2D( "p2_rec_v1_pid", ";theta;centrality", 140, 0.2, 1.6, 12, 0.0, 60 );

  p2_tru_v1_all_ = new TProfile2D( "p2_tru_v1_all", ";theta;centrality", 140, 0.2, 1.6, 12, 0.0, 60 );
  p2_rec_v1_all_ = new TProfile2D( "p2_rec_v1_all", ";theta;centrality", 140, 0.2, 1.6, 12, 0.0, 60 );

  p3_dtheta_dphi_dpT_loss_ = new TProfile3D("p3_dtheta_dphi_dpT_loss_",
                                            ";#phi_{1}-#phi_{2} (rad);#theta_{1}-#theta_{2} (rad);p_{T,1}-p_{T,2} (GeV/c)",
                                            150, -1.5, 1.5,
                                            150, -1.5, 1.5,
                                            150, -1.5, 1.5  );
  h2_theta_phi_sector_lost_population_ = new TH2F( "h2_theta_phi_sector_lost_population",
                                                  ";#phi_{1}-#phi_{2} (rad);#theta_{1}-#theta_{2} (rad)",
                                                  50, -0.75, 0.75,
                                                  120, 0.3, 1.5);
  p2_dphi_dtheta_second_efficiency_ = new TProfile2D( "p2_dphi_dtheta_second_efficiency_",
                             ";#phi_{1}-#phi_{2} (rad);#theta_{1}-#theta_{2} (rad)",
                             150, -1.5, 1.5,
                             150, -1.5, 1.5);
  p2_dphi_dtheta_pair_efficiency_ = new TProfile2D( "p2_dphi_dtheta_pair_efficiency_",
                             ";#phi_{1}-#phi_{2} (rad);#theta_{1}-#theta_{2} (rad)",
                             150, -1.5, 1.5,
                             150, -1.5, 1.5);

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

  for (auto track : tracks_->Loop()) {
    auto pid = track.DataT<Particle>()->GetPid();
    auto mass = track.DataT<Particle>()->GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    }
    auto mom4 = track.DataT<Particle>()->Get4MomentumByMass(mass);
    auto chi2 = track[rec_chi2_var].GetVal();
    auto dca_xy = track[rec_dca_xy_var].GetVal();
    auto dca_z = track[rec_dca_z_var].GetVal();
    auto delta_phi = AngleDifference(mom4.Phi(), psi_rp);
    h3_rec_delta_phi_theta_centrality_all_->Fill(delta_phi, mom4.Theta(), centrality);
    p2_rec_v1_all_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
    h2_rec_theta_centrality_all_->Fill(mom4.Theta(), centrality);
    if( chi2 > 100.0 )
      continue;
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    if( pid != reference_pdg_code_ )
      continue;
    h3_rec_delta_phi_theta_centrality_pid_->Fill(delta_phi, mom4.Theta(), centrality);
    p2_rec_v1_pid_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
  }
}

void Yield::LoopTruParticles() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  auto psi_rp = (*sim_header_)[GetVar( "sim_header/reaction_plane" )].GetVal();

  auto var_is_primary = GetVar("sim_tracks/is_primary");

  int idx1 =-1;
  for( auto particle : sim_particles_->Loop() ){
    idx1++;
    auto mass = particle.DataT<Particle>()->GetMass();
    auto pid = particle.DataT<Particle>()->GetPid();
    auto is_prim = particle[var_is_primary].GetBool();
    auto mom4 = particle.DataT<Particle>()->Get4MomentumByMass(mass);
    auto phi = mom4.Phi() + M_PI;
    int sector1 = phi / ( M_PI/3.0 );
    auto match1 = rec_sim_matching_->GetMatchInverted(idx1);
    double charge=0.0;
    if( TDatabasePDG::Instance()->GetParticle( pid ) ){
      charge= TDatabasePDG::Instance()->GetParticle( pid )->Charge() / 3.0;
    }
    auto delta_phi = AngleDifference(mom4.Phi(), psi_rp);
    if( fabs(charge) < 0.01 )
      continue;
    if(  mom4.Theta() < 0.3 )
      continue;
    if( mom4.Theta() > 1.5 )
      continue;
//    if( mom4.Pt() < 0.4 )
//      continue;
//    if( !is_prim )
//      continue;
    for( int idx2= 0; idx2 < sim_particles_->size(); idx2++ ){
      if( idx2 == idx1 )
        continue;
      auto particle2 = (*sim_particles_)[idx2];
      auto mom42 = particle2.DataT<Particle>()->Get4MomentumByMass(mass);
      auto phi2 = mom42.Phi() + M_PI;
      int sector2 = phi2 / ( M_PI/3.0 );
      if( sector1 != sector2 )
        continue;
      double charge2=0.0;
      auto pid2 = particle2.DataT<Particle>()->GetPid();
      if( TDatabasePDG::Instance()->GetParticle( pid2 ) ){
        charge2= TDatabasePDG::Instance()->GetParticle( pid2 )->Charge() / 3.0;
      }
      auto is_prim2 = particle2[var_is_primary].GetBool();
      if( fabs(charge2) < 0.01 )
        continue;
      if( mom42.Theta() < 0.3 )
        continue;
      if( mom42.Theta() > 1.5 )
        continue;
//      if( !is_prim2 )
//        continue;
//      if( mom42.Pt() < 0.4 )
//        continue;
      auto dphi = AngleDifference( mom4.Phi(), mom42.Phi() );
      auto dtheta = AngleDifference( mom4.Theta(), mom42.Theta() );
      auto dpT =  mom4.Pt() - mom42.Pt();
      auto match2  = rec_sim_matching_->GetMatchInverted(idx2);
      int n_rec = 0;
      if( match1 >= 0 )
        n_rec++;
      if( match2 >= 0 )
        n_rec++;
      if( match1 == match2 && match1 >= 0 )
        n_rec = 1;

      p3_dtheta_dphi_dpT_loss_->Fill( dphi, dtheta, dpT, 2-n_rec );
      if( n_rec < 2 ){
        double sector_phi = ( sector1*M_PI/3.0 + (sector1+1)*M_PI/3.0 ) / 2 - M_PI;
        h2_theta_phi_sector_lost_population_ -> Fill( mom4.Phi() - sector_phi, mom4.Theta(), 2-n_rec );
        h2_theta_phi_sector_lost_population_ -> Fill( mom42.Phi() - sector_phi, mom42.Theta(), 2-n_rec );
      }
      if( match2 < 0 )
        p2_dphi_dtheta_second_efficiency_->Fill( dphi, dtheta, 0 );
      if( match2 >= 0 )
        p2_dphi_dtheta_second_efficiency_->Fill( dphi, dtheta, 1 );

      int first_and_second = match1 >= 0 && match2 >= 0;
      p2_dphi_dtheta_pair_efficiency_->Fill(dphi, dtheta, first_and_second);
    }
    h3_tru_delta_phi_theta_centrality_all_->Fill(delta_phi, mom4.Theta(), centrality);
    p2_tru_v1_all_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
    h2_tru_theta_centrality_all_->Fill(mom4.Theta(), centrality);
    if( pid!=reference_pdg_code_ )
      continue;
    h3_tru_delta_phi_theta_centrality_pid_->Fill(delta_phi, mom4.Theta(), centrality);
    p2_tru_v1_pid_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
  }
}


void Yield::UserFinish() {
  out_file_->cd();
  h1_centrality_->Write();

//  out_file_->mkdir("efficiency_projections");
//  out_file_->cd("efficiency_projections");

  out_file_->cd();
  h3_tru_delta_phi_theta_centrality_all_->Write();
  h3_tru_delta_phi_theta_centrality_pid_->Write();
  h3_rec_delta_phi_theta_centrality_all_->Write();
  h3_rec_delta_phi_theta_centrality_pid_->Write();
  h2_rec_theta_centrality_all_->Write();
  h2_tru_theta_centrality_all_->Write();
  p2_tru_v1_pid_->Write();
  p2_rec_v1_pid_->Write();
  p2_tru_v1_all_->Write();
  p2_rec_v1_all_->Write();
  p3_dtheta_dphi_dpT_loss_->Write();
  h2_theta_phi_sector_lost_population_->Write();
  p2_dphi_dtheta_second_efficiency_->Write();
  p2_dphi_dtheta_pair_efficiency_->Write();
  std::cout << "Finished" << std::endl;
}
