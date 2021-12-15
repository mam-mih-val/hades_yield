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
    ("acceptance-protons", value(&str_protons_acceptance_)->default_value(""), "PDG-code of particle")
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

  h2_tru_pid_pT_theta_ = new TH2F( "h2_tru_pid_pT_theta_", "#theta;centrality (%)",
                                  200, 0.0, 2.0,
                                  170, 0.0, 1.7 );

  auto file = TFile::Open( str_protons_acceptance_.c_str(), "read" );
  if( file )
    file->GetObject( "h2_rec_2212_pT_theta_", h2_acceptacne_2212_pT_theta_ );

  h2_rec_all_nprart_theta_phi_in_event_ = new TH2F("h2_rec_all_nprart_theta_phi_in_event_",
                                                   ";theta (rad); phi (rad)",
                                                   1, 0.2, 1.6,
                                                   6, -M_PI, M_PI );

  h3_rec_pid_nprart_theta_phi_pT_in_event_ = new TH3F("h3_rec_pid_nprart_theta_phi_pT_in_event_",
                                                   ";theta (rad); phi (rad); p_{T} (GeV/c)",
                                                   14, 0.2, 1.6,
                                                   6, -M_PI, M_PI,
                                                   10, 0, 2.0 );
  h3_tru_pid_nprart_theta_phi_pT_in_event_ = new TH3F("h3_tru_pid_nprart_theta_phi_pT_in_event_",
                                                   ";theta (rad); phi (rad); p_{T} (GeV/c)",
                                                   14, 0.2, 1.6,
                                                   6, -M_PI, M_PI,
                                                   10, 0, 2.0 );

  p3_rec_pid_efficiency_theta_pT_track_density_ = new TProfile3D("p3_rec_pid_efficiency_theta_pT_track_density_",
                                                                 ";theta (rad); p_{T} (GeV/c); N_{Tracks}/S",
                                                                 14, 0.2, 1.6,
                                                                 10, 0, 2.0,
                                                                 40, 0.0, 40.0);

  auto y_cm = data_header_->GetBeamRapidity();
  beta_cm_ = tanh(y_cm);
  out_file_->cd();
  std::cout << "Initialized" << std::endl;
}

void Yield::UserExec() {
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  h1_centrality_->Fill( centrality );

  h2_rec_all_nprart_theta_phi_in_event_->Reset("ICESM");
  h3_rec_pid_nprart_theta_phi_pT_in_event_->Reset("ICESM");
  h3_tru_pid_nprart_theta_phi_pT_in_event_->Reset("ICESM");

  this->LoopRecTracks();
  this->LoopTruParticles();

  for( int theta_bin=1; theta_bin<h3_tru_pid_nprart_theta_phi_pT_in_event_->GetNbinsX(); theta_bin++ ){
    auto theta_lo = h2_rec_all_nprart_theta_phi_in_event_->GetXaxis()->GetBinLowEdge(1);
    auto theta_hi = h2_rec_all_nprart_theta_phi_in_event_->GetXaxis()->GetBinUpEdge(1);
    auto surface = ConeSideSquare( theta_hi, theta_lo ) / 6.0;
    auto theta = h3_tru_pid_nprart_theta_phi_pT_in_event_->GetXaxis()->GetBinCenter(theta_bin);
    for( int phi_bin=1; phi_bin<h3_tru_pid_nprart_theta_phi_pT_in_event_->GetNbinsY(); phi_bin++ ){
      auto phi = h3_tru_pid_nprart_theta_phi_pT_in_event_->GetYaxis()->GetBinLowEdge(phi_bin);
      auto n_tracks =  h2_rec_all_nprart_theta_phi_in_event_->GetBinContent(1, phi_bin);
      for( int pT_bin=1; pT_bin<h3_tru_pid_nprart_theta_phi_pT_in_event_->GetNbinsZ(); pT_bin++ ){
        auto pT = h3_rec_pid_nprart_theta_phi_pT_in_event_->GetZaxis()->GetBinLowEdge(pT_bin);
        auto n_rec = h3_rec_pid_nprart_theta_phi_pT_in_event_->GetBinContent( theta_bin, phi_bin, pT_bin );
        auto n_tru = h3_tru_pid_nprart_theta_phi_pT_in_event_->GetBinContent( theta_bin, phi_bin, pT_bin );
        if( n_tru < 1.0 )
          continue;
        auto track_density_unbiased = n_tracks - n_rec;
        p3_rec_pid_efficiency_theta_pT_track_density_->Fill( theta, pT, track_density_unbiased, n_rec/n_tru );
      }
    }
  }
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
    h2_rec_all_nprart_theta_phi_in_event_->Fill( mom4.Theta(), mom4.Phi() );
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    if( pid != reference_pdg_code_ )
      continue;
    p2_rec_v1_pid_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
    h2_acceptacne_2212_pT_theta_->Fill( mom4.Pt(), mom4.Theta() );
    h3_rec_pid_nprart_theta_phi_pT_in_event_->Fill( mom4.Theta(), mom4.Phi(), mom4.Pt() );
    if( pid == 2212 && mom4.Pt() < 0.4 )
      continue;
    h3_rec_delta_phi_theta_centrality_pid_->Fill(delta_phi, mom4.Theta(), centrality);
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
    double charge=0.0;
    if( TDatabasePDG::Instance()->GetParticle( pid ) ){
      charge= TDatabasePDG::Instance()->GetParticle( pid )->Charge() / 3.0;
    }
    auto delta_phi = AngleDifference(mom4.Phi(), psi_rp);
    if( fabs(charge) < 0.01 )
      continue;
    if( pid == 2212 ){
      if( h2_acceptacne_2212_pT_theta_ ) {
        auto pT_bin = h2_acceptacne_2212_pT_theta_->GetXaxis()->FindBin(mom4.Pt());
        auto theta_bin = h2_acceptacne_2212_pT_theta_->GetYaxis()->FindBin(mom4.Theta());
        auto n_entries = h2_acceptacne_2212_pT_theta_->GetBinContent(pT_bin, theta_bin);
        if( n_entries < 1.0 )
          continue;
      }
    }
    h3_tru_delta_phi_theta_centrality_all_->Fill(delta_phi, mom4.Theta(), centrality);
    if( !is_prim )
      continue;
    if( pid!=reference_pdg_code_ )
      continue;
    p2_tru_v1_all_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
    h2_tru_theta_centrality_all_->Fill(mom4.Theta(), centrality);
    p2_tru_v1_pid_->Fill( mom4.Theta(), centrality, cos(delta_phi) );
    h2_tru_pid_pT_theta_->Fill( mom4.Pt(), mom4.Theta() );
    h3_tru_pid_nprart_theta_phi_pT_in_event_->Fill( mom4.Theta(), mom4.Phi(), mom4.Pt() );
    if( pid == 2212 && mom4.Pt() < 0.4 )
      continue;
    h3_tru_delta_phi_theta_centrality_pid_->Fill(delta_phi, mom4.Theta(), centrality);
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
  h2_tru_pid_pT_theta_->Write();
  p3_rec_pid_efficiency_theta_pT_track_density_->Write();
  h2_rec_all_nprart_theta_phi_in_event_->Write();
  h3_rec_pid_nprart_theta_phi_pT_in_event_->Write();
  h3_tru_pid_nprart_theta_phi_pT_in_event_->Write();
  std::cout << "Finished" << std::endl;
}
