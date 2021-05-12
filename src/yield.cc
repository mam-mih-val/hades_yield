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
//  SetInputBranchNames({"mdc_vtx_tracks", "event_header", "sim_tracks"});
}

void Yield::UserInit(std::map<std::string, void *> &Map) {
  // extracting event variables
  event_header_ = GetInBranch("event_header");
  centrality_var_ = GetVar("event_header/selected_tof_rpc_hits_centrality");
  // extracting tracks variables
  tracks_ = GetInBranch("mdc_vtx_tracks");
  charge_var_ = GetVar("mdc_vtx_tracks/charge");
  pdg_code_var_ = GetVar("mdc_vtx_tracks/pid");
  dca_xy_var_ = GetVar("mdc_vtx_tracks/dca_xy");
  dca_z_var_ = GetVar("mdc_vtx_tracks/dca_z");
  chi2_var_ = GetVar("mdc_vtx_tracks/chi2");
  dedx_mdc_var_ = GetVar("mdc_vtx_tracks/dEdx");
  mdc2meta_matching_ = static_cast<AnalysisTree::Matching*>(Map.at("mdc_vtx_tracks2meta_hits"));
  // extracting meta hits variables
  meta_hits_ = GetInBranch("meta_hits");
  mass2_var_ = GetVar("meta_hits/mass2");
  beta_var_ = GetVar("meta_hits/beta");
  beta_var_ = GetVar("meta_hits/beta");
  dedx_meta_var_ = GetVar("meta_hits/dEdx");
  is_rpc_hit_var_ = GetVar("meta_hits/is_rpc_hit");

  try {
    sim_particles_ = GetInBranch("sim_tracks");
    is_mc_=true;
  } catch (std::exception&) {
    is_mc_=false;
  }
  if( is_mc_ ) {
    sim_particles_ = GetInBranch("sim_tracks");
    sim_pdg_code_var_ = GetVar("sim_tracks/pid");
    is_primary_var_ = GetVar("sim_tracks/is_primary");
    mdc2sim_matching_ = static_cast<AnalysisTree::Matching*>(Map.at("mdc_vtx_tracks2sim_tracks"));
  }
  centrality_classes_ = new TH1F( "centrality", ";centrality", 20, 0.0, 100.0 );
  int p=0;
  while(p<60){
    // All
    std::string histo_name{ "m2_vs_pq_all_"+std::to_string(p)+"-"+std::to_string(p+5) };
    m2_vs_pq_all_.push_back( new TH2F( histo_name.c_str(), ";p/q [GeV/c]; m^{2} [GeV^{2}/c^{4}]", 300, -2.50, 5.0, 420, -0.5, 10 ) );
    histo_name =  "beta_vs_pq_all_"+std::to_string(p)+"-"+std::to_string(p+5);
    beta_vs_pq_all_.push_back( new TH2F( histo_name.c_str(), ";p/q [GeV/c]; #beta [1/c]", 300, -2.50, 5.0, 550, 0.0, 1.1 ) );
    histo_name =  "dedx_mdc_vs_pq_all_"+std::to_string(p)+"-"+std::to_string(p+5) ;
    dedx_mdc_vs_pq_all_.push_back( new TH2F( histo_name.c_str(), ";p/q [GeV/c]; dE/dx MDC", 300, -2.50, 5.0, 400, 0.0, 20.0 ) );
    histo_name =  "dedx_meta_vs_pq_all_"+std::to_string(p)+"-"+std::to_string(p+5) ;
    dedx_meta_vs_pq_all_.push_back( new TH2F( histo_name.c_str(), ";p/q [GeV/c]; dE/dx META", 300, -2.50, 5.0, 400, 0.0, 20.0 ) );
    // PID-Reco
    histo_name = "m2_vs_pq_pid_reco_"+std::to_string(p)+"-"+std::to_string(p+5);
    m2_vs_pq_pid_reco_.push_back( new TH2F( histo_name.c_str(), ";p/q [GeV/c]; m^{2} [GeV^{2}/c^{4}]", 300, -2.50, 5.0, 420, -0.5, 10 ) );
    histo_name =  "beta_vs_pq_pid_reco_"+std::to_string(p)+"-"+std::to_string(p+5);
    beta_vs_pq_pid_reco_.push_back( new TH2F( histo_name.c_str(), ";p/q [GeV/c]; #beta [1/c]", 300, -2.50, 5.0, 550, 0.0, 1.1 ) );
    histo_name =  "dedx_mdc_vs_pq_pid_reco_"+std::to_string(p)+"-"+std::to_string(p+5) ;
    dedx_mdc_vs_pq_pid_reco_.push_back( new TH2F( histo_name.c_str(), ";p/q [GeV/c]; dE/dx MDC", 300, -2.50, 5.0, 400, 0.0, 20.0 ) );
    histo_name =  "dedx_meta_vs_pq_pid_reco_"+std::to_string(p)+"-"+std::to_string(p+5) ;
    dedx_meta_vs_pq_pid_reco_.push_back( new TH2F( histo_name.c_str(), ";p/q [GeV/c]; dE/dx META", 300, -2.50, 5.0, 400, 0.0, 20.0 ) );
    // Mismatch
    if( is_mc_ ) {
      histo_name = "m2_vs_pq_mismatch_" + std::to_string(p) + "-" + std::to_string(p + 5);
      m2_vs_pq_mismatch_.push_back(new TH2F(histo_name.c_str(), ";p/q [GeV/c]; m^{2} [GeV^{2}/c^{4}]",300, -2.50, 5.0, 420, -0.5, 10));
      histo_name = "beta_vs_pq_mismatch_" + std::to_string(p) + "-" + std::to_string(p + 5);
      beta_vs_pq_mismatch_.push_back(new TH2F(histo_name.c_str(),";p/q [GeV/c]; #beta [1/c]", 300,-2.50, 5.0, 550, 0.0, 1.1));
      histo_name = "dedx_mdc_vs_pq_mismatch_" + std::to_string(p) + "-" +std::to_string(p + 5);
      dedx_mdc_vs_pq_mismatch_.push_back(new TH2F(histo_name.c_str(), ";p/q [GeV/c]; dE/dx MDC", 300, -2.50,5.0, 400, 0.0, 20.0));
      histo_name = "dedx_meta_vs_pq_mismatch_" + std::to_string(p) + "-" +std::to_string(p + 5);
      dedx_meta_vs_pq_mismatch_.push_back(new TH2F(histo_name.c_str(), ";p/q [GeV/c]; dE/dx META", 300, -2.50,5.0, 400, 0.0, 20.0));
    }
    p+=5;
  }
  pt_distribution_reco_ = new TH1F( "pT_distribution_reco", ";p_{T} [GeV/c]; entries", 2000, 0, 2.0 );
  pt_distribution_sim_ = new TH1F( "pT_distribution_sim", ";p_{T} [GeV/c]; entries", 2000, 0, 2.0 );
  rapidity_true_mass_ = new TH2F( "y_tof_vs_y_pdg", ";PDG-mass y_{cm};TOF-mass y_{cm};", 200, -1.0, 1.0, 200, -1.0, 1.0 );
  out_file_->cd();
  std::cout << "Initialized" << std::endl;
}

void Yield::UserExec() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[centrality_var_].GetVal();
  centrality_classes_->Fill( centrality );
  int c_class = (size_t) ( (centrality-2.5)/5.0 );
  if ( centrality < 0 )
    return;
  if ( centrality > 60 )
    return;
  float y_beam = data_header_->GetBeamRapidity();
  int track_idx=-1;
  for (auto track : tracks_->Loop()) {
    track_idx++;
    auto chi2 = track[chi2_var_].GetVal();
    auto dca_xy = track[dca_xy_var_].GetVal();
    auto dca_z = track[dca_z_var_].GetVal();
    if( chi2 > 100.0 )
      continue;
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    auto meta_idx = mdc2meta_matching_->GetMatchDirect(track_idx);
    auto meta_hit = (*meta_hits_)[meta_idx];
    auto pid = track[pdg_code_var_].GetInt();
    auto charge = track[charge_var_].GetInt();
    auto dEdx_mdc = track[dedx_mdc_var_].GetVal();
    auto mass2 = meta_hit[mass2_var_].GetVal();
    auto beta = meta_hit[beta_var_].GetVal();
    auto dEdx_meta = meta_hit[dedx_meta_var_].GetVal();
    auto mass = track.DataT<Particle>()->GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    }
    auto mom4 = track.DataT<Particle>()->Get4MomentumByMass(mass);
    auto p = mom4.P();
    auto y = mom4.Rapidity() - y_beam;
    if( mom4.Pt() > 0.6 )
      continue;
    if ( fabs(y) > 0.25 )
      continue;
    // All particles
    m2_vs_pq_all_.at(c_class)->Fill( p/charge, mass2 );
    beta_vs_pq_all_.at(c_class)->Fill( p/charge, beta );
    dedx_mdc_vs_pq_all_.at(c_class)->Fill( p/charge, dEdx_mdc );
    dedx_meta_vs_pq_all_.at(c_class)->Fill( p/charge, dEdx_meta );
    // PID-reco
    if( pid == 2212 ){
      m2_vs_pq_pid_reco_.at(c_class)->Fill( p/charge, mass2 );
      beta_vs_pq_pid_reco_.at(c_class)->Fill( p/charge, beta );
      dedx_mdc_vs_pq_pid_reco_.at(c_class)->Fill( p/charge, dEdx_mdc );
      dedx_meta_vs_pq_pid_reco_.at(c_class)->Fill( p/charge, dEdx_meta );
    }
    if( is_mc_ ){
      auto gen_idx = mdc2sim_matching_->GetMatchDirect(track_idx);
      if( gen_idx ==  AnalysisTree::UndefValueInt )
        continue;
      auto gen_particle = (*sim_particles_)[gen_idx];
      auto gen_pid = gen_particle[sim_pdg_code_var_].GetInt();
      if( pid == 2212 && !gen_particle[is_primary_var_].GetBool() ) {
        m2_vs_pq_mismatch_.at(c_class)->Fill(p / charge, mass2);
        beta_vs_pq_mismatch_.at(c_class)->Fill(p / charge, beta);
        dedx_mdc_vs_pq_mismatch_.at(c_class)->Fill(p / charge, dEdx_mdc);
        dedx_meta_vs_pq_mismatch_.at(c_class)->Fill(p / charge, dEdx_meta);
      }
    }
  }
}
void Yield::UserFinish() {
  centrality_classes_->Write();
  pt_distribution_reco_->Write();
  pt_distribution_sim_->Write();
  rapidity_true_mass_->Write();
  for( size_t i=0; i<m2_vs_pq_all_.size(); ++i ){
    m2_vs_pq_all_.at(i)->Write();
    beta_vs_pq_all_.at(i)->Write();
    dedx_mdc_vs_pq_all_.at(i)->Write();
    dedx_meta_vs_pq_all_.at(i)->Write();

    m2_vs_pq_pid_reco_.at(i)->Write();
    beta_vs_pq_pid_reco_.at(i)->Write();
    dedx_mdc_vs_pq_pid_reco_.at(i)->Write();
    dedx_meta_vs_pq_pid_reco_.at(i)->Write();

    if( is_mc_ ){
      m2_vs_pq_mismatch_.at(i)->Write();
      beta_vs_pq_mismatch_.at(i)->Write();
      dedx_mdc_vs_pq_mismatch_.at(i)->Write();
      dedx_meta_vs_pq_mismatch_.at(i)->Write();
    }
  }
  std::cout << "Finished" << std::endl;
}
