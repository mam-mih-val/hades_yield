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
  this->InitEfficiency();
  event_header_ = GetInBranch("event_header");
  sim_header_ = GetInBranch("sim_header");
  tracks_ = GetInBranch("mdc_vtx_tracks");
  sim_particles_ = GetInBranch("sim_tracks");

  h1_centrality_ = new TH1F( "centrality", ";TOF+RPC hits centrality (%)", 20, 0.0, 100.0 );
  h1_phi_theta_all_ = new TH2F( "phi_theta_event_by_event", ";#phi;N particles in event", 18, -M_PI, M_PI, 6, 0.3, 1.5 );


  std::vector<double> pt_axis;
  std::vector<double> theta_axis;
  std::vector<double> phi_axis;
  std::vector<double> npart_sector_axis;
  std::vector<double> centrality_axis;

  for( int i=0; i<60; i++ ){npart_sector_axis.push_back(i);}
  for( int i=0; i<37; i++ ){phi_axis.push_back(-M_PI+(i*2*M_PI)/36.0);}
  for( int i=0; i<25; i+=1 ){theta_axis.push_back(0.3+i*0.05); }
  for( int i=0; i<13; i+=1 ){centrality_axis.push_back(i*5.); }

  if( abs(reference_pdg_code_) == 2212 ) {
    pt_axis = {0.0,     0.29375, 0.35625, 0.41875, 0.48125, 0.54375,
               0.61875, 0.70625, 0.81875, 1.01875, 2.0};

  }
  if( abs(reference_pdg_code_) == 211 ) {
    pt_axis = {0,    0.08, 0.105, 0.13,  0.155, 0.18,
               0.21, 0.25, 0.315, 0.535, 1.0};
  }
  h3_rec_theta_pT_phi_ = new TH3F( "h3_rec_theta_pT_phi", ";#theta [rad];p_{T} [GeV/c];#phi [rad]",
                                  theta_axis.size()-1, theta_axis.data(),
                                  pt_axis.size()-1, pt_axis.data(),
                                  phi_axis.size()-1, phi_axis.data());
  h3_tru_theta_pT_phi_ = new TH3F( "h3_tru_theta_pT_phi", ";#theta [rad];p_{T} [GeV/c];#phi [rad]",
                                  theta_axis.size()-1, theta_axis.data(),
                                  pt_axis.size()-1, pt_axis.data(),
                                  phi_axis.size()-1, phi_axis.data());
  p3_theta_pT_npart_sector_ = new TProfile3D( "p3_theta_pT_npart_sector", ";#theta [rad];p_{T} [GeV/c];N particles in sector",
                                             theta_axis.size()-1, theta_axis.data(),
                                             pt_axis.size()-1, pt_axis.data(),
                                             npart_sector_axis.size()-1, npart_sector_axis.data());

  theta_axis.clear();
  size_t n_theta_bins=6;
  for( int i=0; i<n_theta_bins+1; i+=1 ){theta_axis.push_back(0.3+i*1.2/ (double) n_theta_bins); }
  h1_theta_all_ = new TH1F( "h1_theta_all", ";#theta [rad];N particles in event", theta_axis.size()-1, theta_axis.data() );
  h2_rec_pT_theta_ = new TH2F("h2_rec_pT_theta", "p_{T} [GeV/c]; #theta [rad]",
                                  pt_axis.size()-1, pt_axis.data(),
                                  theta_axis.size()-1, theta_axis.data());

  h2_tru_pT_theta_ = new TH2F("h2_tru_pT_theta", ";p_{T} [GeV/c]; #theta [rad]",
                                  pt_axis.size()-1, pt_axis.data(),
                                  theta_axis.size()-1, theta_axis.data());
  p3_pT_theta_npart_theta_ = new TProfile3D( "p3_pT_theta_npart_theta", ";p_{T} [GeV/c]; #theta [rad];N particles theta",
                                            pt_axis.size()-1, pt_axis.data(),
                                            theta_axis.size()-1, theta_axis.data(),
                                            npart_sector_axis.size()-1, npart_sector_axis.data());
  h3_rec_pT_theta_npart_theta_ = new TH3F( "h3_rec_pT_theta_npart_theta", ";p_{T} [GeV/c]; #theta [rad];N particles theta",
                                            pt_axis.size()-1, pt_axis.data(),
                                            theta_axis.size()-1, theta_axis.data(),
                                            npart_sector_axis.size()-1, npart_sector_axis.data());
  h3_tru_pT_theta_npart_theta_ = new TH3F( "h3_tru_pT_theta_npart_theta", ";p_{T} [GeV/c]; #theta [rad];N particles theta",
                                            pt_axis.size()-1, pt_axis.data(),
                                            theta_axis.size()-1, theta_axis.data(),
                                            npart_sector_axis.size()-1, npart_sector_axis.data());
  h3_theta_centrality_npart = new TH3F( "h3_theta_centrality_npart", ";#theta [rad];centrality (%);N particles theta",
                                            theta_axis.size()-1, theta_axis.data(),
                                                    centrality_axis.size()-1, centrality_axis.data(),
                                             npart_sector_axis.size()-1, npart_sector_axis.data());

  auto y_cm = data_header_->GetBeamRapidity();
  std::vector<double> pT_midrapidity;
  beta_cm_ = tanh(y_cm);
  out_file_->cd();
  std::cout << "Initialized" << std::endl;
}

void Yield::UserExec() {
  using AnalysisTree::Particle;
  // Reseting all of the event by event histograms
  h1_phi_theta_all_->Reset();
  h3_rec_theta_pT_phi_->Reset();
  h3_tru_theta_pT_phi_->Reset();

  h1_theta_all_->Reset();
  h2_rec_pT_theta_->Reset();
  h2_tru_pT_theta_->Reset();

  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  h1_centrality_->Fill( centrality );
  float y_beam = data_header_->GetBeamRapidity();
  std::vector<double> pT_midrapidity;
  h1_centrality_->Fill(centrality);

  auto rec_chi2_var = GetVar("mdc_vtx_tracks/chi2");
  auto rec_dca_xy_var = GetVar("mdc_vtx_tracks/dca_xy");
  auto rec_dca_z_var = GetVar("mdc_vtx_tracks/dca_z");
  auto tru_is_primary = GetVar("sim_tracks/is_primary");

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

    if( chi2 > 100.0 )
      continue;
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    if( pid != reference_pdg_code_ )
      continue;
    auto y = mom4.Rapidity() - y_beam;
    h3_rec_theta_pT_phi_->Fill( mom4.Theta(), mom4.Pt(), mom4.Phi() );
    h2_rec_pT_theta_->Fill( mom4.Pt(), mom4.Theta() );
  }
  int n_charged_tracks=0;
  for( auto particle : sim_particles_->Loop() ){
    auto mass = particle.DataT<Particle>()->GetMass();
    auto pid = particle.DataT<Particle>()->GetPid();
    auto is_prim = particle[tru_is_primary].GetBool();
    auto mom4 = particle.DataT<Particle>()->Get4MomentumByMass(mass);
    double charge=0.0;
    if( TDatabasePDG::Instance()->GetParticle( pid ) ){
      charge= TDatabasePDG::Instance()->GetParticle( pid )->Charge() / 3.0;
    }
    if( fabs(charge) < 0.01 )
      continue;
    h1_phi_theta_all_->Fill( mom4.Phi(), mom4.Theta() );
    h1_theta_all_->Fill( mom4.Theta() );
    if( !is_prim )
      continue;
    if( pid!=reference_pdg_code_ )
      continue;
    h3_tru_theta_pT_phi_->Fill( mom4.Theta(), mom4.Pt(), mom4.Phi() );
    h2_tru_pT_theta_->Fill( mom4.Pt(), mom4.Theta() );
  }
  for( int phi_bin=0; phi_bin < 36; phi_bin++ ){
    for( int theta_bin =0; theta_bin < h3_tru_theta_pT_phi_->GetNbinsX();theta_bin++ ){
      auto theta_mean = h3_tru_theta_pT_phi_->GetXaxis()->GetBinCenter(theta_bin +1);
      auto phi_mean = h3_tru_theta_pT_phi_->GetZaxis()->GetBinCenter(phi_bin +1);
      auto theta_bin_all = h1_phi_theta_all_->GetYaxis()->FindBin( theta_mean );
      auto phi_bin_all = h1_phi_theta_all_->GetXaxis()->FindBin( phi_mean );
      double npart_sector = h1_phi_theta_all_->GetBinContent(phi_bin_all, theta_bin_all);
      for (int pT_bin = 0; pT_bin < h3_tru_theta_pT_phi_->GetNbinsY(); ++pT_bin) {
        double pT_mean = h3_tru_theta_pT_phi_->GetYaxis()->GetBinCenter(pT_bin+1);
        double tru_npart_bin = h3_tru_theta_pT_phi_->GetBinContent(theta_bin +1, pT_bin+1, phi_bin +1 );
        if( fabs(tru_npart_bin) < 0.01 )
          continue;
        auto rec_npart_bin = h3_rec_theta_pT_phi_->GetBinContent(theta_bin +1, pT_bin+1, phi_bin +1 );
        auto efficiency = rec_npart_bin / tru_npart_bin;
        p3_theta_pT_npart_sector_->Fill(theta_mean, pT_mean, npart_sector-tru_npart_bin, efficiency );
      }
    }
    for( int pT_bin=0; pT_bin < h2_tru_pT_theta_->GetNbinsX(); pT_bin++ ){
      auto mean_pT = h2_tru_pT_theta_->GetXaxis()->GetBinCenter(pT_bin+1);
      for( int theta_bin=0; theta_bin<h2_tru_pT_theta_->GetNbinsY(); theta_bin++ ){
        auto mean_theta = h2_tru_pT_theta_->GetYaxis()->GetBinCenter(theta_bin+1);
        auto n_all_particles_theta = h1_theta_all_->GetBinContent( theta_bin+1 );
        auto n_tru_particles_theta_pT = h2_tru_pT_theta_->GetBinContent( pT_bin+1, theta_bin+1 );
        if( fabs(n_tru_particles_theta_pT) < 0.01  )
          continue;
        auto n_rec_particles_theta_pT = h2_rec_pT_theta_->GetBinContent( pT_bin+1, theta_bin+1 );
        auto efficiency = n_rec_particles_theta_pT / n_tru_particles_theta_pT;
        p3_pT_theta_npart_theta_->Fill( mean_pT, mean_theta, n_all_particles_theta - n_tru_particles_theta_pT, efficiency );
        h3_rec_pT_theta_npart_theta_->Fill( mean_pT, mean_theta, n_all_particles_theta - n_tru_particles_theta_pT, n_rec_particles_theta_pT );
        h3_tru_pT_theta_npart_theta_->Fill( mean_pT, mean_theta, n_all_particles_theta - n_tru_particles_theta_pT, n_tru_particles_theta_pT );
        h3_theta_centrality_npart->Fill( mean_theta, centrality, n_all_particles_theta - n_tru_particles_theta_pT );
      }
    }
  }

}
void Yield::UserFinish() {
  out_file_->cd();
  h1_centrality_->Write();

//  out_file_->mkdir("efficiency_projections");
//  out_file_->cd("efficiency_projections");

  out_file_->cd();
  p3_theta_pT_npart_sector_->Write();
  p3_pT_theta_npart_theta_->Write();
  h3_theta_centrality_npart->Write();
  h3_rec_pT_theta_npart_theta_->Write();
  h3_tru_pT_theta_npart_theta_->Write();
  std::cout << "Finished" << std::endl;
}
int Yield::CalculateNumberOfChargedTracks(std::vector<double> theta_range) {
  using AnalysisTree::Particle;
  int n_particles=0;
  for( auto particle : sim_particles_->Loop() ){
    auto pid = particle.DataT<Particle>()->GetPid();
    double charge=0.0;
    if( TDatabasePDG::Instance()->GetParticle( pid ) ){
      charge= TDatabasePDG::Instance()->GetParticle( pid )->Charge() / 3.0;
    }
    if( fabs(charge) < 0.01 )
      continue;
    auto mass = particle.DataT<Particle>()->GetMass();
    auto mom4 = particle.DataT<Particle>()->Get4MomentumByMass(mass);
    if( theta_range.front() < mom4.Theta() && mom4.Theta() < theta_range.back() )
      n_particles++;
  }
  return n_particles;
}
void Yield::InitEfficiency() {

}
