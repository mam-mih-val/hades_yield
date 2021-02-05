//
// Created by mikhail on 2/4/21.
//

#include <TFile.h>
#include <TH2F.h>
#include <stdexcept>

int main(int n_args, char** args){
  if( n_args < 3 )
    throw std::runtime_error( "Minimum 2 argumet expected." );
  std::string input_file{ args[1] };
  std::string output_file{ args[2] };
  auto file_in = TFile::Open( input_file.c_str() );
  std::vector<TH2F*> yields;
  TH1F* centrality;
  file_in->GetObject("centrality",centrality);
  int p=0;
  TH2F* histo;
  int bin=1;
  while(p<40){
    std::string histo_name{ "yield_"+std::to_string(p)+"-"+std::to_string(p+5) };
    file_in->GetObject(histo_name.c_str(), histo);
    auto events = centrality->GetBinContent(bin);
    histo->Scale(1.0/events);
    yields.push_back(histo);
    p+=5;
  }
  std::vector<TH2F*> efficiencies;
  auto file_efficiencies = TFile::Open( "/home/mikhail/hades_yield/efficiency/efficiency_protons.root" );
  p=2;
  while(p<40){
    std::string histo_name{ "efficiency_"+std::to_string(p) };
    file_efficiencies->GetObject(histo_name.c_str(), histo);
    efficiencies.push_back(histo);
    p+=5;
  }
  for( size_t i=0; i<std::size(yields); ++i ){
    yields.at(i)->Divide(efficiencies.at(i));
  }

  std::vector<TH2F*> symmetric_yields;
  std::vector<TH2F*> acceptance;

  float y_axis[16];
  for(int j=0; j<16; ++j){ y_axis[j]=-0.75f+0.1f* (float) j; }
  float pt_axis[]={0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
  p=0;
  for( auto histo : yields ){
    std::string histo_name{ "symmetric_"+std::to_string(p)+"-"+std::to_string(p+5) };
    symmetric_yields.push_back( new TH2F( histo_name.c_str(), ";y_{cm};pT [GeV/c]", 15, y_axis, 10, pt_axis ) );
    for( auto y_bin=1; y_bin <=15; y_bin++ ){
      for( auto pt_bin=1; pt_bin<=10; pt_bin++){
        auto content = histo->GetBinContent( y_bin, pt_bin );
        auto error = histo->GetBinError( y_bin, pt_bin );
        symmetric_yields.back()->SetBinContent(16-y_bin, pt_bin, content);
        symmetric_yields.back()->SetBinError(16-y_bin, pt_bin, error);
      }
    }
    symmetric_yields.back()->SetEntries(histo->GetEntries());
    histo_name = "acceptacne_"+std::to_string(p)+"-"+std::to_string(p+5);
    acceptance.push_back( (TH2F*) histo->Clone(histo_name.c_str()) );
    acceptance.back()->Add( symmetric_yields.back(), -1.0 );
    p+=5;
  }


  auto file_out = TFile::Open(output_file.c_str(), "recreate");
  file_out->mkdir("yields");
  file_out->cd("yields");
  for( auto histo : yields ){
    histo->Write();
  }
  file_out->cd("/");
  file_out->mkdir("symmetric");
  file_out->cd("symmetric");
  for( auto histo : symmetric_yields ){
    histo->Write();
  }
  file_out->cd("/");
  file_out->mkdir("acceptance");
  file_out->cd("acceptance");
  for( auto histo : acceptance ){
    histo->Write();
  }
  file_out->Close();
  return 0;
}