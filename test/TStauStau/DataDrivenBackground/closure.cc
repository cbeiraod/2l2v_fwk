#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TDirectory.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);

void closure()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptLogy();

  TFile plotter("/home/cms/cbeiraod/local-area/DDBkgFRPlots/plotter.root", "READ");
  TFile LIP_ffile("FROut.root", "READ");
  TFile LIP_pfile("PROut.root", "READ");

  TCanvas c1("c1", "c1", 800, 600);
  double IPM_flat_f = 0.4282;
  double IPM_flat_p = 0.7659;
  double LIP_flat_f = 0.4551;
  double LIP_flat_p = 0.7763;
  TFile outfile("closure.root", "RECREATE");

  std::vector<std::string> channels;
  channels.push_back("selected");
  channels.push_back("leadingE");
  channels.push_back("leadingMu");

  std::vector<std::string> plots;
  plots.push_back("MET");
  plots.push_back("TauMET");
  plots.push_back("etaSelectedTau");
  plots.push_back("ptSelectedTau");
  plots.push_back("ptSelectedTauExtended");
//  plots.push_back("ptetaSelectedTau");
//  plots.push_back("ptetaSelectedTauExtended");
//  plots.push_back("ptabsetaSelectedTau");
//  plots.push_back("ptabsetaSelectedTauExtended");
  plots.push_back("varetaSelectedTau");
  plots.push_back("varptSelectedTau");
  plots.push_back("varptSelectedTauExtended");
//  plots.push_back("varptetaSelectedTau");
//  plots.push_back("varptetaSelectedTauExtended");
//  plots.push_back("varptabsetaSelectedTau");
//  plots.push_back("varptabsetaSelectedTauExtended");

  std::stringstream buffer;


  if(plotter.Get("data") == NULL)
  {
    std::cout << " Unable to find data directory. Please fix it." << std::endl;
    return;
  }

  for(std::vector<std::string>::iterator channel = channels.begin(); channel != channels.end(); ++channel)
  {
    std::cout << "  Doing: " << *channel << std::endl;
    for(std::vector<std::string>::iterator plot = plots.begin(); plot!= plots.end(); ++plot)
    {
      std::cout << "    Getting plot: " << *plot << std::endl;

      TH1* Loose = NULL;
      TH1* Tight = NULL;
      Loose = static_cast<TH1*>(plotter.Get(("data/" + *channel + "_" + *plot).c_str()));
      Tight = static_cast<TH1*>(plotter.Get(("data/" + *channel + "_" + *plot + "Tight").c_str()));

      if(Loose == NULL)
      {
        std::cout << "      Unable to find the plot of " << *plot << " for the channel " << *channel << " in data. Continuing..." << std::endl;
        continue;
      }
      if(Tight == NULL)
      {
        std::cout << "      The tight selection for the plot of " << *plot << " for the channel " << *channel << " in data could not be found. Jumping this one" << std::endl;
        continue;
      }

      std::string baseName = *channel + "_" + *plot;

      baseName = ReplaceAll(baseName, " ", "");
      baseName = ReplaceAll(baseName, "+", "");

      Loose->SetName(baseName.c_str());
      Tight->SetName((baseName+"_tight").c_str());
      Loose->Write();
      Tight->Write();

      Loose->Add(Tight, -1.0);


      TH1* IPMestimate_flat = static_cast<TH1*>(Loose->Clone((baseName+"_IPMestimateFlat").c_str()));
      TH1* tmp = static_cast<TH1*>(Tight->Clone("tmp1"));
      double param1 = IPM_flat_f*IPM_flat_p/(IPM_flat_p-IPM_flat_f);
      double param2 = IPM_flat_f*(1-IPM_flat_p)/(IPM_flat_p-IPM_flat_f);
      IPMestimate_flat->Scale(param1);
      tmp->Scale(param2);
      IPMestimate_flat->Add(tmp, -1.0);
      IPMestimate_flat->Write();
      IPMestimate_flat->Draw("hist");
      Tight->Draw("same");
      c1.SaveAs((baseName+"_IPMclosureFlat.png").c_str());
      delete tmp;

      TH1* LIPestimate_flat = static_cast<TH1*>(Loose->Clone((baseName+"_LIPestimateFlat").c_str()));
      tmp = static_cast<TH1*>(Tight->Clone("tmp2"));
      param1 = LIP_flat_f*LIP_flat_p/(LIP_flat_p-LIP_flat_f);
      param2 = LIP_flat_f*(1-LIP_flat_p)/(LIP_flat_p-LIP_flat_f);
      LIPestimate_flat->Scale(param1);
      tmp->Scale(param2);
      LIPestimate_flat->Add(tmp, -1.0);
      LIPestimate_flat->Write();
      LIPestimate_flat->Draw("hist");
      Tight->Draw("same");
      c1.SaveAs((baseName+"_LIPclosureFlat.png").c_str());
      delete tmp;


      TH1* LIPestimate      = static_cast<TH1*>(Loose->Clone((baseName+"_LIPestimate").c_str()));
      tmp = static_cast<TH1*>(Tight->Clone("tmp3"));
      TH1* LIP_f = static_cast<TH1*>(LIP_ffile.Get(("data/data_" + *channel + "_" + *plot + "_FR").c_str()));
      TH1* LIP_p = static_cast<TH1*>(LIP_pfile.Get(("Z #rightarrow ll/Zrightarrowll_" + *channel + "_" + *plot + "_FR").c_str()));
      TH1* denominator = static_cast<TH1*>(LIP_p->Clone("denom"));
      denominator->Add(LIP_f, -1.0);
      LIP_p->Multiply(LIP_f);
      LIP_f->Add(LIP_p,-1.0);
      LIP_p->Divide(denominator); // param1 =     fp/(p-f)
      LIP_f->Divide(denominator); // param2 = f(1-p)/(p-f)
      LIPestimate->Multiply(LIP_p);
      tmp->Multiply(LIP_f);
      LIPestimate->Add(tmp, -1.0);
      LIPestimate->Write();
      LIPestimate->Draw("hist");
      Tight->Draw("same");
      c1.SaveAs((baseName+"_LIPclosure.png").c_str());
      delete tmp;
      delete denominator;
      delete LIP_p;
      delete LIP_f;
    }
  }

  std::cout << buffer.str();

//  TH1D* MetLoose = static_cast<TH1D*>(plotter.Get("W + Jets/selected_MET"));
//  TH1D* MetTight = static_cast<TH1D*>(plotter.Get("W + Jets/selected_METTight"));

//  MetLoose->Write();
//  MetTight->Write();

//  TH1D* METFR = static_cast<TH1D*>(MetTight->Clone("selected_MET_FR"));
//  METFR->Divide(MetLoose);

//  METFR->Write();
}

std::string ReplaceAll(std::string str, const std::string& from, const std::string& to)
{
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos)
  {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
  }
  return str;
}
