#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TDirectory.h"

#include <iostream>
#include <string>
#include <vector>

std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);

void computeFR()
{
  TFile plotter("/home/cms/cbeiraod/local-area/DDBkgPlots/plotter.root", "READ");
  TCanvas c1("c1", "c1", 800, 600);
  TFile outfile("testOut.root", "RECREATE");

  std::vector<std::string> processes;
  processes.push_back("W + Jets");
  processes.push_back("data");

  std::vector<std::string> channels;
  channels.push_back("selected");
  channels.push_back("etau");
  channels.push_back("mutau");

  std::vector<std::string> plots;
  plots.push_back("MET");
  plots.push_back("etaSelectedTau");
  plots.push_back("ptSelectedTau");
  plots.push_back("ptSelectedTauExtended");


  for(std::vector<std::string>::iterator process = processes.begin(); process != processes.end(); ++process)
  {
    std::cout << "Processing fake rate for: " << *process << std::endl;
    if(plotter.Get(process->c_str()) == NULL)
    {
      std::cout << " Unable to find " << *process << " directory. Skipping it." << std::endl;
      continue;
    }

    outfile.mkdir(process->c_str());
    outfile.cd(process->c_str());

    for(std::vector<std::string>::iterator channel = channels.begin(); channel != channels.end(); ++channel)
    {
      std::cout << "  Doing: " << *channel << std::endl;
      for(std::vector<std::string>::iterator plot = plots.begin(); plot!= plots.end(); ++plot)
      {
        std::cout << "    Getting plot: " << *plot << std::endl;

        TH1* Loose = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_" + *plot).c_str()));
        TH1* Tight = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_" + *plot + "Tight").c_str()));

        if(Loose == NULL)
        {
          std::cout << "      Unable to find the plot of " << *plot << " for the channel " << *channel << " in the process " << *process << ". Continuing..." << std::endl;
          continue;
        }
        if(Tight == NULL)
        {
          std::cout << "      The tight selection for the plot of " << *plot << " for the channel " << *channel << " in the process " << *process << " could not be found. Jumping this one" << std::endl;
          continue;
        }

        std::string baseName = *process + "_" + *channel + "_" + *plot;

        baseName = ReplaceAll(baseName, " ", "");
        baseName = ReplaceAll(baseName, "+", "");

        Loose->SetName(baseName.c_str());
        Tight->SetName((baseName+"_tight").c_str());
        Loose->Write();
        Tight->Write();

        TH1* FR = static_cast<TH1*>(Tight->Clone((baseName+"_FR").c_str()));
        FR->Divide(Loose);
        FR->Write();
        FR->Draw();
        c1.SaveAs((baseName+"_FR.png").c_str());
      }
    }
  }

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
