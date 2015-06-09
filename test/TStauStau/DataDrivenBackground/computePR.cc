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

void computePR()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile plotter("/home/cms/cbeiraod/local-area/DDBkgPlots/plotter.root", "READ");
  TCanvas c1("c1", "c1", 800, 600);
  TFile outfile("PROut.root", "RECREATE");

  std::vector<std::string> processes;
  processes.push_back("Z #rightarrow ll");
//  processes.push_back("W + Jets");
//  processes.push_back("data");
//  processes.push_back("data-Z");

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
  plots.push_back("cosPhiSelectedTau");
  plots.push_back("ptetaSelectedTau");
  plots.push_back("ptetaSelectedTauExtended");
  plots.push_back("ptabsetaSelectedTau");
  plots.push_back("ptabsetaSelectedTauExtended");
  plots.push_back("varetaSelectedTau");
  plots.push_back("varptSelectedTau");
  plots.push_back("varptSelectedTauExtended");
  plots.push_back("varptetaSelectedTau");
  plots.push_back("varptetaSelectedTauExtended");
  plots.push_back("varptabsetaSelectedTau");
  plots.push_back("varptabsetaSelectedTauExtended");
  plots.push_back("etaSelectedLep");
  plots.push_back("ptSelectedLep");
  plots.push_back("ptSelectedLepExtended");
  plots.push_back("cosPhiSelectedLep");

  stringstream buffer;


  for(std::vector<std::string>::iterator process = processes.begin(); process != processes.end(); ++process)
  {
    std::cout << "Processing fake rate for: " << *process << std::endl;
    if(plotter.Get(process->c_str()) == NULL && *process != "data-Z")
    {
      std::cout << " Unable to find " << *process << " directory. Skipping it." << std::endl;
      continue;
    }

    outfile.mkdir(process->c_str());
    outfile.cd(process->c_str());

    double numerator = 0, denominator = 0;
    double numeratorUnc = 0, denominatorUnc = 0;

    for(std::vector<std::string>::iterator channel = channels.begin(); channel != channels.end(); ++channel)
    {
      std::cout << "  Doing: " << *channel << std::endl;
      for(std::vector<std::string>::iterator plot = plots.begin(); plot!= plots.end(); ++plot)
      {
        std::cout << "    Getting plot: " << *plot << std::endl;

        TH1* Loose = NULL;
        TH1* Tight = NULL;
        if(*process == "data-Z")
        {
          Loose = static_cast<TH1*>(plotter.Get(("data/" + *channel + "_" + *plot).c_str()));
          Tight = static_cast<TH1*>(plotter.Get(("data/" + *channel + "_" + *plot + "Tight").c_str()));
          Loose->Sumw2();
          Tight->Sumw2();

          if(Loose != NULL && Tight != NULL)
          {
            TH1* temp = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/" + *channel + "_" + *plot).c_str()));
            temp->Sumw2();
            if(temp != NULL)
              Loose->Add(temp, -1);
            else
              Loose = NULL;

            temp = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/" + *channel + "_" + *plot + "Tight").c_str()));
            temp->Sumw2();
            if(temp != NULL)
              Tight->Add(temp, -1);
            else
              Tight = NULL;
          }
        }
        else
        {
          Loose = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_" + *plot).c_str()));
          Tight = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_" + *plot + "Tight").c_str()));
          Loose->Sumw2();
          Tight->Sumw2();
        }

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
        baseName = ReplaceAll(baseName, "#", "");

        Loose->SetName(baseName.c_str());
        Tight->SetName((baseName+"_tight").c_str());
        Loose->Write();
        Tight->Write();

        TH1* FR = static_cast<TH1*>(Tight->Clone((baseName+"_FR").c_str()));
        FR->GetYaxis()->SetTitle("FR");
        FR->Divide(Loose);
        FR->Write();
        FR->Draw();
        c1.SaveAs((baseName+"_FR.png").c_str());

        //numerator = Tight->Integral();
        //denominator = Loose->Integral();
        TAxis* fXaxis = Tight->GetXaxis();
        numerator = Tight->IntegralAndError(fXaxis->GetFirst(), fXaxis->GetLast(), numeratorUnc);
        TAxis* gXaxis = Loose->GetXaxis();
        denominator = Loose->IntegralAndError(gXaxis->GetFirst(), gXaxis->GetLast(), denominatorUnc);
      }

      buffer << *process << ":" << *channel << std::endl;
      buffer << "Numerator: " << numerator  << "+-" << numeratorUnc << "; Denominator: " << denominator << "+-" << denominatorUnc << std::endl;
      double uncertainty = numeratorUnc/static_cast<double>(denominator) * numeratorUnc/static_cast<double>(denominator);
      uncertainty += denominatorUnc*numerator/(denominator * denominator) * denominatorUnc*numerator/(denominator * denominator);
      uncertainty = sqrt(uncertainty);
      buffer << "Ratio: " << numerator/static_cast<double>(denominator) << "\\pm" << uncertainty << std::endl << std::endl;
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
