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

void overlapPlots()
{
  TFile hists("testOut.root", "READ");
  TCanvas c1("c1", "c1", 800, 600);

  std::vector<std::string> processes;
  processes.push_back("W + Jets");
  processes.push_back("data");

  std::vector<std::string> channels;
  channels.push_back("etau");
  channels.push_back("mutau");

  std::vector<std::string> plots;
  plots.push_back("MET");
  plots.push_back("etaSelectedTau");
  plots.push_back("ptSelectedTau");
  plots.push_back("ptSelectedTauExtended");


  for(std::vector<std::string>::iterator process = processes.begin(); process != processes.end(); ++process)
  {
    for(std::vector<std::string>::iterator plot = plots.begin(); plot != plots.end(); ++plot)
    {
      int count = 0;
      for(std::vector<std::string>::iterator channel = channels.begin(); channel != channels.end(); ++channel)
      {
        std::string baseName = *process + "_" + *channel + "_" + *plot;

        baseName = ReplaceAll(baseName, " ", "");
        baseName = ReplaceAll(baseName, "+", "");

        TH1* thisHist = static_cast<TH1*>(hists.Get((*process + "/" + baseName + "_FR").c_str()));
        if(thisHist == NULL)
          continue;

        thisHist->SetNameTitle(channel->c_str(), channel->c_str());

        switch(count%6)
        {
          case 0:
            thisHist->SetMarkerColor(kRed);
            break;
          case 1:
            thisHist->SetMarkerColor(kBlue);
            break;
          case 2:
            thisHist->SetMarkerColor(kGreen);
            break;
          case 3:
            thisHist->SetMarkerColor(kYellow);
            break;
          case 4:
            thisHist->SetMarkerColor(kMagenta);
            break;
          case 5:
            thisHist->SetMarkerColor(kCyan);
            break;
        }

        if(count == 0)
          thisHist->Draw();
        else
          thisHist->Draw("same");

        ++count;
      }

      c1.BuildLegend();
      std::string baseName = *process + "_" + *plot;
      baseName = ReplaceAll(baseName, " ", "");
      baseName = ReplaceAll(baseName, "+", "");
      c1.SaveAs((baseName+".png").c_str());
    }
  }

  for(std::vector<std::string>::iterator channel = channels.begin(); channel != channels.end(); ++channel)
  {
    for(std::vector<std::string>::iterator plot = plots.begin(); plot != plots.end(); ++plot)
    {
      int count = 0;
      for(std::vector<std::string>::iterator process = processes.begin(); process != processes.end(); ++process)
      {
        std::string baseName = *process + "_" + *channel + "_" + *plot;

        baseName = ReplaceAll(baseName, " ", "");
        baseName = ReplaceAll(baseName, "+", "");

        TH1* thisHist = static_cast<TH1*>(hists.Get((*process + "/" + baseName + "_FR").c_str()));
        if(thisHist == NULL)
          continue;

        thisHist->SetNameTitle(process->c_str(), process->c_str());

        switch(count%6)
        {
          case 0:
            thisHist->SetMarkerColor(kRed);
            break;
          case 1:
            thisHist->SetMarkerColor(kBlue);
            break;
          case 2:
            thisHist->SetMarkerColor(kGreen);
            break;
          case 3:
            thisHist->SetMarkerColor(kYellow);
            break;
          case 4:
            thisHist->SetMarkerColor(kMagenta);
            break;
          case 5:
            thisHist->SetMarkerColor(kCyan);
            break;
        }

        if(count == 0)
          thisHist->Draw();
        else
          thisHist->Draw("same");

        ++count;
      }

      c1.BuildLegend();
      std::string baseName = *channel + "_" + *plot;
      baseName = ReplaceAll(baseName, " ", "");
      baseName = ReplaceAll(baseName, "+", "");
      c1.SaveAs((baseName+".png").c_str());
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
