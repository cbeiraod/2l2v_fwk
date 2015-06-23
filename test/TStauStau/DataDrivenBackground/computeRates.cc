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

void computeRates()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile plotter("/home/cms/cbeiraod/local-area/DDBkgPlots/plotter.root", "READ");
  TCanvas c1("c1", "c1", 800, 600);
  TFile outfile("RatesOut.root", "RECREATE");

  std::vector<std::string> saveToFile;
  saveToFile.push_back("data-Zprompt_InvMET_OS_ptSelectedTau");
  saveToFile.push_back("data-Zprompt_InvMET_OS_etaSelectedTau");
  saveToFile.push_back("Zrightarrowll_InvMET_OS_ptSelectedTau");
  saveToFile.push_back("Zrightarrowll_InvMET_OS_etaSelectedTau");

  std::vector<std::string> processes;
  processes.push_back("W + Jets");
  processes.push_back("data");
  processes.push_back("data-Z");
  processes.push_back("data-prompt");
  processes.push_back("data-Zprompt");
  processes.push_back("Z #rightarrow ll");

  std::vector<std::string> channels;
//  channels.push_back("selected");
//  channels.push_back("leadingE");
//  channels.push_back("leadingMu");
  channels.push_back("OS");
  channels.push_back("SS");
  channels.push_back("InvMET");
  channels.push_back("MET_OS");
  channels.push_back("MET_SS");
  channels.push_back("InvMET_OS");
  channels.push_back("InvMET_SS");
  channels.push_back("mm");
  channels.push_back("pp");
  channels.push_back("Gluon");
  channels.push_back("Quark");
  channels.push_back("OS_Gluon");
  channels.push_back("OS_Quark");
  channels.push_back("SS_Gluon");
  channels.push_back("SS_Quark");
  channels.push_back("chargeSymmetric");
  channels.push_back("OS_Prompt");
  channels.push_back("SS_Prompt");
  channels.push_back("Prompt");
  channels.push_back("1Prong");
  channels.push_back("3Prong");
  channels.push_back("OS_1Prong");
  channels.push_back("OS_3Prong");
  channels.push_back("SS_1Prong");
  channels.push_back("SS_3Prong");
//  channels.push_back("OS_leadingE");
//  channels.push_back("OS_leadingMu");
//  channels.push_back("SS_leadingE");
//  channels.push_back("SS_leadingMu");

  std::vector<std::string> plots;
//  plots.push_back("MET");
//  plots.push_back("TauMET");
  plots.push_back("etaSelectedTau");
  plots.push_back("ptSelectedTau");
  plots.push_back("ptSelectedTauExtended");
//  plots.push_back("cosPhiSelectedTau");
  plots.push_back("ptetaSelectedTau");
  plots.push_back("ptetaSelectedTauExtended");
  plots.push_back("ptabsetaSelectedTau");
  plots.push_back("ptabsetaSelectedTauExtended");
//  plots.push_back("varetaSelectedTau");
//  plots.push_back("varptSelectedTau");
//  plots.push_back("varptSelectedTauExtended");
//  plots.push_back("varptetaSelectedTau");
//  plots.push_back("varptetaSelectedTauExtended");
//  plots.push_back("varptabsetaSelectedTau");
//  plots.push_back("varptabsetaSelectedTauExtended");
  plots.push_back("etaSelectedLep");
  plots.push_back("ptSelectedLep");
  plots.push_back("ptSelectedLepExtended");
//  plots.push_back("cosPhiSelectedLep");

  std::stringstream buffer;


  for(std::vector<std::string>::iterator process = processes.begin(); process != processes.end(); ++process)
  {
    std::cout << "Processing fake rate for: " << *process << std::endl;
    if(plotter.Get(process->c_str()) == NULL && *process != "data-Z" && *process != "data-prompt" && *process != "data-Zprompt")
    {
      std::cout << " Unable to find " << *process << " directory. Skipping it." << std::endl;
      continue;
    }

    outfile.mkdir(process->c_str());
    outfile.cd(process->c_str());

    for(std::vector<std::string>::iterator channel = channels.begin(); channel != channels.end(); ++channel)
    {
      double numerator = 0, denominator = 0;
      double numeratorUnc = 0, denominatorUnc = 0;

      std::cout << "  Doing: " << *channel << std::endl;
      for(std::vector<std::string>::iterator plot = plots.begin(); plot!= plots.end(); ++plot)
      {
        std::cout << "    Getting plot: " << *plot << std::endl;

        TH1* Loose = NULL;
        TH1* Tight = NULL;
        if(*process == "data-Z" || *process == "data-prompt" || *process == "data-Zprompt")
        {
          if(*channel != "chargeSymmetric")
          {
            Loose = static_cast<TH1*>(plotter.Get(("data/" + *channel + "_" + *plot).c_str()));
            Tight = static_cast<TH1*>(plotter.Get(("data/" + *channel + "_" + *plot + "Tight").c_str()));
          }
          else
          {
            Loose = static_cast<TH1*>(plotter.Get(("data/lepPlus_" + *plot).c_str()));
            Tight = static_cast<TH1*>(plotter.Get(("data/lepPlus_" + *plot + "Tight").c_str()));

            TH1* temp = static_cast<TH1*>(plotter.Get(("data/lepMinus_" + *plot).c_str()));
            if(Loose == NULL || temp == NULL)
              Loose = NULL;
            else
              Loose->Add(temp, -1);

            temp = static_cast<TH1*>(plotter.Get(("data/lepMinus_" + *plot + "Tight").c_str()));
            if(Tight == NULL || temp == NULL)
              Tight = NULL;
            else
              Tight->Add(temp, -1);
          }
        }
        else
        {
          if(*channel != "chargeSymmetric")
          {
            Loose = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_" + *plot).c_str()));
            Tight = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_" + *plot + "Tight").c_str()));
          }
          else
          {
            Loose = static_cast<TH1*>(plotter.Get((*process + "/lepPlus_" + *plot).c_str()));
            Tight = static_cast<TH1*>(plotter.Get((*process + "/lepPlus_" + *plot + "Tight").c_str()));

            TH1* temp = static_cast<TH1*>(plotter.Get((*process + "/lepMinus_" + *plot).c_str()));
            if(Loose == NULL || temp == NULL)
              Loose = NULL;
            else
              Loose->Add(temp, -1);

            temp = static_cast<TH1*>(plotter.Get((*process + "/lepMinus_" + *plot + "Tight").c_str()));
            if(Tight == NULL || temp == NULL)
              Tight = NULL;
            else
              Tight->Add(temp, -1);
          }
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

        if(Loose == NULL || Tight == NULL)
        {
          std::cout << "      Unable to find " << *plot << " plot for the " << *channel << " channel (" << "data/" + *channel + "_" + *plot  << " or " << "data/" + *channel + "_" + *plot + "Tight" << ")." << std::endl;
          continue;
        }


        if(Loose != NULL && Tight != NULL)
        {
          if(*process == "data-Z")
          {
            TH1* temp = NULL;
            if(*channel != "chargeSymmetric")
            {
              temp = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/" + *channel + "_" + *plot).c_str()));
            }
            else
            {
              temp = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/lepPlus_" + *plot).c_str()));
              TH1* temp2 = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/lepMinus_" + *plot).c_str()));
              if(temp != NULL && temp2 != NULL)
                temp->Add(temp2, -1);
              else
                temp = NULL;
            }
            if(temp != NULL)
              Loose->Add(temp, -1);
            else
              Loose = NULL;

            if(*channel != "chargeSymmetric")
            {
              temp = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/" + *channel + "_" + *plot + "Tight").c_str()));
            }
            else
            {
              temp = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/lepPlus_" + *plot + "Tight").c_str()));
              TH1* temp2 = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/lepMinus_" + *plot + "Tight").c_str()));
              if(temp != NULL && temp2 != NULL)
                temp->Add(temp2, -1);
              else
                temp = NULL;
            }
            if(temp != NULL)
              Tight->Add(temp, -1);
            else
              Tight = NULL;
          }

          if(*process == "data-prompt" || *process == "data-Zprompt")
          {
            std::vector<std::string> toRemove;
            toRemove.push_back("Z #rightarrow ll");
            toRemove.push_back("W + Jets");
            toRemove.push_back("VV-VVV");
            toRemove.push_back("t#bar{t}");
            toRemove.push_back("Single top");
            TH1* temp = NULL;

            for(std::vector<std::string>::iterator removing = toRemove.begin(); removing != toRemove.end(); ++removing)
            {
              if(*channel == "chargeSymmetric")
              {
                Loose = NULL;
                Tight = NULL;
                break;
              }


              std::string prompt = "_Prompt_";
              if(*process == "data-Zprompt" && *removing == "Z #rightarrow ll")
                prompt = "_";

              temp = static_cast<TH1*>(plotter.Get((*removing + "/" + *channel + prompt + *plot).c_str()));
              if(temp != NULL)
                Loose->Add(temp, -1);

              temp = static_cast<TH1*>(plotter.Get((*removing + "/" + *channel + prompt + *plot + "Tight").c_str()));
              if(temp != NULL)
                Tight->Add(temp, -1);
            }
          }
        }

        if(Loose == NULL)
        {
          std::cout << "      Problem when subtracting the plot of " << *plot << " for the channel " << *channel << " in the process " << *process << ". Continuing..." << std::endl;
          continue;
        }
        if(Tight == NULL)
        {
          std::cout << "      Problem when subtracting the tight plot of " << *plot << " for the channel " << *channel << " in the process " << *process << ". Continuing..." << std::endl;
          continue;
        }

        std::string baseName = *process + "_" + *channel + "_" + *plot;

        baseName = ReplaceAll(baseName, " ", "");
        baseName = ReplaceAll(baseName, "+", "");
        baseName = ReplaceAll(baseName, "#", "");

        Loose->SetName(baseName.c_str());
        Tight->SetName((baseName+"_tight").c_str());
        bool save = false;
        for(std::vector<std::string>::iterator name = saveToFile.begin(); name != saveToFile.end(); ++name)
        {
          if(*name == baseName)
           save = true;
        }
        if(saveToFile.size() == 0 || save)
        {
          Loose->Write();
          Tight->Write();
        }

        TH1* FR = static_cast<TH1*>(Tight->Clone((baseName+"_FR").c_str()));
        FR->GetYaxis()->SetTitle("FR");
        FR->Divide(Loose);
        if(!FR->InheritsFrom("TH2"))
          FR->SetAxisRange(0.35,0.6,"Y");
        else
          FR->SetAxisRange(0.35,0.6,"Z");
        if(saveToFile.size() == 0 || save)
          FR->Write();
        FR->Draw();
        c1.SaveAs((baseName+"_FR.png").c_str());

        TAxis* fXaxis = Tight->GetXaxis();
        numerator = Tight->IntegralAndError(fXaxis->GetFirst(), fXaxis->GetLast(), numeratorUnc);
        TAxis* gXaxis = Loose->GetXaxis();
        denominator = Loose->IntegralAndError(gXaxis->GetFirst(), gXaxis->GetLast(), denominatorUnc);
      }

      if(numerator != 0 && denominator != 0)
      {
        buffer << *process << ":" << *channel << std::endl;
        buffer << "Numerator: " << numerator  << "+-" << numeratorUnc << "; Denominator: " << denominator << "+-" << denominatorUnc << std::endl;
        double uncertainty = numeratorUnc/static_cast<double>(denominator) * numeratorUnc/static_cast<double>(denominator);
        uncertainty += denominatorUnc*numerator/(denominator * denominator) * denominatorUnc*numerator/(denominator * denominator);
        uncertainty = sqrt(uncertainty);
        buffer << "Ratio: " << numerator/static_cast<double>(denominator) << "\\pm" << uncertainty << std::endl << std::endl;
      }
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
