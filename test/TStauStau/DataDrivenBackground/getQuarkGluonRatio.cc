#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "../../../interface/doubleWithUncertainty.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);

void getQuarkGluonRatio()
{
  gSystem->Load("libUserCodellvv_fwk");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile plotter("/home/cms/cbeiraod/local-area/DDBkgPlots/plotter.root", "READ");
  TCanvas c1("c1", "c1", 800, 600);
  TFile outfile("RatioOut.root", "RECREATE");

  std::vector<std::string> processes;
  processes.push_back("W + Jets");
  processes.push_back("Z #rightarrow ll");
  processes.push_back("QCD");
  processes.push_back("#gamma + Jets");
  processes.push_back("VV-VVV");
  processes.push_back("t#bar{t}");
  processes.push_back("Single top");
  processes.push_back("total");
  processes.push_back("total-Z");

  std::vector<std::string> channels;
/*  channels.push_back("OS");
  channels.push_back("SS");
  channels.push_back("chargeSymmetric");
  channels.push_back("InvMET");// */
  channels.push_back("MET_OS");
//  channels.push_back("MET_SS");// */
  channels.push_back("InvMET_OS");
/*  channels.push_back("InvMET_SS");
  channels.push_back("mm");
  channels.push_back("pp");// */

  std::stringstream buffer;

  for(std::vector<std::string>::iterator process = processes.begin(); process != processes.end(); ++process)
  {
    std::cout << "Processing ratio for: " << *process << std::endl;
    if(plotter.Get(process->c_str()) == NULL && *process != "total" && *process != "total-Z")
    {
      std::cout << " Unable to find " << *process << " directory. Skipping it." << std::endl;
      continue;
    }

    outfile.mkdir(process->c_str());
    outfile.cd(process->c_str());

    for(std::vector<std::string>::iterator channel = channels.begin(); channel != channels.end(); ++channel)
    {
      std::cout << "  Doing: " << *channel << std::endl;
      TH1* Loose = NULL;
      if(*process != "total" && *process != "total-Z")
      {
        if(*channel != "chargeSymmetric")
        {
          Loose = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_tauTypes").c_str()));
        }
        else
        {
          Loose = static_cast<TH1*>(plotter.Get((*process + "/lepPlus_tauTypes").c_str()));
          TH1* temp = static_cast<TH1*>(plotter.Get((*process + "/lepMinus_tauTypes").c_str()));
          if(Loose == NULL || temp == NULL)
            Loose = NULL;
          else
            Loose->Add(temp, -1);
        }
      }
      else
      {
        std::vector<std::string> toAdd;
        toAdd.push_back("Single top");
        toAdd.push_back("t#bar{t}");
        toAdd.push_back("VV/VVV");
        toAdd.push_back("#gamma + Jets");
        toAdd.push_back("QCD");
        if(*process == "total")
          toAdd.push_back("Z #rightarrow ll");
        toAdd.push_back("W + Jets");

        if(*channel != "chargeSymmetric")
        {
          for(std::vector<std::string>::iterator adding = toAdd.begin(); adding != toAdd.end(); ++adding)
          {
            TH1* temp = static_cast<TH1*>(plotter.Get((*adding + "/" + *channel + "_tauTypes").c_str()));
            if(temp != NULL)
            {
              if(Loose == NULL)
              {
                Loose = temp;
              }
              else
              {
                Loose->Add(temp, 1);
              }
            }
          }
        }
        else
        {
          Loose = NULL;
        }
      }

      bool skip = false;
      if(Loose == NULL)
      {
        std::cout << "    Unable to find tauTypes plot for the channel " << *channel << " in the process " << *process << ". Continuing..." << std::endl;
        skip = true;
      }
      if(skip)
        continue;

      std::string baseName = *process + "_" + *channel;

      baseName = ReplaceAll(baseName, " ", "");
      baseName = ReplaceAll(baseName, "+", "");
      baseName = ReplaceAll(baseName, "#", "");

      Loose->SetName(baseName.c_str());
      Loose->Write();

      buffer << *process << ":" << *channel << "\n";

      doubleUnc prompt;
      doubleUnc wrongSign;
      doubleUnc quark;
      doubleUnc gluon;
      doubleUnc unknown;
      doubleUnc total;

      prompt       = Loose->GetBinContent(2);
      prompt.setUncertainty(Loose->GetBinError(2));
      wrongSign    = Loose->GetBinContent(3);
      wrongSign.setUncertainty(Loose->GetBinError(3));
      gluon        = Loose->GetBinContent(4);
      gluon.setUncertainty(Loose->GetBinError(4));
      quark        = Loose->GetBinContent(5);
      quark.setUncertainty(Loose->GetBinError(5));
      unknown      = Loose->GetBinContent(6);
      unknown.setUncertainty(Loose->GetBinError(6));
      total = quark + gluon + unknown;

      buffer << "  Loose:\n";
      buffer << "    promp      - " << prompt << "\n";
      buffer << "    wrong sign - " << wrongSign << "\n";
      buffer << "    ______________________________________________\n";
      buffer << "    quark      - " << quark << "\n";
      buffer << "    gluon      - " << gluon << "\n";
      buffer << "    unknown    - " << unknown << "\n";
      buffer << "    total      - " << total << "\n";
      buffer << "\n";

      if(total.value() != 0)
      {
        buffer << "    quark component: " << quark/total << "\n";
        buffer << "    gluon component: " << gluon/total << "\n";
        buffer << "    unknown component: " << unknown/total << "\n";
        buffer << "\n";
      }

      if(gluon.value() != 0)
      {
        buffer << "    quark/gluon: " << quark/gluon << "\n";
      }
      if(quark.value() != 0)
      {
        buffer << "    gluon/quark: " << gluon/quark << "\n";
      }
      buffer << "\n\n";

/*        TH1* FR = static_cast<TH1*>(Tight->Clone((baseName+"_FR").c_str()));
        FR->GetYaxis()->SetTitle("FR");
        FR->Divide(Loose);
//        if(!FR->InheritsFrom("TH2"))
          FR->SetAxisRange(0.35,0.6,"Y");
        FR->Write();
        FR->Draw();
        c1.SaveAs((baseName+"_FR.png").c_str());

        TAxis* fXaxis = Tight->GetXaxis();
        numerator = Tight->IntegralAndError(fXaxis->GetFirst(), fXaxis->GetLast(), numeratorUnc);
        TAxis* gXaxis = Loose->GetXaxis();
        denominator = Loose->IntegralAndError(gXaxis->GetFirst(), gXaxis->GetLast(), denominatorUnc);

      buffer << *process << ":" << *channel << std::endl;
      buffer << "Numerator: " << numerator  << "+-" << numeratorUnc << "; Denominator: " << denominator << "+-" << denominatorUnc << std::endl;
      double uncertainty = numeratorUnc/static_cast<double>(denominator) * numeratorUnc/static_cast<double>(denominator);
      uncertainty += denominatorUnc*numerator/(denominator * denominator) * denominatorUnc*numerator/(denominator * denominator);
      uncertainty = sqrt(uncertainty);
      buffer << "Ratio: " << numerator/static_cast<double>(denominator) << "\\pm" << uncertainty << std::endl << std::endl;



      for(std::vector<std::string>::iterator plot = plots.begin(); plot!= plots.end(); ++plot)
      {
        std::cout << "    Getting plot: " << *plot << std::endl;

        TH1* Loose = NULL;
        TH1* Tight = NULL;
        if(*process == "data-Z" || *process == "data-prompt")
        {
          Loose = static_cast<TH1*>(plotter.Get(("data/" + *channel + "_" + *plot).c_str()));
          Tight = static_cast<TH1*>(plotter.Get(("data/" + *channel + "_" + *plot + "Tight").c_str()));
        }
        else
        {
          Loose = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_" + *plot).c_str()));
          Tight = static_cast<TH1*>(plotter.Get((*process + "/" + *channel + "_" + *plot + "Tight").c_str()));
        }

        if(Loose != NULL && Tight != NULL)
        {
          if(*process == "data-Z")
          {
            TH1* temp = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/" + *channel + "_" + *plot).c_str()));
//            temp->Sumw2();
            if(temp != NULL)
              Loose->Add(temp, -1);
            else
              Loose = NULL;

            temp = static_cast<TH1*>(plotter.Get(("Z #rightarrow ll/" + *channel + "_" + *plot + "Tight").c_str()));
//            temp->Sumw2();
            if(temp != NULL)
              Tight->Add(temp, -1);
            else
              Tight = NULL;
          }
          if(*process == "data-prompt")
          {
            std::vector<std::string> toRemove;
            TH1* temp = NULL;

            for(std::vector<std::string>::iterator removing = toRemove.begin(); removing != toRemove.end(); ++removing)
            {
              temp = static_cast<TH1*>(plotter.Get((*removing + "/" + *channel + "_Prompt_" + *plot).c_str()));
              if(temp != NULL)
                Loose->Add(temp, -1);
//              else
//                Loose = NULL;

              temp = static_cast<TH1*>(plotter.Get((*removing + "/" + *channel + "_Prompt_" + *plot + "Tight").c_str()));
              if(temp != NULL)
                Tight->Add(temp, -1);
//              else
//                Tight = NULL;

              if(Loose == NULL || Tight == NULL)
                break;
            }
          }
        }

        std::string baseName = *process + "_" + *channel + "_" + *plot;

        baseName = ReplaceAll(baseName, " ", "");
        baseName = ReplaceAll(baseName, "+", "");

        Loose->SetName(baseName.c_str());
        Tight->SetName((baseName+"_tight").c_str());
        Loose->Write();
        Tight->Write();

        TH1* FR = static_cast<TH1*>(Tight->Clone((baseName+"_FR").c_str()));
        FR->GetYaxis()->SetTitle("FR");
        FR->Divide(Loose);
//        if(!FR->InheritsFrom("TH2"))
          FR->SetAxisRange(0.35,0.6,"Y");
        FR->Write();
        FR->Draw();
        c1.SaveAs((baseName+"_FR.png").c_str());

        TAxis* fXaxis = Tight->GetXaxis();
        numerator = Tight->IntegralAndError(fXaxis->GetFirst(), fXaxis->GetLast(), numeratorUnc);
        TAxis* gXaxis = Loose->GetXaxis();
        denominator = Loose->IntegralAndError(gXaxis->GetFirst(), gXaxis->GetLast(), denominatorUnc);
      }// */

    }
  }

  std::cout << buffer.str() << std::flush;
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
