#include "TROOT.h"
#include "TFile.h"
//#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <iostream>
#include <sstream>
#include <string>

bool pairCompare(const std::pair<int, double>& first, const std::pair<int, double>& second)
{
  return first.second < second.second;
}

int CompareLimits()
{
  std::vector<std::string> files;
  files.push_back("~/Documents/CMS/Stau/UpperLimits/deltaM30.root");
  files.push_back("~/Documents/CMS/Stau/UpperLimits/deltaM70.root");
  files.push_back("~/Documents/CMS/Stau/UpperLimits/deltaM120.root");
//  files.push_back("~/Documents/CMS/Stau/UpperLimits/deltaM180.root");
//  files.push_back("~/Documents/CMS/Stau/UpperLimits/deltaM240.root");
  
  std::vector<TH2D*> hists;
  for(int i = 0; i < files.size(); ++i)
  {
    TFile inFile(files[i].c_str(), "READ");
    gROOT->cd();
    
    if(inFile.IsOpen())
    {
      std::cout << "Input file is open" << std::endl;
      
      TH2D* hist = static_cast<TH2D*>(inFile.Get("twodplot"));
      if(hist != NULL)
      {
        std::stringstream tmp;
        std::string tmpstr;
        tmp << "Selection";
        tmp << i;
        tmp >> tmpstr;
        hists.push_back(static_cast<TH2D*>(hist->Clone(tmpstr.c_str())));
      }
      
      inFile.Close();
    }
  }
  
  gROOT->cd();
  std::cout << "Got the histograms" << std::endl;
  
  if(hists.size() != 0)
  {
    TH2D* newHist = static_cast<TH2D*>(hists[0]->Clone("template"));
    newHist->Clear();
    
    TH2D* bestSelection = static_cast<TH2D*>(newHist->Clone("BestSelection"));
    bestSelection->SetTitle("Best Selection;M_{Stau};M_{Neut};Selection");
    TH2D* ratio2Second = static_cast<TH2D*>(newHist->Clone("RatioToSecond"));
    ratio2Second->SetTitle("Ratio 2nd best selection To 1st Best;M_{Stau};M_{Neut};Ratio");
    
    TH2D* proposedVs2 = static_cast<TH2D*>(newHist->Clone("3SelVs2Sel"));
    proposedVs2->SetTitle("3 SR vs 2 SR;M_{Stau};M_{Neut};Ratio");
    
    std::cout << "Trying to get the bins" << std::endl;
    int nbinsx = newHist->GetNbinsX();
    int nbinsy = newHist->GetNbinsY();
    std::cout << "Got the bins" << std::endl;
    
    for(int binx = 1; binx <= nbinsx; ++binx)
    {
      for(int biny = 1; biny <= nbinsy; ++biny)
      {
        std::vector<std::pair<int, double> > limits;
        for(int i = 0; i < hists.size(); ++i)
        {
          double binContent = hists[i]->GetBinContent(binx, biny);
          if(binContent != 0)
          {
            limits.push_back(std::pair<int, double>(i+1, binContent));
          }
        }
        
        if(limits.size() != 0)
        {
          std::sort(limits.begin(), limits.end(), pairCompare);
          
          bestSelection->SetBinContent(binx, biny, limits[0].first);
          if(limits.size() > 1)
            ratio2Second->SetBinContent(binx, biny, limits[1].second/limits[0].second);
        
          bool contentSet = false;
          double content = 0;
          double MStau = (binx+4)*10;
          double MNeut = (biny-1)*10;
          double TripleSRLimit = 0;
          double DoubleSRLimit = 0;

          if(MStau-MNeut <= 70)
          {
            TripleSRLimit = hists[0]->GetBinContent(binx, biny);
          }
          else
          {
            if(MStau-MNeut <= 160)
            {
              TripleSRLimit = hists[1]->GetBinContent(binx, biny);
            }
            else
            {
              TripleSRLimit = hists[2]->GetBinContent(binx, biny);
            }
          }
        
          if(MStau-MNeut <= 70)
          {
            DoubleSRLimit = hists[0]->GetBinContent(binx, biny);
          }
          else
          {
            DoubleSRLimit = hists[2]->GetBinContent(binx, biny);
          }
      
          if(TripleSRLimit != 0 && DoubleSRLimit != 0)
          {
            contentSet = true;
            content = DoubleSRLimit/TripleSRLimit;
            std::cout << "MStau: " << MStau << "; MNeut: " << MNeut << "; 3SR: " << TripleSRLimit << "; 2SR: " << DoubleSRLimit << "; Limit ratio: " << content << std::endl;
          }

          if(contentSet)
          {
//            proposedVs2->SetBinContent(binx, biny, MStau-MNeut);
            proposedVs2->SetBinContent(binx, biny, content);
          }
        }
      }
    }
    
    TCanvas c1("c1", "c1", 800, 600);
    
    bestSelection->Draw("colz");
    c1.SaveAs("bestSelection.png");
    
    proposedVs2->Draw("colz");
    c1.SaveAs("proposedVs2.png");
    
    c1.SetLogz();
    ratio2Second->Draw("colz");
    ratio2Second->GetZaxis()->SetRangeUser(1, 100);
    c1.SaveAs("ratio2Second.png");
  }
  
  return 0;
}