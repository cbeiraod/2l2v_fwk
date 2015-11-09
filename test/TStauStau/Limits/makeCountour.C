#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

double crossSection(double stauM, double neutM);

void makeCountour(std::string directory="./Results/")
{
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetMarkerStyle(15);
  gStyle->SetMarkerSize(0.25);
  gStyle->SetTextFont(42);
  gStyle->SetMarkerColor(37);


  int minStauM = 50;
  int maxStauM = 500;
  int stepStauM = 10;

  int minNeutM = 0;
  int maxNeutM = 480;
  int stepNeutM = 10;

//  maxStauM = 300;
//  maxNeutM = 300;

  TCanvas c1("c1", "c1", 800, 600);
  TFile outFile("TStauStau_limit.root","RECREATE");

  TH2D histo2D("twodplot", "twodplot", (maxStauM-minStauM)/stepStauM + 1, minStauM - stepStauM/2, maxStauM + stepStauM/2, (maxNeutM-minNeutM)/stepNeutM + 1, minNeutM - stepNeutM/2, maxNeutM + stepNeutM/2);
  TH2D sigStrength("SignalStrength", "r", (maxStauM-minStauM)/stepStauM + 1, minStauM - stepStauM/2, maxStauM + stepStauM/2, (maxNeutM-minNeutM)/stepNeutM + 1, minNeutM - stepNeutM/2, maxNeutM + stepNeutM/2);

  for(int stauM = minStauM; stauM <= maxStauM; stauM += stepStauM)
  {
    for(int neutM = minNeutM; neutM <= maxNeutM; neutM += stepNeutM)
    {
      if(neutM == 0)
        neutM = 1;

      std::stringstream temp;
      temp << directory;
      temp << "higgsCombineS" << stauM << "-N" << neutM << ".Asymptotic.mH120.root";

      std::string fileName;
      temp >> fileName;
      std::cout << fileName << std::endl;

      TFile file(fileName.c_str(), "READ");

//      std::cout << "Opened file" << std::endl;

      if(file.IsOpen())
      {
//        std::cout << "Getting tree" << std::endl;
        TTree *treelimit = (TTree*)file.Get("limit");
        if(treelimit == NULL)
          continue;
//        std::cout << "Got tree" << std::endl;

        TH1D* obs = new TH1D("obs","",100,0,100);
        treelimit->Draw("limit>>obs", "quantileExpected==-1");
        TH1D* expectedm2 = new TH1D("expectedm2","",100,0,100);
        treelimit->Draw("limit>>expectedm2", "quantileExpected>0.02 && quantileExpected<0.03");
        TH1D* expectedm1 = new TH1D("expectedm1","",100,0,100);
        treelimit->Draw("limit>>expectedm1", "quantileExpected>0.15 && quantileExpected<0.16");
        TH1D* expected = new TH1D("expected","",100,0,100);
        treelimit->Draw("limit>>expected", "quantileExpected==0.5");
        TH1D* expectedp1 = new TH1D("expectedp1","",100,0,100);
        treelimit->Draw("limit>>expectedp1", "quantileExpected>0.83 && quantileExpected<0.84");
        TH1D* expectedp2 = new TH1D("expectedp2","",100,0,100);
        treelimit->Draw("limit>>expectedp2", "quantileExpected>0.97 && quantileExpected<0.98");// */

//        std::cout << "Drew everything" << std::endl;

        double limit = expected->GetMean();
        if(limit <= 0)
          limit = 2e8;
        else
        {
//          limit /= 30;
        }

        std::cout << stauM << ", " << neutM << ": " << limit << std::endl;

        histo2D.Fill(stauM, neutM, limit);
        sigStrength.Fill(stauM, neutM, limit/crossSection(stauM, neutM));

        file.Close();
      }
      else
      {
//        histo2D.Fill(stauM, neutM, 500);
      }

      if(neutM == 1)
        neutM = 0;
    }
  }

  outFile.cd();

  double contours[2];
  contours[0] = 2;
  contours[1] = 10000;

  int colors[2] = {2,4}; //red, blue,black
//  gStyle->SetPalette(2,colors);

//  histo2D.SetContour(2);
//  histo2D.SetContourLevel(0,1); //value for your first level
//  histo2D.SetContourLevel(1,1e8); //non-existing high level
//  histo2D.SetFillColor(2);
//  histo2D.Draw("cont3");

  histo2D.Write("twodplot");
  sigStrength.Write("signalStrength");

  TH2D *histo2D2=(TH2D*) histo2D.Clone("histo2Dclone");
//  histo2D2->SetContour(2);
//  histo2D2->SetContourLevel(0,-1e3); //non existing low level
//  histo2D2->SetContourLevel(1,1); //value for your second level
//  histo2D2->SetFillColor(4);
  c1.SetLogz();
  histo2D2->SetTitle("Upper Limit Cross Section;M_{Stau};M_{Neut};pb");
  histo2D2->Draw("colz");
  histo2D2->GetZaxis()->SetRangeUser(0.03, 1000);
//  histo2D2->Draw("cont3 same");

  gStyle->SetTextFont(42);
  c1.SaveAs("contour_1overmu.png");
  c1.SaveAs("contour_1overmu.C");

  sigStrength.Draw("colz");
  c1.SaveAs("signal_strength.png");
  c1.SaveAs("signal_strength.C");

  outFile.Write();
  outFile.Close();
}

double crossSection(double stauM, double neutM)
{
  double a = 0.2979;
  double b = 17.626;
  double c = 67.632;
  double d = 3.463;
  return a / (1 + std::pow((stauM - b) / c, d));
}
