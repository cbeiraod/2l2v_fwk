#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1D.h"

//#include "FWCore/FWLite/interface/AutoLibraryLoader.h"


int makeStauXSecLUT()
{
  gSystem->Load("libFWCoreFWLite");
//  AutoLibraryLoader::enable();
  
  TFile file("../../../data/TStauStau/StauCrossSections.root", "RECREATE");
  
  TH1D* xsec      = new TH1D("xsec", "crossSection;M_{Stau};#sigma", 500/5 - 1, 5, 500);
  TH1D* xsec_UP   = static_cast<TH1D*>(xsec->Clone("xsec_UP"));
  TH1D* xsec_DOWN = static_cast<TH1D*>(xsec->Clone("xsec_DOWN"));
  std::vector<double> xVals, yVals;
  
  double SF = 1;
  {
    double a = 0.2979;
    double b = 17.626;
    double c = 67.632;
    double d = 3.463;
    double val = (a / (1 + std::pow((100 - b) / c, d)));
    SF = 0.126 / val;
  }
  
  for(double mass = 20; mass <= 500; mass += 5)
  {
    // Cross sections from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVsleptonsleptonCMS
    // For MStau < 100, points taken from  http://arxiv.org/abs/1204.2379
    //   Mstau == 100 => 0.1
    //   Mstau == 125 => 0.05
    //   Mstau == 145 => 0.03
    //   Mstau == 195 => 0.01
    //   Mstau == 240 => 0.005
    //   Mstau == 275 => 0.003
    //   Mstau == 300 => 0.002
    //   Mstau == 360 => 0.001
    //   Mstau == 425 => 0.0005

    if(mass < 100)
    {
      double a = 0.2979;
      double b = 17.626;
      double c = 67.632;
      double d = 3.463;
      double val = (a / (1 + std::pow((mass - b) / c, d)));
      val *= SF;
      
      xsec->Fill(mass, val);
      xsec_UP->Fill(mass, val*1.03);
      xsec_DOWN->Fill(mass, val*0.97);
    }
    else
    {
      if(mass < 105)
      {
        xVals.push_back(mass);
        yVals.push_back(0.126);
        xsec->Fill(mass, 0.126);
        xsec_UP->Fill(mass, 0.126 + 0.003);
        xsec_DOWN->Fill(mass, 0.126 - 0.003);
      }
      else if(mass < 110)
      {
        xVals.push_back(mass);
        yVals.push_back(0.1050811);
        xsec->Fill(mass, 0.1050811);
        xsec_UP->Fill(mass, 0.1050811 + 0.002081139);
        xsec_DOWN->Fill(mass, 0.1050811 - 0.002081139);
      }
      else if(mass < 115)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0879);
        xsec->Fill(mass, 0.0879);
        xsec_UP->Fill(mass, 0.0879 + 0.002);
        xsec_DOWN->Fill(mass, 0.0879 - 0.002);
      }
      else if(mass < 120)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0742);
        xsec->Fill(mass, 0.0742);
        xsec_UP->Fill(mass, 0.0742 + 0.001803996);
        xsec_DOWN->Fill(mass, 0.0742 - 0.001803996);
      }
      else if(mass < 125)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0631);
        xsec->Fill(mass, 0.0631);
        xsec_UP->Fill(mass, 0.0631 + 0.0015);
        xsec_DOWN->Fill(mass, 0.0631 - 0.0015);
      }
      else if(mass < 130)
      {
        xVals.push_back(mass);
        yVals.push_back(0.054);
        xsec->Fill(mass, 0.054);
        xsec_UP->Fill(mass, 0.054 + 0.0013);
        xsec_DOWN->Fill(mass, 0.054 - 0.0013);
      }
      else if(mass < 135)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0464);
        xsec->Fill(mass, 0.0464);
        xsec_UP->Fill(mass, 0.0464 + 0.0012);
        xsec_DOWN->Fill(mass, 0.0464 - 0.0012);
      }
      else if(mass < 140)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0401);
        xsec->Fill(mass, 0.0401);
        xsec_UP->Fill(mass, 0.0401 + 0.001);
        xsec_DOWN->Fill(mass, 0.0401 - 0.001);
      }
      else if(mass < 145)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0349);
        xsec->Fill(mass, 0.0349);
        xsec_UP->Fill(mass, 0.0349 + 0.0008089499);
        xsec_DOWN->Fill(mass, 0.0349 - 0.0008089499);
      }
      else if(mass < 150)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0304);
        xsec->Fill(mass, 0.0304);
        xsec_UP->Fill(mass, 0.0304 + 0.0007);
        xsec_DOWN->Fill(mass, 0.0304 - 0.0007);
      }
      else if(mass < 155)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0266);
        xsec->Fill(mass, 0.0266);
        xsec_UP->Fill(mass, 0.0266 + 0.0007);
        xsec_DOWN->Fill(mass, 0.0266 - 0.0007);
      }
      else if(mass < 160)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0234);
        xsec->Fill(mass, 0.0234);
        xsec_UP->Fill(mass, 0.0234 + 0.0005141984);
        xsec_DOWN->Fill(mass, 0.0234 - 0.0005141984);
      }
      else if(mass < 165)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0206);
        xsec->Fill(mass, 0.0206);
        xsec_UP->Fill(mass, 0.0206 + 0.0005384877);
        xsec_DOWN->Fill(mass, 0.0206 - 0.0005384877);
      }
      else if(mass < 170)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0182);
        xsec->Fill(mass, 0.0182);
        xsec_UP->Fill(mass, 0.0182 + 0.0005391333);
        xsec_DOWN->Fill(mass, 0.0182 - 0.0005391333);
      }
      else if(mass < 175)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0162);
        xsec->Fill(mass, 0.0162);
        xsec_UP->Fill(mass, 0.0162 + 0.0005247970);
        xsec_DOWN->Fill(mass, 0.0162 - 0.0005247970);
      }
      else if(mass < 180)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0144);
        xsec->Fill(mass, 0.0144);
        xsec_UP->Fill(mass, 0.0144 + 0.0004350736);
        xsec_DOWN->Fill(mass, 0.0144 - 0.0004350736);
      }
      else if(mass < 185)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0128);
        xsec->Fill(mass, 0.0128);
        xsec_UP->Fill(mass, 0.0128 + 0.0004392250);
        xsec_DOWN->Fill(mass, 0.0128 - 0.0004392250);
      }
      else if(mass < 190)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0115);
        xsec->Fill(mass, 0.0115);
        xsec_UP->Fill(mass, 0.0115 + 0.0003610438);
        xsec_DOWN->Fill(mass, 0.0115 - 0.0003610438);
      }
      else if(mass < 195)
      {
        xVals.push_back(mass);
        yVals.push_back(0.0103);
        xsec->Fill(mass, 0.0103);
        xsec_UP->Fill(mass, 0.0103 + 0.0002339203);
        xsec_DOWN->Fill(mass, 0.0103 - 0.0002339203);
      }
      else if(mass < 200)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00922);
        xsec->Fill(mass, 0.00922);
        xsec_UP->Fill(mass, 0.00922 + 0.0002352341);
        xsec_DOWN->Fill(mass, 0.00922 - 0.0002352341);
      }
      else if(mass < 205)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00830);
        xsec->Fill(mass, 0.00830);
        xsec_UP->Fill(mass, 0.00830 + 0.0002034879);
        xsec_DOWN->Fill(mass, 0.00830 - 0.0002034879);
      }
      else if(mass < 210)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00749);
        xsec->Fill(mass, 0.00749);
        xsec_UP->Fill(mass, 0.00749 + 0.0001848360);
        xsec_DOWN->Fill(mass, 0.00749 - 0.0001848360);
      }
      else if(mass < 215)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00676);
        xsec->Fill(mass, 0.00676);
        xsec_UP->Fill(mass, 0.00676 + 0.0001717110);
        xsec_DOWN->Fill(mass, 0.00676 - 0.0001717110);
      }
      else if(mass < 220)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00612);
        xsec->Fill(mass, 0.00612);
        xsec_UP->Fill(mass, 0.00612 + 0.0001519336);
        xsec_DOWN->Fill(mass, 0.00612 - 0.0001519336);
      }
      else if(mass < 225)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00555);
        xsec->Fill(mass, 0.00555);
        xsec_UP->Fill(mass, 0.00555 + 0.0001420703);
        xsec_DOWN->Fill(mass, 0.00555 - 0.0001420703);
      }
      else if(mass < 230)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00504);
        xsec->Fill(mass, 0.00504);
        xsec_UP->Fill(mass, 0.00504 + 0.0001224106);
        xsec_DOWN->Fill(mass, 0.00504 - 0.0001224106);
      }
      else if(mass < 235)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00458);
        xsec->Fill(mass, 0.00458);
        xsec_UP->Fill(mass, 0.00458 + 0.0001224182);
        xsec_DOWN->Fill(mass, 0.00458 - 0.0001224182);
      }
      else if(mass < 240)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00418);
        xsec->Fill(mass, 0.00418);
        xsec_UP->Fill(mass, 0.00418 + 0.0001030006);
        xsec_DOWN->Fill(mass, 0.00418 - 0.0001030006);
      }
      else if(mass < 245)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00381);
        xsec->Fill(mass, 0.00381);
        xsec_UP->Fill(mass, 0.00381 + 0.00009315606);
        xsec_DOWN->Fill(mass, 0.00381 - 0.00009315606);
      }
      else if(mass < 250)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00348);
        xsec->Fill(mass, 0.00348);
        xsec_UP->Fill(mass, 0.00348 + 0.00008352312);
        xsec_DOWN->Fill(mass, 0.00348 - 0.00008352312);
      }
      else if(mass < 255)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00318);
        xsec->Fill(mass, 0.00318);
        xsec_UP->Fill(mass, 0.00318 + 0.00008089570);
        xsec_DOWN->Fill(mass, 0.00318 - 0.00008089570);
      }
      else if(mass < 260)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00291);
        xsec->Fill(mass, 0.00291);
        xsec_UP->Fill(mass, 0.00291 + 0.00008102447);
        xsec_DOWN->Fill(mass, 0.00291 - 0.00008102447);
      }
      else if(mass < 265)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00267);
        xsec->Fill(mass, 0.00267);
        xsec_UP->Fill(mass, 0.00267 + 0.00007103856);
        xsec_DOWN->Fill(mass, 0.00267 - 0.00007103856);
      }
      else if(mass < 270)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00245);
        xsec->Fill(mass, 0.00245);
        xsec_UP->Fill(mass, 0.00245 + 0.00006119814);
        xsec_DOWN->Fill(mass, 0.00245 - 0.00006119814);
      }
      else if(mass < 275)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00225);
        xsec->Fill(mass, 0.00225);
        xsec_UP->Fill(mass, 0.00225 + 0.00006119922);
        xsec_DOWN->Fill(mass, 0.00225 - 0.00006119922);
      }
      else if(mass < 280)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00207);
        xsec->Fill(mass, 0.00207);
        xsec_UP->Fill(mass, 0.00207 + 0.00005143164);
        xsec_DOWN->Fill(mass, 0.00207 - 0.00005143164);
      }
      else if(mass < 285)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00191);
        xsec->Fill(mass, 0.00191);
        xsec_UP->Fill(mass, 0.00191 + 0.00005143055);
        xsec_DOWN->Fill(mass, 0.00191 - 0.00005143055);
      }
      else if(mass < 290)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00176);
        xsec->Fill(mass, 0.00176);
        xsec_UP->Fill(mass, 0.00176 + 0.00004178052);
        xsec_DOWN->Fill(mass, 0.00176 - 0.00004178052);
      }
      else if(mass < 295)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00162);
        xsec->Fill(mass, 0.00162);
        xsec_UP->Fill(mass, 0.00162 + 0.00004002940);
        xsec_DOWN->Fill(mass, 0.00162 - 0.00004002940);
      }
      else if(mass < 300)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00150);
        xsec->Fill(mass, 0.00150);
        xsec_UP->Fill(mass, 0.00150 + 0.00004181893);
        xsec_DOWN->Fill(mass, 0.00150 - 0.00004181893);
      }
      else if(mass < 305)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00139);
        xsec->Fill(mass, 0.00139);
        xsec_UP->Fill(mass, 0.00139 + 0.00005588316);
        xsec_DOWN->Fill(mass, 0.00139 - 0.00005588316);
      }
      else if(mass < 310)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00128);
        xsec->Fill(mass, 0.00128);
        xsec_UP->Fill(mass, 0.00128 + 0.00003005128);
        xsec_DOWN->Fill(mass, 0.00128 - 0.00003005128);
      }
      else if(mass < 315)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00119);
        xsec->Fill(mass, 0.00119);
        xsec_UP->Fill(mass, 0.00119 + 0.00003235425);
        xsec_DOWN->Fill(mass, 0.00119 - 0.00003235425);
      }
      else if(mass < 320)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00110);
        xsec->Fill(mass, 0.00110);
        xsec_UP->Fill(mass, 0.00110 + 0.00003237770);
        xsec_DOWN->Fill(mass, 0.00110 - 0.00003237770);
      }
      else if(mass < 325)
      {
        xVals.push_back(mass);
        yVals.push_back(0.00102);
        xsec->Fill(mass, 0.00102);
        xsec_UP->Fill(mass, 0.00102 + 0.00003252475);
        xsec_DOWN->Fill(mass, 0.00102 - 0.00003252475);
      }
      else if(mass < 330)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000948);
        xsec->Fill(mass, 0.000948);
        xsec_UP->Fill(mass, 0.000948 + 0.00002449507);
        xsec_DOWN->Fill(mass, 0.000948 - 0.00002449507);
      }
      else if(mass < 335)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000880);
        xsec->Fill(mass, 0.000880);
        xsec_UP->Fill(mass, 0.000880 + 0.00002340214);
        xsec_DOWN->Fill(mass, 0.000880 - 0.00002340214);
      }
      else if(mass < 340)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000818);
        xsec->Fill(mass, 0.000818);
        xsec_UP->Fill(mass, 0.000818 + 0.00002252870);
        xsec_DOWN->Fill(mass, 0.000818 - 0.00002252870);
      }
      else if(mass < 345)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000761);
        xsec->Fill(mass, 0.000761);
        xsec_UP->Fill(mass, 0.000761 + 0.00002891792);
        xsec_DOWN->Fill(mass, 0.000761 - 0.00002891792);
      }
      else if(mass < 350)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000708);
        xsec->Fill(mass, 0.000708);
        xsec_UP->Fill(mass, 0.000708 + 0.00001862963);
        xsec_DOWN->Fill(mass, 0.000708 - 0.00001862963);
      }
      else if(mass < 355)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000660);
        xsec->Fill(mass, 0.000660);
        xsec_UP->Fill(mass, 0.000660 + 0.00001708859);
        xsec_DOWN->Fill(mass, 0.000660 - 0.00001708859);
      }
      else if(mass < 360)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000615);
        xsec->Fill(mass, 0.000615);
        xsec_UP->Fill(mass, 0.000615 + 0.00001615603);
        xsec_DOWN->Fill(mass, 0.000615 - 0.00001615603);
      }
      else if(mass < 365)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000573);
        xsec->Fill(mass, 0.000573);
        xsec_UP->Fill(mass, 0.000573 + 0.00001581577);
        xsec_DOWN->Fill(mass, 0.000573 - 0.00001581577);
      }
      else if(mass < 370)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000535);
        xsec->Fill(mass, 0.000535);
        xsec_UP->Fill(mass, 0.000535 + 0.00001386849);
        xsec_DOWN->Fill(mass, 0.000535 - 0.00001386849);
      }
      else if(mass < 375)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000499);
        xsec->Fill(mass, 0.000499);
        xsec_UP->Fill(mass, 0.000499 + 0.00001349490);
        xsec_DOWN->Fill(mass, 0.000499 - 0.00001349490);
      }
      else if(mass < 380)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000466);
        xsec->Fill(mass, 0.000466);
        xsec_UP->Fill(mass, 0.000466 + 0.00001253196);
        xsec_DOWN->Fill(mass, 0.000466 - 0.00001253196);
      }
      else if(mass < 385)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000436);
        xsec->Fill(mass, 0.000436);
        xsec_UP->Fill(mass, 0.000436 + 0.00001157656);
        xsec_DOWN->Fill(mass, 0.000436 - 0.00001157656);
      }
      else if(mass < 390)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000408);
        xsec->Fill(mass, 0.000408);
        xsec_UP->Fill(mass, 0.000408 + 0.00001109518);
        xsec_DOWN->Fill(mass, 0.000408 - 0.00001109518);
      }
      else if(mass < 395)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000381);
        xsec->Fill(mass, 0.000381);
        xsec_UP->Fill(mass, 0.000381 + 0.00001063177);
        xsec_DOWN->Fill(mass, 0.000381 - 0.00001063177);
      }
      else if(mass < 400)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000357);
        xsec->Fill(mass, 0.000357);
        xsec_UP->Fill(mass, 0.000357 + 0.000009697550);
        xsec_DOWN->Fill(mass, 0.000357 - 0.000009697550);
      }
      else if(mass < 405)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000334);
        xsec->Fill(mass, 0.000334);
        xsec_UP->Fill(mass, 0.000334 + 0.000009321515);
        xsec_DOWN->Fill(mass, 0.000334 - 0.000009321515);
      }
      else if(mass < 410)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000313);
        xsec->Fill(mass, 0.000313);
        xsec_UP->Fill(mass, 0.000313 + 0.000008381153);
        xsec_DOWN->Fill(mass, 0.000313 - 0.000008381153);
      }
      else if(mass < 415)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000294);
        xsec->Fill(mass, 0.000294);
        xsec_UP->Fill(mass, 0.000294 + 0.000007871547);
        xsec_DOWN->Fill(mass, 0.000294 - 0.000007871547);
      }
      else if(mass < 420)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000275);
        xsec->Fill(mass, 0.000275);
        xsec_UP->Fill(mass, 0.000275 + 0.000007400067);
        xsec_DOWN->Fill(mass, 0.000275 - 0.000007400067);
      }
      else if(mass < 425)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000258);
        xsec->Fill(mass, 0.000258);
        xsec_UP->Fill(mass, 0.000258 + 0.000007400054);
        xsec_DOWN->Fill(mass, 0.000258 - 0.000007400054);
      }
      else if(mass < 430)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000242);
        xsec->Fill(mass, 0.000242);
        xsec_UP->Fill(mass, 0.000242 + 0.000007400057);
        xsec_DOWN->Fill(mass, 0.000242 - 0.000007400057);
      }
      else if(mass < 435)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000228);
        xsec->Fill(mass, 0.000228);
        xsec_UP->Fill(mass, 0.000228 + 0.000005546330);
        xsec_DOWN->Fill(mass, 0.000228 - 0.000005546330);
      }
      else if(mass < 440)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000214);
        xsec->Fill(mass, 0.000214);
        xsec_UP->Fill(mass, 0.000214 + 0.000005546401);
        xsec_DOWN->Fill(mass, 0.000214 - 0.000005546401);
      }
      else if(mass < 445)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000201);
        xsec->Fill(mass, 0.000201);
        xsec_UP->Fill(mass, 0.000201 + 0.000005612918);
        xsec_DOWN->Fill(mass, 0.000201 - 0.000005612918);
      }
      else if(mass < 450)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000189);
        xsec->Fill(mass, 0.000189);
        xsec_UP->Fill(mass, 0.000189 + 0.000006973767);
        xsec_DOWN->Fill(mass, 0.000189 - 0.000006973767);
      }
      else if(mass < 455)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000178);
        xsec->Fill(mass, 0.000178);
        xsec_UP->Fill(mass, 0.000178 + 0.000004671806);
        xsec_DOWN->Fill(mass, 0.000178 - 0.000004671806);
      }
      else if(mass < 460)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000167);
        xsec->Fill(mass, 0.000167);
        xsec_UP->Fill(mass, 0.000167 + 0.000004181925);
        xsec_DOWN->Fill(mass, 0.000167 - 0.000004181925);
      }
      else if(mass < 465)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000157);
        xsec->Fill(mass, 0.000157);
        xsec_UP->Fill(mass, 0.000157 + 0.000004250135);
        xsec_DOWN->Fill(mass, 0.000157 - 0.000004250135);
      }
      else if(mass < 470)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000148);
        xsec->Fill(mass, 0.000148);
        xsec_UP->Fill(mass, 0.000148 + 0.000004224096);
        xsec_DOWN->Fill(mass, 0.000148 - 0.000004224096);
      }
      else if(mass < 475)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000139);
        xsec->Fill(mass, 0.000139);
        xsec_UP->Fill(mass, 0.000139 + 0.000004183248);
        xsec_DOWN->Fill(mass, 0.000139 - 0.000004183248);
      }
      else if(mass < 480)
      {
        xVals.push_back(mass);
        yVals.push_back(0.000131);
        xsec->Fill(mass, 0.000131);
        xsec_UP->Fill(mass, 0.000131 + 0.000004178115);
        xsec_DOWN->Fill(mass, 0.000131 - 0.000004178115);
      }
      else
      {
        xVals.push_back(mass);
        yVals.push_back(0.000124);
        xsec->Fill(mass, 0.000124);
        xsec_UP->Fill(mass, 0.000124 + 0.000003233469);
        xsec_DOWN->Fill(mass, 0.000124 - 0.000003233469);
      }
    }
  }
  
  xsec->Write();
  xsec_UP->Write();
  xsec_DOWN->Write();
  
//  file.Close();

  return 0;
}
