#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"

int makeDotC(std::string file)
{
  int lastindex = file.find_last_of(".");
  string rawname = file.substr(0, lastindex);

  std::cout << rawname << std::endl;


  TFile rootFile(file.c_str());

  TCanvas* c1 = static_cast<TCanvas*>(rootFile.Get("c1"));

  c1->Draw();
  c1->SaveAs((rawname+".C").c_str());
//  c1->SaveAs("test.C");

  return 0;
}
