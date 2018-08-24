// Standard Library
#include <iostream>
//#include <fstream>
//#include "boost/algorithm/string.hpp"
//#include "boost/filesystem.hpp"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <cstdio>
//#include <memory>
//#include <stdexcept>
#include <string>
//#include <array>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
using namespace std;

/**
 *  main function
 * Arguments are <root file> <config file (optional)>
 */
int main(int argc, char **argv)
{
  string fName="/Users/cpatrick/uclnemo3/generatorstudy/nd150/run_1/Nd150_2nubb_1E6_flsim_1.root";
  TFile *fIn = new TFile(fName.c_str());
  TTree *tree = (TTree*) fIn->Get("SimData");
  string outname=fName.substr(0,fName.length()-5)+"_trimmed.root";
  cout<<outname<<endl;
  TFile *fOut = new TFile(outname.c_str(),"RECREATE");
  TTree *trimmed = tree->CopyTree("trueparticle.kinenergy[0]+trueparticle.kinenergy[1] > 2");
  cout<<"Fraction of entries kept: "<<(double)trimmed->GetEntries()/(double)tree->GetEntries()<<endl;
  trimmed->Write("",TObject::kOverwrite);
  fOut->Close();
}
