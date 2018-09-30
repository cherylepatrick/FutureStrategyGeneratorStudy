// Standard Library
#include <iostream>
#include <fstream>
//#include "boost/algorithm/string.hpp"
//#include "boost/filesystem.hpp"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstdio>
#include <string>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMath.h"

using namespace std;

string treeName="SimData";


enum ISOTOPE  {SE82, ND150};

// Everything for 82Se and 150Nd
string ISOTOPE_LATEX[2] = {"^{82}Se","^{150}Nd"};
string ISOTOPE_NAME[2] = {"Se82","Nd150"};
double Qbb[2] = {2.99,3.368};


string FILES2NU[2]={"/unix/nemo3/users/cpatrick/generatorstudy/se82/2nubb/sim_2nubb_trimmed_10M.root","/unix/nemo3/users/cpatrick/generatorstudy/nd150/sim_2nubb_150nd_trimmed_5M.root"};
string FILES0NU[2]={"/unix/nemo3/users/cpatrick/generatorstudy/se82/0nubb/Se82_0nubb_1E8_flsim_1.root","/unix/nemo3/users/cpatrick/generatorstudy/nd150/nd150_0nu/1e6/run_1/Nd150_2nubb_1E6_flsim_1.root"}; // name is misleading but it is actually 0nubb


double Smear(double energy, double smearCoefficient);
TH1D * makeSmearedHistogram(ISOTOPE isotope, bool is2nu, double resolutionAt1MeV);
