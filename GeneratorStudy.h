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
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TDictionary.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TLimit.h"
#include "TConfidenceLevel.h"
#include "TLimitDataSource.h"
#include "TError.h"
#include "TGraph.h"

using namespace std;

string treeName="SimData";


enum ISOTOPE  {SE82, ND150};

// Everything for 82Se and 150Nd
string ISOTOPE_LATEX[2] = {"^{82}Se","^{150}Nd"};
string ISOTOPE_NAME[2] = {"Se82","Nd150"};
double Qbb[2] = {2.99,3.368};
//string FILES2NU[2]={"/Users/cpatrick/uclnemo3/generatorstudy/se82/2nubb/sim_2nubb_trimmed_10M.root","/Users/cpatrick/uclnemo3/generatorstudy/nd150/run_1/Nd150_2nubb_1E6_flsim_1_trimmed.root"};
//string FILES0NU[2]={"/Users/cpatrick/uclnemo3/generatorstudy/se82/0nubb/Se82_0nubb_1E8_flsim_1.root"," /Users/cpatrick/uclnemo3/generatorstudy/nd150/nd150_0nu/1e5/run_1/Nd150_0nubb_1E6_flsim_1.root"};
string SMEARED_HISTO_FILE[2]={"smeared_hists_Se82.root","smeared_hists_Nd150.root"};
int TOTAL_2NU_EVENTS[2]={1000000,1000000};
//double FRAC_OVER_2MEV[2]={1,0.0982};
int ATOMIC_MASS[2]={82,150}; //Selenium 82, Neodymium 150
double HALFLIFE2NU[2]={10.07e19,9.1e18}; // 2nubb halflife in years
int TLIMIT_EXPERIMENTS=50000; // Number of pesudoexperiments to run for limit calculation
double DESIRED_CONFIDENCE=0.003; // 0.003 = 3 sigma (99.7% confident)
// How long do we need to run to reach this 0nubb halflife sensitivity (years)?
// The goal is something that corresponds to 1e28 in the big future experiments (which use different isotopes)
// The numbers are based on the average of the calculations from Laurent's list in docdb 4680 tables 5-8
double SENSITIVITY_LEGEND_Se = 2.75e27;// if LEGEND saw 1e28 year halflife in 76Ge
double SENSITIVITY_LEGEND_Nd = 1.81e27;// if LEGEND saw 1e28 year halflife in 76Ge
double SENSITIVITY_NEXO_Se = 6.70e27;// if nEXO saw 1e28 year halflife in 136Xe
double SENSITIVITY_NEXO_Nd = 3.72e27;// if nEXO saw 1e28 year halflife in 136Xe
double AVOGADRO = 6.022140e23;


int main(int argc, char **argv);

double ExpectedLimitSigEvts(double ConfidenceLevel, TH1D* h_signal, TH1D* h_background, TH1D* h_data );
TGraph *ScaledClone(TGraph *graph, double scale);
TGraph* MakeExposureGraph(string experimentText,ISOTOPE isotope,double desiredHalflife);
double GetExposure(TH1D *hist2nu, TH1D *hist0nu, ISOTOPE isotope, double desiredHalflife);
double Get2nuEventsForExposure(double exposure, ISOTOPE isotope);
double Get0nuEventsForExposure(double exposure, ISOTOPE isotope, double desiredHalflife);
double GetExposureFrom2nuEvents(double events, ISOTOPE isotope);
int GetNextSmearing(std::vector<double> smearings, double previousSmearing);
double GetExposureFrom0nuEvents(double events, ISOTOPE isotope, double halflife);
