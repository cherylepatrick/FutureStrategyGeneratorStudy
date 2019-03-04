
#include "SmearHistogram.h"

// All this does is take a smeared dataset, histogram it, and write it to a file
// You can choose the binning.

/**
 *  main function
 *  Arguments are <Se or Nd> <smear percent> <2 or 0 (nu)> <input root file>  <output root file> <number of bins> <min energy in MeV> <max energy in MeV>
 */
string smearedtree="SmearedData";
int main(int argc, char **argv)
{

  gStyle->SetOptStat(0);
  string helptext="<Se or Nd> <smear percent> <2 or 0 (nu)> <input ROOT file>  <output root file> <number of bins> <min energy in MeV> <max energy in MeV>";
  
  double smearing=0;
  ISOTOPE isotope;
  double minEnergy=2.1;
  double maxEnergy=4.;
  double nBins=100;
  
  if (argc>1)
  {
    string arg1=argv[1];
    if ((std::toupper(arg1[0])=='S') && (std::toupper(arg1[1])=='E'))isotope=SE82;
    else if ((std::toupper(arg1[0])=='N') && (std::toupper(arg1[1])=='D'))isotope=ND150;
    else
    {
      cout<<"Not a valid isotope, Se or Nd"<<endl;
      cout<<helptext<<endl;
      return 1;
    }
    
  }
  else
  {
    cout<<helptext<<endl;
    return 1;
  }
  
  if (argc>2)
  {
    try{
      smearing = atof(argv[2]);
    }
    catch(const std::exception&)
    {
      cout<<"Smear percentage is not a number"<<endl;
      cout<<helptext<<endl;
      return 1;
    }
  }
  else
  {
    cout<<"Smear percentage needed e.g. 3.5 "<<endl;
    cout<<helptext<<endl;
    return 1;
  }
  
  bool is2nu=true;
  
  if (argc>3)
  {
    string arg3=argv[3];
    if (arg3=="0") is2nu=false;
    else if (arg3=="2") is2nu=true;
    else
    {
      cout<<"2 or 0 needed for number of neutrinos"<<endl;
      cout<<helptext<<endl;
      return 1;
    }
  }
  else
  {
    cout<<"Number of neutrinos needed - 2 or 0"<<endl;
    cout<<helptext<<endl;
    return 1;
  }
  
  cout<<"Smearing "<<(is2nu?"2":"0")<<"nubb histograms for "<<ISOTOPE_NAME[isotope]<<" by "<<smearing<<"%"<<endl;
  
  string inFileName="";
  bool inFileOK=false;
  TFile *inFile;
  TTree *tree ;
  if (argc>4)
  {
    inFileOK=true; // Until proven otherwise!
    string arg4=argv[4];
    if (arg4.find(".root") + 5 == arg4.length() && arg4.length() > 4)
    {
      inFileName=arg4;
      try{
        inFile = new TFile(inFileName.c_str());
        if (inFile->IsZombie())
        {
          inFileOK=false;
        }
        else{
          try{
            if (! inFile->GetListOfKeys()->Contains(smearedtree.c_str()))
                  {
                    cout<<"No tree named "<<smearedtree<<" in "<<inFileName<<endl;
                    inFileOK=false;
                  }
                else
                {
                  tree = (TTree*) inFile->Get(smearedtree.c_str());
                  if (tree->IsZombie())
                  {
                    inFileOK=false;
                    cout<<"No tree named "<<smearedtree<<" in "<<inFileName<<endl;
                  }
                }
          }
          catch(const std::exception&)
          {
            cout<<"No tree named "<<smearedtree<<" in "<<inFileName<<endl;
            inFileOK=false;
          }
        }
      }
      catch (const std::exception&)
      {
        cout<<"Could not open "<<inFileName<<endl;
        inFileOK=false;
      }

    }
    else
    {
      inFileOK=false;
    }

  }
  
  if (!inFileOK)
  {
    cout<<"Invalid input ROOT file "<<inFileName<<endl;
    cout<<helptext<<endl;
    return 1;
  }
  cout<<"inFile "<<inFileName<<" is ok? "<<inFileOK<<endl;
}

