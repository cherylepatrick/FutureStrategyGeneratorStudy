
#include "SmearHistogram.h"

// All this does is take a smeared dataset, histogram it, and write it to a file
// You can choose the binning.

/**
 *  main function
 *  Arguments are <Se or Nd> <smear percent> <2 or 0 (nu)> <input root file>  <output root file> <number of bins> <min energy in MeV> <max energy in MeV>
 */
int main(int argc, char **argv)
{

  gStyle->SetOptStat(0);
  string helptext="<Se or Nd> <smear percent> <2 or 0 (nu)> <input root file>  <output root file> <number of bins> <min energy in MeV> <max energy in MeV>";
  
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
      cout<<helptext<<endl;
      return 1;
    }
  }
  else
  {
    cout<<helptext<<endl;
    return 1;
  }
  
  cout<<"Smearing histograms for "<<ISOTOPE_NAME[isotope]<<" by "<<smearing<<"%"<<endl;
}

