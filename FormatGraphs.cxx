#include "GeneratorStudy.h"


/**
 *  main function
 * Arguments are <root file> <config file (optional)>
 */
int main(int argc, char **argv)
{
  gStyle->SetOptStat(0);

  string graphfilename="";

  if (argc < 2)
  {
    cout<<"Usage: "<<argv[0]<<" <ROOT file containing graphs> <summary title>"<<endl;
    return -1;
  }
  
  graphfilename=argv[1];
  if (! (graphfilename.find(".root") + 5 == graphfilename.length() && graphfilename.length() > 4))
  {
    cout<<"Not a root file "<<graphfilename<<endl;
    return -1;
  }
  
  string title=graphfilename;
  
  if (argc > 2)
    title=argv[2];
  
  TFile *graphFile = new TFile(graphfilename.c_str());
  if (graphFile->IsZombie())
  {
    cout<<"No valid ROOT file given. Filename was "<<graphfilename<<endl;
    return -1;
  }

  string plotfilename = graphfilename;
  plotfilename.replace(graphfilename.length()-5,5,".png");
  cout<<"Processing plots from "<<graphfilename<<" ("<<title<<") and writing to "<<plotfilename<<endl;
  return 0;
}
