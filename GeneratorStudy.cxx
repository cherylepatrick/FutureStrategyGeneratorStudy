#include "GeneratorStudy.h"


/**
 *  main function
 * Arguments are <root file> <config file (optional)>
 */
int main(int argc, char **argv)
{
  analyze(ND150);
  return 0;
}

void analyze(ISOTOPE isotope)
{
  // Open the file and extract the simulated data tree
  string fName=FILES[isotope];
  TFile *f = new TFile(fName.c_str());
  TTree *tree = (TTree*) f->Get(treeName.c_str());
  double maxEnergy=Qbb[isotope];
  TH1D *h2nubb = new TH1D("h2nubb",(ISOTOPE_LATEX[isotope]).c_str(),300,2,maxEnergy*1.1);
  
  vector<double> *electronEnergy = 0;
  tree->SetBranchAddress("trueparticle.kinenergy", &electronEnergy);
  // Loop the entries
  int nEntries = tree->GetEntries();
 // nEntries=50000; // ######## TAKE THIS OUT
  for (int i=0;i<nEntries;i++)
  {
    tree->GetEntry(i);
    double totalEnergy=electronEnergy->at(0)+electronEnergy->at(1);
    h2nubb->Fill(totalEnergy);
  }

  TCanvas *c = new TCanvas (("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),600,900);
  h2nubb->Draw("HIST");
  c->SaveAs(("totalEnergy_"+ISOTOPE_NAME[isotope]+".png").c_str());
  delete c;
  return ;

}
