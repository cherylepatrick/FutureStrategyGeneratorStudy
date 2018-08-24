#include "GeneratorStudy.h"
TRandom3 *trand;

/**
 *  main function
 * Arguments are <root file> <config file (optional)>
 */
int main(int argc, char **argv)
{
  trand = new TRandom3();
  Analyze(ND150);
  return 0;
}

void Analyze(ISOTOPE isotope)
{
  // Open the file and extract the simulated data tree
  string fName=FILES[isotope];
  TFile *f = new TFile(fName.c_str());
  TTree *tree = (TTree*) f->Get(treeName.c_str());
  double maxEnergy=Qbb[isotope];
  TH1D *h2nubb = new TH1D("h2nubb",(ISOTOPE_LATEX[isotope]).c_str(),300,2.1,maxEnergy*1.1);
  TH1D *s2nubb = new TH1D("smeared2nubb",(ISOTOPE_LATEX[isotope]).c_str(),300,2.1,maxEnergy*1.1);
  
  vector<double> *electronEnergy = 0;
  tree->SetBranchAddress("trueparticle.kinenergy", &electronEnergy);
  // Loop the entries
  int nEntries = tree->GetEntries();
  //nEntries=20; // ######## TAKE THIS OUT
  for (int i=0;i<nEntries;i++)
  {
    tree->GetEntry(i);
    double totalEnergy=electronEnergy->at(0)+electronEnergy->at(1);
    h2nubb->Fill(totalEnergy);
    double smeared1 = Smear(electronEnergy->at(0), 0.01);
    double smeared2 = Smear(electronEnergy->at(1), 0.01);
    s2nubb->Fill(smeared1+smeared2);
  }

  TCanvas *c = new TCanvas (("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),600,900);
  h2nubb->Draw("HIST");
  s2nubb->SetLineColor(kRed);
  s2nubb->Draw("HIST SAME");

  c->SaveAs(("totalEnergy_"+ISOTOPE_NAME[isotope]+".png").c_str());
  delete c;
  return ;

}

double Smear(double energy, double smearCoefficient) // coefficient is fractional smear at 1 MeV
{
  // sigma = k sqrt E so for 1% at 1 MeV, k is 0.01
  double sigma = smearCoefficient * TMath::Sqrt(energy);
  // Pick a random number from a Gaussian distribution with the mean as the true energy and the sigma as calculated
  double smeared= trand->Gaus(energy,sigma);
  return smeared;
  
}
