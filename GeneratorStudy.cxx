#include "GeneratorStudy.h"
TRandom3 *trand;

/**
 *  main function
 * Arguments are <root file> <config file (optional)>
 */
int main(int argc, char **argv)
{
  trand = new TRandom3();
  gStyle->SetOptStat(0);
  
  Analyze(ND150);
  return 0;
}

void Analyze(ISOTOPE isotope)
{
  TH1D *smeared2nu = makeSmearedHistogram(isotope,true);
 // TH1D *smeared0nu = makeSmearedHistogram(isotope,false);

  return ;

}

TH1D * makeSmearedHistogram(ISOTOPE isotope, bool is2nu)
{
  string fName;
  if (is2nu)fName=FILES2NU[isotope]; else fName=FILES0NU[isotope];
  TFile *f = new TFile(fName.c_str());
  TTree *tree = (TTree*) f->Get(treeName.c_str());
  double maxEnergy=Qbb[isotope];
  TH1D *htrue = new TH1D("htrue",(ISOTOPE_LATEX[isotope]).c_str(),300,2.1,maxEnergy*1.1);
  TH1D *hsmeared = new TH1D("hsmeared",(ISOTOPE_LATEX[isotope]).c_str(),300,2.1,maxEnergy*1.1);
  
  vector<double> *electronEnergy = 0;
  tree->SetBranchAddress("trueparticle.kinenergy", &electronEnergy);
  // Loop the entries
  int nEntries = tree->GetEntries();

  for (int i=0;i<nEntries;i++)
  {
    tree->GetEntry(i);
    double totalEnergy=electronEnergy->at(0)+electronEnergy->at(1);
    htrue->Fill(totalEnergy); // Do we actually want to do anything with this?
    double smeared1 = Smear(electronEnergy->at(0), RESOLUTION_AT_1MeV);
    double smeared2 = Smear(electronEnergy->at(1), RESOLUTION_AT_1MeV);
    hsmeared->Fill(smeared1+smeared2);
  }
  TCanvas *c = new TCanvas (("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),("Smeared energy: "+ISOTOPE_NAME[isotope]).c_str(),900,600);
  htrue->Draw("HIST");
  hsmeared->SetLineColor(kRed);
  hsmeared->Draw("HIST SAME");
  
  // Add a legend
  TLegend* legend = new TLegend(0.75,0.8,0.9,0.9);
  legend->AddEntry(htrue, "True", "lep");
  legend->AddEntry(hsmeared,"Smeared", "fl");
  legend->Draw();
  
  // Save a PNG
  string title="energySmear0nu_"+ISOTOPE_NAME[isotope]+".png";
  if (is2nu)title="energySmear2nu_"+ISOTOPE_NAME[isotope]+".png";
  c->SaveAs(title.c_str());
  delete c;
  return hsmeared;
}

double Smear(double energy, double smearCoefficient) // coefficient is fractional smear at 1 MeV
{
  // sigma = k sqrt E so for 1% at 1 MeV, k is 0.01
  double sigma = smearCoefficient * TMath::Sqrt(energy);
  // Pick a random number from a Gaussian distribution with the mean as the true energy and the sigma as calculated
  double smeared= trand->Gaus(energy,sigma);
  return smeared;
  
}
