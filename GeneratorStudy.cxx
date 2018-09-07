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
  
  Analyze(ND150, 0.01); // 1% energy resolution at 1MeV
  Analyze(SE82, 0.01); // 1% energy resolution at 1MeV
  return 0;
}

void Analyze(ISOTOPE isotope, double resolutionAt1MeV)
{
  TH1D *smeared2nu = makeSmearedHistogram(isotope,true,resolutionAt1MeV);
  TH1D *smeared0nu = makeSmearedHistogram(isotope,false, resolutionAt1MeV);
  
  // Plot the two together
  TCanvas *c = new TCanvas (("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),("Energy: "+ISOTOPE_NAME[isotope]).c_str(),900,600);
  
  smeared2nu->SetLineColor(kBlack);
  smeared2nu->Draw("HIST");
  smeared0nu->Scale(smeared2nu->Integral()/smeared0nu->Integral() * 0.1); // Kind of an arbitrary scale for the 0nu sample, just for show
  smeared0nu->SetLineColor(kRed);
  smeared0nu->Draw("HIST SAME");
  
  // Add a legend
  TLegend* legend = new TLegend(0.75,0.8,0.9,0.9);
  legend->AddEntry(smeared2nu, "2#nu#beta#beta", "l");
  legend->AddEntry(smeared0nu,"0#nu#beta#beta", "l");
  legend->Draw();
  
  // Save a PNG
  string title="smearedComparison_"+ISOTOPE_NAME[isotope]+".png";
  c->SaveAs(title.c_str());
  delete c;
  return ;

}

TH1D * FAKESmearedHistogram(ISOTOPE isotope, double resolutionAt1MeV) // Just put this here til we get the 0nu simulation
{
  
  TH1D *hfake = new TH1D("hfake",(ISOTOPE_LATEX[isotope]).c_str(),300,2.1,Qbb[isotope]*1.1);
  for (int i=0;i<10000;i++)
  {
    double fakeSmeared=Smear(Qbb[isotope],resolutionAt1MeV);
    hfake->Fill(fakeSmeared);
  }
  
  TCanvas *c = new TCanvas (("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),("Smeared energy: "+ISOTOPE_NAME[isotope]).c_str(),900,600);
  hfake->Draw("HIST");
  

  // Save a PNG
  string title="energySmear0nu_"+ISOTOPE_NAME[isotope]+".png";
  c->SaveAs(title.c_str());
  delete c;
  return hfake;
}

TH1D * makeSmearedHistogram(ISOTOPE isotope, bool is2nu, double resolutionAt1MeV)
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
    double smeared1 = Smear(electronEnergy->at(0), resolutionAt1MeV);
    double smeared2 = Smear(electronEnergy->at(1), resolutionAt1MeV);
    hsmeared->Fill(smeared1+smeared2);
  }
  TCanvas *c = new TCanvas (("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),("Smeared energy: "+ISOTOPE_NAME[isotope]).c_str(),900,600);
  htrue->Draw("HIST");
  hsmeared->SetLineColor(kRed);
  hsmeared->Draw("HIST SAME");
  
  // Add a legend
  TLegend* legend = new TLegend(0.75,0.8,0.9,0.9);
  legend->AddEntry(htrue, "True", "l");
  legend->AddEntry(hsmeared,"Smeared", "l");
  legend->Draw();
  
  // Save a PNG
  string title="energySmear0nu_"+ISOTOPE_NAME[isotope]+".png";
  if (is2nu)title="energySmear2nu_"+ISOTOPE_NAME[isotope]+".png";
  c->SaveAs(title.c_str());
  delete c;
  hsmeared->Sumw2();
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
