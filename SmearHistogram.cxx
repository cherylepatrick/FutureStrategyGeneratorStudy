//#include "GeneratorStudy.h"
#include "SmearHistogram.h"
TRandom3 *trand;

/**
 *  main function
 * Arguments are <root file>
 */
int main(int argc, char **argv)
{
  
  trand = new TRandom3();
  gStyle->SetOptStat(0);
  string helptext="<Se or Nd> <smear percent> <2nu root file> <0nu root file> <output root file>";
  
  double smearing=0;
  ISOTOPE isotope;
  
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

  string inFile2nu="";
  if (argc>3)
  {
    string arg3=argv[3];
    if (arg3.find(".root") + 5 == arg3.length() && arg3.length() > 4)
    {
      inFile2nu=arg3;
      
      cout<<"Input 2nufile is "<<inFile2nu<<endl;
    }
  }
  if (inFile2nu.length()==0)
  {
    cout<<"Loading 2nu from default files ";
  }
  
  
  string inFile0nu="";
  if (argc>4)
  {
    string arg4=argv[4];
    if (arg4.find(".root") + 5 == arg4.length() && arg4.length() > 4)
    {
      inFile0nu=arg4;
      
      cout<<"Input 0nu file is "+inFile0nu<<endl;
    }
  }
  if (inFile0nu.length()==0)
  {
    cout<<"Loading 0nu from default files ";
  }

  string filename="";
  
  if (argc>5)
  {
    string arg5=argv[5];
    if (arg5.find(".root") + 5 == arg5.length() && arg5.length() > 4)
    {
      cout<<"Saving to file ";
      filename=arg5;
    }
  }
  if (filename.length()==0)
  {
    cout<<"Saving to default file ";
    filename = "smeared_hists_"+ISOTOPE_NAME[isotope]+".root";
  }
  cout<<filename<<endl;
  string smeartext=Form("_%f",smearing*0.01);
  int pos=smeartext.find(".");
  if (pos>=0)
    smeartext.replace(pos,1,"_");

  // 2nubb smeared histogram
  TH1D *hist2nu=makeSmearedHistogram(isotope, true, smearing*0.01, inFile2nu);
  hist2nu->SetName(("smeared_2nu_"+ISOTOPE_NAME[isotope]+smeartext).c_str());
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  hist2nu->Write("",TObject::kOverwrite);
  
  // 0nubb smeared histogram
  TH1D *hist0nu=makeSmearedHistogram(isotope, false, smearing*0.01, inFile0nu);
  hist0nu->SetName(("smeared_0nu_"+ISOTOPE_NAME[isotope]+smeartext).c_str());
  f->cd();
  hist0nu->Write("",TObject::kOverwrite);
  f->Close();
  return 0;
}


TH1D * makeSmearedHistogram(ISOTOPE isotope, bool is2nu, double resolutionAt1MeV, string fName)
{
  if (fName.length()==0)
  {
    if (is2nu)fName=FILES2NU[isotope]; else fName=FILES0NU[isotope];
  }
  TFile *f = new TFile(fName.c_str());
  TTree *tree = (TTree*) f->Get(treeName.c_str());
  double maxEnergy=Qbb[isotope];
  TH1D *htrue = new TH1D("htrue",(ISOTOPE_LATEX[isotope]+(is2nu?" 2#nu#beta#beta":" 0#nu#beta#beta")+Form(" %.1f%% smearing",100*resolutionAt1MeV)).c_str(),500,2.1,maxEnergy*1.4);
  TH1D *hsmeared = new TH1D("hsmeared",(ISOTOPE_LATEX[isotope]+(is2nu?" 2#nu#beta#beta":" 0#nu#beta#beta")+Form(" %.1f%% smearing",100*resolutionAt1MeV)).c_str(),500,2.1,maxEnergy*1.4);
  vector<double> *electronEnergy = 0;
  tree->SetBranchAddress("trueparticle.kinenergy", &electronEnergy);
  // Loop the entries
  int nEntries = tree->GetEntries();
  cout<<"number of entries "<<nEntries<<endl;
  
  // Prepare to write the smeared values to an output root file
  string smearName=Form("_smear_%.2f.root",resolutionAt1MeV);
  int thePos=smearName.find(".");
  if (thePos>=0)
    smearName.replace(thePos,1,"_");
  string outFilename=fName;
  outFilename.replace(outFilename.length()-5,5,smearName);
  
  TFile *outFile = new TFile (outFilename.c_str(),"RECREATE","Smeared data");
  outFile->cd();
  TTree *newTree = new TTree("SmearedData","SmearedData");
  newTree->SetDirectory(outFile);
  vector<double> *smearedEnergy = 0;
  
  newTree->Branch("smeared.kinenergy",&smearedEnergy);
  for (int i=0;i<nEntries;i++)
  //for (int i=0;i<20;i++) // Just to make it run faster
  {
    tree->GetEntry(i);
    double totalEnergy=electronEnergy->at(0)+electronEnergy->at(1);
    htrue->Fill(totalEnergy); // Do we actually want to do anything with this?
    double smeared1 = Smear(electronEnergy->at(0), resolutionAt1MeV);
    double smeared2 = Smear(electronEnergy->at(1), resolutionAt1MeV);
    hsmeared->Fill(smeared1+smeared2);
    smearedEnergy->clear();
    smearedEnergy->push_back(smeared1);
    smearedEnergy->push_back(smeared2);
    newTree->Fill();
  }
  
  outFile->cd();
  newTree->Write();
  outFile->Close();
  
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
  string smeartext=Form("_%f",resolutionAt1MeV);
  int pos=smeartext.find(".");
  if (pos>=0)
    smeartext.replace(pos,1,"_");
  string title="energySmear0nu_"+ISOTOPE_NAME[isotope]+smeartext+".png";
  if (is2nu)title="energySmear2nu_"+ISOTOPE_NAME[isotope]+smeartext+".png";
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

