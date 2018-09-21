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
  
  TGraph * sigevents_se82 = SigEventsVsResolution(SE82);

  GetExposure(sigevents_se82, "LEGEND", SE82, SENSITIVITY_LEGEND_Se);
  GetExposure(sigevents_se82, "nEXO", SE82, SENSITIVITY_NEXO_Se);
  
  return 0;
}

TGraph* GetExposure(TGraph *sigevents, string compExperiment, ISOTOPE isotope, double desiredSensitivity)
{
  
    // From previous sensitivity calculations:
  //double totalTLimitSensitivity= (zeroNuEfficiency/totalExpectedSignalEventLimit) * ((se82IsotopeMass*1000 * AVOGADRO)/se82MolarMass) * TMath::Log(2) * exposureYears;
  
  // Rearrange to get exposure in kg-years needed to reach desired sensitivity
  // We want to scale our graph by the number of kg years divided by the number of signal events
  
  double scale = desiredSensitivity *  ATOMIC_MASS[isotope]/ (AVOGADRO  * 1000. * TMath::Log(2)); // Necessary exposure increases with sensitivity and how many events are needed. Decreases with the number of atoms of isotope per kg. Log 2 is because we want a half-life not a decay constant. 1000 is because a mole of material is (atomic mass in grams)
  
  TGraph *exposure=ScaledClone(sigevents,scale);
  exposure->GetYaxis()->SetTitle("Exposure (kg.years)");
  string title="Exposure for "+ISOTOPE_LATEX[isotope]+" if "+compExperiment+" sees something";
  exposure->SetTitle(title.c_str());
  string pngTitle="exposure_"+ISOTOPE_NAME[isotope]+"_"+compExperiment+".png";
  TCanvas *c1=new TCanvas("c1","c1",900,600);
  exposure->Draw();
  c1->SaveAs(pngTitle.c_str());
  delete c1;
  return exposure;
}

TGraph *ScaledClone(TGraph *input, double scale)
{
  std::vector<double>x;
  std::vector<double>y;
  for (int i=0;i<input->GetN();i++)
  {
    x.push_back(input->GetX()[i]);
    y.push_back(input->GetY()[i] * scale);
  }
  
  TGraph *output= new TGraph(x.size(),&x[0],&y[0]);
  output->SetTitle(input->GetTitle());
  output->GetXaxis()->SetTitle(input->GetXaxis()->GetTitle());
  output->GetYaxis()->SetTitle(input->GetYaxis()->GetTitle());
  return output;
}

TGraph* SigEventsVsResolution(ISOTOPE isotope)
{

  std::vector<double>resolutions;
  resolutions.push_back(0.01);
  resolutions.push_back(0.02);
  resolutions.push_back(0.03);
  resolutions.push_back(0.04);
  resolutions.push_back(0.05);

  std::vector<double>sigEvents;
  for (int i=0;i<resolutions.size();i++)
  {
    sigEvents.push_back(SigEventLimit(isotope, resolutions.at(i)));
  }
  TGraph *eventsGraph= new TGraph(resolutions.size(), &resolutions[0], &sigEvents[0]);
  for (int i=0;i<resolutions.size();i++)
  {
    cout<<resolutions.at(i)<<":"<<sigEvents.at(i)<<endl;
  }
  
  eventsGraph->SetTitle(("Signal event limit "+ISOTOPE_LATEX[isotope]).c_str());
  eventsGraph->GetXaxis()->SetTitle("Fractional resolution at 1MeV");
  eventsGraph->GetYaxis()->SetTitle("Signal event limit");
  eventsGraph->Print("ALL");
  TCanvas *c=new TCanvas("c","c",900,600);
  //TGraph * sigevents_se82 = SigEventsVsResolution(SE82);
  eventsGraph->Draw();
  c->SaveAs(("signal_events_"+ISOTOPE_NAME[isotope]+".png").c_str());
  delete c;
  return eventsGraph;
}

double SigEventLimit(ISOTOPE isotope, double resolutionAt1MeV)
{
  TH1D *smeared2nu = makeSmearedHistogram(isotope,true,resolutionAt1MeV);
  TH1D *smeared0nu = makeSmearedHistogram(isotope,false, resolutionAt1MeV);
  smeared0nu->GetSumw2();
  smeared2nu->GetSumw2();
  // Plot the two together
  TCanvas *c = new TCanvas (("totalEnergy_"+ISOTOPE_NAME[isotope]).c_str(),("Energy: "+ISOTOPE_NAME[isotope]).c_str(),900,600);
  
  smeared2nu->SetLineColor(kBlack);
  smeared2nu->Draw("HIST");
  smeared0nu->Scale(smeared2nu->Integral()/smeared0nu->Integral() * 0.1); // Kind of an arbitrary scale for the 0nu sample, just for show
  //smeared0nu->Scale(smeared2nu->Integral()/smeared0nu->Integral() * 1e-5); // Approximate scaling - 0nu halflives are about 1e27, 2nu are about 1e20, but we only have ~1% of 2nu events. But it looks like the scale doesn't matter?
  smeared0nu->SetLineColor(kRed);
  smeared0nu->Draw("HIST SAME");
  
  // Add a legend
  TLegend* legend = new TLegend(0.75,0.8,0.9,0.9);
  legend->AddEntry(smeared2nu, "2#nu#beta#beta", "l");
  legend->AddEntry(smeared0nu,"0#nu#beta#beta", "l");
  legend->Draw();
  
  // Save a PNG
  string title="smearedComparison_"+ISOTOPE_NAME[isotope]+Form("_%.0f.png",resolutionAt1MeV*100);
  c->SaveAs(title.c_str());

  // Calculate signal event limit
  
  double signalEvents=ExpectedLimitSigEvts(DESIRED_CONFIDENCE, smeared0nu, smeared2nu, smeared2nu ); // "Data" is background
  cout<<signalEvents<<" signal events needed"<<endl;

  delete c;
//  return kgYears;
  return signalEvents;
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
  //for (int i=0;i<20000;i++) // Just to make it run faster
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

//double EstimateBackgroundEvents(double backgroundEfficiency, double isotopeMass, double molarMass, double halfLife)
//{
//  // Get the number of atoms you start with
//  double nSourceAtoms=AVOGADRO * (isotopeMass*1000)/molarMass; //Avogadro is grams
//  // The exponent is around 10^-20, it's too small for TMath::Exp to deal with
//  // Can we just go for a Taylor expansion for 1-e^-x where x is v small?
//  // 1( - e^-x) ~ x so...
//  double totalDecays=nSourceAtoms * (TMath::Log(2) * exposureYears/halfLife);
//  // Multiply by the efficiency and that is the amount of background events you expect to see
//  double events=totalDecays * backgroundEfficiency;
//  //cout<<totalDecays<<" backgrounds, of which we see "<<events<<endl;
//  return events;
//}

// Adapted from James Mott's LimitCalculationFunctions, thanks James!
double ExpectedLimitSigEvts(double ConfidenceLevel, TH1D* h_signal, TH1D* h_background, TH1D* h_data ) {
  
  // The idea is to calculate a confidence level with different signal strengths, to indicate
  // that we have seen no signal when we pass it a background-only sample
  // We will stop when we find the signal strength that gives the confidence level we pass in
  // We will start with by scaling the signal up and down from what we think is a sensible answer
  // One of these will give too high a confidence level, so we can start lowering it
  // (If the 0nubb halflife were short, the expected events would be high, so if we saw none
  // we would be very confident that the signal was not there)
  // The other limit will give too low a confidence level so we must raise it
  // (We expected so few events that we can't tell if there's signal or not, it would be hidden by background)
  
  

  // These numbers are multipliers that will set us to 0.1 signal events and
  // 1000 signal events respectively when we scale the plots
  double low_bound = 0.1/h_signal->Integral();
  double high_bound = 1000.0/h_signal->Integral();
  
  // Start with the small signal of 0.1 signal events for our low bound
  // And 1k for the high bound
  // We can't tell if we have signal or not with this tiny signal strength
  TH1D* null_hyp_signal = (TH1D*) h_signal->Clone("null_hyp_signal"); null_hyp_signal->Scale(low_bound);
  // We could easily tell if we had signal as big as this
  TH1D* disc_hyp_signal = (TH1D*) h_signal->Clone("disc_hyp_signal"); disc_hyp_signal->Scale(high_bound);
  
  // Set up a data source using the small signal
  // Data is the same as background, as we are setting a limit - assume we have no signal
  TLimitDataSource* mydatasource = new TLimitDataSource(null_hyp_signal, h_background, h_data);
  // Calculate a limit using CLs method to get a low-bound confidence level
  TConfidenceLevel* myconfidence = TLimit::ComputeLimit(mydatasource, TLIMIT_EXPERIMENTS);
  double low_bound_cl = myconfidence->CLs(); // This should be below our desired confidence level
  delete mydatasource;
  
  // Now do the same with the scaled-up signal to get a high-bound confidence level
  // We hope that our desired confidence level (input to function) is between these
  mydatasource = new TLimitDataSource(disc_hyp_signal, h_background, h_data);
  myconfidence = TLimit::ComputeLimit(mydatasource, TLIMIT_EXPERIMENTS);
  double high_bound_cl = myconfidence->CLs(); // This should be above our desired confidence level
  // Confidence level for signal is confidence for (signal+bg) / confidence for just background
  delete mydatasource;
  
  // Now we are going to try different numbers of signal events, between those two bounds
  double accuracy = 0.01;
  double this_cl = 0;
  double this_val = 0; // Scale factor for this test - between the high and low bounds
  
  // Now we are going to close in those bounds until the difference between the
  // low and high bound integrals is less than 1% of our initial signal integral, and
  // the desired confidence level falls in between them
  while  (fabs(high_bound - low_bound) * h_signal->Integral() > accuracy) {
    // Pick a new number scale factor between the low and high bounds (nearer to the low) to try next
    // remember that low_bound and high_bound (and thus this_val) are not the number of events -
    // they are a scale factor for the histogram i.e. number of events divided by original signal integral
    this_val = low_bound+(high_bound - low_bound)/3.;
    TH1D* this_signal = (TH1D*) h_signal->Clone("test_signal");
    this_signal->Scale(this_val);
    
    // Calculate a confidence level...
    mydatasource = new TLimitDataSource(this_signal, h_background, h_data);
    myconfidence = TLimit::ComputeLimit(mydatasource, TLIMIT_EXPERIMENTS);
    this_cl = myconfidence->GetExpectedCLs_b(); // Get a new confidence level
    // This is a different calculation than we did before!
    // The documentation for TLimit does not say what it does, but my guess is
    // This is(expected Confidence Level for the signal plus background hypothesis if there is only background)
    // divided by (expected Confidence Level for the background only if there is only background.)
    // So it's a confidence level for signal if there is only background
    if (this_cl > ConfidenceLevel) {
      low_bound = this_val; // If it's higher than our desired confidence level, move the low bound up
      low_bound_cl = this_cl;
    } else { // If it's lower than our desired level, we must move the high bound down
      high_bound = this_val;
      high_bound_cl = this_cl;
    }
    delete mydatasource;
    delete this_signal;
    delete myconfidence;
  }
  // Tidy up
  delete null_hyp_signal;
  delete disc_hyp_signal;
  
  // Number of events is the scale factor * the original integral
  return h_signal->Integral() * this_val;
}
void Renormalize(ISOTOPE isotope, TH1D* hist)
{
  // I'm not sure if we actually NEED the 2nu halflife if we aren't looking at backgrounds? It doesn't seem to matter if I scale it. And we don't  know how long we will run for, so we can't turn it to a number of events.
}
