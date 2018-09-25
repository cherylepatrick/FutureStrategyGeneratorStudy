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
  TFile *outfile = new TFile("exposure_graphs.root","UPDATE");
  TGraph *g=MakeExposureGraph("LEGEND average",SE82,SENSITIVITY_LEGEND_Se);
  
  outfile->cd();
  g->Write("",TObject::kOverwrite);
  return 0;
}

TGraph* MakeExposureGraph(string experimentText,ISOTOPE isotope,double desiredHalflife)
{
  cout<<"Experiment: "<<experimentText<<" isotope "<<ISOTOPE_NAME[isotope]<<endl;
  vector<TH1D*> smeared2nuPlots;
  vector<TH1D*> smeared0nuPlots;

  // Load the plots from the smeared plots file
  
  TFile *f=new TFile(SMEARED_HISTO_FILE[isotope].c_str());
  TList* list = f->GetListOfKeys() ;
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;
  
  while (( key = (TKey*)next() )) {
    obj = key->ReadObj() ;
    if ( obj->InheritsFrom("TH1")) // It's a histogram
    {
      try {
        // If it is a 2nu histogram, write it to smeared2nuplots
        string histname=obj->GetName();
        if (histname.find("smeared_2nu_")==0)
        {
          string name0nu=histname.replace(8,1,"0");
          if (list->Contains(name0nu.c_str()))
          {
            TH1D *hist0nu=(TH1D*)f->Get(name0nu.c_str());
            smeared2nuPlots.push_back((TH1D*)obj);
            smeared0nuPlots.push_back(hist0nu);
            cout<<"Loading histograms "<<histname<<" and "<<name0nu<<endl;
          }
          else
          {
            cout<<"Couldn't get corresponding 0nu histogram for "<<obj->GetName()<<endl;
          }
          // Get the corresponding 0nu histogram, write it to smeared0nuplots
        }
      } catch (exception &e) {
        cout<<"Couldn't get corresponding 2nu and 0nu histograms for "<<obj->GetName()<<endl;
      }
    }
  }

  if (smeared2nuPlots.size()==0 || smeared2nuPlots.size() != smeared0nuPlots.size())
  {
    cout<<"Unable to load smeared histograms, quitting"<<endl;
    return (TGraph*) NULL;
  }

  
  // Calculate the exposure needed for each smearing
  vector<double> smearings;
  vector<double> exposures;
//  for (int i=0;i<1;i++)
  for (int i=0;i<smeared2nuPlots.size();i++)
  {
    // Find out the smearing percentage
    string plotname=smeared2nuPlots.at(i)->GetName();
    int pos=plotname.find(ISOTOPE_NAME[isotope]) + ISOTOPE_NAME[isotope].length()+1;
    string fractionString=plotname.substr(pos);
    pos=fractionString.find("_");
    fractionString=fractionString.replace(pos,1,".");
    double smearing=100. * atof(fractionString.c_str());
    
    double exposure=GetExposure(smeared2nuPlots.at(i),smeared0nuPlots.at(i),isotope,desiredHalflife);
    cout<<"The best exposure is "<<exposure<<" for "<<smearing<<"% smearing"<<endl;
    smearings.push_back(smearing);
    exposures.push_back(exposure);
  }
  
  // Unfortunately for a TGraph we need to sort these. Ugh.
  vector<double> sortedSmearings;
  vector<double> sortedExposures;

  // This is not efficient but I don't think that really matters as they will be small vectors. Might change mind if they get big
  double currentSmear=0;
  for (int i=0;i<smearings.size();i++)
  {
    int nextSmearing = GetNextSmearing(smearings, currentSmear);
    sortedSmearings.push_back(smearings.at(nextSmearing));
    sortedExposures.push_back(exposures.at(nextSmearing));
    currentSmear=smearings.at(nextSmearing);
  }

  // Plot smearing vs exposure
  TGraph *exposureGraph = new TGraph(sortedSmearings.size(),&sortedSmearings[0], &sortedExposures[0]);
  exposureGraph->SetName((ISOTOPE_NAME[isotope]+"_"+experimentText).c_str());
  exposureGraph->SetTitle((ISOTOPE_LATEX[isotope]+": "+experimentText).c_str());
  exposureGraph->GetXaxis()->SetTitle("Percentage resolution at 1MeV");
  exposureGraph->GetYaxis()->SetTitle("Required exposure (kg.year)");

  TCanvas *c=new TCanvas("c","c",900,600);
  exposureGraph->Draw();
  c->SaveAs(("exposures_"+experimentText+"_"+ISOTOPE_NAME[isotope]+".png").c_str());
  delete c;
  exposureGraph->SetTitle(experimentText.c_str());
  f->Close();
  return exposureGraph;
}

// Return the position of the lowest number in the vector that is bigger than the input
int GetNextSmearing(std::vector<double> smearings, double previousSmearing)
{
  int pos=-9999;
  double lowest=999999999;
  for (int i=0;i<smearings.size();i++)
  {
    if (smearings.at(i) > previousSmearing && smearings.at(i) < lowest)
    {
      lowest = smearings.at(i);
      pos=i;
    }
  }
  return pos;
}

double GetExposure(TH1D *hist2nu, TH1D *hist0nu, ISOTOPE isotope, double desiredHalflife)
{
  // ######### take this out!!!!!!!!!
 // return 3500 + trand->Gaus(0,500);
  // Here are the steps:
  // We have a desired 0nubb halflife. That corresponds to an expected number of signal events
  // depending on the exposure (mass x time)
  // The exposure is proportional to the amount of 2nubb, as we know the 2nubb halflife
  // Find the amount of 2nubb that gives us the number of signal events corresponding to our
  // desired halflife
  
  if( hist2nu->GetSumw2N() == 0 )hist2nu->Sumw2();
  if( hist0nu->GetSumw2N() == 0 )hist0nu->Sumw2();

  // If we are worried our limits might not be wide enough, we could check here
//  TH1D *low2nu = (TH1D*)hist2nu->Clone();
//  TH1D *high2nu = (TH1D*)hist2nu->Clone();

  double lowExposure=1000;
  double highExposure=50000;
  double low2nuEvents=Get2nuEventsForExposure(lowExposure, isotope);
  double high2nuEvents=Get2nuEventsForExposure(highExposure, isotope);
  
//  double low0nuEvents=Get0nuEventsForExposure(lowExposure, isotope, desiredHalflife);
//  double high0nuEvents=Get0nuEventsForExposure(highExposure, isotope, desiredHalflife);
  
  double this2nuEvents=0;
  double this0nuEvents=0;
  double thisExposure=0;
  double signalEventLimit=-9999;
  double accuracy=0.1;
  
  // Try different exposures to find one where the signal event limit
  // matches the number of 0nubb events we expect to get
  while  (fabs(signalEventLimit - this0nuEvents) > accuracy)
  {
    this2nuEvents = low2nuEvents + (high2nuEvents - low2nuEvents)/3;
    thisExposure= GetExposureFrom2nuEvents(this2nuEvents, isotope);
    // This how many 0nu events we would get at our desired halflife
    this0nuEvents=Get0nuEventsForExposure(thisExposure,isotope, desiredHalflife);
    // Now see how the signal event limit compares to this number of events
    // to find out whether we could detect them or not
    TH1D *this2nu=(TH1D*)hist2nu->Clone();
    this2nu->Scale(this2nuEvents / TOTAL_2NU_EVENTS[isotope]);
    double signalEventLimit=ExpectedLimitSigEvts(DESIRED_CONFIDENCE, hist0nu, this2nu, this2nu );
    cout<<"**** Exposure "<<thisExposure<<" would have "<<this0nuEvents<<" 0nu and the limit is "<<signalEventLimit<<endl;
    if (fabs(signalEventLimit - this0nuEvents) > accuracy)
    {
      // If the limit is higher than the halflife, we won't see it and we need more events
      // If it's lower, we can get away with less exposure
      if (signalEventLimit < this0nuEvents)
        high2nuEvents = this2nuEvents;
      else
        low2nuEvents =  this2nuEvents;
    }
    else
      return thisExposure;
  }
  return thisExposure;

}

double Get2nuEventsForExposure(double exposure, ISOTOPE isotope)
{
  // Number of atoms times exposure time / 2nubb decay constant
  return exposure * (AVOGADRO * 1000 / ATOMIC_MASS[isotope]) * (TMath::Log(2) / HALFLIFE2NU[isotope]); // 1000 is because we have kg years exposure, but atomic mass is grams per mole
}

double Get0nuEventsForExposure(double exposure, ISOTOPE isotope, double desiredHalflife)
{
  // Number of atoms times exposure time / 2nubb decay constant
  return exposure * (AVOGADRO * 1000 / ATOMIC_MASS[isotope]) * (TMath::Log(2) / desiredHalflife); // 1000 is because we have kg years exposure, but atomic mass is grams per mole
}

double GetExposureFrom2nuEvents(double events, ISOTOPE isotope)
{
  // Inverse of the other calculation: events divided by number of atoms per kg, times decay constant
  return events / ( (AVOGADRO * 1000 / ATOMIC_MASS[isotope]) * (TMath::Log(2) / HALFLIFE2NU[isotope]) );
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
