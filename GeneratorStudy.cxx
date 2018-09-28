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

  string graphname="";
  double halflife=0;
  ISOTOPE isotope=SE82;
  string isotopetext ="";
  string outputfile = "exposure_graphs.root";
  if (argc < 2)
  {
    cout<<"Usage: "<<argv[0]<<" -i <se or nd> -s <sensitivity in years> -x <experiment name> -o <outputfile>"<<endl;
    return -1;
  }
  int flag=0;
  while ((flag = getopt (argc, argv, "i:s:x:o:")) != -1)
  {
    switch (flag)
    {
      case 'i':
        if ((std::toupper(optarg[0])=='S') && (std::toupper(optarg[1])=='E'))isotope=SE82;
        else if ((std::toupper(optarg[0])=='N') && (std::toupper(optarg[1])=='D'))isotope=ND150;
        else
        {
          cout<<"Invalid isotope, pick Se or Nd"<<endl;
          return 1;
        }
        break;
      case 's':
        halflife = atof(optarg);
        break;
      case 'x':
        graphname = optarg;
        break;
      case 'o':
        outputfile = optarg;
        break;
      case '?':
        if (optopt == 'i' || optopt == 's' || optopt == 'x' || optopt == 'o' )
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        cout<<"Usage: "<<argv[0]<<" -i <se or nd> -s <sensitivity in years> -x <experiment name> -o <outputfile>"<<endl;
        return 1;
      default:
        abort ();
    }
  }
  
  if (graphname=="" || halflife<=0)
  {
    cout<<"Usage: "<<argv[0]<<" -i <se or nd> -s <sensitivity in years> -x <experiment name> -o <outputfile>"<<endl;
    return 1;
  }
  cout<<"Writing graph "<<graphname<<" for "<<ISOTOPE_NAME[isotope]<<" with halflife "<<halflife<<" years to "<<outputfile<<endl;
  TGraph *g=MakeExposureGraph(graphname,isotope,halflife);
  
  TFile *outfile = new TFile(outputfile.c_str(),"UPDATE");
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
        string name0nu=histname;
        if (histname.find("smeared_2nu_")==0)
        {
          name0nu.replace(8,1,"0");
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

  double lowExposure=500;
  double highExposure=2000000;
  double low2nuEvents=Get2nuEventsForExposure(lowExposure, isotope);
  double high2nuEvents=Get2nuEventsForExposure(highExposure, isotope);
  
  // Get the signal event limit for those two
  TH1D *high2nu=(TH1D*)hist2nu->Clone();
  high2nu->Scale(high2nuEvents / TOTAL_2NU_EVENTS[isotope]);
  double highEventLimit=ExpectedLimitSigEvts(DESIRED_CONFIDENCE, hist0nu, high2nu, high2nu );
  
  TH1D *low2nu=(TH1D*)hist2nu->Clone();
  low2nu->Scale(low2nuEvents / TOTAL_2NU_EVENTS[isotope]);
  double lowEventLimit=ExpectedLimitSigEvts(DESIRED_CONFIDENCE, hist0nu, low2nu, low2nu );
  
  // The signal event limits don't seem to change THAT Much with exposure, so these can give us new starting points for our search
  lowExposure=GetExposureFrom0nuEvents(lowEventLimit, isotope, desiredHalflife) * 0.5;
  highExposure=GetExposureFrom0nuEvents(highEventLimit, isotope, desiredHalflife) * 2.;
  
  low2nuEvents=Get2nuEventsForExposure(lowExposure, isotope);
  high2nuEvents=Get2nuEventsForExposure(highExposure, isotope);
  cout<<hist2nu->GetName()<<endl;
  cout<<"Starting from "<<lowExposure<<" kg.years (half of "<<lowEventLimit<<" events) to "<<highExposure<<" kg.years ( twice "<<highEventLimit<<" events)"<<endl;

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
    cout<<"Total 2 nu events: "<<this2nuEvents<<endl;
    int binnumber=this2nu->FindBin(2.8);
    cout<<"Of which "<<this2nu->Integral(binnumber,301)<<" are above 2.8MeV"<<endl;
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

double GetExposureFrom0nuEvents(double events, ISOTOPE isotope, double halflife)
{
  // Inverse of the other calculation: events divided by number of atoms per kg, times decay constant
  return events / ( (AVOGADRO * 1000 / ATOMIC_MASS[isotope]) * (TMath::Log(2) / halflife) );
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

