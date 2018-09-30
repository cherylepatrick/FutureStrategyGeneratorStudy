#include "GeneratorStudy.h"
#include "FormatGraphs.h"

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
  std::vector<TGraph*> graphs;
  
  TList* list = graphFile->GetListOfKeys() ;
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;
  while (( key = (TKey*)next() )) {
    obj = key->ReadObj() ;
    if ( obj->InheritsFrom("TGraph")) // It's a graph
    {
      graphs.push_back((TGraph*)obj);
    }
  }
  double max=0;
  int biggestGraph=0;
  std::vector<int> colors = LoadColors();
  TLegend *legend=new TLegend(0.1,0.5,0.48,0.9);
  for (int i=0;i<graphs.size();i++)
  {
    TGraph *graph=graphs.at(i);
    cout<<graph->GetName()<<": "<<graph->GetHistogram()->GetMaximum()<<endl;
    if (graph->GetHistogram()->GetMaximum() > max)
    {
      max=graph->GetHistogram()->GetMaximum();
      biggestGraph=i;
    }
    graph->SetLineColor(colors.at(i));
    graph->SetLineWidth(2);
    graph->GetYaxis()->SetRangeUser(0,max);
    graph->GetYaxis()->SetTitleOffset(1.5);
    graph->GetXaxis()->SetTitle("Percent resolution (sigma) at 1 MeV");
    string legendName = graph->GetName(); // Strip isotope prefix
    int underscorepos = legendName.find("_");
    legendName = legendName.substr(underscorepos + 1);
    legend->AddEntry(graph,legendName.c_str(),"l");
    graph->SetTitle(title.c_str());
  }
  cout<<graphs.at(biggestGraph)->GetName()<<" has the biggest max: "<<max<<endl;
  TCanvas *c = new TCanvas("c",title.c_str(),900,600);
  
  graphs.at(biggestGraph)->Draw();
  
  legend->Draw();
  for (int i=0;i<graphs.size();i++)
  {
    if (i!=biggestGraph)
    {
      graphs.at(i)->Draw("SAME");
    }
  }
  c->SetTitle(title.c_str());
  c->SaveAs(plotfilename.c_str());
  return 0;
}


std::vector<int>LoadColors()
{
  std::vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue+2);
  colors.push_back(kGreen+3);
  colors.push_back(kMagenta+3);
  colors.push_back(kOrange+1);
  colors.push_back(kCyan-2);
  colors.push_back(kPink+10);
  colors.push_back(kViolet-8);
  colors.push_back(kOrange+3);
  colors.push_back(kRed+3);
  colors.push_back(kBlue-7);
  colors.push_back(kSpring-1);
  colors.push_back(kGray+2);
  colors.push_back(kPink-9);
  return colors;
}
