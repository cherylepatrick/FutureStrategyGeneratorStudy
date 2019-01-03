# FutureStrategyGeneratorStudy
To look at what we would need to do to cover the inverted hierarchy range.

Input - simulated 2nubb and 0nubb files. Generator level only - detector simulation is not really needed (but I don't know how to turn it off)

**TrimTree** -  selects only events with total KE > 2 MeV. Better solution is just to use simulation where this cut is already made.

**PruningModule** Falaise reconstruction module. Run it on a simulated file to prune it so the only branch is trueparticle.kinenergy, and output in .ROOT format. *Should we also add four-vectors just in case?*

**SmearHistogram** - applies a smearing of the individual kinetic energies by a specified percentage, and then histograms the smeared energies. *A better thing for it to do would be to save the smeared energies, so we could re-plot them with different binnings without re-doing the smearing.*

<Se or Nd> <smear percent> <2nu root file> <0nu root file> <output root file>

The Se and Nd are just for setting the axis scale and legend/title etc on the output plots. They don't affect the actual files it uses.

**GeneratorStudy** -It will try to work out, from the smeared distributions, how much exposure you would need to achieve the specified sensitivity for the given experiment (NEXO or LEGeND: the experiment name just goes into the naming; the sensitivity must be manually entered). The input file will contain an ensemble of smeared histograms, for each resolution, with a particular naming. It then writes to the specified output file the resulting graph of smearing vs exposure.

cout<<"Usage: "<<argv[0]<<" -i <se or nd> -s <sensitivity in years> -x <experiment name> -o <outputfile>"<<endl;

The input file is the smeared histograms, and it is getting the name of it from a hardcoded name per isotope. *That's scrappy but anyway - it shouldn't be taking smeared histograms, it should be taking the smeared energies and histogramming them on the fly.*

**FormatGraphs** Input is the root file output from the GeneratorStudy, with multiple graphs of smearing vs exposure (with a specific naming convention), and a general title. This program will plot all those graphs neatly on the same axes and put the title above.

