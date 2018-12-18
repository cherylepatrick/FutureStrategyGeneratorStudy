# FutureStrategyGeneratorStudy
To look at what we would need to do to cover the inverted hierarchy range.

Input - simulated 2nubb and 0nubb files. Generator level only - detector simulation is not really needed (but I don't know how to turn it off)

**TrimTree** -  selects only events with total KE > 2 MeV. Better solution is just to use simulation where this cut is already made.

**SmearHistogram** - applies a smearing of the individual kinetic energies by a specified percentage, and then histograms the smeared energies. *A better thing for it to do would be to save the smeared energies, so we could re-plot them with different binnings without re-doing the smearing.*

Arguments are Se or Nd, percentage to smear by, and an output ROOT filename. Se and Nd then lead to hard-coded input file names. *This is bad news, we should be passing in an input filename*

**GeneratorStudy** - takes an input isotope Se or Nd *it really should take filenames* which points to  hard coded filenames for 2nubb and 0nubb smeared distributions; sensitivity in years, experiment name, and output filename

It will try to work out, from the smeared distributions, how much exposure you would need to achieve the specified sensitivity for the given experiment (NEXO or LEGeND: the experiment name just goes into the naming; the sensitivity must be manually entered). The input file will contain an ensemble of smeared histograms, for each resolution, with a particular naming. It then writes to the specified output file the resulting graph of smearing vs exposure.

**FormatGraphs** Input is the root file output from the GeneratorStudy, with multiple graphs of smearing vs exposure (with a specific naming convention), and a general title. This program will plot all those graphs neatly on the same axes and put the title above.

