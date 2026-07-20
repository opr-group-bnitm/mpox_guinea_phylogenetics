# Mpox phylogenies

This repository contains scripts and data to reproduce the results presented in our manuscript "Emergence of mpox in Guinea: sporadic clade IIa cases and a clade IIb upsurge".

The analysis for clade IIa and clade IIb genomes found in guinea are in two separate folders (IIa and IIb). The directory both_clades contains scripts to create the figures in the manuscript. 

## Installation of dependencies

You can install all dependencies with `conda install -f env.yaml`. The calculations were performed using this environment on a MacBook ProIntel Core i9 computer with macOS Sequoia 15.7.7. The environment may have to be adapted to run on another machine.

## Computing alignments and phylogenies

Most of the analysis can be reproduced by simply running the scripts according to the prefixes (step1, ...). To re-run the beast analysis, read the next section.

## Beast configuration and output

Our configuration files were created with the tool beauti. You can find them in results/beast_in.
If you want to create these files yourself, run step3_prepare_beast_alignment.py, load it into beauti and set the parameters as described in the manuscript. You can then run step5_run_beast.sh. Then use Tracer to make sure to all trajacteries converged. Finally, merge the trajectories and create a maximum clade credibility tree with step6_combine_logs.sh.

## Results

To not overload this repository, most output files are not included. However, the files required to create the plots are included (trees) and the combined beast log to inspect the outputs.
Beware, that if you recompute the results, these files will be overwritten.

## Plotting

The plots in the paper can be reproduced with the jupyter notebook plot.ipynb. The relevant trees and metadata are provided in this repository (see previous section).
