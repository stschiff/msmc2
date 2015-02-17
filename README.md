*Disclaimer: This software is still under development and mainly for internal use only at the moment.*

This program implements MSMC2, a method to infer population size history and population separation history from whole genome sequencing data.

To install, you need a modern D compiler, for example [DMD](http://dlang.org/download.html).
You also need the GNU Scientific library [GSL](http://www.gnu.org/software/gsl/) installed.

You can then install the program by typing `make` in the directory.
The resulting executable will be in the `build/release` subdirectory.

You have to manually adjust the Makefile if you have the GSL in a non-trivial location on your platform.

    Usage: msmc2 [options] <datafiles>
      Options:
        -i, --maxIterations=<size_t> :      number of EM-iterations [default=20]
        -o, --outFilePrefix=<string> :      file prefix to use for all output files
        -r, --rhoOverMu=<double> :          ratio of recombination over mutation rate (default: 0.25)
        -t, --nrThreads=<size_t> :          nr of threads to use (defaults to nr of CPUs)
        -p, --timeSegmentPattern=<string> : pattern of fixed time segments [default=1*4+25*2+1*4+1*6]
        -R, --fixedRecombination :          keep recombination rate fixed (rarely needed in MSMC2)
        -I, --indices:                      indices (comma-separated) of alleles in the data file to run over
        -P, --subpopLabels:                 comma-separated list of 0s and 1s to indicate subpopulations. If given, 
                                            estimate coalescence rates only across populations.
        -s, --skipAmbiguous:                skip sites with ambiguous phasing. Recommended for gene flow analysis

minimum command line:

    build/release/msmc2 -o <out_prefix> <input_chr1> <input_chr2> ...

