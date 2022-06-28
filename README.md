This program implements MSMC2, a method to infer population size history and population separation history from whole genome sequencing data. For a general guide, please refer to version 1: https://github.com/stschiff/msmc.

For reference see the following two publications:

* Schiffels, Stephan, and Ke Wang. 2020. “MSMC and MSMC2: The Multiple Sequentially Markovian Coalescent.” In Statistical Population Genomics, edited by Julien Y. Dutheil, 2090:147–66. Methods in Molecular Biology. New York, NY: Springer US.
* Wang, Ke, Iain Mathieson, Jared O’Connell, and Stephan Schiffels. 2020. “Tracking Human Population Structure through Time from Whole Genome Sequences.” PLoS Genetics 16 (3): e1008552.

Binaries are attached to the releases (under the "Releases" tab within github).

To build yourself, you need a modern D compiler, for example [DMD](http://dlang.org/download.html).
You also need the GNU Scientific library [GSL](http://www.gnu.org/software/gsl/) installed.

You can then install the program by typing `make` in the directory.
The resulting executable will be in the `build/release` subdirectory.

You have to manually adjust the Makefile if you have the GSL in a non-trivial location on your platform.

Options:

    -i, --maxIterations=<size_t> :      number of EM-iterations [default=20]
    -o, --outFilePrefix=<string> :      file prefix to use for all output files
    -r, --rhoOverMu=<double> :          initial ratio of recombination over mutation rate (default:
                                        0.25)
    -t, --nrThreads=<size_t> :          nr of threads to use (defaults to nr of CPUs)
    -p, --timeSegmentPattern=<string> : pattern of fixed time segments [default=1*2+25*1+1*2+1*3]
    -R, --fixedRecombination :          keep recombination rate fixed (rarely needed in MSMC2)
    -I, --pairIndices:                  this can be given in two flavors. First, you can enter a
                                        single comma-separated list like this "-I 0,1,4,5".
                                        In this case, the program will
                                        run over all pairs of haplotypes within this set of
                                        indices. This is useful for running on multiple phased
                                        diploid genomes sampled from one population.
                                        In the second flavor, you can give a list of pairs, like
                                        this: "-I 0-1,2-3,4-5". In this case, the
                                        program will run only those specified pairs. This can be
                                        used to run on a number of unphased genomes, to avoid pairs
                                        of haplotypes from different individuals. This should also
                                        be used to indicate a cross-population run, where you want
                                        to run the program only over pairs of haplotypes across
                                        population boundaries. So with two phased genomes, one from
                                        each population you'd run "-I 0-2,0-3,1-2,1-3", and the
                                        program would run only those four pairs of haplotypes.
                                        Note that if you do not use this parameter altogether,
                                        MSMC2 will run on all pairs of input haplotypes.
    -s, --skipAmbiguous:                skip sites with ambiguous phasing. Recommended for cross
                                        population analysis
    --quantileBoundaries:               use quantile boundaries, as in MSMC. To fully replicate
                                        MSMC's time intervals, combine this with -p 10*1+15*2

minimum command line:

    build/release/msmc2 -o <out_prefix> <input_chr1> <input_chr2> ...

#Brief Guide
For population size estimation, you can basically use `msmc2` just as you used `msmc`, using the same input files. Note there are some important changes to `msmc`: First, it is not anymore necessary to fix the recombination size for more than two haplotypes. So I would run without the `-R` switch (which was recommended in `msmc`). Second, the time patterning has changed. We now use the same patterning as in (Li and Durbin, 2011), however with more and more resolution in recent times depending on the number of haplotypes you use.

For cross-population analysis, things are quite different now. While in `msmc`, you could estimate all three coalescence rate functions (two within and one across populations) in one go, now you have to make three runs, corresponding to the three estimates. For example, say you have four diploid individuals, two from each population, you should generate a combined input file with eight haplotypes (see msmc and msmc-tools repositories), and then start three runs:

1. `build/release/msmc2 -I 0,1,2,3 -o within1_msmc <input_chr1> <input_chr2> ...`
2. `build/release/msmc2 -I 4,5,6,7 -o within2_msmc <input_chr1> <input_chr2> ...`
3. `build/release/msmc2 -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -o across_msmc <input_chr1> <input_chr2> ...`

The first two runs just estimate coalescence rates within each population. The third run selects only the 16 cross-population pairs to estimate the cross-coalescence rate.

Finally, you can generate a combined msmc output file using the combineCrossCoal.py script from the msmc-tools repository:
    
    ./combineCrossCoal.py across_msmc.final.txt within1_msmc.final.txt within2_msmc.final.txt > combined12_msmc.final.txt

 This script handles the fact that the three coalescence rate estimates result in different time patternings, using interpolation of the within-coalescence rates.

While this approach seems to be more cumbersome at first, you get a speed-up because you can separate the run into three parallel runs, and you can reuse the population size estimates within each population. 

Note also that with MSMC2 there is no numerical problem anymore for large number of haplotypes. Of course, the complexity still increases quite dramatically, and I don't think it's feasible to go beyond 16 haplotypes (even that will be a huge run), but the model itself is perfectly accurate up to arbitrary numbers of haplotypes. Of course, there may also be constraints by the phasing quality. With limited phasing you will get problems in very recent times, so a large number of haplotypes may not be helpful without good phasing.

