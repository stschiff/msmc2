/* Copyright (c) 2012,2013 Genome Research Ltd.
 *
 * Author: Stephan Schiffels <stephan.schiffels@sanger.ac.uk>
 *
 * This file is part of msmc.
 * msmc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
import std.stdio;
import std.math;
import std.string;
import std.conv;
import std.getopt;
import std.parallelism;
import std.algorithm;
import std.array;
import std.file;
import std.typecons;
import std.regex;
import std.exception;
import core.stdc.stdlib;
import std.range;
import model.data;
import model.psmc_model;
import expectation_step;
import maximization_step;
import logger;
import model.time_intervals;
import core.memory;

auto maxIterations = 20UL;
double mutationRate;
double recombinationRate;
// auto timeSegmentPattern = [4UL, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 6];
auto timeSegmentPattern = [2UL, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3];
uint nrThreads;
auto verbose = false;
string outFilePrefix;
auto memory = false;
auto fixedRecombination = false;
bool skipAmbiguous = false;
string[] inputFileNames, treeFileNames;
SegSite_t[][] inputData;
size_t hmmStrideWidth = 1000;
double[] lambdaVec;
size_t nrTimeSegments;
size_t[2][] pairIndices;
string logFileName, loopFileName, finalFileName;
double time_factor = 1.0;
bool quantileBoundaries = false;

immutable versionString = "2.1.0";


auto helpString = format("This is version %s. Usage: msmc2 [options] <datafiles>
  Options:
    -i, --maxIterations=<size_t> :      number of EM-iterations [default=20]
    -o, --outFilePrefix=<string> :      file prefix to use for all output files
    -r, --rhoOverMu=<double> :          initial ratio of recombination over mutation rate (default: 
                                        0.25)
    -t, --nrThreads=<size_t> :          nr of threads to use (defaults to nr of CPUs)
    -p, --timeSegmentPattern=<string> : pattern of fixed time segments [default=1*2+25*1+1*2+1*3]
    -R, --fixedRecombination :          keep recombination rate fixed (rarely needed in MSMC2)
    -I, --pairIndices:                  this can be given in two flavors. First, you can enter a 
                                        single comma-separated list like this \"-I 0,1,4,5\".
                                        In this case, the program will 
                                        run over all pairs of haplotypes within this set of 
                                        indices. This is useful for running on multiple phased 
                                        diploid genomes sampled from one population.
                                        In the second flavor, you can give a list of pairs, like 
                                        this: \"-I 0-1,2-3,4-5\". In this case, the 
                                        program will run only those specified pairs. This can be 
                                        used to run on a number of unphased genomes, to avoid pairs 
                                        of haplotypes from different individuals. This should also 
                                        be used to indicate a cross-population run, where you want 
                                        to run the program only over pairs of haplotypes across 
                                        population boundaries. So with two phased genomes, one from 
                                        each population you'd run \"-I 0-2,0-3,1-2,1-3\", and the 
                                        program would run only those four pairs of haplotypes.
                                        Note that if you do not use this parameter altogether, 
                                        MSMC2 will run on all pairs of input haplotypes.
    -s, --skipAmbiguous:                skip sites with ambiguous phasing. Recommended for cross 
                                        population analysis
    --quantileBoundaries:               use quantile boundaries, as in MSMC. To fully replicate 
                                        MSMC's time intervals, combine this with -p 10*1+15*2",     versionString);

void main(string[] args) {
  try {
    parseCommandLine(args);
  }
  catch(Exception e) {
    stderr.writeln("error in parsing command line: ", e);
    exit(0);
  }
  run();
}

void parseCommandLine(string[] args) {
  
  void displayHelpMessageAndExit() {
    stderr.writeln(helpString);
    exit(0);
  }
  void handleTimeSegmentPatternString(string option, string patternString) {
    enforce(match(patternString, r"^\d+\*\d+[\+\d+\*\d+]*"), text("illegal timeSegmentPattern: ", patternString));
    timeSegmentPattern.length = 0;
    foreach(product; std.string.split(patternString, "+")) {
      auto pair = array(map!"to!size_t(a)"(std.string.split(product, "*")));
      foreach(i; 0 .. pair[0]) {
        timeSegmentPattern ~= pair[1];
      }
    }
  }

  void handleLambdaVecString(string option, string lambdaString) {
    enforce(match(lambdaString, r"^[\d.]+[,[\d.]+]+"), text("illegal array string: ", lambdaString));
    lambdaVec = std.string.split(lambdaString, ",").map!"to!double(a)"().array();
  }
  
  void handleTreeFileNames(string option, string value) {
    treeFileNames = std.string.split(value, ",").array();
  }
  
  void handleIndices(string option, string value) {
    try {
      auto hapIndices = std.string.split(value, ",").map!"a.to!size_t()"().array();
      foreach(i; hapIndices[0..$-1])
        foreach(j; hapIndices[i..$])
          pairIndices ~= [i, j];
    }
    catch(ConvException e) {
      auto pairIndexSubStrings = std.string.split(value, ",");
      foreach(subStr; pairIndexSubStrings) {
        auto indexPair = std.string.split(subStr, "-").map!"a.to!size_t()"().array();
        enforce(indexPair.length == 2, "cannot parse pairIndices");
        pairIndices ~= [indexPair[0], indexPair[1]];
      }
    }
  }

  if(args.length == 1) {
    displayHelpMessageAndExit();
  }

  auto rhoOverMu = 0.25;
  getopt(args,
      std.getopt.config.caseSensitive,
      "maxIterations|i", &maxIterations,
      "theta|m", &mutationRate,
      "rhoOverMu|r", &rhoOverMu,
      "timeSegmentPattern|p", &handleTimeSegmentPatternString,
      "nrThreads|t", &nrThreads,
      "verbose", &verbose,
      "outFilePrefix|o", &outFilePrefix,
      "pairIndices|I", &handleIndices,
      "skipAmbiguous|s", &skipAmbiguous,
      "help|h", &displayHelpMessageAndExit,
      "hmmStrideWidth", &hmmStrideWidth,
      "fixedRecombination", &fixedRecombination,
      "initialLambdaVec", &handleLambdaVecString,
      "treeFileNames", &handleTreeFileNames,
      "time_factor", &time_factor,
      "quantileBoundaries", &quantileBoundaries
  );
  if(nrThreads)
    std.parallelism.defaultPoolThreads(nrThreads);
  // auto reserved = GC.reserve(20_000_000_000);
  // stderr.writefln("reserved %s bytes for Garbage Collector", reserved);
  enforce(args.length > 1, "need at least one input file");
  enforce(hmmStrideWidth > 0, "hmmStrideWidth must be positive");
  inputFileNames = args[1 .. $];
  if(pairIndices.length == 0) {
    auto nrHaplotypes = getNrHaplotypesFromFile(inputFileNames[0]);
    foreach(i; 0..(nrHaplotypes - 1))
      foreach(j; (i+1)..nrHaplotypes)
        pairIndices ~= [i, j];
  }
  inputData = readDataFromFiles(inputFileNames, pairIndices, skipAmbiguous);
  if(isNaN(mutationRate)) {
    stderr.write("estimating mutation rate: ");
    mutationRate = getTheta(inputData) / 2.0;
    stderr.writeln(mutationRate);
  }
  recombinationRate = mutationRate * rhoOverMu;
  nrTimeSegments = timeSegmentPattern.reduce!"a+b"();
  if(lambdaVec.length > 0) {
    // this is necessary because we read in a scaled lambdaVec.
    lambdaVec[] *= mutationRate;
    enforce(lambdaVec.length == nrTimeSegments, "initialLambdaVec must have correct length");
  }
  enforce(treeFileNames.length == 0 || treeFileNames.length == inputFileNames.length);
  
  logFileName = outFilePrefix ~ ".log";
  loopFileName = outFilePrefix ~ ".loop.txt";
  finalFileName = outFilePrefix ~ ".final.txt";
  logger.logFile = File(logFileName, "w");
  
  printGlobalParams();
}

void printGlobalParams() {
  logInfo(format("Version:                       %s\n", versionString));
  logInfo(format("input files:                   %s\n", inputFileNames));
  logInfo(format("maxIterations:                 %s\n", maxIterations));
  logInfo(format("mutationRate:                  %s\n", mutationRate));
  logInfo(format("recombinationRate:             %s\n", recombinationRate));
  logInfo(format("timeSegmentPattern:            %s\n", timeSegmentPattern));
  logInfo(format("nrThreads:                     %s\n", nrThreads == 0 ? totalCPUs : nrThreads));
  logInfo(format("outFilePrefix:                 %s\n", outFilePrefix));
  logInfo(format("hmmStrideWidth:                %s\n", hmmStrideWidth));
  logInfo(format("fixedRecombination:            %s\n", fixedRecombination));
  logInfo(format("initialLambdaVec:              %s\n", lambdaVec));
  logInfo(format("skipAmbiguous:                 %s\n", skipAmbiguous));
  logInfo(format("pairIndices:                   %s\n", pairIndices));
  logInfo(format("logging information written to %s\n", logFileName));
  logInfo(format("loop information written to    %s\n", loopFileName));
  logInfo(format("final results written to       %s\n", finalFileName));
  logInfo(format("time factor:                   %s\n", time_factor));
  if(verbose)
    logInfo(format("transition matrices written to %s.loop_*.expectationMatrix.txt\n", outFilePrefix));
}

void run() {
  PSMCmodel params;
  auto nrPairs = pairIndices.length;
  auto time_constant = time_factor * 0.1 / to!double(nrPairs);
  auto timeIntervals = TimeIntervals.standardIntervals(nrTimeSegments, time_constant);
  if(quantileBoundaries) {
    time_constant = time_factor / to!double(nrPairs);
    auto b = TimeIntervals.getQuantileBoundaries(nrTimeSegments, time_constant);
    timeIntervals = new TimeIntervals(b);
  }
  if(lambdaVec.length == 0)
    params = new PSMCmodel(mutationRate, recombinationRate, timeIntervals);
  else
    params = new PSMCmodel(mutationRate, recombinationRate, lambdaVec, timeIntervals);
  
  auto nrFiles = inputData.length;

  auto f = File(loopFileName, "w");
  f.close();
  
  foreach(iteration; 0 .. maxIterations) {
    logInfo(format("[%s/%s] Baumwelch iteration\n", iteration + 1, maxIterations));
    auto expectationResult = getExpectation(inputData, params, hmmStrideWidth, 1000);
    auto transitions = expectationResult[0];
    auto emissions = expectationResult[1];
    auto logLikelihood = expectationResult[2];
    
    auto newParams = getMaximization(transitions, emissions, params, timeSegmentPattern, fixedRecombination);
    params = newParams;
    // auto lambdaCI = getLambdaCI(transitions, emissions, params, timeSegmentPattern);
    // auto recCI = getRecombinationCI(transitions, emissions, params, timeSegmentPattern);
    
    // printLoop(loopFileName, params, logLikelihood, recCI, lambdaCI);
    printLoop(loopFileName, params, logLikelihood);
    if(verbose) {
      auto filename = outFilePrefix ~ format(".loop_%s.expectationMatrix.txt", iteration);
      printMatrix(filename, transitions, emissions);
    }
  }
  
  printFinal(finalFileName, params);
}

SegSite_t[][] readDataFromFiles(string[] filenames, size_t[2][] indexPairs, bool skipAmbiguous) {
    SegSite_t[][] ret;
    
    GC.disable();
    foreach(i, filename; filenames) {
        auto data = readSegSites(filename, indexPairs, skipAmbiguous);
        logInfo(format("read %s SNPs from file %s, using indices %s\n", data[0].length, filename, 
          indexPairs));
        ret ~= data;
        if(i % 10 == 0) {
            GC.enable();
            GC.collect();
            GC.minimize();
            GC.disable();
        }
    }
    GC.enable();
    GC.collect();
    GC.minimize();
    return ret;
}

void printMatrix(string filename, double[][] transitions, double[][2] emissions) {
  auto f = File(filename, "w");
  foreach(a; 0 .. transitions.length) {
    f.writeln(transitions[a].map!"text(a)"().join("\t"));
  }
  f.writeln(emissions[0].map!"text(a)"().join("\t"));
  f.writeln(emissions[1].map!"text(a)"().join("\t"));
}

// void printLoop(string filename, PSMCmodel params, double logLikelihood, double[2] recCI, double[2][] lambdaCI) {
//   auto f = File(filename, "a");
//   f.writef("%s,%s,%s\t%.2f\t%s\t", params.recombinationRate, recCI[0], recCI[1], logLikelihood,
//            params.timeIntervals.boundaries.map!(a => text(a * params.mutationRate)).join(",").array());
//   foreach(i; 0 .. lambdaCI.length) {
//     f.writef("%s,%s,%s", params.lambdaVec[i] / mutationRate, lambdaCI[i][0] / mutationRate, lambdaCI[i][1] / mutationRate);
//     if(i < lambdaCI.length - 1)
//       f.write(",");
//   }
//   f.write("\n");
// }

void printLoop(string filename, PSMCmodel params, double logLikelihood) {
  auto f = File(filename, "a");
  f.writef("%s\t%.2f\t%s\t%s\n", params.recombinationRate, logLikelihood, 
           params.timeIntervals.boundaries.map!(a => text(a * params.mutationRate)).join(",").array(),
           params.lambdaVec.map!(a => text(a / params.mutationRate)).join(",").array());
}

void printFinal(string filename, PSMCmodel params) { 
  auto f = File(filename, "w");
  f.writeln("time_index\tleft_time_boundary\tright_time_boundary\tlambda");
  foreach(i; 0 .. params.nrStates) {
    auto left = params.timeIntervals.leftBoundary(i);
    auto right = params.timeIntervals.rightBoundary(i);
    f.writefln("%s\t%s\t%s\t%s", i, left * mutationRate, right * mutationRate, params.lambdaVec[i] / mutationRate);
  }
}
