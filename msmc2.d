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
import std.c.stdlib;
import std.range;
import model.data;
import model.psmc_model;
import expectation_step;
import maximization_step;
import logger;
import model.time_intervals;

auto maxIterations = 20UL;
double mutationRate;
double recombinationRate;
auto timeSegmentPattern = [1UL, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
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
size_t[] indices;
string logFileName, loopFileName, finalFileName;
double time_factor = 0.1;


auto helpString = "Usage: msmc [options] <datafiles>
  Options:
    -i, --maxIterations=<size_t> : number of EM-iterations [default=20]
    -o, --outFilePrefix=<string> : file prefix to use for all output files
    -m, --mutationRate=<double> : mutation rate, scaled by 2N. In case of more than two haplotypes, this needs to be 
          the same as was used in running \"msmc branchlength\".
    -r, --recombinationRate=<double> : recombination rate, scaled by 2N, to begin with
          [by default set to mutationRate / 4]. Recombination rate inference does not work very well for 
          more than two haplotypes. Using the -R option is recommended for more than 2 haplotypes.
    -t, --nrThreads=<size_t> : nr of threads to use (defaults to nr of CPUs)
    -p, --timeSegmentPattern=<string> : pattern of fixed time segments [default=10*1+15*2]
    -R, --fixedRecombination : keep recombination rate fixed [recommended, but not set by default]
    -v, --verbose: write out the expected number of transition matrices (into a separate file)
    -I, --indices: indices (comma-separated) of alleles in the data file to run over
    --skipAmbiguous: skip sites with ambiguous phasing. Recommended for gene flow analysis
    --hmmStrideWidth <int> : stride width to traverse the data in the expectation step [default=1000]
    --initialLambdaVec <str> : comma-separated string of lambda-values to start with. This can be used to
      continue a previous run by copying the values in the last row and the third column of the corresponding
      *.loop file
    --time_factor <double>: factor to reduce or extent the time intervals";

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
    indices = std.string.split(value, ",").map!"a.to!size_t()"().array();
  }
  
  if(args.length == 1) {
    displayHelpMessageAndExit();
  }

  getopt(args,
      std.getopt.config.caseSensitive,
      "maxIterations|i", &maxIterations,
      "mutationRate|m", &mutationRate,
      "recombinationRate|r", &recombinationRate,
      "timeSegmentPattern|p", &handleTimeSegmentPatternString,
      "nrThreads|t", &nrThreads,
      "verbose|v", &verbose,
      "outFilePrefix|o", &outFilePrefix,
      "indices|I", &handleIndices,
      "skipAmbiguous", &skipAmbiguous,
      "help|h", &displayHelpMessageAndExit,
      "hmmStrideWidth", &hmmStrideWidth,
      "fixedRecombination|R", &fixedRecombination,
      "initialLambdaVec", &handleLambdaVecString,
      "treeFileNames", &handleTreeFileNames,
      "time_factor", &time_factor
  );
  if(nrThreads)
    std.parallelism.defaultPoolThreads(nrThreads);
  enforce(args.length > 1, "need at least one input file");
  enforce(hmmStrideWidth > 0, "hmmStrideWidth must be positive");
  inputFileNames = args[1 .. $];
  if(indices.length == 0) {
    auto nrHaplotypes = getNrHaplotypesFromFile(inputFileNames[0]);
    indices = iota(nrHaplotypes).array();
  }
  inputData = readDataFromFiles(inputFileNames, indices, skipAmbiguous);
  if(isNaN(mutationRate)) {
    stderr.write("estimating mutation rate: ");
    mutationRate = getTheta(inputData) / 2.0;
    stderr.writeln(mutationRate);
  }
  if(isNaN(recombinationRate))
    recombinationRate = mutationRate / 4.0;
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
  logInfo(format("input files:         %s\n", inputFileNames));
  logInfo(format("maxIterations:       %s\n", maxIterations));
  logInfo(format("mutationRate:        %s\n", mutationRate));
  logInfo(format("recombinationRate:   %s\n", recombinationRate));
  logInfo(format("timeSegmentPattern:  %s\n", timeSegmentPattern));
  logInfo(format("nrThreads:           %s\n", nrThreads == 0 ? totalCPUs : nrThreads));
  logInfo(format("verbose:             %s\n", verbose));
  logInfo(format("outFilePrefix:       %s\n", outFilePrefix));
  logInfo(format("hmmStrideWidth:      %s\n", hmmStrideWidth));
  logInfo(format("fixedRecombination:  %s\n", fixedRecombination));
  logInfo(format("initialLambdaVec:    %s\n", lambdaVec));
  logInfo(format("skipAmbiguous:       %s\n", skipAmbiguous));
  logInfo(format("indices:             %s\n", indices));
  logInfo(format("logging information written to %s\n", logFileName));
  logInfo(format("loop information written to %s\n", loopFileName));
  logInfo(format("final results written to %s\n", finalFileName));
  logInfo(format("time factor:         %s\n", time_factor));
  if(verbose)
    logInfo(format("transition matrices written to %s.loop_*.expectationMatrix.txt\n", outFilePrefix));
}

void run() {
  PSMCmodel params;
  auto timeIntervals = TimeIntervals.standardIntervals(nrTimeSegments, time_factor);
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
    auto lambdaCI = getLambdaCI(transitions, emissions, params, timeSegmentPattern);
    auto recCI = getRecombinationCI(transitions, emissions, params, timeSegmentPattern);
    
    printLoop(loopFileName, params, logLikelihood, recCI, lambdaCI);
    if(verbose) {
      auto filename = outFilePrefix ~ format(".loop_%s.expectationMatrix.txt", iteration);
      printMatrix(filename, transitions, emissions);
    }
  }
  
  printFinal(finalFileName, params);
}

SegSite_t[][] readDataFromFiles(string[] filenames, size_t[] indices, bool skipAmbiguous) {
  SegSite_t[][] ret;
  foreach(filename; filenames) {
    foreach(i; 0 .. indices.length - 1) {
      foreach(j; i + 1 .. indices.length) {
        size_t[2] ind_pair = [indices[i], indices[j]];
        auto data = readSegSites(filename, ind_pair, skipAmbiguous);
        logInfo(format("read %s SNPs from file %s, using indices %s\n", data.length, filename, ind_pair));
        ret ~= data;
      }
    }
  }
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

void printLoop(string filename, PSMCmodel params, double logLikelihood, double[2] recCI, double[2][] lambdaCI) {
  auto f = File(filename, "a");
  f.writef("%s,%s,%s\t%.2f\t%s\t", params.recombinationRate, recCI[0], recCI[1], logLikelihood, 
           params.timeIntervals.boundaries.map!(a => text(a * params.mutationRate)).join(",").array());
  foreach(i; 0 .. lambdaCI.length) {
    f.writef("%s,%s,%s", params.lambdaVec[i] / mutationRate, lambdaCI[i][0] / mutationRate, lambdaCI[i][1] / mutationRate);
    if(i < lambdaCI.length - 1)
      f.write(",");
  }
  f.write("\n");
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
