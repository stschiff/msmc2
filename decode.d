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
import std.getopt;
import std.exception;
import std.c.stdlib;
import std.algorithm;
import std.parallelism;
import std.array;
import std.conv;
import model.psmc_hmm;
import model.data;
import model.psmc_model;
import model.gsl_matrix_vector;
import model.propagation_core;

double mutationRate, recombinationRate;
size_t nrTimeSegments=64;
size_t stride=1000;
string inputFileName;
size_t nrHaplotypes;
uint nrThreads;
size_t[2] indices = [0, 1];
double[] lambdaVec;

void main(string[] args) {
  try {
    parseCommandlineArgs(args);
  }
  catch (Exception e) {
    stderr.writeln(e.msg);
    displayHelpMessage();
    exit(0);
  }
  run();
}

void parseCommandlineArgs(string[] args) {
  
  void handleIndices(string option, string value) {
    indices = std.string.split(value, ",").map!"a.to!size_t()"().array();
  }
  
  void handleLambdaVecString(string option, string lambdaString) {
    enforce(match(lambdaString, r"^[\d.]+[,[\d.]+]+"), text("illegal array string: ", lambdaString));
    lambdaVec = std.string.split(lambdaString, ",").map!"to!double(a)"().array();
  }
  
  getopt(args,
         std.getopt.config.caseSensitive,
         "mutationRate|m", &mutationRate,
         "recombinationRate|r", &recombinationRate,
         "nrTimeSegments|t", &nrTimeSegments,
         "nrThreads", &nrThreads,
         "indices|I", &handleIndices,
         "lambdaVec|l", &handleLambdaVecString,
         "stride|s", &stride);
  if(nrThreads)
    std.parallelism.defaultPoolThreads(nrThreads);
  enforce(args.length == 2, "need exactly one input file");
  enforce(indices.length == 2, "need exactly two indices");
  inputFileName = args[1];
  nrHaplotypes = getNrHaplotypesFromFile(inputFileName);
  enforce(mutationRate > 0, "need positive mutationRate");
  enforce(recombinationRate > 0, "need positive recombinationRate");
  if(lambdaVec.length == 0) {
    lambdaVec = new double[nrTimeSegments];
    lambdaVec[] = 1.0;
  }
  else {
    nrTimeSegments = lambdaVec.length;
    enforce(nrTimeSegments > 0, "number of time segments is zero");
  }
}

void displayHelpMessage() {
  stderr.writeln("Usage: msmc decode [options] <inputFile>
Options:
-m, --mutationRate <double>
-r, --recombinationRate <double>
-t, --nrTimeSegments <int> (ignored when lambdaVec is given)
-l, --lambdaVec <list> (comma-separated list of scaled inverse population sizes 1/2N(t). By default all set to 1)
-I, --indices <int,int> (two comma separated numbers to select the two haplotypes from the data file, by default: 0,1)
--nrThreads <int> : nr of threads, defaults to nr of CPUs
-s, --stride <int>: stride width in basepairs [default=1000]");
}

void run() {
  auto hmm = makeStandardHmm();
  decodeWithHmm(hmm);
}

PSMC_hmm makeStandardHmm() {
  auto model = new PSMCmodel(mutationRate, recombinationRate, lambdaVec);

  stderr.writeln("generating propagation core");
  auto propagationCore = new PropagationCore(model, 1000);
  
  auto segsites = readSegSites(inputFileName, [indices], false);

  stderr.writeln("generating Hidden Markov Model");
  return new PSMC_hmm(propagationCore, segsites[0]);
}

void decodeWithHmm(PSMC_hmm hmm) {
  hmm.runForward();
  auto forwardState = hmm.propagationCore.newForwardState();
  auto backwardState = hmm.propagationCore.newBackwardState();
  
  double[][] posteriors;
  for(size_t pos = hmm.segsites[$ - 1].pos; pos > hmm.segsites[0].pos && pos <= hmm.segsites[$ - 1].pos; pos -= stride)
  {
    hmm.getForwardState(forwardState, pos);
    hmm.getBackwardState(backwardState, pos);
    auto posterior = new double[nrTimeSegments];
    foreach(i; 0 .. nrTimeSegments)
      posterior[i] = gsl_vector_get(forwardState.vec, i) * gsl_vector_get(backwardState.vec, i);
    auto norm = posterior.reduce!"a+b"();
    posterior[] /= norm;
    posteriors ~= posterior;
  }
  
  foreach_reverse(posterior; posteriors) {
    foreach(p; posterior)
      write(p, " ");
    writeln("");
  } 
}
