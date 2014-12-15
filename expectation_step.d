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
 
import std.typecons;
import std.stdio;
import std.string;
import std.algorithm;
import std.exception;
import std.conv;
import std.parallelism;
import std.math;
import core.memory;
import model.propagation_core;
import model.psmc_model;
import model.psmc_hmm;
import model.data;
import logger;

alias Tuple!(double[][], double[][2], double) ExpectationResult_t;

ExpectationResult_t getExpectation(in SegSite_t[][] inputData, PSMCmodel psmc, size_t hmmStrideWidth, size_t maxDistance)
{
  PropagationCore propagationCore;
  propagationCore = new PropagationCore(psmc, maxDistance);
  
  auto transitions = new double[][](psmc.nrStates, psmc.nrStates);
  double[][2] emissions = [new double[psmc.nrStates], new double[psmc.nrStates]];
  foreach(a; 0 .. psmc.nrStates) {
    transitions[a][] = 0.0;
    emissions[0][] = 0.0;
    emissions[1][] = 0.0;
  }
  auto logLikelihood = 0.0;
  
  auto cnt = 0;
  GC.disable();
  foreach(i, data; taskPool.parallel(inputData, 1)) {
    logInfo(format("\r  * [%s/%s] Expectation Step", ++cnt, inputData.length));
    auto result = singleChromosomeExpectation(data, hmmStrideWidth, propagationCore);
    emissions[0][] += result[1][0][];
    emissions[1][] += result[1][1][];
    foreach(a; 0 .. psmc.nrStates)
      transitions[a][] += result[0][a][];
    logLikelihood += result[2];
    if(i % 500 == 0) {
        GC.enable();
        GC.collect();
        GC.minimize();
        GC.disable();
    }
  }
  GC.enable();
  GC.collect();
  GC.minimize();
  logInfo(format(", log likelihood: %s", logLikelihood));
  logInfo("\n");
  
  return tuple(transitions, emissions, logLikelihood);
}

ExpectationResult_t singleChromosomeExpectation(in SegSite_t[] data, size_t hmmStrideWidth, in PropagationCore propagationCore)
{
  auto psmc_hmm =  new PSMC_hmm(propagationCore, data);

  psmc_hmm.runForward();
  auto exp = psmc_hmm.runBackward(hmmStrideWidth);
  auto logL = psmc_hmm.logLikelihood();
  return tuple(exp[0], exp[1], logL);
}

unittest {
  writeln("test expectation step");
  auto psmc = new PSMCmodel(0.01, 0.001, 4);
  auto fileName = "model/hmm_testData.txt";
  auto data = readSegSites(fileName, [0UL, 1], false);
  auto hmmStrideWidth = 100UL;
  
  auto allData = [data, data, data];
  
  std.parallelism.defaultPoolThreads(1U);
  auto resultSingleThreaded = getExpectation(allData, psmc, hmmStrideWidth, 100UL);
  std.parallelism.defaultPoolThreads(2U);
  auto resultMultiThreaded = getExpectation(allData, psmc, hmmStrideWidth, 100UL);
   
  foreach(a; 0 .. psmc.nrStates) {
    assert(all!"a>0.0"(resultSingleThreaded[0][a]));
    assert(all!"a>0.0"(resultMultiThreaded[0][a]));
    assert(resultSingleThreaded[1][0][a] > 0.0);
    assert(resultMultiThreaded[1][1][a] > 0.0);
  }
}
