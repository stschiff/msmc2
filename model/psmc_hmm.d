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
 
module model.psmc_hmm;
import std.stdio;
import std.math;
import std.string;
import std.conv;
import std.algorithm;
import std.concurrency;
import std.typecons;
import std.random;
import std.exception;
import model.data;
import model.propagation_core;
import model.stateVec;
import model.stateVecAllocator;
import model.gsl_matrix_vector;

SegSite_t[] chop_segsites(in SegSite_t[] segsites, size_t maxDistance) {
  SegSite_t[] ret;
  size_t lastPos = 0;
  foreach(segsite; segsites) {
    while(segsite.pos - lastPos > maxDistance) {
      ret ~= new SegSite_t(lastPos + maxDistance, min(segsite.obs[0], 1));
      lastPos += maxDistance;
    }
    ret ~= segsite.dup;
    lastPos = segsite.pos;
  }
  return ret;
}

unittest {
  writeln("testing chop_sites");
  auto data = [
    new SegSite_t(400, [2]),
    new SegSite_t(3600, [0]), // missing data
    new SegSite_t(5000, [2])
  ];
  auto ret = chop_segsites(data, 1000);
  assert(ret[0].pos == 400);
  assert(ret[0].obs == [2]);
  assert(ret[1].pos == 1400);
  assert(ret[1].obs == [0]);
  assert(ret[3].pos == 3400);
  assert(ret[3].obs == [0]);
  assert(ret[5].pos == 4600);
  assert(ret[5].obs == [1]); // homozygous
  assert(ret[6].pos == 5000);
  assert(ret[6].obs == [2]);
}

class PSMC_hmm {

  size_t L;
  size_t maxDistance;
  size_t indexCache;
  size_t hmmStrideWidth;
  size_t currentBackwardIndex;
  State_t currentBackwardState, nextBackwardState;
  
  const SegSite_t[] segsites;
  double[] scalingFactors;
  const PropagationCore propagationCore;
  StateVecAllocator stateVecAllocator;
  State_t[] forwardStates;
  State_t expectationForwardDummy, expectationBackwardDummy, getBackwardStateDummy;
  State_t runBackwardDummy;
  bool have_run_forward;
  
  this(in PropagationCore propagationCore, in SegSite_t[] segsites) {
    this.propagationCore = propagationCore;
    this.maxDistance = propagationCore.maxDistance;
    this.segsites = chop_segsites(segsites, maxDistance);
    this.L = this.segsites.length;
    this.hmmStrideWidth = hmmStrideWidth;
    
    scalingFactors = new double[L];
    scalingFactors[] = 0.0;
    
    auto stateSize = propagationCore.forwardStateSize;
    stateVecAllocator = new StateVecAllocator(L * stateSize);
    forwardStates = new State_t[L];
    foreach(i; 0 .. L)
      forwardStates[i] = propagationCore.newForwardState(stateVecAllocator);

    expectationForwardDummy = propagationCore.newForwardState();
    expectationBackwardDummy = propagationCore.newBackwardState();
    getBackwardStateDummy = propagationCore.newBackwardState();
    runBackwardDummy = propagationCore.newBackwardState();
    currentBackwardState = propagationCore.newBackwardState();
    nextBackwardState = propagationCore.newBackwardState();
    currentBackwardIndex = L - 1;
  }

  double logLikelihood() const {
    return scalingFactors.map!log().reduce!"a+b"();
  }
  
  void runForward()
  out {
    foreach(i; 0 .. L)
      assert(scalingFactors[i] > 0);
  }
  body {
    enforce(!have_run_forward);
    propagationCore.initialState(forwardStates[0]);
    scalingFactors[0] = forwardStates[0].norm;
    forwardStates[0].scale(1.0 / scalingFactors[0]);

    auto forwardDummyVec = propagationCore.newForwardState();
    foreach(index; 1 .. L) {
      
      if(segsites[index].pos == segsites[index - 1].pos + 1) {
        propagationCore.propagateSingleForward(forwardStates[index - 1],
            forwardStates[index], segsites[index - 1], segsites[index]);
      }
      else {
        auto dummy_site = getSegSite(segsites[index].pos - 1);
        propagationCore.propagateMultiForward(forwardStates[index - 1], 
            forwardDummyVec, segsites[index - 1], dummy_site);
        propagationCore.propagateSingleForward(forwardDummyVec, forwardStates[index],
            dummy_site, segsites[index]);
      }
      scalingFactors[index] = forwardStates[index].norm;
      assert(scalingFactors[index] > 0.0, text(scalingFactors));
      forwardStates[index].scale(1.0 / scalingFactors[index]);
    }
    have_run_forward = true;
  }

  Tuple!(double[][], double[][2]) runBackward(size_t hmmStrideWidth=1000) {
    enforce(have_run_forward);

    auto nrStates = propagationCore.getPSMC.nrStates;

    auto transitions = new double[][](nrStates, nrStates);
    double[][2] emissions = [new double[nrStates], new double[nrStates]];
    foreach(i; 0 .. nrStates)
      transitions[i][] = 0.0;
    foreach(i; 0 .. 2)
      emissions[i][] = 0.0;
    
    currentBackwardIndex = L - 1;
    auto transitionsDummy = new double[][](nrStates, nrStates);
    double[][2] emissionsDummy = [new double[nrStates], new double[nrStates]];
    for(size_t pos = segsites[$ - 1].pos; pos > segsites[0].pos && pos <= segsites[$ - 1].pos; pos -= hmmStrideWidth) {
      getForwardState(expectationForwardDummy, pos - 1);
      getBackwardState(expectationBackwardDummy, pos);
      auto site = getSegSite(pos);
    
      propagationCore.getTransitionExpectation(expectationForwardDummy, expectationBackwardDummy, site, transitionsDummy);
      
      getForwardState(expectationForwardDummy, pos);
      propagationCore.getEmissionExpectation(expectationForwardDummy, expectationBackwardDummy, site, emissionsDummy);

      foreach(i; 0 .. nrStates)
        transitions[i][] += transitionsDummy[i][];
      foreach(i; 0 .. 2)
        emissions[i][] += emissionsDummy[i][];
    }

    return tuple(transitions, emissions);
  }
    
  void getForwardState(State_t s, size_t pos)
  in {
    assert(pos >= segsites[0].pos);
    assert(pos <= segsites[$ - 1].pos);
  }
  body {
    auto index = getRightIndexAtPos(pos);
    if(pos == segsites[index].pos) {
      forwardStates[index].copy_into(s);
    }
    else {
      auto site = getSegSite(pos);
      propagationCore.propagateMultiForward(forwardStates[index - 1], s, segsites[index - 1], site);
    }
  }

  void getBackwardState(State_t s, size_t pos)
  in {
    assert(pos >= segsites[0].pos);
    assert(pos <= segsites[$ - 1].pos);
  }
  body {
    auto index = getRightIndexAtPos(pos);
    auto site = getSegSite(pos);
    if(pos == segsites[index].pos) {
      getBackwardStateAtIndex(index).copy_into(s);
    }
    else {
      if(pos == segsites[index].pos - 1) {
        propagationCore.propagateSingleBackward(getBackwardStateAtIndex(index), s, segsites[index], site);
      }
      else {
        auto dummy_site = getSegSite(segsites[index].pos - 1);
        propagationCore.propagateSingleBackward(getBackwardStateAtIndex(index), getBackwardStateDummy, segsites[index], dummy_site);
        propagationCore.propagateMultiBackward(getBackwardStateDummy, s, dummy_site, site);
      }
    }
  }
  
  private SegSite_t getSegSite(size_t pos) {
    auto index = getRightIndexAtPos(pos);
    if(segsites[index].pos == pos)
      return segsites[index].dup;
    else
      return new SegSite_t(pos, min(segsites[index].obs[0], 1));
  }
  
  private size_t getRightIndexAtPos(size_t pos)
  in {
    assert(pos <= segsites[L - 1].pos);
  }
  out(result) {
    assert(segsites[result].pos >= pos);
    if(result > 0) {
      assert(segsites[result - 1].pos < pos);
    }
  }
  body {
    while(segsites[indexCache].pos < pos) {
      ++indexCache;
    }
    if(indexCache > 0) {
      while(segsites[indexCache - 1].pos >= pos) {
        if(--indexCache == 0)
            break;
      }
    }
    return indexCache;
  }

  private State_t getBackwardStateAtIndex(size_t index)
  in {
    assert(have_run_forward);
    assert(index < L);
  }
  body {
    if(index == L - 1) {
      assert(scalingFactors[L - 1] > 0, text(scalingFactors[L - 1]));
      propagationCore.setState(currentBackwardState, 1.0 / scalingFactors[L - 1], segsites[L - 1]);
      currentBackwardIndex = L - 1;
    }
    else {
      assert(index <= currentBackwardIndex, text([index, L]));
      while(index < currentBackwardIndex) {
        if(segsites[currentBackwardIndex].pos == segsites[currentBackwardIndex - 1].pos + 1) {
          propagationCore.propagateSingleBackward(currentBackwardState, nextBackwardState,
              segsites[currentBackwardIndex], segsites[currentBackwardIndex - 1]);
        }
        else {
          auto dummy_site = getSegSite(segsites[currentBackwardIndex].pos - 1);
          propagationCore.propagateSingleBackward(currentBackwardState, runBackwardDummy,
              segsites[currentBackwardIndex], dummy_site);
          propagationCore.propagateMultiBackward(runBackwardDummy, nextBackwardState,
              dummy_site, segsites[currentBackwardIndex - 1]);
        }
        nextBackwardState.copy_into(currentBackwardState);
        currentBackwardState.scale(1.0 / scalingFactors[currentBackwardIndex - 1]);
        currentBackwardIndex -= 1;
      }
    }
    return currentBackwardState;
  }

}

unittest {
  writeln("testing MSMC_hmm");
  import model.psmc_model;
  
  auto T = 10UL;
  auto params = new PSMCmodel(0.01, 0.001, T);
  auto lvl = 1.0e-8;
  
  auto propagationCore = new PropagationCore(params, 100);

  auto nrS = propagationCore.getPSMC.nrStates;
  
  auto data = readSegSites("model/hmm_testData.txt", [[0UL, 1]], false)[0];
  
  auto psmc_hmm = new PSMC_hmm(propagationCore, data);
  psmc_hmm.runForward();
  
  auto L = psmc_hmm.L;
  
  for(auto pos = L - 1; pos >= 0 && pos < L; --pos) {
    auto sum = 0.0;
    foreach(a; 0 .. T)
      sum += gsl_vector_get(psmc_hmm.forwardStates[pos].vec, a) *
        gsl_vector_get(psmc_hmm.getBackwardStateAtIndex(pos).vec, a) * 
        psmc_hmm.scalingFactors[pos];
    
    assert(approxEqual(sum, 1.0, lvl, 0.0), text(sum));
  }
}


