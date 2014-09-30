/* Copyright (c) 2012,2013 Genome Research Ltd.
 *
 * Author: Stephan Schiffels <stephan.schiffels@sanger.ac.uk>
 *
 * This file is part of psmc.
 * psmc is free software: you can redistribute it and/or modify it under
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
 
module model.propagation_core;
import std.stdio;
import std.algorithm;
import std.conv;
import std.string;
import std.exception;
import model.data;
import model.propagation_core;
import model.psmc_model;
import model.stateVec;
import model.stateVecAllocator;

class PropagationCore {
  
  double[][][] forwardPropagators, backwardPropagators;
  double[][][] forwardPropagatorsMissing, backwardPropagatorsMissing;
  
  double[][] emissionProbs; // first index: second index: obs
  double[][] transitionMatrix;
  
  const PSMCmodel psmc;
  
  this(in PSMCmodel psmc, size_t maxDistance) {
    enforce(maxDistance > 0);
    this.psmc = psmc;

    forwardPropagators = new double[][][](maxDistance, psmc.nrStates, psmc.nrStates);
    backwardPropagators = new double[][][](maxDistance, psmc.nrStates, psmc.nrStates);
    forwardPropagatorsMissing = new double[][][](maxDistance, psmc.nrStates, psmc.nrStates);
    backwardPropagatorsMissing = new double[][][](maxDistance, psmc.nrStates, psmc.nrStates);
    emissionProbs = new double[][](3, psmc.nrStates);
    transitionMatrix = new double[][](psmc.nrStates, psmc.nrStates);
      
    foreach(i; 0 .. 3)
      foreach(a; 0 .. psmc.nrStates)
        emissionProbs[i][a] = psmc.emissionProb(i, a);
    
    foreach(a; 0 .. psmc.nrStates)
      foreach(b; 0 .. psmc.nrStates)
        transitionMatrix[a][b] = psmc.transitionProb(a, b);
      
    computeForwardPropagators(forwardPropagators, false, maxDistance);
    computeBackwardPropagators(backwardPropagators, false, maxDistance);

    computeForwardPropagators(forwardPropagatorsMissing, true, maxDistance);
    computeBackwardPropagators(backwardPropagatorsMissing, true, maxDistance);
    
  }
  
  private void computeForwardPropagators(double[][][] ret, bool missing_data, size_t maxDistance) const
  {
    foreach(a; 0 .. psmc.nrStates) {
      double e = missing_data ? 1.0 : emissionProbs[1][a];
      foreach(b; 0 .. psmc.nrStates)
        ret[0][a][b] = transitionMatrix[a][b] * e;
    }

    foreach(distance; 1 .. maxDistance) {
      foreach(a; 0 .. psmc.nrStates){ 
        foreach(b; 0 .. psmc.nrStates) {
          ret[distance][a][b] = 0.0;
          foreach(c; 0 .. psmc.nrStates)
            ret[distance][a][b] += ret[distance - 1][c][b] * ret[0][a][c];
        }
      }
    }
  }
  
  private void computeBackwardPropagators(double[][][] ret, bool missing_data, size_t maxDistance) const
  {
    foreach(a; 0 .. psmc.nrStates) {
      double e = missing_data ? 1.0 : emissionProbs[1][a];
      foreach(b; 0 .. psmc.nrStates)
        ret[0][a][b] = transitionMatrix[a][b] * e;
    }

    foreach(distance; 1 .. maxDistance) {
      foreach(a; 0 .. psmc.nrStates){ 
        foreach(b; 0 .. psmc.nrStates) {
          ret[distance][a][b] = 0.0;
          foreach(c; 0 .. psmc.nrStates)
            ret[distance][a][b] += ret[distance - 1][a][c] * ret[0][c][b];
        }
      }
    }
  }
  
  private double fullE(in SegSite_t segsite, size_t a) const {
    double ret = 0.0;
    foreach(o; segsite.obs) {
      ret += emissionProbs[o][a];
    }
    ret /= cast(double)segsite.obs.length;
    return ret;
  }
  
  void propagateSingleForward(in State_t from, State_t to, 
        in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  body {
    to.setZero();
    
    foreach(a; 0 .. psmc.nrStates) {
      auto sum = 0.0;
      foreach(b; 0 .. psmc.nrStates) {
        sum += from.vec[b] * transitionMatrix[a][b];
      }
      to.vec[a] = fullE(to_segsite, a) * sum;
    }
  }
  
  void propagateSingleBackward(in State_t to, State_t from,
            in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  body {
    foreach(b; 0 .. psmc.nrStates) {
      auto sum = 0.0;
      foreach(a; 0 .. psmc.nrStates) {
        sum += to.vec[a] * fullE(to_segsite, a) * transitionMatrix[a][b];
      }
      from.vec[b] = sum;
    }
  }
  
  void propagateMultiForward(in State_t from, State_t to,
        in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  body {
    auto dist = to_segsite.pos - from_segsite.pos;
    foreach(a; 0 .. psmc.nrStates) {
      if(to_segsite.obs[0] == 0) {
        auto prop = forwardPropagatorsMissing[dist - 1];
        auto sum = 0.0;
        foreach(b; 0 .. psmc.nrStates) {
          sum += from.vec[b] * prop[a][b];
        }
        to.vec[a] = sum;
      }
      else {
        auto prop = forwardPropagators[dist - 1];
        auto sum = 0.0;
        foreach(b; 0 .. psmc.nrStates) {
          sum += from.vec[b] * prop[a][b];
        }
        to.vec[a] = sum;
      }
    }
  }
  
  void propagateMultiBackward(in State_t to, State_t from,
        in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  body {  
    auto dist = to_segsite.pos - from_segsite.pos;
    foreach(b; 0 .. psmc.nrStates) {
      if(to_segsite.obs[0] == 0) {
        auto prop = backwardPropagatorsMissing[dist - 1];
        auto sum = 0.0;
        foreach(a; 0 .. psmc.nrStates) {
          sum += to.vec[a] * prop[a][b];
        }
        from.vec[b] = sum;
      }
      else {
        auto prop = backwardPropagators[dist - 1];
        auto sum = 0.0;
        foreach(a; 0 .. psmc.nrStates) {
          sum += to.vec[a] * prop[a][b];
        }
        from.vec[b] = sum;
      }
    }

  }
  
  const(PSMCmodel) getPSMC() const {
    return psmc;
  }
  
  @property size_t forwardStateSize() const {
    return psmc.nrStates;
  }

  @property size_t backwardStateSize() const {
    return psmc.nrStates;
  }
  
  State_t newForwardState() const {
    return new State_t(psmc.nrStates);
  }

  State_t newBackwardState() const {
    return new State_t(psmc.nrStates);
  }

  State_t newForwardState(StateVecAllocator stateAllocator) const {
    return new State_t(psmc.nrStates, stateAllocator);
  }

  State_t newBackwardState(StateVecAllocator stateAllocator) const {
    return new State_t(psmc.nrStates, stateAllocator);
  }
  
  void initialState(State_t s) const {
    foreach(a; 0 .. psmc.nrStates) {
      auto val = psmc.equilibriumProb(a);
      s.vec[a] = val;
    }
  }
  
  void setState(State_t s, double x, in SegSite_t segsite) const {
    foreach(a; 0 .. psmc.nrStates)
      s.vec[a] = x;
  }
  
  void getTransitionExpectation(State_t f, State_t b, in SegSite_t to_segsite, double[][] eMat) const
  {
    foreach(a; 0 .. psmc.nrStates)
      foreach(b_; 0 .. psmc.nrStates)
        eMat[a][b_] = f.vec[b_] * transitionMatrix[a][b_] * b.vec[a] * fullE(to_segsite, a);
  }

  void getEmissionExpectation(State_t f, State_t b, in SegSite_t to_segsite, double[][2] eMat) const
  {
    eMat[0][] = 0.0;
    eMat[1][] = 0.0;
    foreach(o; to_segsite.obs) {
      if(o > 0) {
        auto norm = 0.0;
        foreach(a; 0 .. psmc.nrStates) {
          auto n = f.vec[a] * b.vec[a] / cast(double)(to_segsite.obs.length);
          eMat[o - 1][a] = n;
          norm += n;
        }
        foreach(a; 0 .. psmc.nrStates)
          eMat[o - 1][a] /= norm;
      }
    }
  }
  
  @property size_t maxDistance() const {
    return cast(size_t)forwardPropagators.length;
  }
  
}

unittest {
  import std.math;
  writeln("testing propagationCoreFast and propagationCoreNaive propagateForward ");
  auto psmc = new PSMCmodel(0.01, 0.001, 10);
  auto lvl = 1.0e-8;
  auto dist = 10;
  auto propagationCore = new PropagationCore(psmc, dist);
  
  auto f = propagationCore.newForwardState();
  auto fNext = propagationCore.newForwardState();
  auto dummy_site = new SegSite_t(1, 1);
  propagationCore.setState(f, 1.0, dummy_site);
  foreach(i; 0 .. dist) {
    auto left_site = new SegSite_t(1 + i, 1);
    auto right_site = new SegSite_t(1 + i + 1, 1);
    propagationCore.propagateSingleForward(f, fNext, left_site, right_site);
    fNext.copy_into(f);
  }
  auto fSingles = propagationCore.newForwardState();
  f.copy_into(fSingles);
  propagationCore.setState(f, 1.0, dummy_site);
  auto left_site = new SegSite_t(1, 1);
  auto right_site = new SegSite_t(1 + dist, 1);
  propagationCore.propagateMultiForward(f, fNext, left_site, right_site);
  auto fSinglesA = fSingles.vec;
  auto fNextA = fNext.vec;
  foreach(aij; 0 .. psmc.nrStates) {
    assert(
        approxEqual(fNextA[aij], fSinglesA[aij], lvl, 0.0),
        text(fNext, " ", fSingles)
    );
  }
}
  
unittest {
  import std.math;
  writeln("testing propagationCoreFast and propagationCoreNaive propagateBackward ");
  import model.propagation_core;
  auto psmc = new PSMCmodel(0.01, 0.001, 10);
  auto lvl = 1.0e-8;
  auto dist = 10;
  auto propagationCore = new PropagationCore(psmc, dist);
  
  auto b = propagationCore.newBackwardState();
  auto bNext = propagationCore.newBackwardState();
  auto dummy_site = new SegSite_t(dist + 1, 1);
  propagationCore.setState(b, 1.0, dummy_site);
  foreach(i; 0 .. dist) {
    auto left_site = new SegSite_t(dist - i, 1);
    auto right_site = new SegSite_t(dist + 1 - i, 1);
    propagationCore.propagateSingleBackward(b, bNext, right_site, left_site);
    bNext.copy_into(b);
  }
  auto bSingles = propagationCore.newBackwardState();
  b.copy_into(bSingles);
  propagationCore.setState(b, 1.0, dummy_site);
  auto left_site = new SegSite_t(1, 1);
  auto right_site = new SegSite_t(dist + 1, 1);
  propagationCore.propagateMultiBackward(b, bNext, right_site, left_site);
  auto bSinglesA = bSingles.vec;
  auto bNextA = bNext.vec;
  foreach(aij; 0 .. psmc.nrStates) {
    assert(
        approxEqual(bNextA[aij], bSinglesA[aij], lvl, 0.0),
        text(bNext, ", ", bSingles)
    );
  }
}
