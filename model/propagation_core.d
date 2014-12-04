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
import model.gsl_matrix_vector;

class PropagationCore {
  
  gsl_matrix*[] forwardPropagators, backwardPropagators;
  gsl_matrix*[] forwardPropagatorsMissing, backwardPropagatorsMissing;
  
  double[][] emissionProbs; // first index: second index: obs
  gsl_matrix* transitionMatrix;
  
  const PSMCmodel psmc;
  
  this(in PSMCmodel psmc, size_t maxDistance) {
    enforce(maxDistance > 0);
    this.psmc = psmc;

    forwardPropagators = new gsl_matrix*[maxDistance];
    backwardPropagators = new gsl_matrix*[maxDistance];
    forwardPropagatorsMissing = new gsl_matrix*[maxDistance];
    backwardPropagatorsMissing = new gsl_matrix*[maxDistance];
    emissionProbs = new double[][](3, psmc.nrStates);
    transitionMatrix = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
      
    foreach(i; 0 .. 3)
      foreach(a; 0 .. psmc.nrStates)
        emissionProbs[i][a] = psmc.emissionProb(i, a);
    
    foreach(a; 0 .. psmc.nrStates)
      foreach(b; 0 .. psmc.nrStates)
        gsl_matrix_set(transitionMatrix, a, b, psmc.transitionProb(a, b));
      
    computeForwardPropagators(forwardPropagators, false, maxDistance);
    computeBackwardPropagators(backwardPropagators, false, maxDistance);

    computeForwardPropagators(forwardPropagatorsMissing, true, maxDistance);
    computeBackwardPropagators(backwardPropagatorsMissing, true, maxDistance);
    
  }
  ~this() {
    gsl_matrix_free(transitionMatrix);
    foreach(i; 0 .. maxDistance) {
        gsl_matrix_free(forwardPropagators[i]);
        gsl_matrix_free(forwardPropagatorsMissing[i]);
        gsl_matrix_free(backwardPropagators[i]);
        gsl_matrix_free(backwardPropagatorsMissing[i]);
    }
  }
  
  private void computeForwardPropagators(gsl_matrix*[] ret, bool missing_data, size_t maxDistance) const
  {
    ret[0] = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
    foreach(a; 0 .. psmc.nrStates) {
      double e = missing_data ? 1.0 : emissionProbs[1][a];
      foreach(b; 0 .. psmc.nrStates) {
        auto val = gsl_matrix_get(transitionMatrix, a, b) * e;
        gsl_matrix_set(ret[0], a, b, val);
      }
    }

    foreach(distance; 1 .. maxDistance) {
      ret[distance] = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
      gsl_blas_dgemm_checked(CBLAS_TRANSPOSE_t.CblasNoTrans, CBLAS_TRANSPOSE_t.CblasNoTrans,
                     1.0, ret[0], ret[distance - 1], 0.0, ret[distance]);
    }
  }
  
  private void computeBackwardPropagators(gsl_matrix*[] ret, bool missing_data, size_t maxDistance) const
  {
    ret[0] = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
    foreach(a; 0 .. psmc.nrStates) {
      double e = missing_data ? 1.0 : emissionProbs[1][a];
      foreach(b; 0 .. psmc.nrStates) {
        auto val = gsl_matrix_get(transitionMatrix, a, b) * e;
        gsl_matrix_set(ret[0], a, b, val);
      }
    }

    foreach(distance; 1 .. maxDistance) {
      ret[distance] = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
      gsl_blas_dgemm(CBLAS_TRANSPOSE_t.CblasNoTrans, CBLAS_TRANSPOSE_t.CblasNoTrans,
                     1.0, ret[distance - 1], ret[0], 0.0, ret[distance]);
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
    gsl_blas_dgemv(CBLAS_TRANSPOSE_t.CblasNoTrans, 1.0, transitionMatrix, from.vec, 0.0, to.vec);
    foreach(a; 0 .. psmc.nrStates) {
        auto val = gsl_vector_get(to.vec, a) * fullE(to_segsite, a);
        gsl_vector_set(to.vec, a, val);
    }
  }
  
  void propagateSingleBackward(in State_t to, State_t from,
            in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  body {
    gsl_blas_dgemv(CBLAS_TRANSPOSE_t.CblasTrans, 1.0, transitionMatrix, to.vecE, 0.0, from.vec);
    foreach(a; 0 .. psmc.nrStates) {
        auto val = gsl_vector_get(from.vec, a) * fullE(from_segsite, a);
        gsl_vector_set(from.vecE, a, val);
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
      if(to_segsite.obs[0] == 0) {
        auto prop = forwardPropagatorsMissing[dist - 1];
        gsl_blas_dgemv(CBLAS_TRANSPOSE_t.CblasNoTrans, 1.0, prop, from.vec, 0.0, to.vec);
      }
      else {
        auto prop = forwardPropagators[dist - 1];
        gsl_blas_dgemv(CBLAS_TRANSPOSE_t.CblasNoTrans, 1.0, prop, from.vec, 0.0, to.vec);
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
      if(to_segsite.obs[0] == 0) {
        auto prop = backwardPropagatorsMissing[dist - 1];
        gsl_blas_dgemv(CBLAS_TRANSPOSE_t.CblasTrans, 1.0, prop, to.vec, 0.0, from.vec);
      }
      else {
        auto prop = backwardPropagators[dist - 1];
        gsl_blas_dgemv(CBLAS_TRANSPOSE_t.CblasTrans, 1.0, prop, to.vec, 0.0, from.vec);
      }
      foreach(a; 0 .. psmc.nrStates) {
        auto val = gsl_vector_get(from.vec, a) * fullE(from_segsite, a);
        gsl_vector_set(from.vecE, a, val);
      }
  }
  
  const(PSMCmodel) getPSMC() const {
    return psmc;
  }
  
  @property size_t forwardStateSize() const {
    return psmc.nrStates;
  }

  @property size_t backwardStateSize() const {
    return 2 * psmc.nrStates;
  }
  
  State_t newForwardState() const {
    return new State_t(psmc.nrStates, false);
  }

  State_t newBackwardState() const {
    return new State_t(psmc.nrStates, true);
  }

  State_t newForwardState(StateVecAllocator stateAllocator) const {
    return new State_t(psmc.nrStates, false, stateAllocator);
  }

  State_t newBackwardState(StateVecAllocator stateAllocator) const {
    return new State_t(psmc.nrStates, true, stateAllocator);
  }
  
  void initialState(State_t s) const {
    foreach(a; 0 .. psmc.nrStates) {
      auto val = psmc.equilibriumProb(a);
      gsl_vector_set(s.vec, a, val);
    }
  }
  
  void setState(State_t s, double x, in SegSite_t segsite) const 
  body {
    foreach(a; 0 .. psmc.nrStates) {
      gsl_vector_set(s.vec, a, x);
      if(s.isBwd) {
          auto val = x * fullE(segsite, a);
          gsl_vector_set(s.vecE, a, val);
      }
    }
  }
  
  void getTransitionExpectation(State_t f, State_t b, in SegSite_t to_segsite, double[][] eMat) const
  {
    foreach(a; 0 .. psmc.nrStates)
      foreach(b_; 0 .. psmc.nrStates)
        eMat[a][b_] = gsl_vector_get(f.vec, b_) * gsl_matrix_get(transitionMatrix, a, b_) *
                      gsl_vector_get(b.vec, a) * fullE(to_segsite, a);
  }

  void getEmissionExpectation(State_t f, State_t b, in SegSite_t to_segsite, double[][2] eMat) const
  {
    eMat[0][] = 0.0;
    eMat[1][] = 0.0;
    foreach(o; to_segsite.obs) {
      if(o > 0) {
        auto norm = 0.0;
        foreach(a; 0 .. psmc.nrStates) {
          auto n = gsl_vector_get(f.vec, a) * gsl_vector_get(b.vec, a) / cast(double)(to_segsite.obs.length);
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
  writeln("testing propagateForward ");
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
  foreach(a; 0 .. psmc.nrStates) {
    assert(
        approxEqual(gsl_vector_get(fNextA, a), gsl_vector_get(fSinglesA, a), lvl, 0.0),
        text(fNext, " ", fSingles)
    );
  }
}
  
unittest {
  import std.math;
  writeln("testing propagateBackward ");
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
  foreach(a; 0 .. psmc.nrStates) {
    assert(
        approxEqual(gsl_vector_get(bNextA, a), gsl_vector_get(bSinglesA, a), lvl, 0.0),
        text(bNext, ", ", bSingles)
    );
  }
}
