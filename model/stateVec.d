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

module model.stateVec;
import std.conv;
import std.algorithm;
import model.stateVecAllocator;
import model.gsl_matrix_vector;
import std.range;

class State_t {
  gsl_vector_view view;
  gsl_vector_view viewE;
  gsl_vector* vec;
  gsl_vector* vecE;
  bool isBwd;
  
  this(size_t nrS, bool isBwd) { 
    this.isBwd = isBwd;
    vec = gsl_vector_alloc(nrS);
    view = gsl_vector_subvector(vec, 0, nrS);
    if(isBwd) {
        vecE = gsl_vector_alloc(nrS);
        viewE = gsl_vector_subvector(vecE, 0, nrS);
    }
  }
  
  this(size_t nrS, bool isBwd, StateVecAllocator stateAllocator) {
    this.isBwd = isBwd;
    view = stateAllocator.allocate(nrS);
    vec = &view.vector;
    if(isBwd) {
        viewE = stateAllocator.allocate(nrS);
        vecE = &viewE.vector;
    }
  }

  //~this() {
  //  gsl_vector_free(vec);
  //  if(isBwd)
  //      gsl_vector_free(vecE);
  //}
  
  @property double norm() const {
     return iota(vec.size).map!(i => gsl_vector_get(vec, i)).reduce!"a+b"();
  }
  
  void scale(double x) {
    gsl_vector_scale(vec, x);
    if(isBwd)
        gsl_vector_scale(vecE, x);
  }
  
  void copy_into(State_t dest) {
    gsl_vector_memcpy(dest.vec, vec);
    if(isBwd)
        gsl_vector_memcpy(dest.vecE, vecE);
  }
  
  override string toString() const {
    auto body_ = iota(vec.size).map!(i => text(gsl_vector_get(vec, i))).joiner(",").array();
    return to!string('[' ~ body_ ~ ']');
  }

  void setZero() {
    gsl_vector_set_zero(vec);
    if(isBwd)
        gsl_vector_set_zero(vecE);
  }

}
