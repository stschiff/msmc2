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

class State_t {
  double[] vec;
  
  this(size_t nrS) {
    vec = new double[nrS];
  }
  
  this(size_t nrS, StateVecAllocator stateAllocator) {
    vec = stateAllocator.allocate(nrS);
  }
  
  @property double norm() const {
     return reduce!"a+b"(vec);
  }
  
  void scale(double x) {
    vec[] *= x;
  }
  
  void copy_into(State_t dest) {
    dest.vec[] = vec[];
  }
  
  override string toString() const {
    return text(vec);
  }

  void setZero() {
    vec[] = 0.0;
  }

}
