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

module model.transition_rate;
import std.math;
import std.conv;
import std.stdio;
import std.exception;
import std.parallelism;
import std.range;
import model.time_intervals;

class TransitionRate {
  double rho;
  const TimeIntervals timeIntervals;
  const double[] lambdaVec;
  size_t nrTimeIntervals;

  private double[][] transitionProbabilities;
  
  this(in TimeIntervals timeIntervals, double rho, in double[] lambdaVec) {
    enforce(rho > 0.0, "need positive recombination rate");
    this.timeIntervals = timeIntervals;
    this.nrTimeIntervals = timeIntervals.nrIntervals;
    this.rho = rho;
    this.lambdaVec = lambdaVec;
    fillTransitionProbabilitiesParallel();
    // fillTransitionProbabilitiesSingleThread();
  }
  
  private double integrateLambda(double from, double to, size_t fromIndex, size_t toIndex, double lambdaFac=1.0) const {
    double sum = 0.0;
    if(fromIndex == toIndex) {
      return exp(-(to - from) * lambdaFac * lambdaVec[toIndex]);
    }
    foreach(kappa; fromIndex + 1 .. toIndex) {
      sum += lambdaFac * lambdaVec[kappa] * timeIntervals.delta(kappa);
    }
    
    double ret = exp(-(timeIntervals.rightBoundary(fromIndex) - from) *
               lambdaFac * lambdaVec[fromIndex] - sum - (to - timeIntervals.leftBoundary(toIndex)) * lambdaFac * lambdaVec[toIndex]);
    return ret;
  }
  
  private void fillTransitionProbabilitiesSingleThread() {
    transitionProbabilities = new double[][](nrTimeIntervals, nrTimeIntervals);
    foreach(b; 0 .. nrTimeIntervals) {
      auto sum = 0.0;
      foreach(a; 0 .. nrTimeIntervals) {
        if(a != b) {
          transitionProbabilities[a][b] = transitionProbabilityOffDiagonal(a, b);
          sum += transitionProbabilities[a][b];
        }
      }
      transitionProbabilities[b][b] = 1.0 - sum;
    }
  }

  private void fillTransitionProbabilitiesParallel() {
    transitionProbabilities = new double[][](nrTimeIntervals, nrTimeIntervals);
    foreach(b; taskPool.parallel(iota(nrTimeIntervals))) {
      auto sum = 0.0;
      foreach(a; 0 .. nrTimeIntervals) {
        if(a != b) {
          transitionProbabilities[a][b] = transitionProbabilityOffDiagonal(a, b);
          sum += transitionProbabilities[a][b];
        }
      }
      transitionProbabilities[b][b] = 1.0 - sum;
    }
  }
  
  private double transitionProbabilityOffDiagonal(size_t a, size_t b) const {
    if(a < b) {
      return q2IntegralSmaller(a, b);
    }
    if(a > b) {
      return q2IntegralGreater(a, b);
    }
    assert(false);
  }
  
  private double q2IntegralSmaller(size_t a, size_t b) const
    in {
      assert(a < b);
    }
  body {
    auto meanTime = timeIntervals.meanTimeWithLambda(b, lambdaVec[b]);
    double integ = (1.0 - exp(-(timeIntervals.delta(a) * 2.0 * lambdaVec[a]))) / (2.0 * lambdaVec[a]);
    double sum = 0.0;
    foreach(g; 0 .. a) {
      sum += 2.0 * (1.0 - exp(-timeIntervals.delta(g) * 2.0 * lambdaVec[g])) *
        integrateLambda(timeIntervals.rightBoundary(g), timeIntervals.leftBoundary(a), g + 1, a, 2.0) / (2.0 * lambdaVec[g]) * integ;
    }

    sum += 2.0 * (timeIntervals.delta(a) - integ) / (2.0 * lambdaVec[a]);
    
    double ret = (1.0 - exp(-rho * 2 * meanTime)) / (meanTime * 2) * lambdaVec[a] * sum;
    return ret;
  }

  private double q2IntegralGreater(size_t a, size_t b) const
    in {
      assert(a > b);
    }
  body {
    auto meanTime = timeIntervals.meanTimeWithLambda(b, lambdaVec[b]);
    double integ = integrateLambda(meanTime, timeIntervals.leftBoundary(a), b, a) / 
                   lambdaVec[a] * (1.0 - exp(-timeIntervals.delta(a) * lambdaVec[a]));
    double sum = 0.0;
    foreach(g; 0 .. b) {
      sum += 2.0 * (1.0 - exp(-2.0 * lambdaVec[g] * timeIntervals.delta(g))) / (2.0 * lambdaVec[g]) * 
             integrateLambda(timeIntervals.rightBoundary(g), meanTime, g + 1, b, 2.0);
    }
    sum += 2.0 * (1.0 - exp(-2.0 * lambdaVec[b] * (meanTime - timeIntervals.leftBoundary(b)))) / (2.0 * lambdaVec[b]);
      
    return integ * (1.0 - exp(-rho * 2.0 * meanTime)) / (meanTime * 2.0) * lambdaVec[a] * sum;
  }
  
  double transitionProbability(size_t a, size_t b) const {
    return transitionProbabilities[a][b];
  }
  
  double equilibriumProbability(size_t a) const {
    return integrateLambda(0.0, timeIntervals.leftBoundary(a), 0, a) *
        (1.0 - exp(-timeIntervals.delta(a) * lambdaVec[a]));
  }
}

unittest {
  writeln("test transitionRate.equilibriumProbability");
  auto T = 4UL;
  auto r = 0.01;
  auto lambdaVec = [1.0, 0.1, 1, 2];
  auto lvl = 1.0e-8;
  
  auto timeIntervals = TimeIntervals.getQuantileBoundaries(T, 1.0);
  auto transitionRate = new TransitionRate(new TimeIntervals(timeIntervals), r, lambdaVec);
  
  auto sum = 0.0;
  foreach(a; 0 .. T)
    sum += transitionRate.equilibriumProbability(a);
  assert(approxEqual(sum, 1.0, lvl, 0.0), text(sum, " should be 1.0"));
}

unittest {  
  writeln("test transitionRate.transitionProbability");
  auto T = 4UL;
  auto r = 0.01;
  auto lambdaVec = [1, 0.1, 1, 2];
  auto lvl = 1.0e-8;
  
  auto timeIntervals = TimeIntervals.getQuantileBoundaries(T, 1.0);
  auto transitionRate = new TransitionRate(new TimeIntervals(timeIntervals), r, lambdaVec);
  
  foreach(b; 0 .. T) {
    auto sum = 0.0;
    foreach(a; 0 .. T)
      sum += transitionRate.transitionProbability(a, b);
    assert(approxEqual(sum, 1.0, lvl, 0.0), text("transition prob from state ", b, ": ", sum, " should be 1.0"));
  }
}

unittest {
  writeln("test transitionRate.transitionProbability values");
  auto msmc_vals = [[0.999886,1.63405e-05,1.63378e-05,1.63345e-05,1.63303e-05,1.63246e-05,1.63157e-05,1.62937e-05],
                    [3.03282e-05,0.999727,9.29926e-05,9.29737e-05,9.29499e-05,9.29172e-05,9.28668e-05,9.27412e-05],
                    [3.53829e-05,0.00010808,0.999633,0.000189957,0.000189909,0.000189842,0.000189739,0.000189482],
                    [1.74971e-05,5.34459e-05,9.27996e-05,0.99939,0.0001761,0.000176038,0.000175942,0.000175704],
                    [7.77647e-06,2.37537e-05,4.12443e-05,7.76406e-05,0.999035,0.000163879,0.00016379,0.000163568],
                    [9.27006e-07,2.8316e-06,4.91658e-06,9.25526e-06,1.94749e-05,0.998333,4.24982e-05,4.24408e-05],
                    [2.89998e-06,8.85817e-06,1.53807e-05,2.89535e-05,6.09238e-05,0.000133011,0.997586,0.000257459],
                    [1.95024e-05,5.95714e-05,0.000103436,0.000194713,0.000409714,0.000894505,0.00173237,0.999062]];
  auto mu = 0.001;
  auto rho = 0.001;
  auto lambdaVec = [1.0, 2.0, 3.0, 2.0, 1.0, 0.1, 0.2, 1.0];
  auto T = 8UL;
  auto timeIntervals = TimeIntervals.getQuantileBoundaries(T, 1.0);
  auto transitionRate = new TransitionRate(new TimeIntervals(timeIntervals), rho, lambdaVec);

  foreach(a; 0 .. T) {
    foreach(b; 0 .. T) {
      writef("%s,", transitionRate.transitionProbability(a, b));
    }
    write("\n");
  }

}