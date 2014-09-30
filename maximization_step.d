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
 
import std.math;
import std.stdio;
import std.random;
import std.exception;
import std.algorithm;
import model.psmc_model;
import powell;
import logger;

double delta_x = 1.0e-4;

PSMCmodel getMaximization(double[][] transitions, double[][2] emissions, PSMCmodel params,
                          in size_t[] timeSegmentPattern, bool fixedRecombination)
{
  auto minFunc = new MinFunc(transitions, emissions, params, timeSegmentPattern, fixedRecombination);

  auto powell = new Powell!MinFunc(minFunc);
  auto x = minFunc.initialValues();
  auto startVal = minFunc(x);
  auto xNew = powell.minimize(x);
  auto endVal = minFunc(xNew);
  logInfo(format(", Q-function before: %s, after:%s\n", startVal, endVal));
  return minFunc.makeParamsFromVec(xNew);
}

double[2][] getLambdaCI(double[][] transitions, double[][2] emissions, PSMCmodel params, in size_t[] timeSegmentPattern) {

  auto minFunc = new MinFunc(transitions, emissions, params, timeSegmentPattern, false);

  auto basicLL = minFunc.logLikelihood(params);
  
  auto lambdaVec = params.lambdaVec.dup;
  auto ret = new double[2][lambdaVec.length];
  
  auto indexOffset = 0;
  foreach(nrI; timeSegmentPattern) {
    auto newLambdaVec = lambdaVec.dup;
    foreach(i; 0 .. nrI)
      newLambdaVec[indexOffset + i] *= exp(-delta_x);
    auto newParams = new PSMCmodel(params.mutationRate, params.recombinationRate, newLambdaVec, params.timeIntervals);
    auto newLLsmaller = minFunc.logLikelihood(newParams);

    newLambdaVec = lambdaVec.dup;
    foreach(i; 0 .. nrI)
      newLambdaVec[indexOffset + i] *= exp(delta_x);
    newParams = new PSMCmodel(params.mutationRate, params.recombinationRate, newLambdaVec, params.timeIntervals);
    auto newLLgreater = minFunc.logLikelihood(newParams);
    
    auto limitsCI = getCIlimits(basicLL, 0.001, newLLsmaller, newLLgreater);
    foreach(i; 0 .. nrI) {
      ret[indexOffset + i][0] = lambdaVec[indexOffset + i] * exp(limitsCI[0]);
      ret[indexOffset + i][1] = lambdaVec[indexOffset + i] * exp(limitsCI[1]);
    }
    
    indexOffset += nrI;
  }
  return ret;
}

double[2] getRecombinationCI(double[][] transitions, double[][2] emissions, PSMCmodel params, in size_t[] timeSegmentPattern) {
  auto minFunc = new MinFunc(transitions, emissions, params, timeSegmentPattern, false);

  auto basicLL = minFunc.logLikelihood(params);
  
  auto rec = params.recombinationRate;
  
  auto newRec = rec * exp(-delta_x);
  auto newParams = new PSMCmodel(params.mutationRate, newRec, params.lambdaVec, params.timeIntervals);
  auto newLLsmaller = minFunc.logLikelihood(newParams);

  newRec = rec * exp(delta_x);
  newParams = new PSMCmodel(params.mutationRate, newRec, params.lambdaVec, params.timeIntervals);
  auto newLLgreater = minFunc.logLikelihood(newParams);
  
  auto CI = getCIlimits(basicLL, delta_x, newLLsmaller, newLLgreater);
  
  double[2] ret;
  ret[0] = rec * exp(CI[0]);
  ret[1] = rec * exp(CI[1]);
  
  return ret;
}

double[2] getCIlimits(double basicLL, double delta_x, double yl, double yh) {
  auto xl = -delta_x;
  auto xh = delta_x;
  
  auto denom = (xl*xl * xh - xl * xh*xh);
  auto a = (yl * xh - xl * yh) / denom;
  auto b = (xl*xl * yh - yl * xh*xh) / denom;
  
  auto mu = -b / (2.0 * a);
  auto sigma = sqrt(-1.0 / (2.0 * a));
  
  double[2] ret;
  ret[0] = mu - 2.0 * sigma;
  ret[1] = mu + 2.0 * sigma;
  return ret;
}

class MinFunc {
  
  PSMCmodel initialParams;
  const size_t[] timeSegmentPattern;
  size_t nrParams;
  const double[][] transitions;
  const double[][2] emissions;
  bool fixedRecombination;
  
  this(in double[][] transitions, in double[][2] emissions, PSMCmodel initialParams,
       in size_t[] timeSegmentPattern, bool fixedRecombination)
  {
    this.initialParams = initialParams;
    this.timeSegmentPattern = timeSegmentPattern;
    this.transitions = transitions;
    this.emissions = emissions;
    this.fixedRecombination = fixedRecombination;
    nrParams = cast(size_t)timeSegmentPattern.length;
    if(!fixedRecombination)
      nrParams += 1;
  }
  
  double opCall(in double[] x) {
    PSMCmodel newParams = makeParamsFromVec(x);
    return -logLikelihood(newParams);
  };
  
  double[] initialValues()
  out(x) {
    assert(x.length == nrParams);
  }
  body {
    auto x = getXfromLambdaVec(initialParams.lambdaVec);
    if(!fixedRecombination)
      x ~= log(initialParams.recombinationRate);
    return x;
  }
  
  double[] getXfromLambdaVec(in double[] lambdaVec)
  out(x) {
    assert(x.length == timeSegmentPattern.length);
  }
  body {
    double[] ret;
    size_t lIndex = 0;
    foreach(nrIntervalsInSegment; timeSegmentPattern) {
      ret ~= log(lambdaVec[lIndex]);
      lIndex += nrIntervalsInSegment;
    }
    return ret;
  }
  
  PSMCmodel makeParamsFromVec(in double[] x) {
    auto lambdaVec = getLambdaVecFromX(x);
    auto recombinationRate = fixedRecombination ? initialParams.recombinationRate : getRecombinationRateFromX(x);
    return new PSMCmodel(initialParams.mutationRate, recombinationRate, lambdaVec, initialParams.timeIntervals);
  }
    
  double[] getLambdaVecFromX(in double[] x)
  in {
    assert(x.length == nrParams);
  }
  body {
    auto lambdaVec = initialParams.lambdaVec.dup;
    auto timeIndex = 0U;
    foreach(segmentIndex, nrIntervalsInSegment; timeSegmentPattern) {
      foreach(intervalIndex; 0 .. nrIntervalsInSegment) {
        auto xIndex = segmentIndex;
        lambdaVec[timeIndex] = exp(x[xIndex]);
        timeIndex += 1;
      }
    }
    return lambdaVec;
  }
  

  double getRecombinationRateFromX(in double[] x)
  in {
    assert(!fixedRecombination);
  }
  body {
    return exp(x[$ - 1]);
  }

  double logLikelihood(PSMCmodel params) {
    double ret = 0.0;
    foreach(a; 0 .. initialParams.nrStates) {
      foreach(b; 0 .. initialParams.nrStates) {
        ret += transitions[a][b] * log(params.transitionProb(a, b));
      }
      ret += emissions[0][a] * log(params.emissionProb(1, a));
      ret += emissions[1][a] * log(params.emissionProb(2, a));
    }
    return ret;
  }

}

unittest {
  writeln("test minfunc.getLambdaFromX");
  import std.conv;
  
  auto lambdaVec = [1.0, 1, 4, 4];
  auto params = new PSMCmodel(0.01, 0.001, lambdaVec, 4);
  auto transitions = new double[][](params.nrStates, params.nrStates);
  double[][2] emissions = [new double[params.nrStates], new double[params.nrStates]];
  auto timeSegmentPattern = [2UL, 2];
  
  auto minFunc = new MinFunc(transitions, emissions, params, timeSegmentPattern, false);
  auto rho = 0.001;
  auto x = minFunc.getXfromLambdaVec(lambdaVec);
  x ~= log(rho);
  auto lambdaFromX = minFunc.getLambdaVecFromX(x);
  auto rhoFromX = minFunc.getRecombinationRateFromX(x);
  foreach(i; 0 .. lambdaVec.length)
    assert(approxEqual(lambdaFromX[i], lambdaVec[i], 1.0e-8, 0.0), text(lambdaFromX[i], " ", lambdaVec[i]));
  assert(approxEqual(rhoFromX, rho, 1.0e-8, 0.0));

  minFunc = new MinFunc(transitions, emissions, params, timeSegmentPattern, true);
  x = minFunc.getXfromLambdaVec(lambdaVec);
  lambdaFromX = minFunc.getLambdaVecFromX(x);
  foreach(i; 0 .. lambdaVec.length)
    assert(approxEqual(lambdaFromX[i], lambdaVec[i], 1.0e-8, 0.0), text(lambdaFromX[i], " ", lambdaVec[i]));
}
