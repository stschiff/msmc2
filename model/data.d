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

module model.data; 
import std.stdio;
import std.string;
import std.conv;
import std.algorithm;
import std.math;
import core.stdc.stdlib;
import std.regex : match, regex, ctRegex;
import std.exception;
import std.range;
import model.time_intervals;

class SegSite_t {
  size_t pos; // rightMost position in the given segment
  size_t[] obs; // is 0 for missing, 1 for hom, and 2 for het, can have multiple obs for ambiguous phasing
  
  this(size_t pos, in size_t[] obs) {
    this.pos = pos;
    this.obs = obs.dup;
  }
  
  this(size_t pos, size_t obs) {
    this.pos = pos;
    this.obs = [obs];
  }
  
  @property SegSite_t dup() const {
    return new SegSite_t(pos, obs.dup);
  }
  
  override string toString() const {
    return text("Segsite(", pos, ", ", obs, ")");
  }
}

void checkDataLine(const char[] line) {
  auto r = regex(r"^\w+\s\d+\s\d+\s[ACTG01\?,]+$");
  enforce(match(line, r));
}

unittest {
  assertThrown(checkDataLine("1 20 5 AACC,AACA 2.44"));
  assertNotThrown(checkDataLine("1 20 5 AACC"));
  assertThrown(checkDataLine("4 5 2"));
  assertNotThrown(checkDataLine("1 10 5 ACC"));
  assertThrown(checkDataLine("1 20 5 AGGSSXX"));
}

size_t getNrHaplotypesFromFile(string filename) {
  auto _file = File(filename, "r");
  scope(exit) _file.close();
  auto line = _file.readln();
  line = line.strip();
  checkDataLine(line);
  auto fields = line.strip().split();
  if(fields.length < 4)
    return 2;
  else {
    auto splitted = fields[3].split(",");
    return cast(size_t)splitted[0].length;
  }
}

unittest {
  auto tmp = File("/tmp/nrHaplotypesTest.txt", "w");
  tmp.writeln("1 10 5 ACC,CCA");
  tmp.close();
  assert(getNrHaplotypesFromFile("/tmp/nrHaplotypesTest.txt") == 3);
}

SegSite_t[][] readSegSites(string filename, size_t[2][] indices, bool skipAmbiguous) {
  // format: chr position nr_calledSites [alleles] -> tab separated
  // [alleles]: comma-separated for ambiguous phasing
  // if no alleles are given, assume M=2 and "01"
  // returns data for each pair of indices
  
  auto ret = new SegSite_t[][indices.length];

  auto f = File(filename, "r");
  long lastPos = -1;

  foreach(line; f.byLine()) {
    // checkDataLine(line.strip());
    auto fields = line.strip().split();
    auto pos = to!size_t(fields[1]);
    auto nrCalledSites = to!size_t(fields[2]);
    auto raw_allele_strings = split(fields[3], ",");
    if(lastPos == -1) {
      lastPos = pos - nrCalledSites;
    }
    
    enforce(nrCalledSites <= pos - lastPos);
    enforce(nrCalledSites > 0, "nr of called sites must be positive!");
    
    foreach(i, ind; indices) {
     if(has_missing_data(raw_allele_strings, ind)) {
        if(nrCalledSites < pos - lastPos) { // missing data
          ret[i] ~= new SegSite_t(pos - nrCalledSites, 0UL);
        }
        if(nrCalledSites > 1)
          ret[i] ~= new SegSite_t(pos - 1, 1UL);
        ret[i] ~= new SegSite_t(pos, 0UL);
        lastPos = pos;
      }
      else {
        size_t[] allele_indices;
        foreach(allele_string; raw_allele_strings) {
          auto selected_allele_string = [allele_string[ind[0]], allele_string[ind[1]]];
          allele_indices ~= selected_allele_string[0] == selected_allele_string[1] ? 1UL : 2UL;
        }
        if(nrCalledSites < pos - lastPos) { // missing data
          ret[i] ~= new SegSite_t(pos - nrCalledSites, 0UL);
        }
        if(allele_indices.length > 1)
            allele_indices = allele_indices.uniq().array();
        if(skipAmbiguous && allele_indices.length > 1)
          ret[i] ~= new SegSite_t(pos, 0UL);
        else
          ret[i] ~= new SegSite_t(pos, allele_indices);
      }
   }
    
   lastPos = pos;
  } 
  return ret;
}

unittest {
  writeln("test readSegSites");
  auto tmp_file = File("/tmp/msmc_data_unittest.tmp", "w");
  tmp_file.writeln("1 1000000 42 AACC");
  tmp_file.writeln("1 1000004 2 ACCG");
  tmp_file.writeln("1 1000008 3 ACC?,ATTA");
  tmp_file.writeln("1 1000012 4 ACCG,TTGA");
  tmp_file.close();

  auto segsites = readSegSites("/tmp/msmc_data_unittest.tmp", [[0UL, 1]], false)[0];
  assert(segsites[0].pos == 1000000 && segsites[0].obs == [1]);
  assert(segsites[1].pos == 1000002 && segsites[1].obs == [0]);
  assert(segsites[3].pos == 1000005 && segsites[3].obs == [0]);
  assert(segsites[4].pos == 1000008 && segsites[4].obs == [2]);
  assert(segsites[5].pos == 1000012 && segsites[5].obs == [2, 1]);
  
}

unittest {
  writeln("test readSegSites for M=2");
  auto tmp_file = File("/tmp/msmc_data_unittest.tmp", "w");
  tmp_file.writeln("1 1000000 42 AC");
  tmp_file.writeln("1 1000004 2 CC");
  tmp_file.writeln("1 1000008 3 C?,AT");
  tmp_file.writeln("1 1000012 4 AA,TT");
  tmp_file.close();

  auto segsites = readSegSites("/tmp/msmc_data_unittest.tmp", [[0UL, 1]], false)[0];
  assert(segsites[0].pos == 1000000 && segsites[0].obs == [2]);
  assert(segsites[1].pos == 1000002 && segsites[1].obs == [0]);
  assert(segsites[3].pos == 1000005 && segsites[3].obs == [0]);
  assert(segsites[4].pos == 1000007 && segsites[4].obs == [1]);
  assert(segsites[5].pos == 1000008 && segsites[5].obs == [0]);
  assert(segsites[6].pos == 1000012 && segsites[6].obs == [1]);
  
}

bool has_missing_data(in char[][] raw_allele_strings, size_t[2] indices) {
    auto is_missing = false;
    foreach(raw_allele_string; raw_allele_strings) {
      foreach(i; indices) {
       if(!canFind("ACTG01", raw_allele_string[i])) {
          return true;
        }
      }
    }
    return false;
}
 
double getTheta(in SegSite_t[][] data) {
  size_t nr_hets; 
  size_t called_sites;
  foreach(d; data) {
    size_t lastPos = 0;
    foreach(dd; d) {
      if(dd.obs[0] > 0) {
        if(lastPos > 0)
          called_sites += dd.pos - lastPos;
        if(dd.obs.any!(o => o > 1))
          nr_hets += 1;
      }
      lastPos = dd.pos;
    }
  }
  
  return cast(double)nr_hets / cast(double)called_sites;
}
