#include "SpectraSTSearchTaskStats.hpp"
#include "SpectraSTLog.hpp"
#include <sstream>

/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.06 


Copyright (C) 2006 Henry Lam

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA

Henry Lam
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
hlam@systemsbiology.org

*/

/* Class: SpectraSTSearchTaskStats
 * 
 * Class to manage search task stats
 */

extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog* g_log;

// constructor
SpectraSTSearchTaskStats::SpectraSTSearchTaskStats(bool doDecoyAnalysis) :
  m_doDecoyAnalysis(doDecoyAnalysis),
  m_numScans(0),
  m_numNotSelected(0),
  m_numMissing(0),
  m_numMS1(0),
  m_numFailedFilter(0),
  m_numSearched(0),
  m_numLikelyGood(0),
  m_numTopHits(6, 0),
  m_numBadTopHits(6, 0),
  m_numGoodTopHits(6, 0),
  m_numDecoyTopHits(6, 0),
  m_numDecoyBadTopHits(6, 0),
  m_numDecoyGoodTopHits(6, 0),
  m_numLowerHits(11, 0),
  m_numDecoyLowerHits(11, 0),
  m_numTotalLowerHits(0),
  m_numTotalDecoyLowerHits(0) {
//  m_numTotalLowerSingletonHits(0), 
//  m_numTotalDecoyLowerSingletonHits(0) {
}

// destructor
SpectraSTSearchTaskStats::~SpectraSTSearchTaskStats() {

}

string SpectraSTSearchTaskStats::getScanStatsStr() {
  
  stringstream ss;
  ss << "(Max " << m_numScans << " scans; " << m_numSearched << " searched, ";
  ss << m_numLikelyGood << " likely good; ";
  if (m_numNotSelected > 0) {
    ss << m_numNotSelected << " not selected; ";
  }
  ss << m_numFailedFilter << " failed filter; " << m_numMissing << " missing; " << m_numMS1 << " MS1)";

  return (ss.str());
}
    
// processSearch - take the search results and store some interesting stats. TO BE EXPANDED
void SpectraSTSearchTaskStats::processSearchResult(SpectraSTSearch* s) {
  
  if (s->isNoMatch()) {
    return;
  }
  
  if (s->isLikelyGood()) {
    m_numLikelyGood++;
  }    
  
  SpectraSTLibEntry* entry = s->m_candidates[0]->getEntry();
  // unsigned int numPeaks = entry->getPeakList()->getNumPeaks();
  // unsigned int numAssignedPeaks = entry->getPeakList()->getNumAssignedPeaks();
  // double fracAssigned = (double)(numAssignedPeaks) / (double)(numPeaks);
  
  int ch = entry->getCharge() - 1;
  //unsigned int ch = (unsigned int)(fracAssigned * 5 - 0.01) ;
  if (ch < 0) ch = 0;
  if (ch > 4) ch = 4; // treat all 5+ charges in same bin
     

  m_numTopHits[ch]++;
  if (s->isLikelyGood()) {
    m_numGoodTopHits[ch]++;
  } else if (s->isLikelyBad()) {
    m_numBadTopHits[ch]++;
  }
 
  if (!m_doDecoyAnalysis) return;
 
  // Decoy analysis
  
  if (s->isDecoy(1)) {
    m_numDecoyTopHits[ch]++;
    if (s->isLikelyGood()) {
      m_numDecoyGoodTopHits[ch]++;
    } else if (s->isLikelyBad()) {
      m_numDecoyBadTopHits[ch]++;
    }
    
    map<string, int>::iterator found = m_decoyCounts.find(entry->getName());
    if (found != m_decoyCounts.end()) {
        found->second++;
    } else {
        m_decoyCounts[entry->getName()] = 1;
    }
  }
  
  for (unsigned int rank = 2; rank <= 10 && rank <= (unsigned int)(s->m_candidates.size()); rank++) {
    bool isSingleton = s->isSingleton(rank);
    m_numTotalLowerHits++;
//    if (isSingleton) {
//      m_numTotalLowerSingletonHits++;
//    }
    m_numLowerHits[rank]++;
    if (s->isDecoy(rank)) {
      m_numDecoyLowerHits[rank]++;
      m_numTotalDecoyLowerHits++; 
//      if (isSingleton) {
//	m_numTotalDecoyLowerSingletonHits++;
//      }
    }
  }
   
}

// logStats - write the stats in the log file. TO BE EXPANDED
void SpectraSTSearchTaskStats::logStats() {
  
  stringstream bss;
  
  bss << "Breakdown: ";
  bss.precision(0);
  for (int ch = 0; ch <= 4; ch++) {
    bss << '+' << ch + 1 << " = " << m_numTopHits[ch] << " ; ";
  }
  g_log->log("SEARCH STATS", bss.str());
  
  if (!m_doDecoyAnalysis) return;
   
  // decoy analysis
  if (m_numTotalDecoyLowerHits > 0) { // if there is any decoys at all
    stringstream decoyss;
    decoyss.precision(4);
    decoyss << "Decoy analysis:";
    
    for (unsigned int rank = 2; rank <= 10; rank++) {
      if (m_numLowerHits[rank] > 0) {
        decoyss << " D" << rank << " = " << fixed << ((double)m_numDecoyLowerHits[rank] / (double)m_numLowerHits[rank]) << ";";
      }
    }
    
    for (int ch = 0; ch <= 4; ch++) {
      if (m_numBadTopHits[ch] > 0) {
        decoyss << " D1bad(+" << ch + 1 << ") = " << fixed << ((double)m_numDecoyBadTopHits[ch] / (double)m_numBadTopHits[ch]) << ";";
      }
      if (m_numGoodTopHits[ch] > 0) {
        decoyss << " D1good(+" << ch + 1 << ") = " << fixed << ((double)m_numDecoyGoodTopHits[ch] / (double)m_numGoodTopHits[ch]) << ";";
      }
      if (m_numTopHits[ch] > 0) {
        decoyss << " D1(+" << ch + 1 << ") = " << fixed << ((double)m_numDecoyTopHits[ch] / (double)m_numTopHits[ch]) << ";";
      }
    }
    
    if (m_numTotalLowerHits > 0) {
      decoyss << " D2-10 = " << fixed << ((double)m_numTotalDecoyLowerHits / (double)m_numTotalLowerHits) << ";";
 //     if (m_numTotalLowerSingletonHits > 0) {
//	decoyss << " D2-10(SINGLETON) = " << fixed << ((double)m_numTotalDecoyLowerSingletonHits / (double)m_numTotalLowerSingletonHits) << ";";
//      }
    }

    g_log->log("DECOY STATS", decoyss.str());
    
    for (map<string, int>::iterator i = m_decoyCounts.begin(); i != m_decoyCounts.end(); i++) {
      if (i->second > 10) {
	stringstream badDecoyss;
	badDecoyss << "Frequent decoy hit: " << i->first << " (" << i->second << " times)";
	g_log->log("DECOY STATS", badDecoyss.str());
      }
    }
    
  }
  
}
    
SpectraSTSearchTaskStats* SpectraSTSearchTaskStats::aggregateStats(vector<SpectraSTSearchTaskStats*>& statsArray) {
   
  if (statsArray.empty()) {
    return (new SpectraSTSearchTaskStats(false)); // still return an object with zeros
  }
  
  bool doDecoyAnalysis = statsArray[0]->m_doDecoyAnalysis;
  // create new SpectraSTSearchTaskStats object to hold the sums of the stats
  SpectraSTSearchTaskStats* aggregate = new SpectraSTSearchTaskStats(doDecoyAnalysis);
  
  for (vector<SpectraSTSearchTaskStats*>::iterator i = statsArray.begin(); i != statsArray.end(); i++) {
    
    aggregate->m_numScans += (*i)->m_numScans;
    aggregate->m_numNotSelected += (*i)->m_numNotSelected;
    aggregate->m_numMissing += (*i)->m_numMissing;
    aggregate->m_numMS1 += (*i)->m_numMS1;
    aggregate->m_numFailedFilter += (*i)->m_numFailedFilter;
    aggregate->m_numSearched += (*i)->m_numSearched;
    aggregate->m_numLikelyGood += (*i)->m_numLikelyGood;
    
    for (int ch = 0; ch <= 5; ch++) {
      aggregate->m_numTopHits[ch] += (*i)->m_numTopHits[ch];
      aggregate->m_numBadTopHits[ch] += (*i)->m_numBadTopHits[ch];
      aggregate->m_numGoodTopHits[ch] += (*i)->m_numGoodTopHits[ch];

      aggregate->m_numDecoyTopHits[ch] += (*i)->m_numDecoyTopHits[ch];
      aggregate->m_numDecoyBadTopHits[ch] += (*i)->m_numDecoyBadTopHits[ch];
      aggregate->m_numDecoyGoodTopHits[ch] += (*i)->m_numDecoyGoodTopHits[ch];
    }
    
    for (int rank = 0; rank <= 10; rank++) {
      aggregate->m_numLowerHits[rank] += (*i)->m_numLowerHits[rank];
      aggregate->m_numDecoyLowerHits[rank] += (*i)->m_numDecoyLowerHits[rank];
    }
    
    aggregate->m_numTotalLowerHits += (*i)->m_numTotalLowerHits;
    aggregate->m_numTotalDecoyLowerHits += (*i)->m_numTotalDecoyLowerHits;
  
    for (map<string, int>::iterator decoy = (*i)->m_decoyCounts.begin(); decoy != (*i)->m_decoyCounts.end(); decoy++) {
      map<string, int>::iterator found = aggregate->m_decoyCounts.find(decoy->first);  
 
      if (found != aggregate->m_decoyCounts.end()) {
        found->second += decoy->second;
      } else {
        aggregate->m_decoyCounts[decoy->first] = decoy->second;
      }
    }
      
  }
  
  return (aggregate);
  
}
