#ifndef SPECTRASTMZXMLSEARCHTASK_HPP_
#define SPECTRASTMZXMLSEARCHTASK_HPP_

#include "SpectraSTSearchTask.hpp"
#include "SpectraSTSearchParams.hpp"
#include "SpectraSTLib.hpp"
#include "FileUtils.hpp"

#ifdef STANDALONE_LINUX
#include "SpectraST_cramp.hpp"
#else
#include "Parsers/mzParser/cramp.hpp"
#endif

#include <string>
#include <map>
#include <vector>

#ifdef __MINGW__
#define MSVC
#endif

#ifdef MSVC
#include "windows.h"
#else
#include <pthread.h>
#endif

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


/* Class: SpectraSTMzXMLSearchTask
 * 
 * Subclass of SpectraSTSearchTask that handles.mzXML file types. Note that because of
 * the availability of a random access index at the end of the mzXML, conveniently one can
 * sort the spectra by ascending precursor m/z and search them in that order to take advantage
 * of caching.
 * 
 */
 
class SpectraSTMzXMLSearchTask : public SpectraSTSearchTask {
  
public:
  SpectraSTMzXMLSearchTask(vector<string>& searchFileNames, SpectraSTSearchParams& params, SpectraSTLib* lib);
  virtual ~SpectraSTMzXMLSearchTask();
  virtual void search();
  
    
#ifdef MSVC
  static DWORD WINAPI runSearchThread(LPVOID threadArg);
#else
  static void* runSearchThread(void* threadArg);
#endif

private:
  
  // m_files - keeps the open files as tuples (fileName, cRamp*)
  vector<pair<FileName, cRamp*> > m_files;
  
  // m_scans - all the scans to be searched, as (fileIndex, scanInfo*). 
  vector<pair<unsigned int, rampScanInfo*> > m_scans;
  
  // m_batchBoundaries - the index boundaries of the batches. Say m_batchBoundaries contains 0, 9, 19, then
  // there are two batches, files 0 to 8, and files 9 to 18. The last number in m_batchBoundaries is the
  // total number of files (== m_searchFileNames.size())
  vector<unsigned int> m_batchBoundaries;
  
  // private method for sorting queries by precursor m/z before search 
  void prepareSortedSearch(unsigned int batch);
  
  // private method for searching one file (threadIndex = -1 means single-threaded)
  void searchOneFile(unsigned int fileIndex, int threadIndex);
  
  // private method for searching one query
  void searchOneScan(unsigned int fileIndex, rampScanInfo* scanInfo);
  
  // comparator method for sorting
  static bool sortRampScanInfoPtrsByPrecursorMzAsc(pair<unsigned int, rampScanInfo*> a, pair<unsigned int, rampScanInfo*> b);
  
  bool m_isMzData;
  
   // Fingerprinting
  vector<float> m_searchFingerprintIndexedByLibID; 
  vector<float> m_searchFingerprintIndexedBySearchID; // Index starts with zero. If dot IsLikelyGood, the value is the matched libID; else, the value is (number of lib entry + 1)
  vector<float> m_bootstrapSupport;  
  
#ifdef MSVC
  HANDLE m_fingerprintMutex;
#else
  pthread_mutex_t* m_fingerprintMutex;
#endif
  
  void printFingerprintingSummary();
  void calcFingerprint(bool doBootstrapping, ofstream& fout);
  void updateFingerprint(SpectraSTSearch* s);
  void doBootstrapping(); 
  // END Fingerprinting
};

#endif /*SPECTRASTMZXMLSEARCHTASK_HPP_*/
