#include "SpectraSTMspSearchTask.hpp"
#include "SpectraSTSearch.hpp"
#include "SpectraSTQuery.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"

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

/* Class: SpectraSTMspSearchTask
 * 
 * Subclass of SpectraSTSearchTask that handles .msp file types
 * 
 */
 
extern bool g_quiet;
extern bool g_verbose;
extern SpectraSTLog* g_log;

// constructor
SpectraSTMspSearchTask::SpectraSTMspSearchTask(vector<string>& searchFileNames, SpectraSTSearchParams& params, SpectraSTLib* lib) :
  SpectraSTSearchTask(searchFileNames, params, lib) {
}

// destructor
SpectraSTMspSearchTask::~SpectraSTMspSearchTask()
{
}

// search - run the searches
void SpectraSTMspSearchTask::search() {
  
  unsigned int numThreads = (unsigned int)(m_params.numThreadsUsed);
  
  if (numThreads > 1) {

    cout << "Multi-threaded search: Using " << numThreads << " threads." << endl;

    // Multi-threaded search   
    struct threadData* threadDataArray = new struct threadData[numThreads];
    
    for (unsigned int ti = 0; ti < numThreads; ti++) {
      threadDataArray[ti].searchTaskPtr = this;
      threadDataArray[ti].threadIndex = ti;
      for (unsigned int j = ti; j < m_searchFileNames.size(); j += numThreads) {
	threadDataArray[ti].fileIndices.push_back(j);
      }
    }      
    
#ifdef MSVC
      DWORD *pId = new DWORD[numThreads];
      HANDLE *threads = new HANDLE[numThreads];
#else    
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      void *status;
	
      pthread_t* threads = new pthread_t[numThreads];
#endif
      
    for (unsigned int ti = 0; ti < numThreads; ti++) {
      
#ifdef MSVC
      int returnCode = ti + 1;
      threads[ti] = CreateThread(NULL, 0, runSearchThread, (void*)&threadDataArray[ti], 0, NULL);
      if (!threads[ti]) {
	returnCode = 0;
      }
#else
      int returnCode = pthread_create(&threads[ti], &attr, runSearchThread, (void*)(&(threadDataArray[ti]))); 
#endif
      
      if (returnCode != 0) {
	stringstream msg;
	msg << "Cannot spawn new thread #" << ti << "; return code from pthread_create() is " << returnCode;
	g_log->error("MSP SEARCH", msg.str());
	g_log->crash();
      }

    }
    
#ifndef MSVC      
    pthread_attr_destroy(&attr);
#endif
    
    for (unsigned int ti = 0; ti < numThreads; ti++) { 
      
#ifdef MSVC
      int returnCode = WaitForSingleObject(threads[ti],INFINITE);
#else
      int returnCode = pthread_join(threads[ti], &status);
#endif
      
      if (returnCode != 0) {
	stringstream msg;
#ifdef MSVC
	msg << "Cannot join thread #" << ti << "; return code from WaitForSingleObject() is " << returnCode;
#else
	msg << "Cannot join thread #" << ti << "; return code from pthread_join() is " << returnCode;
#endif
	g_log->error("MSP SEARCH", msg.str());
	g_log->crash();
      }
	//printf("Main: completed join with thread %ld having a status of %ld \n",t,(long)status);
      
    }

    delete[] threadDataArray;
    delete[] threads;
    
  } else {
    
    // single-threaded search
      for (unsigned int n = 0; n < (unsigned int)m_searchFileNames.size(); n++) {
      searchOneFile(n, -1);
    }
  
  }
  
  
  logSearchStats("MSP SEARCH");
}

#ifdef MSVC
DWORD WINAPI SpectraSTMspSearchTask::runSearchThread(LPVOID threadArg) {
#else
void* SpectraSTMspSearchTask::runSearchThread(void* threadArg) {
#endif  
  
  struct threadData* threadData = (struct threadData*)threadArg;
  
  for (vector<unsigned int>::iterator i = threadData->fileIndices.begin(); i != threadData->fileIndices.end(); i++) {
  
    threadData->searchTaskPtr->searchOneFile(*i, (int)(threadData->threadIndex));
  
  }
  
  long ti = (long)(threadData->threadIndex);
  
#ifdef MSVC
  ExitThread(0);
#else
  pthread_exit((void*)ti); 
#endif
}


// searchOneFile - search one msp file
void SpectraSTMspSearchTask::searchOneFile(unsigned int fileIndex, int threadIndex) {
  
  string searchFileName(m_searchFileNames[fileIndex]);
  SpectraSTSearchTaskStats* stats = m_searchTaskStats[fileIndex];
  
  ifstream fin;
  if (!myFileOpen(fin, searchFileName)) {
    g_log->error("MSP SEARCH", "Cannot open MSP file \"" + searchFileName + " for reading query spectra. File Skipped.");
    return;
  }	
  
  m_outputs[fileIndex]->openFile();
  m_outputs[fileIndex]->printHeader();
  
  ProgressCount pc(!g_quiet && !g_verbose, 200);
  
  if (threadIndex == -1) {
    // single-threaded

    stringstream msg;
    msg << "Searching \"" << searchFileName << "\" " << "(" << fileIndex + 1 << " of " << m_searchFileNames.size() << ")";
    pc.start(msg.str());	
  
  } else {
    // multi-threaded, difficult to track % progress, just mark the start of a file here
    
    stringstream msg;
    msg << "Starting search of \"" << searchFileName << "\" by thread #" << threadIndex;
    cout << msg.str() << endl; // TODO: Is this thread-safe? Not really but no big deal.

  }
  
  string line("");
  string::size_type pos = 0;
  
  
  while (true) {
    
    if (line == "_EOF_") {
      // no more record
      if (threadIndex == -1) {
	pc.done();
      } else {
	// multi-threaded. Mark the end of a file
	stringstream msg;
	msg << "Finishing search of \"" << m_searchFileNames[fileIndex] << "\" by thread #" << threadIndex;
	cout << msg.str() << endl; 
      }
      
      m_outputs[fileIndex]->printFooter();
      m_outputs[fileIndex]->closeFile();
      return;
    }
    
    // skips over all lines until the line with Name: (will stop when it reaches either "Name:" or the end-of-file) 
    if (line.compare(0, 6, "Name: ") != 0) {	
      while (nextLine(fin, line, "Name: ", ""));
      if (line == "_EOF_") {
	// no more record
	if (threadIndex == -1) {
	  pc.done();
	} else {
	  // multi-threaded. Mark the end of a file
	  stringstream msg;
	  msg << "Finishing search of \"" << m_searchFileNames[fileIndex] << "\" by thread #" << threadIndex;
	  cout << msg.str() << endl; 
	}
	
        m_outputs[fileIndex]->printFooter();
        m_outputs[fileIndex]->closeFile();
	return;
      }
    }
    
    
    // line should now start with "Name: "
    string name = nextToken(line, 5, pos, "\r\n");
    
    bool selected;	
    if (!m_searchAll) {
      selected = isInSelectedList(name);
      // not to mess up the parsing, the following fields will still be 
      // parsed even though the record is not selected to be searched.
      // the search simply will not be initiated.
    } else {
      selected = true;
    }
    
    
    double mw = 0.0;
    double precursorMz = 0.0;
    string comments("");
    int charge = 0;
    
    // read the rest of the header fields until Num peaks:
    while (nextLine(fin, line, "Num peaks: ", "")) {
      if (line.compare(0, 3, "MW:") == 0) {
	mw = atof((nextToken(line, 3, pos, "\r\n")).c_str());				
      } else if (line.compare(0, 8, "Comment:") == 0) {
	comments = nextToken(line, 8, pos, "\r\n");
      } else if (line.compare(0, 12, "PrecursorMZ:") == 0) {
	precursorMz = atof((nextToken(line, 12, pos, "\r\n")).c_str());
      } else if (line.compare(0, 7, "Charge:") == 0) {
	charge = atoi((nextToken(line, 7, pos, "\r\n")).c_str());		
      } else if (line.compare(0, 9, "NumPeaks:") == 0) {
	// reach here because this .msp probably is converted from a .sptxt,
	// hack the line (pad one char) so that the number of peaks will be correctly parsed later
	line = " " + line;
	break;		
      } else if (line.compare(0, 7, "Status:") == 0 || 
		 line.compare(0, 9, "FullName:") == 0 || 
		 line.compare(0, 6, "LibID:") == 0) {
	// reach here because this .msp probably is converted from a .sptxt,
	// just ignore this line	
      } else if (line == "_EOF_" || line.compare(0, 5, "Name:") == 0) {
	// reach the end unexpectedly, or see another name field before the Num peaks field. 
	// ignore this incomplete record, and return
	cerr << "\nBadly formatted .msp file! Library creation truncated." << endl;	
	return;
      } else {
	cerr << "Unrecognized header field. Ignored." << endl;
      }
    }
    
    if (line == "_EOF_") {
      // no "Num peaks:" field. ignore this incomplete record, and return
      cerr << "\nBadly formatted .msp file! Library creation truncated." << endl;	
      m_outputs[fileIndex]->printFooter();
      m_outputs[fileIndex]->closeFile();

      return;
    }
    
    
    if (precursorMz < 0.0001) {
      // Precursor m/z not specified as a header field, try to find it in the comments
      string::size_type parentPos = comments.find("Parent", 0);
      if (parentPos != string::npos) {
	precursorMz = atof((nextToken(comments, parentPos + 7, parentPos, " \t\r\n")).c_str());	
      } else {
	if (charge != 0) {
	  // Charge is specified. Can calculate precursor m/z
	  precursorMz = (mw + charge * 1.00728) / (double)charge;	
	} else {
	  // nothing we can do, maybe the mw is really the precursor m/z
	  precursorMz = mw;
	}
      }
    }
    
    // line now starts with "Num peaks: "
    int numPeaks = atoi((nextToken(line, 10, pos, "\r\n")).c_str());
    
    SpectraSTPeakList* peakList = new SpectraSTPeakList(precursorMz, charge, numPeaks);
    peakList->setNoiseFilterThreshold(m_params.filterRemovePeakIntensityThreshold);
    
    // hack -- if the query file is an .sptxt, then put the binaryFileOffset as
    // part of the query name. This allows the HTML output file to link to the spectrum
    string::size_type binaryFileOffsetPos = comments.find("BinaryFileOffset", 0);
    if (binaryFileOffsetPos != string::npos) {
      string binaryFileOffset = nextToken(comments, binaryFileOffsetPos + 17, binaryFileOffsetPos, " \t\r\n");
      name += "." + binaryFileOffset + "." + binaryFileOffset + ".0";
    }
    
    SpectraSTQuery* query = new SpectraSTQuery(name, precursorMz, charge, comments, peakList);
    
    while (nextLine(fin, line, "Name: ", "")) { // will stop when it reaches the next "Name:" or end-of-file
      
      // here are the peaks		
      double mz = atof((nextToken(line, 0, pos)).c_str());
      // pos1 now stores the position of the first space after the mz
      float intensity = (float)(atof((nextToken(line, pos, pos)).c_str()));
      // pos2 now stores the position of the first space after the intensity
      
      // annotation has quotes around it, remove them by adding the quote char to the skipover and delimiter strings passed into nextToken
      string annotation = nextToken(line, pos, pos, "\"\r\n", "\"\t");
      
      string info("");
      
      string::size_type spacePos = annotation.find_first_of(" \t", 0);
      string::size_type dummyPos = 0;
      if (spacePos != string::npos) {
	info = nextToken(annotation, spacePos + 1, dummyPos, "\r\n", " \t");
	annotation = annotation.substr(0, spacePos);
      }
      // annotation will get an empty string if there's no annotation
      peakList->insertForSearch(mz, intensity, annotation);
      
      
    }	
    
    if (!selected) {
      // not selected to be searched, or bad spectra, ignore
      delete query;
      stats->m_numNotSelected++;
      
    } else if (!peakList->passFilter(m_params)) {
      
      delete query;
      stats->m_numFailedFilter++;
      
    } else {
      
//      if (m_params.filterMaxPeaksUsed <= 50) {
//        double retained = peakList->simplify(m_params.filterMaxPeaksUsed, m_params.filterMaxDynamicRange, true, true);
//      } else {
//        double retained = peakList->simplify(m_params.filterMaxPeaksUsed, m_params.filterMaxDynamicRange);
//      }
        
      FileName fn;
      parseFileName(searchFileName, fn);
      
      // create the search based on what is read, then search
      SpectraSTSearch* s = new SpectraSTSearch(query, m_params, m_outputs[fileIndex]);
      s->search(m_lib);
      
      stats->m_numSearched++; // counting searches in all msp files
      m_searchTaskStats[fileIndex]->processSearchResult(s);
      
      if (threadIndex == -1) pc.increment();
      
      // print search result
      s->print();
      delete s;
    }
    
  }
  m_outputs[fileIndex]->printFooter();
  m_outputs[fileIndex]->closeFile();

}

