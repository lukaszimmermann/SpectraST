#include "SpectraSTMgfSearchTask.hpp"
#include "SpectraSTSearch.hpp"
#include "SpectraSTQuery.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"

#include <sstream>

extern bool g_quiet;
extern bool g_verbose;
extern SpectraSTLog* g_log;

// constructor
SpectraSTMgfSearchTask::SpectraSTMgfSearchTask(vector<string>& searchFileNames, SpectraSTSearchParams& params, SpectraSTLib* lib) :
  SpectraSTSearchTask(searchFileNames, params, lib) {
}

// destructor
SpectraSTMgfSearchTask::~SpectraSTMgfSearchTask() {
}

// search - run the searches
void SpectraSTMgfSearchTask::search() {
  
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
	g_log->error("MGF SEARCH", msg.str());
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
	g_log->error("MGF SEARCH", msg.str());
	g_log->crash();
      }
	//printf("Main: completed join with thread %ld having a status of %ld \n",t,(long)status);
      
    }

    delete[] threadDataArray;
    delete[] threads;
    
  } else {
  
    // Single-threaded search
    for (unsigned int n = 0; n < (unsigned int)m_searchFileNames.size(); n++) {
      searchOneFile(n, -1);
    }
    
  }
  
  logSearchStats("MGF SEARCH");
}

#ifdef MSVC
DWORD WINAPI SpectraSTMgfSearchTask::runSearchThread(LPVOID threadArg) {
#else
void* SpectraSTMgfSearchTask::runSearchThread(void* threadArg) {
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

// searchOneFile - search one mgf file
void SpectraSTMgfSearchTask::searchOneFile(unsigned int fileIndex, int threadIndex) {
  
  string searchFileName(m_searchFileNames[fileIndex]);
  SpectraSTSearchTaskStats* stats = m_searchTaskStats[fileIndex];
  
  ifstream fin;
  if (!myFileOpen(fin, searchFileName)) {
    g_log->error("MGF SEARCH", "Cannot open MGF file \"" + searchFileName + " for reading query spectra. File Skipped.");
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
    if (line.compare(0, 10, "BEGIN IONS") != 0) {
      while (nextLine(fin, line, "BEGIN IONS", ""));
      if (line == "_EOF_") {
	// no more record
	if (threadIndex == -1) {
	  pc.done();
	} else {
	  // multi-threaded. Mark the end of a file
	  if (g_verbose) {
	    stringstream msg;
	    msg << "Finishing search of \"" << m_searchFileNames[fileIndex] << "\" by thread #" << threadIndex;
	    cout << msg.str() << endl; 
	  }
	}

	m_outputs[fileIndex]->printFooter();
	m_outputs[fileIndex]->closeFile();
	return;
      }
    }
    
    double precursorMz = 0.0;
    string title("");
    string comments("");
    int charge = 0;
    
    
    // read the rest of the header fields until TITLE:
    while (nextLine(fin, line, "END IONS", "")) {
      
      // cerr << line << endl;
      
      if (line.find('=') == string::npos) {
	// no more headers
	break; 
      }

      if (line.compare(0, 6, "TITLE=") == 0) {
	title = nextToken(line, 6, pos, " \t\r\n"); 
      } else if (line.compare(0, 8, "PEPMASS=") == 0) {
	precursorMz = atof((nextToken(line, 8, pos," \r\n")).c_str());
	double precursorIntensity = atof((nextToken(line, pos, pos," \r\n")).c_str());
      } else if (line.compare(0, 7, "CHARGE=") == 0) {
	charge = atoi((nextToken(line, 7, pos, " \t\r\n", "+")).c_str());
      } else if (line == "_EOF_" || line.compare(0, 10, "BEGIN IONS") == 0) {
	cerr << "\nBadly formatted .mgf file! Search task truncated." << endl;
	m_outputs[fileIndex]->printFooter();
	m_outputs[fileIndex]->closeFile();
	return;

      } else {
	// some other header, put the whole thing in comments
	comments += nextToken(line, 0, pos, " \t\r\n") + " ";
      }
    }

    bool hasNoPeaks = false;
    if (line.find("END IONS") != string::npos) {
      // no peaks?
      hasNoPeaks = true;
    }

    bool selected;
    if (!m_searchAll) {
      selected = isInSelectedList(title);
      // not to mess up the parsing, the following fields will still be
      // parsed even though the record is not selected to be searched.
      // the search simply will not be initiated.
    } else {
      selected = true;
    }

    bool hasNoPrecursorMz = false;
    if (precursorMz < 0.0001) {
     // no mass, this can't be searched
      hasNoPrecursorMz = true;
    }

    if (hasNoPeaks) {
      
	stats->m_numFailedFilter++;
	
    } else {

      int numPeaks = 0;
      
      SpectraSTPeakList* peakList = new SpectraSTPeakList(precursorMz, charge);
      peakList->setNoiseFilterThreshold(m_params.filterRemovePeakIntensityThreshold);
      
      SpectraSTQuery* query = new SpectraSTQuery(title, precursorMz, charge, comments, peakList);
      
      // line should contain the first peak
      
      do { // will stop when it reaches the next "END IONS" or end-of-file
	
	// here are the peaks
	double mz = atof((nextToken(line, 0, pos)).c_str());
	// pos1 now stores the position of the first space after the mz
	float intensity = atof((nextToken(line, pos, pos)).c_str());
	// pos2 now stores the position of the first space after the intensity
	peakList->insertForSearch(mz, intensity, "");
	
      } while (nextLine(fin, line, "END IONS", ""));
    
      if (!selected) {
	
	// not selected to be searched, or bad spectra, ignore
	delete query;
	stats->m_numNotSelected++;
	
      } else if (hasNoPrecursorMz) {
	
	delete query;
	stats->m_numFailedFilter++;
	
      } else if (!peakList->passFilter(m_params)) {
	
	delete query;
	stats->m_numFailedFilter++;
	  
      } else {
    
	FileName fn;
	parseFileName(searchFileName, fn);
	
	// create the search based on what is read, then search
	SpectraSTSearch* s = new SpectraSTSearch(query, m_params, m_outputs[fileIndex]);
	s->search(m_lib);
	
	stats->m_numSearched++; // counting searches in all mgf files
	m_searchTaskStats[fileIndex]->processSearchResult(s);
	
	if (threadIndex == -1) pc.increment();

	// print search result
	s->print();
	delete (s);
      }
    }

  }
  m_outputs[fileIndex]->printFooter();
  m_outputs[fileIndex]->closeFile();

}


