#include "SpectraSTLog.hpp"
#include <iostream>
#include <stdlib.h>

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

/* Class: SpectraSTLog
 * 
 * Manages a log file that takes care of error/warning reporting and various logging. 
 * 
 */

// Constructor
SpectraSTLog::SpectraSTLog(string logFileName) :
  m_logFileName(logFileName),
  m_fout(),
  m_numError(0),
  m_numWarning(0),
  m_good(false),
  m_errors(),
  m_mutex(NULL),
  m_warnings() {

  m_fout.open(logFileName.c_str(), ios::app);
  if (m_fout.good()) {
    m_good = true;
  } else {
    m_numError++;
    m_errors.push_back("LOGGING: Cannot open log file \"" + logFileName + "\" for writing. No logging for current session.");  
  }  

#ifdef MSVC
  m_mutex = CreateMutex( 
			    NULL,              // default security attributes
			    FALSE,             // initially not owned
			    NULL);             // unnamed mutex
  
  if (m_mutex == NULL) 
    {
      printf("CreateMutex error: %d\n", (int)GetLastError());
      exit(1);
    }
#else
    m_mutex = new pthread_mutex_t();
    pthread_mutex_init(m_mutex, NULL);
#endif

}

// Destructor
SpectraSTLog::~SpectraSTLog() {
  m_fout.close();
  
  if (m_mutex) {
#ifdef MSVC
    CloseHandle(m_mutex);
#else
    pthread_mutex_destroy(m_mutex);
    delete (m_mutex);
#endif
  }
}

// log - logs a general message
void SpectraSTLog::log(string msg) {
  if (m_good) {
    if (m_mutex) {
#ifdef MSVC
      WaitForSingleObject( m_mutex,    // handle to mutex
			   INFINITE);      //#include "windows.h"
#else
      pthread_mutex_lock(m_mutex);
#endif      
    }
    m_fout << msg << endl;   
    if (m_mutex) {
#ifdef MSVC
      ReleaseMutex( m_mutex );      //#include "windows.h"
#else
      pthread_mutex_unlock( m_mutex );//#include <pthread.h>
#endif
    }
  }
}

// log - logs a message with a tag
void SpectraSTLog::log(string tag, string msg) {
  if (m_good) { 
    if (m_mutex) {
#ifdef MSVC
      WaitForSingleObject( m_mutex,    // handle to mutex
			   INFINITE);      //#include "windows.h"
#else
      pthread_mutex_lock(m_mutex);
#endif 
    }
    m_fout << tag << ": " << msg << endl;
    if (m_mutex) {
#ifdef MSVC
      ReleaseMutex( m_mutex );      //#include "windows.h"
#else
      pthread_mutex_unlock( m_mutex );//#include <pthread.h>
#endif
    }
  }
}

// log - logs a message with a tag and a stamp (e.g. a time stamp)
void SpectraSTLog::log(string tag, string msg, string stamp) {
  if (m_good) {
    if (m_mutex) {
#ifdef MSVC
      WaitForSingleObject( m_mutex,    // handle to mutex
			   INFINITE);      //#include "windows.h"
#else
      pthread_mutex_lock(m_mutex);
#endif      
    }
    m_fout << tag << ": (" << stamp << ") " << msg << endl; 
    if (m_mutex) {
#ifdef MSVC
      ReleaseMutex( m_mutex );      //#include "windows.h"
#else
      pthread_mutex_unlock( m_mutex );//#include <pthread.h>
#endif
    }
  }
}

// error - logs an error message. Difference between error() and log() is that 
// 1. errors will be prefixed with the word "ERROR" in the log entry.
// 2. errors will be recorded and all error messages can be dumped later on
void SpectraSTLog::error(string tag, string msg) {
  if (m_mutex) {
#ifdef MSVC
    WaitForSingleObject( m_mutex,    // handle to mutex
			 INFINITE);      //#include "windows.h"
#else
    pthread_mutex_lock(m_mutex);
#endif      
  }
  if (m_good) {
    m_fout << "ERROR " << tag << ": " << msg << endl;
  }
  m_errors.push_back(tag + ": " + msg);
  m_numError++;
  
  if (m_mutex) {
#ifdef MSVC
    ReleaseMutex( m_mutex );      //#include "windows.h"
#else
    pthread_mutex_unlock( m_mutex );//#include <pthread.h>
#endif
  }
}

// warning - logs a warning message.
void SpectraSTLog::warning(string tag, string msg) {
  if (m_mutex) {
#ifdef MSVC
    WaitForSingleObject( m_mutex,    // handle to mutex
			 INFINITE);      //#include "windows.h"
#else
    pthread_mutex_lock(m_mutex);
#endif      
  }
  if (m_good) {
    m_fout << "WARNING " << tag << ": " << msg << endl;
  }
  m_warnings.push_back(tag + ": " + msg);
  m_numWarning++;
  if (m_mutex) {
#ifdef MSVC
    ReleaseMutex( m_mutex );      //#include "windows.h"
#else
    pthread_mutex_unlock( m_mutex );//#include <pthread.h>
#endif
  }
}

// printErrors - dumps all errors to console.
void SpectraSTLog::printErrors() {
  if (m_mutex) {
#ifdef MSVC
    WaitForSingleObject( m_mutex,    // handle to mutex
			 INFINITE);      //#include "windows.h"
#else
    pthread_mutex_lock(m_mutex);
#endif   
  }
  for (vector<string>::iterator i = m_errors.begin(); i != m_errors.end(); i++) {
    cout << (*i) << endl;
  }
  if (m_mutex) {
#ifdef MSVC
    ReleaseMutex( m_mutex );      //#include "windows.h"
#else
    pthread_mutex_unlock( m_mutex );//#include <pthread.h>
#endif
  }
}

// printWarnings - dumps all warnings to console.
void SpectraSTLog::printWarnings() {
  if (m_mutex) {
#ifdef MSVC
    WaitForSingleObject( m_mutex,    // handle to mutex
			 INFINITE);      //#include "windows.h"
#else
    pthread_mutex_lock(m_mutex);
#endif
  }   
  for (vector<string>::iterator i = m_warnings.begin(); i != m_warnings.end(); i++) {
    cout << "WARNING -- " << (*i) << endl;
  }
  if (m_mutex) {
#ifdef MSVC
    ReleaseMutex( m_mutex );      //#include "windows.h"
#else
    pthread_mutex_unlock( m_mutex );//#include <pthread.h>
#endif
  }
}

// crash - called when a fatal enough error has occurred and the caller has no intention of carrying on.
// will dump all error messages before exiting.
// Just a nicer way of crashing.
void SpectraSTLog::crash() {
  if (m_mutex) {
#ifdef MSVC
    WaitForSingleObject( m_mutex,    // handle to mutex
			 INFINITE);      //#include "windows.h"
#else
    pthread_mutex_lock(m_mutex);
#endif  
  } 
  cerr << "\t==== FATAL ERROR. Exiting immediately. ====" << endl;
  cerr << "\tError trace :" << endl;
  for (vector<string>::iterator i = m_errors.begin(); i != m_errors.end(); i++) {
    cerr << '\t' << (*i) << endl;
  }
  cerr << "\t===========================================" << endl;
  if (m_mutex) {
#ifdef MSVC
    ReleaseMutex( m_mutex );      //#include "windows.h"
#else
    pthread_mutex_unlock( m_mutex );//#include <pthread.h>
#endif
  }
  exit (1);

}
