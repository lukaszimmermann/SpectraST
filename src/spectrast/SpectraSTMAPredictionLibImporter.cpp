#include "SpectraSTMAPredictionLibImporter.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"
#include <iostream>
#include <sstream>
#include <stdlib.h>


/*
 * Program : SpectraSTLog
 * Author: HU Yingwei
 * Date: 2013.01.13
 */

extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog* g_log;




//constructor -will open the first file
SpectraSTMAPredictionLibImporter::SpectraSTMAPredictionLibImporter(vector<string>& impFileNames, SpectraSTLib* lib, SpectraSTCreateParams& params):
  SpectraSTLibImporter(impFileNames, lib, params){
}

//desctructor
SpectraSTMAPredictionLibImporter::~SpectraSTMAPredictionLibImporter()
{

}

//import - prints the preamble, then loops over all files and import them one by one
void SpectraSTMAPredictionLibImporter::import(){
  for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
    string fullName(*i);
    makeFullPath(fullName);
    string quoted("\"" + fullName + "\"");
    string desc = m_params.constructDescrStr(quoted,".map");
    m_preamble.push_back(desc);// add full source filename to library preamble "IMPORT From: ..."
    
    m_lib->writePreamble(m_preamble);
    
    for (vector<string>::iterator i = m_impFileNames.begin(); i!= m_impFileNames.end(); i++) {
      readFromFile(*i);
    }
  }  
}

//readFromCurFile - reads one .map file
void SpectraSTMAPredictionLibImporter::readFromFile(string& impFileName) {
  
  ifstream fin;
  
  if (!myFileOpen(fin,impFileName)) {
    g_log->error("MAP IMPORT","Cannot open " + impFileName + "for reading map-format (MassAnalyzer 2.1 MS2 Prediction");
    return;
  }
  
  g_log->log("MAP IMPORT", "Importing .map file \"" + impFileName + "\".");
  
  string line("");
  string::size_type pos = 0;
  
  if (g_verbose){
    cout << "\nImporting spectra from .map (MassAnalyzer 2.1 MS2 Prediction) file..." << endl;
  }
  
  //start the progress count
  ProgressCount pc(!g_quiet && !g_verbose, 1000);
  pc.start("\nImporting spectra from .map (MassAnalyzer 2.1 MS2 Prediction) file");
  
  int absLineNum  = 0;
  int relLineNum = 0;
  int totalEntriesNum = 0;
  string peptide("");
  string description("");
  int charge;
  string fragmentType(""); 
  int fragmentTypeID = 0;
  int collisionEnergy = 0;
  double reactionTime = 0;//second(unit)
  double isolationWidth = 0; 
  int instrumentModelID = 0;
  string instrumentModel("");
  double resolution = 0;// at m/z 400
  double scanLowerBound = 0;
  double scanUpperBound = 0;
  int peakNum = 0;
  vector<pair<double, double> > peakList;
  
  while (nextLine(fin, line)) {
    absLineNum++;
    //cout << "DEBUG:LINENUM="<< absLineNum << endl;
    if (absLineNum == 1) {
      totalEntriesNum = atoi(line.c_str());
      //cout << "DEBUG:MAPImporter: Total entry number:" << totalEntriesNum << endl;
    } else {
      relLineNum++;
      //cout << "DEBUG:MAPImporter relLINE=" << relLineNum << "\t" << line << endl;
      if (relLineNum == 1) {
	peptide = line;
 	// cout << "DEBUG:MAPImporter1:peptide:" << peptide << endl;
      } else if (relLineNum == 2) {
	description = line;
	//cout << "DEBUG:MAPImporter:2:description:" << description << endl;
      } else if (relLineNum == 3) {
	charge = atoi(line.c_str());
	//cout << "DEBUG:MAPImporter:3:charge:" << charge << endl;
      } else if (relLineNum == 4) {
	fragmentTypeID = atoi(line.c_str());
	if (fragmentTypeID == 1) {
	  fragmentType = "CID";
	} else if (fragmentTypeID == 2) {
	  fragmentType = "ETD";
	} else if (fragmentTypeID == 3) {
	  fragmentType = "ETD-SA";
	}
	//cout << "DEBUG:MAPImporter:4:fragmentType:" << fragmentType << endl;
      } else if (relLineNum == 5) {
	collisionEnergy = atoi(line.c_str());
	//cout << "DEBUG:MAPImporter:5:collisionEnergy:" << collisionEnergy << endl; 
      } else if (relLineNum == 6) {
	reactionTime = atof(line.c_str());
	//cout << "DEBUG:MAPImporter:6:reactionTime:" << reactionTime << endl;
      } else if (relLineNum == 7) {
	isolationWidth = atof(line.c_str());
	//cout << "DEBUG:MAPImporter:7:isolationWidth:" << isolationWidth << endl;
      } else if (relLineNum == 8) {
	instrumentModelID = atoi(line.c_str());
	if (instrumentModelID == 1) {
	  instrumentModel = "LCQ";
	} else if (instrumentModelID == 2) {
	  instrumentModel = "LTQ";
	} else if (instrumentModelID == 3) {
	  instrumentModel = "Orbitrap";
	} else if (instrumentModelID == 4) {
	  instrumentModel = "LTQFT";
	} else if (instrumentModelID == 5) {
	  instrumentModel = "QTOF";
	}
	//cout << "DEBUG:MAPImporter:8:instrumentModel:" << instrumentModel << endl;
      } else if (relLineNum == 9) {
	resolution = atof(line.c_str());
	//cout << "DEBUG:MAPImporter:9:resolution" << resolution << endl; 
      } else if (relLineNum == 10) {
	string::size_type pos = 0;
	string lower = nextToken(line, pos, pos);
	string upper = nextToken(line, pos + 1, pos);
	scanLowerBound = atof(lower.c_str());
	scanUpperBound = atof(upper.c_str());
	//cout << "DEBUG:MAPImporter:10:scanLowerBound:" << scanLowerBound << endl;
	//cout << "DEBUG:MAPImporter:10:scanUpperBound:" << scanUpperBound << endl;
      } else if (relLineNum == 11) {
	peakNum = atoi(line.c_str());
	//cout << "DEBUG:MAPImporter:11:peakNum:" << peakNum << endl;
      } else if (relLineNum >= 12) {
	const char* str = line.c_str();
	if ((isupper(str[0])) || (line.compare(0,5,"n[43]") == 0)) {
	  
	  //cerr << "debug:MassAnalyzerPredictionImporter::readfromfile::last_peptide:" << peptide << endl;
	  //cerr << "debug:MassAnalyzerPredictionImporter::readfromfile::cur_line" << line << endl;
	  //cerr << "debug:MassAnalyzerPredictionImporter::readfromfile::str[0]" << str[0] << endl;	  
	  //cerr << "debug:MassAnalyzerPredictionImporter::readfromfile::isupper(str[0])" << isupper(str[0]) << endl;
	  
	  makeLibEntry(peptide, description, charge, fragmentType, collisionEnergy, reactionTime, isolationWidth, instrumentModel, resolution, scanLowerBound, scanUpperBound, peakNum, peakList);
	  //cerr << ">>TEST1" << endl;
	  relLineNum = 1;
	  peakList.clear();
	  pc.increment();//counter steps one
	  peptide = line;//a new entry start
	  //cout << "DEBUG:MAPImporter:1:peptide:" << peptide << endl;
	} else {
	  string::size_type pos = 0;
	  string mz = nextToken(line, pos, pos);
	  //cerr << "pos1=" << pos << endl; 
	  double mzVal = atof(mz.c_str());
	  string intensity = nextToken(line, pos + 1, pos);
	  if (intensity.empty()) {
	    continue;
	  }
	  //cerr << "pos2=" << pos << endl;
	  double intensityVal = atof(intensity.c_str());
	  //cout << "DEBUG:MAPImporter:mz:" << mzVal << ",intensity:" << intensityVal << endl;  
	  peakList.push_back(pair<double, double>(mzVal, intensityVal));
	  //cerr << "debug:MassAnalyzerPredictionImporter::readfromfile:showPeaks:";
	  //cerr << "mz:" << mzVal << " " << "in:" << intensityVal << endl;
	}
      }
    }
    
  } //for while
  //Now _EOF_ of the input file, complete the last entry
  makeLibEntry(peptide, description, charge,
	       fragmentType, collisionEnergy, reactionTime, isolationWidth, 
	       instrumentModel, resolution, scanLowerBound, scanUpperBound, 
	       peakNum, peakList);
  peakList.clear();
  pc.increment();//counter steps one
  
}

void SpectraSTMAPredictionLibImporter::makeLibEntry(string peptide, string description, int charge, string fragmentType, int collisionEnergy,
				double reactionTime, double isolationWidth, string instrumentModel, double resolution,
				double scanLowerBound, double scanUpperBound, int peakNum, vector<pair<double,double> >& peaks) {
  //cout << "makeLibEntry" << endl;
  double mw = 0.0;
  double precursorMz = 0.0;
  string comments("Spec=MassAnalyzer");
  
  string sequence("");
  string::size_type pos = 0;
  sequence = nextToken(peptide,pos,pos,"(");
  //cout << "DEBUG:MAPImporter: sequence:" << sequence << endl;
  vector<string> mods;
  string mod;
  while(!((mod = nextToken(peptide,pos+1,pos,")","(")).empty())){
    //cout << "mod:" << mod << endl;
    
    mods.push_back(mod);
  }
  
  //In MassAnalyzer, U: C[160], O: M[147], J: C[161]
  for (int i = 0; i < sequence.length(); i++) {
  if (sequence[i] == 'J') {
      stringstream ss;
      ss << "C" << (i + 1) << "+Carboxymethyl";
      mods.push_back(ss.str());
//       peptide[i] = 'C';
      sequence[i] = 'C';
    } else if (peptide[i] == 'U') {
      stringstream ss;
      ss << "C" << (i + 1) << "+Carbamidomethyl";
      mods.push_back(ss.str());
//       peptide[i] = 'C';
      sequence[i] = 'C';
    } else if (peptide[i] == 'O') {
       stringstream ss;
      ss << "M" << (i + 1) << "+Oxidation";
      mods.push_back(ss.str());
//       peptide[i] = 'M';
      sequence[i] = 'M';
    }    
  }

  string mspModStr("");
  if (!makeMspModStr(mods, mspModStr)) {
    // something wrong with parsing this mods
    // skip this entry entirely
    return;
  }
    
  //cout << "debug:makeLibEntry:mspModStr=" << mspModStr << endl;
  //note: create peptide
  //g_log->log("DEBUG:"+ mspModStr);
  Peptide* p = createPeptide(sequence, charge, mspModStr, "", "MAP");
  // Peptide* p = new Peptide(sequence, charge, mspModStr);
  precursorMz = p->monoisotopicMZ();
  //create peakList
  //SpectraSTPeakList* peakList = new SpectraSTPeakList(precursorMz,charge,peakNum);
  
  //create libEntry
  SpectraSTLibEntry* entry = new SpectraSTLibEntry(p, comments, "NORMAL", NULL);
  vector<pair<double,double> >::iterator it;
  
  //entry->getPeakList()->setNoiseFilterThreshold(m_params.rawSpectraNoiseThreshold);
  
  for (it = peaks.begin(); it != peaks.end(); it ++) {
    double mz = (*it).first;
    float intensity = float((*it).second);
    //cerr << "debug:MassAnalyzerPredictionImporter::makeLibentry:";
    //cerr << "mz=" << mz << "\t" << "in=" << intensity << endl; 
    entry->getPeakList()->insert(mz, intensity, "", "");
  }
  
  entry->getPeakList()->normalizeTo(10000, m_params.rawSpectraMaxDynamicRange);
  //add annotation for each peak of fragment
  
  entry->getPeakList()->annotate();
  //if (m_params.flattenAllPeaks){
  //  entry->getPeakList()->flattenAllPeaks();
  // }
  
  //cout << "insert entry" << endl;
  if (insertOneEntry(entry, "MAP")) {
    m_count++;
  }

  //cout << "del entry" << endl;
  delete (entry);
  
}
 
bool SpectraSTMAPredictionLibImporter::makeMspModStr(vector<string>& mods, string& mspModStr) {
  
  if (mods.size() == 0) {
    mspModStr = "0";
    return (true);
  }
  
  vector<string> mspMods;
  for(vector<string>::iterator it = mods.begin(); it != mods.end(); it++) {
    string::size_type pos = 0;
    string aa = nextToken(*it, pos, pos, "0123456789");
    string site = nextToken(*it, pos, pos, "+-");
    string pn = (*it).substr(pos, 1);
    string modType = nextToken(*it, pos + 1, pos);//modType (e.g.GnFB) or mod mass shift (e.g. 18)
    //check the aa--amino acid
    //cerr << "##" << aa << "\t" << site << "\t" << pn << "\t" << modType << endl;
    if (aa.length() != 1) {     
      string errorMsg = "Cannot figure out correct modification for " + (*it);
      g_log->error("MAP IMPORT",errorMsg);
      return (false);
    } 
    
    //check the site
    bool siteIsDigit = true;
    for (int i = 0; i < site.size(); i++) {
      if ((site.at(i) > '9') || (site.at(i) < '0')) {
	siteIsDigit = false;
	break;
      }
    }
    if (!siteIsDigit) {
      string errorMsg = "Cannnot figure out correct site for " + (*it);
      g_log->error("MAP IMPORT", errorMsg);
      return (false);
    }
    
    int correctSite = atoi(site.c_str()) - 1;
    stringstream sitess;
    sitess << correctSite;
    
    //check mod type or mod mass shift
    bool modTypeIsDigit = true;
    // for (int i = 0; i < modType.size(); i++){
    if ((modType.at(0) > '9' ) || (modType.at(0) < '0')) {
      modTypeIsDigit = false;
      string modName = modType;
      if (modType == "Methylation") {
	modName = "Methyl";
      } else if (modType == "Carbamidomethyl") {
	modName = "Carbamidomethyl";
      }

      mspMods.push_back(sitess.str() + "," + aa + "," + modName);
      continue;

    } 

    // }
    if (modTypeIsDigit) {
      int digitModType = int(atof(modType.c_str()) + 0.5);
      if (digitModType == 42) {
	if (site == "1") {
	  mspMods.push_back("-1," + aa + ",Acetyl");
	  continue;
	} else {
	  int correctSite = atoi(site.c_str()) - 1;
	  stringstream sitess;
	  sitess << correctSite;
	  mspMods.push_back(sitess.str() + "," + aa + ",Acetyl");
	  continue;
	}
      } else {
	string modName = "";
	if (digitModType == 1) {
	  modName = "Deamidated";
	} else if (digitModType == 14) {
	  modName = "Methyl";
	} else if (digitModType == 28) {
	  modName = "Dimethyl";
	} else if (digitModType == 80) {
	  modName = "Phospho";
	} else if (digitModType == 56) {
	  modName = "Propionyl";
	} else if (digitModType == 70) {
	  modName = "Methyl_Propionyl";
	}
	if (modName != "") {
	  mspMods.push_back(sitess.str() + "," + aa + "," + modName);
	  continue;
	}
      }
      //This is the case for unknown modification mass
      string errorMsg = "Unknown modification, not supported yet: " + (*it);
      g_log->error("MAP IMPORT", errorMsg);
      return (false);
    }
    
    //This is the case for known named modification
    mspMods.push_back(sitess.str() + "," + aa + "," + modType);
  }
  
  int modNum = mspMods.size();
  stringstream mspModss;
  mspModss << modNum ;
  for(vector<string>::iterator it= mspMods.begin(); it != mspMods.end(); it++) {
    mspModss << "/" << (*it) ;
  }
  mspModStr = mspModss.str();
  return (true);
}


