#include "Analyte.hpp"

#include "FileUtils.hpp"

#include <sstream>
#include <string>
#include <vector>
#include <map>

/*

Library       : Analyte
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

/* Class: Analyte
 * 
 */

// static members initialization. NOTE that the static method Peptide::defaultTable() 
// will be called the first one a Peptide object is constructed to set up these tables. If
// the user wishes to have non-default tables, he can supply these tables before calling any methods. 
// (See defaultTable() for an explanation) 
map<string, double>* Analyte::elementAverageMassTable = NULL;
map<string, double>* Analyte::groupAverageMassTable = NULL;
map<string, double>* Analyte::elementMonoisotopicMassTable = NULL;
map<string, double>* Analyte::groupMonoisotopicMassTable = NULL;

// CONSTRUCTORS, DESTRUCTOR AND ASSIGNMENT OPERATOR
	
Analyte::Analyte() {
		
}

// Destructor
Analyte::~Analyte() {
	
}

void Analyte::defaultTables() {
	
  Analyte::elementAverageMassTable = new map<string, double>();
  
  (*elementAverageMassTable)["C"] = 12.0107;
  (*elementAverageMassTable)["H"] = 1.00794;
  (*elementAverageMassTable)["O"] = 15.9994;
  (*elementAverageMassTable)["N"] = 14.0067;
  (*elementAverageMassTable)["+"] = 1.00739;  // proton

  Analyte::groupAverageMassTable = new map<string,double>();
  
  (*groupAverageMassTable)["H2O"] = 18.01524;
  (*groupAverageMassTable)["NH3"] = 17.0305; 
 	
  Analyte::elementMonoisotopicMassTable = new map<string, double>();
  
  (*elementMonoisotopicMassTable)["C"] = 12.00000;
  (*elementMonoisotopicMassTable)["H"] = 1.007825035;
  (*elementMonoisotopicMassTable)["O"] = 15.99491463;
  (*elementMonoisotopicMassTable)["N"] = 14.003074;
  (*elementMonoisotopicMassTable)["+"] = 1.00728;  // proton
  
  Analyte::groupMonoisotopicMassTable = new map<string,double>();
  
  (*groupMonoisotopicMassTable)["H2O"] = 18.01056;
  (*groupMonoisotopicMassTable)["NH3"] = 17.026549; 
   
}

// deleteTables - frees the memory associated with the tables. 
// (Best to call at the end of any program using the Peptide class)
void Analyte::deleteTables() {
	
  if (elementAverageMassTable) {
    delete (elementAverageMassTable);
  }
  elementAverageMassTable = NULL;

  if (groupAverageMassTable) {
    delete (groupAverageMassTable);
  }
  groupAverageMassTable = NULL;
  
  if (elementMonoisotopicMassTable) {
    delete (elementMonoisotopicMassTable);
  }
  elementMonoisotopicMassTable = NULL;

  if (groupMonoisotopicMassTable) {
    delete (groupMonoisotopicMassTable);
  }
  groupMonoisotopicMassTable = NULL;
}

FragmentIon::FragmentIon(double mz, string& annotation, unsigned int curIsotope) :
  m_ion(""),
  m_mz(mz),
  m_pos(0),
  m_charge(0),
  m_loss(0),
  m_prominence(0),
  m_mzDiff(0.0),
  m_isotope(0),
  m_bracket(false),
  m_assigned(false) {
  
  if (annotation.empty() || annotation[0] == '?') {
    return;
  }
  
  string::size_type pos = 0;
  
  if (annotation[0] == '[') {
    m_bracket = true;
  }
  
  m_ion = nextToken(annotation, 0, pos, "/], \t\r\n", "[");
  m_mzDiff = atof(nextToken(annotation, pos, pos, ", \t\r\n").c_str());
  
  if (m_ion[0] == 'I') {
    m_charge = 1; // only charge 1?
    m_isotope = 0; // no higher isotope?
    return;
  }
  
  char ionType = m_ion[0];
  m_pos = 0;
  
  if ((ionType >= 'a' && ionType <= 'c') ||  (ionType >= 'x' && ionType <= 'z')) {
    m_pos = atoi(nextToken(m_ion, 1, pos, "^i-+/]*, \t\r\n", " \t\r\n").c_str());
  }
  m_loss = 0;
  m_charge = 1;
    
  if (m_ion.find('i') == string::npos) {
    m_isotope = 0;
  } else {
    m_isotope = curIsotope + 1;
  }
    
  if (pos < m_ion.length() && (m_ion[pos] == '-' || m_ion[pos] == '+')) {
    // neutral loss
    m_loss = atoi(nextToken(m_ion, pos + 1, pos, "^i/]*, \t\r\n", " \t\r\n").c_str());
      
    if (m_ion[pos] == '+') m_loss *= -1; 
  }
    
  bool isNISTisotope = false;
  if (pos < m_ion.length() && (m_ion[pos] == 'i' || m_ion[pos] == '*')) { // NIST format, i before charge
    pos++;
  } 
    
  if (pos < m_ion.length() && (m_ion[pos] == '^')) {
    m_charge = atoi(nextToken(m_ion, pos + 1, pos, "i/]*, \t\r\n", " \t\r\n").c_str());
  }
    
  if (pos < m_ion.length() && (m_ion[pos] == 'i')) { // SpectraST format, i after charge
    pos++;
  }
   
}

// FragmentIon - creates a FragmentIon struct from the passed in information
FragmentIon::FragmentIon(string ionType, int pos, int loss, double mz, int ch, unsigned int prominence, unsigned int isotope, double mzDiff, bool bracket) :
  m_ion(""),
  m_pos(pos),
  m_loss(loss),
  m_mz(mz),
  m_charge(ch),
  m_prominence(prominence),
  m_isotope(isotope),
  m_mzDiff(mzDiff),
  m_bracket(bracket),
  m_assigned(false) {
  
  stringstream ss;
  ss << ionType;
  
  if (pos > 0) {
    ss << pos;
  }
  
  if (loss > 0) {
    ss << '-' << loss;
  } else if (loss < 0) {
    ss << '+' << -loss;
  }

  if (ch != 1) {
    ss << '^' << ch;
  }

  m_ion = ss.str();
  
}

string FragmentIon::getAnnotation() {
  stringstream ss;
  ss << m_ion << '/';
  ss.precision(2);
  ss << fixed << m_mzDiff;
  
  if (m_bracket) return ("[" + ss.str() + "]");
  
  return (ss.str());
}

// sortFragmentIonsByProminence - comparison function used by sort() to sort fragment ions by prominence
bool FragmentIon::sortFragmentIonPtrsByProminence(FragmentIon* a, FragmentIon* b) {

  if (a->m_prominence > b->m_prominence) {
    return (true);
  } else if (a->m_prominence < b->m_prominence) {
    return (false);
  } else {
  
    // in case of a tie, annotate with ion of smaller charge first
    if (a->m_charge < b->m_charge) {
      return (true);
    } else if (a->m_charge > b->m_charge) {
      return (false);
    } else {
      // if still tied, sort by ion string
      // this should have the effect that smaller losses will go before larger losses
      return (a->m_ion < b->m_ion);
    }
    
  }

}

