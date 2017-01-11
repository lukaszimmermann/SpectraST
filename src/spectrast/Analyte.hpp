#ifndef ANALYTE_HPP_
#define ANALYTE_HPP_

#include <string>
#include <vector>
#include <map>
#include <set>

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

using namespace std;

class FragmentIon {
  
public:  
  

  FragmentIon(double mz, string& annotation, unsigned int curIsotope = 0);  
  FragmentIon(string ionType, int position, int loss, double mz, int ch, unsigned int prom, unsigned int isotope = 0, double mzDiff = 0.0, bool bracket = false);

  string m_ion;
  int m_pos;
  double m_mz;
  int m_charge;
  unsigned int m_isotope;
  unsigned int m_prominence;
  double m_mzDiff;
  bool m_bracket;
  int m_loss;
  bool m_assigned;
  
  string getAnnotation();
  
  static bool sortFragmentIonPtrsByProminence(FragmentIon* a, FragmentIon* b); 
};

class Analyte {
	
public:

  // constructors, destructors and assignment operator
  Analyte();
  virtual ~Analyte();

  virtual string getName() const = 0;

  // methods to calculate masses and mass-to-charge ratios
  virtual double averageNeutralM() const = 0;
  virtual double monoisotopicNeutralM() const = 0;
  
  // fields

  // static methods to manage and access the mass tables
  static void defaultTables();
  static void deleteTables();  

  // static members - pointers to the tables of amino acids, modifications and their masses, etc - see below	
  static map<string, double>* elementAverageMassTable;
  static map<string, double>* elementMonoisotopicMassTable;
  static map<string, double>* groupAverageMassTable; 
  static map<string, double>* groupMonoisotopicMassTable;
  	
protected:
   
  static double calcApproximateAverageMass(double monoisotopicMass);
	

	
};



#endif /*ANALYTE_HPP_*/
