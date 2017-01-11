#ifndef GLYCAN_HPP_
#define GLYCAN_HPP_

#include "Analyte.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;


class Monosaccharide {

public:
  Monosaccharide(char symbol, string comp);
  ~Monosaccharide();
  double getResidueMonoisotopicMass() { return m_residueMonoisotopicMass; }
  double getResidueAverageMass() { return m_residueAverageMass; }
  string getCompCode();

private:
  char m_symbol;
  string m_compCode;
  double m_monoMass;
  double m_averageMass;
  double m_residueMonoisotopicMass;
  double m_residueAverageMass;
  map<char,int>* m_compTable;
  bool parseComp(const string& comp);
  double calcMonoMass();
  double calcAverageMass();

};

class GlycanFragment {
  
public:
  
  GlycanFragment(unsigned int rid);
  ~GlycanFragment();

  unsigned int m_rid;
  vector<pair<unsigned int,unsigned int>* >* m_brokenBonds;
  map<unsigned int, vector<unsigned int>* >* m_fwt;// forward table to store composition
  string smileName;
  void dere();
  bool exists(unsigned int id);
  void printFWT();

private:
  void init();

};

class Glycan : public Analyte {
  
public:
  
  Glycan();
  Glycan(Glycan& g);
  Glycan(string name,string format = "MassAnalyzer");
  ~Glycan();
  
  Glycan& operator=(Glycan& s);
  
  string m_name;
  unsigned int m_uid;
  
  static void defaultTables();
  static void deleteTables();

  static vector<Monosaccharide*>* monosaccharideList;
  static map<char,unsigned int>* monosaccharideSimpleIndex;
  static vector<Glycan*>* glycanList;
  static map<string,unsigned int>* glycanNameIndex;
  static bool sortLenAsc(const string& s1, const string& s2);
  static double getMonosaccharideResidueMonoisotopicMass(char monosaccharide);
  static double getMonosaccharideResidueAverageMass(char monosaccharide);
  static bool isKnownGlycan(string glycan);
  static double getGlycanMonoisotopicMass(string glycan);
  static double getGlycanAverageMass(string glycan);
  
  // This is a pre-made glycan objec, owned by the glycan class (static), and shared by any glycopeptide with this glycan.
  // It is declared const and cannot be changed by the caller. But the caller has to be careful NOT to delete it!
  static const Glycan* getGlycanPtr(string glycan); 

  string getName() const;
  int countMonosaccharides(string monosaccharide) const;
  
  double averageNeutralM() const { return m_averageMass; }
  double monoisotopicNeutralM() const { return m_monoisotopicMass; }

  map<string,double >* m_yfs;// y fragments
  map<string,double >* m_bfs;// b fragments

  string getComposition() const;
  void print() const;
  void getGlycosidicFragmentIons(vector<FragmentIon*>& ions, unsigned int prominence, double glycopeptideMass, int glycopeptideCharge) const;

private:
  double m_monoisotopicMass;
  double m_averageMass;
  map<unsigned int, vector<unsigned int>* >* m_fwt;// forward table to store composition
  map<unsigned int, char>* m_idx; // index for all ids (id=>symbol of monosaccharide)
  map<unsigned int, unsigned int>* m_pidx; // index for finding parent node
  //map<string, Linkage*> lidx; // index for linkages, for future
  
  

  void init();

  void getMAComp(const string& name, map<string,unsigned int>* maComp);
  void parseMAComp(map<string,unsigned int>* maComp);
  void initMAComp(map<string, unsigned int>* maComp);
  void updateMAComp(const string& optstr, const string& valstr, map<string,unsigned int>* maComp);

  unsigned int getNewID();
  void buildCore();
  void generateFragments();
  
  unsigned int addMonosaccharide(char monosaccharide,unsigned int pid);
  void deleteMonosaccharide(unsigned int id);
  void deleteMonosaccharideAsParent(unsigned int id);
  void deleteMonosaccharideAsChild(unsigned int id);
  void moveMonosaccharideToB(unsigned int id, GlycanFragment* bf, 
			     map<unsigned int, vector<unsigned int>* >* fwt);

  unsigned int chooseBranch(unsigned int mode=1);
  string fwt2smile(GlycanFragment* gf, unsigned int id);

  void breakBond(unsigned int p, unsigned int c, GlycanFragment* yf, GlycanFragment* bf, GlycanFragment* pf =NULL);
  void uniqFragments(vector<string>* gfs, map<string,double>* ufs);
  double calcGlycanFragmentMonoMass(string smileName);

  void copyFWT(map<unsigned int, vector<unsigned int>* >* ofwt,map<unsigned int,vector<unsigned int>* >* tfwt);
  double calcMonoisotopicMass();
  double calcAverageMass();
};


#endif /* GLYCAN_HPP_ */
