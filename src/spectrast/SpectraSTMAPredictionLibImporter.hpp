#ifndef SPECTRAST_MAPREDICTIONLIBIMPORTER_HPP_
#define SPECTRAST_MAPREDICTIONLIBIMPORTER_HPP_

#include "SpectraSTLibImporter.hpp"

class SpectraSTMAPredictionLibImporter: public SpectraSTLibImporter {
  
public:
  
  SpectraSTMAPredictionLibImporter(vector<string>& impFileNames, SpectraSTLib* lib, SpectraSTCreateParams& params);
  virtual ~SpectraSTMAPredictionLibImporter();
  virtual void import();

private:
  void readFromFile(string& impFileName);

  void makeLibEntry(string peptide,string description, int charge, string fragmentType, int collisionEnergy, 
				  double reactionTime, double isolationWidth, string instrumentModel, double resolution, 
				  double scanLowerBound, double scanUpperBound, int peakNum, vector< pair<double,double> >& peaks);
  
  bool makeMspModStr(vector<string>& mods, string& mspModStr);
  
};

#endif /*SPECTRASTMASSANALYZERPRECITIMPORTER_HPP_*/



