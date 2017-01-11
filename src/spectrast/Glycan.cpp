#include "Glycan.hpp"
#include <boost/concept_check.hpp>

vector<Monosaccharide*>* Glycan::monosaccharideList = NULL;
map<char,unsigned int>* Glycan::monosaccharideSimpleIndex = NULL;
vector<Glycan*>* Glycan::glycanList = NULL;
map<string,unsigned int>* Glycan::glycanNameIndex = NULL;


Monosaccharide::Monosaccharide(char symbol, string comp){
  m_symbol = symbol;
  m_compTable = new map<char,int>();
  //cout << "get symbol" << symbol <<"\t" << comp<< endl;
  if (parseComp(comp)){
    //cout << "Pass Parsecomp:";
    m_monoMass = calcMonoMass();
    m_averageMass = calcAverageMass();
    m_residueMonoisotopicMass = m_monoMass - 18.01056;//minus water
    m_residueAverageMass = m_averageMass - 18.01524;
    //cout << getMonoMass() << endl;
  }

  if ((symbol == 'G')||(symbol == 'M')||(symbol == 'L')){
    m_compCode = "HEX";
  }else if ( symbol == 'F' ){
    m_compCode = "dHEX";
  }else if ( symbol == 'S' ){
    m_compCode = "NANA";
  }else if ( symbol == 'Y' ){
    m_compCode = "HEXNAc";
  }else if ( symbol == 'Z' ){
    m_compCode = "NGNA";
  }else {
    m_compCode = "UNK";
  }
}

Monosaccharide::~Monosaccharide(){
  if (m_compTable) delete m_compTable;
}

string Monosaccharide::getCompCode(){
  return m_compCode;
}

bool Monosaccharide::parseComp(const string& comp){
  bool success = true;
  if (m_compTable){
    delete m_compTable;
    m_compTable = new map<char,int>();
  }
  stringstream opt, val;
  for ( int i = 0; i < comp.size(); i ++){
    if ((comp[i]>='A')&&(comp[i]<='Z')){
      if (i != 0){
	if (val.str().size() == 0){
	  val << "1";
	}
	string optstr = opt.str();
	(*m_compTable)[(optstr[0])] = atoi(val.str().c_str());
	opt.str("");
	val.str("");
      }
      opt << comp[i];
    }else if ((comp[i]>='0')&&(comp[i]<='9')){
      val << comp[i];
    }
  }
  if (val.str().size() == 0){
    val << "1";
  }
  string optstr = opt.str();
  (*m_compTable)[optstr[0]] = atoi(val.str().c_str());
  return success;
}

double Monosaccharide::calcMonoMass(){
  map<char, int>::iterator it;
  double sum = 0.0;
  for (it = m_compTable->begin(); it != m_compTable->end(); it++) {
    string elem("");
    elem += it->first;
    int n = it->second;
     
    map<string, double>::iterator found = Analyte::elementMonoisotopicMassTable->find(elem);
    if (found != Analyte::elementMonoisotopicMassTable->end()) {
      sum += found->second * n;
    }
  }
  return sum;
}

double Monosaccharide::calcAverageMass(){
  map<char,int>::iterator it;
  double sum = 0.0;
  for ( it = m_compTable->begin(); it != m_compTable->end(); it++) {
    string elem("");
    elem += it->first;
    int n = it->second;
    
    map<string, double>::iterator found = Analyte::elementAverageMassTable->find(elem);
    if (found != Analyte::elementAverageMassTable->end()) {
      sum += found->second * n;
    }
  }
  return sum;
}



GlycanFragment::GlycanFragment(unsigned int rid){
  m_rid = rid;
  init();
}

GlycanFragment::~GlycanFragment(){
  if (m_fwt){
    map<unsigned int, vector<unsigned int>* >::iterator ifwt;
    for ( ifwt = m_fwt->begin(); ifwt != m_fwt->end(); ifwt ++){
      if (ifwt->second){
	delete (ifwt->second);
      }
    }
    delete m_fwt;
  }
  
  if (m_brokenBonds){
    vector<pair<unsigned int, unsigned int>* >::iterator ib;
    for ( ib = m_brokenBonds->begin(); ib != m_brokenBonds->end(); ib ++){
      delete (*ib);
    }
    delete m_brokenBonds;
  }

}

void GlycanFragment::printFWT(){
  cout << "fragment rid:" << m_rid << endl;
  cout << "fragment name:" << smileName << endl;
  map<unsigned int, vector<unsigned int>* >::iterator ifwt;
  for (ifwt = m_fwt->begin(); ifwt != m_fwt->end(); ifwt++){
    vector<unsigned int>::iterator ic;
    cout << ifwt->first << ":" ;
    for ( ic = ifwt->second->begin(); ic != ifwt->second->end(); ic ++){
      cout << *ic << "," ;
    }
    cout << endl;
  }
}

bool GlycanFragment::exists(unsigned int id){
  if (m_fwt){
    map<unsigned int, vector<unsigned int>* >::iterator ifwt = m_fwt->find(id);
    if (ifwt != m_fwt->end()){
      return true;
    }
  }
  return false;
}

void GlycanFragment::init(){
  m_fwt = new map<unsigned int, vector<unsigned int>* >();
  m_brokenBonds = new vector<pair<unsigned int, unsigned int>* >();
}

void Glycan::copyFWT(map<unsigned int, vector<unsigned int>* >* ofwt, map<unsigned int, vector<unsigned int>* >* tfwt){
  map<unsigned int, vector<unsigned int>* >::iterator io;
  for ( io = ofwt->begin(); io != ofwt->end(); io++){
    unsigned int pid = io->first;
    (*tfwt)[pid] = new vector<unsigned int>();
    map<unsigned int, vector<unsigned int>* >::iterator it = tfwt->find(pid);
    vector<unsigned int>* vpc = io->second;
    vector<unsigned int>::iterator ic;
    for (ic = vpc->begin(); ic != vpc->end(); ic ++){
      it->second->push_back(*ic);
    }
  }
}

void GlycanFragment::dere(){
  vector<unsigned int> keys;
  map<unsigned int, vector<unsigned int>* >::iterator ifwt;
  for ( ifwt = m_fwt->begin(); ifwt != m_fwt->end(); ifwt ++){
    keys.push_back(ifwt->first);
  }
  if (keys.size() > 1){
  
    vector<unsigned int>::iterator ik;
    for ( ik = keys.begin(); ik != keys.end(); ik++){
      ifwt = m_fwt->find(*ik);
      if (ifwt->second == NULL){
	m_fwt->erase(ifwt);
      }else if (ifwt->second->size() == 0){
	m_fwt->erase(ifwt);
      }
    }
  }
}


Glycan::Glycan() {
  
}

Glycan::~Glycan() {
  if (m_fwt) delete m_fwt;
  if (m_idx) delete m_idx;
  if (m_pidx) delete m_pidx;
  if (m_yfs) delete m_yfs;
  if (m_bfs) delete m_bfs;
}

string Glycan::getName() const {
  return (m_name);
}

Glycan::Glycan(string name,string format){
  if (format == "MassAnalyzer"){
    m_name = name;
    //cout << name << endl;
    //cout << "getMAComp" << endl;
    map<string,unsigned int>* maComp = new map<string, unsigned int>();
    getMAComp(name,maComp);
    if (0){//for debug
      map<string,unsigned int>::iterator icp;
      for (icp = maComp->begin(); icp != maComp->end(); icp ++){
	cout << icp->first << "\t" << icp->second << endl;
      }
    }
    //cout << "parseMAComp" << endl;
    parseMAComp(maComp);
    //cout << "generateFragment" << endl;
    generateFragments();
    delete maComp;
    m_monoisotopicMass = calcMonoisotopicMass();
    m_averageMass = calcAverageMass();
  }
}


void Glycan::initMAComp(map<string, unsigned int>* maComp){
  (*maComp)["A"] = 0;
  (*maComp)["F"] = 0;
  (*maComp)["B"] = 0;
  (*maComp)["Gn"] = 0;
  (*maComp)["S"] = 0;
  (*maComp)["Sg"] = 0;
  (*maComp)["G"] = 0;
  (*maComp)["M"] = 0;
}

int Glycan::countMonosaccharides(string monosaccharide) const {
  
  char m;
  int sum = 0;
  if (monosaccharide == "Gn") {
    m = 'Y';
  } else if (monosaccharide == "F") {
    m = 'F';
  } else if (monosaccharide == "M") {
    m = 'M';
  } else if (monosaccharide == "Ga") {
    m = 'L';
  } else if (monosaccharide == "S") {
    m = 'S';
  } else if (monosaccharide == "Sg") {
    m = 'Z';
  }
  map<unsigned int, char>::iterator i;
  for (i = m_idx->begin(); i != m_idx->end(); i++) {
    if (i->second == m) {
      sum++;
    }
  }
  return sum;
}

void Glycan::getMAComp(const string& name, map<string, unsigned int>* maComp) {
  initMAComp(maComp);
  if (name == "Gn"){
    (*maComp)["Gn"] = 1;
  }else if (name == "GnF"){
    (*maComp)["Gn"] = 1;
    (*maComp)["F"] = 1;
  }else{
    stringstream opt,val;
    string optstr,valstr;
    for ( int i = 0; i < name.size(); i ++){
      if ((name[i]>='a')&&(name[i]<='z')){
	opt << name[i];
      }else if ((name[i]>='A')&&(name[i] <='Z')){
	if (i != 0){
	  optstr = opt.str();
	  valstr = val.str();
	  updateMAComp(optstr,valstr,maComp);
	  opt.str("");
	  val.str("");
	}
	opt << name[i];
      }else if ((name[i]>='0')&&(name[i] <='9')){
	val << name[i];
      }
    }
    optstr = opt.str();
    valstr = val.str();
    updateMAComp(optstr,valstr,maComp);
    opt.str("");
    val.str("");
  }
}

void Glycan::updateMAComp(const string& optstr,const string& valstr, map<string,unsigned int>* maComp) {
  unsigned int valint = 1;
  if (valstr.size() > 0){
    valint = (unsigned int)atoi(valstr.c_str());
  }
  map<string,unsigned int>::iterator found = maComp->find(optstr);
  if (found != maComp->end()){
    found->second = valint;
  }
}

void Glycan::buildCore() {
  (*m_fwt)[1] = new vector<unsigned int>();
  (*m_fwt)[1]->push_back(2);
  (*m_fwt)[2] = new vector<unsigned int>();
  (*m_fwt)[2]->push_back(3);
  (*m_fwt)[3] = new vector<unsigned int>();
  (*m_fwt)[3]->push_back(4);
  (*m_fwt)[3]->push_back(5);
  
  (*m_idx)[1] = 'Y';
  (*m_idx)[2] = 'Y';
  (*m_idx)[3] = 'M';
  (*m_idx)[4] = 'M';
  (*m_idx)[5] = 'M';
  
  map<unsigned int, vector<unsigned int>* >::iterator it;
  for ( it = m_fwt->begin(); it!= m_fwt->end(); it++) {
    unsigned int pid = it->first;
    vector<unsigned int>* vpc = it->second;
    for (int i = 0; i < vpc->size(); i++) {
      (*m_pidx)[(*vpc)[i]] = pid;
    }
  }
}

void Glycan::parseMAComp(map<string, unsigned int>* maComp) {
  //cout << "init" << endl;
  init();
  
  //cout << "buildCore" << endl;
  buildCore();

  //cout << "buildBranches" << endl;
  if ((*maComp)["A"] > 0){
    // complex or hybrid
    if (((*maComp)["M"] > 0)||((*maComp)["A"] < 2)) {
      //hybrid
      //cout << "GlycanType:Hybrid" << endl;
      if ((*maComp)["Sg"] > 0){
	for (int i = 0; i < (*maComp)["Sg"]; i++) {
	  unsigned int pid = chooseBranch();
	  pid = addMonosaccharide('Y',pid);
	  pid = addMonosaccharide('L',pid);
	  pid = addMonosaccharide('Z',pid);
	}
      }
      if ((*maComp)["S"] > 0){
	for (int i = 0; i < (*maComp)["S"]; i++) {
	  unsigned int pid = chooseBranch();
	  pid = addMonosaccharide('Y',pid);
	  pid = addMonosaccharide('L',pid);
	  pid = addMonosaccharide('S',pid);
	}
      }
      if ((*maComp)["G"] > 0){
	for (int i = 0; i < (*maComp)["G"]; i++) {
	  unsigned int pid = chooseBranch();
	  pid = addMonosaccharide('Y',pid);
	  pid = addMonosaccharide('L',pid);
	}
      }

      int numGn = (*maComp)["A"]-(*maComp)["S"]-(*maComp)["Sg"]-(*maComp)["G"];
      if (numGn > 0) {
	for ( int i = 0; i < numGn; i++) {
	  unsigned int pid = 4;//chooseBranch(1);
	  pid = addMonosaccharide('Y', pid);
	}
      }
      
      if ((*maComp)["M"] == 2) {
	unsigned int id = 5;
	deleteMonosaccharide(id);
      } else if ((*maComp)["M"]>3) {
	for (int i = 3; i < (*maComp)["M"]; i++) {
	  unsigned int pid = 5;
	  addMonosaccharide('M',pid);
	}
      }
    } else {
      //complex
      //cout << "GlycanType:Complex" << endl;
      if ((*maComp)["Sg"] > 0) {
	for (int i = 0; i < (*maComp)["Sg"]; i++) {
	  unsigned int pid = chooseBranch();
	  pid = addMonosaccharide('Y', pid);
	  pid = addMonosaccharide('L', pid);
	  pid = addMonosaccharide('Z', pid);
	}
      }
      if ((*maComp)["S"] > 0){
	for (int i = 0; i < (*maComp)["S"]; i++) {
	  unsigned int pid = chooseBranch();
	  pid = addMonosaccharide('Y', pid);
	  pid = addMonosaccharide('L', pid);
	  pid = addMonosaccharide('S', pid);
	}
      }
      if ((*maComp)["G"]>0){
	for (int i = 0; i < (*maComp)["G"]; i++) {
	  unsigned int pid = chooseBranch();
	  pid = addMonosaccharide('Y', pid);
	  pid = addMonosaccharide('L', pid);
	}
      }
      
      int numGn = (*maComp)["A"] - (*maComp)["S"] - (*maComp)["Sg"] - (*maComp)["G"];
      if (numGn > 0) {
	for (int i = 0; i < numGn; i++) {
	  unsigned int pid = chooseBranch(2);
	  pid = addMonosaccharide('Y', pid);
	}
      }
    }
  } else {
    if ((*maComp)["M"] > 3) {
      //cout << "GlycanType:HighMan" << endl;
      for (int i = 4; i <= (*maComp)["M"]; i++) {
	if (i <= 6) {
	  unsigned int pid = chooseBranch();
	  pid = addMonosaccharide('M', pid);
	} else {
	  unsigned int pid = i - 1;
	  pid = addMonosaccharide('M', pid);
	}
      }
    } else {
      //cout << "GlycanType:Core" << endl;
      if ((*maComp)["Gn"] == 1) {
	unsigned int id = 2;
	deleteMonosaccharide(id);
      } else if ((*maComp)["M"] == 0) {
	unsigned int id = 3;
	deleteMonosaccharide(id);
      } else if ((*maComp)["M"] == 1) {
	unsigned int id = 4;
	deleteMonosaccharide(id);
	id = 5;
	deleteMonosaccharide(id);
      } else if ((*maComp)["M"] == 2) {
	unsigned int id = 5;
	deleteMonosaccharide(id);
      } else if ((*maComp)["M"] == 3) {
	// do nothing
      }
    }
  }

  //cout << "buildF" << endl;
  
  if ((*maComp)["F"] == 1) {
    unsigned int pid = 1;
    addMonosaccharide('F', pid);
  }
  //cout << "buildB" << endl;
  if ((*maComp)["B"] == 1) {
    unsigned int pid = 3;
    addMonosaccharide('Y', pid);
  }

}

unsigned int Glycan::chooseBranch(unsigned int i) {

  int atn4Num = 0;
  int atn5Num = 0;
  unsigned int id = 0;
  map<unsigned int, vector<unsigned int>* >::iterator ifwt = m_fwt->find(4);
  if (ifwt != m_fwt->end()) {
    if (ifwt->second){
	atn4Num = ifwt->second->size();
    }
  }
  
  ifwt = m_fwt->find(5);
  if (ifwt != m_fwt->end()) {
    if (ifwt->second) {
      atn5Num = ifwt->second->size();
    }
  }
  
  if (i == 1) {
    if (atn4Num >= 2) {
      id = 5;
    } else {
      id = 4;
    }
  } else if (i == 2) {
    if (atn4Num > atn5Num) {
      id = 5;
    } else {
      id = 4;
    }
  }
  return (id);
}

unsigned int Glycan::addMonosaccharide(char monosaccharide, unsigned int pid) {
  unsigned int id = getNewID();
  map<unsigned int, vector<unsigned int>* >::iterator ifwt  = m_fwt->find(pid);
  if (ifwt != m_fwt->end()) {
    vector<unsigned int>* vpc = ifwt->second;
    vpc->push_back(id);
  } else {
    (*m_fwt)[pid] = new vector<unsigned int>();
    (*m_fwt)[pid]->push_back(id);
  }
  (*m_idx)[id] = monosaccharide;
  (*m_pidx)[id] = pid;
  return (id);
}

void Glycan::deleteMonosaccharide(unsigned int id) {
  deleteMonosaccharideAsParent(id);
  deleteMonosaccharideAsChild(id);
}

void Glycan::deleteMonosaccharideAsChild(unsigned int id) {
  map<unsigned int, unsigned int>::iterator ip = m_pidx->find(id);
  unsigned int pid = 0;
  if (ip != m_pidx->end()) {
    pid = ip->second;
  } else {
    return;
  }
  map<unsigned int, vector<unsigned int>* >::iterator ifwt = m_fwt->find(pid);
  int i = 0;
  for (i = 0; i < ifwt->second->size(); i++) {
    if (ifwt->second->at(i) == id) {
      break;
    }
  }
  ifwt->second->erase(ifwt->second->begin() + i);
  m_pidx->erase(id);
}

void Glycan::deleteMonosaccharideAsParent(unsigned int id) {
  //cout << "delete id:" << id << endl;
  
  map<unsigned int, vector<unsigned int>* >::iterator ifwt = m_fwt->find(id);
  if (0) {//for debug
    for (ifwt = m_fwt->begin(); ifwt!= m_fwt->end(); ifwt++) {
      cout << ifwt->first << ":" ;
      vector<unsigned int>::iterator i1;
      for (i1 = ifwt->second->begin(); i1 != ifwt->second->end(); i1++) {
	cout << *i1 << ",";
      }
      cout << endl;
    }
    ifwt = m_fwt->find(id);
  }
  // delete all children of id
  if (ifwt != m_fwt->end()) {
    vector<unsigned int>* vpc = ifwt->second;
    vector<unsigned int>::iterator ic;
    //cout << "test child" << endl;
    for (ic = vpc->begin(); ic != vpc->end(); ic ++) {
      //cout << "child:" << *ic << endl;
      //cout << "child size:" << vpc->size()<< endl;
      deleteMonosaccharideAsParent(*ic);
      // delete its record in parent index
      m_pidx->erase(*ic);
    }
    delete vpc;
    // delete its column
    m_fwt->erase(id);
  }
  // delete id from its parent column
  /*
  map<unsigned int, unsigned int>::iterator ip = m_pidx->find(id);
  unsigned int pid;
  if (ip != m_pidx->end()){
    pid = ip->second;
  }else{
    cout << "ERROR! can't find parent id of " << id << endl;
    exit(1);
  }
  ifwt = m_fwt->find(pid);
  if (ifwt != m_fwt->end()){
    vector<unsigned int>* vpc = ifwt->second;
    int i = 0;
    for ( i = 0; i < vpc->size(); i++){
      if ((*vpc)[i] == id){
	break;
      }
    }
    vpc->erase(vpc->begin()+i);
  }
  */
  
  // delete its record in index
  m_idx->erase(id);
}


unsigned int Glycan::getNewID() {
  int nodeNum = m_idx->size();
  unsigned int id = (unsigned int)(nodeNum + 1);
  map<unsigned int,char>::iterator found = m_idx->find(id);
  while (found != m_idx->end()) {
    id++;
    found = m_idx->find(id);
  }
  return (id);
}

void Glycan::generateFragments() {
  
  vector<pair<unsigned int, unsigned int> > bonds;
  map<unsigned int, vector<unsigned int>* >::iterator ifwt;
  for (ifwt = m_fwt->begin(); ifwt != m_fwt->end(); ifwt++) {
    unsigned int pid = ifwt->first;
    vector<unsigned int>* vpc = ifwt->second;
    vector<unsigned int>::iterator ic;
    for (ic = vpc->begin(); ic != vpc->end(); ic++) {
      bonds.push_back(make_pair(pid, *ic));
    }
  }
  if(0) {
    vector< pair<unsigned int, unsigned int> >::iterator ib;
    for (ib = bonds.begin(); ib != bonds.end(); ib++) {
      cout << (*ib).first << "," << (*ib).second << endl;
    }
  }

  vector<string>* yfs =  new vector<string>();
  vector<string>* bfs = new vector<string>();

  for ( int i = 0; i < bonds.size(); i++) {
    unsigned int p1 = bonds[i].first;
    unsigned int c1 = bonds[i].second;
    //cout << "LEVEL1:" << p1 << "," << c1 << endl;
    GlycanFragment* yf = new GlycanFragment(1);
    GlycanFragment* bf = new GlycanFragment(c1);
    breakBond(p1, c1, yf, bf);
    yf->smileName = fwt2smile(yf, yf->m_rid);
    bf->smileName = fwt2smile(bf, bf->m_rid);
    yfs->push_back(yf->smileName);
    bfs->push_back(bf->smileName);
    if (0) {
      cout << "y ions:" << endl;
      yf->printFWT();
      cout << "b ions:" << endl;
      bf->printFWT();
    }
    for (int j = i + 1; j < bonds.size(); j++) {
      unsigned int p2 = bonds[j].first;
      unsigned int c2 = bonds[j].second;
      //cout << "LEVEL2:" << p2 << "," << c2 << endl;
      if (yf->exists(p2)) {
	//cout << "create yf1, yf2" << endl;
	// break y
	GlycanFragment* yf1 = new GlycanFragment(1);
	GlycanFragment* yf2 = new GlycanFragment(c2);
	breakBond(p2, c2, yf1, yf2, yf);
	yf1->smileName = fwt2smile(yf1, yf1->m_rid);
	yf2->smileName = fwt2smile(yf2, yf2->m_rid);
	yfs->push_back(yf1->smileName);
	bfs->push_back(yf2->smileName);
	if (0) {
	  yf1->printFWT();
	  yf2->printFWT();
	}
	delete yf1;
	delete yf2;
      } else if (bf->exists(p2)) {
	//cout << "create bf1,bf2" << endl;
	// break b
	GlycanFragment* bf1 = new GlycanFragment(c1);
	GlycanFragment* bf2 = new GlycanFragment(c2);
	breakBond(p2, c2, bf1, bf2, bf);
	bf1->smileName = fwt2smile(bf1, bf1->m_rid);
	bf2->smileName = fwt2smile(bf2, bf2->m_rid);
	bfs->push_back(bf1->smileName);
	bfs->push_back(bf2->smileName);
	if (0) {
	  bf1->printFWT();
	  bf2->printFWT();
	}
	delete bf1;
	delete bf2;
      } else {
	cout << "ERROR! can't find bond " << p2 << "," << c2 << endl;
	exit(1);
      }      
    }
    delete yf;
    delete bf;
  }
  
  uniqFragments(bfs, m_bfs);
  uniqFragments(yfs, m_yfs);
  delete bfs;
  delete yfs;
  
}


void Glycan::uniqFragments(vector<string>* gfs, map<string, double>* ufs) {
  vector<string>::iterator it;
  for (it = gfs->begin(); it != gfs->end(); it++) {
    map<string, double>::iterator found = ufs->find(*it);
    if (found == ufs->end()) {
      (*ufs)[*it] = calcGlycanFragmentMonoMass(*it);
    }
  }
}

double Glycan::calcGlycanFragmentMonoMass(string smileName) {
  double sum = 0;
  for (int i = 0 ; i < smileName.size(); i++) {
    char m = smileName[i];
    if ((m >= 'A') && (m <= 'Z')) {
      map<char, unsigned int>::iterator im = monosaccharideSimpleIndex->find(m);
      if (im != monosaccharideSimpleIndex->end()) {
	unsigned int mid = im->second; 
	sum += (*monosaccharideList)[mid]->getResidueMonoisotopicMass();
      } else {
	cout << "ERROR! can't find monosaccharide :" << m << endl;
	exit(1);
      }
    }
  }
  return sum;
}

string Glycan::fwt2smile(GlycanFragment* gf, unsigned int id) {
  //cout << "fwt2smile:" << id << endl;
  stringstream gss;
  gss.str("");
  map<unsigned int , char>::iterator it = m_idx->find(id);
  if (it != m_idx->end()){
    gss << it->second ;
  }else{
    cout << "ERROR! can't find root id:" << id << " in glycan index table" << endl;
    exit(1);
  }

  map<unsigned int, vector<unsigned int>* >::iterator ifwt = gf->m_fwt->find(id);
  if (ifwt!=gf->m_fwt->end()){
    vector<unsigned int>* vpc = ifwt->second;
    if (vpc == NULL){
      return gss.str();
    }else if (vpc->size() == 0){
      return gss.str();
    }else if (vpc->size() == 1){
      unsigned int cid = (*vpc)[0];
      return gss.str() + fwt2smile(gf,cid);
    }else{
      vector<string> sc;
      for ( int i =0; i < vpc->size(); i++){
	sc.push_back(fwt2smile(gf,(*vpc)[i]));
      }
      sort(sc.begin(),sc.end(), Glycan::sortLenAsc);
      for ( int i = 0; i < sc.size() -1; i ++){
	gss << "(" << sc[i] << ")"; 
      }
      gss << sc[sc.size()-1];
      return gss.str();
    }
  }else{
    //cout << "ERROR! don't exist root in glycan fragment, root id: " << id << endl; 
    return gss.str();
  }
}

bool Glycan::sortLenAsc(const string& s1, const string& s2){
  int l1 = 0;
  int l2 = 0;
  for ( int i = 0; i < s1.size(); i ++){
    if ((s1[i]=='(')||(s1[i]==')')){
      continue;
    }
    l1 ++;
  }
  for ( int i = 0; i < s2.size(); i ++){
    if ((s2[i]=='(')||(s2[i]==')')){
      continue;
    }
    l2 ++;
  }
  return (l1<l2);
}

void Glycan::breakBond(unsigned int pid, unsigned int cid, GlycanFragment* yf, GlycanFragment* bf, GlycanFragment* pf){
  map<unsigned int,vector<unsigned int>* >* tmpFWT = new map<unsigned int, vector<unsigned int>* >();
  if (pf == NULL){
    copyFWT(m_fwt,tmpFWT);
  }else{
    copyFWT(pf->m_fwt,tmpFWT);
  }
  if(0){
    map<unsigned int, vector<unsigned int>* >::iterator i1 ;
    for (i1 = tmpFWT->begin(); i1 != tmpFWT->end(); i1 ++){
      cout << i1->first << ":" ;
      vector<unsigned int>::iterator i2;
      for ( i2 = i1->second->begin(); i2 != i1->second->end(); i2 ++){
	cout << (*i2) << "," ;
      }
      cout << endl;
    }
  }
  map<unsigned int, vector<unsigned int>* >::iterator ifwt = tmpFWT->find(pid);
  if (ifwt != tmpFWT->end()){
    moveMonosaccharideToB(cid,bf,tmpFWT);
    vector<unsigned int>* vpc = ifwt->second;
    int i = 0;
    for ( i = 0; i < vpc->size(); i ++){
      if ((*vpc)[i] == cid) break;
    }
    vpc->erase(vpc->begin() + i);
  }
  if (pf == NULL){
    yf->m_rid = 1;
    bf->m_rid = cid;
  }else{
    yf->m_rid = pf->m_rid;
    bf->m_rid = cid;
  }
  copyFWT(tmpFWT,yf->m_fwt);
  yf->dere();
  bf->dere();
  delete tmpFWT;// more delete in future
}

void Glycan::moveMonosaccharideToB(unsigned int id, GlycanFragment* bf, map<unsigned int, vector<unsigned int>* >* fwt){
  //cout << "move " << id << " To B" << endl;
  map<unsigned int, vector<unsigned int>* >::iterator ifwt = fwt->find(id);
  if (ifwt != fwt->end()){
    vector<unsigned int>* vpc = ifwt->second;
    for ( int i = 0; i < vpc->size(); i++){
      moveMonosaccharideToB((*vpc)[i],bf,fwt);
    }
    fwt->erase(id);
    (*(bf->m_fwt))[id]=vpc;
  }else{
    (*(bf->m_fwt))[id] = new vector<unsigned int>();
  }
  //cout << "moved" << endl;
}

void Glycan::init(){
  m_fwt = new map<unsigned int, vector<unsigned int>* >();
  m_idx = new map<unsigned int, char>();
  m_pidx = new map<unsigned int, unsigned int>();
  m_yfs = new map<string, double>();
  m_bfs = new map<string, double>();
}

double Glycan::getMonosaccharideResidueMonoisotopicMass(char monosaccharide) {
  vector<Monosaccharide*>::iterator im;
  map<char, unsigned int>::iterator isi = monosaccharideSimpleIndex->find(monosaccharide);
  unsigned int gid;
  if (isi != monosaccharideSimpleIndex->end()){
    gid = isi->second;
  } else {
    cout << "ERROR! can't find id of :" << monosaccharide << endl;
    exit(1);
  }
  //cout << gid << endl;
  return ((*monosaccharideList)[gid])->getResidueMonoisotopicMass();
}

double Glycan::getMonosaccharideResidueAverageMass(char monosaccharide) {
  vector<Monosaccharide*>::iterator im;
  map<char, unsigned int>::iterator isi = monosaccharideSimpleIndex->find(monosaccharide);
  unsigned int gid;
  if (isi != monosaccharideSimpleIndex->end()){
    gid = isi->second;
  }else{
    cout << "ERROR! can't find id of :" << monosaccharide << endl;
    exit(1);
  }
  //cout << gid << endl;
  return ((*monosaccharideList)[gid])->getResidueAverageMass();
}


bool Glycan::isKnownGlycan(string glycan) {
  
  if (glycanNameIndex->find(glycan) != glycanNameIndex->end()) {
    return (true);
  } else {
    return (false);
  }
}

double Glycan::getGlycanMonoisotopicMass(string glycan) {
  
  map<string, unsigned int>::const_iterator found = glycanNameIndex->find(glycan);
  if (found == glycanNameIndex->end()) {
    return (0.0);
  }
  unsigned int id = found->second;
  return ((*glycanList)[id]->monoisotopicNeutralM());
}

double Glycan::getGlycanAverageMass(string glycan) {
  map<string, unsigned int>::const_iterator found = glycanNameIndex->find(glycan);
  if (found == glycanNameIndex->end()) {
    return (0.0);
  }
  unsigned int id = found->second;
  return ((*glycanList)[id]->averageNeutralM());
}

double Glycan::calcMonoisotopicMass() {
  double sum = 0.0;
  if (m_fwt && m_idx && m_pidx) {
    map<unsigned int, char>::iterator id;
    for (id = m_idx->begin(); id != m_idx->end(); id++) {
      sum += getMonosaccharideResidueMonoisotopicMass(id->second);
    }
  }
  //cout << "mono=" << sum << endl;
  return sum;
}

double Glycan::calcAverageMass() {
  double sum = 0.0;
  if (m_fwt && m_idx && m_pidx) {
    map<unsigned int, char>::iterator id;
    for (id = m_idx->begin(); id != m_idx->end(); id ++) {
      sum += getMonosaccharideResidueAverageMass(id->second);
      //cout << getMonosaccharideResidueAverageMass << endl;
      //cout << id->second << endl;
    }
  }
  //cout << "ave=" << sum << endl;
  return sum;
}

const Glycan* Glycan::getGlycanPtr(string glycan) {
  
  map<string, unsigned int>::iterator git = glycanNameIndex->find(glycan);
  if (git != glycanNameIndex->end()) {
    unsigned int id = git->second;
    const Glycan* g = (const Glycan*)((*glycanList)[id]);
    return (g);
  } else {
    return (NULL);
  }
}

/*
Glycan* Glycan::getGlycanInIndexPtr(string glycan) {
  map<string,unsigned int>::iterator git = glycanNameIndex->find(glycan);
  if (git != glycanNameIndex->end()){
    unsigned int id = git->second;
    Glycan* g = (*glycanList)[id];
    return g;
  }else{
    cout << "ERROR! Glycan::getGlycanInIndexPtr: can't find glycan:" << glycan << endl;
    exit(1);
  }
}
*/

void Glycan::defaultTables(){
 
  // cout << "Glycan::monosaccharideList" << endl;

  Glycan::monosaccharideList = new vector<Monosaccharide*>();
  Glycan::monosaccharideSimpleIndex = new map<char,unsigned int>();

  Monosaccharide* G = new Monosaccharide('G',"C6H12O6");//Glc

  Monosaccharide* M = new Monosaccharide('M',"C6H12O6");//Man

  Monosaccharide* L = new Monosaccharide('L',"C6H12O6");//Gal

  Monosaccharide* F = new Monosaccharide('F',"C6H12O5");//Fuc

  Monosaccharide* S = new Monosaccharide('S',"C11H19NO9");//NANA

  Monosaccharide* Z = new Monosaccharide('Z',"C11H19NO10");//NGNA

  Monosaccharide* Y = new Monosaccharide('Y',"C8O6NH15");//GlcNAc

  monosaccharideList->push_back(G); (*monosaccharideSimpleIndex)['G'] = 0;
  monosaccharideList->push_back(M); (*monosaccharideSimpleIndex)['M'] = 1;
  monosaccharideList->push_back(L); (*monosaccharideSimpleIndex)['L'] = 2;
  monosaccharideList->push_back(F); (*monosaccharideSimpleIndex)['F'] = 3;
  monosaccharideList->push_back(S); (*monosaccharideSimpleIndex)['S'] = 4;
  monosaccharideList->push_back(Y); (*monosaccharideSimpleIndex)['Y'] = 5;
  monosaccharideList->push_back(Z); (*monosaccharideSimpleIndex)['Z'] = 6;

  //cout << "Glycan::GlycanList" << endl;
  Glycan::glycanList = new vector<Glycan*>();
  Glycan::glycanNameIndex = new map<string,unsigned int>();
  vector<string> names;
  names.push_back("A2G0");
  names.push_back("A2G0B");
  names.push_back("A2G0F");
  names.push_back("A2G0FB");
  names.push_back("A3G0");
  names.push_back("A3G0B");
  names.push_back("A3G0F");
  names.push_back("A3G0FB");
  names.push_back("A4G0");
  names.push_back("A4G0B");
  names.push_back("A4G0F");
  names.push_back("A4G0FB");
  names.push_back("A2G1");
  names.push_back("A2G1B");
  names.push_back("A2G1F");
  names.push_back("A2G1FB");
  names.push_back("A2G2");
  names.push_back("A2G2B");
  names.push_back("A2G2F");
  names.push_back("A2G2FB");
  names.push_back("A3G1");
  names.push_back("A3G1B");
  names.push_back("A3G1F");
  names.push_back("A3G1FB");
  names.push_back("A4G1");
  names.push_back("A4G1B");
  names.push_back("A4G1F");
  names.push_back("A4G1FB");
  names.push_back("A3G2");
  names.push_back("A3G2B");
  names.push_back("A3G2F");
  names.push_back("A3G2FB");
  names.push_back("A4G2");
  names.push_back("A4G2B");
  names.push_back("A4G2F");
  names.push_back("A4G2FB");
  names.push_back("A3G3");
  names.push_back("A3G3B");
  names.push_back("A3G3F");
  names.push_back("A3G3FB");
  names.push_back("A4G3");
  names.push_back("A4G3B");
  names.push_back("A4G3F");
  names.push_back("A4G3FB");
  names.push_back("A4G4");
  names.push_back("A4G4B");
  names.push_back("A4G4F");
  names.push_back("A4G4FB");
  names.push_back("A2S1G0");
  names.push_back("A2S1G0B");
  names.push_back("A2S1G0F");
  names.push_back("A2S1G0FB");
  names.push_back("A2S1G1");
  names.push_back("A2S1G1B");
  names.push_back("A2S1G1F");
  names.push_back("A2S1G1FB");
  names.push_back("A2S2");
  names.push_back("A2S2B");
  names.push_back("A2S2F");
  names.push_back("A2S2FB");
  names.push_back("A3S1G0");
  names.push_back("A3S1G0B");
  names.push_back("A3S1G0F");
  names.push_back("A3S1G0FB");
  names.push_back("A4S1G0");
  names.push_back("A4S1G0B");
  names.push_back("A4S1G0F");
  names.push_back("A4S1G0FB");
  names.push_back("A3S1G1");
  names.push_back("A3S1G1B");
  names.push_back("A3S1G1F");
  names.push_back("A3S1G1FB");
  names.push_back("A4S1G1");
  names.push_back("A4S1G1B");
  names.push_back("A4S1G1F");
  names.push_back("A4S1G1FB");
  names.push_back("A3S1G2");
  names.push_back("A3S1G2B");
  names.push_back("A3S1G2F");
  names.push_back("A3S1G2FB");
  names.push_back("A4S1G2");
  names.push_back("A4S1G2B");
  names.push_back("A4S1G2F");
  names.push_back("A4S1G2FB");
  names.push_back("A4S1G3");
  names.push_back("A4S1G3B");
  names.push_back("A4S1G3F");
  names.push_back("A4S1G3FB");
  names.push_back("A3S2G0");
  names.push_back("A3S2G0B");
  names.push_back("A3S2G0F");
  names.push_back("A3S2G0FB");
  names.push_back("A4S2G0");
  names.push_back("A4S2G0B");
  names.push_back("A4S2G0F");
  names.push_back("A4S2G0FB");
  names.push_back("A3S2G1");
  names.push_back("A3S2G1B");
  names.push_back("A3S2G1F");
  names.push_back("A3S2G1FB");
  names.push_back("A4S2G1");
  names.push_back("A4S2G1B");
  names.push_back("A4S2G1F");
  names.push_back("A4S2G1FB");
  names.push_back("A4S2G2");
  names.push_back("A4S2G2B");
  names.push_back("A4S2G2F");
  names.push_back("A4S2G2FB");
  names.push_back("A3S3");
  names.push_back("A3S3B");
  names.push_back("A3S3F");
  names.push_back("A3S3FB");
  names.push_back("A4S3G0");
  names.push_back("A4S3G0B");
  names.push_back("A4S3G0F");
  names.push_back("A4S3G0FB");
  names.push_back("A4S3G1");
  names.push_back("A4S3G1B");
  names.push_back("A4S3G1F");
  names.push_back("A4S3G1FB");
  names.push_back("A4S4");
  names.push_back("A4S4B");
  names.push_back("A4S4F");
  names.push_back("A4S4FB");
  names.push_back("A2Sg1G0");
  names.push_back("A2Sg1G0B");
  names.push_back("A2Sg1G0F");
  names.push_back("A2Sg1G0FB");
  names.push_back("A2Sg1G1");
  names.push_back("A2Sg1G1B");
  names.push_back("A2Sg1G1F");
  names.push_back("A2Sg1G1FB");
  names.push_back("A2Sg1S1");
  names.push_back("A2Sg1S1B");
  names.push_back("A2Sg1S1F");
  names.push_back("A2Sg1S1FB");
  names.push_back("A2Sg2");
  names.push_back("A2Sg2B");
  names.push_back("A2Sg2F");
  names.push_back("A2Sg2FB");
  names.push_back("A3Sg1G0");
  names.push_back("A3Sg1G0B");
  names.push_back("A3Sg1G0F");
  names.push_back("A3Sg1G0FB");
  names.push_back("A4Sg1G0");
  names.push_back("A4Sg1G0B");
  names.push_back("A4Sg1G0F");
  names.push_back("A4Sg1G0FB");
  names.push_back("A3Sg1G1");
  names.push_back("A3Sg1G1B");
  names.push_back("A3Sg1G1F");
  names.push_back("A3Sg1G1FB");
  names.push_back("A4Sg1G1");
  names.push_back("A4Sg1G1B");
  names.push_back("A4Sg1G1F");
  names.push_back("A4Sg1G1FB");
  names.push_back("A3Sg1G2");
  names.push_back("A3Sg1G2B");
  names.push_back("A3Sg1G2F");
  names.push_back("A3Sg1G2FB");
  names.push_back("A4Sg1G2");
  names.push_back("A4Sg1G2B");
  names.push_back("A4Sg1G2F");
  names.push_back("A4Sg1G2FB");
  names.push_back("A4Sg1G3");
  names.push_back("A4Sg1G3B");
  names.push_back("A4Sg1G3F");
  names.push_back("A4Sg1G3FB");
  names.push_back("A3Sg1S1G0");
  names.push_back("A3Sg1S1G0B");
  names.push_back("A3Sg1S1G0F");
  names.push_back("A3Sg1S1G0FB");
  names.push_back("A4Sg1S1G0");
  names.push_back("A4Sg1S1G0B");
  names.push_back("A4Sg1S1G0F");
  names.push_back("A4Sg1S1G0FB");
  names.push_back("A3Sg1S1G1");
  names.push_back("A3Sg1S1G1B");
  names.push_back("A3Sg1S1G1F");
  names.push_back("A3Sg1S1G1FB");
  names.push_back("A4Sg1S1G1");
  names.push_back("A4Sg1S1G1B");
  names.push_back("A4Sg1S1G1F");
  names.push_back("A4Sg1S1G1FB");
  names.push_back("A4Sg1S1G2");
  names.push_back("A4Sg1S1G2B");
  names.push_back("A4Sg1S1G2F");
  names.push_back("A4Sg1S1G2FB");
  names.push_back("A3Sg1S2");
  names.push_back("A3Sg1S2B");
  names.push_back("A3Sg1S2F");
  names.push_back("A3Sg1S2FB");
  names.push_back("A4Sg1S2G0");
  names.push_back("A4Sg1S2G0B");
  names.push_back("A4Sg1S2G0F");
  names.push_back("A4Sg1S2G0FB");
  names.push_back("A4Sg1S2G1");
  names.push_back("A4Sg1S2G1B");
  names.push_back("A4Sg1S2G1F");
  names.push_back("A4Sg1S2G1FB");
  names.push_back("A4Sg1S3");
  names.push_back("A4Sg1S3B");
  names.push_back("A4Sg1S3F");
  names.push_back("A4Sg1S3FB");
  names.push_back("A3Sg2G0");
  names.push_back("A3Sg2G0B");
  names.push_back("A3Sg2G0F");
  names.push_back("A3Sg2G0FB");
  names.push_back("A4Sg2G0");
  names.push_back("A4Sg2G0B");
  names.push_back("A4Sg2G0F");
  names.push_back("A4Sg2G0FB");
  names.push_back("A3Sg2G1");
  names.push_back("A3Sg2G1B");
  names.push_back("A3Sg2G1F");
  names.push_back("A3Sg2G1FB");
  names.push_back("A4Sg2G1");
  names.push_back("A4Sg2G1B");
  names.push_back("A4Sg2G1F");
  names.push_back("A4Sg2G1FB");
  names.push_back("A4Sg2G2");
  names.push_back("A4Sg2G2B");
  names.push_back("A4Sg2G2F");
  names.push_back("A4Sg2G2FB");
  names.push_back("A3Sg2S1");
  names.push_back("A3Sg2S1B");
  names.push_back("A3Sg2S1F");
  names.push_back("A3Sg2S1FB");
  names.push_back("A4Sg2S1G0");
  names.push_back("A4Sg2S1G0B");
  names.push_back("A4Sg2S1G0F");
  names.push_back("A4Sg2S1G0FB");
  names.push_back("A4Sg2S1G1");
  names.push_back("A4Sg2S1G1B");
  names.push_back("A4Sg2S1G1F");
  names.push_back("A4Sg2S1G1FB");
  names.push_back("A4Sg2S2");
  names.push_back("A4Sg2S2B");
  names.push_back("A4Sg2S2F");
  names.push_back("A4Sg2S2FB");
  names.push_back("A3Sg3");
  names.push_back("A3Sg3B");
  names.push_back("A3Sg3F");
  names.push_back("A3Sg3FB");
  names.push_back("A4Sg3G0");
  names.push_back("A4Sg3G0B");
  names.push_back("A4Sg3G0F");
  names.push_back("A4Sg3G0FB");
  names.push_back("A4Sg3G1");
  names.push_back("A4Sg3G1B");
  names.push_back("A4Sg3G1F");
  names.push_back("A4Sg3G1FB");
  names.push_back("A4Sg3S1");
  names.push_back("A4Sg3S1B");
  names.push_back("A4Sg3S1F");
  names.push_back("A4Sg3S1FB");
  names.push_back("A4Sg4");
  names.push_back("A4Sg4B");
  names.push_back("A4Sg4F");
  names.push_back("A4Sg4FB");
  names.push_back("A1G0M2");
  names.push_back("A1G0M2B");
  names.push_back("A1G0M2F");
  names.push_back("A1G0M2FB");
  names.push_back("A1G0");
  names.push_back("A1G0B");
  names.push_back("A1G0F");
  names.push_back("A1G0FB");
  names.push_back("A1G0M4");
  names.push_back("A1G0M4B");
  names.push_back("A1G0M4F");
  names.push_back("A1G0M4FB");
  names.push_back("A1G0M5");
  names.push_back("A1G0M5B");
  names.push_back("A1G0M5F");
  names.push_back("A1G0M5FB");
  names.push_back("A2G0M2");
  names.push_back("A2G0M2B");
  names.push_back("A2G0M2F");
  names.push_back("A2G0M2FB");
  names.push_back("A2G0M3");
  names.push_back("A2G0M3B");
  names.push_back("A2G0M3F");
  names.push_back("A2G0M3FB");
  names.push_back("A2G0M4");
  names.push_back("A2G0M4B");
  names.push_back("A2G0M4F");
  names.push_back("A2G0M4FB");
  names.push_back("A2G0M5");
  names.push_back("A2G0M5B");
  names.push_back("A2G0M5F");
  names.push_back("A2G0M5FB");
  names.push_back("A1G1M2");
  names.push_back("A1G1M2B");
  names.push_back("A1G1M2F");
  names.push_back("A1G1M2FB");
  names.push_back("A1G1");
  names.push_back("A1G1B");
  names.push_back("A1G1F");
  names.push_back("A1G1FB");
  names.push_back("A1G1M4");
  names.push_back("A1G1M4B");
  names.push_back("A1G1M4F");
  names.push_back("A1G1M4FB");
  names.push_back("A1G1M5");
  names.push_back("A1G1M5B");
  names.push_back("A1G1M5F");
  names.push_back("A1G1M5FB");
  names.push_back("A2G1M2");
  names.push_back("A2G1M2B");
  names.push_back("A2G1M2F");
  names.push_back("A2G1M2FB");
  names.push_back("A2G1M3");
  names.push_back("A2G1M3B");
  names.push_back("A2G1M3F");
  names.push_back("A2G1M3FB");
  names.push_back("A2G1M4");
  names.push_back("A2G1M4B");
  names.push_back("A2G1M4F");
  names.push_back("A2G1M4FB");
  names.push_back("A2G1M5");
  names.push_back("A2G1M5B");
  names.push_back("A2G1M5F");
  names.push_back("A2G1M5FB");
  names.push_back("A2G2M2");
  names.push_back("A2G2M2B");
  names.push_back("A2G2M2F");
  names.push_back("A2G2M2FB");
  names.push_back("A2G2M3");
  names.push_back("A2G2M3B");
  names.push_back("A2G2M3F");
  names.push_back("A2G2M3FB");
  names.push_back("A2G2M4");
  names.push_back("A2G2M4B");
  names.push_back("A2G2M4F");
  names.push_back("A2G2M4FB");
  names.push_back("A2G2M5");
  names.push_back("A2G2M5B");
  names.push_back("A2G2M5F");
  names.push_back("A2G2M5FB");
  names.push_back("A1S1M2");
  names.push_back("A1S1M2B");
  names.push_back("A1S1M2F");
  names.push_back("A1S1M2FB");
  names.push_back("A1S1");
  names.push_back("A1S1B");
  names.push_back("A1S1F");
  names.push_back("A1S1FB");
  names.push_back("A1S1M4");
  names.push_back("A1S1M4B");
  names.push_back("A1S1M4F");
  names.push_back("A1S1M4FB");
  names.push_back("A1S1M5");
  names.push_back("A1S1M5B");
  names.push_back("A1S1M5F");
  names.push_back("A1S1M5FB");
  names.push_back("A2S1G0M2");
  names.push_back("A2S1G0M2B");
  names.push_back("A2S1G0M2F");
  names.push_back("A2S1G0M2FB");
  names.push_back("A2S1G0M3");
  names.push_back("A2S1G0M3B");
  names.push_back("A2S1G0M3F");
  names.push_back("A2S1G0M3FB");
  names.push_back("A2S1G0M4");
  names.push_back("A2S1G0M4B");
  names.push_back("A2S1G0M4F");
  names.push_back("A2S1G0M4FB");
  names.push_back("A2S1G0M5");
  names.push_back("A2S1G0M5B");
  names.push_back("A2S1G0M5F");
  names.push_back("A2S1G0M5FB");
  names.push_back("A2S1G1M2");
  names.push_back("A2S1G1M2B");
  names.push_back("A2S1G1M2F");
  names.push_back("A2S1G1M2FB");
  names.push_back("A2S1G1M3");
  names.push_back("A2S1G1M3B");
  names.push_back("A2S1G1M3F");
  names.push_back("A2S1G1M3FB");
  names.push_back("A2S1G1M4");
  names.push_back("A2S1G1M4B");
  names.push_back("A2S1G1M4F");
  names.push_back("A2S1G1M4FB");
  names.push_back("A2S1G1M5");
  names.push_back("A2S1G1M5B");
  names.push_back("A2S1G1M5F");
  names.push_back("A2S1G1M5FB");
  names.push_back("A2S2M2");
  names.push_back("A2S2M2B");
  names.push_back("A2S2M2F");
  names.push_back("A2S2M2FB");
  names.push_back("A2S2M3");
  names.push_back("A2S2M3B");
  names.push_back("A2S2M3F");
  names.push_back("A2S2M3FB");
  names.push_back("A2S2M4");
  names.push_back("A2S2M4B");
  names.push_back("A2S2M4F");
  names.push_back("A2S2M4FB");
  names.push_back("A2S2M5");
  names.push_back("A2S2M5B");
  names.push_back("A2S2M5F");
  names.push_back("A2S2M5FB");
  names.push_back("A1Sg1M2");
  names.push_back("A1Sg1M2B");
  names.push_back("A1Sg1M2F");
  names.push_back("A1Sg1M2FB");
  names.push_back("A1Sg1");
  names.push_back("A1Sg1B");
  names.push_back("A1Sg1F");
  names.push_back("A1Sg1FB");
  names.push_back("A1Sg1M4");
  names.push_back("A1Sg1M4B");
  names.push_back("A1Sg1M4F");
  names.push_back("A1Sg1M4FB");
  names.push_back("A1Sg1M5");
  names.push_back("A1Sg1M5B");
  names.push_back("A1Sg1M5F");
  names.push_back("A1Sg1M5FB");
  names.push_back("A2Sg1G0M2");
  names.push_back("A2Sg1G0M2B");
  names.push_back("A2Sg1G0M2F");
  names.push_back("A2Sg1G0M2FB");
  names.push_back("A2Sg1G0M3");
  names.push_back("A2Sg1G0M3B");
  names.push_back("A2Sg1G0M3F");
  names.push_back("A2Sg1G0M3FB");
  names.push_back("A2Sg1G0M4");
  names.push_back("A2Sg1G0M4B");
  names.push_back("A2Sg1G0M4F");
  names.push_back("A2Sg1G0M4FB");
  names.push_back("A2Sg1G0M5");
  names.push_back("A2Sg1G0M5B");
  names.push_back("A2Sg1G0M5F");
  names.push_back("A2Sg1G0M5FB");
  names.push_back("A2Sg1G1M2");
  names.push_back("A2Sg1G1M2B");
  names.push_back("A2Sg1G1M2F");
  names.push_back("A2Sg1G1M2FB");
  names.push_back("A2Sg1G1M3");
  names.push_back("A2Sg1G1M3B");
  names.push_back("A2Sg1G1M3F");
  names.push_back("A2Sg1G1M3FB");
  names.push_back("A2Sg1G1M4");
  names.push_back("A2Sg1G1M4B");
  names.push_back("A2Sg1G1M4F");
  names.push_back("A2Sg1G1M4FB");
  names.push_back("A2Sg1G1M5");
  names.push_back("A2Sg1G1M5B");
  names.push_back("A2Sg1G1M5F");
  names.push_back("A2Sg1G1M5FB");
  names.push_back("A2Sg1S1M2");
  names.push_back("A2Sg1S1M2B");
  names.push_back("A2Sg1S1M2F");
  names.push_back("A2Sg1S1M2FB");
  names.push_back("A2Sg1S1M3");
  names.push_back("A2Sg1S1M3B");
  names.push_back("A2Sg1S1M3F");
  names.push_back("A2Sg1S1M3FB");
  names.push_back("A2Sg1S1M4");
  names.push_back("A2Sg1S1M4B");
  names.push_back("A2Sg1S1M4F");
  names.push_back("A2Sg1S1M4FB");
  names.push_back("A2Sg1S1M5");
  names.push_back("A2Sg1S1M5B");
  names.push_back("A2Sg1S1M5F");
  names.push_back("A2Sg1S1M5FB");
  names.push_back("A2Sg2M2");
  names.push_back("A2Sg2M2B");
  names.push_back("A2Sg2M2F");
  names.push_back("A2Sg2M2FB");
  names.push_back("A2Sg2M3");
  names.push_back("A2Sg2M3B");
  names.push_back("A2Sg2M3F");
  names.push_back("A2Sg2M3FB");
  names.push_back("A2Sg2M4");
  names.push_back("A2Sg2M4B");
  names.push_back("A2Sg2M4F");
  names.push_back("A2Sg2M4FB");
  names.push_back("A2Sg2M5");
  names.push_back("A2Sg2M5B");
  names.push_back("A2Sg2M5F");
  names.push_back("A2Sg2M5FB");
  names.push_back("M4");
  names.push_back("M4B");
  names.push_back("M4F");
  names.push_back("M4FB");
  names.push_back("M5");
  names.push_back("M5B");
  names.push_back("M5F");
  names.push_back("M5FB");
  names.push_back("M6");
  names.push_back("M6B");
  names.push_back("M6F");
  names.push_back("M6FB");
  names.push_back("M7");
  names.push_back("M7B");
  names.push_back("M7F");
  names.push_back("M7FB");
  names.push_back("M8");
  names.push_back("M8B");
  names.push_back("M8F");
  names.push_back("M8FB");
  names.push_back("M9");
  names.push_back("M9B");
  names.push_back("M9F");
  names.push_back("M9FB");
  names.push_back("Gn");
  names.push_back("GnF");
  names.push_back("M0");
  names.push_back("M0F");
  names.push_back("M1");
  names.push_back("M1B");
  names.push_back("M1F");
  names.push_back("M1FB");
  names.push_back("M2");
  names.push_back("M2B");
  names.push_back("M2F");
  names.push_back("M2FB");
  names.push_back("M3");
  names.push_back("M3B");
  names.push_back("M3F");
  names.push_back("M3FB");
  //omp_lock_t glock;
  //omp_init_lock(&glock);
  //#pragma omp parallel for  
  for ( int i = 0; i < names.size(); i++){
    string g = names[i];
    //cout << "create glycan:" << g << endl;
    Glycan* glycan = new Glycan(g);
    //omp_set_lock(&glock);
    glycanList->push_back(glycan); 
    (*glycanNameIndex)[g] = i;
    //omp_unset_lock(&glock);
    
  }
  //omp_destroy_lock(&glock);

  if(0){
    vector<Glycan*>::iterator ig;
    for ( ig = Glycan::glycanList->begin(); ig != Glycan::glycanList->end(); ig ++){
      cout << ">" ;
      (*ig)->print();
    }
  }
}

void Glycan::print() const {
  
  cout << "GLYCAN:" << getName() << "," << monoisotopicNeutralM() << endl;
  cout << "ID INDEX:" << endl;
  map<unsigned int, char>::iterator i1;
  for ( i1 = m_idx->begin(); i1 != m_idx->end(); i1 ++){
    cout << i1->first << "," << i1->second << endl;
  }
  cout << "Parent INDEX:" << endl;
  map<unsigned int, unsigned int>::iterator i2;
  for ( i2 = m_pidx->begin(); i2 != m_pidx->end(); i2 ++){
    cout << i2->first << "," << i2->second << endl;
  }
  
  cout << "FWT:" << endl;
  map<unsigned int, vector<unsigned int>* >::iterator i3;
  for ( i3 = m_fwt->begin(); i3 != m_fwt->end(); i3 ++){
    vector<unsigned int>* vpc = i3->second;
    cout << i3->first << ":";
    for (int i = 0; i < vpc->size(); i ++){
      cout << (*vpc)[i] << ",";
    }
    cout << endl;
  }

  cout << "Fragments(GY) " << m_yfs->size() << endl;
  map<string,double>::iterator i4;
  for ( i4 = m_yfs->begin(); i4 != m_yfs->end(); i4++){
    cout << i4->first << "," << i4->second << endl;
  }

  cout << "Fragments(GB) " << m_bfs->size() << endl;
  map<string,double>::iterator i5;
  for ( i5 =m_bfs->begin(); i5 != m_bfs->end(); i5++){
    cout << i5->first << "," << i5->second << endl;
  }
  cout << "END GLYCAN" << endl;
}

void Glycan::deleteTables(){
 
  if (monosaccharideList) {
    delete monosaccharideList;
  }
  monosaccharideList = NULL;

  if (monosaccharideSimpleIndex) {
    delete monosaccharideSimpleIndex;
  }
  monosaccharideSimpleIndex = NULL;

  if (glycanList) {
    delete glycanList;
  }
  glycanList = NULL ;
  if (glycanNameIndex) {
    delete glycanNameIndex;
  }
  glycanNameIndex = NULL;
}

string Glycan::getComposition() const {
  if (m_idx){
    map<string,int> composition;
    for (map< unsigned int, char>::iterator i = m_idx->begin(); i != m_idx->end(); i++) {
      if (monosaccharideSimpleIndex->find(i->second) != monosaccharideSimpleIndex->end()) {
	unsigned int id = (*monosaccharideSimpleIndex)[i->second];
	string comp = (*monosaccharideList)[id]->getCompCode();
	if (composition.find(comp) != composition.end()) {
	  composition[comp]++;
	} else {
	  composition[comp] = 1;
	}
      }
    }
    
    stringstream css;
    for (map<string,int>::iterator i = composition.begin(); i != composition.end(); i++) {
      css << i->first << i->second;
    }
    
    return (css.str());
    
  } else {
    
    return ("UNK");
  
  }
}

// Get all glycosidic fragments. For glycopeptides, Y ions will include the peptide. Set glycopeptideMass = 0 for fragments of just the glycan.
void Glycan::getGlycosidicFragmentIons(vector<FragmentIon*>& ions, unsigned int prominence, double glycopeptideMass, int glycopeptideCharge) const {
  
  double protonMass = (*Analyte::elementMonoisotopicMassTable)["+"];
  
  map<string,double>::iterator it;
  
  for (it = m_bfs->begin(); it != m_bfs->end(); it++) {

    int ch = 1; // only consider charge +1 for now
    double mz = (it->second + protonMass) / (double)ch;
    
    ions.push_back(new FragmentIon("GB:"+(it->first), 0, 0, mz, ch, prominence));

  }
  
  double peptideMass = 0.0;
  if (glycopeptideMass > 0.00001) peptideMass = glycopeptideMass - monoisotopicNeutralM();
 
  for (it = m_yfs->begin(); it != m_yfs->end(); it++) {
    int ch;
     for (ch = 1; ch <= glycopeptideCharge; ch++) {
      double mz = (it->second + peptideMass + protonMass * (double)ch) / (double)ch;
      ions.push_back(new FragmentIon("GY:"+(it->first), 0, 0, mz, ch, prominence));
      //cout << "ADD GY" << endl;
    }
  }
 
  // TODO: Other ions?
      /*
    ions.push_back(new FragmentIon("OX:F",0,0,147,1,prom));
    ions.push_back(new FragmentIon("OX:Hex",0,0,163,1,prom));
    ions.push_back(new FragmentIon("OX:NANA",0,0,292,1,prom));
    ions.push_back(new FragmentIon("OX:HexNAc",0,0,204,1,prom));
    ions.push_back(new FragmentIon("OX:Hex+HexNAc",0,0,366,1,prom));
    ions.push_back(new FragmentIon("OX:Hex+NANA",0,0,454,1,prom));
    ions.push_back(new FragmentIon("OX:Hex+HexNAc+NANA",0,0,657,1,prom));
    */
    // fragmentation only in peptide
    // generateFragmentIonsCID(ions);
    // fragmentation in peptide and both glycan
}
