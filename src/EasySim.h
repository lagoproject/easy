//----------------------------------------------------------------------
/*
  EASYSIM Program - IPN Orsay since December 2002
  
  File EasySim.h 

 */
//----------------------------------------------------------------------

#ifndef EASYSIM_H
#define EASYSIM_H

using namespace std;

#include <string>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#ifndef CALIBONLY
#include "Shower_ROOT.h"
#include "ShowerGrnd_ROOT.h"
#include "PartGrnd_ROOT.h"
#include "Header_ROOT.h"
#include "HeaderAires_ROOT.h"
#endif
#include "EasySimConfig.h"

inline bool string_contains(string line, string s) {
    // ligb++ equivalent: return line.find(s);
    return line.find(s)==string::npos?false:true;
}

inline string string_before(string line, string seperator) {
    // ligb++ equivalent: return line.before(seperator);
    return line.substr(0,line.find(seperator));
}

inline string string_after(string line, string seperator) {
    // ligb++ equivalent: return line.before(seperator);
    return line.substr(line.find(seperator)+1);
}

#endif
