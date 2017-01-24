#pragma once

#include <vector> 
#include "mace/CosSwitch.h"

using namespace std;
using namespace MOASP;
using namespace MACE;

class CoordinationPairValue 
{
public:
  CoordinationPairValue () {};
  CoordinationPairValue (const double & rmin,
			 const double & rmax);
  ~CoordinationPairValue () {}
public:
  void reinit (const double & rmin,
	       const double & rmax);
  double rcut () const {return rc;}
  double evaluate (const vector<double > & dist) const;
  unsigned valueDim() const {return interValueDim();}
private:
  unsigned interValueDim() const {return 1;}
  unsigned interDerivDim() const {return 3;}
  mutable CosSwitch<double> cs;
  double rc;
  double rc2;
  double prefactor;
  // buffers for CosSwitch
  mutable double r1;
  mutable double cs_v;
  mutable double cs_d;
};    

CoordinationPairValue::
CoordinationPairValue (const double & rmin,
		       const double & rmax)
{
  reinit (rmin, rmax);
}

void
CoordinationPairValue::
reinit (const double & rmin,
	const double & rmax)
{
  cs.reinit (rmin, rmax);
  rc = rmax;
  rc2 = rc * rc;
  // register CosSwitch buffer
  PointerArray<double> ai, av, ad;
  ai.push_back (&r1);
  av.push_back (&cs_v);
  ad.push_back (&cs_d);
  cs.registerData (ai, av, ad);
}

double
CoordinationPairValue::
evaluate (const vector<double > & dist) const
{  
  const double & dx (dist[0]);
  const double & dy (dist[1]);
  const double & dz (dist[2]);

  double r2 = dx * dx + dy * dy + dz * dz;
  if (r2 > rc2) {
    return 0.;
  }
  double rinv  = 1./sqrtf(r2);
  r1 = r2 * rinv;
  
  cs.calculate();
  
  return cs_v;
}




