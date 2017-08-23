#pragma once

#include <vector> 
#include <cmath>
#include "mace/CosHill.h"

using namespace std;
using namespace MOASP;
using namespace MACE;

class CoordinationHillPairValue 
{
public:
  CoordinationHillPairValue () {};
  CoordinationHillPairValue (const double & r0,
			 const double & r1,
			 const double & r2,
			 const double & r3);
  ~CoordinationHillPairValue () {}
public:
  void reinit (const double & r0,
	       const double & r1,
	       const double & r2,
	       const double & r3);
  double rcut () const {return rc;}
  double evaluate (const vector<double > & dist) const;
  unsigned valueDim() const {return interValueDim();}
private:
  unsigned interValueDim() const {return 1;}
  unsigned interDerivDim() const {return 3;}
  mutable CosHill<double> cs;
  double rc;
  double rc2;
  double prefactor;
  // buffers for CosSwitch
  mutable double r1;
  mutable double cs_v;
  mutable double cs_d;
};    

CoordinationHillPairValue::
CoordinationHillPairValue (const double & r0,
		       const double & r1,
		       const double & r2,
		       const double & r3)
{
  reinit (r0, r1, r2, r3);
}

void
CoordinationHillPairValue::
reinit (const double & r0,
	const double & r1_,
	const double & r2,
	const double & r3)
{
  cs.reinit (r0, r1_, r2, r3);
  rc = r3;
  rc2 = rc * rc;
  // register CosSwitch buffer
  PointerArray<double> ai, av, ad;
  ai.push_back (&r1);
  av.push_back (&cs_v);
  ad.push_back (&cs_d);
  cs.registerData (ai, av, ad);
}

double
CoordinationHillPairValue::
evaluate (const vector<double > & dist) const
{  
  const double & dx (dist[0]);
  const double & dy (dist[1]);
  const double & dz (dist[2]);

  double r2 = dx * dx + dy * dy + dz * dz;
  if (r2 > rc2) {
    return 0.;
  }
  double rinv  = 1./sqrt(r2);
  r1 = r2 * rinv;
  
  cs.calculate();
  
  return cs_v;
}




