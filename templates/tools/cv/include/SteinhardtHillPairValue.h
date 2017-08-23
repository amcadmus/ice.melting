#pragma once

#include "mace/CosHill.h"
#include "mace/SphericalHarmonics.h"

template<unsigned LL>
class SteinhardtHillPairValue 
{
public:
  SteinhardtHillPairValue () {};
  SteinhardtHillPairValue (const double & r0,
			   const double & r1,
			   const double & r2,
			   const double & r3);
  ~SteinhardtHillPairValue () {}
public:
  void reinit (const double & r0,
	       const double & r1,
	       const double & r2,
	       const double & r3);
  double rcut () const {return rc;}
  void evaluate (vector<double > & value,
		 const vector<double > & dist) const;
  unsigned valueDim() const {return interValueDim();}
private:
  unsigned interValueDim() const {return 2*(LL+1);}
  unsigned interDerivDim() const {return 6*(LL+1);}
  mutable CosHill<double> cs;
  mutable SphericalHarmonics<double, LL> sh;
  double rc;
  double rc2;
  double prefactor;
  // buffers for CosSwitch
  mutable double r1;
  mutable double cs_v;
  mutable double cs_d;
  // buffers for SphericalHarmonics
  mutable vector<double > sh_dist;
  mutable vector<double > sh_v;
  mutable vector<double > sh_d;  
};    

template<unsigned LL> 
SteinhardtHillPairValue<LL> ::
SteinhardtHillPairValue (const double & r0,
			 const double & r1,
			 const double & r2,
			 const double & r3)
{
  reinit (r0, r1, r2, r3);
}

template<unsigned LL> 
void
SteinhardtHillPairValue<LL> ::
reinit (const double & r0,
	const double & r1_,
	const double & r2,
	const double & r3)
{
  cs.reinit (r0, r1_, r2, r3);
  rc = r3;
  rc2 = rc * rc;
  prefactor = 1. * sqrt(4*M_PI / (2. * double(LL) + 1.));  
  // register CosSwitch buffer
  PointerArray<double> ai, av, ad;
  ai.push_back (&r1);
  av.push_back (&cs_v);
  ad.push_back (&cs_d);
  cs.registerData (ai, av, ad);
  // register SphericalHarmonics buffer
  sh_dist.resize (sh.inputDim());
  sh_v.resize (sh.valueDim());
  sh_d.resize (sh.derivDim());  
  PointerArray<double> bi, bv, bd;
  for (unsigned ii = 0; ii < sh_dist.size(); ++ii) bi.push_back(&sh_dist[ii]);
  for (unsigned ii = 0; ii < sh_v.size(); ++ii) bv.push_back(&sh_v[ii]);
  for (unsigned ii = 0; ii < sh_d.size(); ++ii) bd.push_back(&sh_d[ii]);
  sh.registerData (bi, bv, bd);
}

template<unsigned LL> 
void
SteinhardtHillPairValue<LL> ::
evaluate (vector<double > & value,
	  const vector<double > & dist) const
{
  assert (sh_dist.size() == dist.size());
  for (unsigned ii = 0; ii < sh_dist.size(); ++ii) sh_dist[ii] = dist[ii];
  value.resize(interValueDim());
  fill (value.begin(), value.end(), 0.);
  
  const double & dx (dist[0]);
  const double & dy (dist[1]);
  const double & dz (dist[2]);

  double r2 = dx * dx + dy * dy + dz * dz;
  if (r2 > rc2) {
    return;
  }
  double rinv  = 1./sqrtf(r2);
  r1 = r2 * rinv;
  
  cs.calculate();
  sh.calculate();
  
  double tfx = cs_d * dx * rinv;
  double tfy = cs_d * dy * rinv;
  double tfz = cs_d * dz * rinv;

  for (unsigned ii = 0; ii < interValueDim(); ++ii){
    double vsh = (sh_v[ii]);
    sh_v[ii] = sh_v[ii] * cs_v * prefactor;
    sh_d[ii*3+0] = ( sh_d[ii*3+0] * cs_v + vsh * tfx ) * prefactor;
    sh_d[ii*3+1] = ( sh_d[ii*3+1] * cs_v + vsh * tfy ) * prefactor;
    sh_d[ii*3+2] = ( sh_d[ii*3+2] * cs_v + vsh * tfz ) * prefactor;
  }

  for (unsigned ii = 0; ii < interValueDim(); ++ii) value[ii] = sh_v[ii];
}




