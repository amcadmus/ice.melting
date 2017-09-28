#include "Histogram.h"
#include <iostream>
#include <cassert>

Histogram::
Histogram (const unsigned nbin,
	   const double lower,
	   const double upper)
{
  reinit (nbin, lower, upper);
}

void
Histogram::
reinit (const unsigned nbin_,
	const double lower_,
	const double upper_)
{
  setNumbBins (nbin_);
  setBoundary (lower_, upper_);
  
  hh = (upper - lower) / double(n_bin);
  
  value.clear();
  value.resize (n_bin, 0.0);
}

void
Histogram::
setBoundary (const double lower_,
	     const double upper_)
{
  upper = upper_;
  lower = lower_;
  if (lower >  upper) {
    swap (upper, lower);
  }
}

void
Histogram::
processData (const vector<double > & vec) 
{
  if (lower == upper){
    lower = vec[0];
    upper = vec[0];
    for (unsigned ii = 1; ii < vec.size(); ++ii){
      if (vec[ii] < lower) lower = vec[ii];
      if (vec[ii] > upper) upper = vec[ii] + 1e-6 * vec[ii];
    }
    reinit (n_bin, lower, upper);
  }
  
  double sum = 0;
  for (unsigned ii = 0; ii < vec.size(); ++ii){
    int index = (vec[ii] - lower) / hh;
    if (index >= 0 && index < int(n_bin)){
      value[index] += 1.;
      sum += 1;
    }
  }
  if (sum != 0){
    for (unsigned ii = 0; ii < n_bin; ++ii){
      value[ii] = value[ii] / sum / hh;
    }
  }
}

void
Histogram::
getResult (vector<double > & xx,
	   vector<double > & vv)
{
  xx.resize (n_bin);
  vv.resize (n_bin);

  for (unsigned ii = 0; ii < n_bin; ++ii){
    xx[ii] = lower + (ii + 0.5) * hh;
    vv[ii] = value[ii];
  }
}


Histogram2D::
Histogram2D (const unsigned nbinx_,
	const unsigned nbiny_,
	const double lowerx_,
	const double lowery_,
	const double upperx_,
	const double uppery_)
{
  reinit (nbinx_, nbiny_, lowerx_, lowery_, upperx_, uppery_);
}

void
Histogram2D::
reinit (const unsigned nbinx_,
	const unsigned nbiny_,
	const double lowerx_,
	const double lowery_,
	const double upperx_,
	const double uppery_)
{
  setNumbBins (nbinx_, nbiny_);
  setBoundary (lowerx_, lowery_, upperx_, uppery_);
  
  hx = (upperx - lowerx) / double(nbinx);
  hy = (uppery - lowery) / double(nbiny);
  
  value.clear();
  value.resize (nbinx, vector<double > (nbiny, 0.));
}

void
Histogram2D::
setBoundary (const double lowerx_,
	     const double lowery_,
	     const double upperx_,
	     const double uppery_)
{
  upperx = upperx_;
  lowerx = lowerx_;
  if (lowerx >  upperx) {
    swap (upperx, lowerx);
  }

  uppery = uppery_;
  lowery = lowery_;
  if (lowery >  uppery) {
    swap (uppery, lowery);
  }
}

void
Histogram2D::
processData (const vector<double > & vecx, 
	     const vector<double > & vecy) 
{
  assert (vecx.size() == vecy.size());
  bool need_reinit = false;
  if (lowerx == upperx){
    lowerx = vecx[0];
    upperx = vecx[0];
    for (unsigned ii = 1; ii < vecx.size(); ++ii){
      if (vecx[ii] < lowerx) lowerx = vecx[ii];
      if (vecx[ii] > upperx) upperx = vecx[ii] + 1e-6 * vecx[ii];
    }
    need_reinit = true;
  }
  if (lowery == uppery){
    lowery = vecy[0];
    uppery = vecy[0];
    for (unsigned ii = 1; ii < vecy.size(); ++ii){
      if (vecy[ii] < lowery) lowery = vecy[ii];
      if (vecy[ii] > uppery) uppery = vecy[ii] + 1e-6 * vecy[ii];
    }
    need_reinit = true;
  }
  if (need_reinit) reinit (nbinx, nbiny, lowerx, lowery, upperx, uppery);
  
  double sum = 0;
  for (unsigned ii = 0; ii < vecx.size(); ++ii){
    int index_x = (vecx[ii] - lowerx) / hx;
    int index_y = (vecy[ii] - lowery) / hy;
    if (index_x >= 0 && index_x < int(nbinx) && 
	index_y >= 0 && index_y < int(nbiny) ){
      value[index_x][index_y] += 1.;
      sum += 1;
    }
  }
  if (sum != 0){
    for (unsigned ii = 0; ii < nbinx; ++ii){
      for (unsigned jj = 0; jj < nbiny; ++jj){
	value[ii][jj] = value[ii][jj] / sum / hx / hy;
      }
    }
  }
}

void
Histogram2D::
getResult (vector<double > & xx,
	   vector<double > & yy,
	   vector<vector<double > > & vv)
{
  xx.resize (nbinx);
  yy.resize (nbiny);
  vv.clear ();
  vv.resize (nbinx, vector<double> (nbiny, 0.));

  for (unsigned ii = 0; ii < nbinx; ++ii){
    xx[ii] = lowerx + (ii + 0.5) * hx;
  }
  for (unsigned ii = 0; ii < nbiny; ++ii){
    yy[ii] = lowery + (ii + 0.5) * hy;
  }
  vv = value;
}
