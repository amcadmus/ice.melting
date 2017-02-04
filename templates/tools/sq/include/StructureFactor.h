#pragma once

#include <complex>
#include <vector>
#include <cassert>

using namespace std;

class StructureFactor 
{
public:
  StructureFactor (const unsigned * KK,
		   const unsigned nmol, 
		   const unsigned func_numb_threads);
  void reinit (const unsigned * KK,
	       const unsigned nmol, 
	       const unsigned func_numb_threads);  
  // void getVecQ (vector<vector<double > > & q);
  void getNormQ (vector<double > & q,
		 const vector<double> & box);
  void compute (vector<complex<double > > & sq,
		const vector<vector<double> > & coms,
		const vector<double> & box);
private:
  void normalize_coord (vector<double > & coord,
			const vector<vector<double> > & coms,
			const vector<double> & box);
  unsigned KK[3];
  unsigned NK;
  unsigned nmol;
  unsigned func_numb_threads;
}
    ;

StructureFactor::
StructureFactor (const unsigned * KK,
		 const unsigned nmol, 
		 const unsigned func_numb_threads)
{
  reinit (KK, nmol, func_numb_threads);
}

void
StructureFactor::
reinit (const unsigned * KK_,
	const unsigned nmol_, 
	const unsigned func_numb_threads_)
{
  for (unsigned dd = 0; dd < 3; ++dd) KK[dd] = KK_[dd];
  NK = KK[0] * KK[1] * KK[2];
  nmol = nmol_;
  func_numb_threads = func_numb_threads_;
}


void
StructureFactor::
getNormQ (vector<double > & q,
	  const vector<double> & box)
{
  vector<double > boxi(3);
  for (unsigned dd = 0; dd < 3; ++dd) boxi[dd] = 1./box[dd];
  q.resize (NK);
  vector<double >::iterator iter = q.begin();
  for (unsigned ii = 0; ii < KK[0]; ++ii){
    for (unsigned jj = 0; jj < KK[1]; ++jj){
      for (unsigned kk = 0; kk < KK[2]; ++kk){
	double qq[3];
	qq[0] = ii * 2 * M_PI * boxi[0];
	qq[1] = jj * 2 * M_PI * boxi[1];
	qq[2] = kk * 2 * M_PI * boxi[2];
	*(iter ++) = sqrt (qq[0] * qq[0] + qq[1] * qq[1] + qq[2] * qq[2]);
      }
    }
  }
}

void
StructureFactor::
normalize_coord (vector<double > & coord,
		 const vector<vector<double> > & coms,
		 const vector<double> & box)
{
  coord.resize (3 * nmol);
  assert (coms.size() == nmol);
  vector<double > boxi(3);
  for (unsigned dd = 0; dd < 3; ++dd) boxi[dd] = 1./box[dd];
  vector<double >::iterator iter = coord.begin();
  for (unsigned ii = 0; ii < coms.size(); ++ii){
    for (unsigned dd = 0; dd < 3; ++dd){
      *iter = coms[ii][dd] * boxi[dd] ;
      if      (*iter <  0) *iter += 1;
      else if (*iter >= 1) *iter -= 1;
      *iter -= 0.5;
      iter ++;
    }
  }
}

void
StructureFactor::
compute (vector<complex<double > > & sq,
	 const vector<vector<double> > & coms,
	 const vector<double> & box)
{
  vector<double > coord;
  normalize_coord (coord, coms, box);
  
  sq.resize (NK);
  for (unsigned ii = 0; ii < NK; ++ii){
    sq[ii].real(0);
    sq[ii].imag(0);
  }
  
  vector<complex<double > >::iterator iter = sq.begin();

// #pragma omp parallel for collapse(3) num_threads (func_numb_threads) 
  for (unsigned ii = 0; ii < KK[0]; ++ii){
    for (unsigned jj = 0; jj < KK[1]; ++jj){
      for (unsigned kk = 0; kk < KK[2]; ++kk){
	double sum0, sum1;
	sum0 = sum1 = 0.;
#pragma omp parallel for num_threads (func_numb_threads) reduction (+:sum0,sum1)
	for (unsigned pp = 0; pp < nmol; ++pp){
	  double tmp = 2. * M_PI * (ii * coord[3*pp+0] + jj * coord[3*pp+1] + kk * coord[3*pp+2]);
	  sum0 += cos(tmp);
	  sum1 += sin(tmp);
	}
	iter->real( iter->real() + sum0 );
	iter->imag( iter->imag() + sum1 );
	iter ++;
      }
    }
  }
}


