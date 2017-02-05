#pragma once

#include <complex>
#include <vector>
#include <cassert>

using namespace std;

#define NFFT_PRECISION_DOUBLE
#include "nfft3mp.h"

class StructureFactor 
{
public:
  StructureFactor (const unsigned * KK,
		   const unsigned nmol, 
		   const unsigned func_numb_threads);
  ~StructureFactor ();
  // void getVecQ (vector<vector<double > > & q);
  void getNormQ (vector<double > & q,
		 const vector<double> & box);
  void compute (vector<complex<double > > & sq,
		const vector<vector<double> > & coms,
		const vector<double> & box);
private:
  void init (const unsigned * KK,
	       const unsigned nmol, 
	       const unsigned func_numb_threads);  
  void normalize_coord (vector<double > & coord,
			const vector<vector<double> > & coms,
			const vector<double> & box);
  unsigned KK[3];
  unsigned NK;
  unsigned nmol;
  unsigned func_numb_threads;
  NFFT(plan) p;
}
    ;

StructureFactor::
StructureFactor (const unsigned * KK,
		 const unsigned nmol, 
		 const unsigned func_numb_threads)
{
  init (KK, nmol, func_numb_threads);
}

StructureFactor::
~StructureFactor ()
{
  NFFT(finalize)(&p);
}

void
StructureFactor::
init (const unsigned * KK_,
	const unsigned nmol_, 
	const unsigned func_numb_threads_)
{
  for (unsigned dd = 0; dd < 3; ++dd) KK[dd] = KK_[dd];
  NK = KK[0] * KK[1] * KK[2];
  nmol = nmol_;
  func_numb_threads = func_numb_threads_;

  NFFT(init_3d)(&p, KK[0], KK[1], KK[2], nmol);  
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);
  if (p.flags & PRE_ONE_PSI) NFFT(precompute_one_psi)(&p);

  // cout << p.M_total << endl;
  // cout << p.N_total << endl;
  const char *error_str;
  error_str = NFFT(check)(&p);
  if (error_str != 0)
  {
    printf("Error in nfft module: %s\n", error_str);
    exit (1);
  }
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
  for (int ii = 0; ii < int(KK[0]); ++ii){
    for (int jj = 0; jj < int(KK[1]); ++jj){
      for (int kk = 0; kk < int(KK[2]); ++kk){
	double qq[3];
	// qq[0] = ii * 2 * M_PI * boxi[0];
	// qq[1] = jj * 2 * M_PI * boxi[1];
	// qq[2] = kk * 2 * M_PI * boxi[2];
	qq[0] = (ii-int(KK[0]/2)) * 2 * M_PI * boxi[0];
	qq[1] = (jj-int(KK[1]/2)) * 2 * M_PI * boxi[1];
	qq[2] = (kk-int(KK[2]/2)) * 2 * M_PI * boxi[2];
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
  // for (unsigned ii = 0; ii < nmol; ++ii){
  //   cout << coord[3*ii] << " " 
  // 	 << coord[3*ii+1] << " " 
  // 	 << coord[3*ii+2] << " " 
  // 	 << endl;
  // }

  assert (NK == p.N_total);
  assert (nmol == p.M_total);

  for (unsigned ii = 0; ii < 3 * nmol; ++ii){
    p.x[ii] = coord[ii];
  }
  for (unsigned ii = 0; ii < nmol; ++ii){
    p.f[ii][0] = 1.;
    p.f[ii][1] = 0.;
  }

  NFFT(adjoint)(&p);

  sq.resize (NK);
  for (unsigned ii = 0; ii < NK; ++ii){
    sq[ii].real (p.f_hat[ii][0]);
    sq[ii].imag (p.f_hat[ii][1]);
  }
}


// void
// StructureFactor::
// compute (vector<complex<double > > & sq,
// 	 const vector<vector<double> > & coms,
// 	 const vector<double> & box)
// {
//   vector<double > coord;
//   normalize_coord (coord, coms, box);
  
//   sq.resize (NK);
//   for (unsigned ii = 0; ii < NK; ++ii){
//     sq[ii].real(0);
//     sq[ii].imag(0);
//   }
  
//   vector<complex<double > >::iterator iter = sq.begin();

// // #pragma omp parallel for collapse(3) num_threads (func_numb_threads) 
//   for (unsigned ii = 0; ii < KK[0]; ++ii){
//     for (unsigned jj = 0; jj < KK[1]; ++jj){
//       for (unsigned kk = 0; kk < KK[2]; ++kk){
// 	double sum0, sum1;
// 	sum0 = sum1 = 0.;
// #pragma omp parallel for num_threads (func_numb_threads) reduction (+:sum0,sum1)
// 	for (unsigned pp = 0; pp < nmol; ++pp){
// 	  double tmp = 2. * M_PI * (ii * coord[3*pp+0] + jj * coord[3*pp+1] + kk * coord[3*pp+2]);
// 	  sum0 += cos(tmp);
// 	  sum1 += sin(tmp);
// 	}
// 	iter->real( iter->real() + sum0 );
// 	iter->imag( iter->imag() + sum1 );
// 	iter ++;
//       }
//     }
//   }
// }


