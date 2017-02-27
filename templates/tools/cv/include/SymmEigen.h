#pragma once

#include <vector>

using namespace std;

class SymmEigen 
{
public:
  SymmEigen (const int size,
	     const char uplo = 'U');
  void reinit (const int size,
	       const char uplo = 'U');
  void solve (double * eigen_value,
	      double * eigen_vector,
	      const double * matrix);
  void solve (double * eigen_value,
	      const double * matrix);
private:
  void check() const;
  int nn;
  char job;
  char uplo;
  vector<double > buff, ww, work;
  int lda, lwork, info;
}
    ;


