#include "SymmEigen.h"

#include <iostream> 
#include <stdlib.h>
#include <string.h>

extern "C" {
    void dsyev_ (char *, char *, int *, double *, int *, double *, double *, int *, int *);
}

SymmEigen::
SymmEigen (const int size,
	   const char uplo)
{
  reinit (size, uplo);
}

void
SymmEigen::
check () const
{
  if (info != 0){
    cerr << "error happens in SymmEigen" << endl;
    exit(1);
  }
}

void
SymmEigen::
reinit (const int size,
	const char uplo_)
{
  nn = size;
  job = 'V';
  if (uplo_ == 'U') {
    uplo = 'L';
  }
  else {
    uplo = 'U';
  }
  lda = nn;
  ww.resize(nn);
  buff.resize(nn*nn);

  double tmpwork[1] = {0};
  lwork = -1;  
  dsyev_ (&job, &uplo, &nn, &buff[0], &lda, &ww[0], tmpwork, &lwork, &info);
  check ();

  lwork = int(tmpwork[0]);
  work.resize (lwork);
}

void
SymmEigen::
solve (double * eigen_value,
       const double * matrix)
{
  memcpy (&buff[0], matrix, sizeof(double) * nn * nn);
  job = 'N';
  dsyev_ (&job, &uplo, &nn, &buff[0], &lda, &ww[0], &work[0], &lwork, &info);
  memcpy (eigen_value,  &ww[0], sizeof(double) * nn);
}

void
SymmEigen::
solve (double * eigen_value,
       double * eigen_vector,
       const double * matrix)
{
  memcpy (&buff[0], matrix, sizeof(double) * nn * nn);
  job = 'V';
  dsyev_ (&job, &uplo, &nn, &buff[0], &lda, &ww[0], &work[0], &lwork, &info);
  memcpy (eigen_value,  &ww[0], sizeof(double) * nn);
  for (int ii = 0; ii < nn; ++ii){
    for (int jj = 0; jj < nn; ++jj){
      eigen_vector[ii*nn+jj] = buff[jj*lda+ii];
    }
  }  
  // memcpy (eigen_vector, &buff[0], sizeof(double) * nn * nn);  
}



