#pragma once

#include <cassert>
#include <armadillo>

#include "CellList.h"
#include "EigenDist.h"

using namespace std;

class EigenDistAnalysisBase 
{
public:
  virtual ~EigenDistAnalysisBase () {};
  virtual void deposite (const CellList & clist,
			 const vector<double> box,
			 const vector<vector<double > > & coms,
			 const bool do_avg) = 0;
  virtual void average () = 0;
  virtual vector<double > getStepMole () const = 0;
  virtual vector<double > getAvgMole  () const = 0;
};


template <typename MMatrixAssembler>
class EigenDistAnalysis : public EigenDistAnalysisBase
{
public :
  EigenDistAnalysis (const double rmax,
		     const int func_numb_threads = 1);
  virtual ~EigenDistAnalysis () {};
  void reinit (const double rmax,
	       const int func_numb_threads = 1);
  virtual void deposite (const CellList & clist,
			 const vector<double> box,
			 const vector<vector<double > > & coms,
			 const bool do_avg);
  virtual void average ();
public:
  virtual vector<double > getStepMole () const {return step_value;}
  virtual vector<double > getAvgMole  () const {return avg_value;}
private:
  void deposite_mol (const CellList & clist,
		     const vector<double> box,
		     const vector<vector<double > > & waters);
  void deposite_avg (const CellList & clist,
  		     const vector<double> box,
  		     const vector<vector<double > > & waters);
  void computeMolValue (vector<double > & mol_value,
			vector<int    > & mol_coord,
			const CellList & clist,
			const vector<double> box,
			const vector<vector<double > > & waters) ;
  void avgMolValue (vector<double > & avg_mol_q,
  		    const CellList & clist,
  		    const vector<double> box,
  		    const vector<vector<double > > & waters,
  		    const vector<double > & mol_q) ;
  double comp_eig (const double * eig,
		   const unsigned n_eig);
  void  process_ref ();
private :
  vector<double > avg_value;
  vector<double > step_value;
  vector<int    > step_coord;
  int numb_step;
private:
  double rmax;
  // vector <double> ref_eig;
  vector<vector<double > > ref_eig;
  int func_numb_threads;
}
    ;

