#pragma once

#include <cassert>
#include <armadillo>

#include "CellList.h"

using namespace std;

class EigenDistAnalysis 
{
public :
  EigenDistAnalysis (const double rmax,
		     const vector<double > & ref_eig,
		     const int func_numb_threads = 1);
  ~EigenDistAnalysis () {};
  void reinit (const double rmax,
	       const vector<double > & ref_eig,
	       const int func_numb_threads = 1);
  void deposite (const CellList & clist,
		 const vector<double> box,
		 const vector<vector<double > > & coms,
		 const bool do_avg);
  void deposite_mol (const CellList & clist,
		     const vector<double> box,
		     const vector<vector<double > > & waters);
  void deposite_avg (const CellList & clist,
  		     const vector<double> box,
  		     const vector<vector<double > > & waters);
  void average ();
public:
  double getStepQ () const {return step_Q;};
  vector<double > getStepMole () const {return step_value;}
  vector<double > getAvgMole  () const {return avg_value;}
private:
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
  double dist2 (const vector<double> & a1,
		const vector<double> & a2,
		const vector<double> & box);
  double comp_eig (const double * eig,
		   const unsigned n_eig);
  void  process_ref ();
private :
  vector<double > avg_value;
  vector<double > step_value;
  vector<int    > step_coord;
  double step_Q;
  int numb_step;
private:
  double rmax;
  // vector <double> ref_eig;
  vector<vector<double > > ref_eig;
  int func_numb_threads;
}
    ;

