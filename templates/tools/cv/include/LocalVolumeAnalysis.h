#pragma once

#include <cassert>
#include <vector>
#include "CellList.h"
#include "mace/CosHill.h"

using namespace std;
using namespace MOASP;
using namespace MACE;

class LocalVolumeAnalysis
{
public :
  LocalVolumeAnalysis (const double r0,
		       const double r1,
		       const double r2,
		       const double r3,
		       const int func_numb_threads = 1);
  ~LocalVolumeAnalysis () {};
  void reinit (const double r0,
	       const double r1,
	       const double r2,
	       const double r3,
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
  vector<double > getStepMole () const {return step_value;}
  vector<double > getAvgMole  () const {return avg_value;}
private:
  void computeMolValue (vector<double > & mol_value,
			vector<double > & mol_coord,
			const CellList & clist,
			const vector<double> box,
			const vector<vector<double > > & waters);  
  void avgMolValue (vector<double > & avg_mol_q,
		    const CellList & clist,
		    const vector<double> box,
		    const vector<vector<double > > & waters,
		    const vector<double > & mol_q) ;
private :
  vector<double > avg_value;
  vector<double > step_value;
  vector<double > step_coord;
  int numb_step;
private:
  // LocalVolumePairValue<LL>	spv;
  // CoordinationPairValue		cpv;
  double r0, r1, r2, r3;
  int func_numb_threads;
}
    ;

