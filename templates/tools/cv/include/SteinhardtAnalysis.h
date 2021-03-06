#pragma once

#include <cassert>
#include "CoordinationPairValue.h"
#include "SteinhardtPairValue.h"

using namespace std;

class SteinhardtAnalysisBase 
{
public:
  virtual void deposite (const CellList & clist,
			 const vector<double> box,
			 const vector<vector<double > > & coms,
			 const bool do_avg = false) = 0;
  virtual ~SteinhardtAnalysisBase () {};
  virtual void average () = 0;
  virtual double getStepQ () const = 0;
  virtual vector<double > getStepMole () const = 0;
  virtual vector<double > getAvgMole () const = 0;
};

template<unsigned LL>
class SteinhardtAnalysis : public SteinhardtAnalysisBase
{
public :
  SteinhardtAnalysis (const double rmin,
		      const double rmax,
		      const int func_numb_threads = 1);
  virtual ~SteinhardtAnalysis () {};
  void reinit (const double rmin,
	       const double rmax,
	       const int func_numb_threads = 1);
  virtual void deposite (const CellList & clist,
			 const vector<double> box,
			 const vector<vector<double > > & coms,
			 const bool do_avg);
  void deposite_mol (const CellList & clist,
		     const vector<double> box,
		     const vector<vector<double > > & waters);
  void deposite_avg (const CellList & clist,
		     const vector<double> box,
		     const vector<vector<double > > & waters);
  virtual void average ();
public:
  virtual double getStepQ () const {return step_Q;};
  virtual vector<double > getStepMole () const {return step_value;}
  virtual vector<double > getAvgMole  () const {return avg_value;}
private:
  void computeMolValue (vector<vector<double > > & mol_value,
			vector<double > & mol_coord,
			const CellList & clist,
			const vector<double> box,
			const vector<vector<double > > & waters);  
  void avgMolValue (vector<vector<double > > & avg_mol_q,
		    const CellList & clist,
		    const vector<double> box,
		    const vector<vector<double > > & waters,
		    const vector<vector<double > > & mol_q) ;
private :
  vector<double > avg_value;
  vector<double > step_value;
  double step_Q;
  int numb_step;
private:
  // SteinhardtPairValue<LL>	spv;
  // CoordinationPairValue		cpv;
  double rmin, rmax;
  int func_numb_threads;
}
    ;

template<unsigned LL>
SteinhardtAnalysis<LL> :: 
SteinhardtAnalysis (const double rmin,
		    const double rmax,
		    const int func_numb_threads_)
{
  reinit (rmin, rmax, func_numb_threads_);
}

template<unsigned LL>
void
SteinhardtAnalysis<LL> :: 
reinit (const double rmin_,
	const double rmax_,
	const int func_numb_threads_)
{
  rmin = rmin_;
  rmax = rmax_;
  // cpv.reinit (rmin, rmax);
  // spv.reinit (rmin, rmax);
  numb_step = 0;
  avg_value.clear();
  step_value.clear();
  step_Q = 0;
  func_numb_threads = func_numb_threads_;
}

template<unsigned LL>
void
SteinhardtAnalysis<LL> :: 
average ()
{
  if (numb_step == 0) return;
  for (unsigned ii = 0; ii < avg_value.size(); ++ii){
    avg_value[ii] /= double (numb_step);
  }
}

// template<unsigned LL>
// void
// SteinhardtAnalysis<LL> :: 
// localAvg (const vector<vector<double > > & avg_mol_value
// 	  const CellList & clist,
// 	  const vector<double> box,
// 	  const vector<vector<double > > & waters,
// 	  const vector<vector<double > > & mol_value) 
// {
// }

template<unsigned LL>
void
SteinhardtAnalysis<LL> :: 
computeMolValue (vector<vector<double > > & mol_value,
		 vector<double > & mol_coord,
		 const CellList & clist,
		 const vector<double> box,
		 const vector<vector<double > > & waters) 
{
  double rup = rmax;
  int xiter = rup / clist.getCellSize().x;
  if (xiter * clist.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist.getCellSize().y;
  if (yiter * clist.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist.getCellSize().z;
  if (ziter * clist.getCellSize().z < rup) ziter ++;
  assert (xiter * clist.getCellSize().x >= rup);
  assert (yiter * clist.getCellSize().y >= rup);
  assert (ziter * clist.getCellSize().z >= rup);

  // thread bufferes
  unsigned numb_water = waters.size();
  SteinhardtPairValue<LL> tmp_spv (rmin, rmax);
  vector<vector<vector<double > > > th_mol_value (func_numb_threads);
  vector<vector<double > >	    th_mol_coord (func_numb_threads);
  for (unsigned ii = 0; ii < th_mol_value.size(); ++ii){
    th_mol_value[ii].resize (numb_water);
    for (unsigned jj = 0; jj < th_mol_value[ii].size(); ++jj){
      th_mol_value[ii][jj].resize (tmp_spv.valueDim());
      fill (th_mol_value[ii][jj].begin(), th_mol_value[ii][jj].end(), 0.);
    }
    th_mol_coord[ii].resize (numb_water);
    fill (th_mol_coord[ii].begin(), th_mol_coord[ii].end(), 0.);
  }

  // loop
  IntVectorType nCell = clist.getNumCell();
  unsigned cellIndexUpper = unsigned(nCell.x * nCell.y * nCell.z);
#pragma omp parallel for num_threads (func_numb_threads) 
  for (int tt = 0; tt < func_numb_threads; ++tt){ 
    SteinhardtPairValue<LL>	spv (rmin, rmax);
    CoordinationPairValue	cpv (rmin, rmax);
    for (unsigned iCellIndex = tt;
	 iCellIndex < cellIndexUpper;
	 iCellIndex += func_numb_threads){
      const vector<unsigned> & iCellList (clist.getList()[iCellIndex]);
      vector<unsigned > neighborCellIndex =
	  clist.neighboringCellIndex (iCellIndex, IntVectorType (xiter, yiter, ziter));
      // loop of all neighboring cells of i
      for (unsigned iNeighborCellIndex = 0;
	   iNeighborCellIndex < neighborCellIndex.size();
	   ++iNeighborCellIndex){
	unsigned jCellIndex = neighborCellIndex[iNeighborCellIndex];
	const vector<unsigned> & jCellList (clist.getList()[jCellIndex]);
	bool sameCell (iCellIndex == jCellIndex);
	for (unsigned ii = 0; ii < iCellList.size(); ++ii){
	  int i_index = iCellList[ii];
	  for (unsigned jj = 0; jj < jCellList.size(); ++jj){
	    if (sameCell && ii == jj) continue;	    
	    int j_index = jCellList[jj];
	    if (i_index >= j_index) continue;
	    vector<double > io(3); 
	    for (int dd = 0; dd < 3; ++dd) io[dd] = waters[i_index][dd];
	    vector<double > jo(3);
	    for (int dd = 0; dd < 3; ++dd) jo[dd] = waters[j_index][dd];
	    vector<double > diff (3);
	    for (int dd = 0; dd < 3; ++dd) diff[dd] = jo[dd] - io[dd];
	    vector<int > shift(3, 0);
	    for (int dd = 0; dd < 3; ++dd){
	      if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
	      else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
	    }
	    for (int dd = 0; dd < 3; ++dd){
	      diff[dd] += box[dd] * shift[dd];
	    }
	    vector<double > pair_value;
	    spv.evaluate (pair_value, diff);
	    double pair_coord = cpv.evaluate (diff);
	    for (unsigned kk = 0; kk < spv.valueDim(); ++kk){
	      th_mol_value[tt][i_index][kk] += 0.5 * pair_value[kk];
	      th_mol_value[tt][j_index][kk] += 0.5 * pair_value[kk];
	    }
	    th_mol_coord[tt][i_index] += 0.5 * pair_coord;
	    th_mol_coord[tt][j_index] += 0.5 * pair_coord;
	  }
	}
      }
    }
  }

  mol_value.resize (numb_water);
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned jj = 0; jj < mol_value.size(); ++jj){
    mol_value[jj].resize (tmp_spv.valueDim());
    fill (mol_value[jj].begin(), mol_value[jj].end(), 0.);
  }
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned jj = 0; jj < mol_value.size(); ++jj){
    for (unsigned kk = 0; kk < mol_value[jj].size(); ++kk){
      for (int tt = 0; tt < func_numb_threads; ++tt){
	mol_value[jj][kk] += th_mol_value[tt][jj][kk];
      }
    }
  }
  mol_coord.resize (numb_water); 
  fill (mol_coord.begin(), mol_coord.end(), 0.);
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned jj = 0; jj < mol_coord.size(); ++jj){
    for (int tt = 0; tt < func_numb_threads; ++tt){
      mol_coord[jj] += th_mol_coord[tt][jj];
    }
  }  
}

template<unsigned LL>
void
SteinhardtAnalysis<LL> :: 
avgMolValue (vector<vector<double > > & avg_mol_q,
	     const CellList & clist,
	     const vector<double> box,
	     const vector<vector<double > > & waters,
	     const vector<vector<double > > & mol_q) 
{
  double rup = rmax;
  double rup2 = rmax * rmax;
  int xiter = rup / clist.getCellSize().x;
  if (xiter * clist.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist.getCellSize().y;
  if (yiter * clist.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist.getCellSize().z;
  if (ziter * clist.getCellSize().z < rup) ziter ++;
  assert (xiter * clist.getCellSize().x >= rup);
  assert (yiter * clist.getCellSize().y >= rup);
  assert (ziter * clist.getCellSize().z >= rup);

  unsigned numb_water = waters.size();
  unsigned ldof = mol_q[0].size();
  vector<unsigned > count_add (numb_water, 0);
  avg_mol_q.resize(numb_water);
  for (unsigned ii = 0; ii < avg_mol_q.size(); ++ii){
    avg_mol_q[ii].resize (ldof);
    fill (avg_mol_q[ii].begin(), avg_mol_q[ii].end(), 0.);
  }
  // loop
  IntVectorType nCell = clist.getNumCell();
  unsigned cellIndexUpper = unsigned(nCell.x * nCell.y * nCell.z);
#pragma omp parallel for num_threads (func_numb_threads) 
  for (int tt = 0; tt < func_numb_threads; ++tt){ 
    for (unsigned iCellIndex = tt;
	 iCellIndex < cellIndexUpper;
	 iCellIndex += func_numb_threads){
      const vector<unsigned> & iCellList (clist.getList()[iCellIndex]);
      vector<unsigned > neighborCellIndex =
	  clist.neighboringCellIndex (iCellIndex, IntVectorType (xiter, yiter, ziter));
      // loop of all neighboring cells of i
      for (unsigned iNeighborCellIndex = 0;
	   iNeighborCellIndex < neighborCellIndex.size();
	   ++iNeighborCellIndex){
	unsigned jCellIndex = neighborCellIndex[iNeighborCellIndex];
	const vector<unsigned> & jCellList (clist.getList()[jCellIndex]);
	for (unsigned ii = 0; ii < iCellList.size(); ++ii){
	  int i_index = iCellList[ii];
	  for (unsigned jj = 0; jj < jCellList.size(); ++jj){
	    // if (sameCell && ii == jj) continue;	    
	    int j_index = jCellList[jj];
	    const vector<double > & io(waters[i_index]);
	    const vector<double > & jo(waters[j_index]);	    
	    vector<double > diff(3);
	    for (int dd = 0; dd < 3; ++dd) diff[dd] = jo[dd] - io[dd];
	    vector<int > shift(3, 0);
	    for (int dd = 0; dd < 3; ++dd){
	      if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
	      else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
	    }
	    for (int dd = 0; dd < 3; ++dd){
	      diff[dd] += box[dd] * shift[dd];
	    }
	    double r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
	    if (r2 < rup2){
	      // deposite to i
	      for (unsigned kk = 0; kk < avg_mol_q[i_index].size(); ++kk){
		avg_mol_q[i_index][kk] += mol_q[j_index][kk];
	      }
	      count_add[i_index] ++;
	    }
	  }
	}
      }
    }
  }

#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned ii = 0; ii < avg_mol_q.size(); ++ii){
    for (unsigned kk = 0; kk < avg_mol_q[ii].size(); ++kk){
      avg_mol_q[ii][kk] /= double(count_add[ii]);
    }
  }

  // for (unsigned ii = 0; ii < avg_mol_q.size(); ++ii){
  //   cout << ii << " " << count_add[ii] << endl;
  // }
}

template<unsigned LL>
void
SteinhardtAnalysis<LL> :: 
deposite (const CellList & clist,
	  const vector<double> box,
	  const vector<vector<double > > & waters, 
	  const bool do_avg)
{  
  if (do_avg){
    deposite_avg (clist, box, waters);
  }
  else {
    deposite_mol (clist, box, waters);
  }
}

template<unsigned LL>
void
SteinhardtAnalysis<LL> :: 
deposite_mol (const CellList & clist,
	      const vector<double> box,
	      const vector<vector<double > > & waters)
{  
  unsigned numb_water = waters.size();
  if (step_value.size() != numb_water) {
    assert (numb_step == 0);
    step_value.resize (numb_water);
  }
  if (avg_value.size() != numb_water){
    assert (numb_step == 0);
    avg_value.clear();
    avg_value.resize(numb_water);
    fill(avg_value.begin(), avg_value.end(), 0.);
  }
  fill (step_value.begin(), step_value.end(), 0.);

  vector<vector<double > > mol_value;
  vector<double > mol_coord;
  computeMolValue (mol_value, mol_coord, clist, box, waters);
  
  SteinhardtPairValue<LL> tmp_spv (rmin, rmax);
  vector<double > sum_value (tmp_spv.valueDim(), 0.);
  for (unsigned kk = 0; kk < sum_value.size(); ++kk){
#pragma omp parallel for num_threads (func_numb_threads) 
    for (unsigned jj = 0; jj < mol_value.size(); ++jj){
      sum_value[kk] += mol_value[jj][kk];
    }
  }  
  double sum_coord (0);
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned jj = 0; jj < mol_coord.size(); ++jj){
    sum_coord += mol_coord[jj];
  }  
  double l2_norm (0);
  for (unsigned kk = 0; kk < sum_value.size(); ++kk){
    double pref = 2.; if (kk == 0 || kk == 1) pref = 1.;
    l2_norm += pref * sum_value[kk] * sum_value[kk];
  }
  l2_norm = sqrt(l2_norm);
  // cout << "step " << numb_step << " l2 " << l2_norm / sum_coord << endl;
  step_Q = l2_norm / sum_coord;

#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned ii = 0; ii < numb_water; ++ii){
    step_value[ii] = 0.;
    if (mol_coord[ii] == 0){
      cerr << "coordination number of mol " << ii << " is 0, some thing maybe wrong" << endl;
      continue;
    }
    for (unsigned kk = 0; kk < mol_value[ii].size(); ++kk){
      double tmp = mol_value[ii][kk] / mol_coord[ii];
      double pref = 2.; if (kk == 0 || kk == 1) pref = 1.;
      step_value[ii] += pref * tmp * tmp;
    }
    step_value[ii] = sqrt(step_value[ii]);
  }
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned ii = 0; ii < numb_water; ++ii){
    avg_value[ii] += step_value[ii];
  }

  numb_step ++;
}

template<unsigned LL>
void
SteinhardtAnalysis<LL> :: 
deposite_avg (const CellList & clist,
	      const vector<double> box,
	      const vector<vector<double > > & waters)
{  
  unsigned numb_water = waters.size();
  if (step_value.size() != numb_water) {
    assert (numb_step == 0);
    step_value.resize (numb_water);
  }
  if (avg_value.size() != numb_water){
    assert (numb_step == 0);
    avg_value.clear();
    avg_value.resize(numb_water);
    fill(avg_value.begin(), avg_value.end(), 0.);
  }
  fill (step_value.begin(), step_value.end(), 0.);

  vector<vector<double > > mol_value;
  vector<double > mol_coord;
  computeMolValue (mol_value, mol_coord, clist, box, waters);
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned ii = 0; ii < numb_water; ++ii){
    if (mol_coord[ii] == 0){
      cerr << "coordination number of mol " << ii << " is 0, some thing maybe wrong" << endl;
      continue;
    }
    for (unsigned kk = 0; kk < mol_value[ii].size(); ++kk){
      mol_value[ii][kk] = mol_value[ii][kk] / mol_coord[ii];
    }
  }
  // for (unsigned ii = 0; ii < mol_value.size(); ++ii){
  //   cout << ii << " a " << mol_coord[ii] << endl;;
  // }
  vector<vector<double > > avg_mol_q ;
  avgMolValue (avg_mol_q, clist, box, waters, mol_value);  

#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned ii = 0; ii < numb_water; ++ii){
    step_value[ii] = 0.;
    for (unsigned kk = 0; kk < mol_value[ii].size(); ++kk){
      double pref = 2.; if (kk == 0 || kk == 1) pref = 1.;
      step_value[ii] += pref * avg_mol_q[ii][kk] * avg_mol_q[ii][kk];
    }
    step_value[ii] = sqrt(step_value[ii]);
  }
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned ii = 0; ii < numb_water; ++ii){
    avg_value[ii] += step_value[ii];
  }

  numb_step ++;
}


