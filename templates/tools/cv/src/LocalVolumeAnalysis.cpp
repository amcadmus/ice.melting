#include "LocalVolumeAnalysis.h"

LocalVolumeAnalysis :: 
LocalVolumeAnalysis (const double r0,
		     const double r1,
		     const double r2,
		     const double r3,
		     const int func_numb_threads_)
{
  reinit (r0, r1, r2, r3, func_numb_threads_);
}


void
LocalVolumeAnalysis :: 
reinit (const double r0_,
	const double r1_,
	const double r2_,
	const double r3_,
	const int func_numb_threads_)
{
  r0 = r0_;
  r1 = r1_;
  r2 = r2_;
  r3 = r3_;
  numb_step = 0;
  avg_value.clear();
  step_value.clear();
  func_numb_threads = func_numb_threads_;
}


void
LocalVolumeAnalysis :: 
average ()
{
  if (numb_step == 0) return;
  for (unsigned ii = 0; ii < avg_value.size(); ++ii){
    avg_value[ii] /= double (numb_step);
  }
}

void
LocalVolumeAnalysis :: 
computeMolValue (vector<double > & mol_value,
		 vector<double > & mol_coord,
		 const CellList & clist,
		 const vector<double> box,
		 const vector<vector<double > > & waters) 
{
  double rup = r3;
  double rup2 = rup * rup;
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
  vector<vector<double > > th_mol_value (func_numb_threads);
  vector<vector<double > > th_mol_coord (func_numb_threads);
  for (unsigned ii = 0; ii < th_mol_value.size(); ++ii){
    th_mol_value[ii].resize (numb_water);
    th_mol_coord[ii].resize (numb_water);
    fill (th_mol_value[ii].begin(), th_mol_value[ii].end(), 0.);
    fill (th_mol_coord[ii].begin(), th_mol_coord[ii].end(), 0.);
  }

  // loop
  IntVectorType nCell = clist.getNumCell();
  unsigned cellIndexUpper = unsigned(nCell.x * nCell.y * nCell.z);
#pragma omp parallel for num_threads (func_numb_threads) 
  for (int tt = 0; tt < func_numb_threads; ++tt){ 
    CosHill<double > hill;
    double rr, hill_v, hill_d;
    hill.reinit (r0, r1, r2, r3);
    rr = 0;
    hill_v = 0;
    hill_d = 0;
    PointerArray<double> ai, av, ad;
    ai.push_back (&rr);
    av.push_back (&hill_v);
    ad.push_back (&hill_d);
    hill.registerData (ai, av, ad);    
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
	    // vector<double > io(3); 
	    // for (int dd = 0; dd < 3; ++dd) io[dd] = waters[i_index][dd];
	    // vector<double > jo(3);
	    // for (int dd = 0; dd < 3; ++dd) jo[dd] = waters[j_index][dd];
	    // vector<double > diff (3);
	    // for (int dd = 0; dd < 3; ++dd) diff[dd] = jo[dd] - io[dd];
	    // vector<int > shift(3, 0);
	    // for (int dd = 0; dd < 3; ++dd){
	    //   if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
	    //   else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
	    // }
	    // for (int dd = 0; dd < 3; ++dd){
	    //   diff[dd] += box[dd] * shift[dd];
	    // }
	    double r2 = dist2 (waters[i_index], waters[j_index], box);
	    if (r2 < rup2) {
	      rr = sqrt(r2);
	      hill.calculate ();
	      th_mol_value[tt][i_index] += hill_v * rr;
	      th_mol_value[tt][j_index] += hill_v * rr;
	      th_mol_coord[tt][i_index] += hill_v;
	      th_mol_coord[tt][j_index] += hill_v;
	    }
	  }
	}
      }
    }
  }

  mol_value.resize (numb_water);
  fill (mol_value.begin(), mol_value.end(), 0.);
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned jj = 0; jj < mol_value.size(); ++jj){
    for (int tt = 0; tt < func_numb_threads; ++tt){
      mol_value[jj] += th_mol_value[tt][jj];
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
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned jj = 0; jj < mol_coord.size(); ++jj){
    mol_value[jj] /= mol_coord[jj];
  }
}


void
LocalVolumeAnalysis :: 
avgMolValue (vector<double > & avg_mol_q,
	     const CellList & clist,
	     const vector<double> box,
	     const vector<vector<double > > & waters,
	     const vector<double > & mol_q) 
{
  double rup = r3;
  double rup2 = r3 * r3;
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
  vector<unsigned > count_add (numb_water, 0);
  avg_mol_q.resize(numb_water);
  fill (avg_mol_q.begin(), avg_mol_q.end(), 0.);

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
	      avg_mol_q[i_index] += mol_q[j_index];
	      count_add[i_index] ++;
	    }
	  }
	}
      }
    }
  }

  for (unsigned ii = 0; ii < avg_mol_q.size(); ++ii){
    avg_mol_q[ii] /= double(count_add[ii]);
  }
}

void
LocalVolumeAnalysis :: 
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


void
LocalVolumeAnalysis::
deposite_mol (const CellList & clist,
	      const vector<double> box,
	      const vector<vector<double > > & waters)
{  
  unsigned numb_water = waters.size();
  if (avg_value.size() != numb_water){
    assert (numb_step == 0);
    avg_value.clear();
    avg_value.resize(numb_water);
    fill(avg_value.begin(), avg_value.end(), 0.);
  }

  computeMolValue (step_value, step_coord, clist, box, waters);

  for (unsigned ii = 0; ii < numb_water; ++ii){
    avg_value[ii] += step_value[ii];
  }  

  numb_step ++;
}

void
LocalVolumeAnalysis::
deposite_avg (const CellList & clist,
	      const vector<double> box,
	      const vector<vector<double > > & waters)
{  
  unsigned numb_water = waters.size();
  if (avg_value.size() != numb_water){
    assert (numb_step == 0);
    avg_value.clear();
    avg_value.resize(numb_water);
    fill(avg_value.begin(), avg_value.end(), 0.);
  }

  vector<double > mol_value;
  computeMolValue (mol_value, step_coord, clist, box, waters);
  avgMolValue (step_value, clist, box, waters, mol_value);

  for (unsigned ii = 0; ii < numb_water; ++ii){
    avg_value[ii] += step_value[ii];
  }  

  numb_step ++;
}
