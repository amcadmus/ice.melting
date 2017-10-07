#pragma once

#include <cassert>

using namespace std;

template<typename HydrogenBond>
class HydrogenBondAnalysis 
{
public :
  HydrogenBondAnalysis (const typename HydrogenBond::Parameters & param,
			const unsigned numb_water,
			const int func_numb_threads = 1);
  void reinit (const unsigned numb_water,
	       const int func_numb_threads = 1);
  void deposite (const CellList & clist,
		 const vector<double> box,
		 const vector<vector<vector<double > > > & waters);
  void average ();
  void computeBondList (vector<vector<int > > & bond_list,
			const CellList & clist,
			const vector<double> box,
			const vector<vector<vector<double > > > & waters) const;
  void findDLPair (vector<pair<int, int > > & dpair,
		   vector<pair<int, int > > & lpair,
		   const CellList & clist,
		   const vector<double> box,
		   const vector<vector<vector<double > > > & waters) const;
public :
  vector<double > avg_count_acc;
  vector<double > avg_count_don;
  vector<int > step_count_acc;
  vector<int > step_count_don;
  int numb_step;
private:
  HydrogenBond hydrogen_bond;
  int func_numb_threads;
}
    ;

template<typename HydrogenBond>
HydrogenBondAnalysis<HydrogenBond> :: 
HydrogenBondAnalysis (const typename HydrogenBond::Parameters & param,
		      const unsigned numb_water, 
		      const int func_numb_threads_) 
    : hydrogen_bond (param)
{
  reinit (numb_water, func_numb_threads_);
}

template<typename HydrogenBond>
void
HydrogenBondAnalysis<HydrogenBond> :: 
reinit (const unsigned numb_water, 
	const int func_numb_threads_) 
{
  numb_step = 0;
  avg_count_acc.clear();
  avg_count_don.clear();
  avg_count_acc.resize (numb_water, 0.);
  avg_count_don.resize (numb_water, 0.);
  step_count_acc.clear();
  step_count_don.clear();
  step_count_acc.resize (numb_water, 0);
  step_count_don.resize (numb_water, 0);
  func_numb_threads = func_numb_threads_;
}

template<typename HydrogenBond>
void
HydrogenBondAnalysis<HydrogenBond> :: 
average ()
{
  if (numb_step == 0) return;
  for (unsigned ii = 0; ii < avg_count_don.size(); ++ii){
    avg_count_don[ii] /= double (numb_step);
    avg_count_acc[ii] /= double (numb_step);
  }
}

inline 
double 
dist (const vector<double> & x, 
      const vector<double> & y)
{
  vector<double > diff(3);
  for (int dd = 0; dd < 3; ++dd) diff[dd] = x[dd] - y[dd];
  double dr = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
  return sqrt(dr);
}

template<typename HydrogenBond>
void
HydrogenBondAnalysis<HydrogenBond> :: 
deposite (const CellList & clist,
	  const vector<double> box,
	  const vector<vector<vector<double > > > & waters)
{
  fill (step_count_acc.begin(), step_count_acc.end(), 0);
  fill (step_count_don.begin(), step_count_don.end(), 0);

  double rup = hydrogen_bond.getRcut();
  int xiter = rup / clist.getCellSize().x;
  if (xiter * clist.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist.getCellSize().y;
  if (yiter * clist.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist.getCellSize().z;
  if (ziter * clist.getCellSize().z < rup) ziter ++;
  assert (xiter * clist.getCellSize().x >= rup);
  assert (yiter * clist.getCellSize().y >= rup);
  assert (ziter * clist.getCellSize().z >= rup);

  IntVectorType nCell = clist.getNumCell();
  
  unsigned numb_water = waters.size();
  assert (numb_water == avg_count_acc.size());
  assert (numb_water == avg_count_don.size());
  assert (numb_water == step_count_acc.size());
  assert (numb_water == step_count_don.size());
  
  vector<vector<int > > thd_count_acc (func_numb_threads);
  vector<vector<int > > thd_count_don (func_numb_threads);
  for (int tt = 0; tt < func_numb_threads; ++tt) {
    thd_count_acc[tt].resize (numb_water, 0);
    thd_count_don[tt].resize (numb_water, 0);
  }
  
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
	bool sameCell (iCellIndex == jCellIndex);
	for (unsigned ii = 0; ii < iCellList.size(); ++ii){
	  int i_index = iCellList[ii];
	  for (unsigned jj = 0; jj < jCellList.size(); ++jj){
	    if (sameCell && ii == jj) continue;	    
	    int j_index = jCellList[jj];
	    if (i_index >= j_index) continue;
	    vector<double > io(3); 
	    for (int dd = 0; dd < 3; ++dd) io[dd] = waters[i_index][0][dd];
	    vector<double > jo(3);
	    for (int dd = 0; dd < 3; ++dd) jo[dd] = waters[j_index][0][dd];
	    vector<double > diff (3);
	    for (int dd = 0; dd < 3; ++dd) diff[dd] = jo[dd] - io[dd];
	    // if (dist(waters[i_index][0], waters[i_index][1]) > 0.12) {
	    //   cout << "wrong water" << endl;
	    // }
	    // if (dist(waters[i_index][0], waters[i_index][2]) > 0.12) {
	    //   cout << "wrong water" << endl;
	    // }
	    vector<int > shift(3, 0);
	    for (int dd = 0; dd < 3; ++dd){
	      if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
	      else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
	    }
	    vector<vector<double > > j_water = waters[j_index];
	    for (int kk = 0; kk < 3; ++kk){
	      for (int dd = 0; dd < 3; ++dd){
		j_water[kk][dd] += box[dd] * shift[dd];
	      }
	    }
	    // if (
	    if (hydrogen_bond (waters[i_index], j_water)){
	      thd_count_don[tt][i_index] ++;
	      thd_count_acc[tt][j_index] ++;
	    }
	    if (hydrogen_bond (j_water, waters[i_index])){
	      thd_count_acc[tt][i_index] ++;
	      thd_count_don[tt][j_index] ++;
	    }
	  }
	}
      }
    }
  }
  for (int tt = 0; tt < func_numb_threads; ++tt){
    for (unsigned ii = 0; ii < numb_water; ++ii){
      step_count_acc[ii] += thd_count_acc[tt][ii];
      step_count_don[ii] += thd_count_don[tt][ii];
    }
  }
  for (unsigned ii = 0; ii < numb_water; ++ii){
    avg_count_acc[ii] += step_count_acc[ii];
    avg_count_don[ii] += step_count_don[ii];
    // cout << " ii " << ii 
    // 	 << " acc " << step_count_acc[ii]
    // 	 << " don " << step_count_don[ii]
    // 	 << endl;
  }
  numb_step ++;
}


template<typename HydrogenBond>
void
HydrogenBondAnalysis<HydrogenBond> :: 
computeBondList (vector<vector<int > > & bond_list,
		 const CellList & clist,
		 const vector<double> box,
		 const vector<vector<vector<double > > > & waters) const
{  
  double rup = hydrogen_bond.getRcut();
  int xiter = rup / clist.getCellSize().x;
  if (xiter * clist.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist.getCellSize().y;
  if (yiter * clist.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist.getCellSize().z;
  if (ziter * clist.getCellSize().z < rup) ziter ++;
  assert (xiter * clist.getCellSize().x >= rup);
  assert (yiter * clist.getCellSize().y >= rup);
  assert (ziter * clist.getCellSize().z >= rup);

  IntVectorType nCell = clist.getNumCell();
  
  unsigned numb_water = waters.size();
  assert (numb_water == avg_count_acc.size());
  assert (numb_water == avg_count_don.size());
  assert (numb_water == step_count_acc.size());
  assert (numb_water == step_count_don.size());
  
  bond_list.clear();
  bond_list.resize (numb_water);
  
  unsigned cellIndexUpper = unsigned(nCell.x * nCell.y * nCell.z);

  for (unsigned iCellIndex = 0;
       iCellIndex < cellIndexUpper;
       iCellIndex += 1){
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
	  for (int dd = 0; dd < 3; ++dd) io[dd] = waters[i_index][0][dd];
	  vector<double > jo(3);
	  for (int dd = 0; dd < 3; ++dd) jo[dd] = waters[j_index][0][dd];
	  vector<double > diff (3);
	  for (int dd = 0; dd < 3; ++dd) diff[dd] = jo[dd] - io[dd];
	  vector<int > shift(3, 0);
	  for (int dd = 0; dd < 3; ++dd){
	    if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
	    else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
	  }
	  vector<vector<double > > j_water = waters[j_index];
	  for (int kk = 0; kk < 3; ++kk){
	    for (int dd = 0; dd < 3; ++dd){
	      j_water[kk][dd] += box[dd] * shift[dd];
	    }
	  }
	  // un-directional
	  if (hydrogen_bond (waters[i_index], j_water) ||
	      hydrogen_bond (j_water, waters[i_index]) ){
	    bond_list[i_index].push_back (j_index);
	    bond_list[j_index].push_back (i_index);
	  }
	  // // directional donnor -> acceptor
	  // if (hydrogen_bond (waters[i_index], j_water)){
	  //   bond_list[i_index].push_back (j_index);	    
	  // }
	  // if (hydrogen_bond (j_water, waters[i_index])){
	  //   bond_list[j_index].push_back (i_index);
	  // }
	}
      }
    }
  }

  for (unsigned ii = 0; ii < bond_list.size(); ++ii){
    sort (bond_list[ii].begin(), bond_list[ii].end());
  }
}



template<typename HydrogenBond>
void
HydrogenBondAnalysis<HydrogenBond> :: 
findDLPair (vector<pair<int, int > > & dpair,
	    vector<pair<int, int > > & lpair,
	    const CellList & clist,
	    const vector<double> box,
	    const vector<vector<vector<double > > > & waters) const
{
  HydrogenBond_Geo_1::Parameters param (hydrogen_bond.getRcut(), 60);
  HydrogenBond_Geo_1 hb_rule_1 (param);

  dpair.clear();
  lpair.clear();
  
  double rup = hydrogen_bond.getRcut();
  int xiter = rup / clist.getCellSize().x;
  if (xiter * clist.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist.getCellSize().y;
  if (yiter * clist.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist.getCellSize().z;
  if (ziter * clist.getCellSize().z < rup) ziter ++;
  assert (xiter * clist.getCellSize().x >= rup);
  assert (yiter * clist.getCellSize().y >= rup);
  assert (ziter * clist.getCellSize().z >= rup);

  IntVectorType nCell = clist.getNumCell();
  
  unsigned numb_water = waters.size();
  assert (numb_water == avg_count_acc.size());
  assert (numb_water == avg_count_don.size());
  assert (numb_water == step_count_acc.size());
  assert (numb_water == step_count_don.size());
  
  unsigned cellIndexUpper = unsigned(nCell.x * nCell.y * nCell.z);

  vector<int > count_acc (numb_water, 0);
  vector<int > count_don (numb_water, 0);  

  for (unsigned iCellIndex = 0;
       iCellIndex < cellIndexUpper;
       iCellIndex += 1){
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
	  for (int dd = 0; dd < 3; ++dd) io[dd] = waters[i_index][0][dd];
	  vector<double > jo(3);
	  for (int dd = 0; dd < 3; ++dd) jo[dd] = waters[j_index][0][dd];
	  vector<double > diff (3);
	  for (int dd = 0; dd < 3; ++dd) diff[dd] = jo[dd] - io[dd];
	  vector<int > shift(3, 0);
	  for (int dd = 0; dd < 3; ++dd){
	    if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
	    else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
	  }
	  vector<vector<double > > j_water = waters[j_index];
	  for (int kk = 0; kk < 3; ++kk){
	    for (int dd = 0; dd < 3; ++dd){
	      j_water[kk][dd] += box[dd] * shift[dd];
	    }
	  }
	  // if (
	  if (hydrogen_bond (waters[i_index], j_water)){
	    count_don[i_index] ++;
	    count_acc[j_index] ++;
	  }
	  if (hydrogen_bond (j_water, waters[i_index])){
	    count_acc[i_index] ++;
	    count_don[j_index] ++;
	  }
	}
      }
    }
  }


  for (unsigned iCellIndex = 0;
       iCellIndex < cellIndexUpper;
       iCellIndex += 1){
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
	  for (int dd = 0; dd < 3; ++dd) io[dd] = waters[i_index][0][dd];
	  vector<double > jo(3);
	  for (int dd = 0; dd < 3; ++dd) jo[dd] = waters[j_index][0][dd];
	  vector<double > diff (3);
	  for (int dd = 0; dd < 3; ++dd) diff[dd] = jo[dd] - io[dd];
	  vector<int > shift(3, 0);
	  for (int dd = 0; dd < 3; ++dd){
	    if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
	    else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
	  }
	  vector<vector<double > > j_water = waters[j_index];
	  for (int kk = 0; kk < 3; ++kk){
	    for (int dd = 0; dd < 3; ++dd){
	      j_water[kk][dd] += box[dd] * shift[dd];
	    }
	  }
	  const vector<vector<double > > & i_water = waters[i_index];
	  // un-directional
	  if ( hb_rule_1 (i_water, j_water) 
	       && hb_rule_1 (j_water, i_water)	       
	       && count_don[i_index] >= 1
	       && count_don[j_index] >= 1
	       && count_acc[i_index] == 2
	       && count_acc[j_index] == 2 
	      ) {
	    // cout << i_index << " " 
	    // 	 << j_index << " " 
	    // 	 << count_don[i_index] << " " 
	    // 	 << count_don[j_index] << " " 
	    // 	 << count_acc[i_index] << " " 
	    // 	 << count_acc[j_index] << " "  
	    // 	 <<  endl;
	    dpair.push_back (pair<int, int> (i_index, j_index));
	    // cout << i_index << " " << j_index << endl;
	  }
	  for (int dd = 0; dd < 3; ++dd) diff[dd] = i_water[0][dd] - j_water[0][dd];
	  double rr = 0;
	  for (int dd = 0; dd < 3; ++dd) rr += diff[dd] * diff[dd];
	  rr = sqrt(rr);
	  if (rr < hydrogen_bond.getRcut() 
	      && (! hydrogen_bond (i_water, j_water))
	      && (! hydrogen_bond (j_water, i_water))
	      && count_don[i_index] == 2
	      && count_don[j_index] == 2
	      && count_acc[i_index] == 1
	      && count_acc[j_index] == 1
	      ){
	    // cout << count_don[i_index] << " " 
	    // 	 << count_don[j_index] << " " 
	    // 	 << count_acc[i_index] << " " 
	    // 	 << count_acc[j_index] << " "  
	    // 	 <<  endl;
	    lpair.push_back (pair<int, int> (i_index, j_index));
	  }   
	}
      }
    }
  }
}

