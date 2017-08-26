#include <cassert>
#include "Clusters.h"
#include "CellList.h"

Clusters::
Clusters (const double rcut_) 
{
  reinit (rcut_);
}

void 
Clusters::
reinit (const double rcut_) 
{
  clear();
  rcut = rcut_;
  rcut2 = rcut * rcut;
}

void
Clusters::
clear () 
{
  clstr_list.clear();
  atom_clstr_idx.clear();
}

void 
Clusters::
init_null (const int natoms)
{
  clear();
  // fill(atom_clstr_idx.begin(), atom_clstr_idx.end(), -1);
  for (int ii = 0; ii < natoms; ++ii){
    atom_clstr_idx.push_back (ii);
    clstr_list.push_back (vector<int> (1, ii));
  }
}

// void 
// Clusters::
// init_null (const int natoms)
// {
//   clstr_list.clear();
//   atom_clstr_idx.resize(natoms);
//   fill(atom_clstr_idx.begin(), atom_clstr_idx.end(), -1);
// }

void 
Clusters::
merge_clstr (const int i0, 
	     const int i1)
{
  clstr_list[i0].insert (clstr_list[i0].end(), clstr_list[i1].begin(), clstr_list[i1].end());
  for (unsigned ii = 0; ii < clstr_list[i1].size(); ++ii){
    int atom_idx = clstr_list[i1][ii];
    atom_clstr_idx[atom_idx] = i0;
  }
  clstr_list[i1].clear();
}

void 
Clusters::
add_clstr (const int aidx)
{
  atom_clstr_idx[aidx] = clstr_list.size();
  clstr_list.push_back (vector<int> (1, aidx));
}

void 
Clusters::
add_to_clstr (const int cidx, 
	      const int aidx)
{
  atom_clstr_idx[aidx] = cidx;
  clstr_list[cidx].push_back (aidx);
}

const int &
Clusters::
get_cluster_idx (const int aidx) const
{
  return atom_clstr_idx[aidx];
}

const vector<int > & 
Clusters::
get_atoms (const int cidx) const 
{
  return clstr_list[cidx];
}

void 
Clusters::
analysis (const vector<vector<double > > & atoms, 
	  const vector<double > & box)
{
  init_null (atoms.size());
  
  VectorType box_vec (box[0], box[1], box[2]);  
  CellList clist (atoms.size(), box_vec, rcut);
  clist.rebuild (atoms);
  
  double rup = rcut;
  int xiter = rup / clist.getCellSize().x;
  if (xiter * clist.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist.getCellSize().y;
  if (yiter * clist.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist.getCellSize().z;
  if (ziter * clist.getCellSize().z < rup) ziter ++;
  assert (xiter * clist.getCellSize().x >= rup);
  assert (yiter * clist.getCellSize().y >= rup);
  assert (ziter * clist.getCellSize().z >= rup);
  const vector<vector<unsigned > > & clist_data = clist.getList();
  for (unsigned icell = 0; icell < clist_data.size(); ++icell){
    const vector<unsigned > & ilist = clist.getList()[icell];
    vector<unsigned> nei_list_idx = 
	clist.neighboringCellIndex (icell, IntVectorType (xiter, yiter, ziter));
    for (unsigned jidx = 0; jidx < nei_list_idx.size(); ++jidx){
      unsigned jcell = nei_list_idx[jidx];
      const vector<unsigned > & jlist = clist.getList()[jcell];
      bool sameCell (icell == jcell);

      for (unsigned ii = 0; ii < ilist.size(); ++ii){
	int i_index = ilist[ii];
	// i belongs to a cluster, merge neighbors
	for (unsigned jj = 0; jj < jlist.size(); ++jj){
	  if (sameCell && ii == jj) continue;	    
	  int j_index = jlist[jj];
	  vector<double > diff (3);
	  for (int dd = 0; dd < 3; ++dd) {
	    diff[dd] = atoms[i_index][dd] - atoms[j_index][dd];
	  }
	  for (int dd = 0; dd < 3; ++dd){
	    if      (diff[dd] < -.5 * box[dd]) diff[dd] += box[dd];
	    else if (diff[dd] >= .5 * box[dd]) diff[dd] -= box[dd];
	  }
	  double r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
	  if (r2 < rcut2) {
	    if (atom_clstr_idx[i_index] != atom_clstr_idx[j_index]){
	      merge_clstr (atom_clstr_idx[i_index], atom_clstr_idx[j_index]);
	    }
	  }
	}
      }
    }
  }
}


// void 
// Clusters::
// analysis (const vector<vector<double > > & atoms, 
// 	  const vector<double > & box)
// {
//   init_null (atoms.size());
  
//   VectorType box_vec (box[0], box[1], box[2]);  
//   CellList clist (atoms.size(), box_vec, rcut);
//   clist.rebuild (atoms);
  
//   double rup = rcut;
//   int xiter = rup / clist.getCellSize().x;
//   if (xiter * clist.getCellSize().x < rup) xiter ++;
//   int yiter = rup / clist.getCellSize().y;
//   if (yiter * clist.getCellSize().y < rup) yiter ++;
//   int ziter = rup / clist.getCellSize().z;
//   if (ziter * clist.getCellSize().z < rup) ziter ++;
//   assert (xiter * clist.getCellSize().x >= rup);
//   assert (yiter * clist.getCellSize().y >= rup);
//   assert (ziter * clist.getCellSize().z >= rup);
//   const vector<vector<unsigned > > & clist_data = clist.getList();
//   for (unsigned icell = 0; icell < clist_data.size(); ++icell){
//     const vector<unsigned > & ilist = clist.getList()[icell];
//     vector<unsigned> nei_list_idx = 
// 	clist.neighboringCellIndex (icell, IntVectorType (xiter, yiter, ziter));
//     for (unsigned jidx = 0; jidx < nei_list_idx.size(); ++jidx){
//       unsigned jcell = nei_list_idx[jidx];
//       const vector<unsigned > & jlist = clist.getList()[jcell];
//       bool sameCell (icell == jcell);

//       for (unsigned ii = 0; ii < ilist.size(); ++ii){
// 	int i_index = ilist[ii];
// 	// search neighbors, to see if i belongs to a cluster
// 	if (atom_clstr_idx[i_index] == -1) {
// 	  for (unsigned jj = 0; jj < jlist.size(); ++jj){
// 	    if (sameCell && ii == jj) continue;	    
// 	    int j_index = jlist[jj];
// 	    vector<double > diff (3);
// 	    for (int dd = 0; dd < 3; ++dd) {
// 	      diff[dd] = atoms[i_index][dd] - atoms[j_index][dd];
// 	    }
// 	    for (int dd = 0; dd < 3; ++dd){
// 	      if      (diff[dd] < -.5 * box[dd]) diff[dd] += box[dd];
// 	      else if (diff[dd] >= .5 * box[dd]) diff[dd] -= box[dd];
// 	    }
// 	    double r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
// 	    if (r2 < rcut2) {
// 	      if (atom_clstr_idx[j_index] != -1){
// 		add_to_clstr (atom_clstr_idx[j_index], i_index);
// 		break;
// 	      }
// 	    }
// 	  }
// 	}
// 	// all neighbors are -1, create a new cluster
// 	if (atom_clstr_idx[i_index] == -1) {
// 	  add_clstr (i_index);
// 	}
// 	// i belongs to a cluster, merge neighbors
// 	for (unsigned jj = 0; jj < jlist.size(); ++jj){
// 	  if (sameCell && ii == jj) continue;	    
// 	  int j_index = jlist[jj];
// 	  vector<double > diff (3);
// 	  for (int dd = 0; dd < 3; ++dd) {
// 	    diff[dd] = atoms[i_index][dd] - atoms[j_index][dd];
// 	  }
// 	  for (int dd = 0; dd < 3; ++dd){
// 	    if      (diff[dd] < -.5 * box[dd]) diff[dd] += box[dd];
// 	    else if (diff[dd] >= .5 * box[dd]) diff[dd] -= box[dd];
// 	  }
// 	  double r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
// 	  if (r2 < rcut2) {
// 	    if (atom_clstr_idx[j_index] != -1){
// 	      if (atom_clstr_idx[i_index] != atom_clstr_idx[j_index]){
// 		merge_clstr (atom_clstr_idx[i_index], atom_clstr_idx[j_index]);
// 	      }
// 	    }
// 	    else {
// 	      add_to_clstr (atom_clstr_idx[i_index], j_index);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
// }


