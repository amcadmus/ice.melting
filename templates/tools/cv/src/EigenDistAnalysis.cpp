#include "EigenDistAnalysis.h"

enum {
  global_coord_ref_size = 16,
  global_eig_ref_mat_size = 16
};
const double global_eig_ref_mat[global_eig_ref_mat_size][global_eig_ref_mat_size] = 
{
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{-0.22138,0.22138,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{-0.22138,-0.22138,0.44276,0,0,0,0,0,0,0,0,0,0,0,0,0},
{-0.22138,-0.22138,-0.22138,0.66414,0,0,0,0,0,0,0,0,0,0,0,0},
{-0.37224,-0.22138,-0.22138,-0.074489,0.88949,0,0,0,0,0,0,0,0,0,0,0},
{-0.41547,-0.26652,-0.22138,-0.17623,-0.039016,1.1186,0,0,0,0,0,0,0,0,0,0},
{-0.41832,-0.36658,-0.24787,-0.20962,-0.14177,0.10637,1.2778,0,0,0,0,0,0,0,0,0},
{-0.43338,-0.39622,-0.26901,-0.2431,-0.20443,-0.087705,0.14893,1.4849,0,0,0,0,0,0,0,0},
{-0.43386,-0.41011,-0.35555,-0.26722,-0.21609,-0.19852,0.061628,0.17295,1.6468,0,0,0,0,0,0,0},
{-0.43387,-0.43386,-0.37242,-0.29782,-0.21609,-0.21609,-0.19804,0.17294,0.17296,1.8223,0,0,0,0,0,0},
{-0.45092,-0.43387,-0.38998,-0.29782,-0.24895,-0.21609,-0.21271,-0.15971,0.17296,0.23037,2.0067,0,0,0,0,0},
{-0.4607,-0.44256,-0.40068,-0.29782,-0.28449,-0.23579,-0.21465,-0.20769,-0.11032,0.20418,0.25449,2.196,0,0,0,0},
{-0.4607,-0.4607,-0.41036,-0.29782,-0.28449,-0.28449,-0.22973,-0.20769,-0.20768,-0.054852,0.25449,0.2545,2.3895,0,0,0},
{-0.46111,-0.4607,-0.41927,-0.36708,-0.29782,-0.28449,-0.28169,-0.20789,-0.20768,-0.17275,0.075389,0.2545,0.28639,2.5442,0,0},
{-0.46111,-0.46109,-0.43798,-0.38832,-0.30807,-0.29782,-0.28379,-0.23422,-0.20773,-0.20339,-0.168,0.17163,0.28501,0.28701,2.7079,0},
{-0.46112,-0.46111,-0.45875,-0.39445,-0.30807,-0.30807,-0.29782,-0.23422,-0.23421,-0.20339,-0.20338,-0.16614,0.27732,0.28701,0.28703,2.8794}
};


const double global_coord_ref [global_coord_ref_size][3] = {
  {  0,0,2.68425},
{-2.1917,1.2654,-0.89475},
{2.1917,1.2654,-0.89475},
{-0.0000,-2.5308,-0.89475},
{-2.1917,3.7962,-0.0000},
{-4.3835,0.0000,-0.0000},
{4.3835,0.0000,-0.0000},
{2.1917,3.7962,-0.0000},
{-2.1917,-3.7962,-0.0000},
{2.1917,-3.7962,-0.0000},
{-2.1917,1.2654,-3.5791},
{2.1917,1.2654,-3.5791},
{-0.0000,-2.5308,-3.5791},
{-2.1917,1.2654,3.5791},
{2.1917,1.2654,3.5791},
{-0.0000,-2.5308,3.5791},
};


EigenDistAnalysis:: 
EigenDistAnalysis (const double rmax,
		   const vector<double > & ref_eig,
		   const int func_numb_threads_)
{
  reinit (rmax, ref_eig, func_numb_threads_);
  process_ref ();
}

void
EigenDistAnalysis:: 
reinit (const double rmax_,
	const vector<double > & ref_eig_,
	const int func_numb_threads_)
{
  rmax = rmax_;
  // ref_eig = ref_eig_;
  numb_step = 0;
  avg_value.clear();
  step_value.clear();
  step_coord.clear();
  step_Q = 0;
  func_numb_threads = func_numb_threads_;
}

void
EigenDistAnalysis:: 
process_ref ()
{
  // ref_eig.resize (global_eig_ref_mat_size);
  // for (unsigned ii = 0; ii < global_eig_ref_mat_size; ++ii) {
  //   for (unsigned jj = 0; jj < global_eig_ref_mat_size; ++jj) {
  //     ref_eig[ii].push_back (global_eig_ref_mat[ii][jj] * 10.);
  //   }
  // }
  vector<vector<double > > coord_ref (global_coord_ref_size+1, vector<double> (3,1.));
  for (int ii = 0; ii < global_coord_ref_size; ++ii){
    for (int dd = 0; dd < 3; ++dd){
      coord_ref[ii+1][dd] = global_coord_ref[ii][dd] * 0.1 * 1.0305 + coord_ref[0][dd];
    }
  }  
  vector<double > box (3, 5.);
  ref_eig.resize (global_eig_ref_mat_size, vector<double > (global_eig_ref_mat_size, 0));
  for (int ii = 0; ii < global_eig_ref_mat_size; ++ii){
    vector<int > nlist(ii+1);
    // ii + 1 neighbors
    for (int jj = 0; jj < ii+1; ++jj){
      nlist[jj] = jj+1;
    }
    arma::mat VV = arma::zeros<arma::mat>(ii+1, ii+1);
    assemble_trait_mat_dot (VV, coord_ref, box, 0, nlist);

    arma::vec eigs = arma::zeros<arma::vec>(ii+1);
    eig_sym(eigs, VV);
    // record the result
    for (int jj = 0; jj < ii+1; ++jj){
      ref_eig[ii][jj] = eigs(jj);
    }    
  }
}

inline double 
EigenDistAnalysis:: 
dist2 (const vector<double> & a1,
       const vector<double> & a2,
       const vector<double> & box)
{
  vector<double > diff (3);
  for (int dd = 0; dd < 3; ++dd) diff[dd] = a1[dd] - a2[dd];
  vector<int > shift(3, 0);
  for (int dd = 0; dd < 3; ++dd){
    if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
    else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
  }
  for (int dd = 0; dd < 3; ++dd){
    diff[dd] += box[dd] * shift[dd];
  }
  return diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
}

void 
EigenDistAnalysis::
assemble_trait_mat_dist (arma::mat & VV,
			 const vector<vector<double > > & waters,
			 const vector<double > & box,
			 const int & i_index,
			 const vector<int > & i_neigh_index)
{
  VV.zeros();
  for (unsigned j1 = 1; j1 < i_neigh_index.size(); ++j1){
    int j1_index = i_neigh_index[j1];
    for (unsigned j2 = 0; j2 < j1; ++j2){
      int j2_index = i_neigh_index[j2];
      double rij = sqrt (dist2 (waters[j1_index], waters[j2_index], box) );
      VV(j1, j2) = 1./rij;
    }
  }
  VV = VV + VV.t();  
}
		    
static inline
void
diff ( vector<double> & diff,
	const vector<double> & a1,
	const vector<double> & a2,
	const vector<double> & box)
{
  for (int dd = 0; dd < 3; ++dd) diff[dd] = a1[dd] - a2[dd];
  vector<int > shift(3, 0);
  for (int dd = 0; dd < 3; ++dd){
    if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
    else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
  }
  for (int dd = 0; dd < 3; ++dd){
    diff[dd] += box[dd] * shift[dd];
  }  
}

static inline 
double 
dot (const vector<double> & a1,
     const vector<double> & a2)
{
  return a1[0] * a2[0] + a1[1] * a2[1] + a1[2] * a2[2];
}


static inline 
double 
norm2 (const vector<double> & diff)
{
  return dot (diff, diff);
}

void 
EigenDistAnalysis::
assemble_trait_mat_angle (arma::mat & VV,
			  const vector<vector<double > > & waters,
			  const vector<double > & box,
			  const int & i_index,
			  const vector<int > & i_neigh_index)
{
  for (unsigned j1 = 0; j1 < i_neigh_index.size(); ++j1){
    VV(j1, j1) = 1.;
  }
  
  vector<double > diff1(3), diff2(3);
  for (unsigned j1 = 1; j1 < i_neigh_index.size(); ++j1){
    int j1_index = i_neigh_index[j1];
    diff (diff1, waters[j1_index], waters[i_index], box) ;    
    // for (unsigned kk = 0; kk < 3; ++kk) cout << diff1[kk] << " " ;
    // cout << endl;
    double r1 = sqrt(norm2 (diff1));
    for (unsigned j2 = 0; j2 < j1; ++j2){
      int j2_index = i_neigh_index[j2];
      diff (diff2, waters[j2_index], waters[i_index], box) ;    
      // for (unsigned kk = 0; kk < 3; ++kk) cout << diff2[kk] << " " ;
      // cout << endl;
      double r2 = sqrt(norm2 (diff2));
      VV(j1, j2) = dot (diff1, diff2) / r1 / r2;
      VV(j2, j1) = VV(j1, j2);
    }
  }
}
		    
void 
EigenDistAnalysis::
assemble_trait_mat_dot (arma::mat & VV,
			  const vector<vector<double > > & waters,
			  const vector<double > & box,
			  const int & i_index,
			  const vector<int > & i_neigh_index)
{
  for (unsigned j1 = 0; j1 < i_neigh_index.size(); ++j1){
    VV(j1, j1) = 0.;
  }
  
  vector<double > diff1(3), diff2(3);
  for (unsigned j1 = 1; j1 < i_neigh_index.size(); ++j1){
    int j1_index = i_neigh_index[j1];
    diff (diff1, waters[j1_index], waters[i_index], box) ;    
    for (unsigned j2 = 0; j2 < j1; ++j2){
      int j2_index = i_neigh_index[j2];
      diff (diff2, waters[j2_index], waters[i_index], box) ;    
      VV(j1, j2) = dot (diff1, diff2);
      VV(j2, j1) = VV(j1, j2);
    }
  }
}
		    

double 
EigenDistAnalysis:: 
comp_eig (const double * eig,
	  const unsigned n_eig)
{
  if (n_eig > global_eig_ref_mat_size){
    cerr << "matrix size larger than the known refs, exit" << endl;
    exit(1);
  }
  if (n_eig == 0) {
    cerr << "matrix of size 0, the cut-off may be too small" << endl;
  }
  double diff = 0;
  for (unsigned ii = 0; ii < n_eig; ++ii){
    // double tmp = eig[ii] / eig[n_eig-1] - ref_eig[n_eig-1][ii] / ref_eig[n_eig-1][n_eig-1];
    double tmp = eig[ii] - ref_eig[n_eig-1][ii];
    diff += tmp * tmp;
  }
  return sqrt(diff);
} 

void
EigenDistAnalysis:: 
average ()
{
  if (numb_step == 0) return;
  for (unsigned ii = 0; ii < avg_value.size(); ++ii){
    avg_value[ii] /= double (numb_step);
  }
}

void
EigenDistAnalysis:: 
computeMolValue (vector<double > & mol_value,
		 vector<int    > & mol_coord,
		 const CellList & clist,
		 const vector<double> box,
		 const vector<vector<double > > & waters) 
{
  double rup = rmax;
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
  mol_value.resize (numb_water);
  mol_coord.resize (numb_water);
  fill (mol_value.begin(), mol_value.end(), 0.);
  fill (mol_coord.begin(), mol_coord.end(), 0 );

  // loop
  IntVectorType nCell = clist.getNumCell();
  unsigned cellIndexUpper = unsigned(nCell.x * nCell.y * nCell.z);
  vector<vector<int > > i_neigh_index (numb_water);
  for (unsigned ii = 0; ii < i_neigh_index.size(); ++ii){
    i_neigh_index[ii].reserve (32);
  }

// #pragma omp parallel for num_threads (func_numb_threads) 
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
	  // find neighbors of i_index
	  for (unsigned jj = 0; jj < jCellList.size(); ++jj){
	    if (sameCell && ii == jj) continue;	    
	    int j_index = jCellList[jj];
	    double r2 = dist2 (waters[i_index], waters[j_index], box);
	    if (r2 < rup2){
	      i_neigh_index[i_index].push_back (j_index);
	    }
	  }
	}
      }
    }
  }

  for (unsigned i_index = 0; i_index < numb_water; ++i_index){    
    // assemble matrix for i_index
    unsigned nearest_n = i_neigh_index[i_index].size();
    arma::mat VV = arma::zeros<arma::mat>(nearest_n, nearest_n);
    arma::vec eigs = arma::zeros<arma::vec>(nearest_n);
    assemble_trait_mat_dot (VV, waters, box, i_index, i_neigh_index[i_index]);
    // compute eigen value for VV
    eig_sym(eigs, VV);
    // record the result
    mol_coord[i_index] = i_neigh_index[i_index].size();
    mol_value[i_index] = comp_eig (&eigs(0), mol_coord[i_index]);
    // for (unsigned kk = 0; kk < mol_coord[i_index]; ++kk){
    //   cout << eigs(kk) <<  " " ;
    // }
    // cout << endl;
    // cout << i_index << " " << mol_coord[i_index] << " "   << mol_value[i_index]<< endl;
  }
}

void
EigenDistAnalysis:: 
avgMolValue (vector<double > & avg_mol_q,
	     const CellList & clist,
	     const vector<double> box,
	     const vector<vector<double > > & waters,
	     const vector<double > & mol_q) 
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
  vector<unsigned > count_add (numb_water, 0);
  avg_mol_q.resize(numb_water);
  fill (avg_mol_q.begin(), avg_mol_q.end(), 0.);

  // loop
  IntVectorType nCell = clist.getNumCell();
  unsigned cellIndexUpper = unsigned(nCell.x * nCell.y * nCell.z);
// #pragma omp parallel for num_threads (func_numb_threads) 
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
EigenDistAnalysis:: 
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
EigenDistAnalysis:: 
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
EigenDistAnalysis:: 
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
