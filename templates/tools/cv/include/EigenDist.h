#pragma once

#include "Defines.h"
#include <vector>
#include <armadillo>

using namespace std;

class assemble_trait_mat_dist {
public: 
  void operator () (vector<double > & VV,
		    const vector<vector<double > > & waters,
		    const vector<double > & box,
		    const int & i_index,
		    const vector<int > & i_neigh_index)
      {
	unsigned mat_size = i_neigh_index.size();
	VV.resize (mat_size * mat_size);
	fill (VV.begin(), VV.end(), 0.);
	// only assemble lower part
	for (unsigned j1 = 1; j1 < i_neigh_index.size(); ++j1){
	  int j1_index = i_neigh_index[j1];
	  for (unsigned j2 = 0; j2 < j1; ++j2){
	    int j2_index = i_neigh_index[j2];
	    double rij = sqrt (dist2 (waters[j1_index], waters[j2_index], box) );
	    VV[j1*mat_size+j2] = 1./rij;
	  }
	}
	// VV = VV + VV.t();  
      }
};

    
class assemble_trait_mat_angle {
public:
  void operator () (vector<double > & VV,
		    const vector<vector<double > > & waters,
		    const vector<double > & box,
		    const int & i_index,
		    const vector<int > & i_neigh_index)
      {
	unsigned mat_size = i_neigh_index.size();
	VV.resize (mat_size * mat_size);
	fill (VV.begin(), VV.end(), 0.);
	for (unsigned j1 = 0; j1 < i_neigh_index.size(); ++j1){
	  VV[j1*mat_size+j1] = 1.;
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
	    VV[j1*mat_size+j2] = dot (diff1, diff2) / r1 / r2;
	    // VV(j2, j1) = VV(j1, j2);
	  }
	}
      }
};


class assemble_trait_mat_dot {
public:
  void operator () (vector<double > & VV,
		    const vector<vector<double > > & waters,
		    const vector<double > & box,
		    const int & i_index,
		    const vector<int > & i_neigh_index)
      {
	unsigned mat_size = i_neigh_index.size();
	VV.resize (mat_size * mat_size);
	fill (VV.begin(), VV.end(), 0.);
	for (unsigned j1 = 0; j1 < i_neigh_index.size(); ++j1){
	  VV[j1*mat_size+j1] = 0.;
	}
  
	vector<double > diff1(3), diff2(3);
	for (unsigned j1 = 1; j1 < i_neigh_index.size(); ++j1){
	  int j1_index = i_neigh_index[j1];
	  diff (diff1, waters[j1_index], waters[i_index], box) ;    
	  for (unsigned j2 = 0; j2 < j1; ++j2){
	    int j2_index = i_neigh_index[j2];
	    diff (diff2, waters[j2_index], waters[i_index], box) ;    
	    VV[j1*mat_size+j2] = dot (diff1, diff2);
	    // VV(j2, j1) = VV(j1, j2);
	  }
	}
      }
};

