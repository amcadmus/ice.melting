#pragma once

#include <vector>

using namespace std;

class Clusters 
{
public:
  Clusters (const double rcut);
  void reinit (const double rcut);
  void analysis (const vector<vector<double > > & atoms, 
		 const vector<double > & box);
  const vector<vector<int > > & get_cluster () const {return clstr_list;}
private: 
  void clear();
  void init_null (const int natoms);
private:
  void merge_clstr (const int i0, 
		    const int i1);
  void add_clstr (const int aidx);
  void add_to_clstr (const int cidx, 
		     const int aidx);
  const int & get_cluster_idx (const int aidx) const;
  const vector<int > & get_atoms (const int cidx) const ;
public:
  double rcut, rcut2;
  vector<vector<int > > clstr_list;
  vector<int > atom_clstr_idx;
}
    ;
