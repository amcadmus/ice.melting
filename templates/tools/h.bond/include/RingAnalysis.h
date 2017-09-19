#pragma once

#include <iostream>
#include <vector>
using namespace std;

class RingAnalysis 
{
public:
  void compute (vector<vector<vector<int > > > & ring_list,
		const vector<vector<int > > & hbond_list, 
		const int depth = 10);
private:
  void compute (vector<vector<int > > & ring_list,
		const int i_idx,
		const int j_idx,
		const int pj_idx,
		const vector<vector<int > > & hbond_list, 
		const int depth);
  void cut (vector<vector<int > > & ring_list);
  void remove_dupl (vector<vector<int > > & ring_list);
}
    ;


void
RingAnalysis::
compute (vector<vector<vector<int > > > & ring_list,
	 const vector<vector<int > > & hbond_list, 
	 const int depth)
{
  int natoms = hbond_list.size();
  ring_list.clear();
  ring_list.resize (natoms);
  
  for (int ii = 0; ii < natoms; ++ii) {
    // if (ii % 10 == 0) cout << ii << " \r " ;
    // cout.flush();
    compute (ring_list[ii], 
	     ii,
	     ii,
	     -1,
	     hbond_list,
	     depth);
  }
  // cout << endl;

  for (unsigned ii = 0; ii < ring_list.size(); ++ii){
    cut (ring_list[ii]);
    remove_dupl (ring_list[ii]);
  }  
}


void
RingAnalysis::
compute (vector<vector<int > > & ring_list,
	 const int i_idx,
	 const int j_idx, 
	 const int parent_j_idx,	 
	 const vector<vector<int > > & hbond_list,
	 const int depth)
{
  if (depth == 1){
    ring_list.clear();
    const vector<int > & bd_list (hbond_list[j_idx]);
    if (find (bd_list.begin(), bd_list.end(), i_idx) != bd_list.end()){
      ring_list.push_back (vector<int> (1, j_idx));
    }
  }
  else {
    ring_list.clear();
    for (unsigned ii = 0; ii < hbond_list[j_idx].size(); ++ii){
      const int & tmp_k_idx = hbond_list[j_idx][ii];
      if (parent_j_idx >= 0 && parent_j_idx == tmp_k_idx) continue;
      if (tmp_k_idx == i_idx){
	ring_list.push_back (vector<int> (1, j_idx));
      }
      else {
	vector<vector<int > > tmp_ring_list;
	compute (tmp_ring_list,
		 i_idx,
		 tmp_k_idx,
		 j_idx,
		 // -1,
		 hbond_list,
		 depth - 1);
	for (unsigned jj = 0; jj < tmp_ring_list.size(); ++jj){
	  tmp_ring_list[jj].insert (tmp_ring_list[jj].begin(), j_idx);
	}
	ring_list.insert (ring_list.end(), tmp_ring_list.begin(), tmp_ring_list.end());
      }
    }
  }
}



static
bool
check_tag (const vector<int > & tag)
{
  for (unsigned ii = 0; ii < tag.size(); ++ii){
    if (tag[ii] == 0) return true;
  }
  return false;
}

static
bool
same (const vector<int > & r0, 
      const vector<int > & r1)
{
  assert (r0.size() <= r1.size());
  int count_b = 0;
  while (r0[count_b] == r1[count_b]) {
    count_b ++;
    if (count_b >= int(r0.size())) break;
  }
  int count_e = 0;
  int r0_back = r0.size() - 1;
  int r1_back = r1.size() - 1;
  while (r0[r0_back - count_e] == r1[r1_back - count_e]) {
    count_e ++;
    if (r0_back - count_e < 0) break;
  }
  if (count_e + count_b >= 3) return true;
  return false;
}

static
bool
same_bidir (const vector<int > & r0, 
	    const vector<int > & r1) 
{
  if (same (r0, r1)) {
    return true;
  }
  else {
    vector<int > r1_ (r1);
    reverse (r1_.begin() + 1, r1_.end());
    if (same(r0, r1_)) return true;
  }
  return false;
}

static
bool
identical (const vector<int > & r0, 
	   const vector<int > & r1)
{
  vector<int > r1_ (r1);
  reverse (r1_.begin() + 1, r1_.end());
  return (r0 == r1) || (r0 == r1_);
}


void
RingAnalysis::
cut (vector<vector<int > > & ring_list)
{
  vector<vector<int > > orig (ring_list);
  ring_list.clear();
  vector<int > tag (orig.size(), 0);
  
  while (check_tag (tag)){
    vector<pair<int, int> > sorting;
    for (unsigned ii = 0; ii < tag.size(); ++ii){
      if (tag[ii] == 0){
	sorting.push_back (pair<int, int> (orig[ii].size(), ii));
      }
    }
    sort (sorting.begin(), sorting.end());
    int ref_idx = sorting[0].second;
    tag[ref_idx] = 1;
    ring_list.push_back (orig[ref_idx]);
    for (unsigned ii = 0; ii < tag.size(); ++ii){
      if (tag[ii] == 0){
	if (same_bidir (orig[ref_idx], orig[ii])) {
	  assert (orig[ii].size() >= orig[ref_idx].size());
	  // if (orig[ii].size() == orig[ref_idx].size()){
	  //   ring_list.push_back (orig[ii]);
	  // }
	  // tag[ii] = 1;
	  if (orig[ii].size() > orig[ref_idx].size()){
	    tag[ii] = 1;
	  }
	}
	if (identical (orig[ref_idx], orig[ii])) {
	  tag[ii] = 1;
	}
      }
    }
  }
}

static 
bool
is_dupl (const vector<int> & rlist)
{
  vector<int> tmp (rlist);
  sort (tmp.begin(), tmp.end());
  vector<int>::iterator iter = unique (tmp.begin(), tmp.end());
  return iter != tmp.end();
}

void
RingAnalysis::
remove_dupl (vector<vector<int > > & ring_list)
{
  bool remove = false;
  vector<int > tag (ring_list.size());
  for (unsigned ii = 0; ii < ring_list.size(); ++ii){
    if (is_dupl(ring_list[ii])){
	// for (unsigned kk = 0; kk < ring_list[ii].size(); ++kk){
	//   cout << ring_list[ii][kk] << " ";
	// }
	// cout << endl;
	tag[ii] = 1;
	remove = true;
    }
  }
  if (remove == true) {
    vector<vector<int > > orig (ring_list);
    ring_list.clear();
    for (unsigned ii = 0; ii < tag.size(); ++ii){
      if (tag[ii] == 0) ring_list.push_back (orig[ii]);
    }
    // cout << orig.size() << "  ->  " << ring_list.size() << endl;
  }
}

