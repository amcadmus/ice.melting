#pragma once

#include <vector> 
#include <cassert>

using namespace std;

class Displacement 
{
public:
  Displacement ();
  void reinit ();
  void deposite (const vector<vector<double > > & coms, 
		 const vector<double > & box);
  const vector<double > & get_displ () const {return displacement;}
private:
  void set_first (const vector<vector<double > > & coms);
  void remove_pbc (vector<vector<double > > & curr_coms,
		   const vector<vector<double > > & prev_coms,
		   const vector<double > & box);
  void remove_sys_shift (vector<vector<double > > & new_coms,
			 const vector<vector<double > > & old_coms);
  void compute_center (vector<double > & center,
		       const vector<vector<double > > & coms);
  void compute_displacement (const vector<vector<double > > & curr_coms);
  vector<vector<double > > first_frame;
  vector<vector<double > > prev_frame;
  vector<double > displacement;
  vector<double > center0;
  unsigned nmol;
  int nframe;
}
    ;

Displacement::
Displacement ()
{
  reinit ();
}

void 
Displacement ::
reinit ()
{
  nframe = -1;
  nmol = 0;
}

void
Displacement::
set_first (const vector<vector<double > > & coms) 
{
  first_frame = coms;
  compute_center (center0, first_frame);
  prev_frame = first_frame;
  nmol = first_frame.size();
  displacement.resize (nmol, 0.);
  fill (displacement.begin(), displacement.end(), 0.);
  nframe ++;
}


void 
Displacement::
remove_pbc (vector<vector<double > > & curr_coms,
	    const vector<vector<double > > & prev_coms,
	    const vector<double > & box)
{
  assert (prev_coms.size() == curr_coms.size());
  assert (prev_coms.size() == nmol);
  
  for (unsigned ii = 0; ii < prev_coms.size(); ++ii){
    for (unsigned dd = 0; dd < 3; ++dd){
      while (curr_coms[ii][dd] - prev_coms[ii][dd] >   0.5 * box[dd]) {
	curr_coms[ii][dd] -= box[dd];
      }
      while (curr_coms[ii][dd] - prev_coms[ii][dd] <= -0.5 * box[dd]) {
	curr_coms[ii][dd] += box[dd];
      }
    }
  }
}

void 
Displacement::
compute_center (vector<double > & center,
		const vector<vector<double > > & coms)
{
  center.resize(3);
  for (int dd = 0; dd < 3; ++dd) center[dd] = 0.;
  for (unsigned ii = 0; ii < coms.size(); ++ii){
    for (int dd = 0; dd < 3; ++dd) {
      center[dd] += coms[ii][dd];
    }
  }
  for (int dd = 0; dd < 3; ++dd) center[dd] /= double(coms.size());
  // cout << center[0] << " " << center[1] << " " << center[2] << endl;
}

void 
Displacement::
remove_sys_shift (vector<vector<double > > & new_coms,
		  const vector<vector<double > > & old_coms)
{
  assert (new_coms.size() == old_coms.size());
  vector<double > center(3);
  compute_center (center, old_coms);
  vector<double > shift(3);  
  for (int dd = 0; dd < 3; ++dd) shift[dd] = center[dd] - center0[dd];
  for (unsigned ii = 0; ii < old_coms.size(); ++ii){
    for (int dd = 0; dd < 3; ++dd) {
      new_coms[ii][dd] = old_coms[ii][dd] - shift[dd];
    }
  }
}

void
Displacement::
compute_displacement (const vector<vector<double > > & curr_coms)
{
  assert (curr_coms.size() == nmol);

  for (unsigned ii = 0; ii < curr_coms.size(); ++ii){
    vector<double> diff(3);
    for (unsigned dd = 0; dd < 3; ++dd){
      double real_posi = curr_coms[ii][dd];
      diff[dd] = real_posi - first_frame[ii][dd];
    }
    displacement[ii] = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
  }
}

void 
Displacement::
deposite (const vector<vector<double > > & coms_, 
	  const vector<double > & box)
{
  if (nframe == -1){
    set_first (coms_);
  }
  else {
    vector<vector<double > > curr_coms (coms_);
    vector<vector<double > > new_curr_coms (coms_);
    remove_pbc (curr_coms, prev_frame, box);
    remove_sys_shift (new_curr_coms, curr_coms);
    prev_frame = curr_coms;
    compute_displacement (new_curr_coms);
  }
}

