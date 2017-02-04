#pragma once

#include <vector> 
#include <cassert>

using namespace std;

class Displacement 
{
public:
  Displacement (const vector<double > & box);
  void reinit (const vector<double > & box);
  void deposite (const vector<vector<double > > & coms);
  const vector<double > & get_displ () const {return displacement;}
private:
  void set_first (const vector<vector<double > > & coms);
  void compute_shift (const vector<vector<double > > & prev_coms,
		      const vector<vector<double > > & curr_coms);
  void compute_displacement (const vector<vector<double > > & curr_coms);
  vector<vector<double > > first_frame;
  vector<vector<double > > prev_frame;
  vector<vector<int > > shifts;
  vector<double > displacement;
  unsigned nmol;
  int nframe;
  vector<double > box;
}
    ;

Displacement::
Displacement (const vector<double > & box)
{
  reinit (box);
}

void 
Displacement ::
reinit (const vector<double > & box_)
{
  nframe = -1;
  nmol = 0;
  box = box_;
}

void
Displacement::
set_first (const vector<vector<double > > & coms) 
{
  first_frame = coms;
  nmol = first_frame.size();
  prev_frame = first_frame;
  shifts.resize (nmol);
  for (unsigned ii = 0; ii < nmol; ++ii){
    shifts[ii] = vector<int> (3, 0);
  }
  displacement.resize (nmol, 0.);
  nframe ++;
}

void 
Displacement::
compute_shift (const vector<vector<double > > & prev_coms,
	       const vector<vector<double > > & curr_coms)
{
  assert (prev_coms.size() == curr_coms.size());
  assert (prev_coms.size() == shifts.size());
  assert (prev_coms.size() == nmol);
  
  for (unsigned ii = 0; ii < prev_coms.size(); ++ii){
    for (unsigned dd = 0; dd < 3; ++dd){
      double diff = curr_coms[ii][dd] - prev_coms[ii][dd];
      if      (diff >  0.5 * box[dd]) shifts[ii][dd] -= 1;
      else if (diff < -0.5 * box[dd]) shifts[ii][dd] += 1;
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
      double real_posi = curr_coms[ii][dd] + box[dd] * shifts[ii][dd];
      diff[dd] = real_posi - first_frame[ii][dd];
    }
    displacement[ii] = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
  }
}

void 
Displacement::
deposite (const vector<vector<double > > & coms)
{
  if (nframe == -1){
    set_first (coms);
  }
  else {
    compute_shift (prev_frame, coms);
    prev_frame = coms;
    compute_displacement (coms);
  }
  
}

