#pragma once

class RingStats 
{
public:
  RingStats (const int nmols, 
	     const int max_ring = 10);
  // void mol_compute (vector<vector<int > > & count,
  // 		    const vector<vector<vector<int > > > & ring_list);
  void sys_deposite (const vector<vector<int > > & ring_list);
  const vector<double > & get_frame_histo () const {return frame_value;}
  void print_head  (ofstream & fout) const;
  void print_frame (ofstream & fout, const double time) const;
private:
  vector<double > frame_value;
  vector<double > avg_value;
  int nmols;
  int max_ring;
  int nframes;
  int word_size;
};

RingStats::
RingStats (const int nmols_, const int max_ring_)
    : frame_value (1+max_ring_, 0.),
      avg_value (1+max_ring_, 0.),
      nmols (nmols_),
      max_ring (max_ring_), 
      nframes (0),
      word_size (12)
{
}

void 
RingStats::
sys_deposite (const vector<vector<int > > & ring_list)
{
  fill (frame_value.begin(), frame_value.end(), 0);
  for (unsigned ii = 0; ii < ring_list.size(); ++ii){
    frame_value[ring_list[ii].size()] += 1.;
  }
  for (unsigned ii = 0; ii < frame_value.size(); ++ii){
    frame_value[ii] /= double (nmols);
  }
}

void 
RingStats::
print_head  (ofstream & fout) const
{
  fout << "#       time" ;
  for (unsigned ii = 3; ii < frame_value.size(); ++ii){
    fout << setw (word_size) << ii ;
  }
  fout << endl;
}

void 
RingStats::
print_frame  (ofstream & fout, const double time) const
{
  fout.setf(ios::fixed, ios::floatfield);
  fout.setf(ios::showpoint);
  fout.precision (2);
  fout << setw(12) << time;
  fout.precision (4);
  fout << scientific;
  for (unsigned ii = 3; ii < frame_value.size(); ++ii){
    fout << setw (word_size) << frame_value[ii] ;
  }
  fout << endl;
}



