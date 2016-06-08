#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <cmath>
#include <boost/program_options.hpp>
#include <vector>
#include "StringSplit.h"

#define MaxLineLength 65536

namespace po = boost::program_options;
using namespace std;

template <typename TYPE>
void 
insert (const vector<vector<TYPE> > & vold,
	vector<vector<TYPE> > & vnew,
	const int posii,
	const int posij)
{
  for (unsigned ii = 0; ii < vold.size(); ++ii){
    for (unsigned jj = 0; jj < vold[ii].size(); ++jj){
      vnew[posii + ii][posij + jj] = vold[ii][jj];
    }
  }
}

template <typename TYPE>
void 
insert (const vector<vector<vector<TYPE > > > & vold,
	vector<vector<vector<TYPE > > > & vnew,
	const int posii,
	const int posij,
	const int posik)
{
  for (unsigned ii = 0; ii < vold.size(); ++ii){
    for (unsigned jj = 0; jj < vold[ii].size(); ++jj){
      for (unsigned kk = 0; kk < vold[ii][jj].size(); ++kk){
	vnew[posii + ii][posij + jj][posik + kk] = vold[ii][jj][kk];
      }
    }
  }
}


inline double
dist (const double & c0, const double & c1, const double ref_c0, const double ref_c1)
{
  return sqrt ( (c0 - ref_c0) * (c0 - ref_c0) + 
		(c1 - ref_c1) * (c1 - ref_c1) );
}

inline double
dist (const double & c0, const double & c1, const double & c2, const double ref_c0, const double ref_c1, const double ref_c2)
{
  return sqrt ( (c0 - ref_c0) * (c0 - ref_c0) + 
		(c1 - ref_c1) * (c1 - ref_c1) + 
		(c2 - ref_c2) * (c2 - ref_c2) );
}

template <typename TYPE>
void
make_matrix (const unsigned nbin_c0, 
	     const unsigned nbin_c1, 
	     vector<vector<TYPE> > &  matrix, 
	     const TYPE default_value = 0)
{
  matrix.resize(nbin_c0+1);
  for (unsigned ii = 0; ii < matrix.size(); ++ii){
    matrix[ii].resize (nbin_c1 + 1, default_value);
  }
}

int 
compute_idx (const double value,
	     const double low,
	     const double bin_size)
{
  int shift (0);
  if (value - low < 0) shift = -1;
  return int((value - low) / bin_size) + shift;
}

void process (const vector<double > & time,
	      const vector<double > & c0,
	      const vector<double > & c1,
	      const double bin_size0,
	      const double bin_size1,
	      double & low_c0,
	      double & low_c1,
	      vector<vector<double > > & hit_time,
	      vector<vector<double > > & min_dist)
{
  if (c0.size() == 0) return;
  if (c0.size() != c1.size() || c0.size() != time.size()) {
    cerr << "unmatching c0 and c1 or c0 and time, exit" << endl;
    exit(1);
  }
  low_c0 = int(c0[0] / bin_size0) * bin_size0;
  low_c1 = int(c1[0] / bin_size1) * bin_size1;
  unsigned nbin_c0;
  unsigned nbin_c1;
  nbin_c0 = 1;
  nbin_c1 = 1;
  make_matrix<double> (nbin_c0, nbin_c1, hit_time, -1);
  make_matrix<double> (nbin_c0, nbin_c1, min_dist, 1e24);
  
  for (unsigned ii = 1; ii < c0.size(); ++ii){
    int idx0 = compute_idx (c0[ii], low_c0, bin_size0);
    int idx1 = compute_idx (c1[ii], low_c1, bin_size1);
    if (idx0 < 0 || idx1 < 0 || idx0 >= int(nbin_c0) || idx1 >= int(nbin_c1)){
      int new_idx0, new_idx1, new_nbin_c0, new_nbin_c1;
      int posi0, posi1;
      if (idx0 < 0){
	new_idx0 = 0;
	posi0 = -idx0;
	new_nbin_c0 = nbin_c0 + posi0;
	low_c0 = low_c0 + idx0 * bin_size0;
      }
      else if (idx0 >= int(nbin_c0)){
	new_idx0 = idx0;
	posi0 = 0;
	new_nbin_c0 = new_idx0 + 1;
      }
      else {
	new_idx0 = idx0;
	posi0 = 0;
	new_nbin_c0 = nbin_c0;
      }
      if (idx1 < 0){
	new_idx1 = 0;
	posi1 = -idx1;
	new_nbin_c1 = nbin_c1 + posi1;
	low_c1 = low_c1 + idx1 * bin_size1;
      }
      else if (idx1 >= int(nbin_c1)){
	new_idx1 = idx1;
	posi1 = 0;
	new_nbin_c1 = new_idx1 + 1;
      }
      else {
	new_idx1 = idx1;
	posi1 = 0;
	new_nbin_c1 = nbin_c1;
      }
      vector<vector<double > > new_hit_time;
      vector<vector<double > > new_min_dist;
      make_matrix<double> (new_nbin_c0, new_nbin_c1, new_hit_time, -1);
      insert (hit_time, new_hit_time, posi0, posi1);
      hit_time = new_hit_time;
      make_matrix<double> (new_nbin_c0, new_nbin_c1, new_min_dist, 1e12);
      insert (min_dist, new_min_dist, posi0, posi1);
      min_dist = new_min_dist;
      idx0 = new_idx0;
      idx1 = new_idx1;
      nbin_c0 = new_nbin_c0;
      nbin_c1 = new_nbin_c1;
    }
    double tmp_dist;
    int pt0, pt1;
    double ref_c0, ref_c1;
//////////////////////////////////////////////////
    pt0 = idx0;
    pt1 = idx1;
    ref_c0 = low_c0 + bin_size0 * pt0;
    ref_c1 = low_c1 + bin_size1 * pt1;
    tmp_dist = dist (c0[ii], c1[ii], ref_c0, ref_c1);
    if (tmp_dist < min_dist[pt0][pt1]){
      min_dist[pt0][pt1] = tmp_dist;
      hit_time[pt0][pt1] = time[ii];
    }
//////////////////////////////////////////////////
    pt0 = idx0;
    pt1 = idx1+1;
    ref_c0 = low_c0 + bin_size0 * pt0;
    ref_c1 = low_c1 + bin_size1 * pt1;
    tmp_dist = dist (c0[ii], c1[ii], ref_c0, ref_c1);
    if (tmp_dist < min_dist[pt0][pt1]){
      min_dist[pt0][pt1] = tmp_dist;
      hit_time[pt0][pt1] = time[ii];
    }
//////////////////////////////////////////////////
    pt0 = idx0+1;
    pt1 = idx1;
    ref_c0 = low_c0 + bin_size0 * pt0;
    ref_c1 = low_c1 + bin_size1 * pt1;
    tmp_dist = dist (c0[ii], c1[ii], ref_c0, ref_c1);
    if (tmp_dist < min_dist[pt0][pt1]){
      min_dist[pt0][pt1] = tmp_dist;
      hit_time[pt0][pt1] = time[ii];
    }
//////////////////////////////////////////////////
    pt0 = idx0+1;
    pt1 = idx1+1;
    ref_c0 = low_c0 + bin_size0 * pt0;
    ref_c1 = low_c1 + bin_size1 * pt1;
    tmp_dist = dist (c0[ii], c1[ii], ref_c0, ref_c1);
    if (tmp_dist < min_dist[pt0][pt1]){
      min_dist[pt0][pt1] = tmp_dist;
      hit_time[pt0][pt1] = time[ii];
    }
  }
}

int main(int argc, char * argv[])
{
  std::string ifile, ofile;
  string cols, bins;
  double begin, end;
  
  po::options_description desc ("Calculates the average first hitting time in unit of T.\nAllow options");
  desc.add_options()
      ("help,h", "print this message")
      ("input", po::value<string > (&ifile)->default_value ("COLVAR"), "the input trajectory of q4 and q6.")
      ("begin", po::value<double > (&begin)->default_value (0), "the start time.")
      ("end", po::value<double > (&end)->default_value (0), "the end time.")
      ("columns", po::value<string > (&cols)->default_value (""), "the columns of CV.")
      // ("min-hit,m", po::value<int > (&min_hit)->default_value (1), "minimum number of hitting.")
      ("bin-size,b", po::value<string > (&bins)->default_value (""), "sizes of the bins.")
      ("output", po::value<string > (&ofile)->default_value ("centers.out"), "the output of selected centers.");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  vector<string > words;
  vector<double > binSizes;
  vector<unsigned > columns;
  StringOperation::split (cols, words);
  for (unsigned ii = 0; ii < words.size(); ++ii){
    columns.push_back (atoi(words[ii].c_str()) - 1);
  }
  StringOperation::split (bins, words);
  for (unsigned ii = 0; ii < words.size(); ++ii){
    binSizes.push_back (atof(words[ii].c_str()));
  }
  if (columns.size() == 0) return 0;
  if (binSizes.size() == 1){
    for (unsigned ii = 1; ii < columns.size(); ++ii){
      binSizes.push_back (binSizes[0]);
    }
  }
  if (columns.size() != binSizes.size()){
    cerr << " the number of words for columns and binsizes does not match!, exit!"<< endl;
    return 1;
  }
  
  unsigned maxcol = 0;
  for (unsigned ii = 0; ii < columns.size(); ++ii){
    if (columns[ii] > maxcol) maxcol = columns[ii];
  }

  vector<vector<double > > data;
  data.resize (1 + columns.size());

  char * line = NULL;
  size_t lineSize = 0;
  FILE * fp = fopen (ifile.c_str(), "r");
  int count = 0;
  while (-1 != getline(&line, &lineSize, fp)){
    if (line[0] == '#') continue;
    if (count % 1000 == 0){
      printf ("# reading line %d      \r", count);
      fflush(stdout);
    }
    count ++;
    vector<string > words;
    StringOperation::split (string(line), words);
    if (words.size() < maxcol + 1){
      cerr << "wrong line format at " << line << endl;
      exit (1);
    }
    double time = atof(words[0].c_str());
    if (time < begin) continue;
    if (end != 0 && time > end) break;

    data[0].push_back (time);
    for (unsigned ii = 0; ii < columns.size(); ++ii){
      data[ii+1].push_back ( atof(words[columns[ii]].c_str()) );
    }
  }
  cout << endl;
  fclose (fp);  
  cout << "# finish reading " << endl;
  cout << "# start computing " << endl;

  if (data.size() - 1 == 2){
    double low_c0;
    double low_c1;
    vector<vector<double > > hit_time;
    vector<vector<double > > min_dist;
    process (data[0], 
	     data[1], data[2], 
	     binSizes[0], binSizes[1],
	     low_c0, low_c1,
	     hit_time, min_dist);
    fp = fopen (ofile.c_str(), "w");
    for (unsigned ii = 0; ii < hit_time.size(); ++ii){
      for (unsigned jj = 0; jj < hit_time[ii].size(); ++jj){
	if (hit_time[ii][jj] >= 0){
	  fprintf (fp, "%f %f   %f %f\n", 
		   low_c0 + double(ii * binSizes[0]), 
		   low_c1 + double(jj * binSizes[1]), 
		   hit_time[ii][jj], 
		   min_dist[ii][jj]);
	}
      }
    }
    fclose (fp);
  }
  else {
    cerr << "not implemented!" << endl;
    return 1;
  }
  
  
  return 0;
}
