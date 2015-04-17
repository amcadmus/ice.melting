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

inline double
dist (const double & q4, const double & q6,
      const double ref_q4, const double ref_q6)
{
  return sqrt ( (q4 - ref_q4) * (q4 - ref_q4) + 
		(q6 - ref_q6) * (q6 - ref_q6) );
}

int main(int argc, char * argv[])
{
  std::string ifile, ofile;
  unsigned cq4, cq6;
  int min_hit;
  double binSize;
  double q4low(0), q4up(1);
  double q6low(0), q6up(1);
  double begin, end;
  
  po::options_description desc ("Calculates the average first hitting time in unit of T.\nAllow options");
  desc.add_options()
      ("help,h", "print this message")
      ("input", po::value<string > (&ifile)->default_value ("COLVAR"), "the input trajectory of q4 and q6.")
      ("begin", po::value<double > (&begin)->default_value (0), "the start time.")
      ("end", po::value<double > (&end)->default_value (0), "the end time.")
      ("column-q4", po::value<unsigned > (&cq4)->default_value (2), "the column of Q4.")
      ("column-q6", po::value<unsigned > (&cq6)->default_value (3), "the column of Q6.")
      ("min-hit,m", po::value<int > (&min_hit)->default_value (1), "minimum number of hitting.")
      ("bin-size,b", po::value<double > (&binSize)->default_value (0.01), "size of the bins.")
      ("output", po::value<string > (&ofile)->default_value ("centers.out"), "the output of selected centers.");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  unsigned nq4 = (q4up - q4low + 0.5 * binSize) / binSize;
  unsigned nq6 = (q6up - q6low + 0.5 * binSize) / binSize;
  
  vector<vector<bool > > active (nq4);
  vector<vector<int > > count_hit (nq4);
  vector<vector<double > > min_distance (nq4);
  vector<vector<double > > hit_time (nq4);
  for (unsigned ii = 0; ii < nq4; ++ii){
    active[ii].resize(nq6, false);
    count_hit[ii].resize(nq6, 0);
    min_distance[ii].resize(nq6, 100 * binSize);
    hit_time[ii].resize(nq6, -1);
  }
  
  char *line = NULL;
  size_t lineSize = 0;
  FILE * fp = fopen (ifile.c_str(), "r");
  cq4--;
  cq6--;	// make them c numbering 
  unsigned max46 = (cq4 > cq6) ? cq4 : cq6;
  unsigned count = 0;

  while (-1 != getline(&line, &lineSize, fp)){
    if (line[0] == '#') continue;
    if (count % 1000 == 0){
      printf ("# reading line %d      \r", count);
      fflush(stdout);
    }
    count ++;
    vector<string > words;
    StringOperation::split (string(line), words);
    if (words.size() < max46+1){
      cerr << "wrong line format at " << line << endl;
      exit (1);
    }
    double time = atof(words[0].c_str());
    if (time < begin) continue;
    if (end != 0 && time > end) break;

    double q4 = atof(words[cq4].c_str());
    double q6 = atof(words[cq6].c_str());
    unsigned iq4 = (q4) / binSize;
    unsigned iq6 = (q6) / binSize;
    if (iq4 >= nq4){
      cerr << "q4 " << q4 << " is out of range " << q4low << " " << q4up << endl;
      continue;
    }
    if (iq6 >= nq6){
      cerr << "q6 " << q6 << " is out of range " << q6low << " " << q6up << endl;
      continue;
    }

    active[iq4][iq6] = true;
    active[iq4+1][iq6] = true;
    active[iq4][iq6+1] = true;
    active[iq4+1][iq6+1] = true;
    count_hit[iq4][iq6] ++;
    count_hit[iq4+1][iq6] ++;
    count_hit[iq4][iq6+1] ++;
    count_hit[iq4+1][iq6+1] ++;
    double tmp_dist;
    double ref_q4, ref_q6;
    int pt_q4, pt_q6;
//////////////////////////////////////////
    pt_q4 = iq4;
    pt_q6 = iq6;
    ref_q4 = binSize * pt_q4;
    ref_q6 = binSize * pt_q6;
    tmp_dist = dist (q4, q6, ref_q4, ref_q6);
    if (tmp_dist < min_distance[pt_q4][pt_q6]) {
      min_distance[pt_q4][pt_q6] = tmp_dist;
      hit_time[pt_q4][pt_q6] = time;
    }
//////////////////////////////////////////
    pt_q4 = iq4+1;
    pt_q6 = iq6;
    ref_q4 = binSize * pt_q4;
    ref_q6 = binSize * pt_q6;
    tmp_dist = dist (q4, q6, ref_q4, ref_q6);
    if (tmp_dist < min_distance[pt_q4][pt_q6]) {
      min_distance[pt_q4][pt_q6] = tmp_dist;
      hit_time[pt_q4][pt_q6] = time;
    }
//////////////////////////////////////////
    pt_q4 = iq4;
    pt_q6 = iq6+1;
    ref_q4 = binSize * pt_q4;
    ref_q6 = binSize * pt_q6;
    tmp_dist = dist (q4, q6, ref_q4, ref_q6);
    if (tmp_dist < min_distance[pt_q4][pt_q6]) {
      min_distance[pt_q4][pt_q6] = tmp_dist;
      hit_time[pt_q4][pt_q6] = time;
    }
//////////////////////////////////////////
    pt_q4 = iq4+1;
    pt_q6 = iq6+1;
    ref_q4 = binSize * pt_q4;
    ref_q6 = binSize * pt_q6;
    tmp_dist = dist (q4, q6, ref_q4, ref_q6);
    if (tmp_dist < min_distance[pt_q4][pt_q6]) {
      min_distance[pt_q4][pt_q6] = tmp_dist;
      hit_time[pt_q4][pt_q6] = time;
    }
  }
  cout << endl;
  fclose (fp);
  
  cout << "# finish reading " << endl;
  
  fp = fopen (ofile.c_str(), "w");
  for (unsigned ii = 0; ii < nq4; ++ii){
    for (unsigned jj = 0; jj < nq6; ++jj){
      if (active[ii][jj] && count_hit[ii][jj] >= min_hit){
	fprintf (fp, "%f %f   %f %f\n", double(ii * binSize), double(jj * binSize), hit_time[ii][jj], min_distance[ii][jj]);
      }
    }
  }
  
  return 0;
}
