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
#include "BlockAverage.h"
#include "StringSplit.h"

#define MaxLineLength 65536

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string ifile;
  unsigned ccv;
  double ref_value, evl_value;
  double kappa;
  unsigned numbDataInBlock;
  static const double kb = 0.008314;
  double temperature;
  double begin, end;
  
  po::options_description desc ("Calculates the average first hitting time in unit of T.\nAllow options");
  desc.add_options()
      ("help,h", "print this message")
      ("input", po::value<string > (&ifile)->default_value ("COLVAR"), "the input trajectory of collective variable.")
      ("column-cv", po::value<unsigned > (&ccv)->default_value (2), "the column of collective variable.")
      ("begin", po::value<double > (&begin)->default_value (0), "the start time.")
      ("end", po::value<double > (&end)->default_value (0), "the end time.")
      ("ref-value,r", po::value<double > (&ref_value)->default_value (0), "the reference value of collective variable.")
      ("evl-value,e", po::value<double > (&evl_value)->default_value (0), "the value where the force is evaluated.")
      ("kappa,k", po::value<double > (&kappa)->default_value (1), "the spring constant of the restraint simulation.")
      ("temperature,t", po::value<double > (&temperature)->default_value (300), "the temperature of the system, needed for reweighting.")
      ("numb-data-in-block,b", po::value<unsigned > (&numbDataInBlock)->default_value (100), "the number of data in each block to evaluate the error.");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }
  
  double kbT = temperature * kb ;
  double beta = 1./kbT;
  char *line = NULL;
  size_t lineSize = 0;
  FILE * fp = fopen (ifile.c_str(), "r");
  ccv--;
  unsigned count = 0;
  vector<double > force, vij;
  double vij_min = 0;
  
  while (-1 != getline(&line, &lineSize, fp)){
    if (line[0] == '#') continue;
    if (count % 1000 == 0){
      printf ("# reading line %d      \r", count);
      fflush(stdout);
    }
    count ++;
    vector<string > words;
    StringOperation::split (string(line), words);
    if (words.size() < ccv+1 || words.size() < 1){
      cerr << "wrong line format at " << line << endl;
      exit (1);
    }
    double time = atof(words[0].c_str());
    double cv = atof(words[ccv].c_str());
    if (time < begin) continue;
    if (end != 0 && time > end) break;
    double force_esti = - kappa * (ref_value - cv);
    double vij_esti = 0.5 * (+kappa * (evl_value - cv) * (evl_value - cv)
			     -kappa * (ref_value - cv) * (ref_value - cv) );
    force.push_back (force_esti);
    vij.push_back (vij_esti);
    if (vij_esti < vij_min) vij_min = vij_esti;
  }
  cout << endl;
  fclose (fp);
  cout << "# finish reading " << endl;

  cout << "# post process the vij " << endl;
  for (unsigned ii = 0; ii < vij.size(); ++ii){
    vij[ii] -= vij_min;
  }
  BlockAverage_acc esti_force (numbDataInBlock);
  BlockAverage_acc esti_partition (numbDataInBlock);
  for (unsigned ii = 0; ii < vij.size(); ++ii){
    // esti_force.deposite (force[ii] * (1 + (-beta * (vij[ii] - vij_min))));
    // esti_partition.deposite (1 + (-beta * (vij[ii] - vij_min)));
    esti_force.deposite (force[ii] * exp(-beta * (vij[ii] - vij_min)));
    esti_partition.deposite (exp(-beta * (vij[ii] - vij_min)));
  }  
  esti_force.calculate();
  esti_partition.calculate();

  double v0 = esti_force.getAvg();
  double v1 = esti_partition.getAvg();
  double e0 = esti_force.getAvgError();
  double e1 = esti_partition.getAvgError();
  double ee = sqrt ( e0 / v1 * e0 / v1 +
		     e1 * v0 / v1 / v1 * e1 * v0 / v1 / v1);

  cout << "# force \t error" << endl;
  cout << v0/v1 << "\t" << ee << endl;
  
  return 0;
}
