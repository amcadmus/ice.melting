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
  double ref_value;
  double kappa;
  unsigned numbDataInBlock;
  
  po::options_description desc ("Calculates the average first hitting time in unit of T.\nAllow options");
  desc.add_options()
      ("help,h", "print this message")
      ("input", po::value<string > (&ifile)->default_value ("COLVAR"), "the input trajectory of collective variable.")
      ("column-cv", po::value<unsigned > (&ccv)->default_value (2), "the column of collective variable.")
      ("ref-value,r", po::value<double > (&ref_value)->default_value (0), "the reference value of collective variable.")
      ("kappa,k", po::value<double > (&kappa)->default_value (1), "the spring constant of the restraint simulation.")
      ("numb-data-in-block,b", po::value<unsigned > (&numbDataInBlock)->default_value (100), "the number of data in each block to evaluate the error.");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }
  
  char *line = NULL;
  size_t lineSize = 0;
  FILE * fp = fopen (ifile.c_str(), "r");
  ccv--;
  unsigned count = 0;
  BlockAverage_acc naive_force (numbDataInBlock);

  while (-1 != getline(&line, &lineSize, fp)){
    if (line[0] == '#') continue;
    if (count % 1000 == 0){
      printf ("# reading line %d      \r", count);
      fflush(stdout);
    }
    vector<string > words;
    StringOperation::split (string(line), words);
    if (words.size() < ccv+1){
      cerr << "wrong line format at " << line << endl;
      exit (1);
    }
    double cv = atof(words[ccv].c_str());
    double force_esti = - kappa * (ref_value - cv);
    naive_force.deposite (force_esti);
    count ++;
  }
  cout << endl;
  fclose (fp);
  
  cout << "# finish reading " << endl;

  naive_force.calculate();
  cout << "# force \t error \n" << endl;
  cout << naive_force.getAvg() << "\t" << naive_force.getAvgError() << endl;
  
  return 0;
}
