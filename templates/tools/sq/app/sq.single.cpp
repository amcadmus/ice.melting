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
#include <vector>
#include <cmath>
#include <sys/stat.h>

#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
#include "StructureFactor.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

double
sq_single (const vector<double > qq,
	   const vector<vector<double > > & com)
{
  double sqr=0, sqi=0;
  for (unsigned ii = 0; ii < com.size(); ++ii){
    double tmp = 0;
    for (int dd = 0; dd < 3; ++dd) tmp += qq[dd] * com[ii][dd];
    sqr += cos(tmp);
    sqi += sin(tmp);
  }
  return sqrt(sqr * sqr + sqi * sqi) / double(com.size());
}

void
sq_single_part (vector<double > & part_value,
		const vector<double > qq,
		const vector<vector<double > > & com)
{
  part_value.resize (com.size());
  for (unsigned ii = 0; ii < com.size(); ++ii){
    double sqr=0, sqi=0;    
    cout << ii << endl;
    for (unsigned jj = 0; jj < com.size(); ++jj){
      double diff[3];
      for (int dd = 0; dd < 3; ++dd) diff[dd] = com[ii][dd] - com[jj][dd];
      double tmp = 0;
      for (int dd = 0; dd < 3; ++dd) tmp += qq[dd] * diff[dd];
      sqr += cos(tmp);
      sqi += sin(tmp);
    }
    part_value[ii] = sqrtf(sqr * sqr + sqi * sqi) / double(com.size());
    
    // cout << sqr << " " 
    // 	 << sqi << " " 
    // 	 << sqrt(sqr * sqr + sqi * sqi) / double(com.size()) << endl;
  }
}

void print_step (const string dir,
		 const int step, 
		 const double time, 
		 const vector<double > & step_value)
{
  char fname [1024];
  // sprintf (fname, "/step_%09d", step);
  sprintf (fname, "/step_%09d_time_%.3f", step, time);
  string fpath = dir + string(fname);
  FILE * fp = fopen (fpath.c_str(), "w");
  if (fp == NULL){
    cerr << "cannot open file " << fpath << endl;
    exit (1);
  }
  fprintf (fp, "# mol_idx mole_Q_value\n");
  for (unsigned ii = 0; ii < step_value.size(); ++ii){
    fprintf (fp, "%d %f\n", ii, step_value[ii]);
  }
  fclose (fp);
}

// FILE ** 
void
open_mol_defect (const string dir,
		 const int numb_mol)
{
  // FILE ** fp = (FILE ** ) malloc (sizeof(FILE *) * numb_mol);
  for (int ii = 0; ii < numb_mol; ++ii){
    FILE * fp;
    char name[1024];
    sprintf (name, "%s/mol_%06d", dir.c_str(), ii);
    fp = fopen (name, "w");
    if (fp == NULL){
      cerr << "cannot open file " << string(name) << endl;
      exit (1);
    }
    fclose (fp);
  }
  // return fp;
}

void 
print_mol (const string dir,
	   const int step, 
	   const double time, 
	   const vector<double > & step_value, 
	   const int func_numb_threads)
{
#pragma omp parallel for num_threads (func_numb_threads) 
  for (unsigned ii = 0; ii < step_value.size(); ++ii){
    char name[1024];
    sprintf (name, "%s/mol_%06d", dir.c_str(), ii);
    FILE * fp = fopen (name, "a");
    fprintf (fp, "%09d %.3f %f\n", step, time, step_value[ii]);
    fclose (fp);
  }
}

int main(int argc, char * argv[])
{
  double begin, end;
  vector<double> qq(3);
  string ifile, ofile, odir;
  int func_numb_threads;
  int numb_mol_atom;
  bool p_mol (false);
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b", po::value<double > (&begin)->default_value(0.f), "start time")
      ("end,e",   po::value<double > (&end  )->default_value(0.f), "end   time")
      ("qx", po::value<double > (&qq[0]), "value of q_x")
      ("qy", po::value<double > (&qq[1]), "value of q_y")
      ("qz", po::value<double > (&qq[2]), "value of q_z")
      ("mol-value", "print the Q trajectory for each atom is printed")
      ("numb-mol-atom", po::value<int > (&numb_mol_atom)->default_value(4), "number of sites in the water molecule")
      ("numb-threads,t", po::value<int > (&func_numb_threads)->default_value(1), "number of threads")
      ("input,f",   po::value<string > (&ifile)->default_value ("traj.xtc"), "the input .xtc file")
      ("output,o",  po::value<string > (&ofile)->default_value ("sqs.out"), "the output time average of mol sq")
      ("output-dir",  po::value<string > (&odir)->default_value ("sqs"), "the output directory of mol sq");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    cout << desc<< "\n";
    return 0;
  }
  if (vm.count("mol-value")){
    p_mol = true;
  }
  if ( (!vm.count("qx")) || 
       (!vm.count("qy")) ||
       (!vm.count("qz")) ){
    cerr << "qx, qy, qz should be provided, exit" << endl;
    return 1;
  }

  cout << "###################################################" << endl;
  cout << "# computes the displacement of molecules" << endl;
  cout << "# begin->end: " << begin << " " << end << endl;
  cout << "# numb sites in water: " << numb_mol_atom << endl;
  cout << "# qx: " << qq[0] << endl;
  cout << "# qy: " << qq[1] << endl;
  cout << "# qz: " << qq[2] << endl;
  cout << "# input: " << ifile << endl;
  cout << "# output: " << ofile << endl;
  if (p_mol){
    cout << "# output dir: " << odir << endl;
  }
  cout << "###################################################" << endl;  
  
  XDRFILE *fp;
  int natoms, step;
  float time;
  matrix box;
  rvec * xx;
  float prec = 1000;
  float time_prec = .01;

  char tmpfname[1024];
  strncpy (tmpfname, ifile.c_str(), 1023);
  {
    int c;
    if ((c = read_xtc_natoms (tmpfname, &natoms)) == 0) {
      xx = (rvec *) malloc (sizeof(rvec) * natoms);
    }
    else {
      fprintf (stderr, "error read_xtc_natoms\n");
      exit (1);
    }
  }

  fp = xdrfile_open (ifile.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << ifile << endl;
    exit (1);
  }
  if (read_xtc (fp, natoms, &step, &time, box, xx, &prec) != 0) {
    cerr << "error reading file " << ifile << endl;
    return 1;
  }
  xdrfile_close (fp);
  
  int nmolecules = 0;
  nmolecules = natoms / numb_mol_atom;

  vector<double > mybox(3);
  mybox[0] = box[0][0];
  mybox[1] = box[1][1];
  mybox[2] = box[2][2];  
  vector<vector<double > > coms;
  coms.reserve (nmolecules);
  vector<double > sum_mol_sq (nmolecules);

  int countread = 0;
  int count = 0;
  if (p_mol){
    if (access (odir.c_str(), 0) == -1) {
      cout << "# dir " << odir << " does not exist, create." << endl;
      if (mkdir (odir.c_str(), 0755)){
	cerr << "# dir " << odir << " creation failed." << endl;
	return 1;
      }
    }
  }
  
  fp = xdrfile_open (ifile.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << ifile << endl;
    exit (1);
  }
  // if (p_mol) {
  //   open_mol_defect (odir, nmolecules);
  // }
  while (read_xtc (fp, natoms, &step, &time, box, xx, &prec) == 0){
    if (end != 0.f) {
      if (time < begin - time_prec){
	continue;
      }
      else if (time > end + time_prec) {
	break;
      }	
    }
    else {
      if (time < begin - time_prec) continue;
    }
    if (((countread++)) % 10 == 0){
      // printf ("# load frame at time: %.1f ps\r", time);
      // fflush (stdout);
    }
    
    coms.clear ();
    int nmol = natoms / numb_mol_atom;    
    vector<double > read_com(3, 0.);

    for (int ii = 0; ii < nmol; ++ii){
      for (int dd = 0; dd < 3; ++dd){
	read_com[dd] = xx[ii*numb_mol_atom][dd];
	if      (read_com[dd] <  0          ) read_com[dd] += box[dd][dd];
	else if (read_com[dd] >= box[dd][dd]) read_com[dd] -= box[dd][dd];
      }
      coms.push_back (read_com);
    }
    vector<double >  mybox(3);
    for (unsigned dd = 0; dd < 3; ++dd) mybox[dd] = box[dd][dd];
    
    // vector<double > mol_sq;
    cout << time << "\t " << sq_single (qq, coms) << endl;;
    // if (p_mol){
    //   print_mol (odir, step, time, mol_sq, func_numb_threads);
    // }
    // for (unsigned ii = 0; ii < mol_sq.size(); ++ii){
    //   sum_mol_sq[ii] += mol_sq[ii];
    // }
    count += 1;
  }
  printf ("\n");


  // for (unsigned ii = 0; ii < sum_mol_sq.size(); ++ii){
  //   sum_mol_sq[ii] /= double(count);
  // }
  // FILE *fout = fopen (ofile.c_str(), "w");
  // if (fout == NULL){
  //   cerr << "cannot open file " << ofile << endl;
  //   exit (1);
  // }
  // for (unsigned ii = 0; ii < sum_mol_sq.size(); ++ii){
  //   if (ii == 0) continue;
  //   fprintf (fout, "%06d %f\n", ii, sum_mol_sq[ii]);
  // }
  // fclose (fout);

  free (xx);  
  xdrfile_close (fp);

  return 0;
}
