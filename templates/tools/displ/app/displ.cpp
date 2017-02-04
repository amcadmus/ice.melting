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
#include "Defines.h"
#include "Displacement.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

void align_water (matrix box,
		  vector<vector<double > > & water)
{
  vector<double > oh1(3), oh2(3);
  for (int dd = 0; dd < 3; ++dd){
    oh1[dd] = water[1][dd] - water[0][dd];
    oh2[dd] = water[2][dd] - water[0][dd];    
    double hbox = box[dd][dd] * 0.5;
    if      (oh1[dd] >  hbox) water[1][dd] -= box[dd][dd];
    else if (oh1[dd] < -hbox) water[1][dd] += box[dd][dd];
    if      (oh2[dd] >  hbox) water[2][dd] -= box[dd][dd];
    else if (oh2[dd] < -hbox) water[2][dd] += box[dd][dd];
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
close_mol_defect (const int numb_mol)
{
  // for (int ii = 0; ii < numb_mol; ++ii){
  //   fclose (fp[ii]);
  // }
  // free (fp);
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
  string ifile, ofile, odir;
  int func_numb_threads;
  int numb_mol_atom;
  bool p_detail (false), p_mol(false);
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b", po::value<double > (&begin)->default_value(0.f), "start time")
      ("end,e",   po::value<double > (&end  )->default_value(0.f), "end   time")
      ("detail", "print the displacement of each molecule at each step")
      ("mol-value", "print the displacement trajectory for each atom")
      ("numb-mol-atom", po::value<int > (&numb_mol_atom)->default_value(4), "number of sites in the water molecule")
      ("numb-threads,t", po::value<int > (&func_numb_threads)->default_value(1), "number of threads")
      ("input,f",   po::value<string > (&ifile)->default_value ("traj.xtc"), "the input .xtc file")
      ("output,o",  po::value<string > (&ofile)->default_value ("displ.out"), "the output file")
      ("output-dir",  po::value<string > (&odir)->default_value ("displ"), "the output directory of displacement value information");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    cout << desc<< "\n";
    return 0;
  }
  if (vm.count("detail")){
    p_detail = true;
  }
  if (vm.count("mol-value")){
    p_mol = true;
  }

  cout << "###################################################" << endl;
  cout << "# computes the displacement of molecules" << endl;
  cout << "# begin->end: " << begin << " " << end << endl;
  cout << "# numb sites in water: " << numb_mol_atom << endl;
  cout << "# input: " << ifile << endl;
  cout << "# output: " << ofile << endl;
  if (p_detail){
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
  vector<vector<vector<double > > > waters;
  waters.reserve (nmolecules * 3);

  Displacement displ (mybox);

  int countread = 0;
  if (p_detail || p_mol){
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
  if (p_mol) {
    open_mol_defect (odir, nmolecules);
  }
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
    if (((countread++)) % 100 == 0){
      printf ("# load frame at time: %.1f ps\r", time);
      fflush (stdout);
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

    displ.deposite (coms);

    // fprintf (fout, "%f\t %f\n", time, p_sa->getStepQ());
    if (p_detail){
      print_step (odir, step, time, displ.get_displ());
    }
    if (p_mol){
      print_mol (odir, step, time, displ.get_displ(), func_numb_threads);
    }
  }
  printf ("\n");

  FILE *fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cerr << "cannot open file " << ofile << endl;
    exit (1);
  }
  for (unsigned ii = 0; ii < displ.get_displ().size(); ++ii){
    fprintf (fout, "%06d %f\n", ii, displ.get_displ()[ii]);
  }
  fclose (fout);

  free (xx);  
  xdrfile_close (fp);

  return 0;
}
