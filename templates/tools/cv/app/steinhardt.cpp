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
#include "CellList.h"
#include "SteinhardtAnalysis.h"

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
  double begin, end, cellSize;
  string ifile, ofile, odir;
  int func_numb_threads;
  int numb_mol_atom;
  int order;
  double rmin, rmax;
  bool p_detail (false), p_mol(false), do_avg(false);
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b", po::value<double > (&begin)->default_value(0.f), "start time")
      ("end,e",   po::value<double > (&end  )->default_value(0.f), "end   time")
      ("r-min",   po::value<double > (&rmin)->default_value(0.31), "the min for r switch function")
      ("r-max",   po::value<double > (&rmax)->default_value(0.36), "the max for r switch function")
      ("order,l",   po::value<int > (&order)->default_value(6), "the order of the Steinhardt parameter")      
      ("loc-avg,a", "the local averaged Steinhardt parameter due to Leichner and Dellago")      
      ("detail", "print the Q value of each molecule at each step")
      ("mol-value", "print the Q trajectory for each atom is printed")
      ("numb-mol-atom", po::value<int > (&numb_mol_atom)->default_value(4), "number of sites in the water molecule")
      ("numb-threads,t", po::value<int > (&func_numb_threads)->default_value(1), "number of threads")
      ("input,f",   po::value<string > (&ifile)->default_value ("traj.xtc"), "the input .xtc file")
      ("output,o",  po::value<string > (&ofile)->default_value ("steinhardt.out"), "the output file")
      ("output-dir",  po::value<string > (&odir)->default_value ("steinhardt"), "the output directory of Q value information");
  
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
  if (vm.count("loc-avg")){
    do_avg = true;
  }

  double rcut = rmax;
  cellSize = rcut + 1e-6;
  cout << "###################################################" << endl;
  cout << "# computes the H-bond" << endl;
  cout << "# begin->end: " << begin << " " << end << endl;
  cout << "# numb sites in water: " << numb_mol_atom << endl;
  cout << "# input: " << ifile << endl;
  cout << "# output: " << ofile << endl;
  if (p_detail){
    cout << "# output dir: " << odir << endl;
  }
  cout << "# rmin: " << rmin << endl;
  cout << "# rmax: " << rmax << endl;
  cout << "# order: " << order << endl;
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

  VectorType vbox;
  vbox.x = box[0][0];
  vbox.y = box[1][1];
  vbox.z = box[2][2];
  if (cellSize >= .5 * vbox.x){
    cerr << "the cell size should be less than half of the box size" << endl;
    return 1;
  }
  
  vector<vector<double > > coms;
  coms.reserve (nmolecules);
  vector<vector<vector<double > > > waters;
  waters.reserve (nmolecules * 3);

  CellList clist (nmolecules, vbox, cellSize);
  SteinhardtAnalysisBase * p_sa;
  switch (order){
  case 3: 
      p_sa = new SteinhardtAnalysis<3> (rmin, rmax, func_numb_threads);
      break;
  case 4:
      p_sa = new SteinhardtAnalysis<4> (rmin, rmax, func_numb_threads);
      break;
  case 6:
      p_sa = new SteinhardtAnalysis<6> (rmin, rmax, func_numb_threads);
      break;
  default:
      cerr << " unsupported Steinhardt order " << order << " exit " << endl;
      return 1;
  }

  int countread = 0;
  FILE *fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cerr << "cannot open file " << ofile << endl;
    exit (1);
  }
  fprintf (fout, "# time  tot_numb_donator  tot_numb_acceptor\n");
  if (p_detail){
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
    if (((countread++)) % 10 == 0){
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

    vbox.x = box[0][0];
    vbox.y = box[1][1];
    vbox.z = box[2][2];
    clist.reinit (nmolecules, vbox, cellSize);
    clist.rebuild (coms);
    vector<double > vect_box(3);
    for (int dd = 0; dd < 3; ++dd) vect_box[dd] = box[dd][dd];

    p_sa->deposite (clist, vect_box, coms, do_avg);

    // fprintf (fout, "%f\t %f\n", time, p_sa->getStepQ());
    if (p_detail){
      print_step (odir, step, time, p_sa->getStepMole());
    }
    if (p_mol){
      print_mol (odir, step, time, p_sa->getStepMole(), func_numb_threads);
    }
  }
  printf ("\n");

  p_sa->average();
  vector<double > avg = p_sa->getAvgMole();
  for (unsigned ii = 0; ii < avg.size(); ++ii){
    fprintf (fout, "%06d %f\n", ii, avg[ii]);
  }
  
  free (xx);  
  fclose (fout);
  xdrfile_close (fp);
  delete p_sa;

  return 0;
}
