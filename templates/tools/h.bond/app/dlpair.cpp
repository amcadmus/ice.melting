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
#include "HydrogenBond.h"
#include "HydrogenBondAnalysis.h"

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

void print_step_d (const string dir,
		   const int step, 
		   const double time, 
		   const vector<pair<int,int> > & dpair)
{
  char fname [1024];
  // sprintf (fname, "/step_%09d", step);
  sprintf (fname, "/dpair_%09d_time_%.3f", step, time);
  string fpath = dir + string(fname);
  FILE * fp = fopen (fpath.c_str(), "w");
  if (fp == NULL){
    cerr << "cannot open file " << fpath << endl;
    exit (1);
  }
  fprintf (fp, "# mol_idx numb_donator numb_acceptor\n");
  for (unsigned ii = 0; ii < dpair.size(); ++ii){
    fprintf (fp, "%d %d \n", dpair[ii].first, dpair[ii].second);
  }
  fclose (fp);
}

void print_step_l (const string dir,
		   const int step, 
		   const double time, 
		   const vector<pair<int,int> > & dpair)
{
  char fname [1024];
  // sprintf (fname, "/step_%09d", step);
  sprintf (fname, "/lpair_%09d_time_%.3f", step, time);
  string fpath = dir + string(fname);
  FILE * fp = fopen (fpath.c_str(), "w");
  if (fp == NULL){
    cerr << "cannot open file " << fpath << endl;
    exit (1);
  }
  fprintf (fp, "# mol_idx numb_donator numb_acceptor\n");
  for (unsigned ii = 0; ii < dpair.size(); ++ii){
    fprintf (fp, "%d %d \n", dpair[ii].first, dpair[ii].second);
  }
  fclose (fp);
}

void print_step (const string dir,
		 const int step, 
		 const double time, 
		 const vector<pair<int,int> > & dpair,
		 const vector<pair<int,int> > & lpair)
{
  print_step_d (dir, step, time, dpair);
  print_step_l (dir, step, time, lpair);
}


int main(int argc, char * argv[])
{
  double begin, end, cellSize;
  string ifile, ofile, odir;
  int func_numb_threads;
  int numb_mol_atom;
  double rcut, acut;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b", po::value<double > (&begin)->default_value(0.f), "start time")
      ("end,e",   po::value<double > (&end  )->default_value(0.f), "end   time")
      ("r-cut,r",   po::value<double > (&rcut)->default_value(0.35), "the cut-off of O-O dist")
      ("angle-cut,a",   po::value<double > (&acut)->default_value(30), "the cut-off of O-H .. O angle")
      ("numb-mol-atom", po::value<int > (&numb_mol_atom)->default_value(4), "number of sites in the water molecule")
      ("numb-threads,t", po::value<int > (&func_numb_threads)->default_value(1), "number of threads")
      ("input,f",   po::value<string > (&ifile)->default_value ("traj.xtc"), "the input .xtc file")
      ("output,o",  po::value<string > (&ofile)->default_value ("dlpair.out"), "the output file")
      ("output-dir",  po::value<string > (&odir)->default_value ("dlpair"), "the output directory of H-bond detailed information");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    cout << desc<< "\n";
    return 0;
  }

  cellSize = rcut + 1e-6;
  cout << "###################################################" << endl;
  cout << "# computes the H-bond" << endl;
  cout << "# begin->end: " << begin << " " << end << endl;
  cout << "# numb sites in water: " << numb_mol_atom << endl;
  cout << "# input: " << ifile << endl;
  cout << "# output: " << ofile << endl;
  cout << "# output dir: " << odir << endl;
  cout << "# rcut: " << rcut << endl;
  cout << "# acut: " << acut << endl;
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
  HydrogenBond_Geo_1::Parameters param (rcut, acut);
  HydrogenBondAnalysis<HydrogenBond_Geo_1> hba (param, nmolecules, func_numb_threads);

  int countread = 0;
  FILE *fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cerr << "cannot open file " << ofile << endl;
    exit (1);
  }
  fprintf (fout, "# time  tot_numb_donator  tot_numb_acceptor\n");
  if (access (odir.c_str(), 0) == -1) {
    cout << "# dir " << odir << " does not exist, create." << endl;
    if (mkdir (odir.c_str(), 0755)){
      cerr << "# dir " << odir << " creation failed." << endl;
      return 1;
    }
  }
  
  fp = xdrfile_open (ifile.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << ifile << endl;
    exit (1);
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
    waters.clear ();
    int nmol = natoms / numb_mol_atom;    
    vector<double > read_com(3, 0.);
    vector<vector<double > > read_water(3);
    for (int ii = 0; ii < 3; ++ii) read_water[ii].resize(3, 0.);

    for (int ii = 0; ii < nmol; ++ii){
      for (int dd = 0; dd < 3; ++dd){
	read_com[dd] = xx[ii*numb_mol_atom][dd];
	if      (read_com[dd] <  0          ) read_com[dd] += box[dd][dd];
	else if (read_com[dd] >= box[dd][dd]) read_com[dd] -= box[dd][dd];
      }
      for (int jj = 0; jj < 3; ++jj){
	for (int dd = 0; dd < 3; ++dd){
	  read_water[jj][dd] = xx[ii*numb_mol_atom + jj][dd];
	}
      }
      for (int dd = 0; dd < 3; ++dd){
	if      (read_water[0][dd] <  0          ) read_water[0][dd] += box[dd][dd];
	else if (read_water[0][dd] >= box[dd][dd]) read_water[0][dd] -= box[dd][dd];
      }
      align_water (box, read_water);
      coms.push_back (read_com);
      waters.push_back (read_water);
    }

    vbox.x = box[0][0];
    vbox.y = box[1][1];
    vbox.z = box[2][2];
    clist.reinit (nmolecules, vbox, cellSize);
    clist.rebuild (coms);
    vector<double > vect_box(3);
    for (int dd = 0; dd < 3; ++dd) vect_box[dd] = box[dd][dd];

    vector<pair<int, int> > dpair, lpair;
    hba.findDLPair (dpair, lpair, clist, vect_box, waters);


    fprintf (fout, "%f\t %d \t %d\n", time, int(dpair.size()), int(lpair.size()));
    print_step (odir, step, time, dpair, lpair);
  }
  printf ("\n");
  
  free (xx);  
  fclose (fout);
  xdrfile_close (fp);

  return 0;
}
