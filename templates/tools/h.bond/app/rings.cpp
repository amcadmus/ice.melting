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
#include "RingAnalysis.h"
#include "RingStats.h"

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
		 const vector<vector<int> > & urlist)
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
  for (unsigned ii = 0; ii < urlist.size(); ++ii){
    if (urlist[ii].size() != 6){
      for (unsigned jj = 0; jj < urlist[ii].size(); ++jj){
	fprintf (fp, "%d ", urlist[ii][jj]);
      }
      fprintf (fp, "\n");
    }
  }
  fclose (fp);
  
}


int main(int argc, char * argv[])
{
  double begin, end, cellSize;
  string ifile, ofile, odir;
  int func_numb_threads;
  int numb_mol_atom;
  int max_ring;
  double rcut, acut;
  bool p_detail (false), p_mol(false);
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b", po::value<double > (&begin)->default_value(0.f), "start time")
      ("end,e",   po::value<double > (&end  )->default_value(0.f), "end   time")
      ("r-cut,r",   po::value<double > (&rcut)->default_value(0.35), "the cut-off of O-O dist")
      ("angle-cut,a",   po::value<double > (&acut)->default_value(30), "the cut-off of O-H .. O angle")
      ("max-ring,m",   po::value<int > (&max_ring)->default_value(10), "the maximum size of rings")      
      ("detail", "print the numb of H-bond of each molecule at each step")
      ("mol-defect", "print if the molecule is in defect status. the history for each atom is printed")
      ("numb-mol-atom", po::value<int > (&numb_mol_atom)->default_value(4), "number of sites in the water molecule")
      ("numb-threads,t", po::value<int > (&func_numb_threads)->default_value(1), "number of threads")
      ("input,f",   po::value<string > (&ifile)->default_value ("traj.xtc"), "the input .xtc file")
      ("output,o",  po::value<string > (&ofile)->default_value ("hbring.out"), "the output file")
      ("output-dir",  po::value<string > (&odir)->default_value ("hbring"), "the output directory of H-bond detailed information");
 
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
  if (vm.count("mol-defect")){
    p_mol = true;
  }

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
  RingStats rs (nmolecules, max_ring);

  int countread = 0;
  ofstream fout (ofile.c_str());
  rs.print_head (fout);
  
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
    // open_mol_defect (odir, nmolecules);
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
    if (((countread++)) % 1 == 0){
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

    vector<vector<int > > h_list;
    hba.computeBondList (h_list, clist, vect_box, waters);
    // for (unsigned ii = 0; ii < h_list.size(); ++ii){
    //   cout << ii << " \t " << h_list[ii].size() << " \t ";
    //   for (unsigned jj = 0; jj < h_list[ii].size(); ++jj){
    // 	cout << h_list[ii][jj] << " " ;
    //   }
    //   cout << endl;
    // }
    
    RingAnalysis ra;
    vector<vector<vector<int > > > r_list;
    ra.compute (r_list, h_list, max_ring, func_numb_threads);
    vector<vector<int > > ur_list;
    ra.unique_list (ur_list, r_list);

    rs.sys_deposite (ur_list);
    rs.print_frame (fout, time);

    // for (unsigned jj = 0; jj < r_list[0].size(); ++jj){
    //   cout << "ring " << jj << "  " ;
    //   for (unsigned kk = 0; kk < r_list[0][jj].size(); ++kk){
    // 	cout << r_list[0][jj][kk] << " ";
    //   }
    //   cout << endl;
    // }

    // for (unsigned tt = 0; tt < r_list.size(); ++tt){
    //   cout << " atom " << tt << endl;
    //   for (unsigned jj = 0; jj < r_list[tt].size(); ++jj){
    // 	cout << "ring " << jj << "  " ;
    // 	for (unsigned kk = 0; kk < r_list[tt][jj].size(); ++kk){
    // 	  cout << r_list[tt][jj][kk] << " ";
    // 	}
    // 	cout << endl;
    //   }
    //   // break;
    // }

    // for (unsigned tt = 0; tt < ur_list.size(); ++tt){
    //   cout << tt << " \t " ;
    //   for (unsigned kk = 0; kk < ur_list[tt].size(); ++kk){
    // 	cout << ur_list[tt][kk] << " ";
    //   }
    //   cout << endl;
    // }
    
    if (p_detail){
      print_step (odir, step, time, ur_list);
    }
    if (p_mol){
      // print_mol_defect (odir, step, time, hba.step_count_don, hba.step_count_acc, func_numb_threads);
    }
  }
  printf ("\n");
  
  free (xx);  
  xdrfile_close (fp);

  return 0;
}
