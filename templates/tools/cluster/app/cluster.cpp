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
#include <stdio.h>

#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
#include "Defines.h"
#include "CellList.h"
#include "Clusters.h"
#include "StringSplit.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

void 
read_sel (FILE * fp_sel, 
	  const int step, 
	  vector<int > & selection)
{
  selection.clear();
  char * line_buff (NULL);
  size_t numb = 0;
  getline (&line_buff, &numb, fp_sel);
  vector<string> words;
  StringOperation::split (string(line_buff), words);
  if (step != atoi (words[0].c_str())){
    cerr << "the selection does not match the step, exit" << endl;
    exit (1);
  }
  for (unsigned ii = 1; ii < words.size(); ++ii){
    selection.push_back (atoi(words[ii].c_str()));
  }
  free (line_buff);
}


void 
print_numb_cluster (FILE * fout, 
		    const int step,
		    const float time,
		    const vector<vector<int > > & clstr_list)
{
  int count = 0;
  unsigned max = 0;
  for (unsigned ii = 0; ii < clstr_list.size(); ++ii){
    if (clstr_list[ii].size() != 0) count ++;
    if (clstr_list[ii].size() > max) max = clstr_list[ii].size();
  }
  fprintf (fout, "%09d %d %d\n", step, count, max);
}

void
print_cluster (const string odir, 
	       const int step,
	       const float time,	       
	       const vector<vector<int > > & clstr_list, 
	       const vector<int> & selection)
{
  FILE * fp;
  char name[1024];
  sprintf (name, "%s/step_%09d", odir.c_str(), step);
  fp = fopen (name, "w");
  if (fp == NULL){
    cerr << "cannot open file " << string(name) << endl;
    exit (1);
  }
  
  for (unsigned ii = 0; ii < clstr_list.size(); ++ii){
    if (clstr_list[ii].size() != 0){
      fprintf (fp, "%d ", ii);
      for (unsigned jj = 0; jj < clstr_list[ii].size(); ++jj){
	fprintf (fp, "%d ", selection[clstr_list[ii][jj]]);
      }
      fprintf (fp, "\n");
    }
  }

  fclose (fp);  
}



int main(int argc, char * argv[])
{
  double begin, end;
  string ifile, isfile, ofile, odir;
  int func_numb_threads;
  int numb_mol_atom;
  double rcut;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b", po::value<double > (&begin)->default_value(0.f), "start time")
      ("end,e",   po::value<double > (&end  )->default_value(0.f), "end   time")
      ("rcut,c",   po::value<double > (&rcut)->default_value(0.3), "the cut-off of the cluster.")
      ("numb-mol-atom", po::value<int > (&numb_mol_atom)->default_value(4), "number of sites in the water molecule")
      ("numb-threads,t", po::value<int > (&func_numb_threads)->default_value(1), "number of threads")
      ("input-traj,f",   po::value<string > (&ifile)->default_value ("traj.xtc"), "the input .xtc file")
      ("input-sel,s",   po::value<string > (&isfile)->default_value ("sel.out"), "the input molecular selection")
      ("output,o",  po::value<string > (&ofile)->default_value ("cluster.out"), "the output file")
      ("output-dir",  po::value<string > (&odir)->default_value ("cluster"), "the output directory of the cluster definition of each frame");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    cout << desc<< "\n";
    return 0;
  }

  cout << "###################################################" << endl;
  cout << "# computes the local Steinhardt parameter" << endl;
  cout << "# begin->end: " << begin << " " << end << endl;
  cout << "# numb sites in water: " << numb_mol_atom << endl;
  cout << "# input traj: " << ifile << endl;
  cout << "# input sel : " << isfile << endl;
  cout << "# output: " << ofile << endl;
  cout << "# rcut: " << rcut << endl;
  cout << "###################################################" << endl;  
  
  if (access (odir.c_str(), 0) == -1) {
    cout << "# dir " << odir << " does not exist, create." << endl;
    if (mkdir (odir.c_str(), 0755)){
      cerr << "# dir " << odir << " creation failed." << endl;
      return 1;
    }
  }

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
  
  vector<vector<double > > coms;
  coms.reserve (natoms / numb_mol_atom);
  vector<int> selection;
  Clusters clstr (rcut);  

  int countread = 0;
  FILE *fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cerr << "cannot open file " << ofile << endl;
    exit (1);
  }

  fp = xdrfile_open (ifile.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << ifile << endl;
    exit (1);
  }
  
  FILE * fp_sel = fopen (isfile.c_str(), "r");
  if (fp_sel == NULL){
    cerr << "cannot open file " << isfile << endl;
    exit (1);
  }

  while (read_xtc (fp, natoms, &step, &time, box, xx, &prec) == 0){
    read_sel (fp_sel, step, selection);

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
    for (int ii = 0; ii < nmol; ++ii){
      vector<double > read_com(3, 0.);
      for (int dd = 0; dd < 3; ++dd){
	read_com[dd] = xx[ii*numb_mol_atom][dd];
	if      (read_com[dd] <  0          ) read_com[dd] += box[dd][dd];
	else if (read_com[dd] >= box[dd][dd]) read_com[dd] -= box[dd][dd];
      }
      coms.push_back (read_com);
    }
    vector<vector<double > > sel_coms;
    for (unsigned ii = 0; ii < selection.size(); ++ii){
      sel_coms.push_back (coms[selection[ii]]);
    }

    vector<double > boxsize(3);
    for (int dd = 0; dd < 3; ++dd) boxsize[dd] = box[dd][dd];

    clstr.analysis (sel_coms, boxsize);

    const vector<vector<int > > & clstr_list = clstr.get_cluster();

    print_numb_cluster (fout, step, time, clstr_list);
    print_cluster (odir, step, time, clstr_list, selection);
    
  }
  printf ("\n");

  
  free (xx);  
  fclose (fout);
  fclose (fp_sel);
  xdrfile_close (fp);

  return 0;
}
