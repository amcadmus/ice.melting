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

void 
statistic_sq (vector<double > & stat,
	      vector<unsigned > & count,
	      const double hh,
	      const vector<double > & qq,
	      const vector<double > & sq)
{
  unsigned nbin = stat.size();
  assert (nbin == count.size());
  unsigned nthreads = omp_get_num_threads() ;
  
  vector<vector<unsigned > > thread_count (nthreads);
  vector<vector<double > > thread_stat (nthreads);
#pragma omp parallel for
  for (unsigned tt = 0; tt < nthreads; ++tt){
    thread_stat[tt].resize (nbin,0.);
    thread_count[tt].resize (nbin,0);
  }
  
#pragma omp parallel for
  for (unsigned tt = 0; tt < nthreads; ++tt){
    for (unsigned ii = 0; ii < qq.size(); ++ii){
      int idx = qq[ii] / hh;
      if (idx >= int(nbin) || idx < 0) continue;
      thread_stat[tt][idx] += sq[ii];
      thread_count[tt][idx] ++;
    }
  }

#pragma omp parallel for
  for (unsigned ii = 0; ii < nbin; ++ii){
    for (unsigned tt = 0; tt < nthreads; ++tt){
      stat[ii] += thread_stat[tt][ii];
      count[ii] += thread_count[tt][ii];
    }
  }  
}


int main(int argc, char * argv[])
{
  double begin, end;
  string ifile, ofile, odir;
  int func_numb_threads;
  int numb_mol_atom;
  unsigned KK[3];
  bool p_detail (false);
  double qbin (0.05);
  double qmax (45);
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b", po::value<double > (&begin)->default_value(0.f), "start time")
      ("end,e",   po::value<double > (&end  )->default_value(0.f), "end   time")
      ("detail", "print the displacement of each molecule at each step")
      ("q-bin",   po::value<double > (&qbin )->default_value(0.05), "the bin size of sq computing ")
      ("q-max",   po::value<double > (&qmax )->default_value(45), "the range of print sq, default for water O-H bond ")
      ("kx", po::value<unsigned > (&KK[0])->default_value(12), "number of modes in dir x")
      ("ky", po::value<unsigned > (&KK[1])->default_value(12), "number of modes in dir y")
      ("kz", po::value<unsigned > (&KK[2])->default_value(12), "number of modes in dir z")
      ("numb-mol-atom", po::value<int > (&numb_mol_atom)->default_value(4), "number of sites in the water molecule")
      ("numb-threads,t", po::value<int > (&func_numb_threads)->default_value(1), "number of threads")
      ("input,f",   po::value<string > (&ifile)->default_value ("traj.xtc"), "the input .xtc file")
      ("output,o",  po::value<string > (&ofile)->default_value ("sq.out"), "the output file")
      ("output-dir",  po::value<string > (&odir)->default_value ("sq"), "the output directory of displacement value information");
  
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

  cout << "###################################################" << endl;
  cout << "# computes the displacement of molecules" << endl;
  cout << "# begin->end: " << begin << " " << end << endl;
  cout << "# numb sites in water: " << numb_mol_atom << endl;
  cout << "# qmax: " << qmax << endl;
  cout << "# qbin: " << qbin << endl;
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
  StructureFactor sf (KK, nmolecules, func_numb_threads);
  vector<double > nq;
  vector<complex<double > > sq;
  vector<double > norm_sq;
  vector<pair<double, double > > sq_print;
  vector<double > qstat (qmax / qbin + 1, 0.);
  vector<unsigned > countstat (qmax / qbin + 1, 0);

  int countread = 0;
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
    vector<double >  mybox(3);
    for (unsigned dd = 0; dd < 3; ++dd) mybox[dd] = box[dd][dd];
    
    sf.getNormQ (nq, mybox);
    sf.compute (sq, coms, mybox);
    norm_sq.resize (sq.size());
#pragma omp parallel for
    for (unsigned ii = 0; ii < sq.size(); ++ii){
      norm_sq[ii] = norm(sq[ii]) / double(nmol);
    }
    statistic_sq (qstat, countstat, qbin, nq, norm_sq);
    
    // fprintf (fout, "%f\t %f\n", time, p_sa->getStepQ());
    // if (p_detail){
    //   print_step (odir, step, time, displ.get_displ());
    // }
  }
  printf ("\n");

  for (unsigned ii = 0; ii < qstat.size(); ++ii){
    if (countstat[ii] == 0) continue;
    qstat[ii] = qstat[ii] / double(countstat[ii]);
  }

  FILE *fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cerr << "cannot open file " << ofile << endl;
    exit (1);
  }
  for (unsigned ii = 0; ii < qstat.size(); ++ii){
    if (ii == 0) continue;
    fprintf (fout, "%f %f\n", (ii+0.5) * qbin, qstat[ii]);
  }
  fclose (fout);

  free (xx);  
  xdrfile_close (fp);

  return 0;
}
