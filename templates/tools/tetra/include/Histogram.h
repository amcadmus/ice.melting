#ifndef __Tool_Histogram_h_wanghan__
#define __Tool_Histogram_h_wanghan__

#include <vector>

using namespace std;

class Histogram
{
public:
  Histogram		(const unsigned nbin = 10,
			 const double lower = 0.,
			 const double upper = 0.);
  void reinit		(const unsigned nbin = 10,
			 const double lower = 0.,
			 const double upper = 0.);
public:
  void processData	(const vector<double > & vec) ;
  void getResult	(vector<double > & xx,
			 vector<double > & vv);
private:
  void setBoundary	(const double lower_,
			 const double upper_);
  void setNumbBins	(const unsigned nbin)			{n_bin = nbin;}  
  // void setNumbBlocks	(const unsigned nblock)			{n_block = nblock;}
  double lower;
  double upper;
  unsigned n_bin;
  double hh;
  vector<double > value;
}
    ;


class Histogram2D 
{
public:
  Histogram2D		(const unsigned nbinx = 10,
			 const unsigned nbiny = 10,
			 const double lowerx = 0.,
			 const double lowery = 0.,
			 const double upperx = 0.,
			 const double uppery = 0.);
  void reinit		(const unsigned nbinx = 10,
			 const unsigned nbiny = 10,
			 const double lowerx = 0.,
			 const double lowery = 0.,
			 const double upperx = 0.,
			 const double uppery = 0.);
public:
  void processData	(const vector<double > & vecxx, 
			 const vector<double > & vecyy) ;
  void getResult	(vector<double > & xx,
			 vector<double > & yy,
			 vector<vector<double > > & vv);
private:
  void setBoundary	(const double lowerx,
			 const double lowery,
			 const double upperx,
			 const double uppery);
  void setNumbBins	(const unsigned nbinx_, 
			 const unsigned nbiny_) {nbinx = nbinx_; nbiny = nbiny_;}  
  double lowerx;
  double lowery;
  double upperx;
  double uppery;
  unsigned nbinx, nbiny;
  double hx, hy;
  vector<vector<double> > value;  
};

#endif


