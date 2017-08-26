#ifndef __Defines_h_wanghan__
#define __Defines_h_wanghan__

#include <vector>
using namespace std;

typedef int	Identity;
typedef double	ValueType;

struct IntVectorType
{
  int x, y, z;
  IntVectorType (const int xx=0, const int yy=0, const int zz=0);
};

struct VectorType
{
  ValueType x, y, z;
  VectorType (const ValueType xx=0, const ValueType yy=0, const ValueType zz=0);
};

inline IntVectorType::
IntVectorType (const int xx, const int yy, const int zz)
    : x(xx), y(yy), z(zz)
{
}

inline VectorType::
VectorType (const ValueType xx, const ValueType yy, const ValueType zz)
    : x(xx), y(yy), z(zz)
{
}

inline void
diff ( vector<double> & rdiff,
       const vector<double> & a1,
       const vector<double> & a2,
       const vector<double> & box)
{
  for (int dd = 0; dd < 3; ++dd) rdiff[dd] = a1[dd] - a2[dd];
  vector<int > shift(3, 0);
  for (int dd = 0; dd < 3; ++dd){
    if      (rdiff[dd] < -.5 * box[dd]) shift[dd] += 1;
    else if (rdiff[dd] >= .5 * box[dd]) shift[dd] -= 1;
  }
  for (int dd = 0; dd < 3; ++dd){
    rdiff[dd] += box[dd] * shift[dd];
  }  
}

inline double 
dot (const vector<double> & a1,
     const vector<double> & a2)
{
  return a1[0] * a2[0] + a1[1] * a2[1] + a1[2] * a2[2];
}

inline double 
norm2 (const vector<double> & diff)
{
  return dot (diff, diff);
}

inline double 
dist2 (const vector<double> & a1,
       const vector<double> & a2,
       const vector<double> & box)
{
  vector<double > diff (3);
  for (int dd = 0; dd < 3; ++dd) diff[dd] = a1[dd] - a2[dd];
  vector<int > shift(3, 0);
  for (int dd = 0; dd < 3; ++dd){
    if      (diff[dd] < -.5 * box[dd]) shift[dd] += 1;
    else if (diff[dd] >= .5 * box[dd]) shift[dd] -= 1;
  }
  for (int dd = 0; dd < 3; ++dd){
    diff[dd] += box[dd] * shift[dd];
  }
  return diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
}

#endif
