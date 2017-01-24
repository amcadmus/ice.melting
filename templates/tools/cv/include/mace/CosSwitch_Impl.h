#ifndef __MACE_CosSwitch_Impl_h_wanghan__
#define __MACE_CosSwitch_Impl_h_wanghan__

#include <cmath>

template <typename TYPE>
MOASP::MACE::CosSwitch<TYPE>::
CosSwitch ()
    :
    rmin(0),
    rmax(0)
{
}

template <typename TYPE>
MOASP::MACE::CosSwitch<TYPE>::
CosSwitch (const TYPE & rmin_,
	   const TYPE & rmax_)
    :
    rmin(rmin_),
    rmax(rmax_)
{
}

template <typename TYPE>
void
MOASP::MACE::CosSwitch<TYPE>::
reinit (const TYPE & rmin_,
	const TYPE & rmax_)
{
  rmin = (rmin_);
  rmax = (rmax_);
}

// template <typename TYPE>
// void
// MOASP::MACE::CosSwitch<TYPE>::
// registerData (const TYPE *const p_in,
// 	      TYPE *const p_v,
// 	      TYPE *const p_d)
// {
//   Function<TYPE>::registerData (p_in, p_v, p_d);
// }

template <typename TYPE>
void
MOASP::MACE::CosSwitch<TYPE>::
calculate () const 
{  
  const TYPE & xx (*(Function<TYPE>::getPtrInput()[0]));
  TYPE & vv (*(Function<TYPE>::getPtrValue()[0]));
  TYPE & dd (*(Function<TYPE>::getPtrDeriv()[0]));

  if (xx >= 0){
    if (xx < rmin) {
      dd = 0;
      vv = 1;
    }
    else if (xx < rmax){
      TYPE value = (xx - rmin) / (rmax - rmin) * M_PI;
      dd = -0.5 * sin(value) * M_PI / (rmax - rmin);
      vv = 0.5 * (cos(value) + 1);
    }
    else {
      dd = 0;
      vv = 0;
    }
  }
  else {
    if (xx > -rmin){
      dd = 0;
      vv = 1;
    }
    else if (xx > -rmax){
      TYPE value = (-xx - rmin) / (rmax - rmin) * M_PI;
      dd = 0.5 * sin(value) * M_PI / (rmax - rmin);
      vv = 0.5 * (cos(value) + 1);      
    }
    else {
      dd = 0;
      vv = 0;
    }
  }

}



#endif

