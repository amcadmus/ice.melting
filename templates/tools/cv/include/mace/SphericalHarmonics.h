#ifndef __MACE_SphericalHarmonics_h_wanghan__
#define __MACE_SphericalHarmonics_h_wanghan__

#include "mace/Function.h"

namespace MOASP{
  namespace MACE{

    template<typename TYPE, unsigned LL>
    class SphericalHarmonics : public Function <TYPE> 
    {
  public:
      SphericalHarmonics ();
      virtual ~SphericalHarmonics () {};
  public:
      virtual void calculate () const;
  public:
      virtual unsigned inputDim () const {return 3;}
      virtual unsigned valueDim () const {return 2*(LL+1);}
      virtual unsigned derivDim () const {return 6*(LL+1);}
  private:
      TYPE normalize	[LL+1];
      TYPE coeff_poly	[LL+1];
      TYPE coeff_matrix [LL+1][LL+1];
      void makeCoeffMatrix ();
      TYPE deriv_poly (const unsigned & m, const TYPE & val, TYPE & df ) const;
      void deriv_poly(const TYPE val, TYPE* vf, TYPE* df ) const ;
    };
  }
}

// #ifdef MOASP_INLINE_IMPLEMENTATION
#include "SphericalHarmonics_Impl.h"
// #endif

#endif


