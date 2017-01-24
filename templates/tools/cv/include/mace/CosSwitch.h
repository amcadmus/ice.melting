#ifndef __MACE_CosSwitch_h_wanghan__
#define __MACE_CosSwitch_h_wanghan__

#include "mace/Function.h"

namespace MOASP{
  namespace MACE{

    template<typename TYPE>
    class CosSwitch : public Function <TYPE> 
    {
  public:
      CosSwitch ();
      CosSwitch (const TYPE & rmin,
		 const TYPE & rmax);
      virtual ~CosSwitch () {};
      void reinit (const TYPE & rmin,
		   const TYPE & rmax);
  public:
      virtual void calculate () const;
  public:
      virtual unsigned inputDim () const {return 1;}
      virtual unsigned valueDim () const {return 1;}
      virtual unsigned derivDim () const {return 1;}
  private:
      TYPE rmin, rmax;
    };
  }
}

#ifdef MOASP_INLINE_IMPLEMENTATION
#include "CosSwitch_Impl.h"
#endif

#endif


