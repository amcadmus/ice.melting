#pragma once

#include "mace/Function.h"
#include "mace/CosSwitch.h"

namespace MOASP{
  namespace MACE{

    template<typename TYPE>
    class CosHill : public Function <TYPE> 
    {
  public:
      CosHill ();
      CosHill (const TYPE & r0,
	       const TYPE & r1,
	       const TYPE & r2,
	       const TYPE & r3);
      virtual ~CosHill () {};
      void reinit (const TYPE & r0,
		   const TYPE & r1,
		   const TYPE & r2,
		   const TYPE & r3);
  public:
      virtual void calculate () const;
  public:
      virtual unsigned inputDim () const {return 1;}
      virtual unsigned valueDim () const {return 1;}
      virtual unsigned derivDim () const {return 1;}
  private:
      mutable CosSwitch<TYPE> neg_swith;
      mutable CosSwitch<TYPE> pos_swith;
      mutable TYPE neg_i, neg_v, neg_d;
      mutable TYPE pos_i, pos_v, pos_d;
    };
  }
}

#ifdef MOASP_INLINE_IMPLEMENTATION
#include "CosHill_Impl.h"
#endif
