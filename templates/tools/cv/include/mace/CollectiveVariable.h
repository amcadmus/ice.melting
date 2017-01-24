#ifndef __MACE_CollectiveVariable_h_wanghan__
#define __MACE_CollectiveVariable_h_wanghan__

namespace MOASP{
  namespace MACE {    
    
    class CollectiveVariable 
    {
  public:
      virtual ~CollectiveVariable () {};
  public:
      virtual void check () const = 0;
      virtual unsigned inputDim () const = 0;
      virtual unsigned valueDim () const = 0;
      virtual unsigned derivDim () const = 0;
      // operations.. can be viewed as operations on memory... calcualte both value and deriv
      virtual void calculate () const = 0;
      virtual unsigned computeDerivIndex (const unsigned vidx, const unsigned iidx) const;
      // virtual void value (const ValueIterator<Acceleration> &ibegin,
      // 			  const ValueIterator<Acceleration> &iend,
      // 			  const ValueIterator<Acceleration> &obegin,
      // 			  const ValueIterator<Acceleration> &oend) = 0;
      // virtual void derivative (const ValueIterator<Acceleration> &ibegin,
      // 			       const ValueIterator<Acceleration> &iend,
      // 			       const DerivIterator<Acceleration> &dbegin,
      // 			       const DerivIterator<Acceleration> &dend) = 0;
    };
    
  }
}

inline unsigned
MOASP::MACE::CollectiveVariable::
computeDerivIndex (const unsigned vidx,
		   const unsigned iidx) const
{
  return vidx * inputDim() + iidx;
}

#endif


