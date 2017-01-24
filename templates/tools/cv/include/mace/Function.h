#ifndef __MACE_Function_h_wanghan__
#define __MACE_Function_h_wanghan__

// #include "mpi.h"
#include "PointerArray.h"
#include "mace/CollectiveVariable.h"

namespace MOASP{
  namespace MACE {

    template<typename TYPE>
    class Function : public CollectiveVariable
    {
  public:
      Function () ;
      virtual ~Function () {};
  public:
      void registerData (const PointerArray<TYPE> & in,
			 const PointerArray<TYPE> & v,
			 const PointerArray<TYPE> & d) 
	  { p_in = in; p_v = v; p_d = d;}
      virtual void check () const;
      // void registerMPIComm (const MPI_Comm comm)	{my_comm = comm;}
  public:
      const PointerArray<TYPE> & getPtrInput () const {return p_in;}
      const PointerArray<TYPE> & getPtrValue () const {return p_v;}
      const PointerArray<TYPE> & getPtrDeriv () const {return p_d;}
      // const MPI_Comm &		 getMPIComm  () const {return my_comm;}
  private:
      PointerArray<TYPE> p_in, p_v, p_d;
      // MPI_Comm my_comm;
    };
  }
}

#ifdef MOASP_INLINE_IMPLEMENTATION
#include "Function_Impl.h"
#endif

#endif


