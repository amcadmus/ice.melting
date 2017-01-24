#ifndef __MACE_Function_Impl_h_wanghan__
#define __MACE_Function_Impl_h_wanghan__

template<typename TYPE>
MOASP::MACE::Function<TYPE>::
Function ()
    // : my_comm (MPI_COMM_NULL)
{
}

template<typename TYPE>
void
MOASP::MACE::Function<TYPE>::
check () const
{
#ifdef DEBUG_CHECK_ASSERTIONS
  // assert (dim_d == dim_in * dim_v);
  // assert (p_in != NULL);
  // assert (p_v != NULL);
  // assert (p_d != NULL);
#endif
}


#endif


