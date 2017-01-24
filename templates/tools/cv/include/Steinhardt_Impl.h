#ifndef __MACE_Steinhardt_Impl_h_wanghan__
#define __MACE_Steinhardt_Impl_h_wanghan__

template<typename Acceleration, unsigned LL>
MOASP::MACE::Steinhardt<Acceleration,LL>::
Steinhardt (const typename Acceleration::ValueType & rmin,
	    const typename Acceleration::ValueType & rmax)
{
  reinit (rmin, rmax);
}

template<typename Acceleration, unsigned LL>
void
MOASP::MACE::Steinhardt<Acceleration,LL>::
reinit (const typename Acceleration::ValueType & rmin,
	const typename Acceleration::ValueType & rmax)
{
  cs.reinit (rmin, rmax);
  rc = rmax;
  rc2 = rc * rc;
  PointerArray<typename Acceleration::ValueType> ai, av, ad;
  ai.push_back (&r1);
  av.push_back (&vcs);
  ad.push_back (&dcs);
  cs.registerData (ai, av, ad);
  prefactor = 1. * sqrt(4*M_PI / (2. * double(LL) + 1.));
}

template<typename Acceleration, unsigned LL>
void
MOASP::MACE::Steinhardt<Acceleration,LL>::
registerFunctionData (const PointerArray<typename Acceleration::ValueType>& dist_,
		      const PointerArray<typename Acceleration::ValueType>& value_,
		      const PointerArray<typename Acceleration::ValueType>& dvalue_) const
{
  dist = dist_;
  value = value_;
  dvalue = dvalue_;
  sh.registerData (dist, value, dvalue);
}

template<typename Acceleration, unsigned LL>
void
MOASP::MACE::Steinhardt<Acceleration,LL>::
evaluate () const 
{
  const typename Acceleration::ValueType & dx (*dist[0]);
  const typename Acceleration::ValueType & dy (*dist[1]);
  const typename Acceleration::ValueType & dz (*dist[2]);

  typename Acceleration::ValueType r2 = dx * dx + dy * dy + dz * dz;
  if (r2 > rc2) return;
  typename Acceleration::ValueType rinv  = MOASP::MathUtilities::invsqrt (r2);
  r1 = r2 * rinv;

  cs.calculate();
  sh.calculate();

  typename Acceleration::ValueType tfx = dcs * dx * rinv;
  typename Acceleration::ValueType tfy = dcs * dy * rinv;
  typename Acceleration::ValueType tfz = dcs * dz * rinv;

  for (unsigned ii = 0; ii < interValueDim(); ++ii){
    typename Acceleration::ValueType vsh = *(value[ii]);
    *(value[ii]) = *(value[ii]) * vcs * prefactor;
    *(dvalue[ii*3+0]) = ( *(dvalue[ii*3+0]) * vcs + vsh * tfx ) * prefactor;
    *(dvalue[ii*3+1]) = ( *(dvalue[ii*3+1]) * vcs + vsh * tfy ) * prefactor;
    *(dvalue[ii*3+2]) = ( *(dvalue[ii*3+2]) * vcs + vsh * tfz ) * prefactor;
  }
}



#endif
