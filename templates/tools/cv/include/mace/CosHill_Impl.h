#pragma once

#include <cmath>

template <typename TYPE>
MOASP::MACE::CosHill<TYPE>::
CosHill ()
{
}

template <typename TYPE>
MOASP::MACE::CosHill<TYPE>::
CosHill (const TYPE & r0,
	 const TYPE & r1,
	 const TYPE & r2,
	 const TYPE & r3)
{
  reinit (r0, r1, r2, r3);
}

template <typename TYPE>
void
MOASP::MACE::CosHill<TYPE>::
reinit (const TYPE & r0,
	const TYPE & r1,
	const TYPE & r2,
	const TYPE & r3)
{
  neg_swith.reinit (r0, r1);
  pos_swith.reinit (r2, r3);
  PointerArray<TYPE> p_neg_i, p_neg_v, p_neg_d;
  p_neg_i.push_back (&neg_i);
  p_neg_v.push_back (&neg_v);
  p_neg_d.push_back (&neg_d);
  neg_swith.registerData (p_neg_i, p_neg_v, p_neg_d);
  PointerArray<TYPE> p_pos_i, p_pos_v, p_pos_d;
  p_pos_i.push_back (&pos_i);
  p_pos_v.push_back (&pos_v);
  p_pos_d.push_back (&pos_d);
  pos_swith.registerData (p_pos_i, p_pos_v, p_pos_d);
}

template <typename TYPE>
void
MOASP::MACE::CosHill<TYPE>::
calculate () const 
{  
  const TYPE & xx (*(Function<TYPE>::getPtrInput()[0]));
  TYPE & vv (*(Function<TYPE>::getPtrValue()[0]));
  TYPE & dd (*(Function<TYPE>::getPtrDeriv()[0]));

  neg_i = xx;
  pos_i = xx;
  neg_swith.calculate();
  pos_swith.calculate();

  vv = pos_v - neg_v;
  dd = pos_d - neg_d;
}

