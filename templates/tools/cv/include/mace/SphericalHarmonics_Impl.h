#ifndef __MACE_SphericalHarmonics_Impl_h_wanghan__
#define __MACE_SphericalHarmonics_Impl_h_wanghan__

template<typename TYPE, unsigned LL>
class SphericalHarmonics_ParameterInit
{
public:
  void set (TYPE * norm, TYPE * coeff) 
      {
	cerr << "SphericalHarmonics with l = " << LL << " has not been implemented" << endl;	
      }
};

template<typename TYPE>
class SphericalHarmonics_ParameterInit <TYPE,3>
{
public:
  void set (TYPE * normalize, TYPE * coeff_poly)
      {
	normalize[0] = sqrt( ( 7.0*6.0 ) / (4.0*M_PI*6.0) );
	normalize[1] =-sqrt( ( 7.0*2.0 ) / (4.0*M_PI*24.0) );
	normalize[2] = sqrt( ( 7.0*1.0 ) / (4.0*M_PI*120.0) );
	normalize[3] =-sqrt( ( 7.0*1.0 ) / (4.0*M_PI*720.0) );
	coeff_poly[0] = 0.0;
	coeff_poly[1] =-1.5;
	coeff_poly[2] = 0.0;
	coeff_poly[3] = 2.5;
      }
};

template<typename TYPE>
class SphericalHarmonics_ParameterInit <TYPE,4>
{
public:
  void set (TYPE * normalize, TYPE * coeff_poly)
      {
	normalize[0] = sqrt( ( 9.0*24.0 ) / (4.0*M_PI*24.0) );
	normalize[1] =-sqrt( ( 9.0* 6.0 ) / (4.0*M_PI*120.0) );
	normalize[2] = sqrt( ( 9.0* 2.0 ) / (4.0*M_PI*720.0) );
	normalize[3] =-sqrt( ( 9.0* 1.0 ) / (4.0*M_PI*5040.0) );
	normalize[4] = sqrt( ( 9.0* 1.0 ) / (4.0*M_PI*40320.0) );
	coeff_poly[0] = 0.375;
	coeff_poly[1] = 0.0;
	coeff_poly[2] =-3.75;
	coeff_poly[3] = 0.0;
	coeff_poly[4] = 4.375;
      }
};

template<typename TYPE>
class SphericalHarmonics_ParameterInit <TYPE,6>
{
public:
  void set (TYPE * normalize, TYPE * coeff_poly)
      {
	normalize[0] = sqrt( ( 13.0*720.0 ) / (4.0*M_PI*720.0) );
	normalize[1] =-sqrt( ( 13.0*120.0 ) / (4.0*M_PI*5040.0) );
	normalize[2] = sqrt( ( 13.0* 24.0 ) / (4.0*M_PI*40320.0) );
	normalize[3] =-sqrt( ( 13.0*  6.0 ) / (4.0*M_PI*362880.0) );
	normalize[4] = sqrt( ( 13.0*  2.0 ) / (4.0*M_PI*3628800.0) );
	normalize[5] =-sqrt( ( 13.0*  1.0 ) / (4.0*M_PI*39916800.0) );
	normalize[6] = sqrt( ( 13.0*  1.0 ) / (4.0*M_PI*479001600.0) );
	coeff_poly[0] =- 0.3125;
	coeff_poly[1] =  0.0;
	coeff_poly[2] =  6.5625;
	coeff_poly[3] =  0.0;
	coeff_poly[4] =-19.6875;
	coeff_poly[5] =  0.0;
	coeff_poly[6] = 14.4375;
      }
};



template<typename TYPE, unsigned LL>
MOASP::MACE::SphericalHarmonics<TYPE, LL>::
SphericalHarmonics ()
{
  SphericalHarmonics_ParameterInit<TYPE, LL> init;
  init.set (normalize, coeff_poly);
  makeCoeffMatrix ();
}


template<typename TYPE, unsigned LL>
void
MOASP::MACE::SphericalHarmonics<TYPE, LL>::
makeCoeffMatrix ()
{
  for (unsigned ii = 0; ii < LL+1; ++ii){
    coeff_matrix[0][ii] = coeff_poly[ii];
  }
  for (unsigned ii = 1; ii < LL+1; ++ii){
    for (unsigned jj = 0; jj < LL+1-ii; ++jj){
      coeff_matrix[ii][jj] = coeff_matrix[ii-1][jj+1] * (jj+1);
    }
    for (unsigned jj = LL+1-ii; jj < LL+1; ++jj){
      coeff_matrix[ii][jj] = 0.;
    }
  }
}

// (aa + bb * ii) * (cc + dd * ii)
template<typename TYPE>
inline void
mul_complex (const TYPE & aa, const TYPE & bb,
	     const TYPE & cc, const TYPE & dd,
	     TYPE & o_r, TYPE & o_i) 
{
  o_r = aa * cc - bb * dd;
  o_i = aa * dd + bb * cc;
}

// generate an array of size MM containing the powers up to MM-1
template<typename TYPE, unsigned MM>
inline void
pow_complex (const TYPE & i_real, const TYPE & i_imag,
	     TYPE * o_real, TYPE * o_imag) 
{
  o_real[0] = TYPE(1);
  o_imag[0] = TYPE(0);
  for (unsigned ii = 1; ii < MM; ++ii){
    mul_complex<TYPE> (o_real[ii-1], o_imag[ii-1], i_real, i_imag, o_real[ii], o_imag[ii]);
  }
}

template<typename TYPE, unsigned LL>
inline TYPE
MOASP::MACE::SphericalHarmonics<TYPE, LL>::
deriv_poly( const unsigned& m, const TYPE& val, TYPE& df ) const {
  TYPE fact=1.0;
  for(unsigned j=1;j<=m;++j) fact=fact*j;
  TYPE res=coeff_poly[m]*fact;

  TYPE pow=1.0, xi=val, dxi=1.0; df=0.0;
  for(unsigned i=m+1;i<=LL;++i){
    TYPE fact=1.0;
    for(unsigned j=i-m+1;j<=i;++j) fact=fact*j;
    res=res+coeff_poly[i]*fact*xi;
    df = df + pow*coeff_poly[i]*fact*dxi;
    xi=xi*val;
    dxi=dxi*val;
    pow+=1.0;
  }
  df = df*normalize[m];
  return normalize[m]*res;
}

template<typename TYPE, unsigned LL>
inline void
MOASP::MACE::SphericalHarmonics<TYPE, LL>::
deriv_poly(const TYPE val,
	   TYPE *vf,
	   TYPE *df ) const
{
  TYPE rsh[LL+1];
  // computeRsh<LL> (val, rsh);
  rsh[0] = 1.;
  for (unsigned ii = 1; ii < LL+1; ++ii){
    rsh[ii] = rsh[ii-1] * val;
  }

  for (unsigned ii = 0; ii < LL+1; ++ii){
    vf [ii] = TYPE(0.);
  }
  for (unsigned ii = 0; ii < LL+1; ++ii){
    for (unsigned jj = 0; jj < LL+1-ii; ++jj){
      vf[ii] += coeff_matrix[ii][jj] * rsh[jj];
    }
  }
  df[LL] = 0;
  for (unsigned ii = 0; ii < LL; ++ii){
    df[ii] = vf[ii+1];
  }
  for (unsigned ii = 0; ii < LL+1; ++ii){
    vf[ii] *= normalize[ii];
    df[ii] *= normalize[ii];
  }
}


template<typename TYPE, unsigned LL>
void
MOASP::MACE::SphericalHarmonics<TYPE, LL>::
calculate () const
{  
  const TYPE & dx (*(this->getPtrInput()[0]));
  const TYPE & dy (*(this->getPtrInput()[1]));
  const TYPE & dz (*(this->getPtrInput()[2]));
  PointerArray<TYPE> pv = this->getPtrValue();
  PointerArray<TYPE> pd = this->getPtrDeriv();

  const TYPE r2 = dx * dx + dy * dy + dz * dz;
  const TYPE rinv  = 1./sqrt(r2);
  const TYPE dlen = r2 * rinv;
  const TYPE dleni = TYPE(1.)/dlen;
  const TYPE dlen3 = dlen * dlen * dlen;
  const TYPE dlen3i = TYPE(1.)/dlen3;
  TYPE vdz[3];
  vdz[0] = (-( dz * dlen3i )* dx);
  vdz[1] = (-( dz * dlen3i )* dy);
  vdz[2] = (-( dz * dlen3i )* dz) + (dleni);

  // Calculate Legendre Polynomial
  TYPE mat_poly_ass[LL+1], mat_dpoly_ass[LL+1];
  deriv_poly (dz*dleni, mat_poly_ass, mat_dpoly_ass );
  
  // deal with m == 0
  TYPE myrealvec[3];
  myrealvec[0] = mat_dpoly_ass[0]*vdz[0];
  myrealvec[1] = mat_dpoly_ass[0]*vdz[1];
  myrealvec[2] = mat_dpoly_ass[0]*vdz[2];
  // And store the vector function
  *(pv[0]) = mat_poly_ass[0];
  *(pv[1]) = 0;
  // Accumulate the derivatives
  *(pd[0]) = myrealvec[0];
  *(pd[1]) = myrealvec[1];
  *(pd[2]) = myrealvec[2];
  *(pd[3 + 0]) = 0;
  *(pd[3 + 1]) = 0;
  *(pd[3 + 2]) = 0;

  const TYPE com1_real = dx * dleni;
  const TYPE com1_imag = dy * dleni;
  TYPE powered_real[LL+1];
  TYPE powered_imag[LL+1];  
  pow_complex<TYPE,LL+1> (com1_real, com1_imag, powered_real, powered_imag);
  
  TYPE tmp_real[3], tmp_imag[3];
  tmp_real[0] = ( (1.0*dleni) - (dx*dx)*dlen3i );
  tmp_imag[0] = ( -(dx*dy)*dlen3i );
  tmp_real[1] = tmp_imag[0];
  tmp_imag[1] = ( (1.0*dleni) - (dy*dy)*dlen3i );
  tmp_real[2] = ( -(dx*dz)*dlen3i );
  tmp_imag[2] = ( -(dy*dz)*dlen3i );
  
  for (unsigned mm = 1; mm < LL+1; ++mm){
    const TYPE md (static_cast<const TYPE>(mm));
    const TYPE & real_z ( powered_real[mm] );
    const TYPE & imag_z ( powered_imag[mm] );
    TYPE real_dz[3], imag_dz[3];
    mul_complex (tmp_real[0], tmp_imag[0], powered_real[mm-1], powered_imag[mm-1], real_dz[0], imag_dz[0]);
    mul_complex (tmp_real[1], tmp_imag[1], powered_real[mm-1], powered_imag[mm-1], real_dz[1], imag_dz[1]);
    mul_complex (tmp_real[2], tmp_imag[2], powered_real[mm-1], powered_imag[mm-1], real_dz[2], imag_dz[2]);
    real_dz[0] *= md;
    real_dz[1] *= md;
    real_dz[2] *= md;
    imag_dz[0] *= md;
    imag_dz[1] *= md;
    imag_dz[2] *= md;

    TYPE tq6=mat_poly_ass[mm]*real_z;   // Real part of steinhardt parameter
    TYPE myrealvec[3];
    myrealvec[0] = ( mat_dpoly_ass[mm]*real_z*vdz[0]  + mat_poly_ass[mm]*real_dz[0] ); 
    myrealvec[1] = ( mat_dpoly_ass[mm]*real_z*vdz[1]  + mat_poly_ass[mm]*real_dz[1] ); 
    myrealvec[2] = ( mat_dpoly_ass[mm]*real_z*vdz[2]  + mat_poly_ass[mm]*real_dz[2] ); 
    *(pv[mm*2 + 0]) = tq6;
    *(pd[mm*6 + 0]) = myrealvec[0];
    *(pd[mm*6 + 1]) = myrealvec[1];
    *(pd[mm*6 + 2]) = myrealvec[2];
    
    TYPE itq6=mat_poly_ass[mm]*imag_z;  // Imaginary part of steinhardt parameter
    TYPE myimagvec[3];
    myimagvec[0] = mat_dpoly_ass[mm]*imag_z*vdz[0] + mat_poly_ass[mm]*imag_dz[0];
    myimagvec[1] = mat_dpoly_ass[mm]*imag_z*vdz[1] + mat_poly_ass[mm]*imag_dz[1];
    myimagvec[2] = mat_dpoly_ass[mm]*imag_z*vdz[2] + mat_poly_ass[mm]*imag_dz[2];
    *(pv[mm*2 + 1]) = itq6;
    *(pd[mm*6 + 3]) = myimagvec[0];
    *(pd[mm*6 + 4]) = myimagvec[1];
    *(pd[mm*6 + 5]) = myimagvec[2];
  }
}

#endif

