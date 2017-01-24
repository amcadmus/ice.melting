#ifndef __PointerArray_h_wanghan__
#define __PointerArray_h_wanghan__

#include <stdlib.h>
#include <assert.h>

namespace MOASP {
  template <typename TYPE>
  class PointerArray 
  {
public:
    typedef TYPE * PointerType;
    PointerArray ();
    PointerArray (const unsigned, const PointerType p = NULL);
    PointerArray (const PointerArray<TYPE> & a) ;
    ~PointerArray ();
public:
    void resize (const unsigned, const PointerType p = NULL);
    unsigned size () const {return numb;}
    bool empty () const {return 0 == size();}
    void push_back (const PointerType &);
    void clear ();
    PointerType& operator [] (const unsigned &i) {return buff[i];}
    const PointerType& operator [] (const unsigned &i) const {return buff[i];}
    PointerType& back () ;
    const PointerType& back () const;
    // void set (const unsigned i, PointerType p);
    PointerArray<TYPE> & operator = (const PointerArray<TYPE> & array) {return copy(array);}
    PointerArray<TYPE> & copy (const PointerArray<TYPE> & a);
private:
    PointerType *buff;
    unsigned numb;
  };
}

template<typename TYPE>
MOASP::PointerArray<TYPE>::
PointerArray () :
    buff (NULL), numb(0)
{
}

template<typename TYPE>
MOASP::PointerArray<TYPE>::
PointerArray (const unsigned n, const PointerType p) :
    buff (NULL), numb(0)
{
  resize (n, p);
}

template<typename TYPE>
MOASP::PointerArray<TYPE>::
PointerArray (const PointerArray<TYPE> & a) :
    buff (NULL), numb(0)
{
  copy (a);
}

template<typename TYPE>
MOASP::PointerArray<TYPE>::
~PointerArray () 
{
  clear ();
}

template<typename TYPE>
void
MOASP::PointerArray<TYPE>::
clear () 
{
  free (buff);
  buff = NULL;
  numb = 0;
  // resize(0);
}

template<typename TYPE>
MOASP::PointerArray<TYPE> &
MOASP::PointerArray<TYPE>::
copy (const PointerArray & array)
{
  // identity test!
  if (buff != array.buff){
    // clear();
    resize (array.size());
    for (unsigned ii = 0; ii < size() ; ++ii){
      // set(ii,  array[ii] );
      // this->operator[](ii) = array[ii];
      buff[ii] = array.buff[ii];
    }
  }
  else {
    assert (numb == array.numb);
  }
  return *this;
}

template<typename TYPE>
inline void
MOASP::PointerArray<TYPE>::
resize (const unsigned new_numb, const PointerType ptr)
{
  if (new_numb == numb) return;
  buff = (PointerType*) realloc (buff, sizeof(PointerType) * new_numb);
  for (unsigned ii = numb; ii < new_numb; ++ii){
    buff[ii] = ptr;
  }
  numb = new_numb;
}

template<typename TYPE>
inline void
MOASP::PointerArray<TYPE>::
push_back (const PointerType & ptr)
{
  resize (numb+1);
  buff[numb-1] = ptr;
}

template<typename TYPE>
inline typename MOASP::PointerArray<TYPE>::PointerType&
MOASP::PointerArray<TYPE>::
back () 
{
  assert (size() != 0);
  return buff[size()-1];
}

template<typename TYPE>
inline const typename MOASP::PointerArray<TYPE>::PointerType&
MOASP::PointerArray<TYPE>::
back () const
{
  assert (size() != 0);
  return buff[size()-1];
}


#endif

