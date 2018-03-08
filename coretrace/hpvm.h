
/*******************************************************************************************************
*  Copyright 2018 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  ("Alliance") under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as "SolTrace". Except to comply with the 
*  foregoing, the term "SolTrace", or any confusingly similar designation may not be used to refer to 
*  any modified version of this software or any modified version of the underlying software originally 
*  provided by Alliance without the prior written consent of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/


#ifndef __hpvm_h
#define __hpvm_h

#include <assert.h>


#if defined(_DEBUG) && defined(_MSC_VER) && defined(_WIN32) && !defined(_WIN64)
#define VEC_ASSERT(x) {if(!(x)) _asm{int 0x03}}
#else
#define VEC_ASSERT(X) assert(X)
#endif

template< typename T >
class matrix_t
{
protected:
	T *t_array;
	size_t n_rows, n_cols;
public:

	matrix_t()
	{
		t_array = new T[1];
		n_rows = n_cols = 1;
	}

	matrix_t( const matrix_t &rhs )
	{
		t_array = 0;
		n_rows = n_cols = 0;
		copy( rhs );
	}

	matrix_t(size_t nr, size_t nc)
	{
		n_rows = n_cols = 0;
		t_array = NULL;
		if (nr < 1) nr = 1;
		if (nc < 1) nc = 1;
		resize(nr,nc);
	}
		
	matrix_t(size_t nr, size_t nc, const T &val)
	{
		n_rows = n_cols = 0;
		t_array = NULL;
		if (nr < 1) nr = 1;
		if (nc < 1) nc = 1;
		resize(nr,nc);
		fill(val);
	}


	virtual ~matrix_t()
	{
		if (t_array) delete [] t_array;
	}
		
	void clear()
	{
		if (t_array) delete [] t_array;
		n_rows = n_cols = 1;
		t_array = new T[1];
	}
		
	void copy( const matrix_t &rhs )
	{
		if (this != &rhs)
		{
			resize( rhs.nrows(), rhs.ncols() );
			size_t nn = n_rows*n_cols;
			for (size_t i=0;i<nn;i++)
				t_array[i] = rhs.t_array[i];
		}
	}

	void assign( const T *pvalues, size_t len )
	{
		resize( len );
		if ( n_cols == len && n_rows == 1 )
			for (size_t i=0;i<len;i++)
				t_array[i] = pvalues[i];
	}
		
	void assign( const T *pvalues, size_t nr, size_t nc )
	{
		resize( nr, nc );
		if ( n_rows == nr && n_cols == nc )
		{
			size_t len = nr*nc;
			for (size_t i=0;i<len;i++)
				t_array[i] = pvalues[i];
		}
	}

	matrix_t &operator=(const matrix_t &rhs)
	{
		copy( rhs );
		return *this;
	}
		
	matrix_t &operator=(const T &val)
	{
		resize(1,1);
		t_array[0] = val;
		return *this;
	}
		
	inline operator T()
	{
		return t_array[0];
	}
		
	bool equals( const matrix_t & rhs )
	{
		if (n_rows != rhs.n_rows || n_cols != rhs.n_cols)
			return false;
			
		size_t nn = n_rows*n_cols;
		for (size_t i=0;i<nn;i++)
			if (t_array[i] != rhs.t_array[i])
				return false;
			
		return true;
	}
		
	inline bool is_single()
	{
		return (n_rows == 1 && n_cols == 1);
	}
			
	inline bool is_array()
	{
		return (n_rows == 1);
	}
		
	void fill( const T &val )
	{
		size_t ncells = n_rows*n_cols;
		for (size_t i=0;i<ncells;i++)
			t_array[i] = val;
	}

	void resize(size_t nr, size_t nc)
	{
		if (nr < 1 || nc < 1) return;
		if (nr == n_rows && nc == n_cols) return;
			
		if (t_array) delete [] t_array;
		t_array = new T[ nr * nc ];
		n_rows = nr;
		n_cols = nc;
	}

	void resize_fill(size_t nr, size_t nc, const T &val)
	{
		resize( nr, nc );
		fill( val );
	}

	void resize_preserve( size_t nr, size_t nc, const T &val )
	{
		matrix_t<T> old( *this );
		resize( nr, nc );
		fill( val );
		for( size_t r=0;r<nr && r<old.nrows();r++)
			for( size_t c=0;c<nc && c<old.ncols();c++)
				at(r,c) = old(r,c);
	}
		
	void resize(size_t len)
	{
		resize( 1, len );
	}
		
	void resize_fill(size_t len, const T &val)
	{
		resize_fill( 1, len, val );
	}
		
	inline T &at(size_t i)
	{
#ifdef _DEBUG
		VEC_ASSERT( i >= 0 && i < n_cols );
#endif
		return t_array[i];
	}

	inline const T&at(size_t i) const
	{
#ifdef _DEBUG
		VEC_ASSERT( i >= 0 && i < n_cols );
#endif
		return t_array[i];
	}
		
	inline T &at(size_t r, size_t c)
	{
#ifdef _DEBUG
		VEC_ASSERT( r >= 0 && r < n_rows && c >= 0 && c < n_cols );
#endif
		return t_array[n_cols*r+c];
	}

	inline const T &at(size_t r, size_t c) const
	{
#ifdef _DEBUG
		VEC_ASSERT( r >= 0 && r < n_rows && c >= 0 && c < n_cols );
#endif
		return t_array[n_cols*r+c];
	}
		
	inline T &operator()(size_t r, size_t c)
	{
#ifdef _DEBUG
		VEC_ASSERT( r >= 0 && r < n_rows && c >= 0 && c < n_cols );
#endif
		return t_array[n_cols*r+c];
	}

	inline const T &operator()(size_t r, size_t c) const
	{
#ifdef _DEBUG
		VEC_ASSERT( r >= 0 && r < n_rows && c >= 0 && c < n_cols );
#endif
		return t_array[n_cols*r+c];
	}
		
	T operator[] (size_t i) const
	{
#ifdef _DEBUG
		VEC_ASSERT( i >= 0 && i < n_cols );
#endif
		return t_array[i];
	}
		
	T &operator[] (size_t i)
	{
#ifdef _DEBUG
		VEC_ASSERT( i >= 0 && i < n_cols );
#endif
		return t_array[i];
	}
				
	inline size_t nrows() const
	{
		return n_rows;
	}
		
	inline size_t ncols() const
	{
		return n_cols;
	}
		
	inline size_t ncells() const
	{
		return n_rows*n_cols;
	}
		
	inline size_t membytes() const
	{
		return n_rows*n_cols*sizeof(T);
	}
		
	void size(size_t &nr, size_t &nc) const
	{
		nr = n_rows;
		nc = n_cols;
	}
		
	size_t length() const
	{
		return n_cols;
	}
		
	inline T *data()
	{
		return t_array;
	}

	inline T value() const
	{
		return t_array[0];
	}
};


typedef unsigned long muint;
typedef matrix_t<double> HPM2D;

/*
class HPM2D
{
public:
	HPM2D()
	{
		m_array = 0;
		m_nrows = m_ncols = 0;
	}

	HPM2D(muint nr, muint nc)
	{
		m_array = 0;
		resize(nr,nc);
	}

	HPM2D( const HPM2D &rhs )
	{
		m_array = 0;
		m_nrows = m_ncols = 0;
		copy(rhs);
	}

	HPM2D &operator=( const HPM2D &rhs )
	{
		copy(rhs);
		return *this;
	}

	~HPM2D()
	{
		if (m_array) delete [] m_array;
	}


	void copy( const HPM2D &rhs )
	{
		if (m_array != 0) delete [] m_array;
		m_nrows = rhs.m_nrows;
		m_ncols = rhs.m_ncols;
		m_array = new double[ m_nrows*m_ncols ];
		muint n = m_nrows*m_ncols;
		for ( muint i=0;i<n;i++ )
			m_array[i] = rhs.m_array[i];
	}


	void resize(muint nr, muint nc)
	{
		if (nr == 0 || nc == 0) return;

		if (m_array) delete [] m_array;
		m_array = new double[ nr * nc ];
		m_nrows = nr;
		m_ncols = nc;
		zero();
	}

	void zero()
	{
		muint len = m_nrows*m_ncols;
		for (muint i=0;i<len;i++)
			m_array[i] = 0.0;
	}

	inline muint nrows() const { return m_nrows; }
	inline muint ncols() const { return m_ncols; }
	inline double *data() { return m_array; }

	inline double &at(muint r, muint c)
	{
		return m_array[m_ncols*r+c];
	}

	inline const double &at(muint r, muint c) const
	{
		return m_array[m_ncols*r+c];
	}



private:

	double *m_array;
	muint m_nrows;
	muint m_ncols;
};
*/

#endif
