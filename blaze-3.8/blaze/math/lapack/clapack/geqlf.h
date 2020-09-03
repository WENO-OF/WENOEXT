//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/geqlf.h
//  \brief Header file for the CLAPACK geqlf wrapper functions
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_GEQLF_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_GEQLF_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/Types.h>
#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>


//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if !defined(INTEL_MKL_VERSION)
extern "C" {

void sgeqlf_( blaze::blas_int_t* m, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
              float* tau, float* work, blaze::blas_int_t* lwork, blaze::blas_int_t* info );
void dgeqlf_( blaze::blas_int_t* m, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
              double* tau, double* work, blaze::blas_int_t* lwork, blaze::blas_int_t* info );
void cgeqlf_( blaze::blas_int_t* m, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
              float* tau, float* work, blaze::blas_int_t* lwork, blaze::blas_int_t* info );
void zgeqlf_( blaze::blas_int_t* m, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
              double* tau, double* work, blaze::blas_int_t* lwork, blaze::blas_int_t* info );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK QL DECOMPOSITION FUNCTIONS (GEQLF)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK QL decomposition functions (geqlf) */
//@{
void geqlf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda,
            float* tau, float* work, blas_int_t lwork, blas_int_t* info );

void geqlf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda,
            double* tau, double* work, blas_int_t lwork, blas_int_t* info );

void geqlf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda,
            complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info );

void geqlf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda,
            complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QL decomposition of the given dense single precision column-major
//        matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix QL decomposition of a general \a m-by-\a n single
// precision column-major matrix based on the LAPACK sgeqlf() function. The resulting decomposition
// has the form

                              \f[ A = Q \cdot L, \f]

// where the \c Q is represented as a product of elementary reflectors

               \f[ Q = H(k) ... H(2) H(1) \texttt{, with k = min(m,n).} \f]

// Each H(i) has the form

                      \f[ H(i) = I - tau \cdot v \cdot v^T, \f]

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(m-k+i+1:m) = 0</tt> and
// <tt>v(m-k+i) = 1</tt>. <tt>v(1:m-k+i-1)</tt> is stored on exit in <tt>A(1:m-k+i-1,n-k+i)</tt>,
// and \c tau in \c tau(i). Thus in case \a m >= \a n, the lower triangle of the subarray
// A(m-n+1:m,1:n) contains the \a n-by-\a n lower triangular matrix \c L and in case \a m <= \a n,
// the elements on and below the (\a n-\a m)-th subdiagonal contain the \a m-by-\a n lower
// trapezoidal matrix \c L; the remaining elements in combination with the array \c tau represent
// the orthogonal matrix \c Q as a product of min(\a m,\a n) elementary reflectors.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the sgeqlf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void geqlf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda,
                   float* tau, float* work, blas_int_t lwork, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   sgeqlf_( &m, &n, A, &lda, tau, work, &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QL decomposition of the given dense double precision column-major
//        matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix QL decomposition of a general \a m-by-\a n double
// precision column-major matrix based on the LAPACK sgeqlf() function. The resulting decomposition
// has the form

                              \f[ A = Q \cdot L, \f]

// where the \c Q is represented as a product of elementary reflectors

               \f[ Q = H(k) ... H(2) H(1) \texttt{, with k = min(m,n).} \f]

// Each H(i) has the form

                      \f[ H(i) = I - tau \cdot v \cdot v^T, \f]

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(m-k+i+1:m) = 0</tt> and
// <tt>v(m-k+i) = 1</tt>. <tt>v(1:m-k+i-1)</tt> is stored on exit in <tt>A(1:m-k+i-1,n-k+i)</tt>,
// and \c tau in \c tau(i). Thus in case \a m >= \a n, the lower triangle of the subarray
// A(m-n+1:m,1:n) contains the \a n-by-\a n lower triangular matrix \c L and in case \a m <= \a n,
// the elements on and below the (\a n-\a m)-th subdiagonal contain the \a m-by-\a n lower
// trapezoidal matrix \c L; the remaining elements in combination with the array \c tau represent
// the orthogonal matrix \c Q as a product of min(\a m,\a n) elementary reflectors.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the sgeqlf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void geqlf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda,
                   double* tau, double* work, blas_int_t lwork, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   dgeqlf_( &m, &n, A, &lda, tau, work, &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QL decomposition of the given dense single precision complex
//        column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix QL decomposition of a general \a m-by-\a n single
// precision complex column-major matrix based on the LAPACK sgeqlf() function. The resulting
// decomposition has the form

                              \f[ A = Q \cdot L, \f]

// where the \c Q is represented as a product of elementary reflectors

               \f[ Q = H(k) ... H(2) H(1) \texttt{, with k = min(m,n).} \f]

// Each H(i) has the form

                      \f[ H(i) = I - tau \cdot v \cdot v^T, \f]

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(m-k+i+1:m) = 0</tt> and
// <tt>v(m-k+i) = 1</tt>. <tt>v(1:m-k+i-1)</tt> is stored on exit in <tt>A(1:m-k+i-1,n-k+i)</tt>,
// and \c tau in \c tau(i). Thus in case \a m >= \a n, the lower triangle of the subarray
// A(m-n+1:m,1:n) contains the \a n-by-\a n lower triangular matrix \c L and in case \a m <= \a n,
// the elements on and below the (\a n-\a m)-th subdiagonal contain the \a m-by-\a n lower
// trapezoidal matrix \c L; the remaining elements in combination with the array \c tau represent
// the orthogonal matrix \c Q as a product of min(\a m,\a n) elementary reflectors.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the sgeqlf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void geqlf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda,
                   complex<float>* tau, complex<float>* work, blas_int_t lwork, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   cgeqlf_( &m, &n, reinterpret_cast<ET*>( A ), &lda, reinterpret_cast<ET*>( tau ),
            reinterpret_cast<ET*>( work ), &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QL decomposition of the given dense double precision complex
//        column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix QL decomposition of a general \a m-by-\a n double
// precision complex column-major matrix based on the LAPACK sgeqlf() function. The resulting
// decomposition has the form

                              \f[ A = Q \cdot L, \f]

// where the \c Q is represented as a product of elementary reflectors

               \f[ Q = H(k) ... H(2) H(1) \texttt{, with k = min(m,n).} \f]

// Each H(i) has the form

                      \f[ H(i) = I - tau \cdot v \cdot v^T, \f]

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(m-k+i+1:m) = 0</tt> and
// <tt>v(m-k+i) = 1</tt>. <tt>v(1:m-k+i-1)</tt> is stored on exit in <tt>A(1:m-k+i-1,n-k+i)</tt>,
// and \c tau in \c tau(i). Thus in case \a m >= \a n, the lower triangle of the subarray
// A(m-n+1:m,1:n) contains the \a n-by-\a n lower triangular matrix \c L and in case \a m <= \a n,
// the elements on and below the (\a n-\a m)-th subdiagonal contain the \a m-by-\a n lower
// trapezoidal matrix \c L; the remaining elements in combination with the array \c tau represent
// the orthogonal matrix \c Q as a product of min(\a m,\a n) elementary reflectors.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the sgeqlf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void geqlf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda,
                   complex<double>* tau, complex<double>* work, blas_int_t lwork, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zgeqlf_( &m, &n, reinterpret_cast<ET*>( A ), &lda, reinterpret_cast<ET*>( tau ),
            reinterpret_cast<ET*>( work ), &lwork, info );
}
//*************************************************************************************************

} // namespace blaze

#endif