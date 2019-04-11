/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2010 Alistair Boyle <alistair.js.boyle@gmail.com>
 *
 *     This file is part of Meagre-Crowd.
 *
 *     Meagre-Crowd program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "config.h"

#include <stdlib.h>
//#include <bebop/smc/sparse_matrix.h>

// whether the first index is zero or one for sparse storage
enum matrix_base_t { FIRST_INDEX_ZERO = 0, FIRST_INDEX_ONE = 1 };

/// matrix storage format, either row or column based
// TODO block storage formats BCOO, BCSR, JAD
// TODO remove SM prefixes (bebop conflict)
enum matrix_format_t { 
  INVALID = 0, 
  DROW, ///< dense storage, no sparse optimizations, Dense (rows), store along rows first when serializing matrix, C uses this format naturally
  DCOL, ///< dense storage, no sparse optimizations, Dense (columns), store along columns first when serializing matrix, Fortran uses this format naturally
  SM_COO, ///< coordinate format, valid format request but will never be stored: will always set 'order' field as appropriate and use 'COMPRESSED' in the 'format' field
  SM_CSC, ///< compressed sparse column format, valid format request but will never be stored: will always set 'order' field as appropriate and use 'COMPRESSED' in the 'format' field
  SM_CSR ///< compressed sparse row format, valid format request but will never be stored: will always set 'order' field as appropriate and use 'COMPRESSED' in the 'format' field
};

/// matrix symmetry
//  TODO remove SM prefixes (bebop conflict)
enum matrix_symmetry_t { SM_UNSYMMETRIC = 0, 
  SM_SYMMETRIC, 
  SM_SKEW_SYMMETRIC, 
  SM_HERMITIAN 
}; 

/// where the data is stored (upper or lower triangular); if not MC_STORE_BOTH, then storage locations can be validated
enum matrix_symmetric_storage_t { MC_STORE_BOTH = 0, 
  UPPER_TRIANGULAR, 
  LOWER_TRIANGULAR 
};

/// data type; unpacked: A->i [A->p [j] ... A->p [j]+A->nz[j]-1] vs packed: A->i [A->p [j] ... A->p [j+1]-1]
// TODO I like the way TAUCS used a union to split out the different data types the pointer could hold (null, d, s, c, z)... look at taucs doc/
// TODO remove SM prefixes (bebop conflict)
// TODO support MATLAB complex format where real and imag are split ("ZOMPLEX") (CHOLMOD)
// TODO maybe split REAL/COMPLEX/PATTERN from storage type (float, double, int, long, uint, ulong)
// TODO support various index formats for ii, jj (int, unsigned int, long) (CHOLMOD)
// TODO support sorted OR unsorted flag (if sorted, doesn't need sorting later, mark as unsorted when modifying)
// TODO support "packed" -- nzmax is malloc size, nz is ptr to end-of-row/col (CHOLMOD)
enum matrix_data_type_t { 
  SM_REAL = 0,
  REAL_DOUBLE = 0, 
  REAL_SINGLE = 1, ///< C float, single-precision floating point number
  SM_COMPLEX = 2, 
  COMPLEX_DOUBLE = 2, 
  COMPLEX_SINGLE = 3, ///< two C floats make a complex single-precision floating point number
  SM_PATTERN = 4 };  ///< no data (dd), pattern of non-zero entries is indicated by ii, jj

///Note: initializing via 'matrix_t A = {0};' should get the most common defaults and a valid but empty matrix
typedef struct matrix_t {
  size_t m; ///< rows
  // TODO rename 'nrow'
  size_t n; ///< columns; m=n=0 indicates an empty/invalid matrix, if m=0, then n=0 and vice-versa, so you only need to check one of them
  // TODO rename 'ncol'
  size_t nz; ///< non-zeros, invalid for DENSE
  enum matrix_base_t   base;   ///< FIRST_INDEX_ZERO, FIRST_INDEX_ONE (default zero index, 'one' supports Fortran)
  enum matrix_format_t format; /*!< data storage (meaning varies by format)
   * DENSE: ii and jj are ignored
   * COO: ii=row indices, jj=column indices
   * CSR: ii=per-row ptrs into jj, jj=column indices
   * Note: ii[0]=0 always (w/ FIRST_INDEX_ZERO), ii[m]=nz always,
   * ii[1] points to the first non-zero entry of the second row in both dd and jj,
   * and there are ii[2]-ii[1] entries in the second row
   * CSC: ii=row indices, jj=per-column ptrs into ii
   * Note: swaps the row and column vs. CSR */

  enum matrix_symmetry_t sym;  ///< UNSYMMETRIC, SYMMETRIC, SKEW_SYMMETRIC, HERMITIAN
  enum matrix_symmetric_storage_t location; ///< MC_STORE_BOTH, UPPER_TRIANGULAR, LOWER_TRIANGULAR
  enum matrix_data_type_t data_type; ///< REAL, COMPLEX, etc
  void* dd; ///< data (size of entries defined by 'data_type', is float or double)
  // TODO support single precision (float) and double precision (double) in real and complex (paired) and complex (split)
  // TODO (this can be treated as seperate pointers by those who need to, though it can't be safely free'd as two pointers...)
  // TODO for complex (paired) use { x1, i y1, x2, i y2 }
  unsigned int* ii;
  // TODO support int/uint 8, 16, 32, 64 and conversion between (limited by max rows/columns values - fail if it could be exceeded), allow this to be upsized if the matrix grows, allow these to be converted to match a specific solver, minimize memory footprint for small matrices == faster solutions! (unless copying it around takes longer...)
  unsigned int* jj; ///< max is UINT_MAX (limits.h), at least 4,294,967,295
} matrix_t;

// helper functions to deallocate the components of a matrix_t or the whole thing, assuming it was malloc-ed
/*! \brief Allocates matrix.
 * \return Pointer to the allocated space as matrix_t*.
 */
matrix_t* malloc_matrix();

/*! \brief Frees allocated space of matrix m.
 */
void free_matrix( matrix_t* m );

/*! \brief I don't know what this does.
 */
void clear_matrix( matrix_t* m );

/*! \brief deep copy matrix
 */
// TODO const correctness
matrix_t* copy_matrix( matrix_t* m ); // deep copy

/*! \brief Compare matrices.
 * \return Zero on match.
 */
// TODO const correctness
int cmp_matrix( matrix_t* a, matrix_t* b );

/*! \brief Converts matrix to a new format.
 * \return Returns non-zero on failure.
 */
int convert_matrix( matrix_t* m, enum matrix_format_t f, enum matrix_base_t b );

/*! \brief This function.
 * \return Returns
 */
int convert_matrix_symmetry( 
    matrix_t* m /*! The matrix to convert.*/, 
    enum matrix_symmetric_storage_t loc /*! The matrix storage type of m. */
    );

/*! \brief
 * \return
 */
int detect_matrix_symmetry( matrix_t* m );

/*! \brief Checks if the matrix isn't malformed
 * \return 0: ok, <0: problem found
 */
// TODO const correctness
int validate_matrix( matrix_t* m );

/*! \brief For debugging, will print the matrix to stdout.
 * \return
 */
void printf_matrix( 
    char const *const pre /*! String at the start of each line, like indent and/or matrix name. */, 
    matrix_t* m /*! The matrix. */
    );

// from enum, returns width of ea. value in the matrix in bytes
size_t _data_width( const enum matrix_data_type_t t );


//int results_match( matrix_t* expected_matrix, matrix_t* result_matrix, const double precision );
/*! \details Test result of matrix computations.
 * \return 1: match w/in precision, 0: not matching
 */
// TODO refactor mv to matrix.h, operate on matrix_t objects
// TODO cmp_matrix()
int results_match(matrix_t* expected_matrix, matrix_t* result_matrix, const double precision);

#endif
