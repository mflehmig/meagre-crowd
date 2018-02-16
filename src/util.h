/*
 * Util.h
 *
 *  Created on: 16 Feb 2018
 *      Author: mf
 */

#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_


#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <libgen.h>
#include <string.h>
#include <float.h> // machine epsilon: LDBL_EPSILON, DBL_EPSILON, FLT_EPSILON

#include <mpi.h>
#include <omp.h>

// for getrusage
#include <sys/time.h>
#include <sys/resource.h>

#include "args.h"
//#include "perftimer.h"
//#include "file.h"
#include "matrix.h"
//#include "solvers.h"



//void mpi_sum(void* in, void* inout, int *len, MPI_Datatype *dptr);
void mpi_sum(void* in, void* inout, int *len, MPI_Datatype *dptr);

/// Time in micro seconds.
long long current_timestamp();

int get_mpi_num_procs();

/** \brief Get the number of OpenMP threads to use.
 *
 * If the environment variable OMP_NUM_THREADS is specified (by the user), we take this value. If it is not defined, we
 * check if Meagre-Crowd was compiled with OpenMP support and take the value provided by omp_get_max_threads().
 * Otherwise, we return 1.
 *
 * \remark This (read-only) value does not influence the solver behaviour!
 */
int get_omp_num_threads();

/** \brief Print configuration
 *
 */
void print_verbose_output(struct parse_args* args, matrix_t* A, matrix_t* b, matrix_t* expected, int c_mpi, int c_omp);

#endif /* SRC_UTIL_H_ */
