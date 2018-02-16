/*
 * util.c
 *
 *  Created on: 16 Feb 2018
 *      Author: mf
 */

#include "util.h"

//void mpi_sum(void* in, void* inout, int *len, MPI_Datatype *dptr);
void mpi_sum(void* in, void* inout, int *len, MPI_Datatype *dptr)
{
  assert(*dptr == MPI_DOUBLE);
  double* const din = in;
  double* const dinout = inout;
  for (int i = 0; i < *len; i++) {
    dinout[i] += din[i];
  }
}

/// Time in micro seconds.
long long current_timestamp()
{
    struct timeval te;
    gettimeofday(&te, NULL);
    long long usecs = te.tv_sec * 1000000LL + te.tv_usec;
    return usecs;
}

int get_mpi_num_procs()
{
  int num_procs;
  const char* mpi_procs = getenv("OMPI_COMM_WORLD_SIZE");
  if (mpi_procs == NULL) {  // no env value configured
    num_procs = 1;
  } else {
    num_procs = atoi(mpi_procs);
  }

  return num_procs;
}

/** \brief Get the number of OpenMP threads to use.
 *
 * If the environment variable OMP_NUM_THREADS is specified (by the user), we take this value. If it is not defined, we
 * check if Meagre-Crowd was compiled with OpenMP support and take the value provided by omp_get_max_threads().
 * Otherwise, we return 1.
 *
 * \remark This (read-only) value does not influence the solver behaviour!
 */
int get_omp_num_threads()
{
  int num_threads;
  const char* omp_threads = getenv("OMP_NUM_THREADS");
  if (omp_threads == NULL) {  // no env value configured
    printf("Env. variable OMP_NUM_THREADS is not defined.\n");
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
#else
    num_threads = 1;
#endif
  } else {
    num_threads = atoi(omp_threads);
  }

  return num_threads;
}

/** \brief Print configuration
 *
 */
void print_verbose_output(struct parse_args* args, matrix_t* A, matrix_t* b, matrix_t* expected, int c_mpi, int c_omp) {
  const int false = 0;
  assert(A != NULL);
  int ierr = convert_matrix(A, SM_COO, FIRST_INDEX_ZERO);
  assert(ierr == 0);
  const char* sym, *location, *type;
  switch (A->sym) {
    case SM_UNSYMMETRIC:
      sym = "unsymmetric";
      break;
    case SM_SYMMETRIC:
      sym = "symmetric";
      break;
    case SM_SKEW_SYMMETRIC:
      sym = "skew symmetric";
      break;
    case SM_HERMITIAN:
      sym = "hermitian";
      break;
    default:
      assert(false);  // fell through
  }
  if ((args->verbosity < 2) || (A->sym == SM_UNSYMMETRIC)) {
    location = "";
  }
  else {
    switch (A->location) {
      case UPPER_TRIANGULAR:
        location = " (upper)";
        break;
      case LOWER_TRIANGULAR:
        location = " (lower)";
        break;
      case MC_STORE_BOTH:
        location = "";
        break;  // nothing
      default:
        assert(false);
    }
  }
  switch (A->data_type) {
    case REAL_DOUBLE:
      type = "real";
      break;
    case REAL_SINGLE:
      type = "real (single-precision)";
      break;
    case COMPLEX_DOUBLE:
      type = "complex";
      break;
    case COMPLEX_SINGLE:
      type = "complex (single-precision)";
      break;
    case SM_PATTERN:
      type = "pattern";
      break;
    default:
      assert(false);  // fell through
  }

//      printf("Ax=b: A is %zux%zu, nz=%zu, %s%s, %s, b is %zux%zu, nz=%zu\nsolved with %s on %d core%s, %d thread%s\n",
//             A->m, A->n, A->nz, sym, location, type, b->m, b->n, b->nz, solver2str(args->solver), c_mpi,
//             c_mpi == 1 ? "" : "s", c_omp, c_omp == 1 ? "" : "s");
  printf("\n====================== Configuration ======================\n");
  printf("                A: %zu x %zu, nz=%zu, %s%s, %s\n", A->m, A->n, A->nz, sym, location, type);
  printf("                b: %zu x %zu, nz=%zu\n", b->m, b->n, b->nz);
  if (expected->format != INVALID)
  {
    printf("          Known x: %s\n", args->expected);
    printf("  Comp. precision: %e\n", args->expected_precision);
  }
  if (NULL != args->output)
  {
    printf("       Output (x): %s\n", args->output);
  }
  printf("           Solver: %s \n", solver2str(args->solver));
  printf("    Cores/Threads: %d / %d\n", c_mpi, c_omp);
  printf("      Repetitions: %d\n", args->rep);
  printf("===========================================================\n\n");

  if (args->verbosity >= 2) {
    printf_matrix("  A", A);
  }

  if (args->verbosity >= 2) {  // show the rhs matrix
    printf_matrix("  b", b);
  }
}



