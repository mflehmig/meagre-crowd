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
/* Example program using the C interface to the
 * double precision version of MUMPS, dmumps_c.
 * We solve the system A x = RHS with
 * A = diag(1 2) and RHS = [1 4]ˆT
 * Solution is [1 2]ˆT */


#define TIMEON 1

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
#include "perftimer.h"
#include "file.h"
#include "matrix.h"
#include "solvers.h"
#include "util.h"


int main(int argc, char ** argv)
{
  int retval;
  const int extra_timing = 0;  // allow extra timing: initialization, matrix loading, etc.

  // handle command-line arguments
  struct parse_args* args = calloc(1, sizeof(struct parse_args));
  // TODO should default to appropriate epsilon for solver, may need to be *2 or some larger value given numerical instability? should print out epsilon of solution in verbose mode
  args->expected_precision = 5e-14;  // TODO was DBL_EPSILON=1.11e-16 but not stored with enough digits? // default to machine epsilon for 'double'
  // args->expected_precision = FLT_EPSILON*2;

  assert(args != NULL);  // calloc failure
  if ((retval = parse_args(argc, argv, args)) != EXIT_SUCCESS) {
    free(args);
    return retval;
  }

  // initialize MPI
  perftimer_t* timer = perftimer_malloc();
  if (extra_timing && args->rep == 0) {
    perftimer_inc(timer, "initialization", -1);
    perftimer_adjust_depth(timer, +1);
    perftimer_inc(timer, "MPI init", -1);
  }

  // if this solver is not thread/mpi safe, report an error if its been launched that way
  const int uses_mpi = solver_uses_mpi(args->solver);
  const int requires_mpi = solver_requires_mpi(args->solver);
  const int uses_omp = solver_uses_omp(args->solver);
  const int requires_omp = solver_requires_omp(args->solver);
  int c_mpi = get_mpi_num_procs();
  //  {
  //    const char* mpi_world_size = getenv("OMPI_COMM_WORLD_SIZE");
  //    if (mpi_world_size == NULL) {  // no env value configured
  //      c_mpi = 0;
  //      // printf( "no OMPI_COMM_WORLD_SIZE\n" );
  //    }
  //    else {  // works #ifdef OPEN_MPI -- we're using openMPI rather than lam or mpich...
  //      int ret = sscanf(mpi_world_size, "%d", &c_mpi);
  //      assert(ret == 1);
  //      // printf( "OMPI_COMM_WORLD_SIZE=%d\n", c_mpi );
  //    }
  //  }
  int c_omp = get_omp_num_threads();

  const int is_mpi = (requires_mpi || (uses_mpi && (c_mpi != 0)));
  const int is_omp = (requires_omp || (uses_omp && (c_omp != 0)));
  const int is_single_threaded = (!is_mpi && !is_omp);

  // start up MPI, if we're using it
  const int is_single_threaded_but_mpi = is_single_threaded && (c_mpi > 1);
  if (is_mpi || is_single_threaded_but_mpi) {
    int ierr;
    int provided_threading;
    ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided_threading);
    assert(ierr == 0);
    // assert(provided_threading == MPI_THREAD_MULTIPLE);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &(args->mpi_rank));
    assert(ierr == 0);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &c_mpi);
    assert(ierr == 0);
  }

  // wait till after firing up MPI, so we only put out an error on the rank=0 machine
  if (is_single_threaded_but_mpi) {
    if (args->mpi_rank == 0)
      fprintf( stderr, "error: selected solver (%s) is single threaded but was launched with %d threads\n",
              solver2str(args->solver), c_mpi);
    int ierr = MPI_Finalize();
    assert(ierr == 0);
    return 10;
  }

  // Define the problem on the host
  matrix_t* b = NULL;
  matrix_t* expected = NULL;
  matrix_t* A = NULL;
  matrix_t* rhs = NULL;
  if (args->mpi_rank == 0) {
    // we only load matrices for the zero-rank master process
    b = malloc_matrix();
    expected = malloc_matrix();
    A = malloc_matrix();
    rhs = malloc_matrix();

    if (extra_timing && args->rep == 0) {
      perftimer_adjust_depth(timer, -1);
      perftimer_inc(timer, "input", -1);
      perftimer_adjust_depth(timer, +1);
    }

    if (extra_timing && args->rep == 0)
      perftimer_inc(timer, "load", -1);

    // Load A
    if ((retval = load_matrix(args->input, A)) != 0) {
      return 1;
    }
    assert(validate_matrix(A) == 0);
    unsigned int m = matrix_rows(A); // rows
    // TODO load a rhs from the --input matrix file

    // Load b if provided, else create "random" rhs.
    // allocate an sequentially numbered right-hand side of A.m rows
    // TODO warn if there is already a rhs loaded
    if (args->rhs != NULL) {
      if ((retval = load_matrix(args->rhs, b)) != 0) {
        return 1;
      }
      assert(b->m == m);  // rows must match // TODO nice error (user could load some random matrix file, also testcases)

      // TODO but don't remove this here, since we want to test dense-matrix performance? OR do we need a new option! OR do we decide based on matrix sparsity? (EIT is sparse...)
    }
    else {
      assert(b != NULL);  // malloc failure
      *b = ( matrix_t ) { 0 };
      b->m = m;
      b->n = 1;
      b->nz = m;
      b->format = DCOL;
      b->dd = malloc(m * sizeof(double));
      b->data_type = REAL_DOUBLE;
      {  // initialize right-hand-side (b)
        double* d = b->dd;

        for (unsigned int i = 0; i < m; i++) {
            *d = i;
	          d++;
        }
      }
    }
    assert(validate_matrix(b) == 0);

    // Load x (if provided)
    if (args->expected != NULL) {
      // TODO refactor: this is a cut and paste of the loader for 'b'
      if ((retval = load_matrix(args->expected, expected)) != 0) {
        return 1;
      }
      assert(validate_matrix(expected) == 0);
      assert(expected->m == m);  // rows must match // TODO nice error (user could load some random matrix file, also testcases)
    }

    if (extra_timing && args->rep == 0) {
      perftimer_inc(timer, "solver", -1);
      perftimer_inc(timer, "rhs", -1);
    }

    // verbose output
    if (args->verbosity >= 1)
      print_verbose_output(args, A, b, expected, c_mpi, c_omp);

    if (extra_timing && args->rep == 0)
      perftimer_adjust_depth(timer, -1);
    //destroy_sparse_matrix (A); // TODO can't release it unless we're copying it...

  } // MPI master

  solver_state_t* state = solver_init(args->solver, args->verbosity, args->mpi_rank, timer);  //okay

  // MS: Since I don't understand how this timer works, and I am currently only interested in "time to solution" for
  // every solve call within this loop, I refactored it to a for loop (see below).
//    unsigned int r = 0;
//    do {
//      if (r != 0)
//        perftimer_restart(&timer);
//
//  #if TIMEON
//      long long start = current_timestamp();
//  #endif
//      solver_solve(state, A, b, rhs);  // TODO rhs -> x
//  #if TIMEON
//      long long finish = current_timestamp();
//      printf("Iteration %d\t time %lld\n", r, (finish-start));
//  #endif
//
//      r++;
//    }
//    while (r < args->rep);

  printf("\nStart solving Ax=b %d time(s) ...\n", args->rep);
  // for (int i = 1; i <= args->rep; ++i) { // BOYLE
  for (unsigned int i = 1; i <= args->rep; ++i) {
#if TIMEON
    long long start = current_timestamp();
#endif
    solver_solve(state, A, b, rhs);  // TODO rhs -> x
#if TIMEON
    long long finish = current_timestamp();
    printf("Iteration %d\t time %lld\n", i, (finish-start));
#endif
  }
  printf("Done.\n");

  if (extra_timing && args->rep == 0) {
    perftimer_inc(timer, "clean up", -1);
    perftimer_adjust_depth(timer, -1);
    perftimer_inc(timer, "output", -1);
  }

  // print/store result
  if ((args->mpi_rank == 0) && (args->output != NULL)) {
    if (strncmp(args->output, "-", 2) == 0) {
      printf_matrix("  x", rhs);
    }
    else {
      // TODO test args->output to decide if its an acceptable filename before now
      // TODO really, we shouldn't care what the file name is, just the file format...
      if (save_matrix(rhs, args->output) != 0)
        retval = 12;  // TODO normalize error codes so they are defined in one place (defines?)
    }
  }

  // test result?
  if ((args->mpi_rank == 0) && (expected->format != INVALID)) {
    printf("\nCheck result against expected result...\n");
    perftimer_inc(timer, "test", -1);
    if (results_match(expected, rhs, args->expected_precision)) {
      retval = 0;
      printf("  PASS.\n");
    }
    else {
      retval = 100;
      printf("  FAIL.\n");
    }
  }

  solver_finalize(state);

  // we can report on memory usage per-process
  // RUSAGE_SELF includes the usage for all threads and children
  struct rusage usage;
  int rr = getrusage(RUSAGE_SELF, &usage);
  if ((rr == 0) && (args->verbosity >= 2)) {
    if (args->verbosity >= 4)
      printf(
          "user %ld:%ld, sys: %ld:%ld,  maxrss: %ld, ixrss: %ld, idrss: %ld, isrss: %ld minflt: %ld, majflt: %ld, nswap: %ld, inblock: %ld, outblock: %ld, msgsnd: %ld, msgrcv: %ld, signals: %ld, nvcsw: %ld, nivcsw: %ld\n",
          usage.ru_utime.tv_sec, usage.ru_utime.tv_usec, /* user time used */
          usage.ru_stime.tv_sec, usage.ru_stime.tv_usec, /* system time used */
          usage.ru_maxrss, /* maximum resident set size */
          usage.ru_ixrss, /* integral shared memory size */
          usage.ru_idrss, /* integral unshared data size */
          usage.ru_isrss, /* integral unshared stack size */
          usage.ru_minflt, /* page reclaims */
          usage.ru_majflt, /* page faults */
          usage.ru_nswap, /* swaps */
          usage.ru_inblock, /* block input operations */
          usage.ru_oublock, /* block output operations */
          usage.ru_msgsnd, /* messages sent */
          usage.ru_msgrcv, /* messages received */
          usage.ru_nsignals, /* signals received */
          usage.ru_nvcsw, /* voluntary context switches */
          usage.ru_nivcsw); /* involuntary context switches */
    if (args->verbosity >= 1)
      printf("memory usage: maximum resident=%lgMB, page reclaims=%ld, page faults w/ IO=%ld\n", usage.ru_maxrss / 1e3, /* maximum resident set size */
             usage.ru_minflt, /* page reclaims */
             usage.ru_majflt); /* page faults */
  }

  // show timing info, if requested, to depth N
  if (extra_timing && args->rep == 0) {
    perftimer_adjust_depth(timer, -1);
    perftimer_inc(timer, "finished", -1);
  }

  // collect memory total
  double mem_sum;
  if (is_mpi) {
    double mem_ind = usage.ru_maxrss / 1e3;
    MPI_Op sum_f;
    int ret;
    ret = MPI_Op_create(&mpi_sum, 1, &sum_f);  // addition is a commutative operation (order doesn't matter)
    assert(ret == MPI_SUCCESS);
    ret = MPI_Reduce(&mem_ind, &mem_sum, 1, MPI_DOUBLE, sum_f, 0, MPI_COMM_WORLD);
    assert(ret == MPI_SUCCESS);
  }
  else {
    mem_sum = usage.ru_maxrss / 1e3;
  }

  if (args->mpi_rank == 0) {
    if (args->timing_enabled == 1) {
      printf("status, solver, threads, rows, max. memory (MB), total memory (MB), ");
      perftimer_printf_csv_header(timer, 2);
      printf("%s, %s, %d, %zu, %lg, %lg, ", retval == 100 ? "FAIL" : "PASS", solver2str(args->solver), c_mpi, A->m,
             usage.ru_maxrss / 1e3, mem_sum);
      perftimer_printf_csv_body(timer, 2);
    }
    else if (args->timing_enabled != 0) {
      perftimer_printf(timer, args->timing_enabled - 2);
    }
  }

  // close down MPI
  if (extra_timing && args->rep == 0) {
    perftimer_inc(timer, "clean up", -1);
    perftimer_adjust_depth(timer, +1);

    perftimer_inc(timer, "MPI", -1);
  }

  if (is_mpi) {
    int ierr = MPI_Finalize();
    assert(ierr == 0);
  }

  // clean up matrices
  free_matrix(b);
  free_matrix(expected);
  free_matrix(A);
  free_matrix(rhs);

  perftimer_free(timer);
  free(args);
  return retval;
}
