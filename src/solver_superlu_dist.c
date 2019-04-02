/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2011 Alistair Boyle <alistair.js.boyle@gmail.com>
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

#include "solver_superlu_dist.h"
#include "solvers.h"
#include "matrix.h"
#include "matrix_share.h"
#include "util.h"
#include <superlu_ddefs.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>
#include <stdint.h> // int64_t
#include <unistd.h> // dup -> redirect stdout temporarily

#include <mpi.h>
#include <math.h>

typedef struct {
  superlu_dist_options_t options;
  gridinfo_t grid;
  int active; // is this node active in the superlu grid?
  int rank0; // who is the old rank0 in the new communicator
  SuperMatrix A;
  ScalePermstruct_t scale_permute;
  LUstruct_t lu;
} solve_system_superlu_dist_t;

void solver_init_superlu_dist( solver_state_t* s ) {
  assert( s != NULL );
  solve_system_superlu_dist_t * const p = calloc( 1, sizeof( solve_system_superlu_dist_t ) );
  assert( p != NULL );
  s->specific = p;

  // set default options
  set_default_options_dist(&(p->options));
  // Get ordering from environment variable ORDERING
  p->options.ColPerm = getSuperLUOrdering();
  if((s->mpi_rank == 0) && (s->verbosity >= 3)) {
    p->options.PrintStat = YES;
    print_options_dist(&(p->options));
    print_sp_ienv_dist(&(p->options));
    printf("ISPEC for column permutation): %d\n", p->options.ColPerm);
  }
  else {
    p->options.PrintStat = NO;
  }

  // determine the size of the MPI communicator
  const MPI_Comm initial_comm = MPI_COMM_WORLD; // TODO pass in communicator instead of assuming MPI_COMM_WORLD
  int size;
  int ret = MPI_Comm_size(initial_comm, &size);
  assert(ret == MPI_SUCCESS);

  // determine size: a 2D grid
  int npcol = floor(sqrt(size));
  int nprow = floor(size/npcol);
  superlu_gridinit(initial_comm, nprow, npcol, &(p->grid));
  p->active = (p->grid.iam < nprow * npcol);
  if(s->mpi_rank == 0) {
    assert(p->active); // must have the rank=0 node active or we're broken (A is only loaded on rank=0)
    p->rank0 = p->grid.iam;
  }
  ret = MPI_Bcast(&(p->rank0), 1, MPI_INT, 0, initial_comm);
  assert(ret == MPI_SUCCESS);
}

// TODO split analyze stage into ordering and symbolic factorization stages?
void solver_analyze_superlu_dist( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  solve_system_superlu_dist_t* const p = s->specific;
  assert( p != NULL );
  if( !p->active ) // check if this grid node is active
    return;

  // setup the A matrix, shared globally
  matrix_t* AA;
  if(s->mpi_rank == 0) {
    AA = copy_matrix(A);
  }
  else {
    AA = malloc_matrix();
  }
  matrix_bcast(AA, p->rank0, p->grid.comm);
  dCreate_CompCol_Matrix_dist(&(p->A), AA->m, AA->n, AA->nz,
                              AA->dd, (int*) AA->ii, (int*) AA->jj, SLU_NC, SLU_D, SLU_GE);
  // last 3 enums are: stype=column-wise(no super-nodes), dtype=double, mtype=general);

  // clear the data pointers since these are now held by p->A and
  // release the rest of the AA matrix pointer
  AA->ii = NULL;
  AA->jj = NULL;
  AA->dd = NULL;
  free_matrix(AA);
  AA = NULL;

  // TODO analyze
}

void solver_factorize_superlu_dist( solver_state_t* s, matrix_t* A ) {
  // TODO factorize
}

void solver_evaluate_superlu_dist( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  solve_system_superlu_dist_t* const p = s->specific;
  assert( p != NULL );
  if( !p->active ) // check if this grid node is active
    return;

  // initialize structures
  SuperLUStat_t stat;
  ScalePermstructInit(p->A.nrow, p->A.ncol, &(p->scale_permute));
  LUstructInit(/*p->A.nrow,*/ p->A.ncol, &(p->lu));
  PStatInit(&stat);

  // setup for solver
  int ldb;
  int nrhs;
  if(s->mpi_rank == 0) {
    nrhs = b->n;
    ldb = b->m;
  }
  int ret;
  ret = MPI_Bcast(&nrhs, 1, MPI_INT, p->rank0, p->grid.comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&ldb, 1, MPI_INT, p->rank0, p->grid.comm);
  assert(ret == MPI_SUCCESS);
  // now we know the size of the rhs on all nodes, allocate space
  double* bb = doubleMalloc_dist(nrhs*ldb);
  double* berr = doubleMalloc_dist(nrhs*ldb);
  assert(bb != NULL);
  assert(berr != NULL);
  // share the rhs to all nodes
  if(s->mpi_rank == 0) {
    assert(b->format == DCOL);
    assert(b->data_type == REAL_DOUBLE);
    memcpy(bb, b->dd, sizeof(double)*nrhs*ldb);
  }
  ret = MPI_Bcast(bb, nrhs*ldb, MPI_DOUBLE, p->rank0, p->grid.comm); // TODO broken
  assert(ret == MPI_SUCCESS);

  // call solver
  int info_;
  pdgssvx_ABglobal(&(p->options), &(p->A), &(p->scale_permute), bb, ldb, nrhs, &(p->grid),
                   &(p->lu), berr, &stat, &info_);

  // making info an unsigned int
  unsigned int info;
  assert(info_ >= 0);
  info = (unsigned int)info_;
  // Compare info with A->ncol which is equal to b->nrow. Only master prints this information.
  // TODO If info != 0 the solver failed. Thus return error code and do not compare results.
  if(!s->mpi_rank) {
    if (info == 0)
      printf("pdgssvx_ABglobal() returns info %d which means FINE\n", info);
    if (info > 0 && info <= b->m)
      printf("pdgssvx_ABglobal() returns info %d which means U(%d,%d) is exactly zero\n", info, info, info);
    if (info > 0 && info > b->m)
      printf("pdgssvx_ABglobal() returns info %d which means memory allocation failed\n", info);
  }

  // TODO calculate inf_norm

  // print statistics
  if(s->verbosity >= 3)
    PStatPrint(&(p->options), &stat, &(p->grid));

  // release structures
  ScalePermstructFree(&(p->scale_permute));
  Destroy_LU(p->A.nrow, &(p->grid), &(p->lu));

  // copy the answer into x
  if(s->mpi_rank == 0) {
    clear_matrix(x);
    x->format = DCOL;
    x->m = ldb; // since the A matrix is square, the rows in b match the columns in A, which match the rows in x
    x->n = nrhs; // matches the b's columns
    x->nz = ldb * nrhs;
    x->dd = malloc(sizeof(double)*(x->m)*(x->n));
    assert(x->dd != NULL);
    memcpy(x->dd, bb, sizeof(double)*(x->m)*(x->n));
    assert(validate_matrix(x) == 0);
  }

  // release data
  PStatFree(&stat);
  SUPERLU_FREE(bb);
  SUPERLU_FREE(berr);
}

void solver_finalize_superlu_dist( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_superlu_dist_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {

    // release the A matrix
    if( p->active )
      Destroy_CompCol_Matrix_dist(&(p->A));

    // shutdown the MPI grid for superlu
    superlu_gridexit(&(p->grid));
  }
  free( p );
  s->specific = NULL;
}

colperm_t getSuperLUOrdering()
{
  colperm_t ret = NATURAL;  // Default value.
  printf("SuperLU_dist Ordering: ");
  const char* env_p = getenv("ORDERING");
  // if(const char* env_p = getenv("ORDERING")) {
  if (env_p != NULL)
  {
    // std::string env(env_p);
    //if (env.compare("NATURAL") == 0) {
    if (strcmp(env_p, "NATURAL") == 0)
    {
      printf("NATURAL\n");
      return NATURAL;
    }
    if (strcmp(env_p, "MMD_ATA") == 0)
    {
      printf("MMD_ATA\n");
      return MMD_ATA;
    }
    if (strcmp(env_p, "MMD_AT_PLUS_A") == 0)
    {
      printf("MMD_AT_PLUS_A\n");
      return MMD_AT_PLUS_A;
    }
    if (strcmp(env_p, "COLAMD") == 0)
    {
      printf("COLAMD\n");
      return COLAMD;
    }
    if (strcmp(env_p, "METIS_AT_PLUS_A") == 0)
    {
      printf("METIS_AT_PLUS_A\n");
      return METIS_AT_PLUS_A;
    }
    if (strcmp(env_p, "PARMETIS") == 0)
    {
      printf("PARMETIS\n");
      return PARMETIS;
    }
    // ZOLTAN is only defined in SuperLU, not in SuperLU_MT. But since the user
    // guide does not hold any information about ZOLTAN, we do not use it.
//    if (strcmp(env_p,"ZOLTAN") == 0) {
//      //printf( "ZOLTAN\n";
//      m_warning(E_NULL, "Hqp_IpMatrix::getPermcSpec ZOLTAN is not enabled. "
//                        "Set to COLAMD.");
//    } else if (env.compare("MY_PERMC") == 0) {
//      //printf( "MY_PERMC\n";
//      ret = MY_PERMC;
//    }
    printf("Info: Not a valid Ordering. Use COLAMD.\n");
  }
  printf("NATURAL\n");
  return ret;
}

