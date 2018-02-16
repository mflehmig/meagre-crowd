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
#ifndef _ARGS_H_
#define _ARGS_H_

#include <argp.h>
#include "config.h"
#include "solvers.h"

/** \brief Structure holding the parsed command line arguments.
 *
 */
struct parse_args
{
  char* input;                    ///< Input matrix A
  char* output;                   ///< output solution vector x
  char* rhs;                      ///< right-hand side b
  char* expected;                 ///< Expected solution vector x to compare solution from meagre-crowd against
  double expected_precision;      ///< Expected precision of solution x
  unsigned int timing_enabled;    ///< Enable timing functionality
  unsigned int verbosity;         ///< Verbosity
  unsigned int rep;               ///< Number of repetitions to solve the system
  int mpi_rank;                   ///< Set by meagre-crowd
  int solver;                     ///< Solver to use
};

int parse_args(int argc, char** argv, struct parse_args* args);

#endif
