AC_DEFUN([ACXX_MPI], [
  AC_ARG_VAR(MPICXX,[MPI CXX compiler command])
  AC_CHECK_PROGS(MPICXX, mpicxx, $CPP)
  acxx_mpi_save_CXX="$CXX"
  CXX="$MPICXX"
  AC_SUBST(MPIXX)
])
