AC_DEFUN([ACXX_MPI], [
  AC_ARG_VAR(MPICPP,[MPI CPP compiler command])
  AC_CHECK_PROGS(MPICPP, mpiccxx, $CPP)
  acxx_mpi_save_CPP="$CPP"
  CPP="$MPICPP"
  AC_SUBST(MPICPP)
])
