AC_INIT([lib], [0.0.2])
AC_CONFIG_SRCDIR([swarm/print.c])

AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

AC_CANONICAL_HOST

AC_PROG_CC
AC_PROG_F77
AC_PROG_FC
AC_PROG_INSTALL
AC_PROG_LIBTOOL

AM_PROG_CC_C_O

AC_C_INLINE

if test "$ac_cv_c_inline"; then
  AC_DEFINE(HAVE_INLINE, 1)
  AC_SUBST(HAVE_INLINE)
fi

AC_SUBST(DATAHOME)

AC_CONFIG_FILES([            \
  Makefile                   \
  backus/Makefile            \
  champ/Makefile             \
  cond/Makefile              \
  dmsp/Makefile              \
  efi/Makefile               \
  emm/Makefile               \
  estist/Makefile            \
  gd/Makefile                \
  green/Makefile             \
  grid/Makefile              \
  grobs/Makefile             \
  hwm/Makefile               \
  invert_eej/Makefile        \
  iri/Makefile               \
  julia/Makefile             \
  julia/julia.h              \
  lapack_wrapper/Makefile    \
  lis/Makefile               \
  lls/Makefile               \
  magfit/Makefile            \
  mfield/Makefile            \
  mserv/Makefile             \
  msis/Makefile              \
  pej/Makefile               \
  poltor/Makefile            \
  superlu/Makefile           \
  swarm/Makefile             \
  pca/Makefile               \
  pca/pca.h                  \
  pca3d/Makefile             \
  pca3d/pca3d.h              \
])

AC_OUTPUT
