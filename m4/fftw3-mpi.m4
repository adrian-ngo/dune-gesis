## -*- autoconf -*-
# $Id$
# searches for FFTW3-stuff

# DUNE_PATH_FFTW3_MPI()
#
# shell variables:
#   with_fftw3
#     no, yes
#   with_fftw3_libs
#     empty or something apropriate for LIBS
#   FFTW3_CPPFLAGS
#   FFTW3_LDFLAGS
#   FFTW3_LIBS
#   direct_FFTW3_CPPFLAGS
#   direct_FFTW3_LDFLAGS
#   direct_FFTW3_LIBS
#     same as above, but without variable indirections (e.g. the contents of
#     DUNEMPICPPFLAGS instead of '${DUNEMPICPPFLAGS}')
#   FFTW3_PARALLEL
#     1 or undef
#   HAVE_FFTW3
#     0 or 1
#
# substitutions:
#   FFTW3_CPPFLAGS
#   FFTW3_LDFLAGS
#   FFTW3_LIBS
#
# defines:
#   HAVE_FFTW3
#
# conditionals:
#   FFTW3
AC_DEFUN([DUNE_PATH_FFTW3_MPI],[
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PATH_XTRA])
  AC_REQUIRE([DUNE_MPI])

  AC_ARG_WITH(fftw3,
    AC_HELP_STRING([--with-fftw3=PATH],[directory with FFTW3 inside]),
    # expand tilde / other stuff
    eval with_fftw3=$with_fftw3
  )

  AC_ARG_WITH(fftw_libs,
    [AC_HELP_STRING([--with-fftw3-libs=LIBS],[additional libraries needed to link fftw programs. Those might be needed if your fftw library is static. Possible values are: -lz or -lz -lsz.])],[])

 # store values
 ac_save_CPPFLAGS="$CPPFLAGS"
 ac_save_LDFLAGS="$LDFLAGS"
 ac_save_LIBS="$LIBS"

 # start building variables

 # use special FFTW3-lib-path if it's set
 if test x$with_fftw3 != x ; then
   if test x"$with_fftw3" = xyes; then
     # Search in default locations
     with_fftw3="/usr"
   fi
   # extract absolute path
   if test -d $with_fftw3; then
     eval with_fftw3=`cd $with_fftw3 ; pwd`
   else
     AC_MSG_ERROR([FFTW3-directory $with_fftw3 does not exist])
   fi
   _dune_fftw3_libpath="-L$with_fftw3/lib"
   _dune_fftw3_incpath="-I$with_fftw3/include"
 else
   _dune_fftw3_libpath="-L/usr/lib"
   _dune_fftw3_incpath="-I/usr/include"
 fi

 CPPFLAGS="$CPPFLAGS $DUNEMPICPPFLAGS $_dune_fftw3_incpath"

 direct_FFTW3_CPPFLAGS="$DUNEMPICPPFLAGS $_dune_fftw3_incpath"
 nodep_FFTW3_CPPFLAGS="$DUNEMPICPPFLAGS $_dune_fftw3_incpath"
 FFTW3_CPPFLAGS="$DUNEMPICPPFLAGS $_dune_fftw3_incpath"
 direct_FFTW3_LDFLAGS="$DUNEMPILDFLAGS"
 nodep_FFTW3_LDFLAGS="$DUNEMPILDFLAGS"
 FFTW3_LDFLAGS="$DUNEMPILDFLAGS"
 direct_FFTW3_LIBS="$DUNEMPILIBS $_dune_fftw3_libpath"
 nodep_FFTW3_LIBS="$DUNEMPILIBS $_dune_fftw3_libpath"
 FFTW3_LIBS="$DUNEMPILIBS $_dune_fftw3_libpath"
 FFTW3_PARALLEL=0

 # test for an arbitrary header
 AC_CHECK_HEADER([fftw3-mpi.h],
   [HAVE_FFTW3=1],
   [HAVE_FFTW3=0])

 # Just for the configure check.  In the end, -L has to go into LIBS.
 LDFLAGS="$LDFLAGS $_dune_fftw3_libpath"
 # test for lib
 if test x$HAVE_FFTW3 = x1 ; then
   AC_CHECK_LIB([fftw3], [fftw_execute],
     [
       direct_FFTW3_LIBS="$_dune_fftw3_libpath -lfftw3_mpi -lfftw3 $with_fftw3_libs $direct_FFTW3_LIBS"
       nodep_FFTW3_LIBS="$_dune_fftw3_libpath -lfftw3_mpi -lfftw3 $with_fftw3_libs $nodep_FFTW3_LIBS"
       FFTW3_LIBS="$_dune_fftw3_libpath -lfftw3_mpi -lfftw3 $with_fftw3_libs $FFTW3_LIBS"
     ],
     [HAVE_FFTW3=0], ["$with_fftw3_libs"])
 fi

 # pre-set variable for summary
 with_fftw3="no"

 # did we succeed?
 if test x$HAVE_FFTW3 = x1 ; then
   AC_DEFINE(HAVE_FFTW3, 1, [Define to 1 if fftw3 was found])

   # proudly show in summary
   with_fftw3="yes"
 else
   # clear variables
   direct_FFTW3_CPPFLAGS=
   nodep_FFTW3_CPPFLAGS=
   FFTW3_CPPFLAGS=
   direct_FFTW3_LDFLAGS=
   nodep_FFTW3_LDFLAGS=
   FFTW3_LDFLAGS=
   direct_FFTW3_LIBS=
   nodep_FFTW3_LIBS=
   FFTW3_LIBS=
   FFTW3_PARALLEL=0
 fi

 AC_SUBST([FFTW3_CPPFLAGS])
 AC_SUBST([FFTW3_LDFLAGS])
 AC_SUBST([FFTW3_LIBS])

 # also tell automake
 AM_CONDITIONAL(FFTW3, test x$HAVE_FFTW3 = x1)

 # add to global list
 DUNE_ADD_ALL_PKG([FFTW3], [$nodep_FFTW3_CPPFLAGS],
                  [$nodep_FFTW3_LDFLAGS], [$nodep_FFTW3_LIBS])

 # reset values
 LIBS="$DUNEMPILIBS $ac_save_LIBS"
 LDFLAGS="$DUNEMPILDFLAGS $ac_save_LDFLAGS"
 CPPFLAGS="$DUNEMPICPPFLAGS $ac_save_CPPFLAGS"

 DUNE_ADD_SUMMARY_ENTRY([FFTW3],[$with_fftw3])

])
