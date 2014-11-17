# -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
AC_DEFUN([DUNE_BOOST_PROPERTY_TREE],
[
        AC_REQUIRE([AC_PROG_CC])
        CPPFLAGS_SAVED="$CPPFLAGS"
        if test x"$with_boost" != x"no" ; then
          if test x"$with_boost" != x"yes" ; then
            if test x"$with_boost" != x"" ; then
              # USER-specified path
              BOOST_CPPFLAGS="$with_boost"
            fi
          fi
        fi
        CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
        export CPPFLAGS

        LDFLAGS_SAVED="$LDFLAGS"
        LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
        export LDFLAGS

        AC_CACHE_CHECK(whether the Boost.PropertyTree library is available,
                               dune_cv_boost_property_tree,
                               [AC_LANG_PUSH([C++])
                               AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[ @%:@include <boost/property_tree/ptree.hpp>
                                  ]],
                                  [[
                                      using boost::property_tree::ptree;
                                      ptree pt;
                                      return 0;
                                   ]])],
                               dune_cv_boost_property_tree=yes, dune_cv_boost_property_tree=no)
                               AC_LANG_POP([C++])])
         if test "x$dune_cv_boost_property_tree" = "xyes"; then
               AC_DEFINE(HAVE_BOOST_PROPERTY_TREE,,[define if the Boost.PropertyTree headers are available])
         fi
         CPPFLAGS="$CPPFLAGS_SAVED"
         LDFLAGS="$LDFLAGS_SAVED"
])
