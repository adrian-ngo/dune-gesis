#
# Some variables pointing to the locations of external libs:
#
# Paths to HFD5
set(HDF5_LIBS  -L${HOME}/Software/hdf5-1.8.13/lib -lhdf5)
set(HDF5_INCL  -I${HOME}/Software/hdf5-1.8.13/include -DH5_USE_16_API)
#
# Paths to FFTW3
set(FFTW3_LIBS  -L${HOME}/Software/fftw-3.3.4/lib -lfftw3_mpi -lfftw3)
set(FFTW3_INCL  -I${HOME}/Software/fftw-3.3.4/include)
#
#
set(BOOST_INCL  -I${HOME}/Software/boost_1_47_0)
#
#
##############################################################
#
# applications:
#
##############################################################

set(SUBDIRS )
set(noinst_PROGRAMS  dune_gesis gesis2dDG gesis3dDG gesis2dFEM gesis3dFEM)

# 2D DG version for YASP:
add_executable(gesis2dDG  geo_inversion_2d.cc)
add_dune_mpi_flags(gesis2dDG)
add_dune_superlu_flags(gesis2dDG)
add_dune_boost_flags(gesis2dDG)
set_property(TARGET gesis2dDG APPEND PROPERTY COMPILE_DEFINITIONS "USE_YASP")


set(gesis2dDG_SOURCES  geo_inversion_2d.cc)
set(gesis2dDG_CPPFLAGS  ${AM_CPPFLAGS} -DUSE_YASP
        ${HDF5_INCL}
        ${FFTW3_INCL}
        ${BOOST_INCL}
        ${PARMETIS_CPPFLAGS}
        ${SUPERLU_CPPFLAGS}
        ${DUNEMPICPPFLAGS})
set(gesis2dDG_LDADD
        ${DUNEMPILDFLAGS}
        ${DUNEMPILIBS}
        ${HDF5_LIBS}
        ${FFTW3_LIBS}
        ${PARMETIS_LIBS}
        ${SUPERLU_LIBS}
        ${DUNE_LDFLAGS} ${DUNE_LIBS}
        ${LDADD})
set(gesis2dDG_LDFLAGS  ${AM_LDFLAGS}
        ${SUPERLU_LDFLAGS}
        ${DUNE_LDFLAGS})





# 3D DG version for YASP:
set(gesis3dDG_SOURCES  geo_inversion_3d.cc)
set(gesis3dDG_CPPFLAGS  ${AM_CPPFLAGS} -DUSE_YASP
        ${HDF5_INCL}
        ${FFTW3_INCL}
        ${BOOST_INCL}
        ${PARMETIS_CPPFLAGS}
        ${SUPERLU_CPPFLAGS}
        ${DUNEMPICPPFLAGS})
set(gesis3dDG_LDADD
        ${DUNEMPILDFLAGS}
        ${DUNEMPILIBS}
        ${HDF5_LIBS}
        ${FFTW3_LIBS}
        ${PARMETIS_LIBS}
        ${SUPERLU_LIBS}
        ${DUNE_LDFLAGS} ${DUNE_LIBS}
        ${LDADD})
set(gesis3dDG_LDFLAGS  ${AM_LDFLAGS}
        ${SUPERLU_LDFLAGS}
        ${DUNE_LDFLAGS})
# 2D FEM version for YASP:
set(gesis2dFEM_SOURCES  geo_inversion_2d.cc)
set(gesis2dFEM_CPPFLAGS  ${AM_CPPFLAGS} -DUSE_YASP -DUSE_FEM
        ${HDF5_INCL}
        ${FFTW3_INCL}
        ${BOOST_INCL}
        ${PARMETIS_CPPFLAGS}
        ${SUPERLU_CPPFLAGS}
        ${DUNEMPICPPFLAGS})
set(gesis2dFEM_LDADD
        ${DUNEMPILDFLAGS}
        ${DUNEMPILIBS}
        ${HDF5_LIBS}
        ${FFTW3_LIBS}
        ${PARMETIS_LIBS}
        ${SUPERLU_LIBS}
        ${DUNE_LDFLAGS} ${DUNE_LIBS}
        ${LDADD})
set(gesis2dFEM_LDFLAGS  ${AM_LDFLAGS}
        ${SUPERLU_LDFLAGS}
        ${DUNE_LDFLAGS})
# 3D FEM version for YASP:
set(gesis3dFEM_SOURCES  geo_inversion_3d.cc)
set(gesis3dFEM_CPPFLAGS  ${AM_CPPFLAGS} -DUSE_YASP -DUSE_FEM
        ${HDF5_INCL}
        ${FFTW3_INCL}
        ${BOOST_INCL}
        ${PARMETIS_CPPFLAGS}
        ${SUPERLU_CPPFLAGS}
        ${DUNEMPICPPFLAGS})
set(gesis3dFEM_LDADD
        ${DUNEMPILDFLAGS}
        ${DUNEMPILIBS}
        ${HDF5_LIBS}
        ${FFTW3_LIBS}
        ${PARMETIS_LIBS}
        ${SUPERLU_LIBS}
        ${DUNE_LDFLAGS} ${DUNE_LIBS}
        ${LDADD})
set(gesis3dFEM_LDFLAGS  ${AM_LDFLAGS}
        ${SUPERLU_LDFLAGS}
        ${DUNE_LDFLAGS})
# default app:
add_executable(dune_gesis dune_gesis.cc)
add_dune_alberta_flags(dune_gesis )
add_dune_amiramesh_flags(dune_gesis )
add_dune_ug_flags(dune_gesis )
add_dune_mpi_flags(dune_gesis )
target_link_libraries(dune_gesis ${DUNE_LIBS})
add_dune_alugrid_flags(dune_gesis )
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dune_gesis_ LDADD
#        ${DUNE_LDFLAGS} ${DUNE_LIBS}
#        ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS}
#        ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS}
#        ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS}
#        ${UG_LDFLAGS} ${UG_LIBS}
#        ${DUNEMPILIBS}
#        ${LDADD})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dune_gesis_ LDFLAGS ${AM_LDFLAGS}
#        ${DUNEMPILDFLAGS}
#        ${UG_LDFLAGS}
#        ${AMIRAMESH_LDFLAGS}
#        ${ALBERTA_LDFLAGS}
#        ${ALUGRID_LDFLAGS}
#        ${DUNE_LDFLAGS})
##############################################################
##############################################################
# don't follow the full GNU-standard
# we need automake 1.9
set(AUTOMAKE_OPTIONS  foreign 1.9)
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-geometry=${DUNE_GEOMETRY_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-typetree=${DUNE_TYPETREE_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")
set(EXTRA_DIST  CMakeLists.txt)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
