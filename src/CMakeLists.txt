#
# Some variables pointing to the locations of external libs:
#
# Paths to HFD5
#set(HDF5_LIBS  -L${HOME}/Software/hdf5-1.8.13/lib -lhdf5)
#set(HDF5_INCL  -I${HOME}/Software/hdf5-1.8.13/include -DH5_USE_16_API)
#
# Paths to FFTW3
#set(FFTW3_LIBS  -L${HOME}/Software/fftw-3.3.4/lib -lfftw3_mpi -lfftw3)
#set(FFTW3_INCL  -I${HOME}/Software/fftw-3.3.4/include)
#
#
#set(BOOST_INCL  -I${HOME}/Software/boost_1_47_0)
#
#
##############################################################
#
# applications:
#
##############################################################

# 2D DG version for YASP:
add_executable(gesis2dDG  geo_inversion_2d.cc)
set_property(TARGET gesis2dDG APPEND PROPERTY COMPILE_DEFINITIONS "USE_YASP")
target_link_libraries(gesis2dDG ${DUNE_LIBS})
add_dune_superlu_flags(gesis2dDG)
add_dune_boost_flags(gesis2dDG) 
add_dune_hdf5_flags(gesis2dDG)
add_dune_fftw3_flags(gesis2dDG)
add_dune_mpi_flags(gesis2dDG)


# 2D FEM version for YASP:
add_executable(gesis2dFEM  geo_inversion_2d.cc)
set_property(TARGET gesis2dFEM APPEND PROPERTY COMPILE_DEFINITIONS "USE_YASP")
set_property(TARGET gesis2dFEM APPEND PROPERTY COMPILE_DEFINITIONS "USE_FEM")
target_link_libraries(gesis2dFEM ${DUNE_LIBS})
add_dune_superlu_flags(gesis2dFEM)
add_dune_boost_flags(gesis2dFEM) 
add_dune_hdf5_flags(gesis2dFEM)
add_dune_fftw3_flags(gesis2dFEM)
add_dune_mpi_flags(gesis2dFEM)

