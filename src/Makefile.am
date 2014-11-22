#
# Some variables pointing to the locations of external libs:
#
# 
# The special variable BOOST_INCL is required only due to the faulty boost check 
# in dune-common (2.3.1):
BOOST_INCL = -I$(BOOST_CPPFLAGS)/include
# 
#
##############################################################
#
# applications:
#
##############################################################
SUBDIRS =

noinst_PROGRAMS = \
	gesis2dDG gesis2dFEM \
	gesis3dDG  gesis3dFEM \
	gesisUG2dDG \
	gesisALU3dDG

# 2D DG version for YASP:
gesis2dDG_SOURCES = geo_inversion_2d.cc

gesis2dDG_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_YASP \
	$(DUNEMPICPPFLAGS) \
	$(HDF5_CPPFLAGS) -DH5_USE_16_API \
	$(FFTW3_CPPFLAGS) \
	$(BOOST_INCL) \
	$(PARMETIS_CPPFLAGS) \
	$(SUPERLU_CPPFLAGS)

gesis2dDG_LDADD = \
	$(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS)	\
	$(HDF5_LIBS) \
	$(FFTW3_LIBS) \
	$(PARMETIS_LIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(LDADD)

gesis2dDG_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNE_LDFLAGS)

# 2D DG version for UG:
gesisUG2dDG_SOURCES = ug_2d.cc

gesisUG2dDG_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_UG \
	$(DUNEMPICPPFLAGS) \
	$(HDF5_CPPFLAGS) -DH5_USE_16_API \
	$(FFTW3_CPPFLAGS) \
	$(BOOST_INCL) \
	$(PARMETIS_CPPFLAGS) \
	$(SUPERLU_CPPFLAGS) \
	$(UG_CPPFLAGS)

gesisUG2dDG_LDADD = \
	$(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS)	\
	$(HDF5_LIBS) \
	$(FFTW3_LIBS) \
	$(PARMETIS_LIBS) \
	$(SUPERLU_LIBS) \
	$(UG_LIBS) \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(LDADD)

gesisUG2dDG_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNE_LDFLAGS)


# 3D DG version for YASP:
gesis3dDG_SOURCES = geo_inversion_3d.cc

gesis3dDG_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_YASP \
	$(HDF5_CPPFLAGS) -DH5_USE_16_API \
	$(FFTW3_CPPFLAGS) \
	$(BOOST_INCL) \
	$(PARMETIS_CPPFLAGS) \
	$(SUPERLU_CPPFLAGS) \
	$(DUNEMPICPPFLAGS)

gesis3dDG_LDADD = \
	$(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS)	\
	$(HDF5_LIBS) \
	$(FFTW3_LIBS) \
	$(PARMETIS_LIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(LDADD)

gesis3dDG_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNE_LDFLAGS)

# 3D DG version for ALUGrid:
gesisALU3dDG_SOURCES = alugrid_3d.cc

gesisALU3dDG_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_ALUGRID \
	$(DUNEMPICPPFLAGS) \
	$(HDF5_CPPFLAGS) -DH5_USE_16_API \
	$(FFTW3_CPPFLAGS) \
	$(BOOST_INCL) \
	$(PARMETIS_CPPFLAGS) \
	$(SUPERLU_CPPFLAGS)

gesisALU3dDG_LDADD = \
	$(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS)	\
	$(HDF5_LIBS) \
	$(FFTW3_LIBS) \
	$(PARMETIS_LIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(LDADD)

gesisALU3dDG_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNE_LDFLAGS)


# 2D FEM version for YASP:
gesis2dFEM_SOURCES = geo_inversion_2d.cc

gesis2dFEM_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_YASP -DUSE_FEM \
	$(HDF5_CPPFLAGS) -DH5_USE_16_API \
	$(FFTW3_CPPFLAGS) \
	$(BOOST_INCL) \
	$(PARMETIS_CPPFLAGS) \
	$(SUPERLU_CPPFLAGS) \
	$(DUNEMPICPPFLAGS)

gesis2dFEM_LDADD = \
	$(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS)	\
	$(HDF5_LIBS) \
	$(FFTW3_LIBS) \
	$(PARMETIS_LIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(LDADD)

gesis2dFEM_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNE_LDFLAGS)

# 3D FEM version for YASP:
gesis3dFEM_SOURCES = geo_inversion_3d.cc

gesis3dFEM_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_YASP -DUSE_FEM \
	$(HDF5_CPPFLAGS) -DH5_USE_16_API \
	$(FFTW3_CPPFLAGS) \
	$(BOOST_INCL) \
	$(PARMETIS_CPPFLAGS) \
	$(SUPERLU_CPPFLAGS) \
	$(DUNEMPICPPFLAGS)

gesis3dFEM_LDADD = \
	$(DUNEMPILDFLAGS) \
	$(DUNEMPILIBS)	\
	$(HDF5_LIBS) \
	$(FFTW3_LIBS) \
	$(PARMETIS_LIBS) \
	$(SUPERLU_LIBS) \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(LDADD)

gesis3dFEM_LDFLAGS = $(AM_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNE_LDFLAGS)

# default app:

dune_gesis_SOURCES = dune_gesis.cc

dune_gesis_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS)
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
dune_gesis_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD)
dune_gesis_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)
##############################################################




##############################################################

# don't follow the full GNU-standard
# we need automake 1.9
AUTOMAKE_OPTIONS = foreign 1.9

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-geometry=$(DUNE_GEOMETRY_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-typetree=$(DUNE_TYPETREE_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

EXTRA_DIST = CMakeLists.txt

include $(top_srcdir)/am/global-rules

all: all-recursive
clear:
	-rm -v *.vtu *.log *.h5 *.xmf *.pvd
