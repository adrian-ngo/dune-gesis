#
# Some variables pointing to the locations of external libs:
#
# Paths to HFD5
HDF5_LIBS = -L$(HOME)/Software/hdf5-1.8.13/lib -lhdf5
HDF5_INCL = -I$(HOME)/Software/hdf5-1.8.13/include -DH5_USE_16_API
#
# Paths to FFTW3
FFTW3_LIBS = -L$(HOME)/Software/fftw-3.3.4/lib -lfftw3_mpi -lfftw3
FFTW3_INCL = -I$(HOME)/Software/fftw-3.3.4/include
# 
#
BOOST_INCL = -I$(HOME)/Software/boost_1_47_0
#
#
##############################################################
#
# applications:
#
##############################################################
SUBDIRS =

noinst_PROGRAMS = dune_gesis gesis2dDG

# 2D DG version for YASP:
gesis2dDG_SOURCES = geo_inversion_2d.cc

gesis2dDG_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_YASP \
	$(HDF5_INCL) \
	$(FFTW3_INCL) \
	$(BOOST_INCL) \
	$(PARMETIS_CPPFLAGS) \
	$(SUPERLU_CPPFLAGS) \
	$(DUNEMPICPPFLAGS)

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

