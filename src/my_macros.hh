/***
    This header file contains a lot of macros that are used to
    control the compilation flow in order to switch between
    different implementation options.

    Warning: Do not change these macros unless you are really sure of what you are doing!
    
    A. Ngo (2010-2014)
***/

#define MAJOR_VERSION 1
#define MINOR_VERSION 0
#define PDELAB_VERSION 23
#define BUILD_VERSION 101
// Some examples for switching compiler flags.
//
// release flags:
// CXXFLAGS='-O3 -DNDEBUG'
//
// debug flags:
// CXXFLAGS='-g3 -O0 -Wall -Wno-vla -pedantic -UNDEBUG'
// CXXFLAGS='-g3 -O0 -UNDEBUG'
//
// ************************************************
// The macro pMAX is very important. 
// It is used to set the polynomial degree for 
// the Q(k) or P(k) basis.
// ************************************************
#ifdef USE_FEM
#define pMAX 1  // For FEM, do not change this value.
#else
#define pMAX 1  // For DG, you can set pMax 2 or 3 if you want to try higher order DG.
#endif



// *********************************************************
//
// WARNING: The following macros are best left as they are!
//
// *********************************************************


/**
   These macros are compared to the flag "verbosity" inside the input XML 
   specified by the user.
   The condition (inputdata.verbosity >= VERBOSITY_EQ_SUMMARY)
   is valid only for all values of "verbosity" being at least 4.
   So, if "verbosity"<4 that part will not occur in the output.
   If want to insert a new verbosity level in between, 
   increment all of the following levels by +1.
*/
#define VERBOSITY_SILENT          0
#define VERBOSITY_INVERSION       1
#define VERBOSITY_TIMER_SUMMARY   2
#define VERBOSITY_TIMER_DETAILS   3
#define VERBOSITY_EQ_SUMMARY      4
#define VERBOSITY_EQ_DETAILS      5
#define VERBOSITY_LS_SUMMARY      6
#define VERBOSITY_LS_DETAILS      7
#define VERBOSITY_L_CONTRIB       8
#define VERBOSITY_DEBUG_LEVEL     9


// **************************************************
// The macro PARALLEL should be defined (=activated)
// by default. Different solvers and constraints 
// are being used for the sequential mode!
// **************************************************
#define PARALLEL 


// **************************************************
// The following macros are used to choose among
// a diversity of parallel linear solvers from 
// dune-istl. Do not change them if you are not familiar
// with the PDELab solver backends.
// **************************************************
#ifdef PARALLEL
// Note:
// The solver for the GWE and the TPE can be different, 
// They should be both overlapping!

// Choose linear solvers for FEM and CCFV:
#define USE_OVLP_AMG_PAR_GWE_CG   // for GWE(FEM) : uses AMG
#define USE_OVLP_AMG_PAR_GWE_CCFV // for GWE(CCFV) : uses AMG

// Choose TPE linear solver for SDFEM:
#define USE_OVLP_GMRES_ILU0_TPE_CG
//#define USE_OVLP_AMG_PAR_TPE_CG    // for TPE(SDFEM): uses AMG

// Choose TPE linear solver for DG:
#define USE_OVLP_GMRES_ILU0_TPE_DG

// Choose linear solver for L2-Projection:
#define USE_OVLP_AMG_L2PROJ  // L2 (DG -> CG) : uses AMG


#else // sequential case (This is of minor interest here!):

#define USE_SUPERLU_SEQ       // Sequential SuperLU
#define USE_SUPERLU_SEQ_TPE_DG

#endif




// ************************************************
// These macros should be active only for debugging 
// purposes:
// ************************************************
//#define DEBUG_LOG         // activate more debug logging inside logger output
//#define OUTPUT_TO_MATLAB  // Use dune-istl function to store the sparse stiffness matrix in matlab format.



// *******************************************
// BE CAREFUL!
// These macros should not be touched at all.
// *******************************************
/* 
   PRESERVE_STRUCTURE = true :
   The backend vector is stored in a structured 2D or 3D way to hdf5 
   such that on can directly plot the solution from the hdf5 file
   If too many processes fight for write access on the overlap elements
   this hdf5 output takes too long!
   This would be necessary only when used with parallel-in-parallel mode.

   PRESERVE_STRUCTURE = false:
   The backend vector is stored in as a 1D vector to hdf5.
   This should be much faster because there each process has its own
   block on the hdf5 file.
*/
#define PRESERVE_STRUCTURE true
#define PARALLEL_KSI_JQ  // This macro enables parallel computation of ksi times JQ.
#define USE_ALL_ADJOINTS // This macro should be the default unless you want to save time during the m1 inversion and can work with a result of lower quality.


#define USE_SECURE_MODE // Keep this macro activated. Otherwise, non-ovlp AMG linear solver changes the stiffness matrix.
#define CLEAN_S   // This removes all HDF5 files for the Sensitivity fields after usage and saves your diskspace.
#define CLEAN_JQ  // This removes all HDF5 files for the Sensitivity fields after usage and saves your diskspace.


//#define HOMOGEN  // CENTER_SQUARE // <-- testcase for forward solver
//#define HARMONIC_AVG // makes no difference at all!

//===================================================================
// The compiler flag -DUSE_FEM is used to 
// switch on the Standard Galerkin version and to swith off the DG version.
// By default, the DG version is built.
//===================================================================

#ifdef USE_FEM
#define GroundWaterEquation GWE_Galerkin
#define GEP_Equation GWE_Galerkin
#define TransportEquation TPE_SDFEM
#define GradientVectorField Dune::Gesis::DiscreteGridFunctionDarcy
#else
#define L2ProjectionOfM0
#define GroundWaterEquation GWE_CCFV
#define GEP_Equation GWE_CCFV
#define TransportEquation TPE_DG
#define GradientVectorField Dune::Gesis::RT0FluxDGF
//#define PLOT_GHOST // This must not be used for USE_FEM
#endif



#ifndef USE_FEM
#define NEW_GV   // use the reordered gridview, not implemented for HOMOGEN case
#ifdef USE_YASP
#define USE_DGF_PressureField
#endif
#endif


#define USE_CACHE  // uses local basis cache to speed up matrix assembly do not change this!

// epsilon value used in the code (mostly for comparison of REAL)
#define GEO_EPSILON 1e-6


//delete old output data (LOG and *.VTU)
#define CLEAN_OUTPUT 1 // 1 = clean OUTPUT directory at the beginning
#define CLEAN_BUFFER 0 // 1 = clean BUFFER directory at the beginning


//Define # of threads used for the FFT calculations (comment it out if you don't have the fftw version for distributed memory)
// IMPORTANT: link with '-lfftw3_threads' and use FFTW 3.3 (alpha) or higher
//#define FFT_THREADS 2


// The ESTIMATION_VARIANCE can be de-activated if you are not interested in uncertainty quantification.
#define ESTIMATION_VARIANCE


// Define the data type of the HDF5 files here. 32 or 64 bit!!!!
#define HDF5_DATA_TYPE H5T_IEEE_F64LE  //define for 64 bit machine
//#define HDF5_DATA_TYPE H5T_IEEE_F32LE //define for 32 bit machine


#define CONVECTIVE_FORMULATION // Leave this flag be defined for the SDFEM method to work properly!
#define REGULARIZED_DIRICHLET  // Leave this flag be defined for the SDFEM method to have the Dirichlet boundary function to be continuous.
//

// (This is of minor interest here!):
#if defined USE_NOVLP_AMG_PAR || defined USE_NOVLP_AMG_PAR_TR
#define USE_NOVLP_MODE
#endif
