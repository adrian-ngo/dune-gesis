#define MAJOR_VERSION 1
#define MINOR_VERSION 3
#define PDELAB_VERSION 21
#define BUILD_VERSION 20140730
// release flags:
// CXXFLAGS='-O3 -funroll-loops -DNDEBUG'
// CXXFLAGS='-O3 -funroll-loops -fno-strict-aliasing -DNDEBUG'
// debug flags:
// CXXFLAGS='-g3 -O0 -Wall -Wno-vla -pedantic -UNDEBUG'
// CXXFLAGS='-g -O0 -Wno-unused-result -fno-strict-aliasing -UNDEBUG'

// release build
// CXXFLAGS='-g -O3 -DNDEBUG -Wall -Wno-unused-but-set-variable -Wno-unused-variable -Wno-literal-suffix -Wno-unused-local-typedefs'

//#define FLOATING_POINT_EXEPTION_CHECK

//#define OLES_TEST

//#define BACKUP_REPLAY

#define PARALLEL_KSI_JQ

#define USE_ALL_ADJOINTS // Check: Is this really needed?

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

//#define PRESERVE_STRUCTURE false
#define PRESERVE_STRUCTURE true

// #define USE_SEQ_WRITE_TO_HDF5  // This may speed up hdf5 write of /Sensitivity if too many processors try to write a small hdf5 file.

//#define TEST_MAX
#define TWO_NORM_STRATEGY //  test new error fraction strategy

#define USE_SECURE_MODE

#define CLEAN_S // remove h5 file for Sensitivity field and JQ after usage
#define CLEAN_JQ

// Note: Use these two MACROS alternatively!

//#define HOMOGEN  // CENTER_SQUARE // <-- testcase for forward solver: homogeneous field

#define USE_Qk   // use Qk FEM vs. Pk FEM
#define USE_CUBE // use cube vs. simplex grid elements

//#define HARMONIC_AVG // make no difference at all!

//===================================================================
// The compiler flag -DUSE_FEM is used to 
// switch on the Standard Galerkin version and to swith off the DG version.
// By default, the DG version is built.
//===================================================================

#ifdef USE_FEM
#define pMAX 1
//#define L2ProjectionOfM0
#define GroundWaterEquation GWE_Galerkin
#define GEP_Equation GWE_Galerkin
#define TransportEquation TPE_SDFEM
#define GradientVectorField Dune::PDELab::DiscreteGridFunctionDarcy
#else
#define pMAX 1
#define L2ProjectionOfM0
#define GroundWaterEquation GWE_CCFV
#define GEP_Equation GWE_CCFV
#define TransportEquation TPE_DG
#define GradientVectorField Dune::GeoInversion::RT0FluxDGF
//#define PLOT_GHOST // This must not be used for USE_FEM
#endif


#define SPECIAL_YASP_PARTITION
//#define USE_YASP     (moved to "Makefile.am")
//#define USE_ALUGRID  (moved to "Makefile.am")
//#define USE_UG       (moved to "Makefile.am")


#ifndef USE_FEM
#define NEW_GV   // use the reordered gridview, not implemented for HOMOGEN case
#ifdef USE_YASP
#define USE_DGF_PressureField
#endif
#endif


//#define WELL_FRACTURE_MODEL
//#define USE_NEW_GRIDOPERATORSPACE
#define USE_CACHE  // uses local basis cache to speed up matrix assembly do not change this!


// epsilon value used in the code (mostly for comparison of REAL)
#define GEO_EPSILON 1e-6

// Run in parallel mode or not. It uses other solvers than sequential mode!
#define PARALLEL 

//delete old output data (LOG and *.VTU)
#define CLEAN_OUTPUT 1 // 1 = clean OUTPUT directory at the beginning
#define CLEAN_BUFFER 0 // 1 = clean BUFFER directory at the beginning


//Define # of threads used for the FFT calculations (comment it out if you don't have the fftw version for distributed memory)
// IMPORTANT: link with '-lfftw3_threads' and use FFTW 3.3 (alpha) or higher
//#define FFT_THREADS 2


#define USE_HDF5 // This requires HDF5 version 1.6.6 or 1.8.5 to be installed on your machine! If not available switch this off and MPI_Bcast will be used to distribute the Y-field. That takes much longer.



// Choose appropriate linear solvers from dune-istl:
#ifdef PARALLEL

// Note:
// The solver for the GWE and the TPE can be different, 
// BUT they should be both overlapping or both non-overlapping!!!!!

// choose GWE solvers:
#define USE_OVLP_AMG_PAR_GWE_CCFV // for GWE(CCFV) : uses AMG
#define USE_OVLP_AMG_PAR_GWE_CG   // for GWE(FEM) and for L2 (DG -> CG) : uses AMG
//#define USE_OVLP_PAR_GWE_CG
//#define USE_AMG_PAR         // Parallel Overlapping AMG
//#define USE_NOVLP_AMG_PAR   // Parallel Non-Overlapping AMG
//#define USE_OVLP_PAR        // Parallel Overlapping BICGSTAB

// choose TPE solvers for SDFEM:
//#define USE_OVLP_BCGS_SuperLU_TPE_CG
//#define USE_OVLP_AMG_PAR_TPE_CG    // for TPE(SDFEM): uses AMG
//#define USE_OVLP_BCGS_SSORk_TPE_CG
//#define USE_OVLP_BCGS_ILU0_TPE_CG
#define USE_OVLP_GMRES_ILU0_TPE_CG

// choose TPE solvers for DG:
//#define USE_OVLP_BCGS_ILU0_TPE_DG       // for TPE(DG): BiCGStab ILU0
//#define USE_OVLP_BCGS_SSORk_TPE_DG       // for TPE(DG): BiCGStab SSORk
#define USE_OVLP_GMRES_ILU0_TPE_DG

// choose solver for L2 Projection:
#define USE_OVLP_AMG_L2PROJ
//#define USE_OVLP_BCG_L2PROJ


#else // sequential case:

//#define USE_AMG_SEQ         // Sequential AMG
#define USE_SUPERLU_SEQ       // Sequential SuperLU
#define USE_SUPERLU_SEQ_TPE_DG
//#define USE_BICGSTAB_SEQ    // Sequential BICGSTAB

#endif



#if defined USE_NOVLP_AMG_PAR || defined USE_NOVLP_AMG_PAR_TR
#define USE_NOVLP_MODE
#endif



//#define DEBUG_PLOT

#ifdef DEBUG_PLOT
//#define VTK_PLOT_WELLS
#define VTK_PLOT_PSI_HEAD //vtk output of adjoint state of head: on, otherwise off
#define VTK_PLOT_PSI_transport //vtk output of adjoint state of transport: on, otherwise off
#define VTK_PLOT_S   //vtk output of sensitivities: on, otherwise off
#define VTK_PLOT_JQ //vtk output of JQ: on, otherwise off
#endif //DEBUG_PLOT
//
//#define VTK_PLOT_Y_SMOOTH
#define VTK_PLOT_P_FIELD //vtk output of hydr. head on, otherwise off
#define VTK_PLOT_Q_FIELD //vtk output of darcy flux on, otherwise off
#define VTK_PLOT_RGV
#define VTK_PLOT_C_FIELD //vtk output of concentration on, otherwise off

//#define VTK_PLOT_TRIAL_FIELDS
//#define VTK_PLOT_YTRY  // plot intermediate estimated YFields
//#define DEBUG_PLOT_H_OLD    // See if 'h' gets transferred correctly via HDF5 io.
//#define DEBUG_PLOT_M0_OLD   // See if 'm0' gets transferred correctly via HDF5 io.


// Note: Use these MACROS to switch on/off graphical output
// #define VTK_PLOT_ASCII   //vtk output in ascii mode, otherwise binary mode
#define VTK_PLOT_HEAT_FIELD //vtk output of heat on, otherwise off
#define VTK_PLOT_EL_POTENTIAL_FIELD //vtk output of geoelectrical potential on, otherwise off

// #define VTK_PLOT_KFIELD  // This would plot the exp(logK)-field.
//#define VTK_PLOT_PSI_heat_transport //vtk output of adjoint state of transport: on, otherwise off
//#define VTK_PLOT_PSI_GP // vtk output of adjoint state variables for geoelectrical potential: on, otherwise off


#define PLOT_VTU_CR  // plot the conditional fields as VTU

#define VTK_PLOT_Y_FIELD //vtk output of Yfield on, otherwise off
#define VTK_PLOT_Y_OLD //vtk output of Y_old in each iteration: on, otherwise off

//#define VTK_PLOT_ESTIMATES  //vtk output of head/transport/etc related to the estimated field!

#define ESTIMATION_VARIANCE // the estimation variance of the inversion will be calculated (might be a lot to do!)

//#define DEBUG_LOG
//#define DEBUG_LOG_LEVEL_2 // Use this ONLY for small numbers of elements! Otherwise gigantic logfiles will be produced!


//#define OUTPUT_TO_MATLAB  // Use dune-istl function to store the large sparse stiffness matrix in matlab format
#define VERBOSE_LEVEL 2 // Output level of the ISTL linear solvers: 0=least, 2=most


//#define _USE_FLOAT_

#define LINEAR_PROBLEM // As opposed to the non-linear case for shock-capturing. Then we would need the Newton solver!

// Set this flag to smear the delta peak of the adjoint! Sampling Volume procedure in m0 adjoint disabled!!!!!
//#define ADJ_SMEARING

#define REGULARIZED_DIRICHLET
#define CONVECTIVE_FORMULATION
#define JacVol // switch off numerical differentiation of alpha_vol, makes assembling faster
#define JacBnd // switch off numerical differentiation of alpha_boundary


#define EFFECTIVE_MESHSIZE // for the correct meshsize at flow field angles at 45 degrees

/* 
 * Macros for the SDFEM operator to distinguish between 
 * the stationary scalar convection diffusion problem 
 * and 
 * the stationary solute transport equation (default)
 */
//#define _SCALAR_CDE_


/* 
 * Macros for debugging the Yfield while reading parallel hdf5
 */
//#define _L_FIELD_
//#ifdef DIMENSION3
//#define _EXTRA_DATA_ // This is necessary only for the completeness of the xdmf-format, but not for anything else at the moment.
//#endif

//#define RAVIART_THOMAS
//#define FLOATING_POINT_EXEPTION_CHECK


#ifdef USE_HDF5
// define the data type of the HDF5 files here. 32 or 64 bit!!!!
#define HDF5_DATA_TYPE H5T_IEEE_F64LE  //define for 64 bit machine
//#define HDF5_DATA_TYPE H5T_IEEE_F32LE //define for 32 bit machine
#include "hdf5.h"
#endif






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
