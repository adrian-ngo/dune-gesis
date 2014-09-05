// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_L2PROJECTION_HH
#define DUNE_PDELAB_L2PROJECTION_HH

#include<dune/pdelab/gridoperator/gridoperator.hh>

//#include "general.hh"


namespace Dune{
  namespace GeoInversion{

template<
  typename GFS
  , typename DGF
  , typename IDT
  >
class L2Projection
{
  // extract types from GFS:
  typedef typename GFS::Traits::GridViewType GV;
  typedef typename GFS::Traits::BackendType VBE;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type VCType;

  typedef Dune::GeoInversion::L2ProjectionOperator<DGF> LOP;
  typedef typename GFS::template ConstraintsContainer<REAL>::Type CC;
  
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,REAL,REAL,REAL,CC,CC> GOS;
  
  typedef typename GOS::template MatrixContainer<REAL>::Type GLOBAL_MATRIX;


  // set up solver
#ifdef USE_OVLP_AMG_L2PROJ
  typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GOS> LS;
  //typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOS> LS;
#endif

#ifdef USE_OVLP_BCG_L2PROJ
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS, CC> LS;
  //typedef Dune::PDELab::ISTLBackend_OVLP_CG_SSORk<GFS, CC> LS;
#endif




private:
  GV gv; // Using a copy of the gridview is cheap and safe.
  const GFS& gfs;
  const IDT& inputdata;


  VCType xSolution;

  CC cc; // Take care: CC cc must be declared before GOS gos! Otherwise you will get unprecdictable crashes! (cc: container for transformation)
  LOP lop;
  MBE mbe;
  GOS gos;
  GLOBAL_MATRIX stiffness_matrix;
  UINT ls_verbosity;
  bool bReUseMatrixHierarchy;
  LS linearsolver;
  Dune::PDELab::LinearSolverResult<REAL> ls_result;

  UINT iCounter; // used to count repeated solving of the same test problem

public:
  
  // constructor:
  L2Projection( const GV& gv_
                , const GFS& gfs_
                , const DGF& dgf_
                , const REAL l2_diffusion
                , const IDT& inputdata_ ) 
    : 
    gv( gv_ )
    , gfs( gfs_ )
    , inputdata( inputdata_ )
    , xSolution( gfs, 0.0 )
    , lop( dgf_, l2_diffusion )
    , mbe(9)
    , gos( gfs, cc, gfs, cc, lop, mbe)
    , stiffness_matrix( gos )
    , ls_verbosity( std::max(0, inputdata.verbosity-VERBOSITY_EQ_DETAILS) )
    , bReUseMatrixHierarchy(false)
    

#ifdef USE_OVLP_AMG_L2PROJ
    , linearsolver( gfs, 5000, ls_verbosity, bReUseMatrixHierarchy ) // AMG
#endif

#ifdef USE_OVLP_BCG_L2PROJ
    , linearsolver(gfs, cc, 5000, 5, ls_verbosity)
#endif

  {

    int processor_rank=gos.trialGridFunctionSpace().gridView().comm().rank();
    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
      std::cout << "=== L2-Projection: diffusion_factor = " 
                << l2_diffusion << std::endl;
    
#ifdef USE_OVLP_AMG_L2PROJ
    Dune::Amg::Parameters amg_params;
    amg_params.setAccumulate( Dune::Amg::atOnceAccu );
    amg_params.setDebugLevel( ls_verbosity );
    //linearsolver.setparams(amg_params);
    linearsolver.setParameters(amg_params);
    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
      std::cout << "=== L2-Projection: AMG with atOnceAccu" << std::endl;
#endif

    iCounter = 0;

    //Dune::PDELab::constraints( bctype, gfs, cc); // fill container

    UINT dofs = gfs.globalSize();
    UINT cdofs = cc.size();
    if( inputdata.verbosity>=VERBOSITY_EQ_DETAILS ) {
      logger << "=== L2-Projection:  dofs = " << dofs << std::endl;
      logger << "=== L2-Projection: cdofs = " << cdofs << std::endl;
    }
    
    dofs = gv.comm().sum( dofs );
    cdofs = gv.comm().sum( cdofs );
    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS ) {
      std::cout << "=== L2-Projection:  dofs (all procs) = " << dofs << std::endl;
      std::cout << "=== L2-Projection: cdofs (all procs) = " << cdofs << std::endl;
    }
    
    //Dune::PDELab::set_nonconstrained_dofs( cc, 0.0, xSolution ); // clear interior

    stiffness_matrix = 0.0;
    gos.jacobian( xSolution, stiffness_matrix );
    
  }


  bool solve_forward( VCType& xSolutionOut ) {

    Dune::Timer watch;

    int processor_rank=gos.trialGridFunctionSpace().gridView().comm().rank();

    // assemble residual:
    typedef typename GOS::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
    typedef typename Dune::PDELab::BackendVectorSelector<TrialGridFunctionSpace,REAL>::Type W;
    
    W r(gos.testGridFunctionSpace(), 0.0);

    gos.residual(xSolution, r); // residual is additive

    General::log_elapsed_time( watch.elapsed(),
                               gos.trialGridFunctionSpace().gridView().comm(),
                               inputdata.verbosity,
                               "L2Proj",
                               "Residual assembly" );
    watch.reset();

    // Solve Flow Equation here...
    REAL reduction = 1e-10;
    
    VCType z(gos.trialGridFunctionSpace(), 0.0);
    watch.reset();



#ifdef DEBUG_CODE
    std::vector<REAL> vector_r;
    r.std_copy_to ( vector_r );

    std::string filename;    
    if( bAdjoint )
      filename = "r_Adj_p" + number2string( processor_rank ) + "_i" + number2string( iCounter ) + ".txt";
    else
      filename = "r_Fwd_p" + number2string( processor_rank ) + "_i" + number2string( iCounter ) + ".txt";

    std::ofstream ofs( filename, std::ios::trunc );
    if( ofs.is_open() ){
      for( UINT i=0; i<vector_r.size(); i++ ){
        ofs << i << ": " << vector_r[i] << std::endl;
      }
    }
#endif

    

#ifdef DEBUG_CODE
    std::string A_filename;
    if( bAdjoint )
      A_filename = "A_Adj_p" + number2string( processor_rank ) + "_i" + number2string( iCounter ) + ".txt";
    else
      A_filename = "A_Fwd_p" + number2string( processor_rank ) + "_i" + number2string( iCounter ) + ".txt";
	Dune::writeMatrixToMatlab( stiffness_matrix.base(), A_filename );
#endif


#ifdef USE_NOVLP_AMG_PAR
    // Take care: NOVLP_AMG Solver modifies the stiffness matrix!
    GLOBAL_MATRIX tmpStiffnessMatrix(gos);
    tmpStiffnessMatrix=0.0;
    tmpStiffnessMatrix.base() = stiffness_matrix.base();
    linearsolver.apply( tmpStiffnessMatrix, z, r, reduction ); 
#else
    linearsolver.apply( stiffness_matrix, z, r, reduction );
#endif

    ls_result = linearsolver.result();

#ifdef DEBUG_CODE
    std::string A_filename_postSolve;
    if( bAdjoint )
      A_filename_postSolve 
        = "A_postSolve_Adj_p" + number2string( processor_rank ) + "_i" + number2string( iCounter ) + ".txt";
    else
      A_filename_postSolve 
        = "A_postSolve_Fwd_p" + number2string( processor_rank ) + "_i" + number2string( iCounter ) + ".txt";

    Dune::writeMatrixToMatlab( stiffness_matrix.base(), A_filename_postSolve );
#endif


    xSolution -= z;
    xSolutionOut = xSolution;

    General::log_elapsed_time( watch.elapsed(),
                               gos.trialGridFunctionSpace().gridView().comm(),
                               inputdata.verbosity,
                               "L2Proj",
                               "Linear Solver" );
    return true;
	
  };




  const Dune::PDELab::LinearSolverResult<REAL>& result() const{
    return ls_result;
  }


  std::string show_ls_result(){
    std::stringstream ss;
    std::string msg = ls_result.converged ? "=== L2-Projection: Converged with" : "=== L2-Projection: WARNING: Not converged! "; 

    ss << msg 
       << " T=" << ls_result.elapsed
       << " IT=" << ls_result.iterations
       << " TIT=" << ls_result.elapsed/ls_result.iterations
       << " reduction=" << ls_result.reduction
       << " rate=" << ls_result.conv_rate;
    return ss.str();
  }




  void cutUndershoots( VCType& vc_m ){
    enum{blocksize = VBE::BlockSize};
    for( int i=0; i<vc_m.flatsize(); i++ ){
      if( vc_m[i/blocksize][i%blocksize]<1e-12 )
        vc_m[i/blocksize][i%blocksize]=0.0;
    }
  };

  void cutOvershoots( VCType& vc_m ){
    enum{blocksize = VBE::BlockSize};
    for( int i=0; i<vc_m.flatsize(); i++ ){
      if( vc_m[i/blocksize][i%blocksize]>-1e-12 )
        vc_m[i/blocksize][i%blocksize]=0.0;
    }
  };
  
  void cutOverAndUndershoots( VCType& vc_m,
                              const REAL maxValue, 
                              const REAL minValue
                              ){
    enum{blocksize = VBE::BlockSize};
    for( int i=0; i<vc_m.flatsize(); i++ ){
      
      if( vc_m[i/blocksize][i%blocksize]<minValue )
        vc_m[i/blocksize][i%blocksize]=minValue;

      if( vc_m[i/blocksize][i%blocksize]>maxValue )
        vc_m[i/blocksize][i%blocksize]=maxValue;

    }
  };


  bool isAcceptableUpDown( const VCType& vc_m, 
                           const REAL maxValue=0.0, 
                           const REAL minValue=0.0, 
                           const REAL tolerance=5.0 // in percent
                           ){
    
    REAL minimum = 1E+100;
    REAL maximum = -1E+100;
    Dune::GeoInversion::General::logMinAndMax( vc_m, gv.comm(), minimum, maximum );

    REAL lowerlimit=minimum;
    REAL upperlimit=maximum;
    
    if(maxValue>minValue){
      lowerlimit=-tolerance/100.0*(maxValue-minValue);
      upperlimit=(100.0+tolerance)/100.0*(maxValue-minValue);
    }
    else{
      lowerlimit=-tolerance/100.0*(maximum);
    }
    
    if( minimum >= lowerlimit
        && 
        maximum <= upperlimit
        ){

      if( inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
        std::cout << "Overshoots and undershoots within 5 percent tolerance!" << std::endl;
      return true;
    }
    else
      return false;
  }

  
};


  }
}

#endif// GROUND_WATER_EQUATION_HH



