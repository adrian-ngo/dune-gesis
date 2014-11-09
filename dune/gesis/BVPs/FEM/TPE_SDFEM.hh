#ifndef DUNE_GESIS_TP_SDFEM_HH
#define DUNE_GESIS_TP_SDFEM_HH

//#include "ovlpistlgmresbackend.hh"
#include "TransportOperatorSDFEM.hh"


namespace Dune{
  namespace Gesis{

template<typename GFS
         , typename DGF_DARCY
         , typename SOURCE_TYPE
         , typename TP
         , typename IDT
         , typename SDT>
class TPE_SDFEM
{
  // extract types from GFS:
  typedef typename GFS::Traits::GridViewType GV;
  enum {dim=GV::dimension};
  typedef Dune::FieldVector<CTYPE,dim> DomainType;

  typedef typename GFS::Traits::FiniteElementMapType FEM;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type VCType;

  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<TP> BCType;
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<TP> BCExtension;

  typedef Dune::Gesis::StreamlineDiffusionOperator<TP,FEM,BCType,SOURCE_TYPE,IDT,SDT> LOP;
  typedef typename GFS::template ConstraintsContainer<REAL>::Type CC;

  // Retrieve VBE from GFS:
  typedef typename GFS::Traits::BackendType VBE;
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

#ifdef USE_NOVLP_AMG_PAR
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,REAL,REAL,REAL,CC,CC,true> GOS;
#else
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,REAL,REAL,REAL,CC,CC> GOS;
#endif

  typedef typename GOS::template MatrixContainer<REAL>::Type GLOBAL_MATRIX;



  // sequential solver
#ifdef USE_SUPERLU_SEQ
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
#endif

#ifdef USE_AMG_SEQ
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GOS> LS;
#endif

#ifdef USE_BICGSTAB_SEQ_TPE_CG
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
#endif

  // parallel solver
#ifdef USE_NOVLP_AMG_PAR
  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_AMG_SSOR<GOS> LS;
#endif

#ifdef USE_OVLP_AMG_PAR_TPE_CG
  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOS> LS;
  //typedef Dune::PDELab::ISTLBackend_GMRES_AMG_SSOR<GOS,VCType> LS;
#endif

#ifdef USE_OVLP_BCGS_SSORk_TPE_CG
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
#endif

#ifdef USE_OVLP_BCGS_ILU0_TPE_CG
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0<GFS,CC> LS;
#endif

#ifdef USE_OVLP_BCGS_SuperLU_TPE_CG
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS,CC> LS;
#endif


#ifdef USE_OVLP_GMRES_ILU0_TPE_CG
  typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<GFS,CC> LS;
#endif



private:
  const GFS& gfs;
  const GV& gv;
  VCType xSolution;

  const IDT& inputdata;
  const SDT& setupdata;

  TP& tp;
  BCType bctype;
  BCExtension dirichlet_function;
  SOURCE_TYPE& sourceterm;
  const Passenger::Type passenger;
  const EQ::Mode equationMode;
  bool bReUseMatrixHierarchy;
  LOP lop;
  MBE mbe;
  CC cc;  // Take care: CC cc must be declared before GOS gos! Otherwise you will get unprecdictable crashes! (cc: container for transformation)
  GOS gos;
  GLOBAL_MATRIX stiffness_matrix;
  UINT ls_verbosity;
  UINT istl_max_iter;
  LS linearsolver;
  Dune::PDELab::LinearSolverResult<REAL> ls_result;
  
public:
  
  // constructor:
  TPE_SDFEM( const GFS& gfs_
             , const IDT& inputdata_
             , const SDT& setupdata_
             , const DGF_DARCY& dgf_
             , TP& tp_
             , SOURCE_TYPE& sourceterm_
             , const Passenger::Type passenger_=Passenger::solute
             , const EQ::Mode equationMode_=EQ::forward
             , const bool& bReUseMatrixHierarchy_=false )
    : 
    gfs( gfs_ )
    , gv( gfs.gridView() )
    , xSolution( gfs, 0.0 )
    , inputdata( inputdata_ )
    , setupdata( setupdata_ )
    , tp( tp_ )
    , bctype( gv, tp )
    , dirichlet_function( gv, tp )  // inflow data
    , sourceterm( sourceterm_ )
    , passenger( passenger_ )
    , equationMode( equationMode_ )
    , bReUseMatrixHierarchy( bReUseMatrixHierarchy_ )
    , lop( inputdata,
           setupdata,
           tp,
           bctype,
           sourceterm,
           passenger,
           equationMode
           ) // Q1, P1
    , mbe(9)
    , gos( gfs, cc, gfs, cc, lop, mbe )
    , stiffness_matrix( gos )
    , ls_verbosity( std::max(0, inputdata.verbosity-VERBOSITY_EQ_DETAILS) )
    , istl_max_iter( std::min(5000, inputdata.transport_parameters.istl_max_iter) )
    
#ifdef USE_SUPERLU_SEQ
    , linearsolver( ls_verbosity )
#endif

#ifdef USE_AMG_SEQ
    , linearsolver( istl_max_iter, ls_verbosity, bReUseMatrixHierarchy )
#endif

#ifdef USE_BICGSTAB_SEQ_TPE_CG
    , linearsolver( istl_max_iter, ls_verbosity )
#endif

#ifdef USE_OVLP_AMG_PAR_TPE_CG
    , linearsolver( gfs, istl_max_iter, ls_verbosity, bReUseMatrixHierarchy )
#endif

#ifdef USE_NOVLP_AMG_PAR
    , linearsolver( gfs, istl_max_iter, ls_verbosity, bReUseMatrixHierarchy )
#endif

#ifdef USE_OVLP_BCGS_SSORk_TPE_CG
    , linearsolver( gfs, cc, istl_max_iter, 2, ls_verbosity )
#endif

#ifdef USE_OVLP_BCGS_ILU0_TPE_CG  // ILU0 is more stable than SSORk
    , linearsolver( gfs, cc, istl_max_iter, ls_verbosity )
#endif

#ifdef USE_OVLP_BCGS_SuperLU_TPE_CG
    , linearsolver( gfs, cc, istl_max_iter, ls_verbosity )
#endif

#ifdef USE_OVLP_GMRES_ILU0_TPE_CG
    , linearsolver( gfs, cc, istl_max_iter, ls_verbosity, 100, true )
#endif

  {

    int processor_rank = gos.trialGridFunctionSpace().gridView().comm().rank();

    std::string eqMode = (equationMode==EQ::forward) ? "forward" : "adjoint";
    std::string eqType = (passenger==Passenger::solute) ? "solute" : "heat";
    std::string logMsg = "=== SDFEM: New TPE: " + eqMode + " " + eqType;
    logger << std::endl << logMsg << std::endl;
    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
      std::cout << logMsg << std::endl;

#ifdef USE_OVLP_AMG_PAR_TPE_CG
#if HAVE_PARMETIS
    if(inputdata.transport_parameters.istl_atOnceAccu=="on"){
      Dune::Amg::Parameters amg_params;
      amg_params.setAccumulate( Dune::Amg::atOnceAccu );
      amg_params.setDebugLevel( ls_verbosity );
      linearsolver.setparams(amg_params);
      if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
        std::cout << "=== setAccumulate( Dune::Amg::atOnceAccu )" << std::endl;
    }
#else
    Dune::Amg::Parameters amg_params;
    amg_params.setAccumulate( Dune::Amg::atOnceAccu );
    amg_params.setDebugLevel( ls_verbosity );
    linearsolver.setparams(amg_params);
    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
      std::cout << "=== setAccumulate( Dune::Amg::atOnceAccu )" << std::endl;
#endif
#endif

    Dune::PDELab::constraints( bctype, gfs, cc); // fill container

    UINT dofs = gfs.size();
    UINT cdofs = cc.size();
    if( inputdata.verbosity>=VERBOSITY_EQ_DETAILS ){
      logger << "=== TPE_SDFEM:  dofs = " << dofs << std::endl;
      logger << "=== TPE_SDFEM: cdofs = " << cdofs << std::endl;
    }

    dofs = gv.comm().sum( dofs );
    cdofs = gv.comm().sum( cdofs );
    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS ){
      std::cout << "=== TPE_SDFEM:  dofs (all procs) = " << dofs << std::endl;
      std::cout << "=== TPE_SDFEM: cdofs (all procs) = " << cdofs << std::endl;
    }

    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
      std::cout << "=== TPE_SDFEM: interpolate dirichlet_function" << std::endl;
    Dune::PDELab::interpolate( dirichlet_function, gfs, xSolution );

    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
      std::cout << "=== TPE_SDFEM: set_nonconstrained_dofs( 0.0 ) " << std::endl;
    Dune::PDELab::set_nonconstrained_dofs( cc, 0.0, xSolution ); // clear interior

    stiffness_matrix = 0.0;
    gos.jacobian( xSolution, stiffness_matrix );

#ifdef USE_OVLP_AMG_PAR_TPE_CG
    std::string recycle = bReUseMatrixHierarchy ? "Recycle AMG coarsening hierarchy" : "No Recycling of AMG coarsening hierarchy, maybe not needed here.";
    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
      std::cout << "=== " << recycle << std::endl;
#endif
  }
  
  /*

    bool set_SourceNature( const SourceNature sn )
    {
    lop.set_SourceNature( sn );
    }

  */

  void reset_BoundaryConditions()
  {
    xSolution = 0.0;
  }


  //template<typename VCType>
  void set_rhs( const VCType& vc_solution )
  {
    lop.set_rhs( vc_solution );
  }
  
  template<typename VType,typename UType>
  void set_rhs( const VType& vc_1, const UType& vc_2 )
  {
    lop.set_rhs( vc_1,vc_2 );
  }


  void set_PointSourceLocation( const DomainType& pointsource_global )
  {
    lop.set_PointSourceLocation( pointsource_global );
  }


  void solve_adjoint( const DomainType& pointsource_global, VCType& xSolutionOut )
  {
    assert( equationMode == EQ::adjoint );
    reset_BoundaryConditions();
    set_PointSourceLocation( pointsource_global );
    solve_forward( xSolutionOut );
  }

  template<typename UType>
  void solve_adjoint( const UType& vcM1_adj
                     , VCType& xSolutionOut
                     )
  {
    assert( equationMode == EQ::adjoint );
    reset_BoundaryConditions();
    set_rhs( vcM1_adj );
    solve_forward( xSolutionOut );
  }

  template<typename VType,typename UType>
  void solve_adjoint( const VType& vcphi_adj   // eqn: (150) for k=0
                      , const UType& vc_phi0
                      , VCType& xSolutionOut
                      )
  {
    assert( equationMode == EQ::adjoint );
    reset_BoundaryConditions();
    set_rhs( vcphi_adj, vc_phi0 );
    solve_forward( xSolutionOut );
  }


  bool solve_forward( 
                     VCType& xSolutionOut 
                     , bool bReUseSolution=false // if set to true, re-use xSolutionOut
                      )
  {
    int processor_rank=gos.trialGridFunctionSpace().gridView().comm().rank();

    std::string eqMode = (equationMode==EQ::forward) ? "forward" : "adjoint";
    std::string eqType = (passenger==Passenger::solute) ? "solute" : "heat";
    std::string logMsg = "=== SDFEM: Solving TPE: " + eqMode + " " + eqType;
    logger << std::endl << logMsg << std::endl;
    if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
      std::cout << logMsg << std::endl;

    if( bReUseSolution )
      xSolution = xSolutionOut;
	
    Dune::Timer watch;

    // assemble residual:
    typedef typename GOS::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
    typedef typename Dune::PDELab::BackendVectorSelector<TrialGridFunctionSpace,REAL>::Type W;


    W r(gos.testGridFunctionSpace(), 0.0);

    gos.residual(xSolution, r); // residual is additive

    General::log_elapsed_time( watch.elapsed(),
                               gos.trialGridFunctionSpace().gridView().comm(),
                               inputdata.verbosity,
                               "TPE",
                               "Residual assembly" );
    
    // Solve Transport Equation here...
    REAL reduction = inputdata.transport_parameters.istl_reduction;

    VCType z(gos.trialGridFunctionSpace(), 0.0);
    watch.reset();


    // Example code for matrix I/O:
    /*
    typedef int GlobalId;
    typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;
    Communication mcomm( gv.comm() );
    Dune::storeMatrixMarket( stiffness_matrix.base(), std::string("matrix_TPE_SDFEM"), mcomm );
    Dune::storeMatrixMarket( r.base(), std::string("vector_TPE_SDFEM"), mcomm );
    */

    /*
    // DEBUG Code:
    // plot b-Ax:
    W b(gos.testGridFunctionSpace(), 0.0);
    b = r;
    Dune::Gesis::VTKPlot::output2vtu( gfs,
                                             b,
                                             "m0_residual_0",
                                             "m0_residual_0",
                                             1,
                                             true,
                                             0 );
    */

#ifdef USE_NOVLP_AMG_PAR
    // Take care: NOVLP_AMG Solver modifies the stiffness matrix!
    GLOBAL_MATRIX tmpStiffnessMatrix(gos);
    tmpStiffnessMatrix=0.0;
    tmpStiffnessMatrix.base() = stiffness_matrix.base();
    linearsolver.apply( tmpStiffnessMatrix, z, r, reduction ); 
#else
    linearsolver.apply( stiffness_matrix, z, r, reduction ); // solver makes right hand side consistent
#endif


    /*
    // DEBUG Code:
    // plot end defect r:
    Dune::Gesis::VTKPlot::output2vtu( gfs,
                                             r,
                                             "m0_residual_1",
                                             "m0_residual_1",
                                             1,
                                             true,
                                             0 );

    // plot b-Az:
    W y(gos.testGridFunctionSpace(), 0.0);
    stiffness_matrix.base().mv(z,y);
    b -= y;
    Dune::Gesis::VTKPlot::output2vtu( gfs,
                                             b,
                                             "m0_residual_2",
                                             "m0_residual_2",
                                             1,
                                             true,
                                             0 );
    */

    ls_result = linearsolver.result();


    /*
    //Dune::SeqGS<GLOBAL_MATRIX,VCType,VCType> preconditioner(stiffness_matrix,1,1);       // GS preconditioner
    //Dune::MatrixAdapter<GLOBAL_MATRIX,VCType,VCType> op(stiffness_matrix);   // operator  
    //Dune::LoopSolver<VCType> linearsolver(op,preconditioner,reduction,iter,1); // solver

    // Test GMRes with SSOR:
    Dune::RestartedGMResSolver<VCType> linearsolver(op,
                                                    preconditioner,
                                                    reduction,
                                                    100,
                                                    istl_max_iter,
                                                    ls_verbosity );
    Dune::InverseOperatorResult result;             // stores statistics about solver
    linearsolver.apply(z,r,result);                 // now solve
    ls_result.converged = result.converged;
    ls_result.iterations = result.iterations;
    ls_result.elapsed = result.elapsed;
    ls_result.reduction = result.reduction;
    ls_result.conv_rate = result.conv_rate;
    */

    xSolution -= z;
    xSolutionOut = xSolution;


#ifdef OUTPUT_TO_MATLAB

    logger << "Writing SDFEM stiffness matrix A to 'sdfem_A.dat' for Matlab" << std::endl;
    Dune::writeMatrixToMatlab( stiffness_matrix.base(), "sdfem_A.dat" );

    logger << "Writing SDFEM residual vector r to 'sdfem_r.dat' for Matlab" << std::endl;
    std::ofstream ofs_r( "sdfem_r.dat" );
    ofs_r << r;

    logger << "Writing SDFEM update vector z to 'sdfem_z.dat' for Matlab" << std::endl;
    std::ofstream ofs_z( "sdfem_z.dat" );
    ofs_z << z;

    // In matlab, you can check, if A*z = r is fulfilled.

    logger << "Writing SDFEM solution vector x to 'sdfem_x.dat' for Matlab" << std::endl;
    std::ofstream ofs_xSolution( "sdfem_x.dat" );
    ofs_xSolution << xSolution;

#endif


    General::log_elapsed_time( watch.elapsed(),
                               gos.trialGridFunctionSpace().gridView().comm(),
                               inputdata.verbosity,
                               "TPE",
                               "Linear Solver" );


    return true;
	
  };


  const Dune::PDELab::LinearSolverResult<REAL>& result() const{
    return ls_result;
  }


  std::string show_ls_result(){
    std::stringstream ss;
    std::string msg = ls_result.converged ? "=== TPE_SDFEM: Converged with" : "WARNING: TPE_SDFEM: Not converged! "; 

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
  
  void cutOverAndUndershoots( VCType& vc_m ){
    enum{blocksize = VBE::BlockSize};
    REAL injection_time = setupdata.transport_equation.injection_time;
    REAL concentration = setupdata.transport_equation.boundaries[0].stripes[0].value;
    REAL maxValue = injection_time * concentration;
    REAL minValue=0.0;
    for( int i=0; i<vc_m.flatsize(); i++ ){
      if( vc_m[i/blocksize][i%blocksize]<minValue )
        vc_m[i/blocksize][i%blocksize]=minValue;
      if( vc_m[i/blocksize][i%blocksize]>maxValue )
        vc_m[i/blocksize][i%blocksize]=maxValue;
    }
  };


  bool isAcceptableUpDown( const VCType& vc_m ){
    
    REAL minimum = 1E+100;
    REAL maximum = -1E+100;
    Dune::Gesis::General::logMinAndMax( vc_m, gv.comm(), minimum, maximum );

    REAL injection_time = setupdata.transport_equation.injection_time;
    REAL concentration = setupdata.transport_equation.boundaries[0].stripes[0].value;
    REAL maxValue = injection_time * concentration;
    // REAL minValue=0.0;

    if( std::abs(minimum) < 0.05*maxValue && maximum < 1.05*maxValue ){
      if(gv.comm().rank()==0)
        std::cout << "Overshoots and undershoots within 5 percent tolerance!" << std::endl;
      return true;
    }
    else
      return false;
  }


};

  }
}
#endif// DUNE_GESIS_TP_SDFEM_HH



