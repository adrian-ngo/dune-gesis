#ifndef DUNE_GESIS_TPE_DG_HH
#define DUNE_GESIS_TPE_DG_HH

//#include "ovlpistlgmresbackend.hh"
#include "dune/gesis/BVPs/TransportParameters.hh"
#include "TransportOperatorDG.hh"

namespace Dune {

  namespace Gesis {

    template<
      typename GFS
      , typename DGF_DARCY
      , typename SOURCE_TYPE
      , typename TP
      , typename IDT
      , typename SDT
      >
    class TPE_DG    // TransportEquationDG
    {

      //typedef Dune::PDELab::TransportProblem< GV, REAL, DGF_DARCY, CTransportParameters > TP;

      // extract types from GFS:
      typedef typename GFS::Traits::GridViewType GV;
      enum {dim=GV::dimension};
      typedef Dune::FieldVector<CTYPE,dim> DomainType;

      typedef typename GFS::Traits::FiniteElementMapType FEM;
      typedef typename GFS::Traits::BackendType VBE;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type VCType;

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<TP> BCType;
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<TP> BCExtension;

      typedef Dune::Gesis::ConvectionDiffusionDG<TP,
                                                 FEM,
                                                 SOURCE_TYPE,
                                                 IDT,
                                                 SDT
                                                 > LOP;

      typedef typename GFS::template ConstraintsContainer<REAL>::Type CC;
      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

#ifdef USE_NOVLP_MODE
      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,REAL,REAL,REAL,CC,CC,true> GOS;
#else
      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,REAL,REAL,REAL,CC,CC> GOS;
#endif

      typedef typename GOS::template MatrixContainer<REAL>::Type GLOBAL_MATRIX;



      // sequential solvers
#ifdef USE_SUPERLU_SEQ_TPE_DG
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
#endif

#ifdef USE_AMG_SEQ_TPE_DG
      typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GOS> LS;
#endif

#ifdef USE_BICGSTAB_SEQ_TPE_DG
      typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
#endif

      // parallel solvers

#ifdef USE_OVLP_BCGS_ILU0_TPE_DG
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0<GFS,CC> LS;      // ilu0
#endif

#ifdef USE_OVLP_BCGS_SuperLU_TPE_DG
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS,CC> LS; // superlu
#endif

#ifdef USE_OVLP_BCGS_SSORk_TPE_DG
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;  // ssork
#endif

#ifdef USE_NOVLP_PAR_TPE_DG
      typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GOS> LS; // SSOR with k
      //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS; // <-- no preconditioner
#endif

#ifdef USE_OVLP_GMRES_ILU0_TPE_DG
      typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<GFS,CC> LS;
#endif

#ifdef USE_OVLP_GMRES_SSORk_TPE_DG
      typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_SSORk<GFS,CC> LS;
#endif



    private:
      const GFS& gfs;
      GV gv; // Using a copy of the gridview is cheap and safe.
      VCType xSolution;
      const IDT& inputdata;
      const SDT& setupdata;
      const DGF_DARCY& darcyflux;
      const TP tp;
      const BCType bctype;
      //const BCExtension& dirichlet_function;
      SOURCE_TYPE& sourceterm;

      const Passenger::Type passenger;
      const EQ::Mode equationMode;

      const ConvectionDiffusionDGMethod::Type method;
      const ConvectionDiffusionDGWeights::Type weights;
      const REAL gamma;

      int GP_indicator;
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
      TPE_DG( const GFS& gfs_
              , const IDT& inputdata_
              , const SDT& setupdata_
              , const DGF_DARCY& darcyflux_
              , const TP& tp_
              , SOURCE_TYPE& sourceterm_ // non-const because sourceterm changes with measure location
              , const Passenger::Type passenger_=Passenger::solute
              , const EQ::Mode equationMode_=EQ::forward
              , const ConvectionDiffusionDGMethod::Type method_
              = ConvectionDiffusionDGMethod::SIPG
              , const ConvectionDiffusionDGWeights::Type weights_
              = ConvectionDiffusionDGWeights::weightsOn
              , const REAL gamma_=2.0
              , const int& GP_indicator_=-1
              , const bool& bReUseMatrixHierarchy_=false
              )
      :
      gfs( gfs_ )
        , gv( gfs.gridView() )
        , xSolution( gfs, 0.0 )
        , inputdata( inputdata_ )
        , setupdata( setupdata_ )
        , darcyflux( darcyflux_ )
        , tp( tp_ )
        , bctype( gv, tp )
        , sourceterm( sourceterm_ )
        , passenger( passenger_ )
        , equationMode( equationMode_ )
        , method(method_)
        , weights(weights_)
        , gamma(gamma_)
        , GP_indicator(GP_indicator_)
        , bReUseMatrixHierarchy( false )
        , lop( inputdata,
               setupdata,
               tp,
               sourceterm,
               passenger,
               equationMode,
               method,
               weights,
               gamma
               ) // Q1, P1
        , mbe(9)
        , gos( gfs, cc, gfs, cc, lop, mbe )
        , stiffness_matrix( gos )
        , ls_verbosity( std::max(0, inputdata.verbosity-VERBOSITY_EQ_DETAILS) )
        , istl_max_iter( std::min(5000, inputdata.transport_parameters.istl_max_iter) )

        // sequential solvers:
#ifdef USE_SUPERLU_SEQ_TPE_DG
        , linearsolver( ls_verbosity )
#endif

#ifdef USE_AMG_SEQ_TPE_DG
        , linearsolver( istl_max_iter, ls_verbosity, bReUseMatrixHierarchy )
#endif

#ifdef USE_BICGSTAB_SEQ_TPE_DG
        , linearsolver( istl_max_iter, ls_verbosity )
#endif

        // parallel solvers:

#ifdef USE_OVLP_BCGS_ILU0_TPE_DG
        , linearsolver( gfs, cc, istl_max_iter, ls_verbosity )
#endif

#ifdef USE_OVLP_BCGS_SuperLU_TPE_DG
        , linearsolver( gfs, cc, istl_max_iter, ls_verbosity )
#endif

#ifdef USE_OVLP_BCGS_SSORk_TPE_DG
        , linearsolver( gfs, cc, istl_max_iter, 5, ls_verbosity )
#endif

#ifdef USE_NOVLP_PAR_TPE_DG
        , linearsolver(gfs, istl_max_iter, 2, ls_verbosity)  // SSOR with k=2
        //, linearsolver(gfs, istl_max_iter, ls_verbosity)  // <-- no preconditioner
#endif

#ifdef USE_OVLP_GMRES_ILU0_TPE_DG
        , linearsolver( gfs, cc, istl_max_iter, ls_verbosity, 100, true ) // gemres + ilu0
#endif

#ifdef USE_OVLP_GMRES_SSORk_TPE_DG
        , linearsolver( gfs, cc, istl_max_iter, ls_verbosity, 5, 100, true ) // gemres + ssork
#endif



      {
        int processor_rank=gos.trialGridFunctionSpace().gridView().comm().rank();

        std::string eqMode = (equationMode==EQ::forward) ? "forward" : "adjoint";
        std::string eqType = (passenger==Passenger::solute) ? "solute" : "heat";
        std::string logMsg = "=== DG: New TPE: " + eqMode + " " + eqType;
        logger << std::endl << logMsg << std::endl;
        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << logMsg << std::endl;

        //if(inputdata.verbosity>8){
        //  logger<<"GP_indicator = "<<GP_indicator<<std::endl;
        //}


        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << "=== DG: call constraints assembler" << std::endl;
        // This call is very important to make parallel communication work!
        Dune::PDELab::constraints(bctype,gfs,cc);
        //Dune::PDELab::constraints(dirichlet_function,gfs,cc,false);


        UINT dofs = gfs.size();
        UINT cdofs = cc.size();
        if( inputdata.verbosity>=VERBOSITY_EQ_DETAILS ) {
          logger << "=== DG:  dofs = " << dofs << std::endl;
          logger << "=== DG: cdofs = " << cdofs << std::endl;
        }

        dofs = gv.comm().sum( dofs );
        cdofs = gv.comm().sum( cdofs );
        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS ) {
          std::cout << "=== DG:  dofs (all procs) = " << dofs << std::endl;
          std::cout << "=== DG: cdofs (all procs) = " << cdofs << std::endl;
        }


        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << "=== DG: matrix pattern statistics: "
                    << std::endl
                    << stiffness_matrix.patternStatistics()
                    << std::endl;
        stiffness_matrix = 0.0;


        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << "=== DG: Matrix assembly..." << std::endl;

        gos.jacobian( xSolution, stiffness_matrix );


#ifdef DEBUG_TPE_MATRIX
        if( stiffness_matrix.base().N() < 1000 )
          Dune::writeMatrixToMatlab( stiffness_matrix.base(), "tpe_dg_stiffness_matrix.dat" );
        else
          std::cout << "Stiffness matrix has row size " << stiffness_matrix.base().N() << std::endl;
#endif // DEBUG_TPE_MATRIX

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


      template<typename UType>
      void set_rhs( const UType& vc_solution )
      {
        lop.set_rhs( vc_solution );
      }

      template<typename UType>
      void set_rhs( const UType& vc_1, const UType& vc_2 )
      {
        lop.set_rhs( vc_1,vc_2 );
      }


      void set_PointSourceLocation( const DomainType& pointsource_global )
      {
        lop.set_PointSourceLocation( pointsource_global );
      }


      void solve_adjoint( const DomainType& pointsource_global,
                          VCType& xSolutionOut )
      {

        std::stringstream strlogMsg;
        strlogMsg << "=== DG: Solving adjoint TPE with pointsource at location " << pointsource_global;

        std::string logMsg = strlogMsg.str();
        logger << std::endl << logMsg << std::endl;
        int processor_rank=gos.trialGridFunctionSpace().gridView().comm().rank();
        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << logMsg << std::endl;

        assert( equationMode == EQ::adjoint );
        reset_BoundaryConditions();  // clean up cached solution from last measuring point!
        set_PointSourceLocation( pointsource_global );
        solve_forward( xSolutionOut );
      }

      template<typename UType>
      void solve_adjoint( const UType& vcphi_adj
                          , const UType& vc_phi0
                          , VCType& xSolutionOut )
      {
        assert( equationMode == EQ::adjoint );
        reset_BoundaryConditions(); // clean up cached solution from last measuring point!
        set_rhs( vcphi_adj, vc_phi0);
        solve_forward( xSolutionOut );
      }

      template<typename UType>
      void solve_adjoint( const UType& vcM1_adj
                          , VCType& xSolutionOut )
      {
        assert( equationMode == EQ::adjoint );
        reset_BoundaryConditions(); // clean up cached solution from last measuring point!
        set_rhs( vcM1_adj );
        solve_forward( xSolutionOut );
      }

      bool solve_forward( VCType& xSolutionOut
                          , bool bReUseSolution=false // if set to true, re-use xSolutionOut
                          )
      {
        int processor_rank=gos.trialGridFunctionSpace().gridView().comm().rank();

        std::string eqMode = (equationMode==EQ::forward) ? "forward" : "adjoint";
        std::string eqType = (passenger==Passenger::solute) ? "solute" : "heat";
        std::string logMsg = "=== DG: Solving TPE: " + eqMode + " " + eqType;
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

        watch.reset();
        // Solve Transport Equation here...
        REAL reduction = inputdata.transport_parameters.istl_reduction;

        VCType z(gos.trialGridFunctionSpace(), 0.0);

#ifdef USE_SECURE_MODE
        // Take care: NOVLP_AMG Solver modifies the stiffness matrix!
        GLOBAL_MATRIX tmpStiffnessMatrix(gos);
        tmpStiffnessMatrix=0.0;
        tmpStiffnessMatrix.base() = stiffness_matrix.base();
        linearsolver.apply( tmpStiffnessMatrix, z, r, reduction );
#else
        linearsolver.apply( stiffness_matrix, z, r, reduction ); // solver makes right hand side consistent
#endif

        ls_result = linearsolver.result();


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
        std::string msg = ls_result.converged ? "=== DG: Converged with" : "=== DG: WARNING: Not converged! ";

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


      bool isAcceptableUpDown( const VCType& vc_m,
                               const REAL maxValue=0.0,
                               const REAL minValue=0.0,
                               const REAL tolerance=5.0 // in percent
                               ){

        REAL minimum = 1E+100;
        REAL maximum = -1E+100;
        Dune::Gesis::General::logMinAndMax( vc_m, gv.comm(), minimum, maximum );

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
            std::cout << "Overshoots and undershoots within " << tolerance << " percent tolerance!" << std::endl;
          return true;
        }
        else
          return false;
      }

    };

  }
}

#endif// DUNE_GESIS_TPE_DG_HH
