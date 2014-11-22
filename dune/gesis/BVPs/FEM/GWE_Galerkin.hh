#ifndef DUNE_GESIS_GWE_GALERKIN_HH
#define DUNE_GESIS_GWE_GALERKIN_HH

#include<dune/pdelab/gridoperator/gridoperator.hh>

#include "GroundWaterOperatorGalerkin.hh"

namespace Dune {
  namespace Gesis {

    template<typename GFS
             , typename GWP
             , typename SOURCE_TYPE
             , typename IDT
             , typename SDT
             >
    class GWE_Galerkin
    {
      // extract types from GFS:
      typedef typename GFS::Traits::GridViewType GV;
      enum{dim=GV::dimension};
      typedef Dune::FieldVector<REAL,dim> DomainType;

      typedef typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type VCType;

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<GWP> BCType;
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<GWP> BCExtension;

      typedef Dune::Gesis::GwFlowOperator<GWP,BCType,SOURCE_TYPE,IDT,SDT> LOP;
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


      // set up solver
#ifdef USE_SUPERLU_SEQ
      typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
#endif

#ifdef USE_AMG_SEQ
      typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GOS> LS;
#endif

#ifdef USE_NOVLP_AMG_PAR
      typedef Dune::PDELab::ISTLBackend_NOVLP_CG_AMG_SSOR<GOS> LS;
#endif

#ifdef USE_OVLP_AMG_PAR_GWE_CG
      typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOS> LS;
#endif

#ifdef USE_OVLP_PAR_GWE_CG
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS, CC> LS;
#endif

#ifdef USE_BICGSTAB_SEQ
      typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GOS> LS;
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
#endif


    private:
      const GFS& gfs;
      const GV& gv;

      const IDT& inputdata;
      const SDT& setupdata;

      GWP& gwp;

      SOURCE_TYPE& sourceterm;
      const EQ::Mode equationMode;

      bool bReUseMatrixHierarchy; // flag telling the parallel AMG solver to re-use grid coarsening and the setup of the matrix hierarchy in case of re-use of the stiffness matrix

      VCType xSolution;

      BCType bctype;
      BCExtension dirichlet_function;

      CC cc; // Take care: CC cc must be declared before GOS gos! Otherwise you will get unprecdictable crashes! (cc: container for transformation)
      LOP lop;
      MBE mbe;
      GOS gos;
      GLOBAL_MATRIX stiffness_matrix;
      UINT ls_verbosity;
      LS linearsolver;
      Dune::PDELab::LinearSolverResult<REAL> ls_result;

      UINT iCounter; // used to count repeated solving of the same test problem

    public:

      // constructor:
      GWE_Galerkin( const GFS& gfs_
                    , const IDT& inputdata_
                    , const SDT& setupdata_
                    , GWP& gwp_
                    , SOURCE_TYPE& sourceterm_
                    , const EQ::Mode equationMode_=EQ::forward
                    , const bool& bReUseMatrixHierarchy_=false
                    )
        :
        gfs( gfs_ )
        , gv( gfs.gridView() )
	, inputdata( inputdata_ )
	, setupdata( setupdata_ )
	, gwp( gwp_ )
	, sourceterm( sourceterm_ )
        , equationMode(equationMode_)
        , bReUseMatrixHierarchy( bReUseMatrixHierarchy_ )
	, xSolution( gfs, 0.0 )
        , bctype( gv, gwp )
	, dirichlet_function( gv, gwp )
	, lop( gwp, bctype, sourceterm, inputdata, setupdata, equationMode ) // Q1, P1
        , mbe(9)
        , gos( gfs, cc, gfs, cc, lop, mbe )
	, stiffness_matrix( gos )
        , ls_verbosity( std::max(0, inputdata.verbosity-VERBOSITY_EQ_DETAILS) )

#ifdef USE_SUPERLU_SEQ
        , linearsolver( ls_verbosity )
#endif

#ifdef USE_AMG_SEQ
	, linearsolver( 5000, ls_verbosity, bReUseMatrixHierarchy )
#endif

#ifdef USE_OVLP_AMG_PAR_GWE_CG
        , linearsolver( gfs, 5000, ls_verbosity, bReUseMatrixHierarchy )
#endif

#ifdef USE_NOVLP_AMG_PAR
        , linearsolver( gfs, 5000, ls_verbosity, bReUseMatrixHierarchy )
#endif

#ifdef USE_OVLP_PAR_GWE_CG
        , linearsolver(gfs, cc, 5000, 5, ls_verbosity)
#endif

#ifdef USE_BICGSTAB_SEQ
        , linearsolver( 5000, ls_verbosity )
#endif

      {

        int processor_rank = gos.trialGridFunctionSpace().gridView().comm().rank();

        std::string eqMode = (equationMode==EQ::adjoint) ? "adjoint" : "forward";
        std::string logMsg = "=== FEM: New elliptic equation: " + eqMode;

        if( sourceterm.source_nature != GEOELECTRICAL_POTENTIAL_SOURCE &&
            sourceterm.source_nature != GEP_FUNCTIONAL_SOURCE &&
            sourceterm.source_nature != GEOELECTRICAL_POINT_SOURCE ){
          logMsg += " GWE";
        }
        else {
          logMsg += " GEP";
        }

        logger << std::endl << logMsg << std::endl;
        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << std::endl << logMsg << std::endl;

#ifdef USE_OVLP_AMG_PAR_GWE_CG
        Dune::Amg::Parameters amg_params;
        amg_params.setAccumulate( Dune::Amg::atOnceAccu );
        amg_params.setDebugLevel( ls_verbosity );
        linearsolver.setParameters( amg_params );
        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << "=== FEM: AMG with atOnceAccu" << std::endl;
#endif


        Dune::PDELab::constraints( bctype, gfs, cc); // fill container

        UINT dofs = gfs.size();
        UINT cdofs = cc.size();
        if( inputdata.verbosity>=VERBOSITY_EQ_DETAILS ) {
          logger << "=== FEM:  dofs = " << dofs << std::endl;
          logger << "=== FEM: cdofs = " << cdofs << std::endl;
        }

        dofs = gv.comm().sum( dofs );
        cdofs = gv.comm().sum( cdofs );
        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS ) {
          std::cout << "=== FEM:  dofs (all procs) = " << dofs << std::endl;
          std::cout << "=== FEM: cdofs (all procs) = " << cdofs << std::endl;
        }

        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << "=== FEM: interpolate dirichlet_function" << std::endl;
        Dune::PDELab::interpolate( dirichlet_function, gfs, xSolution );

        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << "=== FEM: set_nonconstrained_dofs( 0.0 ) " << std::endl;

        Dune::PDELab::set_nonconstrained_dofs( cc, 0.0, xSolution ); // clear interior


        stiffness_matrix = 0.0;
        gos.jacobian( xSolution, stiffness_matrix );


#ifdef USE_OVLP_AMG_PAR_GWE_CG
        std::string recycle = bReUseMatrixHierarchy ? "Recycle AMG coarsening hierarchy" : "No Recycling of AMG coarsening hierarchy, maybe not needed here.";
        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << "=== FEM: " << recycle << std::endl;
#endif
      }


      void reset_BoundaryConditions()
      {
        xSolution = 0.0;
      }

      /*
        bool set_SourceNature( const SourceNature sn )
        {
        lop.set_SourceNature( sn );
        }
      */

      template<typename UType>
      void set_rhs( const UType& vc1 )
      {
        lop.set_rhs( vc1 );
      }

      template<typename UType>
      void set_rhs( const UType& vc1, // eqn (66): m_k
                    const UType& vc2  // eqn (66): phi0
                    )
      {
        lop.set_rhs( vc1, vc2 );
      }

      template<typename UType>
      void set_rhs( const UType& vc1, const VCType& vc2 )
      {
        lop.set_rhs( vc1, vc2 );
      }

      template<typename UType>
      void set_rhs( const UType& vc1, const UType& vc2, const UType& vc3, const UType& vc4 )
      {
        lop.set_rhs( vc1, vc2, vc3, vc4 );
      }

      void set_PointSourceLocation( const DomainType& pointsource_global )
      {
        lop.set_PointSourceLocation( pointsource_global );
      }

      void set_PointSourceLocation( const DomainType& ps1_global, const DomainType& ps2_global)
      {
        lop.set_PointSourceLocation( ps1_global,ps2_global );
      }

      // For debugging purposes only:
      // This method is used to repeat solving the same test problem and log the output in r and A to find out differences
      void solve_adjoint( const UINT& iCounter_
                          , const DomainType& pointsource_global
                          , VCType& xSolutionOut )
      {
        iCounter = iCounter_;
        solve_adjoint( pointsource_global, xSolutionOut );
      }

      void solve_adjoint( const DomainType& pointsource_global
                          , VCType& xSolutionOut )
      {
        assert( equationMode==EQ::adjoint );
        set_PointSourceLocation( pointsource_global );
        //set_SourceNature( POINT_SOURCE );
        reset_BoundaryConditions();
        solve_forward( xSolutionOut );
      }

      void solve_adjoint( const DomainType& ps1_global
                          , const DomainType& ps2_global
                          , VCType& xSolutionOut )
      {
        assert( equationMode==EQ::adjoint );
        assert( sourceterm.source_nature==GEOELECTRICAL_POINT_SOURCE);
        set_PointSourceLocation( ps1_global,  ps2_global);
        //set_SourceNature( POINT_SOURCE );
        reset_BoundaryConditions();
        solve_forward( xSolutionOut );
      }


      template<typename UType>
      void solve_adjoint( const UType& vc_m1_adj
                          , VCType& xSolutionOut )
      {
        assert( equationMode==EQ::adjoint );
        reset_BoundaryConditions();
        set_rhs( vc_m1_adj );
        solve_forward( xSolutionOut );
      }


      template<typename UType>
      void solve_adjoint( const UType& vc_m0_adj
                          , const UType& vc_m0
                          , VCType& xSolutionOut )
      {
        assert( equationMode==EQ::adjoint );
        reset_BoundaryConditions();
        set_rhs( vc_m0_adj, vc_m0 );
        solve_forward( xSolutionOut );
      }



      template<typename UType>
      void solve_adjoint( const UType& vc_m0_adj
                          , const UType& vc_m0
                          , const UType& vc_m1_adj
                          , const UType& vc_m1
                          , VCType& xSolutionOut
                          )
      {
        assert( equationMode==EQ::adjoint );
        reset_BoundaryConditions();
        set_rhs( vc_m0_adj, vc_m0, vc_m1_adj, vc_m1 );
        solve_forward( xSolutionOut );
      }


      template<typename UType>
      void solve_forward(
                         const UType& vc_phi0
                         , const VCType& vc_m0
                         , VCType& xSolutionOut
                         )
      {
        reset_BoundaryConditions();
        set_rhs( vc_phi0, vc_m0 );
        solve_forward( xSolutionOut );
      }

      bool solve_forward(
                         VCType& xSolutionOut
                          )
      {

        int processor_rank=gos.trialGridFunctionSpace().gridView().comm().rank();

        std::string eqMode = (equationMode==EQ::adjoint) ? "adjoint" : "forward";
        std::string logMsg = "=== FEM: Solving elliptic equation: " + eqMode;

        if( sourceterm.source_nature != GEOELECTRICAL_POTENTIAL_SOURCE &&
            sourceterm.source_nature != GEP_FUNCTIONAL_SOURCE &&
            sourceterm.source_nature != GEOELECTRICAL_POINT_SOURCE ){
          logMsg += " GWE";
        }
        else {
          logMsg += " GEP";
        }

        logger << std::endl << logMsg << std::endl;
        if( processor_rank==0 && inputdata.verbosity>=VERBOSITY_EQ_DETAILS )
          std::cout << std::endl << logMsg << std::endl;


        // Class for evaluating the boundary function for the velocity ( needed only for RT, = 0 here)
        //V_T<GV, REAL> v(gv);

#ifdef DEBUG_REPLOT
        // Get the conductivity field on the grid for the vtkwriter
        typedef typename GV::Grid GRIDTYPE; // get the grid-type out of the gridview-type
        typedef typename GV::Grid::template Codim < 0 > ::Entity Entity;
        typedef typename GV::Grid::template Codim < 0 > ::EntityPointer EntityPointer;
        typedef typename GV::template Codim < 0 > ::Iterator ElementLeafIterator;
        typedef typename ElementLeafIterator::Entity::Geometry LeafGeometry;

        // make a mapper for codim 0 entities in the leaf grid
        Dune::LeafMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout>
	  mapper(gv.grid()); // get the underlying hierarchical grid out ot the gridview

        const int nGridCells = mapper.size();
        int rank;
	MPI_Comm_rank(gv.comm(),&rank);
	logger << "P" << rank << ": mapper.size() = " << nGridCells << std::endl;
        const typename GV::IndexSet& indexset = gv.indexSet();

        std::vector<REAL> conductivity_field(nGridCells);
        for (ElementLeafIterator it = gv.template begin < 0 > (); it != gv.template end < 0 > (); ++it) {
          int id = indexset.index(*it);
          REAL kk = 0.0;
          Yfieldwrapper.getElementLogConductivity(it, kk);
          conductivity_field[id] = kk;
        }

        std::cout << "DEBUG-LOG: VTK output of k-field again..." << std::endl;
        Dune::VTKWriter<GV> vtkwriterK(gv);
        vtkwriterK.addCellData(conductivity_field, "conductivity");
        vtkwriterK.write("k_field_debug", Dune::VTKOptions::ascii);
#endif // DEBUG_REPLOT




        Dune::Timer watch;

        // assemble residual:
        typedef typename GOS::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
        typedef typename Dune::PDELab::BackendVectorSelector<TrialGridFunctionSpace,REAL>::Type W;

        W r(gos.testGridFunctionSpace(), 0.0);

        gos.residual(xSolution, r); // residual is additive

        General::log_elapsed_time( watch.elapsed(),
                                   gos.trialGridFunctionSpace().gridView().comm(),
                                   inputdata.verbosity,
                                   "GWE",
                                   "Residual Assembly" );

        watch.reset();

        // Solve Flow Equation here...
        REAL reduction = 1e-10;

        VCType z(gos.trialGridFunctionSpace(), 0.0);



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
                                   "GWE",
                                   "Linear Solver" );

        return true;

      };


      const Dune::PDELab::LinearSolverResult<REAL>& result() const{
        return ls_result;
      }


      std::string show_ls_result(){
        std::stringstream ss;
        std::string msg = ls_result.converged ? "=== FEM: Converged with" : "WARNING: Not converged! ";

        ss << msg
           << " T=" << ls_result.elapsed
           << " IT=" << ls_result.iterations
           << " TIT=" << ls_result.elapsed/ls_result.iterations
           << " reduction=" << ls_result.reduction
           << " rate=" << ls_result.conv_rate;
        return ss.str();
      }

    };

  }
}

#endif// DUNE_GESIS_GWE_GALERKIN_HH
