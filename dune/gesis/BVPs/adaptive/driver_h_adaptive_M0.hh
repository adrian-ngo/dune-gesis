#ifndef DUNE_GESIS_ADAPTIVE_DRIVE_HH
#define DUNE_GESIS_ADAPTIVE_DRIVE_HH

#include <dune/gesis/BVPs/MeshPlecletNumber.hh>
#include <dune/gesis/BVPs/sourceterms/source_head_transport.hh>
#include <dune/gesis/BVPs/GroundWaterParameters.hh>
#include <dune/gesis/BVPs/DG/GWE_CCFV.hh>
#include <dune/gesis/BVPs/TransportParameters.hh>
#include <dune/gesis/BVPs/DG/TransportOperatorDG.hh>
#include <dune/gesis/BVPs/DG/TPE_DG.hh>

#include "ee_DG_h_adaptive.hh"  // residual error estimator
#ifdef USE_ALUGRID
//#include "RedistributeDataHandle.hh"  // redistributing vc_h after grid load-balancing
//#include <dune/gesis/BVPs/projectors/CoarseGridP0Datahandle.hh>
#endif
#include "adaptivity.hh"

//#include <dune/gesis/common/differenceAdapters.hh>
//#include "../integrateFunction.hh"
//#include "../sortedplot.hh"

//#include "rgv_adapt.hh"

namespace Dune {
  namespace Gesis {

    template<
      typename GRID,
      typename GFS_GW,
      typename IDT,
      typename YFG,
      typename DIR
      >
    void hAdaptiveLoopForM0( GRID& grid,
                              GFS_GW& gfs_gw,  // <-- changed after refinement and load-balancing!
                              const IDT& inputdata,
                              YFG& yfg, // <-- changed after refinement and load-balancing!
                              const UINT baselevel,
                              const DIR& dir,
                              const Dune::MPIHelper& helper
                              ){

      typedef typename IDT::SDT SDT;
      const SDT& setupdata = inputdata.setups[0];

      Dune::Timer watch;
      //double elapsed_time;

      logger << "Get Gridview for the GWE ..." << std::endl;
      typedef typename GFS_GW::Traits::GridViewType GV_GW;
      enum{ dim = GV_GW::dimension };
      GV_GW gv_gw = gfs_gw.gridView();

      // Setup and solve GWE:
      typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD;
      GWP_FWD gwp_fwd( inputdata, setupdata, yfg );

      typedef PumpingSource<GV_GW,REAL,SDT> PumpingSourceTypeGW;
      PumpingSourceTypeGW source_h( setupdata );

      typedef GroundWaterEquation<GFS_GW,GWP_FWD,PumpingSourceTypeGW,IDT,SDT> GWE_H;

      GWE_H gwe_h( gfs_gw,
                   inputdata,
                   setupdata,
                   gwp_fwd,
                   source_h );

      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;
      VCType_GW vc_head( gfs_gw, 0.0 );
      gwe_h.solve_forward( vc_head );

      if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
        std::cout << gwe_h.show_ls_result() << std::endl;
      General::logMinAndMax( vc_head, gv_gw.comm() );

      VTKPlot::output2vtu( gfs_gw,
                           vc_head,
                           dir.vtudir + "/test_h_orig",
                           "h_orig",
                           -1 );

      typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_DGF;
      DARCY_FLUX_DGF darcyflux_dgf( gwp_fwd,
                                    gfs_gw,
                                    vc_head,
                                    baselevel );

      VTKPlot::output_dgf_to_vtu( gv_gw,
                                  gfs_gw,
                                  darcyflux_dgf,
                                  dir.vtudir + "/test_q_orig",
                                  "q_orig",
                                  1 );


#ifdef NEW_GV

#ifdef USE_DGF_PressureField

      logger << "Using DGF_PressureField ... " << std::endl;
      typedef DGF_PressureField<GFS_GW,VCType_GW> RT0_PF;
      Dune::shared_ptr<RT0_PF> rt0_pressurefield
        = Dune::make_shared<RT0_PF>(gfs_gw,vc_head,baselevel);

#else

      logger << "Using RT0_PressureField ... " << std::endl;
      typedef RT0_PressureField<GV_GW,GFS_GW,GWP_FWD> RT0_PF;
      Dune::shared_ptr<RT0_PF> rt0_pressurefield
        = Dune::make_shared<RT0_PF>( gv_gw,
                                     gwp_fwd,
                                     baselevel );
      rt0_pressurefield->assemble(gfs_gw,vc_head);

#endif

      PressureLikeOrdering comparefunctor;

      typedef typename GRID::LeafGridView OLD_GV_TP;
      //typedef typename GRID::LevelGridView OLD_GV_TP;
      typedef ReorderedGridView<OLD_GV_TP,RT0_PF,PressureLikeOrdering> GV_TP;

      // Grid elements will be sorted according to the appropriate hydraulic head.

#else

      //====================================================================
      // Use default (lexicographical) sorting of grid elements for SDFEM!
      // Resorting would not reduce the number of iterations for the
      // AMG based linear solver.
      //====================================================================
      typedef typename GRID::LeafGridView GV_TP;

#endif


      GV_TP gv_tp( grid.leafGridView(), rt0_pressurefield, comparefunctor );

      //Dune::shared_ptr<GV_TP> pgv_tp
      //  = Dune::make_shared<GV_TP>( grid.leafGridView(),
      //                              //grid.levelView(0),
      //                              rt0_pressurefield,
      //                              comparefunctor );


#ifdef PARALLEL
      // This one must be used for parallel ALUGRID (Rebecca).
      // Ghost cells are treated as overlap cells (overlap=1) for the linear solver.
      typedef Dune::PDELab::P0ParallelGhostConstraints CONSTRAINTS;
#else
      // This one must be used in sequential mode.
      typedef Dune::PDELab::NoConstraints CONSTRAINTS;
#endif // PARALLEL

      CONSTRAINTS con_tp;

      typedef Dune::PDELab::QkDGLocalFiniteElementMap<REAL,REAL,pMAX,dim> FEM_HYPER;
      FEM_HYPER fem_hyper;

      const int blocksize = Dune::QkStuff::QkSize<pMAX,dim>::value;

      if(helper.rank()==0)
        std::cout << "Using QkDG for the hyperbolic PDE with blocksize "
                  << blocksize << std::endl;
      typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking,blocksize> VBE_TP;

      typedef Dune::PDELab::GridFunctionSpace<GV_TP,FEM_HYPER,CONSTRAINTS,VBE_TP> GFS_TP;
      GFS_TP gfs_tp( gv_tp, fem_hyper, con_tp );
      //GFS_TP gfs_tp( *pgv_tp, fem_hyper, con_tp );

      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;


#if defined L2ProjectionOfM0
      Dune::GeometryType gt;
      gt.makeCube(dim);
      typedef Dune::PDELab::P0LocalFiniteElementMap<CTYPE,REAL,dim> P0FEM;
      P0FEM p0fem(gt);
      typedef Dune::PDELab::ISTLVectorBackend<> VBE1;
      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE1;
      MBE1 mbe1(9);

      typedef Dune::PDELab::GridFunctionSpace<GV_TP,P0FEM,Dune::PDELab::NoConstraints,VBE1> P0GFS;
      typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,REAL>::Type U0Type;
#endif


      //#ifdef BACKUP_REPLAY
      VCType_TP vc_m0( gfs_tp, 0.0 );
      //#endif


      // ------------------------
      // Adaptive refinement loop
      // ------------------------

#ifdef USE_UG
      UINT maxsteps = inputdata.domain_data.ug_maxsteps;
#endif

#ifdef USE_ALUGRID
      UINT maxsteps = inputdata.domain_data.alugrid_maxsteps;
#endif


      //std::vector< Dune::shared_ptr<GV_TP> > pool_gv_tp;
      //std::vector< GV_TP > pool_gv_tp;
      //std::vector< std::vector<REAL> > pool_std_vc_m0;


      // some arrays to store results
#ifdef STORE_SOLUTION_HIERARCHY
      std::vector<double> l2e_hierarchy;
#endif

      std::vector<double> h1s;
      std::vector<double> dge;
      std::vector<double> ee;
      std::vector<int> Ndofs;
      std::vector<int> nIterations; // of the linear solver
      std::vector<double> ls_elapsed; // total elapsed time of the linear solver
      std::vector<double> minima;
      std::vector<double> maxima;

#ifdef STORE_SOLUTION_HIERARCHY
      std::vector<std::vector<REAL > > solution_hierarchy(maxsteps);
#endif

      for(UINT step=0; step<maxsteps+1; step++){

        Ndofs.push_back( gfs_tp.globalSize() );
        std::cout << "grid.maxLevel() = " << grid.maxLevel() << std::endl;


        //std::cout << "pgv_tp->size(0) = " << pgv_tp->size(0) << std::endl;
        //pool_gv_tp.push_back( pgv_tp );
        //std::cout << "gv_tp.size(0) = " << gv_tp.size(0) << std::endl;
        //pool_gv_tp.push_back( gv_tp );  // <-- This is not possible!

        if(helper.rank()==0){

          std::cout << std::endl;
          std::cout << "------------------------------" << std::endl;
          std::cout << std::endl;
          std::cout << "Step " << step << std::endl;
          std::cout << std::endl;
          std::cout << "------------------------------" << std::endl;

        }

        logger << std::endl;
        logger << "------------------------------------" << std::endl;
        logger << "Refinement step " << step << std::endl;
        logger << "gfs_tp.globalSize() = " << gfs_tp.globalSize() << std::endl;
        logger << "    grid.maxLevel() = " << grid.maxLevel() << std::endl;
        logger << "------------------------------------" << std::endl;

        if(step<=maxsteps){
          /*
            std::stringstream vtu_Y;
            vtu_Y << dir.vtudir << "/Y_step" << step;
            gwp_fwd.getYfieldgenerator().plot2vtu( gv_gw, vtu_Y.str(), baselevel );

            std::stringstream vtu_h;
            vtu_h  << dir.vtudir << "/h_step" << step;
            Dune::Gesis::VTKPlot::output2vtu( gfs_gw.gridView(),
            gfs_gw,
            vc_h,
            vtu_h.str(),
            "h",
            inputdata.verbosity,
            true,
            std::max(0,pMAX-1) );

            std::stringstream vtu_hLeaf;
            vtu_hLeaf << dir.vtudir << "/h_leaf_0" << step;
            rt0_pressurefield->plot2vtu( gv_tp, vtu_hLeaf.str() );
          */


          if( inputdata.plot_options.vtk_plot_element_ordering ){
            std::stringstream vtu_gv_tp;
            vtu_gv_tp << dir.vtudir << "/gv_tp_orig_" << step;
            VTKPlot::outputGridviewIndexToDGF( gv_tp, vtu_gv_tp.str() );
          }
        }

        // gfs_gw and vc_h updated after each refinement step!


#ifndef TEST_Adjoint
        typedef SolutePumpingSource<GV_TP,REAL,IDT,SDT> SourceTypeTP;
        typedef TransportProblemM0<GV_TP,REAL,DARCY_FLUX_DGF,IDT,SDT> TPM0;
        SourceTypeTP source_m0( inputdata, setupdata );
#else
        typedef TransportProblemAdjoint<GV_TP,REAL,DARCY_FLUX_DGF,IDT,SDT> TPM0;
        typedef PointFunctionSource<GFS_TP,IDT> SourceTypeTP;
        SourceTypeTP source_m0( gfs_tp,
                                inputdata,
                                -1.0, // negative pointsource
                                setupdata.solute_concentration_inversion_data.regularization_factor,
                                setupdata.solute_concentration_inversion_data.bfixedwidth );
#endif

        TPM0 tpm0(darcyflux_dgf,inputdata,setupdata);

        //MeshPlecletNumber<GV_TP,TPM0> meshPeclet( *pgv_tp, tpm0, step );
        MeshPlecletNumber<GV_TP,TPM0> meshPeclet( gv_tp, tpm0, step );
        std::stringstream vtu_Peclet;
        vtu_Peclet << dir.vtudir << "/peclet_step" << step;
        meshPeclet.plot2vtu( vtu_Peclet.str() );
        if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
          std::cout << "Maximal mesh Peclet numbers: Pe_l = "
                    << meshPeclet.maximum_l() << " / Pe_t = "
                    << meshPeclet.maximum_t() << std::endl;
        }


        typedef TransportEquation
          < GFS_TP
            , DARCY_FLUX_DGF
            , SourceTypeTP
            , TPM0
            , IDT
            , SDT
            > TPE_M0;

        TPE_M0 tpe_m0( gfs_tp,
                       inputdata,
                       setupdata,
                       darcyflux_dgf,
                       tpm0,
                       source_m0,
                       Passenger::solute,
#ifndef TEST_Adjoint
                       EQ::forward
#else
                       EQ::adjoint
#endif
                       );

#ifndef BACKUP_REPLAY
        //    VCType_TP vc_m0( gfs_tp, 0.0 );
#else
        if( step > 0 ){

          if( inputdata.plot_options.vtk_plot_m0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
            std::stringstream vtu_m0_file_before_solve;
            vtu_m0_file_before_solve  << dir.vtudir << "/m0_before_step" << step;

            VTKPlot::output2vtu( gfs_tp,
                                 vc_m0,
                                 vtu_m0_file_before_solve.str(),
                                 "m0_orig_before",
                                 pMAX-1 );
          }

        }
#endif


#ifndef TEST_Adjoint
        tpe_m0.solve_forward( vc_m0 );
#else
        Dune::FieldVector<REAL,dim> measure_location( 200.0 );
        measure_location[0] = setupdata.head_inversion_data.mplist.pointdata_vector[0].x;
        measure_location[1] = setupdata.head_inversion_data.mplist.pointdata_vector[0].y;
#ifdef DIMENSION3
        measure_location[2] = setupdata.head_inversion_data.mplist.pointdata_vector[0].z;
#endif
        tpe_m0.solve_adjoint( measure_location, vc_m0 );
#endif

        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << tpe_m0.show_ls_result() << std::endl;

        REAL u_minimum=1e100;
        REAL u_maximum=-1e100;
        General::logMinAndMax( vc_m0, gv_gw.comm(), u_minimum, u_maximum );
        minima.push_back( u_minimum );
        maxima.push_back( u_maximum );

        UINT allDOFs = gfs_tp.globalSize();
        allDOFs = gv_gw.comm().sum( allDOFs );

        if( gv_gw.comm().rank()==0 ){
          std::cout << "Step " << step
                    << " All DOFs = " << allDOFs
                    << " minimum = " << u_minimum
                    << " maximum = " << u_maximum
                    << std::endl;
          std::cout << tpe_m0.show_ls_result()
                    << std::endl;
        }

        //std::vector<REAL> std_vc_m0;
        //General::copy_to_std( vc_m0, std_vc_m0 );
        //pool_std_vc_m0.push_back( std_vc_m0 );

        if( inputdata.plot_options.vtk_plot_m0 ){
          std::stringstream vtu_m0_file;
          vtu_m0_file  << dir.vtudir << "/m0_step" << step;
          VTKPlot::output2vtu( gfs_tp,
                               vc_m0,
                               vtu_m0_file.str(),
                               "m0_orig",
                               pMAX-1 );

        }

        // *********************************
        // Plot over line for gnuplot:
        // *********************************

        /*
          std::stringstream gnuplot_datfile;
          gnuplot_datfile << dir.vtudir << "/p"
          << std::setfill('0') << std::setw(4) << helper.rank()
          << "-lineplot_m0_step" << step << ".dat";

          LinePlot<GV_TP> lplot(gv_tp);
          lplot.setEndPoints( inputdata.domain_data.lineplot.startpoint,
          inputdata.domain_data.lineplot.endpoint );
          lplot.addDGData( gfs_tp,
          vc_m0,
          "m0_adapt" );
          lplot.write( gnuplot_datfile.str() );
          std::cout << "Write to " << gnuplot_datfile.str() << std::endl;
        */

        if( tpe_m0.isAcceptableUpDown( vc_m0, 100.0, 0.0, inputdata.transport_parameters.tolerance ) )
          break;


        P0GFS p0gfs( gv_tp, p0fem );


#ifdef TEST_M1

        typedef FunctionSource<GWP_FWD,
                               GFS_TP,
                               IDT,SDT> FunctionSourceType;

        FunctionSourceType source_m1_m0( gwp_fwd,
                                         gfs_tp,
                                         inputdata,
                                         setupdata );

        typedef TransportProblemM1<GV_TP,REAL,DARCY_FLUX_DGF,IDT,SDT> TPM1;
        TPM1 tpm1(darcyflux_dgf,inputdata,setupdata);

        typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,FunctionSourceType,TPM1,IDT,SDT> TPE_M1_M0;

        TPE_M1_M0 tpe_m1_m0( gfs_tp,
                             inputdata,
                             setupdata,
                             darcyflux_dgf,
                             tpm1,
                             source_m1_m0 );



        // Define sourceterm for the forward transport equation for m1:
        VCType_TP vc_theta_m0 = vc_m0; // theta_m0 := porosity * m0
        // multiply a vector container with the porosity. needed for zonation!

        //General::vc_times_theta(vc_theta_m0,*pgv_tp,inputdata);
        General::vc_times_theta(vc_theta_m0,gv_tp,inputdata);

        tpe_m1_m0.set_rhs( vc_theta_m0 );


        VCType_TP vc_m1_m0( gfs_tp, 0.0 );
        tpe_m1_m0.solve_forward( vc_m1_m0 );

        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << tpe_m1_m0.show_ls_result() << std::endl;
        General::logMinAndMax( vc_m1_m0, gv_gw.comm() );

        std::stringstream vtu_m1m0_file;
        vtu_m1m0_file  << dir.vtudir << "/test_m1m0_step" << step;
        VTKPlot::output2vtu( gfs_tp,
                             vc_m1_m0,
                             vtu_m1m0_file.str(),
                             "m1_orig",
                             pMAX-1 );

        std::vector<REAL> vm0;
        General::copy_to_std( vc_m0, vm0 );
        std::vector<REAL> vm1m0;
        General::copy_to_std( vc_m1_m0, vm1m0 );

        std::vector<REAL> vmat;
        vmat.resize(vc_m0.flatsize());

        for(UINT i=0;i<vc_m0.flatsize();i++)
          {
            if( vm0[i] > 1e-3 )
              vmat[i] = vm1m0[i] / vm0[i];
          }

        VCType_TP vc_MeanArrivalTime( gfs_tp, 0.0 );
        vc_MeanArrivalTime.std_copy_from( vmat );

        std::stringstream vtu_mat_file;
        vtu_mat_file  << dir.vtudir << "/mat_step" << step;
        VTKPlot::output2vtu( gfs_tp,
                             vc_MeanArrivalTime,
                             vtu_mat_file.str(),
                             "MeanArrivalTime",
                             pMAX-1
                             );

#endif // TEST_M1

        if( maxsteps < 1 )
          break;

        typedef DGResidualEstimator<IDT,SDT,TPM0> ESTLOP;
        ESTLOP estlop(inputdata,setupdata,tpm0);

#define Compute_Error_Estimate

#ifdef Compute_Error_Estimate

        logger << std::endl;
        logger << "Calculate Error Estimate..." << std::endl;

        // Estimate Error for the last solution...
        /*
          Dune::GeometryType gt;
          gt.makeCube(dim);
          typedef Dune::PDELab::P0LocalFiniteElementMap<CTYPE,REAL,dim> P0FEM;
          P0FEM p0fem(gt);
          typedef Dune::PDELab::ISTLVectorBackend<> VBE1;
          typedef typename VBE1::MatrixBackend MBE1;

          typedef Dune::PDELab::GridFunctionSpace<GV_TP,P0FEM,Dune::PDELab::NoConstraints,VBE1> P0GFS;
          P0GFS p0gfs(gv_tp,p0fem);
        */
        typedef Dune::PDELab::PowerGridFunctionSpace
          < P0GFS
            , 2
            , VBE1
            , Dune::PDELab::LexicographicOrderingTag
            > P0PowerGFS;
        P0PowerGFS p0powergfs(p0gfs);

        typedef Dune::PDELab::EmptyTransformation NoTrafo;

        typedef Dune::PDELab::GridOperator
          < GFS_TP
            , P0PowerGFS
            , ESTLOP
            , MBE1
            , REAL
            , REAL
            , REAL
            , NoTrafo
            , NoTrafo
            , true // nonoverlapping_mode!!!
            > ESTGO;
        ESTGO estgo(gfs_tp,p0powergfs,estlop,mbe1);

        typedef typename Dune::PDELab::BackendVectorSelector<P0PowerGFS,REAL>::Type U0PowerType;
        U0PowerType x0Power( p0powergfs, 0.0 );
        estgo.residual( vc_m0, x0Power );   // <-- Run over grid cells and evaluate!

        //typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,REAL>::Type U0Type;
        U0Type eta(p0gfs,0.0); // local error indicator
        U0Type rho(p0gfs,0.0); // some other local indicator

        General::extract_two_indicators( p0powergfs, x0Power, p0gfs, eta, rho );

        // Compute global error estimate out of local error indicators

        for (unsigned int i=0; i<eta.N(); i++)
          eta.base()[i] = std::sqrt(eta.base()[i]);  // very important for error_fraction strategy!

        std::stringstream eta_filename;
        eta_filename << dir.vtudir << "/p"
                     << std::setfill('0') << std::setw(4) << helper.rank()
                     << "-sorted_eta_step" << step << ".dat";
        //    SortedPlot sortedplot;
        //sortedplot.output2dat( eta, eta_filename.str() );

        REAL estimated_error = eta.two_norm();
        //REAL estimated_error = std::sqrt( eta.one_norm() );

        ee.push_back(estimated_error);

        std::cout << ">>> Estimated Error = "
                  << estimated_error
                  << std::endl;

        if(true){

          std::vector<REAL> eta_flat;
          General::copy_to_std( eta, eta_flat );

          std::stringstream vtu_m0_eta_file;
          vtu_m0_eta_file  << dir.vtudir << "/eta_M0_step" << step;

          VTKPlot::output_vector_to_vtu( gv_tp,
                                         eta_flat,
                                         vtu_m0_eta_file.str(),
                                         "eta"
                                         );
        }

#endif // Compute_Error_Estimate

        // Prepare for adaptation...
        /*
          typedef Dune::PDELab::L2Projection<GFS_TP,VCType_TP> L2Projection;
          L2Projection l2projection;


          Dune::GeometryType::BasicType basic_geometry_type = Dune::GeometryType::cube;

          #ifdef Compute_Error_Estimate
          typedef Dune::PDELab::myResidualErrorEstimation<VCType_TP,U0Type> REE;
          REE ree(eta);
          #else
          typedef Dune::PDELab::ResidualErrorEstimation
          < GFS_TP
          , VCType_TP
          , ESTLOP
          > REE;
          REE ree(gfs_tp,estlop,basic_geometry_type);
          #endif
        */
        // do refinement

        UINT maxDOFs = 99123456;
        UINT allnewDOFs = gfs_tp.globalSize();
        allnewDOFs = gv_gw.comm().sum( allnewDOFs );

        //std::vector<Dune::shared_ptr<VCType_TP> > solution_hierarchy(maxsteps);

        if( step<maxsteps && allnewDOFs<maxDOFs ){

#ifdef USE_UG
          REAL refinementfraction = inputdata.domain_data.ug_refinementfraction;
          REAL coarseningfraction = inputdata.domain_data.ug_coarseningfraction;
          UINT strategy = inputdata.domain_data.ug_strategy;
#endif
#ifdef USE_ALUGRID
          REAL refinementfraction = inputdata.domain_data.alugrid_refinementfraction;
          REAL coarseningfraction = inputdata.domain_data.alugrid_coarseningfraction;
          UINT strategy = inputdata.domain_data.alugrid_strategy;
#endif
          REAL refine_threshold(0);
          REAL coarsen_threshold(0);


          //error_distribution( eta, 10 );

          if( strategy == 1){
            std::cout << "error fraction strategy was chosen" << std::endl;
            error_fraction( eta,
                            gv_tp,
                            refinementfraction,
                            coarseningfraction,
                            refine_threshold,
                            coarsen_threshold,
                            inputdata.verbosity );
          }
          else{
            std::cout << "element fraction strategy was chosen" << std::endl;
            element_fraction( eta,
                              refinementfraction,
                              coarseningfraction,
                              refine_threshold,
                              coarsen_threshold,
                              inputdata.verbosity );
          }
          mark_grid_rgv( grid,
                         gv_tp,
                         eta,
                         refine_threshold,
                         coarsen_threshold,
                         inputdata.verbosity );

          // prepare the grid for refinement
          grid.preAdapt();

#ifdef BACKUP_REPLAY
          typedef Dune::PDELab::L2Projection<GFS_TP,VCType_TP> Projection;
          Projection projection;
          typedef Dune::PDELab::GridAdaptor<GRID,GFS_TP,VCType_TP,Projection> GA;
          GA grid_adaptor;
          // save u
          typename GA::MapType transferMap1;
          grid_adaptor.backupData(grid,gfs_tp,projection,vc_m0,transferMap1);

          // save intermediate solution hierarchy

#ifdef STORE_SOLUTION_HIERARCHY
          General::copy_to_std( vc_m0, solution_hierarchy[step] );

          std::vector<typename GA::MapType> solutionMap(step+1);
          for(int ii=0;ii<step+1;++ii){
            VCType_TP tmpSolution(gfs_tp,0.0);

            //std::cout << "before backup: solution_hierarchy[" << ii << "].size() = "
            //          << solution_hierarchy[ii].size()
            //          << std::endl;

            tmpSolution.std_copy_from(solution_hierarchy[ii]);

            //std::cout << "before backup: tmpSolution.flatsize() = "
            //          << tmpSolution.flatsize()
            //          << std::endl;

            grid_adaptor.backupData(grid,gfs_tp,projection,tmpSolution,solutionMap[ii]);
          }
#endif// STORE_SOLUTION_HIERARCHY

#endif// BACKUP_REPLAY

          adapt_grid_rgv( grid,
                          gfs_gw,
                          vc_head,
                          gv_tp,
                          gfs_tp,
                          yfg,
                          inputdata,
                          dir,
                          step
#ifdef BACKUP_REPLAY
                          , vc_m0
#endif
                          );


#ifdef STORE_SOLUTION_HIERARCHY
          // reset u
          vc_m0 = VCType_TP(gfs_tp,0.0);
          grid_adaptor.replayData(grid,gfs_tp,projection,vc_m0,transferMap1);


          for(int ii=0;ii<step+1;++ii){
            VCType_TP solution( gfs_tp, 0.0 );
            grid_adaptor.replayData(grid,gfs_tp,projection,solution,solutionMap[ii]);
            General::copy_to_std( solution, solution_hierarchy[ii] );
            //std::cout << "after replay: solution.flatsize() = "
            //          << solution.flatsize()
            //          << std::endl;
            //std::cout << "after replay: solution_hierarchy[" << ii << "].size() = "
            //          << solution_hierarchy[ii].size()
            //          << std::endl;
          }
#endif // STORE_SOLUTION_HIERARCHY

          // clean up
          grid.postAdapt();

        }
        else {
          std::cout << "max. steps or max. DOFs reached at step " << maxsteps
                    << " with #dofs = " <<  gfs_tp.globalSize() << std::endl;
          break;
        }

        //std::cout << "GV step = " << step << std::endl;
        //std::cout << "GV size  = " << pgv_tp->size(0) << std::endl;
        //std::cout << "vc base size = " << vc_m0.base().size() << std::endl;
        //std::cout << "vc flatsize  = " << vc_m0.flatsize() << std::endl;

      } // end of refinement loop




#ifdef STORE_SOLUTION_HIERARCHY
      // dgf for the reference solution:
      typedef Dune::PDELab::DiscreteGridFunction<GFS_TP,VCType_TP> DGF_TP;
      DGF_TP udgf( gfs_tp, vc_m0 );

      // Compare solution hierarchy against reference solution
      std::vector<Dune::shared_ptr<VCType_TP> > pSol(maxsteps);
      for( int ii=0; ii<maxsteps; ++ii ){
        pSol[ii] = Dune::make_shared<VCType_TP>( gfs_tp, 0.0 );
        (pSol[ii])->std_copy_from( solution_hierarchy[ii] );
        DGF_TP udgf_ii( gfs_tp, *(pSol[ii]) );
        DifferenceSquaredAdapter<DGF_TP,DGF_TP> differencesquared(udgf,udgf_ii);
        Dune::FieldVector<REAL,1> l2errorsquared(0.0);
        std::cout << "Calculate l2-error on grid level " << ii << std::endl;
        Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,2*pMAX);
        std::cout << "l2error = " << sqrt(l2errorsquared) << std::endl;
        l2e_hierarchy.push_back(sqrt(l2errorsquared));


        // debug plots:
        if( inputdata.plot_options.vtk_plot_m0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
          std::stringstream vtu_m0_file;
          vtu_m0_file  << dir.vtudir << "/pSol_m0_step" << ii;
          VTKPlot::output2vtu( gfs_tp,
                               *(pSol[ii]),
                               vtu_m0_file.str(),
                               "m0_orig",
                               pMAX-1 );
        }


      }

      std::cout << "Set l2-error on grid level " << maxsteps
                << " to 1E-12." << std::endl;
      l2e_hierarchy.push_back( 1E-12 );

#endif //STORE_SOLUTION_HIERARCHY


#ifdef USE_MIRROR_ADAPT

      typedef Dune::PDELab::L2Projection<GFS_TP,VCType_TP> Projection;
      Projection projection;
      typedef Dune::PDELab::GridAdaptor<GRID,GFS_TP,VCType_TP,Projection> GA;
      GA grid_adaptor;

      for( int ibackstep=grid.maxLevel(); ibackstep>0; --ibackstep ){
        Dune::PDELab::mark_all_coarsen( grid );

        // prepare the grid for refinement
        grid.preAdapt();

        // save u
        typename GA::MapType transferMap1;
        grid_adaptor.backupData(grid,gfs_tp,projection,vc_m0,transferMap1);

        adapt_grid_rgv( grid,
                        gfs_gw,
                        vc_head,
                        gv_tp,
                        gfs_tp,
                        yfg,
                        inputdata,
                        dir,
                        ibackstep
#ifdef BACKUP_REPLAY
                        , vc_m0
#endif
                        );

        // reset u
        vc_m0 = VCType_TP(gfs_tp,0.0);
        grid_adaptor.replayData(grid,gfs_tp,projection,vc_m0,transferMap1);

        // clean up
        grid.postAdapt();
      }

      std::stringstream vtu_m0_mirror;
      vtu_m0_mirror  << dir.vtudir << "/mirror_m0_step0";
      VTKPlot::output2vtu( gfs_tp,
                           vc_m0,
                           vtu_m0_mirror.str(),
                           "m0_mirror",
                           pMAX-1 );

#if defined L2ProjectionOfM0
      //typedef Dune::PDELab::Q22DLocalFiniteElementMap<CTYPE,REAL> FEMCG;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<CTYPE,REAL,dim> FEMCG;
      FEMCG femcg;
      typedef Dune::PDELab::ISTLVectorBackend<> VBE_CG;

      typedef Dune::PDELab::OverlappingConformingDirichletConstraints CGCONSTRAINTS;
      CGCONSTRAINTS cgcons;

      typedef Dune::PDELab::GridFunctionSpace<GV_TP,FEMCG,CGCONSTRAINTS,VBE_CG> GFS_CG;
      GFS_CG gfs_cg( gv_tp, femcg, cgcons );
      std::cout << "gfs_cg.maxLocalSize() = " << gfs_cg.maxLocalSize() << std::endl;
      std::cout << "gfs_cg.size() = " << gfs_cg.size() << std::endl;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;

      REAL l2_diffusion = inputdata.transport_parameters.l2_diffusion;
      if( setupdata.solute_concentration_inversion_data.bfixedwidth ) {
        l2_diffusion = setupdata.solute_concentration_inversion_data.regularization_factor / 8.0;
      }


      L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);
      VCType_CG vc_m0_cg( gfs_cg, 0.0 );

      L2sp.apply( vc_m0, vc_m0_cg, l2_diffusion );
      General::logMinAndMax( vc_m0_cg, gv_tp.comm() );
      VTKPlot::output2vtu( gfs_cg,
                           vc_m0_cg,
                           vtu_m0_mirror.str()+"_cg",
                           "vc_m0_cg",
                           -1 );

#endif

#endif



      /*
      // Evaluation loop:
      typedef Dune::PDELab::DiscreteGridFunction<GFS_TP,VCType_TP> DGF_TP;
      typedef Dune::PDELab::GridFunctionToFunctionAdapter<DGF_TP> FUNC;

      // Reference solution on the finest grid level:
      GFS_TP gfs1_tp( *(pool_gv_tp[ maxsteps ]), // pool_gv_tp[ maxsteps ], //
      fem_hyper,
      con_tp );
      VCType_TP vc1_m0( gfs1_tp, 0.0 );
      vc1_m0.std_copy_from( pool_std_vc_m0[ maxsteps ] );
      DGF_TP udgf1( gfs1_tp, vc1_m0 );

      std::cout << "GV level = " << maxsteps << std::endl;
      std::cout << "GV level = " << grid.maxLevel() << std::endl;
      std::cout << "GV size  = " << pool_gv_tp[ maxsteps ]->size(0) << std::endl;
      std::cout << "vc1 basesize = " << vc1_m0.base().size() << std::endl;
      std::cout << "vc1 flatsize = " << vc1_m0.flatsize() << std::endl;
      std::cout << "std size     = " << pool_std_vc_m0[ maxsteps ].size() << std::endl;


      // Loop over grid hierarchy:
      for( int i=0; i<maxsteps+1; ++i ){

      GFS_TP gfs2_tp( *(pool_gv_tp[ i ]), // pool_gv_tp[ i ], //
      fem_hyper,
      con_tp );
      VCType_TP vc2_m0( gfs2_tp, 0.0 );

      vc2_m0.std_copy_from( pool_std_vc_m0[ i ] );

      std::cout << "GV level = " << i << std::endl;
      std::cout << "GV size  = " << pool_gv_tp[i]->size(0) << std::endl;
      std::cout << "vc basesize = " << vc2_m0.base().size() << std::endl;
      std::cout << "vc flatsize = " << vc2_m0.flatsize() << std::endl;
      std::cout << "std size    = " << pool_std_vc_m0[i].size() << std::endl;

      DGF_TP udgf2( gfs2_tp, vc2_m0 );
      FUNC uf2( udgf2 );

      REAL l2errorsquared(0.0);
      std::cout << "Calculate l2-error on grid level " << i << std::endl;
      Dune::PDELab::L2error_on_different_grid_levels( udgf1, uf2, l2errorsquared, 2 );
      std::cout << "l2error = " << sqrt(l2errorsquared) << std::endl;

      }

      */

      std::cout << "Summary:" << std::endl;
      std::cout << "           N"
        /*
          << "    IT"
          << "     T"
        */
                << "     min"
                << "     max"

                << "        l2"
                << "    l2rate"
                << "   eff(l2)"
        /*
          << "    h1semi"
          << "  h1s-rate"
          << "   eff(h1)"

          << " error(dg)"
          << "  rate(dg)"
          << "   eff(dg)"
        */
                << " estimator"
                << std::endl;
      for(std::size_t i=0; i<Ndofs.size(); i++) {

#ifdef STORE_SOLUTION_HIERARCHY
        REAL rate1=0.0;
        if (i>0)
          rate1 = log( l2e_hierarchy[i] / l2e_hierarchy[i-1] ) / log(0.5);
#endif

        /*
          REAL rate2=0.0;
          if (i>0)
          rate2 = log( h1s[i]/h1s[i-1] ) / log(0.5);

          REAL rate3=0.0;
          if (i>0)
          rate3 = log( dge[i]/dge[i-1] ) / log(0.5);
        */
        std::cout << std::setw(3) << i
                  << std::setw(9) << Ndofs[i]
          /*
            << std::setw(6) << nIterations[i]
            << std::setw(6) << std::setprecision(2) << std::fixed << ls_elapsed[i]
          */
                  << std::setw(8) << std::setprecision(2) << std::fixed << minima[i]
                  << std::setw(8) << std::setprecision(2) << std::fixed << maxima[i]

#ifdef STORE_SOLUTION_HIERARCHY
                  << std::setw(10) << std::setprecision(2) << std::scientific << l2e_hierarchy[i]
                  << std::setw(10) << std::setprecision(2) << std::scientific << rate1
                  << std::setw(10) << std::setprecision(2) << std::scientific << ee[i]/(l2e_hierarchy[i])
#endif

          /*
            << std::setw(10) << std::setprecision(2) << std::scientific << h1s[i]
            << std::setw(10) << std::setprecision(2) << std::scientific << rate2
            << std::setw(10) << std::setprecision(2) << std::scientific << ee[i]/(h1s[i])

            << std::setw(10) << std::setprecision(2) << std::scientific << dge[i]
            << std::setw(10) << std::setprecision(2) << std::scientific << rate3
            << std::setw(10) << std::setprecision(2) << std::scientific << ee[i]/(dge[i])
          */
                  << std::setw(10) << std::setprecision(2) << std::scientific << ee[i]
                  << std::endl;
      }

      return;

    }

  } // Gesis

} // Dune

#endif // DUNE_GESIS_ADAPTIVE_DRIVE_HH
