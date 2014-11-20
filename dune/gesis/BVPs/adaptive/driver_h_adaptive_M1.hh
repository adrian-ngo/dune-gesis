#ifndef DUNE_GESIS_ADAPTIVE_M1_HH
#define DUNE_GESIS_ADAPTIVE_M1_HH

#include <dune/gesis/BVPs/DarcyVelocityCache.hh>
#include <dune/gesis/BVPs/MeshPlecletNumber.hh>
#include <dune/gesis/BVPs/sourceterms/source_head_transport.hh>
#include <dune/gesis/BVPs/GroundWaterParameters.hh>
#include <dune/gesis/BVPs/DG/GWE_CCFV.hh>
#include <dune/gesis/BVPs/TransportParameters.hh>
#include <dune/gesis/BVPs/DG/TransportOperatorDG.hh>
#include <dune/gesis/BVPs/DG/TPE_DG.hh>

#include "ee_DG_h_adaptive.hh"  // residual error estimator
#include "adaptivity.hh"

//#include "../integrateFunction.hh"
//#include "../sortedplot.hh"

#include <dune/gesis/BVPs/totalMass.hh>

namespace Dune {
  namespace Gesis {

    template<
      typename GRID,
      typename GFS_GW,
      typename IDT,
      typename YFG,
      typename MEASLIST,
      typename DIR
      >
    void hAdaptiveLoopForM1( GRID& grid,
                             GFS_GW& gfs_gw,  // <-- changed after refinement and load-balancing!
                             const IDT& inputdata,
                             YFG& yfg, // <-- changed after refinement and load-balancing!
                             MEASLIST& orig_measurements,
                             const UINT baselevel,
                             const DIR& dir,
                             const Dune::MPIHelper& helper
                             ){

      typedef typename IDT::SDT SDT;
      int nSetups = inputdata.setups.size();
      for(int iSetup=0; iSetup<nSetups; iSetup++){

      const SDT& setupdata = inputdata.setups[iSetup]; // We consider only the first setup from the input XML.

      Dune::Timer watch;

      if( helper.rank() == 0 )
        std::cout << "Get Gridview for the GWE ..." << std::endl;

      typedef typename GFS_GW::Traits::GridViewType GV_GW;
      enum{ dim = GV_GW::dimension };
      GV_GW gv_gw = gfs_gw.gridView();

      // Setup and solve GWE:
      typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD;
      GWP_FWD gwp_fwd( inputdata, setupdata, yfg );

      typedef PumpingSource<GV_GW,REAL,SDT> PumpingSourceTypeGW;
      PumpingSourceTypeGW source_h( setupdata );

      typedef GroundWaterEquation<GFS_GW,GWP_FWD,PumpingSourceTypeGW,IDT,SDT> GWE_H;  

      watch.reset();
  
      GWE_H gwe_h( gfs_gw,
                   inputdata,
                   setupdata,
                   gwp_fwd,
                   source_h );

        General::log_elapsed_time( watch.elapsed(),
                                   gv_gw.comm(),
                                   inputdata.verbosity,
                                   "GWE",
                                   "Matrix Pattern + Assembly" );

      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;
      VCType_GW vc_head( gfs_gw, 0.0 );
      gwe_h.solve_forward( vc_head );

      if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
        std::cout << gwe_h.show_ls_result() << std::endl;
      General::logMinAndMax( vc_head, gv_gw.comm() );

      if( inputdata.plot_options.vtk_plot_head ){
        std::stringstream vtu_head;
        vtu_head << dir.vtudir << "/h_orig_s" << iSetup;
        VTKPlot::output2vtu( gfs_gw,
                             vc_head,
                             vtu_head.str(),
                             "h_orig",
                             inputdata.verbosity,
                             true,
                             0 );
      }

      //Take the head measurements!
      if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION )
        std::cout << "Take new measurements of head." << std::endl;
      orig_measurements.take_measurements( 1, vc_head, gfs_gw, helper, iSetup);

      std::stringstream bufferfile_head_measurements;
      bufferfile_head_measurements << dir.bufferdimdir 
                                   << "/head_orig.meas.iSetup" << iSetup 
                                   << ".rank" << helper.rank();

      orig_measurements.store_measurements( 1, iSetup, bufferfile_head_measurements.str() );

      if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
        std::cout << "Head measurements stored to " 
                  << bufferfile_head_measurements.str() << " etc." << std::endl;
        // copy file to DATA:
        std::stringstream datafile_per_Setup;
        datafile_per_Setup << dir.datadimdir 
                           << "/head_orig.meas.iSetup" << iSetup;

        if( inputdata.problem_types.generate_measurement_data ){
          General::systemCall( "cp -v " + bufferfile_head_measurements.str() + " " + datafile_per_Setup.str() );
        }
      }

      //typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_DGF;
      typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_BASE;
      typedef DarcyVelocityCache<GWP_FWD,GFS_GW,DARCY_FLUX_BASE> DARCY_FLUX_DGF;

      DARCY_FLUX_DGF darcyflux_dgf( gwp_fwd,
                                    gfs_gw,
                                    vc_head,
                                    baselevel );

      VTKPlot::output_dgf_to_vtu( gv_gw,
                                  gfs_gw,
                                  darcyflux_dgf.exportDGF(),
                                  dir.vtudir + "/q_orig",
                                  "q_orig",
                                  inputdata.verbosity,
                                  true,
                                  0 );


#ifdef NEW_GV

#ifdef USE_DGF_PressureField
      if( helper.rank() == 0 )
        std::cout << "Using DGF_PressureField ... " << std::endl;
      typedef DGF_PressureField<GFS_GW,VCType_GW> RT0_PF;
      Dune::shared_ptr<RT0_PF> rt0_pressurefield
        = Dune::make_shared<RT0_PF>(gfs_gw,vc_head,baselevel);
#else
      if( helper.rank() == 0 )
        std::cout << "Using RT0_PressureField ... " << std::endl;
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
      // This one must be used for parallel ALUGRID. 
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
      std::vector<double> vminima;
      std::vector<double> vmaxima;

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

        if( gv_gw.comm().rank()==0  && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
          std::cout << "grid-size = " << gv_tp.size(0) << std::endl;
        }

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


        UINT allDOFs = gfs_tp.globalSize();
        allDOFs = gv_gw.comm().sum( allDOFs );
        if( gv_gw.comm().rank()==0 ){
          std::cout << "Step " << step
                    << " All m0 DOFs = " << allDOFs
                    << std::endl;
        }

        typedef TransportEquation
          < GFS_TP
            , DARCY_FLUX_DGF
            , SourceTypeTP
            , TPM0
            , IDT
            , SDT
            > TPE_M0;
  
        watch.reset();
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
                                 inputdata.verbosity,
                                 true,
                                 std::max(0,pMAX-1)  );
          }

        }
#endif

        General::log_elapsed_time( watch.elapsed(),
                                   gv_gw.comm(),
                                   inputdata.verbosity,
                                   "TPE",
                                   "Matrix Pattern + Assembly" );
  
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


        if( gv_gw.comm().rank()==0 ){
          std::cout << tpe_m0.show_ls_result()
                    << std::endl;
        }
        REAL u_minimum = 1e100;
        REAL u_maximum = -1e100;
        std::cout << std::setprecision(2) << std::fixed;
        General::logMinAndMax( vc_m0, gv_gw.comm(), u_minimum, u_maximum );
        minima.push_back( u_minimum );
        maxima.push_back( u_maximum );


        //std::vector<REAL> std_vc_m0;
        //General::copy_to_std( vc_m0, std_vc_m0 );
        //pool_std_vc_m0.push_back( std_vc_m0 );

        if( inputdata.plot_options.vtk_plot_m0 ){
          std::stringstream vtu_m0_file;
          vtu_m0_file  << dir.vtudir << "/m0_s" << iSetup << "_step" << step;
          VTKPlot::output2vtu( gfs_tp,
                               vc_m0,
                               vtu_m0_file.str(),
                               "m0_orig",
                               inputdata.verbosity,
                               true,
                               std::max(0,pMAX-1)  );

        }


        REAL m0dg_negMass(0), m0dg_posMass(0), m0dg_totMass(0);
        GridFunctionTools::totalMass( gfs_tp, vc_m0, m0dg_negMass, m0dg_posMass, m0dg_totMass );
        if( helper.rank() == 0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
          std::cout << "pos./neg. mass (m0) = " 
                    << std::fixed << std::setprecision(5)
                    << m0dg_posMass << " " << m0dg_negMass  << std::endl;
          std::cout << "total mass (m0) = " 
                    << std::fixed << std::setprecision(5)
                    << m0dg_totMass << std::endl;
        }


        //Take the measurements!
        if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION )
          std::cout << "Take new measurements of m0." << std::endl;
        orig_measurements.take_measurements( 2, vc_m0, gfs_tp, helper, iSetup);

        std::stringstream bufferfile_M0_measurements;
        bufferfile_M0_measurements << dir.bufferdimdir 
                                   << "/m0_orig.meas"
                                   << ".iSetup" << iSetup 
                                   << ".step" << step
                                   << ".rank" << helper.rank();

        orig_measurements.store_measurements( 2, iSetup, bufferfile_M0_measurements.str() );


        if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
          std::cout << "m0 measurements stored to " 
                    << bufferfile_M0_measurements.str() << " etc." << std::endl;
          // copy file to DATA:
          std::stringstream datafile_per_Setup;
          datafile_per_Setup << dir.datadimdir 
                             << "/m0_orig.meas"
                             << ".iSetup" << iSetup
                             << ".step" << step;
          if( inputdata.problem_types.generate_measurement_data ){
            General::systemCall( "cp -v " + bufferfile_M0_measurements.str() + " " + datafile_per_Setup.str() );
          }
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



        P0GFS p0gfs( gv_tp, p0fem );


        //#ifdef TEST_M1

        typedef FunctionSource<GWP_FWD, 
                               GFS_TP,
                               IDT,SDT> FunctionSourceType;

        FunctionSourceType source_m1_m0( gwp_fwd, 
                                         gfs_tp,
                                         inputdata,
                                         setupdata );

        typedef TransportProblemM1<GV_TP,REAL,DARCY_FLUX_DGF,IDT,SDT> TPM1;
        TPM1 tpm1(darcyflux_dgf,inputdata,setupdata);

        if( gv_gw.comm().rank()==0  && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
          std::cout << "Step " << step
                    << " All m1 DOFs = " << allDOFs
                    << std::endl;
        }

        typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,FunctionSourceType,TPM1,IDT,SDT> TPE_M1_M0;

        watch.reset();
        TPE_M1_M0 tpe_m1_m0( gfs_tp,
                             inputdata,
                             setupdata,
                             darcyflux_dgf,
                             tpm1,
                             source_m1_m0 );
              
        General::log_elapsed_time( watch.elapsed(),
                                   gv_gw.comm(),
                                   inputdata.verbosity,
                                   "TPE",
                                   "Matrix Pattern + Assembly" );


        // Define sourceterm for the forward transport equation for m1:
        VCType_TP vc_theta_m0 = vc_m0; // theta_m0 := porosity * m0
        // multiply a vector container with the porosity. needed for zonation!

        //General::vc_times_theta(vc_theta_m0,*pgv_tp,inputdata);
        General::vc_times_theta(vc_theta_m0,gv_tp,inputdata);

        tpe_m1_m0.set_rhs( vc_theta_m0 );


        VCType_TP vc_m1_m0( gfs_tp, 0.0 );
        tpe_m1_m0.solve_forward( vc_m1_m0 );


        if( gv_gw.comm().rank()==0  && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
          std::cout << tpe_m1_m0.show_ls_result()
                    << std::endl;
        }
        REAL v_minimum = 1e100;
        REAL v_maximum = -1e100;
        std::cout << std::setprecision(2) << std::scientific;
        General::logMinAndMax( vc_m1_m0, gv_gw.comm(), v_minimum, v_maximum );
        vminima.push_back( v_minimum );
        vmaxima.push_back( v_maximum );





        if( inputdata.plot_options.vtk_plot_m1 ){
          std::stringstream vtu_m1m0_file;
          vtu_m1m0_file  << dir.vtudir << "/m1m0_s" << iSetup << "_step" << step;
          VTKPlot::output2vtu( gfs_tp,
                               vc_m1_m0,
                               vtu_m1m0_file.str(),
                               "m1_orig",
                               inputdata.verbosity,
                               true,
                               std::max(0,pMAX-1)  );
        }

        REAL m1dg_negMass(0), m1dg_posMass(0), m1dg_totMass(0);
        GridFunctionTools::totalMass( gfs_tp, vc_m1_m0, m1dg_negMass, m1dg_posMass, m1dg_totMass );
        if( helper.rank() == 0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
          std::cout << "pos./neg. mass (m1) = " 
                    << std::scientific << std::setprecision(5)
                    << m1dg_posMass << " " << m1dg_negMass  << std::endl;
          std::cout << "total mass (m1) = " 
                    << std::scientific << std::setprecision(5)
                    << m1dg_totMass << std::endl;
        }

        //Take the measurements!
        if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION )
          std::cout << "Take new measurements of m1." << std::endl;
        orig_measurements.take_measurements( 3, vc_m1_m0, gfs_tp, helper, iSetup);



        std::stringstream bufferfile_M1_measurements;
        bufferfile_M1_measurements << dir.bufferdimdir 
                                   << "/m1_orig.meas"
                                   << ".iSetup" << iSetup 
                                   << ".step" << step
                                   << ".rank" << helper.rank();

        if( helper.rank() == 0 )
          std::cout << "Absolute m1 error for Setup #" << iSetup << " is "
                    << orig_measurements.calibrateAbsErrorForM1( iSetup )
                    << std::endl;

        orig_measurements.store_measurements( 3, iSetup, bufferfile_M1_measurements.str() );


        if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
          std::cout << "m1 measurements stored to " 
                    << bufferfile_M1_measurements.str() << " etc." << std::endl;
          // copy file to DATA:
          std::stringstream datafile_per_Setup;
          datafile_per_Setup << dir.datadimdir 
                             << "/m1_orig.meas"
                             << ".iSetup" << iSetup
                             << ".step" << step;
          if( inputdata.problem_types.generate_measurement_data ){
            General::systemCall( "cp -v " + bufferfile_M1_measurements.str() + " " + datafile_per_Setup.str() );
          }
        }

        /*
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
        //vc_MeanArrivalTime.std_copy_from( vmat );
        General::copy_from_std( vmat, vc_MeanArrivalTime );

        std::stringstream vtu_mat_file;
        vtu_mat_file  << dir.vtudir << "/mat_step" << step;
        VTKPlot::output2vtu( gfs_tp, 
                             vc_MeanArrivalTime, 
                             vtu_mat_file.str(), 
                             "MeanArrivalTime", 
                             inputdata.verbosity, 
                             true, 
                             std::max(0,pMAX-1) 
                             );
        */
        //#endif // TEST_M1
  
        if( maxsteps < 1 )
          break;

        //typedef DGResidualEstimator<IDT,SDT,TPM0> ESTLOP;
        //ESTLOP estlop(inputdata,setupdata,tpm0);
        typedef DGResidualEstimator<IDT,SDT,TPM1> ESTLOP;
        ESTLOP estlop(inputdata,setupdata,tpm1);

#define Compute_Error_Estimate

#ifdef Compute_Error_Estimate

        if( helper.rank() == 0 ){
          std::cout << std::endl;
          std::cout << "Calculate Error Estimate..." << std::endl;
        }
        

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

        //estgo.residual( vc_m0, x0Power );   // <-- Run over grid cells and evaluate!
        estgo.residual( vc_m1_m0, x0Power );   // <-- Run over grid cells and evaluate!

        //typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,REAL>::Type U0Type;
        U0Type eta(p0gfs,0.0); // local error indicator
        U0Type rho(p0gfs,0.0); // some other local indicator
        General::extract_two_indicators( p0powergfs, x0Power, p0gfs, eta, rho );
        //General::extract_two_indicators( p0powergfs, x0Power, p0gfs, eta, eta );

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

        std::cout << std::setprecision(2) << std::scientific << std::endl;
        std::cout << ">>> Estimated Error = " 
                  << estimated_error
                  << std::endl
                  << std::endl;


        if(true){
      
          std::vector<REAL> eta_flat;
          General::copy_to_std( eta, eta_flat );

          std::stringstream vtu_m1_eta_file;
          vtu_m1_eta_file  << dir.vtudir << "/eta_M1_step" << step;

          VTKPlot::output_vector_to_vtu( gv_tp,
                                         eta_flat,
                                         vtu_m1_eta_file.str(), 
                                         "eta"
                                         );
          /*
            VTKPlot::output2vtu( p0gfs, 
            eta, 
            vtu_m0_eta_file.str(), 
            "eta", 
            inputdata.verbosity
            //, true
            //, 0
            );
          */
        }

#endif // Compute_Error_Estimate


        if( tpe_m0.isAcceptableUpDown( vc_m0, 100.0, 0.0, inputdata.transport_parameters.tolerance ) )
          break; // Exit the adaptation loop when tolerance level is reached.


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
          //grid_adaptor.backupData(grid,gfs_tp,projection,vc_m0,transferMap1);
          grid_adaptor.backupData(grid,gfs_tp,projection,vc_m1_m0,transferMap1);

          // save intermediate solution hierarchy

#ifdef STORE_SOLUTION_HIERARCHY
          //General::copy_to_std( vc_m0, solution_hierarchy[step] );
          General::copy_to_std( vc_m1_m0, solution_hierarchy[step] );

          std::vector<typename GA::MapType> solutionMap(step+1);
          for(int ii=0;ii<step+1;++ii){
            VCType_TP tmpSolution(gfs_tp,0.0);

            //std::cout << "before backup: solution_hierarchy[" << ii << "].size() = " 
            //          << solution_hierarchy[ii].size()
            //          << std::endl;

            //tmpSolution.std_copy_from(solution_hierarchy[ii]);
            General::copy_from_std( solution_hierarchy[ii], tmpSolution );

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
                          //, vc_m0
                          , vc_m1_m0
#endif
                          );


#ifdef STORE_SOLUTION_HIERARCHY
          // reset u
          //vc_m0 = VCType_TP(gfs_tp,0.0);
          //grid_adaptor.replayData(grid,gfs_tp,projection,vc_m0,transferMap1);
          vc_m1_m0 = VCType_TP(gfs_tp,0.0);
          grid_adaptor.replayData(grid,gfs_tp,projection,vc_m1_m0,transferMap1);

      
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
          std::cout << "STOP: max. steps or max. DOFs reached at step " << maxsteps 
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
      //DGF_TP udgf( gfs_tp, vc_m0 );
      DGF_TP udgf( gfs_tp, vc_m1_m0 );

      // Compare solution hierarchy against reference solution
      std::vector<Dune::shared_ptr<VCType_TP> > pSol(maxsteps);
      for( int ii=0; ii<maxsteps; ++ii ){
        pSol[ii] = Dune::make_shared<VCType_TP>( gfs_tp, 0.0 );

        //(pSol[ii])->std_copy_from( solution_hierarchy[ii] );
        General::copy_from_std( solution_hierarchy[ii], *(pSol[ii]) );

        DGF_TP udgf_ii( gfs_tp, *(pSol[ii]) );
        DifferenceSquaredAdapter<DGF_TP,DGF_TP> differencesquared(udgf,udgf_ii);
        Dune::FieldVector<REAL,1> l2errorsquared(0.0);
        std::cout << "Calculate l2-error on grid level " << ii << std::endl;
        Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,2*pMAX);
        std::cout << "l2error = " << sqrt(l2errorsquared) << std::endl;
        l2e_hierarchy.push_back(sqrt(l2errorsquared));


        // debug plots:
        if( inputdata.plot_options.vtk_plot_m1 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
          std::stringstream vtu_m1_file;
          vtu_m1_file  << dir.vtudir << "/pSol_m1_step" << ii;
          VTKPlot::output2vtu( gfs_tp,
                               *(pSol[ii]),
                               vtu_m1_file.str(),
                               "m1_orig",
                               inputdata.verbosity,
                               true,
                               std::max(0,pMAX-1)  );
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
        Dune::Gesis::mark_all_coarsen( grid );

        // prepare the grid for refinement
        grid.preAdapt();
    
        // save u
        typename GA::MapType transferMap1;
        //grid_adaptor.backupData(grid,gfs_tp,projection,vc_m0,transferMap1);
        grid_adaptor.backupData(grid,gfs_tp,projection,vc_m1_m0,transferMap1);

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
                        //, vc_m0
                        , vc_m1_m0
#endif
                        );

        // reset u
        //vc_m0 = VCType_TP(gfs_tp,0.0);
        vc_m1_m0 = VCType_TP(gfs_tp,0.0);
        //grid_adaptor.replayData(grid,gfs_tp,projection,vc_m0,transferMap1);
        grid_adaptor.replayData(grid,gfs_tp,projection,vc_m1_m0,transferMap1);

        // clean up
        grid.postAdapt();
      }

      std::stringstream vtu_m1_mirror;
      vtu_m1_mirror  << dir.vtudir << "/mirror_m1_step0";
      VTKPlot::output2vtu( gfs_tp,
                           vc_m1_m0,
                           vtu_m1_mirror.str(),
                           "m1_mirror",
                           inputdata.verbosity,
                           true,
                           std::max(0,pMAX-1)  );

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
      //VCType_CG vc_m0_cg( gfs_cg, 0.0 );
      VCType_CG vc_m1_cg( gfs_cg, 0.0 );

      //L2sp.apply( vc_m0, vc_m0_cg, l2_diffusion );
      L2sp.apply( vc_m1_m0, vc_m1_cg, l2_diffusion );

      General::logMinAndMax( vc_m1_cg, gv_tp.comm() );
      VTKPlot::output2vtu( gfs_cg,
                           vc_m1_cg,
                           vtu_m1_mirror.str()+"_cg",
                           "vc_m1_cg",
                           inputdata.verbosity,
                           true,
                           0 );

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


      for( int ibackstep=grid.maxLevel(); ibackstep>0; --ibackstep ){
        Dune::Gesis::mark_all_coarsen( grid );
        grid.adapt();
        gfs_tp.update();
      }

      std::cout << "Summary:" << std::endl;
      std::cout << "           N"
        /*
          << "    IT" 
          << "     T" 
        */
                << "     min" 
                << "     max" 

                << "   est.err" 
                << "      rate" 
        /*
                << "   eff(l2)"
          << "    h1semi" 
          << "  h1s-rate"
          << "   eff(h1)"

          << " error(dg)" 
          << "  rate(dg)" 
          << "   eff(dg)"
                << " estimator"
        */
                << std::endl;
      for(std::size_t i=0; i<Ndofs.size(); i++) {
    
        REAL l2_error_rate=0.0;
        if (i>0)
          l2_error_rate = log( ee[i]/ee[i-1] ) / log(0.5); 

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
                  << std::setw(10) << std::setprecision(2) << std::scientific << l2_error_rate
                  << std::endl;
      }

      } // iSetup loop

      return;

    }

  } // Gesis

} // Dune

#endif // DUNE_GESIS_ADAPTIVE_M1_HH
