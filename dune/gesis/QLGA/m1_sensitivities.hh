/*
 * This function is to calculate the sensitivity of a m1 measurement wrt to lnK.
 *   It calculates also the sensitivity times the covariance matrix (for the current measurement)
 *   and the output.
 * 
 *   "Sensitivity" and "sensitivity times the covariance matrix" are stored on the disk via HDF5
 * 
 * 
 * INPUTS:
 * gv:            grid view
 * gfs:           grid function space ?
 * inputdata:     the class holding all input information
 * dir:           information about IO Paths and file locations!
 * nAllCells:     number of all cells
 * orig_measurements: class holding information about the measuring points
 * Lambdas:       eigenvalues of the extended covariance matrix
 * Y_old:         previous Y field
 * YfieldGenerator_old: needed for the YFieldGenerator wrapper class
 * qorder:        needed for the GWE_ADJ, ...
 * iteration_number: the current iteration number (only for some VTK output needed)
 * helper:        Dune MPI helper
 * CommunicatorPool : vector containing all information about the different communicators! 
 * 
 * 
 * OUTPUTS
 * JX:            value for the multiplication of the sensitivity times the trend matrix (at the current meas location)
 * J_times_Y_old: value for the multiplication of the sensitivity times old Y field (at the current meas location)
 */

#ifndef DUNE_GESIS_M1_SENSITIVITIES_HH
#define DUNE_GESIS_M1_SENSITIVITIES_HH
#include "../common/MyMPIComm.hh"

#include "m1_sensitivity_field.hh"

namespace Dune{
  namespace Gesis{

    // nCellsExt : Dimensions of extended field per zone
    // Lambdas: eigenvalues of the covariance matrix per zone
    // X: zonation matrix -> only on root process 0 (leader of the communicator group) 

    template<typename POOL_GRID,
             typename GFS_GW,
             typename GFS_TP,
             typename GFS_CG,
             typename VCType_GW,
             typename VCType_CG,
             typename MEASLIST,
             typename YFG,
             typename DIR,
             typename IDT,
             typename SDT
             >
    void m1_sensitivities( // input:
                          POOL_GRID& pool_grid,
                          const GFS_GW& gfs_gw,
                          const GFS_TP& gfs_tp,
                          const GFS_CG& gfs_cg,

                          const IDT& inputdata,
                          const SDT& setupdata,
                          const DIR& dir,
                          
                          const MEASLIST& orig_measurements,
                          const MEASLIST& measurements,
                          
                          const Vector<REAL>& Lambdas,
                          const Vector<REAL>& Xzones,
                          
                          const Vector<REAL>& Y_old,// TODO: Check whether this can be exported from YfieldGenerator_old
                          const YFG& YfieldGenerator_old,

                          const VCType_GW& vchead_old,
                          const VCType_CG& vcM0_old_cg,
                          const VCType_CG& vcM1_old_cg,

                          UINT iteration_number,
                          const Dune::MPIHelper& helper,
                          std::vector<MyMPIComm> CommunicatorPool,
                          // output:
                          Vector<REAL>& JX,
                          Vector<REAL>& J_times_Y_old
                           )
    {

      typedef Dune::Gesis::FieldData FD;
      FD fielddata(inputdata);

      if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION )
#ifdef USE_ALL_ADJOINTS
        std::cout << "m1 sensitivities using full set of adjoint equations." << std::endl;
#else
        std::cout << std::endl;
        std::cout << "Note: m1 sensitivities NOT using full set of adjoint equations!!!" << std::endl;
#endif

      logger << "m1_sensitivities:" << std::endl;

#ifdef USE_YASP
      int baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_UG
      int baselevel = inputdata.domain_data.ug_baselevel;
#endif
#ifdef USE_ALUGRID
      int baselevel = inputdata.domain_data.alugrid_baselevel;
#endif

      // Retrieve types and data from parameters:
      typedef typename GFS_GW::Traits::GridViewType GV_GW;
      enum{dim=GV_GW::dimension};
      typedef typename GFS_TP::Traits::GridViewType GV_TP;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;

      GV_GW gv_0 = pool_grid.levelGridView(0);
      GV_GW gv_gw = gfs_gw.gridView();
      GV_TP gv_tp = gfs_tp.gridView(); // re-ordered by vchead_old

      // Total number of all cells required to resolve the parameter field
      const UINT nAllCells = inputdata.domain_data.nCells.componentsproduct();
      

      UINT iSetup = setupdata.index;

 
      logger<<"m1 SENS: precalculations ... "<<std::endl;
      Dune::Timer watch;

  
      // number of m1 measurements!
      UINT nPoints_m1=orig_measurements.nMeasPerSetupPerType(iSetup,3);
  
      //number of available communicators. BE CAREFUL "nComm <= inputdata.maxNComm", because inputdata holds the given input for the pool size and this is a max. value!!!! 
      UINT nComm=CommunicatorPool.size();
  
      // some more needed variables
      Dune::FieldVector<REAL,dim> measure_location(0.0);
      Vector<UINT> read_local_offset;
      Vector<UINT> read_local_count;

      typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD_OLD;
      GWP_FWD_OLD gwp_fwd_old( inputdata, setupdata, YfieldGenerator_old );
      
      typedef GroundwaterAdjointProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_ADJ_OLD;
      GWP_ADJ_OLD gwp_adj_old( inputdata, setupdata, YfieldGenerator_old );
      
      //typedef GradientVectorField<GWP_FWD_OLD,GFS_GW> DARCY_FLUX_DGF;
      typedef GradientVectorField<GWP_FWD_OLD,GFS_GW> DARCY_FLUX_BASE;
      typedef DarcyVelocityCache<GWP_FWD_OLD,GFS_GW,DARCY_FLUX_BASE> DARCY_FLUX_DGF;

      // darcy_flux dgf of the old head
      DARCY_FLUX_DGF darcyflux_dgf( gwp_fwd_old,
                                    gfs_gw,
                                    vchead_old,
                                    baselevel );

      L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);

      /*
       * define m1_adj which is actually identical to m0_adj in "m0_sensitivities.hh"
       */

      // New: Point source with adjustable smearing!
      typedef PointFunctionSource<GFS_TP,IDT> PointSourceTypeTP;
      PointSourceTypeTP source_m1_adj( gfs_tp, 
                                       inputdata, 
                                       -1.0, // negative pointsource
                                       setupdata.solute_concentration_inversion_data.regularization_factor,
                                       setupdata.solute_concentration_inversion_data.bfixedwidth );
      // TODO: How about heat_inversion_data for m1? well_flag?


      typedef Dune::Gesis::TransportProblemAdjoint<GV_TP,
                                                          REAL,
                                                          DARCY_FLUX_DGF,
                                                          IDT,
                                                          SDT> TP_ADJ;

      TP_ADJ tpm1_adj(darcyflux_dgf,inputdata,setupdata);

      typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,PointSourceTypeTP,TP_ADJ,IDT,SDT> 
        TPE_M1_ADJ;

      // Setup stiffness matrix for adjoint TPE, only once for all measuring points:
      const EQ::Mode equationMode = EQ::adjoint;

      watch.reset();
      TPE_M1_ADJ tpe_m1_adj( gfs_tp,
                             inputdata,
                             setupdata,
                             darcyflux_dgf,
                             tpm1_adj,
                             source_m1_adj,
                             Passenger::solute,
                             equationMode );

      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 inputdata.verbosity,
                                 "TPE",
                                 "Matrix Pattern + Assembly, M1 Adjoint" );


      
      /*
       * for m0 adjoint depending on m1 adjoint
       */
#ifdef USE_FEM
      typedef FunctionSource<GWP_ADJ_OLD,GFS_TP,IDT,SDT> FunctionSourceADJ;
      FunctionSourceADJ source_m0adj_m1adj( gwp_adj_old, 
                                            gfs_tp, 
                                            inputdata, 
                                            setupdata, 
                                            0 );
#else
      typedef FunctionSource<GWP_ADJ_OLD,GFS_CG,IDT,SDT> FunctionSourceADJ;
      FunctionSourceADJ source_m0adj_m1adj( gwp_adj_old, 
                                            gfs_cg, 
                                            inputdata, 
                                            setupdata, 
                                            0 );
#endif
      
  
      //if( helper.rank()==0 && inputdata.verbosity > 3 )
      //  std::cout << "Setup TPE_M0ADJ_M1ADJ with function source (m1_sensitivities)" << std::endl;

      typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,FunctionSourceADJ,TP_ADJ,IDT,SDT> TPE_M0ADJ_M1ADJ;

      TP_ADJ tp_m0m1adj(darcyflux_dgf,inputdata,setupdata);

      // matrix assembly for m0 adjoint
      watch.reset();
      TPE_M0ADJ_M1ADJ tpe_m0adj_m1adj( gfs_tp,
                                       inputdata,
                                       setupdata,
                                       darcyflux_dgf,
                                       tp_m0m1adj,
                                       source_m0adj_m1adj,
                                       Passenger::solute,
                                       equationMode );
      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 inputdata.verbosity,
                                 "TPE",
                                 "Matrix Pattern + Assembly, M0M1 Adjoint" );
  
  
      /*
       * for h adjoint
       */
#ifdef USE_FEM
      typedef FunctionSource<GWP_ADJ_OLD,GFS_TP,IDT,SDT>
        FunctionSourceCombiADJ;
      FunctionSourceCombiADJ source_hadj_m0m1( gwp_adj_old,
                                               gfs_tp,
                                               inputdata,
                                               setupdata,
                                               pool_grid.maxLevel()
                                               );
#else
      typedef FunctionSource<GWP_ADJ_OLD,GFS_GW,IDT,SDT>
        FunctionSourceCombiADJ;
      FunctionSourceCombiADJ source_hadj_m0m1( gwp_adj_old, 
                                               gfs_gw, 
                                               inputdata,
                                               setupdata,
                                               baselevel );
#endif


#ifdef USE_ALL_ADJOINTS // GWE
      if( helper.rank()==0 && inputdata.verbosity > VERBOSITY_EQ_SUMMARY )
        std::cout << "Setup GWE_H_ADJ_M0M1 with function source (m1_sensitivities)" << std::endl;

      typedef GroundWaterEquation<GFS_GW,GWP_ADJ_OLD,FunctionSourceCombiADJ,IDT,SDT> 
        GWE_H_ADJ_M0M1;

      // stiffness matrix is assembled here
      watch.reset();
      GWE_H_ADJ_M0M1 gwe_h_adj_m0m1( gfs_gw,
                                   inputdata,
                                   setupdata,
                                   gwp_adj_old,
                                   source_hadj_m0m1,
                                   equationMode,
                                   setupdata.flow_equation.bRecycleMatrixHierarchy );
      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 inputdata.verbosity,
                                 "GWE",
                                 "Matrix Pattern + Assembly, M0M1 Adjoint Head" );

#endif // USE_ALL_ADJOINTS
  
  
      // Loop over all m1 measurements
      for (UINT meas_point=0; meas_point<nPoints_m1;meas_point++){
        watch.reset();
     
        UINT global_meas_id = orig_measurements.get_global(iSetup,3,meas_point);
      
        /*
         * distribute the measurements to the available MPI pools!!!! 
         */
        if(CommunicatorPool[meas_point%nComm].I_am_in()){
          logger<<"working on meas #"<<global_meas_id<<" (in parallel)"<<std::endl;
      

          // actual measurement location
          measure_location[0] = orig_measurements[global_meas_id].x;
          measure_location[1] = orig_measurements[global_meas_id].y;
#ifdef DIMENSION3
          measure_location[2] = orig_measurements[global_meas_id].z;
#endif
          // vector of the sensitivity
          Vector<REAL> Sensitivity(gv_0.size(0),0.0);
          Vector<REAL> iSensitivity;

          // If the corresponding Sensitivity.h5 file exists there is 
          // no need to waste time doing the calculation once again!
          if( !General::bFileExists( dir.Sensitivity_h5file[global_meas_id] ) ){
            /*
             * solve for m1 adjoint
             */
            VCType_TP vc_m1adj( gfs_tp, 0.0 );
            tpe_m1_adj.solve_adjoint( measure_location, vc_m1adj );

            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << tpe_m1_adj.show_ls_result() << std::endl;
            General::logMinAndMax( vc_m1adj, gv_tp.comm() );

            VCType_CG vc_m1adj_cg( gfs_cg, 0.0 );
            L2sp.apply( vc_m1adj, vc_m1adj_cg,
                        inputdata.transport_parameters.l2_diffusion_adjoint );
            General::logMinAndMax( vc_m1adj_cg, gv_tp.comm() );

     
            if( inputdata.plot_options.vtk_plot_adjoint_m1 ){
              std::stringstream vtu_m1_adj;
              vtu_m1_adj << dir.vcM1_adj_prefix
                         << "_s" << iSetup
                         << "_m" << global_meas_id 
                         << "_i" << iteration_number;
              VTKPlot::output2vtu( gfs_cg, 
                                   vc_m1adj_cg, 
                                   vtu_m1_adj.str(), 
                                   "m1_Adjoint", 
                                   -1 );
            }
          
#ifdef USE_ALL_ADJOINTS // m0m1.solveAdjoint
            /*
             * solve for m0 adjoint (given m1 adjoint x porosity)
             */
          
            // Adrian: This is not required!? Ohne geht's deutlich besser!?
            VCType_CG vc_m1adj_theta_cg = vc_m1adj_cg;
            //multiply a vector container with the porosity. needed for zonation!
            Dune::Gesis::General::
              vc_times_theta(vc_m1adj_theta_cg,gv_tp,inputdata);

            VCType_TP vc_m0adj_m1adj( gfs_tp, 0.0 );
            tpe_m0adj_m1adj.solve_adjoint( vc_m1adj_theta_cg, // input
                                           vc_m0adj_m1adj     // output (solution)
                                           );
            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << tpe_m0adj_m1adj.show_ls_result() << std::endl;
            General::logMinAndMax( vc_m0adj_m1adj, gv_tp.comm() );

#endif// USE_ALL_ADJOINTS // m0m1.solveAdjoint

            VCType_CG vc_m0adj_m1adj_cg( gfs_cg, 0.0 );

#ifdef USE_ALL_ADJOINTS // m0m1.solveAdjoint
            L2sp.apply( vc_m0adj_m1adj, vc_m0adj_m1adj_cg );

            General::logMinAndMax( vc_m0adj_m1adj_cg, gv_tp.comm() );

            if( inputdata.plot_options.vtk_plot_adjoint_m1 ){
              std::stringstream vtu_m0adj_m1adj;
              vtu_m0adj_m1adj << dir.vcM0_adj_given_M1_prefix
                              << "_s" << iSetup
                              << "_m" << global_meas_id 
                              << "_i" << iteration_number;

              Dune::Gesis::VTKPlot::output2vtu( gfs_cg, 
                                                vc_m0adj_m1adj_cg,
                                                vtu_m0adj_m1adj.str(), 
                                                "m0adj_m1adj", 
                                                -1 );
            }

#endif// USE_ALL_ADJOINTS // m0m1.solveAdjoint
     
            /*
             * solve for head adjoint (given m0, m1, m0 adjoint, and m1 adjoint)
             */
            VCType_GW vc_hadj_m0m1( gfs_gw, 0.0 );


#ifdef USE_ALL_ADJOINTS // GWE.solveAdjoint

#ifdef USE_FEM
            gwe_h_adj_m0m1.solve_adjoint( vc_m0adj_m1adj_cg
                                          , vcM0_old_cg     // second input
                                          , vc_m1adj_cg     // third input
                                          , vcM1_old_cg      // fourth input
                                          , vc_hadj_m0m1   // output (solution)
                                          );
#else

            CoarseGridP0Projector<GFS_CG,GFS_GW,IDT,dim> sp_ccfv(gfs_cg,gfs_gw,inputdata);

            VCType_GW vc_m0adj_m1adj_gw( gfs_gw, 0.0 );
            sp_ccfv.apply( vc_m0adj_m1adj_cg, 
                           pool_grid.maxLevel(),
                           baselevel, 
                           vc_m0adj_m1adj_gw );
          
            VCType_GW vcM0_old_gw( gfs_gw, 0.0 );
            sp_ccfv.apply( vcM0_old_cg, 
                           pool_grid.maxLevel(),
                           baselevel, 
                           vcM0_old_gw );
          
            VCType_GW vc_m1adj_gw( gfs_gw, 0.0 );
            sp_ccfv.apply( vc_m1adj_cg, 
                           pool_grid.maxLevel(),
                           baselevel, 
                           vc_m1adj_gw );

            VCType_GW vcM1_old_gw( gfs_gw, 0.0 );
            sp_ccfv.apply( vcM1_old_cg, 
                           pool_grid.maxLevel(),
                           baselevel, 
                           vcM1_old_gw );

            gwe_h_adj_m0m1.solve_adjoint( vc_m0adj_m1adj_gw // first input (order is very important here!)
                                          , vcM0_old_gw     // second input
                                          , vc_m1adj_gw     // third input
                                          , vcM1_old_gw      // fourth input
                                          , vc_hadj_m0m1   // output (solution)
                                          );
#endif

            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << gwe_h_adj_m0m1.show_ls_result() << std::endl;

#endif // USE_ALL_ADJOINTS

     
            if( inputdata.plot_options.vtk_plot_adjoint_m1 ){
              std::stringstream vtu_hadj_m0m1;
              vtu_hadj_m0m1 << dir.vchead_adj_given_M0M1_prefix
                            << "_s" << iSetup
                            << "_m" << global_meas_id 
                            << "_i" << iteration_number;
              Dune::Gesis::VTKPlot::output2vtu( gfs_gw, 
                                                vc_hadj_m0m1,
                                                vtu_hadj_m0m1.str(), 
                                                "hadj_m0m1", 
                                                -1 );
            }

     
            logger << "calculate sensitivity ... " << std::endl;

            DARCY_FLUX_BASE grad_hadj_m0m1_dgf( gwp_fwd_old, // dummy place-holder!
                                               gfs_gw, 
                                               vc_hadj_m0m1, 
                                               baselevel, 
                                               true //compute gradient only
                                               );


            get_m1_Sensitivities( gv_0,
                                  gv_gw,
                                  darcyflux_dgf,
                                  grad_hadj_m0m1_dgf,
                                  gfs_cg,
                                  vc_m0adj_m1adj_cg,
                                  vc_m1adj_cg,
                                  vcM0_old_cg,
                                  vcM1_old_cg,
                                  pool_grid.maxLevel(),
                                  baselevel,
                                  2*pMAX,
                                  Sensitivity,
                                  iSensitivity
                                  );


            // writing the sensitivity to HDF5
            // the inversion kernel will read it later

            HDF5Tools::h5g_pWrite( gv_0
                                   , iSensitivity
                                   , dir.Sensitivity_h5file[global_meas_id]
                                   , "/Sensitivity"
                                   , fielddata
                                   , 1
                                   , FEMType::DG
                                   , 0
                                   , false // preserve_structure
                                   );

          }
          else {
            logger << "m0_sensitivities: Read existing Sensitivity field from " 
                   << dir.Sensitivity_h5file[global_meas_id]
                   << std::endl;

            Vector<UINT> local_offset;
            Vector<UINT> local_count;
            HDF5Tools::h5g_pRead( gv_0
                                  , Sensitivity
                                  , dir.Sensitivity_h5file[global_meas_id]
                                  , "/Sensitivity"
                                  , local_offset
                                  , local_count
                                  , fielddata
                                  );
          }

          REAL s1norm = Sensitivity.one_norm();
          s1norm = gv_gw.comm().sum( s1norm );
          if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
            std::cout << General::getDateAndTime()
                      << "====> m" << global_meas_id 
                      << ": " << measure_location 
                      << ", m1 Sensitivity (total) = " << s1norm << std::endl;
          }

          // wait for all processes, only on the current MPI pool
          if(CommunicatorPool[meas_point%nComm].get_size()>1)
            MPI_Barrier(CommunicatorPool[meas_point%nComm].get_comm());
          
          
          //calculate cross_covariance_Parallel!
          logger<<"calculating JQ PARALLEL ...."<<std::endl;
          
          //calculate JQ!
          cross_covariance_JQ( fielddata,
                               Lambdas,
                               dir.Sensitivity_h5file[global_meas_id],
                               dir.JQ_h5file[global_meas_id],
                               CommunicatorPool[meas_point%nComm].get_comm(),
                               global_meas_id,
                               iteration_number,
                               dir
                               );
    
      
          logger<<"working on meas #"<<global_meas_id<<" (parallel part DONE)"<<std::endl;
      
          // sequential part of the function -> Only P0 (in the responsible MPI pool) works now
          //   calculated the needed stuff, e.g., cross covariance
          if(CommunicatorPool[meas_point%nComm].get_rank()==0){ // something to do for P0
            logger<<"working on meas #"<<global_meas_id<<" (sequential part)"<<std::endl;
	
            // if parallel MPI pool then get the sensitivity data from the disk
            if( CommunicatorPool[meas_point%nComm].get_size()>1){
              logger<<"read the sensitivity sequential"<<std::endl; 
              HDF5Tools::h5g_Read( Sensitivity
                                   , dir.Sensitivity_h5file[global_meas_id]
                                   , "/Sensitivity"
                                   , read_local_offset
                                   , read_local_count
                                   );

              if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
                std::cout << "====> m" << global_meas_id 
                          << ": " << measure_location 
                          << ", m1 Sensitivity (seq.read) = " << Sensitivity.one_norm() << std::endl;
              }

            }
        
            for(UINT iCell=0; iCell<nAllCells; iCell++)
              {
                JX[global_meas_id] += Sensitivity[iCell] * Xzones[iCell];//X[iCell];//1.0; //ZonationMatrixX[ iCell ];
                J_times_Y_old[global_meas_id] += Sensitivity[ iCell ] * Y_old[ iCell ];
              }
         
            logger<<"working on meas #"<<global_meas_id<<" (sequential part DONE)"<<std::endl;
          }// END: if(CommunicatorPool[meas_point%nComm].get_rank()==0)   
        }//END: if(CommunicatorPool[meas_point%nComm].I_am_in())

      }//END: for (UINT meas_point=0; meas_point<nPoints_head;meas_point++)

      darcyflux_dgf.emptyCache();

    } // END: m1_sensitivities

  }

}

#endif // DUNE_GESIS_M1_SENSITIVITIES_HH
