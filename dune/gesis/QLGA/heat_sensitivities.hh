/*
 * This function is to calculate the sensitivity of a temperature measurement wrt to lnK.
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
 * nCellsExt:     number of all cells of the extended domain_data
 * orig_measurements: class holding information about the measuring points
 * measurements   information about the current measurements
 * Lambdas:       eigenvalues of the extended covariance matrix
 * Y_old:         previous Y field
 * YfieldGenerator_old: needed for the YFieldGenerator wrapper class
 * qorder:        needed for the GWE_ADJ
 * iSetup:         id for the current setup (important for the boundary conditions)
 * iteration_number: the current iteration number (only for some VTK output needed)
 * helper:        Dune MPI helper
 * CommunicatorPool : vector containing all information about the different communicators! 
 * heat_tomography_setup:  the current heat_tomography setup
 * previous_heat_tomography:   offset value due to previous heat tomography setups
 * 
 * 
 * OUTPUTS
 * JX:            value for the multiplication of the sensitivity times the trend matrix (at the current meas location)
 * J_times_Y_old: value for the multiplication of the sensitivity times old Y field (at the current meas location)
 */

#ifndef DUNE_HEAT_SENSITIVITIES_HH
#define DUNE_HEAT_SENSITIVITIES_HH
#include "../common/MyMPIComm.hh"

#include "m0_sensitivity_field.hh"
#include "m1_sensitivity_field.hh"
#include "AT_sensitivity_field.hh"


namespace Dune {
  namespace GeoInversion {
    
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
    void heat_sensitivities( // input:
                            POOL_GRID& pool_grid,
                            const GFS_GW& gfs_gw,
                            const GFS_TP& gfs_tp,
                            const GFS_CG& gfs_cg,

                            const IDT& inputdata,
                            const SDT& setupdata,
                            const DIR& dir,
                          
                            const MEASLIST& orig_measurements,
                            const MEASLIST& measurements,
                          
                            const std::vector< Vector<UINT> >& nCellsExt,
                            const std::vector< Vector<REAL> >& Lambdas,
                            const std::vector< Vector<REAL> >& Xzones,
                          
                            const Vector<REAL>& Y_old,// TODO: Check whether this can be exported from YfieldGenerator_old
                            const YFG& YfieldGenerator_old,

                            const VCType_GW& vchead_old,
                            //const VCType_CG& vcM0_old_cg,
                            //const VCType_CG& vcM1_old_cg,

                            UINT iteration_number,
                            const Dune::MPIHelper& helper,
                            std::vector<MyMPIComm> CommunicatorPool,
                            // output:
                            std::vector< Vector<REAL> >& JX,
                            Vector<REAL>& J_times_Y_old
                             ) {
      
      /*
       * some needed variables and typedefinitions
       */
      logger<<"heat_sensitivities(...) "<<std::endl;
      int baselevel =0; // TODO: read from inputdata

      // Retrieve types and data from parameters:
      typedef typename GFS_GW::Traits::GridViewType GV_GW;
      enum{dim=GV_GW::dimension};
      typedef typename GFS_TP::Traits::GridViewType GV_TP;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;

      const GV_GW& gv_0 = pool_grid.levelGridView(0);
      const GV_GW& gv_gw = gfs_gw.gridView();
      const GV_TP& gv_tp = gfs_tp.gridView(); // re-ordered by vchead_old

      // Total number of all cells required to resolve the parameter field
      const UINT nAllCells = inputdata.domain_data.nCells.componentsproduct();
      
      // number of zones
      UINT nzones = inputdata.yfield_properties.nz;
      UINT iSetup = setupdata.index;

 
      logger<<"heat AT SENS: precalculations ... "<<std::endl;


      Dune::Timer watch,watch2;
      watch.reset();
  
  
      // number of heat measurements!
      UINT nPoints_AT = orig_measurements.nMeasPerSetupPerType(iSetup,4);
  
      //number of available communicators. BE CAREFUL "nComm <= inputdata.maxNComm", because inputdata holds the given input for the pool size and this is a max. value!!!! 
      UINT nComm=CommunicatorPool.size();
  
      // just some variables
      Dune::FieldVector<REAL,dim> measure_location(0.0);
      Vector<UINT> read_local_count,read_local_offset;

      typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD_OLD;
      GWP_FWD_OLD gwp_fwd_old( inputdata, setupdata, YfieldGenerator_old );

      typedef GroundwaterAdjointProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_ADJ_OLD;
      GWP_ADJ_OLD gwp_adj_old( inputdata, setupdata, YfieldGenerator_old );
      
      typedef GradientVectorField<GWP_FWD_OLD,GFS_GW> DARCY_FLUX_DGF;
      
      // darcy_flux dgf of the old head
      DARCY_FLUX_DGF darcyflux_dgf( gwp_fwd_old, 
                                    gfs_gw, 
                                    vchead_old,
                                    baselevel );

      L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);

      // New: Point source with adjustable smearing!
      typedef PointFunctionSource<GFS_TP,IDT> PointSourceTypeTP;
      PointSourceTypeTP source_m0adj( gfs_tp, 
                                      inputdata, 
                                      -1.0, // negative pointsource
                                       setupdata.heat_inversion_data.regularization_factor,
                                       setupdata.heat_inversion_data.bfixedwidth );
      
      typedef HeatTransportProblemAdjoint<
        GV_TP,
        REAL,
        DARCY_FLUX_DGF,
        IDT,
        SDT> TP_ADJ;
      TP_ADJ tp_m0adj(darcyflux_dgf,inputdata,setupdata);
      typedef TransportEquation<GFS_TP,
                                DARCY_FLUX_DGF,
                                PointSourceTypeTP,
                                TP_ADJ,
                                IDT,
                                SDT> TPE_M0ADJ;


      // Setup stiffness matrix for adjoint TPE, only once for all measuring points:
      
      TPE_M0ADJ tpe_m0adj( gfs_tp,
                           inputdata,
                           setupdata,
                           darcyflux_dgf,
                           tp_m0adj,
                           source_m0adj,
                           Passenger::heat,
                           EQ::adjoint );  // TODO: matrix hierarchy recycling!



      /*
       * for m0 adjoint depending on m1 adjoint
       */
#ifdef USE_FEM
      typedef FunctionSource<GWP_ADJ_OLD,GFS_TP,IDT,SDT> FunctionSourceADJ;
      FunctionSourceADJ source_m0adj_m1adj( gwp_adj_old, 
                                            gfs_tp, 
                                            inputdata, 
                                            setupdata, 
                                            0,
                                            Passenger::heat
                                            );
#else
      typedef FunctionSource<GWP_ADJ_OLD,GFS_CG,IDT,SDT> FunctionSourceADJ;
      FunctionSourceADJ source_m0adj_m1adj( gwp_adj_old, 
                                            gfs_cg, 
                                            inputdata, 
                                            setupdata, 
                                            0,
                                            Passenger::heat
                                            );
#endif






      typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,FunctionSourceADJ,TP_ADJ,IDT,SDT> TPE_M0ADJ_M1ADJ;

      TP_ADJ tp_m0m1adj(darcyflux_dgf,inputdata,setupdata);

      // matrix assembly for m0 adjoint
      TPE_M0ADJ_M1ADJ tpe_m0adj_m1adj( gfs_tp,
                                       inputdata,
                                       setupdata,
                                       darcyflux_dgf,
                                       tp_m0m1adj,
                                       source_m0adj_m1adj,
                                       Passenger::heat,
                                       EQ::adjoint );


      /*
       * define head adjoint
       */
#ifdef USE_FEM
      typedef FunctionSource<GWP_ADJ_OLD,GFS_TP,IDT,SDT>
        FunctionSourceCombiADJ;
      FunctionSourceCombiADJ source_hadj_m( gwp_adj_old,
                                            gfs_tp,
                                            inputdata,
                                            setupdata,
                                            pool_grid.maxLevel()
                                            );
#else
      typedef FunctionSource<GWP_ADJ_OLD,GFS_GW,IDT,SDT>
        FunctionSourceCombiADJ;
      FunctionSourceCombiADJ source_hadj_m( gwp_adj_old, 
                                            gfs_gw, 
                                            inputdata,
                                            setupdata,
                                            baselevel 
                                            );
#endif


      if( helper.rank()==0 && inputdata.verbosity > 3 )
        std::cout << "=== h_adj wrt heat M1/M0" << std::endl;
      
      typedef GroundWaterEquation<GFS_GW,GWP_ADJ_OLD,FunctionSourceCombiADJ,IDT,SDT> 
        GWE_H_ADJ_M;
      
      // stiffness matrix is assembled here
      GWE_H_ADJ_M gwe_hadj_m( gfs_gw,
                              inputdata,
                              setupdata,
                              gwp_adj_old,
                              source_hadj_m,
                              EQ::adjoint,
                              setupdata.flow_equation.bRecycleMatrixHierarchy );


      // Read heatM0 and headM1
      VCType_CG vcM0_old_cg( gfs_cg, 0.0 );
      HDF5Tools::read_BackendVector_from_HDF5( gfs_cg,
                                               inputdata,
                                               dir.vcheatM0_old_h5file[iSetup],
                                               "/vcheatM0_old",
                                               vcM0_old_cg,
                                               true // preserve_structure
                                               );

      VCType_CG vcM1_old_cg( gfs_cg, 0.0 );
      HDF5Tools::read_BackendVector_from_HDF5( gfs_cg,
                                               inputdata,
                                               dir.vcheatM1_old_h5file[iSetup],
                                               "/vcheatM1_old",
                                               vcM1_old_cg,
                                               true // preserve_structure
                                               );


#ifdef VTK_PLOT_PSI_heat_transport_OFF
      std::stringstream vtu_head_old;
      vtu_head_old << dir.vchead_adj_given_heatM0heatM1_prefix
                   << "_old_head"
                   << "_i" <<  iteration_number;
      VTKPlot::output2vtu( gfs_gw, 
                           vchead_old,
                           vtu_head_old.str(), 
                           "vchead_old", 
                           inputdata.verbosity, 
                           true, 
                           0 
                           );

      std::stringstream vtu_m0_old;
      vtu_m0_old << dir.vchead_adj_given_heatM0heatM1_prefix
                 << "_oldM0"
                 << "_i" <<  iteration_number;
      VTKPlot::output2vtu( gfs_cg, 
                           vcM0_old_cg,
                           vtu_m0_old.str(), 
                           "vcM0_old_cg", 
                           inputdata.verbosity, 
                           true, 
                           0 
                           );
      std::stringstream vtu_m1_old;
      vtu_m1_old << dir.vchead_adj_given_heatM0heatM1_prefix
                 << "_oldM1"
                 << "_i" <<  iteration_number;
      VTKPlot::output2vtu( gfs_cg, 
                           vcM1_old_cg,
                           vtu_m1_old.str(), 
                           "vcM1_old_cg", 
                           inputdata.verbosity, 
                           true, 
                           0 
                           );
#endif


      // loop over all heat AT measurements points
      for( UINT meas_point=0; meas_point < nPoints_AT; meas_point++) {
      
        UINT global_meas_id = orig_measurements.get_global( iSetup, 4, meas_point);
        /*
         * distribute the measurements to the available MPI pools!!!! 
         */
        if(CommunicatorPool[meas_point%nComm].I_am_in()){
       
          watch.reset();
          logger<<"working on meas #"<<global_meas_id<<" (in parallel)"<<std::endl;
          
          // actual measurement location
          measure_location[0] = orig_measurements[global_meas_id].x;
          measure_location[1] = orig_measurements[global_meas_id].y;
#ifdef DIMENSION3
          measure_location[2] = orig_measurements[global_meas_id].z;
#endif
          
          // vector of the sensitivity
          //Vector<REAL> Sensitivity(gv_gw.size(0),0.0);
          Vector<REAL> Sensitivity( nAllCells, 0.0);

          // If the corresponding Sensitivity.h5 file exists there is 
          // no need to waste time doing the calculation once again!
          if( !General::bFileExists( dir.Sensitivity_h5file[global_meas_id] ) ){

            //logger<<"calculation of heatM0 sensitivity "<<std::endl;
       
            /*
             * solve for heat m0 adjoint
             */
            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << "=== Solve for heat m0 adjoint with point source at " << measure_location << std::endl;

            VCType_TP vc_m0adj( gfs_tp, 0.0 );
            tpe_m0adj.solve_adjoint( measure_location, vc_m0adj );

            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << tpe_m0adj.show_ls_result() << std::endl;

            General::logMinAndMax( vc_m0adj, gv_tp.comm() );

            VCType_CG vc_m0adj_cg( gfs_cg, 0.0 );
            L2sp.apply( vc_m0adj, vc_m0adj_cg,
                        inputdata.transport_parameters.l2_diffusion_adjoint );
            General::logMinAndMax( vc_m0adj_cg, gv_tp.comm() );
 
#ifdef VTK_PLOT_PSI_heat_transport
            std::stringstream vtu_m0adj;
            vtu_m0adj << dir.vcheatM0_adj_prefix << "_m" << global_meas_id << "_i" << iteration_number;
            VTKPlot::output2vtu( gfs_cg, 
                                 vc_m0adj_cg, 
                                 vtu_m0adj.str(), 
                                 "heat_M0adj", 
                                 inputdata.verbosity, 
                                 true, 
                                 0 );
#endif


            /*
             * solve for adjoint head ( given m0 and m0 adjoint)
             */
            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << "=== Solve for h adjoint wrt heat m0 adjoint with point source at " << measure_location << std::endl;
            VCType_GW vc_hadj_m0( gfs_gw, 0.0 );

#ifdef USE_FEM
            gwe_hadj_m.solve_adjoint( vc_m0adj_cg, // first
                                      vcM0_old_cg, // second 
                                      vc_hadj_m0   // output
                                      );
#else

            VCType_GW vc_m0adj_gw( gfs_gw, 0.0 );
            //std::cout << "gfs_gw.size() = " << gfs_gw.size() << std::endl;

            CoarseGridP0Projector<GFS_CG,GFS_GW,IDT,dim> sp_ccfv(gfs_cg,gfs_gw,inputdata);
            sp_ccfv.apply( vc_m0adj_cg, 
                           pool_grid.maxLevel(),
                           baselevel, 
                           vc_m0adj_gw );
          
            VCType_GW vcM0_old_gw( gfs_gw, 0.0 );
            sp_ccfv.apply( vcM0_old_cg, 
                           pool_grid.maxLevel(),
                           baselevel, 
                           vcM0_old_gw );
          
            gwe_hadj_m.solve_adjoint( vc_m0adj_gw, // first
                                      vcM0_old_gw, // second
                                      vc_hadj_m0   // output
                                      );
#endif

      
#ifdef VTK_PLOT_PSI_heat_transport
            std::stringstream vtu_hadj_m0;
            vtu_hadj_m0 << dir.vchead_adj_given_heatM0_prefix
                        << "_m" << global_meas_id 
                        << "_i" << iteration_number;
            VTKPlot::output2vtu( gfs_gw, 
                                 vc_hadj_m0, 
                                 vtu_hadj_m0.str(), 
                                 "hadj_heatm0", 
                                 inputdata.verbosity, 
                                 true, 
                                 0 );
#endif
      
            // gradient of the head adjoint
            DARCY_FLUX_DGF grad_hadj_m0_dgf( gwp_fwd_old, // dummy place-holder!
                                             gfs_gw, 
                                             vc_hadj_m0, 
                                             baselevel, 
                                             true //switches off Darcy, compute gradient
                                             );
    
            // vector for the sensitivities
            Vector<REAL> heatM0Sensitivity(gv_0.size(0),0.0);
            Vector<REAL> iheatM0Sensitivity;
            get_m0_Sensitivities( gv_0,
                                  gv_gw,
                                  darcyflux_dgf,
                                  grad_hadj_m0_dgf,
                                  gfs_cg,
                                  vc_m0adj_cg,
                                  vcM0_old_cg,
                                  pool_grid.maxLevel(),
                                  baselevel,
                                  2*pMAX,
                                  heatM0Sensitivity,
                                  iheatM0Sensitivity
                                  );
          
            //logger<<"... Done !"<<" took: "<<watch2.elapsed()<<" sec!" <<std::endl;


          
            //logger<<"calculation of heatM1 sensitivity "<<std::endl;
            /*
             * solve for m1 adjoint
             */
            // NOTE: 
            // We can re-use vc_m0adj_cg because m1adjoint == m0adjoint with the very same point source!
            // This serves as a functional RHS for the solution of "m0adjoint given m1adjoint" 
            /*
             * solve m0 adjoint given m1 adjoint (== m0 adjoint)
             */
            VCType_CG vc_m1adj_theta_cg = vc_m0adj_cg;
            General::vc_times_thetaRT(vc_m1adj_theta_cg,gv_tp,inputdata);


            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << "=== Solve for m0 adjoint wrt heat m1 adjoint with point source at " << measure_location << std::endl;
            VCType_TP vc_m0adj_m1adj( gfs_tp, 0.0 );
            tpe_m0adj_m1adj.solve_adjoint( vc_m1adj_theta_cg, // input
                                           vc_m0adj_m1adj     // output
                                           );
            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << tpe_m0adj_m1adj.show_ls_result() << std::endl;
            General::logMinAndMax( vc_m0adj_m1adj, gv_tp.comm() );

            VCType_CG vc_m0adj_m1adj_cg( gfs_cg, 0.0 );
            L2sp.apply( vc_m0adj_m1adj, vc_m0adj_m1adj_cg,
                        inputdata.transport_parameters.l2_diffusion_adjoint );

            General::logMinAndMax( vc_m0adj_m1adj_cg, gv_tp.comm() );

#ifdef VTK_PLOT_PSI_heat_transport
            std::stringstream vtu_m0adj_m1adj;
            vtu_m0adj_m1adj << dir.vcheatM0_adj_given_heatM1_prefix
                            << "_m" << global_meas_id 
                            << "_i" << iteration_number;

            VTKPlot::output2vtu( gfs_cg, 
                                 vc_m0adj_m1adj_cg,
                                 vtu_m0adj_m1adj.str(), 
                                 "heat_M0adj_M1adj", 
                                 inputdata.verbosity, 
                                 true, 
                                 0 
                                 );
#endif




            /*
             * solve for head adjoint (given m0, m1, m0 adjoint, and m1 adjoint)
             */
            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << "=== Solve for h adjoint wrt heat m0/m1 adjoint with point source at " << measure_location << std::endl;
            VCType_GW vc_hadj_m0m1( gfs_gw, 0.0 );

#ifdef USE_FEM
            gwe_hadj_m.solve_adjoint( vc_m0adj_m1adj_cg // first
                                      , vcM0_old_cg     // second
                                      , vc_m0adj_cg     // == vc_m1adj_cg
                                      , vcM1_old_cg     // fourth
                                      , vc_hadj_m0m1    // output
                                      );
#else

            VCType_GW vc_m0adj_m1adj_gw( gfs_gw, 0.0 );
            sp_ccfv.apply( vc_m0adj_m1adj_cg, 
                           pool_grid.maxLevel(),
                           baselevel, 
                           vc_m0adj_m1adj_gw );

            VCType_GW vcM1_old_gw( gfs_gw, 0.0 );
            sp_ccfv.apply( vcM1_old_cg, 
                           pool_grid.maxLevel(),
                           baselevel, 
                           vcM1_old_gw );

            gwe_hadj_m.solve_adjoint( vc_m0adj_m1adj_gw // first
                                      , vcM0_old_gw     // second
                                      , vc_m0adj_gw     // == vc_m1adj_gw
                                      , vcM1_old_gw     // fourth
                                      , vc_hadj_m0m1    // output
                                      );
#endif

            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << gwe_hadj_m.show_ls_result() << std::endl;

      
#ifdef VTK_PLOT_PSI_heat_transport
            std::stringstream vtu_hadj_m0m1;
            vtu_hadj_m0m1 << dir.vchead_adj_given_heatM0heatM1_prefix
                          << "_m" << global_meas_id << "_i" <<  iteration_number;
            VTKPlot::output2vtu( gfs_gw, 
                                 vc_hadj_m0m1,
                                 vtu_hadj_m0m1.str(), 
                                 "hadj_heatm0m1", 
                                 inputdata.verbosity, 
                                 true, 
                                 0
                                 );
#endif
            // watch2.reset();
            // logger << "calculate sensitivity heat m1... " << std::endl;
          
            DARCY_FLUX_DGF grad_hadj_m0m1_dgf( gwp_fwd_old, // dummy place-holder!
                                               gfs_gw, 
                                               vc_hadj_m0m1, 
                                               baselevel, 
                                               true //compute gradient only
                                               );

            // vector of the sensitivity
            Vector<REAL> heatM1Sensitivity(gv_0.size(0),0.0);
            Vector<REAL> iheatM1Sensitivity;
            get_m1_Sensitivities( gv_0,
                                  gv_gw,
                                  darcyflux_dgf,
                                  grad_hadj_m0m1_dgf,
                                  gfs_cg,
                                  vc_m0adj_m1adj_cg,
                                  vc_m0adj_cg, // == vc_m1adj_cg
                                  vcM0_old_cg,
                                  vcM1_old_cg,
                                  pool_grid.maxLevel(),
                                  baselevel,
                                  2*pMAX,
                                  heatM1Sensitivity,
                                  iheatM1Sensitivity
                                  );
          
            // logger<<"... Done !"<<" took: "<<watch2.elapsed()<<" sec!" <<std::endl;

            logger<<"measurements.get_value_heatM0( "<< iSetup <<" , "<<meas_point<<" ) : "
                  << measurements.get_value_heatM0( iSetup, meas_point )<<std::endl;
            logger<<"measurements.get_value_heatM1( "<< iSetup <<" , "<<meas_point<<" ) : "
                  << measurements.get_value_heatM1( iSetup, meas_point )<<std::endl;
          
            AT_sensitivity_field( gv_0
                                  , measurements.get_value_heatM0( iSetup, meas_point )
                                  , heatM0Sensitivity
                                  , measurements.get_value_heatM1( iSetup, meas_point )
                                  , heatM1Sensitivity
                                  , Sensitivity);

          
      
            // writing the sensitivity to disk in HDF5, needed because the inversion kernel will try to read them from the disk
            HDF5Tools::write_parallel_to_HDF5( gv_0
                                               , inputdata
                                               , Sensitivity
                                               , "/Sensitivity"
                                               , inputdata.domain_data.nCells
                                               , dir.Sensitivity_h5file[global_meas_id]
                                               );

          }
          else {
            logger << "heat AT sensitivities: Read existing Sensitivity field from " 
                   << dir.Sensitivity_h5file[global_meas_id]
                   << std::endl;

            Vector<UINT> local_count;
            Vector<UINT> local_offset;
            HDF5Tools::read_parallel_from_HDF5( gv_0
                                                , inputdata
                                                , Sensitivity
                                                , "/Sensitivity"
                                                , local_count
                                                , local_offset
                                                , dir.Sensitivity_h5file[global_meas_id]
                                                );
          }

          REAL s1norm = Sensitivity.one_norm();
          s1norm = gv_gw.comm().sum( s1norm );
          if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
            std::cout << "====> m" << global_meas_id 
                      << ": " << measure_location 
                      << ", heat AT Sensitivity (total) = " << s1norm << std::endl;
          }



          
          // wait for all processes, only on the current MPI pool
          if(CommunicatorPool[meas_point%nComm].get_size()>1)
            MPI_Barrier(CommunicatorPool[meas_point%nComm].get_comm());
      
          //calculate cross_covariance_Parallel!
          watch2.reset();
          logger<<"calculating JQ PARALLEL ...."<<std::endl;
     
          //calculate JQ!
          cross_covariance_JQ( inputdata, 
                               nCellsExt,
                               Lambdas,
                               dir.Sensitivity_h5file[global_meas_id],
                               dir.JQ_h5file[global_meas_id],
                               CommunicatorPool[meas_point%nComm].get_comm(),
                               global_meas_id,
                               iteration_number,
                               dir
                               );
    
          logger<<"..JQ calculation PARALLEL DONE!!! (took: "<<watch2.elapsed()<<" sec)"<<std::endl;
      
          logger<<"working on meas #"<<global_meas_id<<" (parallel part DONE! in "<<watch.elapsed()<<" sec. )"<<std::endl;
      
          // sequential part of the function -> Only P0 (in the responsible MPI pool) works now
          //   calculated the needed stuff, e.g., cross covariance
          if(CommunicatorPool[meas_point%nComm].get_rank()==0){ // something to do for P0
            logger<<"working on meas #"<<global_meas_id<<" (sequential part)"<<std::endl;
            watch.reset();
	
            // if parallel MPI pool then get the sensitivity data from the disk
            if( CommunicatorPool[meas_point%nComm].get_size()>1){
              logger<<"read the sensitivity sequential"<<std::endl; 
              HDF5Tools::read_sequential_from_HDF5( Sensitivity
                                                    , "/Sensitivity"
                                                    , read_local_count
                                                    , read_local_offset
                                                    , dir.Sensitivity_h5file[global_meas_id]
                                                    , inputdata
                                                    );
              if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
                std::cout << "====> m" << global_meas_id 
                          << ": " << measure_location 
                          << ", heat AT Sensitivity (seq.read) = " << Sensitivity.one_norm() << std::endl;
              }

            }
        
            for(UINT iCell=0; iCell<nAllCells; iCell++) {
              for(UINT ii=0; ii<nzones; ii++)
                JX[global_meas_id][ii] += Sensitivity[iCell] * Xzones[ii][iCell];//X[iCell][ii];//1.0; //ZonationMatrixX[ iCell ];
              J_times_Y_old[global_meas_id] += Sensitivity[ iCell ] * Y_old[ iCell ];
            }
         
            logger<<"working on meas #"<<global_meas_id<<" (sequential part DONE!!! in "<<watch.elapsed()<<" sec. )"<<std::endl;
      
          }// END: if(CommunicatorPool[meas_point%nComm].get_rank()==0)  
       
        }//END: if(CommunicatorPool[meas_point%nComm].I_am_in())

      } //END: loop over all meas points


    }// END: void AT_heat_tomography_sensitivities(...)

  } // GeoInversion

} // Dune

#endif // #ifndef DUNE_HEAT_SENSITIVITIES_HH
