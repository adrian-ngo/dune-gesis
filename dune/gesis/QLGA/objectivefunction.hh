// -*- C++ -*-
/* 
 * File:   objectivefunction.hh
 * Author: ngo
 * 
 * modification by Ronnie
 *
 */

#ifndef _OBJECTIVEFUNCTION_HH
#define	_OBJECTIVEFUNCTION_HH

#include <math.h>
#include "../io/inputfile.hh"
#include "../ForwardSimulator_DG.hh"
#include "../io/IO_routines.hh"


namespace Dune {
  namespace GeoInversion {

    template<typename GRID,
             typename GFS_GW,
             typename GFS_TP,
             typename GFS_CG,
             typename DIR,
             typename IDT,
             //typename SDT,
             typename MEASLIST,
             //typename FieldWrapperType
             typename YFG
             >
    double likelihoodfunction( GRID& theGrid, 
                               const GFS_GW& gfs_gw,
                               const DIR& dir,
                               const IDT& inputdata,
                               MEASLIST& orig_measurements,
                               MEASLIST& measurements,
                               const int old,  //flag if Y_old or Y_try?
                               const YFG& log_electricalConductivity,
                               const YFG& log_kappafield,
                               const Dune::MPIHelper& helper,
                               const int it_counter,
                               const int wCounter )
    {
#ifdef USE_YASP
      UINT baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_UG
      int baselevel = inputdata.domain_data.ug_baselevel;
#endif

#ifdef USE_ALUGRID
      int baselevel = inputdata.domain_data.alugrid_baselevel;
#endif

      // Retrieve types and data from parameters:
      typedef typename GFS_GW::Traits::GridViewType GV_GW;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;


      typedef typename GFS_TP::Traits::GridViewType GV_TP;
      typedef typename GFS_TP::Traits::FiniteElementMapType FEM_HYPER;
      typedef typename GFS_TP::Traits::ConstraintsType CON_TP;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;


      typedef typename GFS_CG::Traits::FiniteElementMapType FEM_CG;
      typedef typename GFS_CG::Traits::ConstraintsType CON_CG;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;


      const GV_GW& gv_0 = theGrid.levelGridView( 0 ); // grid for Y_try

#ifndef USE_DGF_PressureField
      GV_GW gv_gw = gfs_gw.gridView(); // must be non-const, maybe updated for rt0_pressurefield after load balancing
#endif

      // Declare variable for the return value:
      double L_try = 0.0;
  
      //important for taking the correct measurements only for the parallel case needed. -> not needed if measurements_try is given in inversion kernel!!!
      //if(helper.size()>1)
      //  measurements.set_value_zero();
  
      // Retrieve values...
      enum{dim=GV_GW::dimension};

      const UINT nMeas=measurements.nMeasurements();
  
      //std::cout << "DEBUG: nMeas = " << nMeas << std::endl;

      // HERE: field generator is used to import data only!
      // typedef FFTFieldGenerator<IDT,REAL,dim> YFG;
      YFG Yfieldgenerator_try( inputdata,dir, helper.getCommunicator() );  
  
      Vector<REAL> local_Y_old;
      Vector<UINT> local_count,local_offset;
  
      // which should e read. Y_old or Y_try?
      // read the needed Y field from disk
      if (old){
        HDF5Tools::
          read_parallel_from_HDF5( gv_0
                                   , inputdata
                                   , local_Y_old
                                   , "/Y_old"
                                   , local_count
                                   , local_offset
                                   , dir.Y_old_h5file
                                   );
      }else{
        HDF5Tools::
          read_parallel_from_HDF5( gv_0
                                   , inputdata
                                   , local_Y_old
                                   , "/Y_try"
                                   , local_count
                                   , local_offset
                                   , dir.Y_try_h5file
                                   );
      }


      if( it_counter == 0 && inputdata.inversion_parameters.max_iter == 0 ){
        // save the estimated field again with preserved structure
        HDF5Tools::write_parallel_to_HDF5( gv_0
                                           , inputdata
                                           , local_Y_old
                                           , "/Y_old"
                                           , inputdata.domain_data.nCells
                                           , dir.Y_est_h5file
                                           , 1
                                           , FEMType::DG
                                           , 0
                                           , false // bWriteOnOverlap
                                           );
      }


      if(helper.size()>1)
        Yfieldgenerator_try.parallel_import_from_local_vector( local_Y_old, local_count, local_offset );
      else
        Yfieldgenerator_try.import_from_vector( local_Y_old );
  
      logger << "objectivefunction: " << std::endl;
      inputdata.loglistOfAllWellCenters();

      logger << "Call setWellConductivities() for Y_try" << std::endl;
      Yfieldgenerator_try.setWellConductivities( gv_0 /* gv_gw */ );


#ifdef VTK_PLOT_YTRY
      std::stringstream vtu_Y_try;
      vtu_Y_try << dir.Y_try_vtu 
                << "_i"<< it_counter 
                << "_w"<< wCounter;
      Yfieldgenerator_try.plot2vtu( gv_0, vtu_Y_try.str(), "Y_try" );
#endif


      typedef typename IDT::SDT SDT;
      typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD;
      typedef GEPotentialForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GEP_FWD;
  
      // typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_DGF;
      typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_BASE;
      typedef DarcyVelocityCache<GWP_FWD,GFS_GW,DARCY_FLUX_BASE> DARCY_FLUX_DGF;
  

      // define forward simulator
      typedef ForwardSimulator_DG<GWP_FWD,
                                  GEP_FWD,
                                  GFS_GW,
                                  GFS_TP,
                                  GFS_CG,
                                  IDT,
                                  SDT,
                                  DIR> FWD_SIM;


      MPI_Barrier(helper.getCommunicator());  
  
      /*
       * take measurements of lnK (if needed)
       */
      if(inputdata.problem_types.lnK_inversion){
        measurements.take_measurements_lnK(Yfieldgenerator_try,gv_0,helper);
      }
      
      // Loop over the setups ...
      for(UINT iSetup=0; iSetup<inputdata.setups.size(); iSetup++) {

        const SDT& setupdata = inputdata.setups[iSetup];

        GWP_FWD gwp_fwd_try( inputdata, setupdata, Yfieldgenerator_try );
        GEP_FWD gep_fwd( inputdata, setupdata, log_electricalConductivity );
        GEP_FWD log_kappawrapper_gw( inputdata, setupdata, log_kappafield ); // equation (66)

        FWD_SIM forward_try( gwp_fwd_try,
                             gep_fwd,
                             log_kappawrapper_gw,
                             gfs_gw,
                             inputdata,
                             setupdata,
                             dir );



        //solve for head and take measurements if needed
        if(( measurements.nMeasPerSetupPerType(iSetup,1) || 
             measurements.nMeasPerSetupPerType(iSetup,2) || 
             measurements.nMeasPerSetupPerType(iSetup,3) || 
             measurements.nMeasPerSetupPerType(iSetup,4) || 
             measurements.nMeasPerSetupPerType(iSetup,5)  ) ){
      
          /*
           *  for head
           */
          //MPI_Barrier(helper.getCommunicator());  
          if( helper.rank()==0 && inputdata.verbosity>1 ){
            std::cout << "=== Solve forward h_try for objective-function" << std::endl;
          }
          VCType_GW vchead_try(gfs_gw, 0.0);
          forward_try.head_simulation(vchead_try);


#ifdef VTK_PLOT_TRIAL_FIELDS
          std::stringstream vtu_head_try;
          vtu_head_try << dir.vtudir << "/h_try"
                       << "_s" << iSetup
                       << "_i" << it_counter 
                       << "_w" << wCounter;
          VTKPlot::output2vtu( gfs_gw, 
                               vchead_try, 
                               vtu_head_try.str(), 
                               "h_try", 
                               inputdata.verbosity, 
                               true, 
                               0
                               );
#endif

          //MPI_Barrier(helper.getCommunicator());  
      
          if(measurements.nMeasPerSetupPerType(iSetup,1))
            measurements.take_measurements( 1, vchead_try, gfs_gw, helper, iSetup);
        
          // save the head
          HDF5Tools::
            write_BackendVector_to_HDF5( gfs_gw
                                         , inputdata
                                         , dir.vchead_old_h5file[iSetup]
                                         , "/vchead_old"
                                         , vchead_try
                                         , PRESERVE_STRUCTURE
                                         , baselevel
                                         );

          if( it_counter == 0 && inputdata.inversion_parameters.max_iter == 0 ){
            // save the head again with preserved structure
            HDF5Tools::
              write_BackendVector_to_HDF5( gfs_gw
                                           , inputdata
                                           , dir.vchead_est_h5file[iSetup]
                                           , "/vchead_old"
                                           , vchead_try
                                           , true // PRESERVE_STRUCTURE
                                           , baselevel
                                           );
          }



          // Calculate the Darcy flux out of the head "vchead_orig"
          DARCY_FLUX_DGF darcyflux_dgf_try( gwp_fwd_try, 
                                            gfs_gw, 
                                            vchead_try,
                                            baselevel );

          // re-order gridview ...
#ifdef NEW_GV

#ifdef USE_DGF_PressureField
          logger << "Using DGF_PressureField ... " << std::endl;
          typedef DGF_PressureField<GFS_GW,VCType_GW> RT0_PF;
#else
          logger << "Using RT0_PressureField ... " << std::endl;
          typedef RT0_PressureField<GV_GW,GFS_GW,GWP_FWD> RT0_PF;
#endif

#ifdef USE_DGF_PressureField
          Dune::shared_ptr<RT0_PF> rt0_pressurefield 
            = Dune::make_shared<RT0_PF>(gfs_gw,vchead_try,baselevel);
#else
          Dune::shared_ptr<RT0_PF> rt0_pressurefield 
            = Dune::make_shared<RT0_PF>(gv_gw,gwp_fwd_try,baselevel);
          rt0_pressurefield->assemble(gfs_gw,vchead_try);
#endif

          PressureLikeOrdering comparefunctor;
          //ReversePressureLikeOrdering comparefunctor;
          GV_TP gv_tp( theGrid.leafGridView(), 
                       rt0_pressurefield,
                       comparefunctor );

#ifdef DEBUG_PLOT_OFF
          std::stringstream vtu_gv_tp;
          vtu_gv_tp << dir.vtudir << "/gv_tp_try"
                    << "_s" << iSetup
                    << "_i" << it_counter 
                    << "_w" << wCounter;
          outputGridviewIndexToDGF( gv_tp, inputdata, vtu_gv_tp.str() );
#endif

#else
          const GV_TP& gv_tp = theGrid.leafGridView(); 
#endif

#ifdef USE_FEM
          FEM_HYPER fem_hyper(gv_tp);
#else
          FEM_HYPER fem_hyper;
#endif
          CON_TP con_tp;
          FEM_CG fem_cg(gv_tp);
          CON_CG con_cg;

          GFS_TP gfs_tp(gv_tp,fem_hyper,con_tp);
          GFS_CG gfs_cg(gv_tp,fem_cg,con_cg );

          CoarseGridP0Projector<GFS_CG,GFS_GW,IDT,dim> sp_ccfv(gfs_cg,gfs_gw,inputdata);


        
          /*
           * for solute and GE
           */
          if(( measurements.nMeasPerSetupPerType(iSetup,2) || 
               measurements.nMeasPerSetupPerType(iSetup,3) || 
               measurements.nMeasPerSetupPerType(iSetup,5)  ) ){

            /*
             * solve for solute M0
             */
            //MPI_Barrier(helper.getCommunicator());  
            if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION ){
              std::cout << "=== Solve forward m0_try for objective-function" << std::endl;
            }
            VCType_TP vcM0_try( gfs_tp, 0.0 );
            VCType_CG vcM0_try_cg( gfs_cg, 0.0 );
            // solve solute transport M0
            forward_try.soluteM0_simulation( darcyflux_dgf_try,
                                             gfs_tp,
                                             gfs_cg,
                                             vcM0_try,
                                             vcM0_try_cg
                                             );


#ifdef VTK_PLOT_TRIAL_FIELDS
            std::stringstream vtu_M0_try;
            vtu_M0_try << dir.vtudir 
                       << "/M0_try_cg"
                       << "_s" << iSetup
                       << "_i" << it_counter 
                       << "_w" << wCounter;
            VTKPlot::output2vtu( gfs_cg, 
                                 vcM0_try_cg, 
                                 vtu_M0_try.str(), 
                                 "m0_try_cg", 
                                 inputdata.verbosity, 
                                 true, 
                                 0 
                                 );
#endif

            if(measurements.nMeasPerSetupPerType(iSetup,2))
              measurements.take_measurements( 2, vcM0_try_cg, gfs_cg, helper, iSetup );

            logger << "DEBUG: write /vcM0_old" << std::endl;
            // save the m0
            HDF5Tools::
              write_BackendVector_to_HDF5( gfs_cg
                                           , inputdata
                                           , dir.vcM0_old_h5file[iSetup]
                                           , "/vcM0_old"
                                           , vcM0_try_cg
                                           , PRESERVE_STRUCTURE
                                           , baselevel
                                           );


            if( it_counter == 0 && inputdata.inversion_parameters.max_iter == 0 ){
              // save m0 again in structured mode
              HDF5Tools::
                write_BackendVector_to_HDF5( gfs_cg
                                             , inputdata
                                             , dir.vcM0_est_h5file[iSetup]
                                             , "/vcM0_old"
                                             , vcM0_try_cg
                                             , true // PRESERVE_STRUCTURE
                                             , baselevel
                                             );
            }
            
            //logger<<"measurements.nMeasPerSetupPerType(iSetup,3) : "<<measurements.nMeasPerSetupPerType(iSetup,3)<<std::endl;
            //logger<<"measurements.nMeasPerSetupPerType(iSetup,5) : "<<measurements.nMeasPerSetupPerType(iSetup,5)<<std::endl;
            
            if( ( measurements.nMeasPerSetupPerType(iSetup,3) || measurements.nMeasPerSetupPerType(iSetup,5)  ) ){
              /*
               * solve for solute M1
               */
              if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION ){
                std::cout << "=== Solve forward m1_try for objective-function" << std::endl;
              }

              VCType_TP vcM1_try( gfs_tp, 0.0 );
              VCType_CG vcM1_try_cg( gfs_cg, 0.0 );
              // solve solute transport M0
              forward_try.soluteM1_simulation(darcyflux_dgf_try,
                                              gfs_tp,
                                              gfs_cg,
                                              vcM0_try_cg,// input
                                              vcM1_try,   // unused output
                                              vcM1_try_cg // output
                                              );

#ifdef VTK_PLOT_TRIAL_FIELDS
              std::stringstream vtu_M1_try;
              vtu_M1_try << dir.vtudir << "/M1_try_cg" 
                         << "_s" << iSetup
                         << "_i" << it_counter 
                         << "_w" << wCounter;
              VTKPlot::output2vtu( gfs_cg,
                                   vcM1_try_cg, 
                                   vtu_M1_try.str(), 
                                   "m1_try_cg", 
                                   inputdata.verbosity, 
                                   true, 
                                   0 
                                   );
#endif
              
              if(measurements.nMeasPerSetupPerType(iSetup,3))
                measurements.take_measurements(3, vcM1_try_cg, gfs_cg, helper, iSetup);
              // save the m0
              HDF5Tools::
                write_BackendVector_to_HDF5( gfs_cg
                                             , inputdata
                                             , dir.vcM1_old_h5file[iSetup]
                                             , "/vcM1_old"
                                             , vcM1_try_cg
                                             , PRESERVE_STRUCTURE
                                             , baselevel
                                             );

              if( it_counter == 0 && inputdata.inversion_parameters.max_iter == 0 ){
                // save m0 again in structured mode
                HDF5Tools::
                  write_BackendVector_to_HDF5( gfs_cg
                                               , inputdata
                                               , dir.vcM1_est_h5file[iSetup]
                                               , "/vcM1_old"
                                               , vcM1_try_cg
                                               , true // PRESERVE_STRUCTURE
                                               , baselevel
                                               );
              }

              /*
               * solve for GE
               */
              if( (measurements.nMeasPerSetupPerType(iSetup,5)  ) ){
                UINT n_GP_config = orig_measurements.nGE_config(iSetup);
                //inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig;
                
                for(UINT iconfig=0; iconfig<n_GP_config; iconfig++){
                    
                  VCType_GW vcphi0(gfs_gw,0.0);
                  VCType_GW vcM0phi_try(gfs_gw,0.0);
                  VCType_GW vcM1phi_try(gfs_gw,0.0);
         
                  // load the phi0 values
                  HDF5Tools::
                    read_BackendVector_from_HDF5( gfs_gw,
                                                  inputdata,
                                                  dir.vcphi0_orig_h5file[iSetup][iconfig],
                                                  "/vcphi0_orig",
                                                  vcphi0
                                                  ); 
           
                  //solve GE moments
                  forward_try.GE_simulation( gfs_cg,   
                                             vcM0_try_cg,
                                             vcM1_try_cg,
                                             vcphi0,
                                             vcM0phi_try,
                                             vcM1phi_try
                                             );

                  //Take measurements !!!
                  logger << "GE: take_measurements_AT of Y_try" << std::endl;
                  measurements.take_measurements_AT( 5,
                                                     vcM0phi_try,
                                                     vcM1phi_try,
                                                     gfs_gw,
                                                     helper,
                                                     iSetup,
                                                     iconfig );

                  HDF5Tools::
                    write_BackendVector_to_HDF5( gfs_gw
                                                 , inputdata
                                                 , dir.vcM0phi_old_h5file[iSetup][iconfig]
                                                 , "/vcM0phi_old"
                                                 , vcM0phi_try
                                                 , PRESERVE_STRUCTURE
                                                 , baselevel
                                                 );

                  HDF5Tools::
                    write_BackendVector_to_HDF5( gfs_gw
                                                 , inputdata
                                                 , dir.vcM1phi_old_h5file[iSetup][iconfig]
                                                 , "/vcM1phi_old"
                                                 , vcM1phi_try
                                                 , PRESERVE_STRUCTURE
                                                 , baselevel
                                                 );

                } // END: for(UINT iconfig=0; iconfig<n_GP_config; iconfig++)
                
              } // END: if( (measurements.nMeasPerSetupPerType(iSetup,5)  ) )

            }  // END: if( ( measurements.nMeasPerSetupPerType(iSetup,3) || measurements.nMeasPerSetupPerType(iSetup,5)  ) )
          } // END: if(( measurements.nMeasPerSetupPerType(iSetup,2) || measurements.nMeasPerSetupPerType(iSetup,3) || measurements.nMeasPerSetupPerType(iSetup,5)  ) )



          if( (measurements.nMeasPerSetupPerType(iSetup,4)  ) ){
            // Calculate the Darcy flux out of the head "vchead_orig"
            //DARCY_FLUX_DGF darcyflux_dgf_try( gfs, vchead_try, Yfieldwrapper_try, inputdata.setups[iSetup]);
            
            VCType_TP vcheatM0_try( gfs_tp, 0.0 );
            VCType_TP vcheatM1_try( gfs_tp, 0.0 );
            VCType_CG vcheatM0_try_cg( gfs_cg, 0.0 );
            VCType_CG vcheatM1_try_cg( gfs_cg, 0.0 );
     
            forward_try.heat_simulation(darcyflux_dgf_try,
                                        gfs_tp,
                                        gfs_cg,
                                        vcheatM0_try,  // unused output
                                        vcheatM1_try,  // unused output
                                        vcheatM0_try_cg, // output
                                        vcheatM1_try_cg  // output
                                        );
            

#ifdef VTK_PLOT_TRIAL_FIELDS
            std::stringstream vtu_heatM0_try;
            vtu_heatM0_try << dir.vtudir << "/heatM0_try"
                           << "_s" << iSetup
                           << "_i" << it_counter 
                           << "_w" << wCounter;
            VTKPlot::output2vtu( gfs_cg, 
                                 vcheatM0_try_cg, 
                                 vtu_heatM0_try.str(), 
                                 "heatm0_try", 
                                 inputdata.verbosity, 
                                 true, 
                                 0
                                 );

            std::stringstream vtu_heatM1_try;
            vtu_heatM1_try << dir.vtudir << "/heatM1_try"
                           << "_s" << iSetup
                           << "_i" << it_counter 
                           << "_w" << wCounter;
            VTKPlot::output2vtu( gfs_cg, 
                                 vcheatM1_try_cg, 
                                 vtu_heatM1_try.str(), 
                                 "heatm1_try", 
                                 inputdata.verbosity, 
                                 true, 
                                 0
                                 );
#endif

            //Take measurements !!!
            logger << "HEAT: take_measurements_AT of Y_try" << std::endl;
            measurements.take_measurements_AT( 4, vcheatM0_try_cg, vcheatM1_try_cg, gfs_cg, helper, iSetup );

            // save the m0
            HDF5Tools::write_BackendVector_to_HDF5( gfs_cg
                                                    , inputdata
                                                    , dir.vcheatM0_old_h5file[iSetup]
                                                    , "/vcheatM0_old"
                                                    , vcheatM0_try_cg
                                                    , PRESERVE_STRUCTURE
                                                    , baselevel
                                                    );

            HDF5Tools::write_BackendVector_to_HDF5( gfs_cg
                                                    , inputdata
                                                    , dir.vcheatM1_old_h5file[iSetup]
                                                    , "/vcheatM1_old"
                                                    , vcheatM1_try_cg
                                                    , PRESERVE_STRUCTURE
                                                    , baselevel
                                                    );

          } // END: if( (measurements.nMeasPerSetupPerType(iSetup,4)  ) )


        } // END: if(( measurements.nMeasPerSetupPerType(iSetup,1) || measurements.nMeasPerSetupPerType(iSetup,2) || measurements.nMeasPerSetupPerType(iSetup,3) || measurements.nMeasPerSetupPerType(iSetup,4) || measurements.nMeasPerSetupPerType(iSetup,5)  ) )
    
      } // end: Loop over the setups
  
      // calculate L_try sequential!
      if(helper.rank()==0){

        Dune::Timer watch;

        for( UINT iMeas=0; iMeas < nMeas; iMeas++ )
          {
            REAL original_value = orig_measurements[ iMeas ].value;
            REAL simulated_value = measurements[ iMeas ].value;

            REAL discrepancy = original_value - simulated_value;

            REAL abs_orig = std::abs( original_value );
            // REAL abs_simu = std::abs( simulated_value );
            //REAL relative_error = std::max( abs_orig, abs_simu ) * orig_measurements[ iMeas ].rel_error;
            REAL relative_error = abs_orig * orig_measurements[ iMeas ].rel_error;

            REAL messfehler = orig_measurements[ iMeas ].abs_error + relative_error;
            //REAL messfehler = std::max( orig_measurements[ iMeas ].abs_error, relative_error );

            REAL contribution = discrepancy / messfehler;
            contribution *= contribution;

            if( inputdata.verbosity >= VERBOSITY_L_CONTRIB ){
              std::cout << "Lm: iMeas=" << iMeas
                        << ", orig=" << original_value
                        << ", simu=" << simulated_value
                        << ", discrepancy=" << discrepancy
                        << ", messfehler=" << messfehler
                        << ", abs_error=" << orig_measurements[ iMeas ].abs_error
                        << ", contribution=" << contribution
                        << std::endl;
            }

            L_try += contribution;

          }

        General::log_elapsed_time( watch.elapsed(),
                                   General::verbosity,
                                   "EVAL",
                                   "Compute L_m" );

      }
      return L_try;
    }








    template<typename GRID,
             typename GFS_GW,
             typename GFS_TP,
             typename GFS_CG,
             typename DIR,
             typename IDT,
             typename YFG,
             typename MEASLIST,
             typename SensCl
             >
    void objectivefunction( GRID& theGrid,
                            const GFS_GW& gfs_gw,
                            const DIR& dir,
                            const IDT& inputdata,
                            MEASLIST& orig_measurements,
                            MEASLIST& measurements,
                            const Vector<REAL>& ksi,
                            const int old,// old=1 uses Y_old, old=0 uses Y_try
                            const YFG& log_electricalConductivity,
                            const YFG& log_kappafield,
                            SensCl& sensitivity,
                            const Dune::MPIHelper& helper,
                            const int it_counter,
                            const int wCounter,
                            REAL& L_m,
                            REAL& L_p
                            ) {

      double Likelihood = likelihoodfunction
        <GRID,
         GFS_GW,
         GFS_TP,
         GFS_CG,
         DIR,
         IDT,
         MEASLIST,
         YFG
         >
        ( theGrid,
          gfs_gw,
          dir,
          inputdata,
          orig_measurements,
          measurements,
          old,
          log_electricalConductivity,
          log_kappafield,
          helper,
          it_counter,
          wCounter
          );

      L_m = Likelihood;

      /*
      if( wCounter > 1 ){
        logger << "objectivefunction: 1 < wCounter = " << wCounter << std::endl;
        logger << "objectivefunction: Do not evaluate L_prior. Do interpolate." << std::endl;
        return;
      }
      */

      if( old ){
        if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION ){
          std::cout << "L_measure = " << Likelihood << std::endl;
        }
        //return Likelihood;
        return;
      }
      else{
        double L_prior=0.0;
        if(helper.rank()==0){

          //std::cout << "DEBUG: objectivefunction ksi_try = "
          //          << ksi << std::endl;

          L_prior = sensitivity.calculate_prior_term( ksi, wCounter );
          /* Note by Adrian:
             Calculated after Nowak's simplified L_p, but be careful:
             The ksi-term is involved in the line-search,
             but the JQ terms are not involved in the line search.
             That's why you cannot expect the new L_p value to converge
             to the old L_p value as weighting tends to zero.
             We do a linear interpolation of Gyy between Gyy_OLD and Gyy_NEW
             inside calculate_prior_term().
          */

          Vector<REAL> tmp(1,L_prior);
          Vector<UINT> dimensions(1,1);

          HDF5Tools::
            write_sequential_to_HDF5_without_DUNE( dimensions,
                                                   tmp,
                                                   "/L_prior",
                                                   dir.L_prior_h5file );

        }
        if(helper.size()>1)
          MPI_Bcast(&(L_prior),1,MPI_DOUBLE,0,helper.getCommunicator());

        L_p = L_prior;

        if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION ){
          // logging only on process 0
          std::cout << "L_measure = " << L_m << std::endl;
          std::cout << "L_prior   = " << L_p << std::endl;
        }

        //return L_prior + Likelihood; // L_prior seems to be buggy!
        //return Likelihood; // --> fwd_50x50x3 test: good!!!
        return;
      }
    }



  } // GeoInversion
}  // Dune
#endif	/* _OBJECTIVEFUNCTION_HH */

