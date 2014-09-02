/*
 * This function is to calculate the sensitivity of a head measurement wrt to lnK.
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
 * Lambdas:       eigenvalues of the extended covariance matrix
 * Y_old:         previous Y field
 * YfieldGenerator_old: needed for the YFieldGenerator wrapper class
 * qorder:        needed for the GWE_ADJ
 * iteration_number: the current iteration number (only for some VTK output needed)
 * helper:        Dune MPI helper
 * CommunicatorPool : vector containing all information about the different communicators! 
 * 
 * 
 * OUTPUTS
 * JX:            value for the multiplication of the sensitivity times the trend matrix (at the current meas location)
 * J_times_Y_old: value for the multiplication of the sensitivity times old Y field (at the current meas location)
 */

#ifndef DUNE_HEAD_SENSITIVITIES_HH
#define DUNE_HEAD_SENSITIVITIES_HH
#include "../common/MyMPIComm.hh"


#include "head_sensitivity_field.hh"

namespace Dune {
  namespace GeoInversion {

// zonation matrix -> only on process 0 (process leader of the communicator group) 
    
    template<typename GV_GW,
             typename GFS_GW,
             typename VCType_GW,
             typename PGV,
             typename DIR,
             typename MEASLIST,
             typename IDT,
             typename SDT,
             typename YFG
             >
    void head_sensitivities( // input:
                            const GV_GW& gv_0
                            , const GFS_GW& gfs_gw
                            , const VCType_GW& vchead_old
                            , const PGV pRootGridView
                            , const DIR& dir
                            , const MEASLIST& orig_measurements
                            , const IDT& inputdata
                            , const SDT& setupdata
                            
                            , const std::vector< Vector<UINT> >& nCellsExt
                            , const std::vector< Vector<REAL> >& Lambdas
                            , const std::vector< Vector<REAL> >& Xzones 
                            , const Vector<REAL>& Y_old
                            , const YFG& YfieldGenerator_old

                            , UINT iteration_number
                            , const Dune::MPIHelper& helper
                            , std::vector<MyMPIComm> CommunicatorPool
                            // output:
                            , std::vector< Vector<REAL> >& JX 
                            , Vector<REAL>& J_times_Y_old
                            
                             )
    {
      logger << "head_sensitivities:" << std::endl;

#ifdef USE_YASP
      UINT baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_UG
      int baselevel = inputdata.domain_data.ug_baselevel;
#endif

#ifdef USE_ALUGRID
      int baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
      int qorder = 2*pMAX;

 
      /*
       * some needed variables and typedefinitions
       */

      //typedef typename GFS_GW::Traits::GridViewType GV_GW;

      enum{dim=GV_GW::dimension};

      GV_GW gv_gw = gfs_gw.gridView();

      const UINT nAllCells = inputdata.domain_data.nCells.componentsproduct();
      UINT iSetup = setupdata.index;

      Dune::Timer watch;
   
      // number of zones
      UINT nzones=inputdata.yfield_properties.nz;
  
      typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD_OLD;
      GWP_FWD_OLD gwp_fwd_old( inputdata, setupdata, YfieldGenerator_old );

      typedef GroundwaterAdjointProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_ADJ_OLD;
      GWP_ADJ_OLD gwp_adj_old( inputdata, setupdata, YfieldGenerator_old );
  

      typedef PointFunctionSource<GFS_GW,IDT> PointSourceTypeGW;
      REAL pointValue = 1;
      PointSourceTypeGW source_h_adj( gfs_gw,
                                      inputdata, 
                                      pointValue,
                                      setupdata.head_inversion_data.regularization_factor,
                                      setupdata.head_inversion_data.bfixedwidth );

      // define the ground water equation for the adjoint -> build the stiffness matrix
      if( helper.rank()==0 && inputdata.verbosity > VERBOSITY_EQ_SUMMARY )
        std::cout << "Setup GWE_ADJ with point source (head sensitivities)" << std::endl;

      typedef GroundWaterEquation<GFS_GW,GWP_ADJ_OLD,PointSourceTypeGW,IDT,SDT> GWE_ADJ;
    
      watch.reset();
      GWE_ADJ gwe_adj( gfs_gw,
                       inputdata,
                       setupdata,
                       gwp_adj_old,
                       source_h_adj,
                       EQ::adjoint,
                       setupdata.flow_equation.bRecycleMatrixHierarchy );

      General::log_elapsed_time( watch.elapsed(),
                                 gv_gw.comm(),
                                 inputdata.verbosity,
                                 "GWE",
                                 "Matrix Pattern + Assembly, Adjoint Head" );

      // number of head measurements!
      UINT nPoints_head=orig_measurements.nMeasPerSetupPerType(iSetup,1);
   
      //number of available communicators. BE CAREFUL "nComm <= inputdata.maxNComm", because inputdata holds the given input for the pool size and this is a max. value!!!! 
      UINT nComm=CommunicatorPool.size();

      /*
       * some more variables and typedefs
       */
      //typedef GradientVectorField<GWP_FWD_OLD,GFS_GW> DARCY_FLUX_DGF;
      typedef GradientVectorField<GWP_FWD_OLD,GFS_GW> DARCY_FLUX_BASE;
      typedef DarcyVelocityCache<GWP_FWD_OLD,GFS_GW,DARCY_FLUX_BASE> DARCY_FLUX_DGF;


      Vector<UINT> read_local_count,read_local_offset;
      Dune::FieldVector<REAL,dim> measure_location(0.0);

      // darcy_flux dgf of the old head
      DARCY_FLUX_DGF darcyflux_dgf( gwp_fwd_old, 
                                    gfs_gw, 
                                    vchead_old,
                                    baselevel
                                    );

   
      /*
       * Loop over all head measurements
       */
      for (UINT meas_point=0; meas_point<nPoints_head;meas_point++){
     
        /*
         * distribute the measurements to the available MPI pools!!!! 
         */
        if(CommunicatorPool[meas_point%nComm].I_am_in()){
       
          UINT global_meas_id= orig_measurements.get_global(iSetup,1,meas_point);

          // actual measurement location
          measure_location[0] = orig_measurements[global_meas_id].x;
          measure_location[1] = orig_measurements[global_meas_id].y;
#ifdef DIMENSION3
          measure_location[2] = orig_measurements[global_meas_id].z;
#endif


          Vector<REAL> Sensitivity( gv_0.size(0), 0.0);
          Vector<REAL> iSensitivity; // i=interior, the size of this vector is to be determined later via an interior_partition loop!
          
          typedef typename GV_GW::Grid::GlobalIdSet::IdType GlobalIdType;
          Vector<GlobalIdType> iGlobalGridIndexVector;

          // If the corresponding Sensitivity.h5 file exists there is 
          // no need to waste time doing the calculation once again!
          if( !General::bFileExists( dir.Sensitivity_h5file[global_meas_id] ) ){

            logger << "head_sensitivities:: working on meas #"<<global_meas_id<<" (in parallel)"<<std::endl;

            //vector for the adjoint and the sensitivity
            VCType_GW vcAdjoint(gfs_gw, 0.0);

            //Solve the Adjoint Flow Problem here for this measuring_location:
            logger << "=== Solve adjoint GWE at location: " << measure_location << std::endl;

            gwe_adj.solve_adjoint( measure_location, vcAdjoint );

#ifdef DEBUG_PLOT
            std::stringstream pointsource_filename;
            pointsource_filename << dir.vchead_adj_prefix
                                 << "_m" << global_meas_id
                                 << "_i" << iteration_number
                                 << "_PS";
            source_h_adj.plot2vtu( pointsource_filename.str() );
#endif // DEBUG_PLOT


            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
              std::cout << gwe_adj.show_ls_result() << std::endl;
      

#ifdef VTK_PLOT_PSI_HEAD
            std::stringstream vtu_h_adj;
            vtu_h_adj << dir.vchead_adj_prefix << "_m" << global_meas_id << "_i" << iteration_number;
            Dune::GeoInversion::VTKPlot::output2vtu( gfs_gw, 
                                                     vcAdjoint,
                                                     vtu_h_adj.str(),
                                                     "head_adj",
                                                     inputdata.verbosity, 
                                                     true, 
                                                     std::max(0,pMAX-1) );
#endif
      

            logger << "head_sensitivities: calculate sensitivity ... " << std::endl;
       
            // gradient dgf of the adjoint
            DARCY_FLUX_BASE grad_h_adj_dgf( gwp_fwd_old, // dummy place-holder!
                                            gfs_gw, 
                                            vcAdjoint, 
                                            baselevel, 
                                            true // <--- This switches off Darcyflux and computes the gradient field!
                                            );
          
#if defined DEBUG_PLOT && defined VTK_PLOT_PSI_HEAD
            Dune::GeoInversion::VTKPlot::
              output_dgf_to_vtu( gv_gw,
                                 gfs_gw,
                                 grad_h_adj_dgf,
                                 vtu_h_adj.str() + "_GradhAdj",
                                 "GradhAdj",
                                 inputdata.verbosity,
                                 true,
                                 0 );

            Dune::GeoInversion::VTKPlot::
              output_dgf_to_vtu( gv_gw,
                                 gfs_gw,
                                 darcyflux_dgf,
                                 vtu_h_adj.str() + "_Darcy",
                                 "Darcy",
                                 inputdata.verbosity,
                                 true,
                                 0 );
#endif
      
            // calculate the sensitivities: the function is defined in this FILE!!
            Dune::GeoInversion::head_sensitivity_field( gv_0
                                                        , gv_gw
                                                        , darcyflux_dgf
                                                        , grad_h_adj_dgf
                                                        , qorder
                                                        , Sensitivity
                                                        , iSensitivity
                                                        , iGlobalGridIndexVector
                                                        );


      


#if defined DEBUG_PLOT && defined VTK_PLOT_PSI_HEAD
            Dune::GeoInversion::VTKPlot::
              output_vector_to_vtu(
                                   gv_0
                                   , Sensitivity
                                   , vtu_h_adj.str() + "_Sens"
                                   , "hSens"
                                   );
#endif


            // writing the sensitivity to HDF5
            // the inversion kernel will read it later
            /*
            HDF5Tools::write_parallel_to_HDF5( gv_0
                                               , inputdata
                                               , Sensitivity
                                               , "/Sensitivity"
                                               , inputdata.domain_data.nCells
                                               , dir.aSensitivity_h5file[global_meas_id]
                                               );
            */

#ifdef USE_SEQ_WRITE_TO_HDF5
            // Alternative way of writing to HDF5:
            // Send data from all processors to the root process
            // and write in seqential mode to the HDF5 file.
            // This may be faster than a parallel write in certain cases.
            HDF5Tools::write_to_HDF5_on_root( gv_0,
                                              pRootGridView,
                                              inputdata,
                                              iSensitivity,
                                              iGlobalGridIndexVector,
                                              dir.Sensitivity_h5file[global_meas_id],
                                              "/Sensitivity"
                                              );
#else
            HDF5Tools::write_parallel_to_HDF5( gv_0
                                               , inputdata
                                               , iSensitivity
                                               , "/Sensitivity"
                                               , inputdata.domain_data.nCells
                                               , dir.Sensitivity_h5file[global_meas_id]
                                               , 1
                                               , FEMType::DG
                                               , 0
                                               , false // preserve_structure
                                               );
#endif

          }
          else {
            logger << "head_sensitivities: Read existing Sensitivity field from " 
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

          REAL s1norm = iSensitivity.one_norm();
          s1norm = gv_gw.comm().sum( s1norm );
          if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
            std::cout << General::getDateAndTime()
                      << "====> m" << global_meas_id
                      << ": " << measure_location
                      << ", head Sensitivity (total) = " << s1norm << std::endl;
          }
      
          // wait for all processes, only on the current MPI pool
          if(CommunicatorPool[meas_point%nComm].get_size()>1)
            MPI_Barrier(CommunicatorPool[meas_point%nComm].get_comm());
      
      
          //calculate cross_covariance_Parallel!
          logger<<"head_sensitivities: calculating JQ PARALLEL ...."<<std::endl;
    
    
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
	
          /*
           *  DOES NOT WORK -> I CAN NOT EXCLUDE THE OVERLAPING ELEMENTS, GIVES THEN WRONG RESULTS!!!
           *                    -> Thus this needs to be sequential in the MPI pool which is responsible for this measurement
           *    
           for(UINT iCell=0; iCell<Sensitivity.size(); iCell++)
           {
           JX[global_meas_id] += Sensitivity[iCell] * 1.0; //ZonationMatrixX[ iCell ];
           J_times_Y_old[global_meas_id] += Sensitivity[ iCell ] * Y_old[ iCell ];
           }
           *
           */
      
          logger<<"head_sensitivities: working on meas #"<<global_meas_id<<" (parallel part) DONE."<<std::endl;
      
#ifdef _TEST_PAR_READ_SENS_Y_OLD_

          // Adrian: This part is just testwise.
          // In order to parallelize this part, even Y_old needs to be parallelized! Too much for the moment!
          // Better things to do: make adaptive forward solver work, add trial solutions to check-pointing, ...

          // LeafIterator on interior partition:
          typedef typename GV_GW::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator InteriorElementIterator;

          const typename GV_GW::IndexSet& is = gv_gw.indexSet();
          for (InteriorElementIterator it = gv_gw.template begin<0,Dune::Interior_Partition>()
                 ; it != gv_gw.template end<0,Dune::Interior_Partition>()
                 ; ++it) {
            
            int iCell = is.index(*it);

            logger << "DEBUG: iCell = " << iCell << std::endl;
            logger << "DEBUG: Y_old[ iCell ] = " << Y_old[ iCell ] << std::endl;
            logger << "DEBUG: Sensitivity[ iCell ] = " << Sensitivity[ iCell ] << std::endl;

            for(UINT ii=0; ii<nzones; ii++)
              JX[global_meas_id][ii] += Sensitivity[iCell] * Xzones[ii][iCell];
            J_times_Y_old[global_meas_id] += Sensitivity[ iCell ] * Y_old[ iCell ];
          }

#else
          // sequential part of the function -> Only P0 (in the responsible MPI pool) works now
          // calculated the needed output
          if(CommunicatorPool[meas_point%nComm].get_rank()==0){ // something to do for P0
            logger<<"head_sensitivities: working on meas #"<<global_meas_id<<" (sequential part)"<<std::endl;
	
            // if parallel MPI pool then get the sensitivity data from the disk
            if( CommunicatorPool[meas_point%nComm].get_size()>1){
              logger<<"head_sensitivities: read the sensitivity sequential"<<std::endl; 
              HDF5Tools::read_sequential_from_HDF5( Sensitivity
                                                    , "/Sensitivity"
                                                    , read_local_count
                                                    , read_local_offset
                                                    , dir.Sensitivity_h5file[global_meas_id]
                                                    , inputdata
                                                    );

              if( inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
                std::cout << "====> m" << global_meas_id 
                          << ": " << measure_location 
                          << ", head Sensitivity (seq.read) = " << Sensitivity.one_norm() << std::endl;
              }
            }
        
            for(UINT iCell=0; iCell<nAllCells; iCell++){
              for(UINT ii=0; ii<nzones; ii++)
                JX[global_meas_id][ii] += Sensitivity[iCell] * Xzones[ii][iCell];
              J_times_Y_old[global_meas_id] += Sensitivity[ iCell ] * Y_old[ iCell ];
            }

            logger << "head_sensitivities: working on meas #" << global_meas_id <<" (sequential part) DONE." <<std::endl;

          }// END: if(CommunicatorPool[meas_point%nComm].get_rank()==0)
#endif
        }//END: if(CommunicatorPool[meas_point%nComm].I_am_in())

      }//END: for (UINT meas_point=0; meas_point<nPoints_head;meas_point++)

      darcyflux_dgf.emptyCache();

    }// END: void head_sensitivities(...)

  } // namespace GeoInversion
} // namespace Dune

#endif // DUNE_HEAD_SENSITIVITIES_HH
