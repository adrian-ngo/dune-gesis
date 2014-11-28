/*
 * This function is to calculate the sensitivity of a m0 measurement wrt to lnK.
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
 * orig_measuring_points: class holding information about the measuring points
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

#ifndef DUNE_GESIS_M0_SENSITIVITIES_HH
#define DUNE_GESIS_M0_SENSITIVITIES_HH
#include "../common/MyMPIComm.hh"


#include "m0_sensitivity_field.hh"


namespace Dune {
  namespace Gesis {


    // nCellsExt : Dimensions of extended field per zone
    // Lambdas   : eigenvalues of the covariance matrix per zone
    // X         : zonation matrix -> only on root process 0 (leader of the communicator group) 
    // 

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
    void m0_sensitivities( // input:
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
                          const VCType_CG& vcM0_old_cg,

                          UINT iteration_number,
                          const Dune::MPIHelper& helper,
                          std::vector<MyMPIComm> CommunicatorPool,
                          // output:
                          std::vector< Vector<REAL> >& JX,
                          Vector<REAL>& J_times_Y_old,
                          UINT& SV_increased 
                           )
    {
      logger << "m0_sensitivities(...)" << std::endl;

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
      enum{dim=GV_GW::dimension};
      typedef typename GFS_TP::Traits::GridViewType GV_TP;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;


      GV_GW gv_0 = pool_grid.levelGridView(0);
      GV_GW gv_gw = gfs_gw.gridView();

      UINT iSetup = setupdata.index;

#ifdef DEBUG_PLOT
      const GV_TP& gv_tp = gfs_tp.gridView();  // re-ordered by vchead_old
      std::stringstream indexrangefile;
      indexrangefile << dir.vcM0_adj_prefix 
                     << "_s" << iSetup
                     << "_i" << iteration_number
                     << "_gv_tp";
      outputGridviewIndexToDGF( gv_tp, inputdata, indexrangefile.str() );
#endif // DEBUG_PLOT

      // Total number of all cells required to resolve the parameter field
      const UINT nAllCells = inputdata.domain_data.nCells.componentsproduct();
      
      // number of zones
      UINT nzones=inputdata.yfield_properties.nz;
 
      logger<<"m0_sensitivities: precalculations ... "<<std::endl;
      Dune::Timer watch;
  
      // number of m0 measurements!
      UINT nPoints_m0=orig_measurements.nMeasPerSetupPerType(iSetup,2);
  
      //number of available communicators. BE CAREFUL "nComm <= inputdata.maxNComm", because inputdata holds the given input for the pool size and this is a max. value!!!! 
      UINT nComm=CommunicatorPool.size();
  
  
      // just some variables
      Dune::FieldVector<REAL,dim> measure_location(0.0);
      Vector<UINT> read_local_count,read_local_offset;
  
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

#ifdef DEBUG_PLOT
      std::stringstream darcyfluxfile;
      darcyfluxfile << dir.vcM0_adj_prefix 
                    << "_s" << iSetup
                    << "_i" << iteration_number
                    << "_qflux";
      VTKPlot::output_dgf_to_vtu( gv_gw, 
                                  gfs_gw, 
                                  darcyflux_dgf, 
                                  darcyfluxfile.str(), 
                                  "q_orig", 
                                  inputdata.verbosity, 
                                  true,
                                  0 );
#endif // DEBUG_PLOT
      
      L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);
      
      /*
       * define m0_adj
       */

      // New: Point source with adjustable smearing!
      typedef PointFunctionSource<GFS_TP,IDT> PointSourceTypeTP;
      PointSourceTypeTP source_m0_adj( gfs_tp, 
                                       inputdata, 
                                       -1.0, // negative pointsource
                                       setupdata.solute_concentration_inversion_data.regularization_factor,
                                       setupdata.solute_concentration_inversion_data.bfixedwidth );
      // TODO: How about heat_inversion_data for m0? heat_flag and well_flag?


      typedef Dune::Gesis::TransportProblemAdjoint<GV_TP,
                                                          REAL,
                                                          DARCY_FLUX_DGF,
                                                          IDT,
                                                          SDT> TP_ADJ;

      TP_ADJ tpm0_adj(darcyflux_dgf,inputdata,setupdata);

#ifdef DEBUG_PLOT
      MeshPlecletNumber<GV_TP,TP_ADJ> meshPeclet( gv_tp, tpm0_adj, 0 );
      std::stringstream vtu_Peclet;
      vtu_Peclet << dir.vtudir << "/peclet_step" << 0;
      meshPeclet.plot2vtu( vtu_Peclet.str() );
      if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
        std::cout << "Maximal mesh Peclet numbers: Pe_l = "
                  << meshPeclet.maximum_l() << " / Pe_t = "
                  << meshPeclet.maximum_t() << std::endl;
      }
#endif


      typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,PointSourceTypeTP,TP_ADJ,IDT,SDT> TPE_M0_ADJ;

      watch.reset();
      // Setup stiffness matrix for adjoint TPE, only once for all measuring points:
      TPE_M0_ADJ tpe_m0_adj( gfs_tp,
                             inputdata,
                             setupdata,
                             darcyflux_dgf,
                             tpm0_adj,
                             source_m0_adj,
                             Passenger::solute,
                             EQ::adjoint ); //TODO: matrix hierarchy recycling!

      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 inputdata.verbosity,
                                 "TPE",
                                 "Matrix Pattern + Assembly, m0 Adjoint" );


     
      /*
       * define head adjoint
       */
      // function source type
#ifdef USE_FEM
      typedef FunctionSource<GWP_ADJ_OLD,GFS_TP,IDT,SDT>
        FunctionSourceCombiADJ;
      FunctionSourceCombiADJ source_hadj_m0( gwp_adj_old,
                                             gfs_tp,
                                             inputdata,
                                             setupdata,
                                             pool_grid.maxLevel()
                                             );
#else
      typedef FunctionSource<GWP_ADJ_OLD,GFS_GW,IDT,SDT>
        FunctionSourceCombiADJ;
      FunctionSourceCombiADJ source_hadj_m0( gwp_adj_old, 
                                             gfs_gw, 
                                             inputdata,
                                             setupdata,
                                             baselevel );
#endif

      typedef GroundWaterEquation<GFS_GW,GWP_ADJ_OLD,FunctionSourceCombiADJ,IDT,SDT> 
        GWE_H_ADJ_M0;


      // stiffness matrix is assembled here, only once for all measuring points!
      watch.reset();
      GWE_H_ADJ_M0 gwe_h_adj_m0( gfs_gw,
                                 inputdata,
                                 setupdata,
                                 gwp_adj_old,
                                 source_hadj_m0,
                                 EQ::adjoint,
                                 setupdata.flow_equation.bRecycleMatrixHierarchy );

      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 inputdata.verbosity,
                                 "GWE",
                                 "Matrix Pattern + Assembly, m0 Adjoint Head" );

      // Loop over all measurements points
      for(UINT meas_point=0; meas_point<nPoints_m0;meas_point++){

        UINT global_meas_id = orig_measurements.get_global(iSetup,2,meas_point);
        
        Vector<INT> SV_procedure(dim-1,0);
        REAL norm_H, norm_H_old=0.0,SV_slope,SV_slope_old=0.0;


        /*
         * distribute the measurements to the available MPI pools!!!! 
         */
        if(CommunicatorPool[meas_point%nComm].I_am_in()){
          
          UINT SV_loop_flag=1;

          while(SV_loop_flag){

            logger << "m0_sensitivities: working on meas #" << global_meas_id << " at location " << measure_location << std::endl;
      

            // actual measurement location
            measure_location[0] = orig_measurements[global_meas_id].x;
            measure_location[1] = orig_measurements[global_meas_id].y;
#ifdef DIMENSION3
            measure_location[2] = orig_measurements[global_meas_id].z;
#endif

            // vector for the sensitivities
            Vector<REAL> Sensitivity( gv_0.size(0), 0.0);
            Vector<REAL> iSensitivity; // i=interior, the size of this vector is to be determined later via an interior_partition loop!

            // If the corresponding Sensitivity.h5 file exists there is 
            // no need to waste time doing the calculation once again!
            if( !General::bFileExists( dir.Sensitivity_h5file[global_meas_id] ) ){      
              /*
               * solve for m0 adjoint
               */
              VCType_TP vc_m0adj( gfs_tp, 0.0 );
#ifdef ADJ_SMEARING_OFF 
              if(SV_procedure[0] 
#ifdef DIMENSION3
                 &&
                 SV_procedure[1]
#endif   
                 ){
                std::vector<REAL> setSV(GV::dimension,1);
                setSV[1]=SV_procedure[0]*setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_step[1];
#ifdef DIMENSION3
                setSV[2]=SV_procedure[1]*setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_step[2];
#endif
                logger<<"m0_sensitivities: set_SV  with: "<<setSV[1];
#ifdef DIMENSION3
                logger<<", "<<setSV[2];
#endif
                logger<<std::endl;
                //source_m0_adj.set_SV(setSV);
              
              }
#endif

              /*

                WARNING: evaluating YFG at some measure_location can be
                dangerous if that measure_location is not in the current
                process partition! This can lead to an out-of-scope index 
                and a segmentation fault.

              REAL k_old = 0.0;
              REAL k_Well_old = 0.0;
#ifdef DIMENSION3
              YfieldGenerator_old.evaluate3d( measure_location, k_old, k_Well_old );
#else
              YfieldGenerator_old.evaluate2d( measure_location, k_old, k_Well_old );
#endif
              logger << "DEBUG: m0_sensitivities: measure_location = "
                     << measure_location
                     << " k_try = " << k_old
                     << " k_Well_try = " << k_Well_old
                     << std::endl;
              */

              // Maybe, we must add option to set setSV here:
              tpe_m0_adj.solve_adjoint( measure_location, vc_m0adj );

              if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
                std::cout << tpe_m0_adj.show_ls_result() << std::endl;
              General::logMinAndMax( vc_m0adj, gv_gw.comm() );

              VCType_CG vc_m0adj_cg( gfs_cg, 0.0 );
              L2sp.apply( vc_m0adj, vc_m0adj_cg,
                          inputdata.transport_parameters.l2_diffusion_adjoint );
              General::logMinAndMax( vc_m0adj_cg, gv_gw.comm() );

              //if(CommunicatorPool[meas_point%nComm].get_size()>1)
              //  MPI_Barrier(CommunicatorPool[meas_point%nComm].get_comm());
#ifdef ADJ_SMEARING_OFF
              if(SV_procedure[0] 
#ifdef DIMENSION3
                 &&
                 SV_procedure[1]
#endif   
                 ){
                // set SV to zero for the next measurement
                std::vector<REAL> setSV(GV::dimension,0);
                //source_m0_adj.set_SV(setSV);
              }
#endif


              if( inputdata.plot_options.vtk_plot_adjoint_m0 ){
                std::stringstream vtu_m0_adj;
                vtu_m0_adj << dir.vcM0_adj_prefix 
                           << "_s" << iSetup
                           << "_m" << global_meas_id 
                           << "_i" << iteration_number;
                VTKPlot::output2vtu( gfs_cg, 
                                     vc_m0adj_cg, 
                                     vtu_m0_adj.str(), 
                                     "m0_Adjoint", 
                                     inputdata.verbosity, 
                                     true, 
                                     0
                                     );
              }
            
              /*
               * solve for adjoint head ( given m0 and m0 adjoint)
               */
              VCType_GW vc_h_adj_m0( gfs_gw, 0.0 );


#ifdef USE_FEM
              gwe_h_adj_m0.solve_adjoint( vc_m0adj_cg,
                                          vcM0_old_cg,    // second input 
                                          vc_h_adj_m0  // output (solution)
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

              gwe_h_adj_m0.solve_adjoint(  vc_m0adj_gw, // first input 
                                           vcM0_old_gw,  // second input 
                                           vc_h_adj_m0   // output (solution)
                                           );
#endif

              if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
                std::cout << gwe_h_adj_m0.show_ls_result() << std::endl;


              if( inputdata.plot_options.vtk_plot_adjoint_m0 ){
                std::stringstream vtu_hadj_m0;
                vtu_hadj_m0 << dir.vchead_adj_given_M0_prefix 
                            << "_s" << iSetup
                            << "_m" << global_meas_id 
                            << "_i" << iteration_number;
                VTKPlot::output2vtu( gfs_gw, 
                                     vc_h_adj_m0, 
                                     vtu_hadj_m0.str(), 
                                     "h_adj_m0", 
                                     inputdata.verbosity, 
                                     true, 
                                     0
                                     );
              }



              // gradient of the head adjoint
              DARCY_FLUX_BASE grad_hadj_m0_dgf( gwp_fwd_old, // dummy place-holder!
                                                gfs_gw, 
                                                vc_h_adj_m0, 
                                                baselevel, 
                                                true //switches off Darcy, compute gradient
                                                );

              logger << "m0_sensitivities: calculate M0 sensitivity ... " << std::endl;

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
                                    Sensitivity,
                                    iSensitivity
                                    );
   
              logger << "====> m" << global_meas_id 
                     << ": " << measure_location 
                     << ", m0 Sensitivity = " 
                     << iSensitivity.one_norm() << std::endl;

      
              // writing the sensitivity to HDF5
              // the inversion kernel will read it later
              /*
                HDF5Tools::h5g_pWrite( gv_0
                , inputdata
                , Sensitivity
                , "/Sensitivity"
                , inputdata.domain_data.nCells
                , dir.aSensitivity_h5file[global_meas_id]
                );
              */

              HDF5Tools::h5g_pWrite( gv_0
                                     , iSensitivity
                                     , dir.Sensitivity_h5file[global_meas_id]
                                     , "/Sensitivity"
                                     , inputdata
                                     , inputdata.domain_data.nCells
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

              Vector<UINT> local_count;
              Vector<UINT> local_offset;
              HDF5Tools::h5g_pRead( gv_0
                                    , Sensitivity
                                    , dir.Sensitivity_h5file[global_meas_id]
                                    , "/Sensitivity"
                                    , local_offset
                                    , local_count
                                    , inputdata
                                    );
            }

            
            REAL s1norm = iSensitivity.one_norm();
            s1norm = gv_gw.comm().sum( s1norm );
            if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
              std::cout << General::getDateAndTime()
                        << "====> m" << global_meas_id 
                        << ": " << measure_location 
                        << ", m0 Sensitivity (total) = " << s1norm << std::endl;
            }


            // wait for all processes, only on the current MPI pool
            if(CommunicatorPool[meas_point%nComm].get_size()>1)
              MPI_Barrier(CommunicatorPool[meas_point%nComm].get_comm());
      
      
            //calculate cross_covariance_Parallel!
            logger<<"m0_sensitivities: calculating JQ PARALLEL ...."<<std::endl;
     
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

            logger<<"m0_sensitivities: working on meas #"<<global_meas_id<<" (parallel part DONE)"<<std::endl;
      
            // sequential part of the function -> Only P0 (in the responsible MPI pool) works now
            //   calculated the needed stuff, e.g., cross covariance
            if(CommunicatorPool[meas_point%nComm].get_rank()==0){ // something to do for P0
              logger<<"m0_sensitivities: working on meas #"<<global_meas_id<<" (sequential part)"<<std::endl;
	
              // if parallel MPI pool then get the sensitivity data from the disk
              if( CommunicatorPool[meas_point%nComm].get_size()>1){
                logger<<"m0_sensitivities: read the sensitivity sequentially"<<std::endl; 
                HDF5Tools::h5g_Read( Sensitivity
                                     , dir.Sensitivity_h5file[global_meas_id]
                                     , "/Sensitivity"
                                     , read_local_offset
                                     , read_local_count
                                     , inputdata
                                     );

                if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
                  std::cout << "====> m" << global_meas_id 
                            << ": " << measure_location 
                            << ", m0 Sensitivity (seq.read) = " << Sensitivity.one_norm() << std::endl;
                }
      

              }
              norm_H = 0.0;
              logger << "norm_H: " << norm_H << std::endl;
        
              //make JX and J_times_Y_old to zero again (needed for the SV procedure, otherwise JX is wrong)
              for(UINT iCell=0; iCell<nAllCells; iCell++){
                for(UINT ii=0; ii<nzones; ii++)
                  JX[global_meas_id][ii]=0.0;
                J_times_Y_old[global_meas_id]=0.0;
              }
        
              for(UINT iCell=0; iCell<nAllCells; iCell++){
                norm_H += std::abs( Sensitivity[iCell] );
                for(UINT ii=0; ii<nzones; ii++)
                  JX[global_meas_id][ii] += Sensitivity[iCell] * Xzones[ii][iCell];//X[iCell][ii];//1.0; //ZonationMatrixX[ iCell ];
                J_times_Y_old[global_meas_id] += Sensitivity[ iCell ] * Y_old[ iCell ];
              }

              logger<<"m0_sensitivities: working on meas #"<<global_meas_id<<" (sequential part DONE)"<<std::endl;
	
              SV_slope=(norm_H-norm_H_old)/2.0;
              if(inputdata.verbosity>1){
                logger<<"norm_H = "<<norm_H<<std::endl;
                logger<<"norm_H_old = "<<norm_H_old<<std::endl;
                logger<<"SV_slope = "<<SV_slope<<std::endl;
                logger<<"SV_slope_old = "<<SV_slope_old<<std::endl;
              }
      
              if(SV_slope_old>SV_slope && SV_slope>1.e-5){ // optimal sampling volume reached
                logger<<"Optimal sampling volume reached at: y-axis ="<< SV_procedure[0]*setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_step[1];
#ifdef DIMENSION3
                logger<<" and z-axis = "<<SV_procedure[1]*setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_step[2];
#endif
                logger<<std::endl;
                SV_loop_flag=0;
              }else{
                if(SV_slope>1.e-3){
                  norm_H_old=norm_H;
                  SV_slope_old=SV_slope;
                }
              }
	
              if(SV_loop_flag){       
                if(setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_max[1]>=(REAL)(SV_procedure[0]+1)*setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_step[1] && fabs(measurements[global_meas_id].value-orig_measurements[global_meas_id].value)> setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_limit && setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_max[1]>0.0){
                  SV_procedure[0]+=1;
                  SV_increased=1;
                }else{
                  SV_loop_flag=0;
                }
#ifdef DIMENSION3
                if(setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_max[2]>=(REAL)(SV_procedure[1]+1)*setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_step[2] && fabs(measurements[global_meas_id].value-orig_measurements[global_meas_id].value)> setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_limit&& setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_max[2]>0.0){
                  SV_procedure[1]+=1; 
                  SV_increased=1;
                  if(SV_loop_flag==0)
                    SV_loop_flag=1;
                }
                if(setupdata.solute_concentration_inversion_data.m0_inversion_data.SV_max[2]<=0.0)
                  SV_loop_flag=0;
#endif
                if(SV_loop_flag==0 && (SV_procedure[0]>0
#ifdef DIMENSION3
                                       || SV_procedure[1]>0
#endif
                                       )){
                  logger<<"Optimal SV NOT reached!!! Using maximum possible!"<<std::endl;
	 
                }

              }
   
            }// END: if(CommunicatorPool[meas_point%nComm].get_rank()==0)   
            //collect SV_procedure and SV_loop_flag on all processors of this communicator group!
            if(CommunicatorPool[meas_point%nComm].get_size()>1){
              MPI_Bcast(&(SV_loop_flag),1,MPI_INT,0,CommunicatorPool[meas_point%nComm].get_comm());
              MPI_Bcast(&(SV_procedure[0]),SV_procedure.size(),MPI_UNSIGNED,0,CommunicatorPool[meas_point%nComm].get_comm());
            }

            // Adrian: switch off SV_loop_flag for the time being
            SV_loop_flag = 0;

          }//END: while(SV_loop_flag)
        }//END: if(CommunicatorPool[meas_point%nComm].I_am_in())
    
      }//END: for (UINT meas_point=0; meas_point<nPoints_head;meas_point++) 

      darcyflux_dgf.emptyCache();
  
    } // END: m0_sensitivities



  }

}

#endif // DUNE_GESIS_M0_SENSITIVITIES_HH
