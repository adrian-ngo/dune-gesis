/*
 * This function is to calculate the sensitivity of geoelectrical potential measurement wrt to lnK.
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
 * log_sigma0     log base electrical conductivity
 * log_kappa      log kappa (linearization field)
 * qorder:        needed for the GWE_ADJ
 * iteration_number: the current iteration number (only for some VTK output needed)
 * helper:        Dune MPI helper
 * CommunicatorPool : vector containing all information about the different communicators! 
 * HT_setup:      the current HT setup
 * previous_HT:   offset value due to previous HT setups
 * 
 * 
 * OUTPUTS
 * JX:            value for the multiplication of the sensitivity times the trend matrix (at the current meas location)
 * J_times_Y_old: value for the multiplication of the sensitivity times old Y field (at the current meas location)
 */

#ifndef DUNE_GESIS_GEOELECTRICAL_POTENTIAL_SENSITIVITIES_HH
#define DUNE_GESIS_GEOELECTRICAL_POTENTIAL_SENSITIVITIES_HH
#include "../common/MyMPIComm.hh"

// for the internal calculation of the sensitivity of m0 and m1, respectively!
#include "m0_sensitivity_field.hh"
#include "m1_sensitivity_field.hh"
#include "AT_sensitivity_field.hh"

//needed for the use of the "sensitivity_head_internal" function!
//#include "head_sensitivities.hh"

namespace Dune {
  namespace Gesis {

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
void geoelectrical_potential_sensitivities( // input:
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

                                           const Vector<REAL>& Y_old,
                                           const YFG& YfieldGenerator_old,
                                           const YFG& log_electricalConductivity,
                                           const YFG& log_kappa,

                                           const VCType_GW& vchead_old,
                                           const VCType_CG& vcM0_old_cg,
                                           const VCType_CG& vcM1_old_cg,

                                           UINT iteration_number,
                                           const Dune::MPIHelper& helper,

                                           std::vector<MyMPIComm> CommunicatorPool,
                                           // output
                                           std::vector< Vector<REAL> >& JX,
                                           Vector<REAL>& J_times_Y_old
                                            ) {

  /*
   * some needed variables and typedefinitions
   */
  logger<<"GEP_sensitivities(...) "<<std::endl;
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

  const GV_GW& gv_0 = pool_grid.levelGridView(0);
  const GV_GW& gv_gw = gfs_gw.gridView();
  const GV_TP& gv_tp = gfs_tp.gridView(); // re-ordered by vchead_old

  // Total number of all cells required to resolve the parameter field
  const UINT nAllCells = inputdata.domain_data.nCells.componentsproduct();
      
  // number of zones
  UINT nzones = inputdata.yfield_properties.nz;
  UINT iSetup = setupdata.index;
      
  logger<<"GEP SENS: precalculations ... "<<std::endl;

  Dune::Timer watch,watch2;
  watch.reset();

  // number of geoelectrical potential measurements measurements for the current setup!
  UINT nGP_config = orig_measurements.nGE_config(iSetup);


  typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD_OLD;
  GWP_FWD_OLD gwp_fwd_old( inputdata, setupdata, YfieldGenerator_old );

  typedef GradientVectorField<GWP_FWD_OLD,GFS_GW> DARCY_FLUX_DGF;
  DARCY_FLUX_DGF darcyflux_dgf( gwp_fwd_old, 
                                gfs_gw, 
                                vchead_old,
                                baselevel );

  if(inputdata.plot_options.vtk_plot_q ){
    VTKPlot::output_dgf_to_vtu( gv_gw, 
                                gfs_gw, 
                                darcyflux_dgf, 
                                dir.q_orig_vtu[iSetup] + "_TEST", 
                                "q_orig", 
                                inputdata.verbosity, 
                                true,
                                0 );
  }

  // some more variables
  Vector<UINT> read_local_count,read_local_offset;
  
  VCType_GW vc_phi_adj(gfs_gw,0.0);
  VCType_GW vcphi0(gfs_gw,0.0);
  
  typedef GeoelectricalPointSource<GV_GW,REAL,IDT> GeoelectricalPointSourceType;
  typedef GEPotentialAdjointProblem<GV_GW,REAL,IDT,SDT,YFG> GEP_ADJ;

  /*
   * Define the geoelectrical potential adjoint equation
   */
  typedef GEP_Equation
    <GFS_GW,
     GEP_ADJ,
     GeoelectricalPointSourceType,
     IDT,
     SDT
     > GEP_EQ_ADJ;

  GeoelectricalPointSourceType source_phi_adj;
  GEP_ADJ gep_adj( inputdata, setupdata, log_electricalConductivity );
  GEP_EQ_ADJ gep_eq_adj( gfs_gw,
                         inputdata,
                         setupdata,
                         gep_adj,
                         source_phi_adj,
                         EQ::adjoint,
                         setupdata.geoelectrical_potential_equation.bRecycleMatrixHierarchy );
  
  typedef TransportProblemAdjoint<GV_TP,
                                  REAL,
                                  DARCY_FLUX_DGF,
                                  IDT,
                                  SDT> TP_ADJ;

  typedef GEPotentialForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GEP_FWD;
  GEP_FWD log_kappawrapper_gw( inputdata, setupdata, log_kappa ); // equation (66), (147): kappa
#ifdef USE_FEM
  typedef GEP_FunctionSource<GEP_FWD,
                             GFS_GW,   //<-- eqn (150): psi_phi
                             GFS_GW,   //<-- eqn (150): phi0
                             IDT,SDT> FunctionSourceType_Kappa;
  FunctionSourceType_Kappa source_m0_adj( log_kappawrapper_gw,
                                          gfs_gw,
                                          gfs_gw,
                                          inputdata,
                                          setupdata );
#else
  typedef GEP_FunctionSource<GEP_FWD,
                             GFS_GW,   //<-- eqn (150): psi_phi
                             GFS_GW,   //<-- eqn (150): phi0
                             IDT,SDT> FunctionSourceType_Kappa;
  FunctionSourceType_Kappa source_m0_adj( log_kappawrapper_gw,
                                          gfs_gw,
                                          gfs_gw,
                                          inputdata,
                                          setupdata );
#endif

  /*
   * Define the transport adjoint equation for m0 and m1 with functional term  phi_adj
   */
  
  typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,FunctionSourceType_Kappa,TP_ADJ,IDT,SDT>
    TPE_M0_ADJ; // eqn (150) with k=0

  TP_ADJ tpm0_adj(darcyflux_dgf,inputdata,setupdata);

  // matrix assembly for m0 adjoint        
  TPE_M0_ADJ tpe_m0_adj( gfs_tp,
                         inputdata,
                         setupdata,
                         darcyflux_dgf,
                         tpm0_adj,
                         source_m0_adj,
                         Passenger::solute,
                         EQ::adjoint ); // TODO: matrix hierarchy recycling!

  
  VCType_TP vcM0_adj_phi_adj(gfs_tp,0.0);
  VCType_CG vcM0_adj_phi_adj_cg(gfs_cg,0.0);
  
  /*
   * Define groundwater flow adjoint equation with functional source term
   */
  
  typedef GroundwaterAdjointProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_ADJ_OLD;
  GWP_ADJ_OLD gwp_adj_old( inputdata, setupdata, YfieldGenerator_old );

#ifdef USE_FEM
  typedef FunctionSource<GWP_ADJ_OLD,GFS_TP,IDT,SDT>
    FunctionSourceCombiADJ;
  FunctionSourceCombiADJ source_FS_Yfield( gwp_adj_old,
                                           gfs_tp,
                                           inputdata,
                                           setupdata,
                                           pool_grid.maxLevel()
                                           );
#else
  typedef FunctionSource<GWP_ADJ_OLD,GFS_GW,IDT,SDT>
    FunctionSourceCombiADJ;
  FunctionSourceCombiADJ source_FS_Yfield( gwp_adj_old, 
                                           gfs_gw, 
                                           inputdata,
                                           setupdata,
                                           baselevel );
#endif


  typedef GroundWaterEquation<GFS_GW,GWP_ADJ_OLD,FunctionSourceCombiADJ,IDT,SDT> 
    GWE_H_ADJ_FS;

  // stiffness matrix for head adjoint is assembled here
  GWE_H_ADJ_FS gwe_h_adj_FS( gfs_gw,
                             inputdata,
                             setupdata,
                             gwp_adj_old,
                             source_FS_Yfield,
                             EQ::adjoint,
                             setupdata.flow_equation.bRecycleMatrixHierarchy );
    
  VCType_GW vc_hadj_m0(gfs_gw,0.0);
  VCType_GW vc_hadj_m0m1(gfs_gw,0.0);
   
   
  /*
   * Define transport equation with for m0 with  m1 adjoint as RHS
   */
  
  typedef FunctionSource<GWP_FWD_OLD,GFS_CG,IDT,SDT> FunctionSourceTypeTP;
  typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,FunctionSourceTypeTP,TP_ADJ,IDT,SDT> 
    TPE_M0M1_ADJ; // eqn (150) with k=1
  
  TP_ADJ tp_m0m1_adj(darcyflux_dgf,inputdata,setupdata);
  FunctionSourceTypeTP source_m0m1_adj( gwp_fwd_old, // dummy placeholder
                                        gfs_cg,
                                        inputdata,
                                        setupdata );

  TPE_M0M1_ADJ tpe_m0m1_adj( gfs_tp,
                             inputdata,
                             setupdata,
                             darcyflux_dgf,
                             tp_m0m1_adj,
                             source_m0m1_adj,
                             Passenger::solute,
                             EQ::adjoint 
                             );

  VCType_TP vc_m0adj_m1adj(gfs_tp,0.0);
  VCType_CG vc_m0adj_m1adj_cg(gfs_cg,0.0);

  logger << "GEP SENS: precalculations took " << watch.elapsed() 
         << " sec. (for setup #"<< iSetup <<" )"
         << std::endl;
  
  // number of available communicators. 
  // BE CAREFUL "nComm <= inputdata.maxNComm", 
  // because inputdata holds the given input for the pool size and this is a max. value!!!! 
  UINT nComm = CommunicatorPool.size();

  for(UINT iConfig=0; iConfig<nGP_config; iConfig++){
    /*
     * load phi0 (undisturbed electrical potential)
     */
    if(CommunicatorPool[iConfig%nComm].I_am_in()){
      HDF5Tools::read_BackendVector_from_HDF5( gfs_gw,
                                               inputdata,
                                               dir.vcphi0_orig_h5file[iSetup][iConfig],
                                               "/vcphi0_orig",
                                               vcphi0
                                               ); 
      
      Dune::FieldVector<REAL,dim> measure_location1(0.0);
      Dune::FieldVector<REAL,dim> measure_location2(0.0);
      /*
       * Loop over all head measurements
       */
      
      UINT nPoints_GP = orig_measurements.nGE_Meas_per_setup_per_config(iSetup,iConfig);

      for( UINT meas_point=0; meas_point<nPoints_GP; meas_point++ ){
          
        UINT global_meas_id = orig_measurements.get_global(iSetup,5,meas_point,iConfig);

        watch.reset();
        logger<<"working on meas #"<<global_meas_id<<" (in parallel)"<<std::endl;
       
        // actual measurement location 1 and 2
        measure_location1[0] = orig_measurements[global_meas_id].x;
        measure_location1[1] = orig_measurements[global_meas_id].y;
        measure_location2[0] = orig_measurements[global_meas_id].x2;
        measure_location2[1] = orig_measurements[global_meas_id].y2;
#ifdef DIMENSION3
        measure_location1[2] = orig_measurements[global_meas_id].z;
        measure_location2[2] = orig_measurements[global_meas_id].z2;
#endif   


        // vector of the sensitivity
        Vector<REAL> Sensitivity(gv_gw.size(0),0.0);

        // If the corresponding Sensitivity.h5 file exists there is 
        // no need to waste time doing the calculation once again!
        if( !General::bFileExists( dir.Sensitivity_h5file[global_meas_id] ) ){


          //logger<<"solving GPE adjoint!"<<std::endl;
          gep_eq_adj.solve_adjoint( measure_location1,
                                    measure_location2,
                                    vc_phi_adj );  // eqn: (144)
       
          if( inputdata.plot_options.vtk_plot_gp_adjoint ){
            std::stringstream vtu_phi_adj;
            vtu_phi_adj << dir.vcphi_adj_prefix << "_m" << global_meas_id << "_i" << iteration_number;
            VTKPlot::output2vtu( gfs_gw,
                                 vc_phi_adj,
                                 vtu_phi_adj.str(),
                                 "phi_adj",
                                 inputdata.verbosity, 
                                 true, 
                                 0 );

            std::stringstream vtu_phi_adj_vcphi0;
            vtu_phi_adj_vcphi0 << dir.vcphi_adj_prefix << "_m" << global_meas_id << "_i" << iteration_number;
            VTKPlot::output2vtu( gfs_gw,
                                 vcphi0,
                                 vtu_phi_adj_vcphi0.str() + "_vcphi0",
                                 "vcphi0",
                                 inputdata.verbosity, 
                                 true, 
                                 0 );

          }

        
          //logger<<"solving for adjoint Transport M0 wrt. phi adjoint!"<<std::endl;
          tpe_m0_adj.solve_adjoint( vc_phi_adj, vcphi0, // inputs
                                    vcM0_adj_phi_adj ); // eqn: (150) with k=0

          L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);
          L2sp.apply( vcM0_adj_phi_adj,     // input
                      vcM0_adj_phi_adj_cg,  // output
                      inputdata.transport_parameters.l2_diffusion_adjoint );


          if( inputdata.plot_options.vtk_plot_gp_adjoint ){
            std::stringstream vtu_M0gp_adj_phi_adj;
            vtu_M0gp_adj_phi_adj << dir.vcM0gp_adj_phi_adj_prefix << "_m" << global_meas_id << "_i" << iteration_number;
            VTKPlot::output2vtu( gfs_cg,
                                 vcM0_adj_phi_adj_cg,
                                 vtu_M0gp_adj_phi_adj.str() + "_cg",
                                 "GEP_M0_adj_phi_adj_cg",
                                 inputdata.verbosity, 
                                 true, 
                                 0 );
          }



#ifdef USE_FEM
          gwe_h_adj_FS.solve_adjoint( vcM0_adj_phi_adj_cg, vcM0_old_cg, vc_hadj_m0 );  // eqn: (153) ?
#else
          CoarseGridP0Projector<GFS_CG,GFS_GW,IDT,dim> sp_ccfv(gfs_cg,gfs_gw,inputdata);
          VCType_GW vcM0_adj_phi_adj_gw(gfs_gw,0.0);
          sp_ccfv.apply( vcM0_adj_phi_adj_cg,
                         pool_grid.maxLevel(),
                         baselevel, 
                         vcM0_adj_phi_adj_gw );
          VCType_GW vcM0_old_gw(gfs_gw,0.0);
          sp_ccfv.apply( vcM0_old_cg,
                         pool_grid.maxLevel(),
                         baselevel, 
                         vcM0_old_gw );
          //logger<<"solving for adjoint head wrt M0 "<<std::endl;
          gwe_h_adj_FS.solve_adjoint( vcM0_adj_phi_adj_gw,
                                      vcM0_old_gw,
                                      vc_hadj_m0 );  // eqn: (153) for k=0
#endif


       

          if( inputdata.plot_options.vtk_plot_gp_adjoint ){
            std::stringstream vtu_head_adj_M0;
            vtu_head_adj_M0 << dir.vcheadgp_adj_M0_prefix << "_m" << global_meas_id << "_i" << iteration_number;
            VTKPlot::output2vtu( gfs_gw,
                                 vc_hadj_m0, 
                                 vtu_head_adj_M0.str(),
                                 "headGP_adj_M0_adj", 
                                 inputdata.verbosity, 
                                 true, 
                                 0 );
          }



          // gradient of the head adjoint
          DARCY_FLUX_DGF grad_hadj_m0_dgf( gwp_fwd_old, // dummy place-holder!
                                           gfs_gw, 
                                           vc_hadj_m0, 
                                           baselevel, 
                                           true //switches off Darcy, compute gradient
                                           );
     
          // vector for the sensitivities
          Vector<REAL> phiM0Sensitivity(gv_0.size(0),0.0);
          Vector<REAL> iphiM0Sensitivity;

          // This function is defined in "m0_sensitivity_field.hh":
          get_m0_Sensitivities( gv_0,
                                gv_gw,
                                darcyflux_dgf,
                                grad_hadj_m0_dgf,
                                gfs_cg,
                                vcM0_adj_phi_adj_cg,
                                vcM0_old_cg,
                                pool_grid.maxLevel(),
                                baselevel,
                                2*pMAX,
                                phiM0Sensitivity,
                                iphiM0Sensitivity
                                );
            

          /*
          std::stringstream vtu_file;
          vtu_file << dir.Sensitivity_vtu_prefix << "_phiM0Sensitivity";
          VTKPlot::output_vector_to_vtu( gv_0
                                         , phiM0Sensitivity
                                         , vtu_file.str()
                                         , "phiM0Sensitivity"
                                         , inputdata.verbosity
                                         );
          */


          VCType_CG vc_theta_phiM1( gfs_cg, 0.0 );
          vc_theta_phiM1 = vcM0_adj_phi_adj_cg; // eqn (150): m1_adj(k=1) == m0_adj(k=0)
          General::vc_times_theta( vc_theta_phiM1, gv_tp, inputdata );
       
          tpe_m0m1_adj.solve_adjoint( vc_theta_phiM1,     // RHS
                                      vc_m0adj_m1adj );   // output
            
          L2sp.apply( vc_m0adj_m1adj,      // input
                      vc_m0adj_m1adj_cg,   // output
                      inputdata.transport_parameters.l2_diffusion_adjoint );


          if( inputdata.plot_options.vtk_plot_gp_adjoint ){
            std::stringstream vtu_gep_m0adj_m1adj;
            vtu_gep_m0adj_m1adj << dir.vcM0gp_adj_M1_prefix << "_m" << global_meas_id << "_i" << iteration_number;
            VTKPlot::output2vtu( gfs_cg,
                                 vc_m0adj_m1adj_cg,
                                 vtu_gep_m0adj_m1adj.str() + "_cg",
                                 "GEP_m0adj_m1adj_cg",
                                 inputdata.verbosity, 
                                 true,
                                 0 );
          }

            
          /*
           * solve for adjoint head (GP) ( given m0/m1 and M0/M1 adjoint)
           */
#ifdef USE_FEM
          // eqn: (153) for k=1:
          gwe_h_adj_FS.solve_adjoint( vc_m0adj_m1adj_cg,        // first input   (parameter is very important here!)
                                      vcM0_old_cg,              // second input
                                      vcM0_adj_phi_adj_cg,      // third input vcM0_adj_phi_adj==vcM1_adj_phi_adj
                                      vcM1_old_cg,              // fourth input
                                      vc_hadj_m0m1              // output (solution)
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
          // eqn: (153) for k=1:
          gwe_h_adj_FS.solve_adjoint( vc_m0adj_m1adj_gw,        // first input   (parameter is very important here!)
                                      vcM0_old_gw,              // second input
                                      vcM0_adj_phi_adj_gw,      // third input vcM0_adj_phi_adj==vcM1_adj_phi_adj
                                      vcM1_old_gw,              // fourth input
                                      vc_hadj_m0m1              // output (solution)
                                      );
#endif


          if( inputdata.plot_options.vtk_plot_gp_adjoint ){
            std::stringstream vtu_head_adj_M0M1;
            vtu_head_adj_M0M1 << dir.vcheadgp_adj_M0M1_prefix << "_m" << global_meas_id << "_i" << iteration_number;
            VTKPlot::output2vtu( gfs_gw,
                                 vc_hadj_m0m1,
                                 vtu_head_adj_M0M1.str(),
                                 "headGP_adj_M0M1_adj",
                                 inputdata.verbosity, 
                                 true, 
                                 0 );
          }


          // dgf: of adjoint head
          DARCY_FLUX_DGF grad_hadj_m0m1_dgf( gwp_fwd_old, // dummy place-holder!
                                             gfs_gw, 
                                             vc_hadj_m0m1, 
                                             baselevel, 
                                             true //compute gradient only
                                             );

          Vector<REAL> phiM1Sensitivity(gv_0.size(0),0.0);
          Vector<REAL> iphiM1Sensitivity;
          get_m1_Sensitivities( gv_0,
                                gv_gw,
                                darcyflux_dgf,
                                grad_hadj_m0m1_dgf,
                                gfs_cg,
                                vc_m0adj_m1adj_cg,
                                vcM0_adj_phi_adj_cg,   // eqn (150): m1_adj(k=1) == m0_adj(k=0)
                                vcM0_old_cg,
                                vcM1_old_cg,
                                pool_grid.maxLevel(),
                                baselevel,
                                2*pMAX,
                                phiM1Sensitivity,
                                iphiM1Sensitivity
                                );


          if( inputdata.plot_options.vtk_plot_sensitivities ){
            std::stringstream vtu_file1;
            vtu_file1 << dir.Sensitivity_vtu_prefix << "_phiM1Sensitivity";
            VTKPlot::output_vector_to_vtu( gv_0
                                           , phiM1Sensitivity
                                           , vtu_file1.str()
                                           , "phiM1Sensitivity"
                                           , inputdata.verbosity
                                           );
          }

            
          // calculate sensitivities of the ratio!
       
          logger<<"measurements.get_value_GEM0( " << iSetup << " , " << iConfig << " , " << meas_point << " ) : "<<measurements.get_value_GEM0( iSetup, iConfig, meas_point ) << std::endl;
          logger<<"measurements.get_value_GEM1( " << iSetup << " , " << iConfig << " , " << meas_point << " ) : "<<measurements.get_value_GEM1( iSetup, iConfig, meas_point ) << std::endl;
          // This function is defined in "AT_sensitivity_field.hh"

          AT_sensitivity_field( gv_0,
                                measurements.get_value_GEM0(iSetup,iConfig,meas_point),
                                phiM0Sensitivity,
                                measurements.get_value_GEM1(iSetup,iConfig,meas_point),
                                phiM1Sensitivity,
                                Sensitivity );

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
          logger << "GEP sensitivities: Read existing Sensitivity field from " 
                 << dir.Sensitivity_h5file[global_meas_id]
                 << std::endl;

          Vector<UINT> local_count;
          Vector<UINT> local_offset;
          HDF5Tools::
            read_parallel_from_HDF5( gv_0
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
                    << " with dipole:"
                    << "  pole1 = " << measure_location1
                    << "  pole2 = " << measure_location2
                    << ", GEP Sensitivity (total) = " << s1norm << std::endl;
        }

          
        // Ask Ronnie: Is this not needed here?
        // wait for all processes, only on the current MPI pool
        //if(CommunicatorPool[meas_point%nComm].get_size()>1)
        //  MPI_Barrier(CommunicatorPool[meas_point%nComm].get_comm());

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
    
        /*
         *  DOES NOT WORK -> I CAN NOT EXCLUDE THE OVERLAPING ELEMENTS, GIVES THEN WRONG RESULTS!!!
         *                    -> Thus this needs to be sequential in the MPI pool which is responsible for this measurement
         *    
         for(UINT iCell=0; iCell<Sensitivity.size(); iCell++)
         {
         JX[meas_point] += Sensitivity[iCell] * 1.0; //ZonationMatrixX[ iCell ];
         J_times_Y_old[meas_point] += Sensitivity[ iCell ] * Y_old[ iCell ];
         }
         *
         */
        logger<<"..JQ calculation PARALLEL DONE!!! (took: "<<watch2.elapsed()<<" sec)"<<std::endl;
            
        logger<<"working on meas #"<<global_meas_id<<" (parallel part DONE! in "<<watch.elapsed()<<" sec. )"<<std::endl;
            
        // sequential part of the function -> Only P0 (in the responsible MPI pool) works now
        //   calculated the needed output
        if(CommunicatorPool[iConfig%nComm].get_rank()==0){ // something to do for P0
          logger<<"working on meas #"<<global_meas_id<<" (sequential part)"<<std::endl;
          watch.reset();
    
          // if parallel MPI pool then get the sensitivity data from the disk
          if( CommunicatorPool[iConfig%nComm].get_size()>1){
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
                        << " with dipole:"
                        << "  pole1 = " << measure_location1
                        << "  pole2 = " << measure_location2
                        << ", GEP Sensitivity (seq.read) = " << Sensitivity.one_norm() << std::endl;
            }
          }

          for(UINT iCell=0; iCell<nAllCells; iCell++){
            for(UINT ii=0; ii<nzones; ii++)
              JX[global_meas_id][ii] += Sensitivity[iCell] * Xzones[ii][iCell];
            J_times_Y_old[global_meas_id] += Sensitivity[ iCell ] * Y_old[ iCell ];
          }

          logger<<"working on meas #"<<global_meas_id<<" (sequential part DONE!!! in "<<watch.elapsed()<<" sec. )"<<std::endl;
        }// END: if(CommunicatorPool[iConfig%nComm].get_rank()==0) 
       
      } //END: for (UINT meas_point=0; meas_point<nPoints_GP;meas_point++)
    }//END: if(CommunicatorPool[iConfig%nComm].I_am_in())
      
  } // END: for(iConfig=0; iConfig<nGP_config; iConfig++)
  
}// END: void GEOELECTRICAL_POTENTIAL_sensitivities(...)


  }
}

#endif // DUNE_GESIS_GEOELECTRICAL_POTENTIAL_SENSITIVITIES_HH
