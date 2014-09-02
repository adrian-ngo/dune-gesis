/* 
 * File:    ForwardSimulator.hh
 * Authors: Ronnie Schwede and Adrian Ngo, 2010-2013
 */

#ifndef FORWARDSIMULATOR_HH
#define FORWARDSIMULATOR_HH


#include <dune/gesis/common/my_gfs_utilities.hh>


// include header file for source terms
#include <dune/gesis/sourceterms/source_head_transport.hh>

// include header files for solving the flow problem (gridview gv_gw)
#include <dune/gesis/BVPs/GroundWaterParameters.hh>
#include <dune/gesis/FEM/GWE_Galerkin.hh>
#include <dune/gesis/DG/GWE_CCFV.hh>

// include header files for reordering grid elements (gridview gv_tp)
#include <dune/gesis/DG/dgf_pressurefield.hh>
#include <dune/gesis/DG/rt0_pressurefield.hh>
#include <dune/gesis/DG/reorderedgridview.hh>

// include header files for solving the transport problem (gridview gv_tp)
#include <dune/gesis/BVPs/TransportParameters.hh>
#include <dune/gesis/BVPs/HeatTransportParameters.hh>
#include <dune/gesis/FEM/TPE_SDFEM.hh>
#include <dune/gesis/DG/TPE_DG.hh>

#include <dune/gesis/BVPs/MeshPlecletNumber.hh>

namespace Dune {

  namespace GeoInversion {

    //
    // The Forward Simulator is required for taking measurements on the synthetic test field
    // and on the simulated trial fields.
    // It serves as a wrapper for the equation + source classes to reduce code.
    //
    // Plan for both uniform and adaptive refine:
    // let the gridviews gv_gw and gv_tp be updated outside and be passed via the 
    // gridfunction spaces gfs_gw and gfs_tp!
    //
    template<typename GWP_FWD,  // <--- YField and BC conditions
             typename GEP_FWD,
             //typename FIELD_GW,
             typename GFS_GW,
             typename GFS_TP,
             typename GFS_CG,
             typename IDT,
             typename SDT,
             typename DIR
             >
    class ForwardSimulator {
  
    public:
      // Retrieve types from gridfunction spaces:
      typedef typename GFS_GW::Traits::GridViewType GV_GW;
      typedef typename GFS_TP::Traits::GridViewType GV_TP;
      enum{dim=GV_GW::dimension};
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;

      // Choose types for the forward GWE:
      typedef PumpingSource<GV_GW,REAL,SDT> PumpingSourceTypeGW;
      typedef GroundWaterEquation<GFS_GW,GWP_FWD,PumpingSourceTypeGW,IDT,SDT> GWE_H;

      // Choose type for the Darcy flux and gradient vector:
      //typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_DGF;
      typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_BASE;
      typedef DarcyVelocityCache<GWP_FWD,GFS_GW,DARCY_FLUX_BASE> DARCY_FLUX_DGF;

      // Choose types for the forward TPE for m0:
      typedef TransportProblemM0<GV_TP,REAL,DARCY_FLUX_DGF,IDT,SDT> TPM0;
      typedef SolutePumpingSource<GV_TP,REAL,IDT,SDT> SolutePumpingSourceType;

      typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,SolutePumpingSourceType,TPM0,IDT,SDT> TPE_M0;
      // Choose types for the forward TPE for m1:
      typedef TransportProblemM1<GV_TP,REAL,DARCY_FLUX_DGF,IDT,SDT> TPM1;
      typedef FunctionSource<GWP_FWD, 
                             GFS_CG, //GFS_TP,
                             IDT,SDT> FunctionSourceType;
      typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,FunctionSourceType,TPM1,IDT,SDT> TPE_M1_M0;


      typedef HeatPumpingSource<GV_TP,REAL,IDT,SDT> HeatPumpingSourceType;

      typedef HeatTransportProblemM0<GV_TP,REAL,DARCY_FLUX_DGF,IDT,SDT> HEAT_TPM0;
      typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,HeatPumpingSourceType,HEAT_TPM0,IDT,SDT> HEAT_TPE_M0;
      
      typedef HeatTransportProblemM1<GV_TP,REAL,DARCY_FLUX_DGF,IDT,SDT> HEAT_TPM1;
      typedef TransportEquation<GFS_TP,DARCY_FLUX_DGF,FunctionSourceType,HEAT_TPM1,IDT,SDT> HEAT_TPE_M1_M0;


      // Geo-electrical potential:
      // Choose types for the forward GEP:
      typedef GeoelectricalPotentialSource<GV_GW,REAL,SDT> DipoleSourceType;
      typedef GEP_Equation<GFS_GW,GEP_FWD,DipoleSourceType,IDT,SDT> GEP_EQ;
      // Choose type for the Darcy flux and gradient vector:
      //typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_DGF;
        
#ifdef USE_FEM
      typedef GEP_FunctionSource<GEP_FWD, 
                                 GFS_CG,   //<-- eqn (66): m_k
                                 GFS_GW,   //<-- eqn (66): phi0
                                 IDT,SDT> FunctionSourceType_Kappa;
#else
      typedef GEP_FunctionSource<GEP_FWD, 
                                 GFS_GW,   //<-- eqn (66): m_k
                                 GFS_GW,   //<-- eqn (66): phi0
                                 IDT,SDT> FunctionSourceType_Kappa;
#endif

      typedef GEP_Equation<GFS_GW,GEP_FWD,FunctionSourceType_Kappa,IDT,SDT> Moment_GEP_EQ;


    private:
      GWP_FWD& gwp_fwd; // Take care: This must be non-const because ConvectionDiffusionBoundaryConditionAdapter expects a non-const parameter class
      GEP_FWD& gep_fwd; // dito.
      //FIELD_GW& log_kappafield_gw;
      GEP_FWD& log_kappafield_gw;
      
      const GFS_GW& gfs_gw;
      const IDT& inputdata;
      const SDT& setupdata;
      const DIR& dir;
      GV_GW gv_gw;
      const int qorder;


    public:
    
      //constructor
      ForwardSimulator< GWP_FWD,
                           GEP_FWD,
                           //FIELD_GW,
                           GFS_GW,
                           GFS_TP,
                           GFS_CG,
                           IDT,
                           SDT,
                           DIR >( GWP_FWD& gwp_fwd_,
                                  GEP_FWD& gep_fwd_,
                                  //FIELD_GW& log_kappafield_gw_,
                                  GEP_FWD& log_kappafield_gw_,
                                  const GFS_GW& gfs_gw_,
                                  const IDT& inputdata_,
                                  const SDT& setupdata_,
                                  const DIR& dir_) 
      : gwp_fwd(gwp_fwd_),
        gep_fwd(gep_fwd_),
        log_kappafield_gw(log_kappafield_gw_),
        gfs_gw(gfs_gw_),
        inputdata(inputdata_),
        setupdata(setupdata_),
        dir(dir_),
        gv_gw(gfs_gw.gridView()),
        qorder(2) {
      }
  

  
      void head_simulation(VCType_GW& vchead) {

        Dune::Timer watch;

        logger<<"Forward Simulation: head"<<std::endl;    
        /*
         *  This solves the forward flow problem for the synthethic Y-Field:
         */
        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== Forward head" << std::endl;
        PumpingSourceTypeGW source_h( setupdata );
        GWE_H gwe_h(gfs_gw,
                    inputdata,
                    setupdata,
                    gwp_fwd,
                    source_h);

        General::log_elapsed_time( watch.elapsed(),
                                   gv_gw.comm(),
                                   inputdata.verbosity,
                                   "GWE",
                                   "Matrix Pattern + Assembly" );

        gwe_h.solve_forward( vchead );

        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << gwe_h.show_ls_result() << std::endl;
        General::logMinAndMax( vchead, gv_gw.comm() );
      }
    

      void soluteM0_simulation( const DARCY_FLUX_DGF & darcyflux_dgf,
                                const GFS_TP& gfs_tp,
                                const GFS_CG& gfs_cg,
                                VCType_TP & vcM0,
                                VCType_CG & vcM0_cg
                                ){
        Dune::Timer watch;

        logger<<"Forward Simulation: solute M0"<<std::endl; 
        GV_TP gv_tp = gfs_tp.gridView();
    
        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== Foward M0" << std::endl;
        SolutePumpingSourceType source_m0( inputdata, setupdata );
        TPM0 tpm0(darcyflux_dgf,inputdata,setupdata);

        //std::stringstream vtu_Peclet;
        //vtu_Peclet << dir.vtudir << "/meshPecletM0";
        MeshPlecletNumber<GV_TP,TPM0> meshPeclet( gv_tp, tpm0 );
        //meshPeclet.plot2vtu( vtu_Peclet.str() );

        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){

          std::cout << "Maximal mesh Peclet numbers: Pe_l = " 
                    << meshPeclet.maximum_l() << " / Pe_t = " 
                    << meshPeclet.maximum_t() << std::endl;

#ifdef USE_FEM
          std::cout << "Maximal deltaSDFEM = " 
                    << meshPeclet.maximum_deltaSDFEM() << std::endl;
#endif
        }

        TPE_M0 tpe_m0( gfs_tp
                       , inputdata
                       , setupdata
                       , darcyflux_dgf
                       , tpm0
                       , source_m0 );

        General::log_elapsed_time( watch.elapsed(),
                                   gv_tp.comm(),
                                   inputdata.verbosity,
                                   "TPE",
                                   "Matrix Pattern + Assembly for M0" );


        tpe_m0.solve_forward( vcM0 );    
        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << tpe_m0.show_ls_result() << std::endl;
        General::logMinAndMax( vcM0, gv_tp.comm() );

        L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);

        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== L2 Projection of Foward M0" << std::endl;
        L2sp.apply( vcM0, vcM0_cg,
                    inputdata.transport_parameters.l2_diffusion );
        General::logMinAndMax( vcM0_cg, gv_tp.comm() );

      }


      // TODO: template<typename GFS, typename VCType>
      void soluteM1_simulation( DARCY_FLUX_DGF &darcyflux_dgf,
                                const GFS_TP& gfs_tp,
                                const GFS_CG& gfs_cg,
                                VCType_CG vcM0, //VCType_TP vcM0, 
                                VCType_TP &vcM1M0,
                                VCType_CG &vcM1M0_cg
                                ){

        Dune::Timer watch;
    
        logger<<"Forward Simulation: solute M1"<<std::endl;
        GV_TP gv_tp = gfs_tp.gridView();

        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== Foward M1" << std::endl;
        FunctionSourceType source_m1_m0( gwp_fwd, 
                                         gfs_cg,
                                         inputdata,
                                         setupdata );    
        TPM1 tpm1(darcyflux_dgf,inputdata,setupdata);
    
        TPE_M1_M0 tpe_m1_m0( gfs_tp,
                             inputdata,
                             setupdata,
                             darcyflux_dgf,
                             tpm1,
                             source_m1_m0 );

        General::log_elapsed_time( watch.elapsed(),
                                   gv_tp.comm(),
                                   inputdata.verbosity,
                                   "TPE",
                                   "Matrix Pattern + Assembly for M1" );
              
        // Define sourceterm for the forward transport equation for m1:
        //VCType_TP vc_theta_m0 = vcM0; // theta_m0 := zone_porosity * m0
        VCType_CG vc_theta_m0 = vcM0; // theta_m0 := zone_porosity * m0
        General::vc_times_theta(vc_theta_m0,gv_tp,inputdata);
        tpe_m1_m0.set_rhs( vc_theta_m0 );
        tpe_m1_m0.solve_forward( vcM1M0 );
        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << tpe_m1_m0.show_ls_result() << std::endl;
        General::logMinAndMax( vcM1M0, gv_tp.comm() );

        L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);

        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== L2 Projection of Foward M1" << std::endl;
        L2sp.apply( vcM1M0, vcM1M0_cg,
                    inputdata.transport_parameters.l2_diffusion );
        General::logMinAndMax( vcM1M0_cg, gv_tp.comm() );

      }



      void heat_simulation( const DARCY_FLUX_DGF & darcyflux_dgf,
                            const GFS_TP& gfs_tp,
                            const GFS_CG& gfs_cg,
                            VCType_TP & vc_heatM0,
                            VCType_TP & vc_heatM1,
                            VCType_CG & vc_heatM0_cg,
                            VCType_CG & vc_heatM1_cg
                            ) {        
        logger<<"Forward Simulation: Heat M0 and M1..."<<std::endl; 
        GV_TP gv_tp = gfs_tp.gridView();

        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== Foward Heat M0" << std::endl;
        HeatPumpingSourceType source_heat_m0(inputdata,setupdata);

        HEAT_TPM0 heat_tpm0(darcyflux_dgf,inputdata,setupdata);
        //std::stringstream vtu_Peclet;
        //vtu_Peclet << dir.vtudir << "/meshPecletHeatM0";
        MeshPlecletNumber<GV_TP,HEAT_TPM0> meshPeclet( gv_tp, heat_tpm0 );
        //meshPeclet.plot2vtu( vtu_Peclet.str() );
        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){

          std::cout << "Maximal mesh Peclet numbers (heat): Pe_l = " 
                    << meshPeclet.maximum_l() << " / Pe_t = " 
                    << meshPeclet.maximum_t() << std::endl;

        }

        HEAT_TPE_M0 heat_tpe_m0( gfs_tp
                                 , inputdata
                                 , setupdata
                                 , darcyflux_dgf
                                 , heat_tpm0
                                 , source_heat_m0 
                                 , Passenger::heat
                                 );
        heat_tpe_m0.solve_forward( vc_heatM0 );
        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << heat_tpe_m0.show_ls_result() << std::endl;
        General::logMinAndMax( vc_heatM0, gv_tp.comm() );
        
        L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);
        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== L2 Projection of Foward heat M0" << std::endl;

        L2sp.apply( vc_heatM0, vc_heatM0_cg,
                    inputdata.transport_parameters.l2_diffusion );
        General::logMinAndMax( vc_heatM0_cg, gv_tp.comm() );


        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== Foward Heat M1" << std::endl;
        FunctionSourceType source_heat_m1_m0( gwp_fwd, 
                                              gfs_cg,
                                              inputdata,
                                              setupdata,
                                              0,// m0 and m1 are on the same grid level
                                              Passenger::heat
                                              );
        HEAT_TPM1 heat_tpm1(darcyflux_dgf,inputdata,setupdata);
        HEAT_TPE_M1_M0 heat_tpe_m1_m0( gfs_tp,
                                       inputdata,
                                       setupdata,
                                       darcyflux_dgf,
                                       heat_tpm1,
                                       source_heat_m1_m0,
                                       Passenger::heat
                                       );

        VCType_CG vc_theta_m0 = vc_heatM0_cg;

        General::vc_times_thetaRT(vc_theta_m0,gv_tp,inputdata);
        heat_tpe_m1_m0.set_rhs( vc_theta_m0 );
        heat_tpe_m1_m0.solve_forward( vc_heatM1 );
        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << heat_tpe_m1_m0.show_ls_result() << std::endl;
        General::logMinAndMax( vc_heatM1, gv_tp.comm() );

        if( gv_tp.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << std::endl << "=== L2 Projection of Foward heat M1" << std::endl;
        L2sp.apply( vc_heatM1, vc_heatM1_cg,
                    inputdata.transport_parameters.l2_diffusion );
        General::logMinAndMax( vc_heatM1_cg, gv_tp.comm() );

      }



      void GE_phi0_simulation( const UINT config, VCType_GW& vc_phi0 ){
        DipoleSourceType source_phi0( setupdata, config );
        GEP_EQ gpe_phi0( gfs_gw,
                         inputdata,
                         setupdata,
                         gep_fwd,
                         source_phi0 );
        gpe_phi0.solve_forward( vc_phi0 );
        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << gpe_phi0.show_ls_result() << std::endl;
        General::logMinAndMax( vc_phi0, gv_gw.comm() );
      }


      void GE_simulation( const GFS_CG& gfs_cg,
                          const VCType_CG& vcM0_cg,
                          const VCType_CG& vcM1_cg,
                          const VCType_GW& vc_phi0,
                          VCType_GW& vc_M0phi,
                          VCType_GW& vc_M1phi
                          ){


#ifdef USE_FEM

        FunctionSourceType_Kappa source_moment( log_kappafield_gw,
                                                gfs_cg,   // eqn (66): m_k  
                                                gfs_gw,   // eqn (66): phi0 
                                                inputdata,
                                                setupdata
                                                //, baselevel 
                                                );

        Moment_GEP_EQ moment_gpe( gfs_gw,
                                  inputdata,
                                  setupdata,
                                  gep_fwd,
                                  source_moment );
        
        moment_gpe.solve_forward( vcM0_cg, // eqn (66): m_k
                                  vc_phi0, // eqn (66): phi0
                                  vc_M0phi // output
                                  );
        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << moment_gpe.show_ls_result() << std::endl;
        General::logMinAndMax( vc_M0phi, gv_gw.comm() );

        moment_gpe.solve_forward( vcM1_cg, // eqn (66): m_k 
                                  vc_phi0, // eqn (66): phi0
                                  vc_M1phi // output        
                                  );
        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << moment_gpe.show_ls_result() << std::endl;
        General::logMinAndMax( vc_M1phi, gv_gw.comm() );

#else

        const GV_TP& gv_tp = gfs_cg.gridView();
        int baselevel=0; // TODO: read from inputfile

        CoarseGridP0Projector<GFS_CG,GFS_GW,IDT,dim> sp_ccfv(gfs_cg,gfs_gw,inputdata);
        VCType_GW vcM0_gw( gfs_gw, 0.0 );
        sp_ccfv.apply( vcM0_cg, 
                       gv_tp.grid().maxLevel(),
                       baselevel, 
                       vcM0_gw );

        VCType_GW vcM1_gw( gfs_gw, 0.0 );
        sp_ccfv.apply( vcM1_cg, 
                       gv_tp.grid().maxLevel(),
                       baselevel, 
                       vcM1_gw );

        FunctionSourceType_Kappa source_moment( log_kappafield_gw,
                                                gfs_gw,    // eqn (66): m_k  
                                                gfs_gw,    // eqn (66): phi0 
                                                inputdata,
                                                setupdata
                                                //, baselevel 
                                                );

        Moment_GEP_EQ moment_gpe( gfs_gw,
                                  inputdata,
                                  setupdata,
                                  gep_fwd,
                                  source_moment );
        
        moment_gpe.solve_forward( vcM0_gw, // eqn (66): m_k 
                                  vc_phi0, // eqn (66): phi0
                                  vc_M0phi // output        
                                  );
        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << moment_gpe.show_ls_result() << std::endl;
        General::logMinAndMax( vc_M0phi, gv_gw.comm() );

        moment_gpe.solve_forward( vcM1_gw, // eqn (66): m_k 
                                  vc_phi0, // eqn (66): phi0
                                  vc_M1phi // output        
                                  );
        if( gv_gw.comm().rank()==0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY )
          std::cout << moment_gpe.show_ls_result() << std::endl;
        General::logMinAndMax( vc_M1phi, gv_gw.comm() );
#endif
      }

    }; // class ForwardSimulator

  }

}
#endif /*FORWARDSIMULATOR_HH*/
