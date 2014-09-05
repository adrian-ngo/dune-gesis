#ifndef DUNE_GESIS_FUNCTION_SOURCE_HH
#define DUNE_GESIS_FUNCTION_SOURCE_HH

namespace Dune {
  namespace Gesis {



    //=========================================================================
    // The template class to define the source term for forward flow simulations
    //========================================================================
    
    template<typename GWP
             , typename GFS
             , typename IDT
             , typename SDT>
    class FunctionSource
      : public SourceTermInterface<
      SourceTermTraits<typename GFS::Traits::GridViewType,REAL>,
      FunctionSource<GWP,GFS,IDT,SDT>
      > 
    {
    private:
      typedef typename GFS::Traits::GridViewType GV;
      enum{dim=GV::dimension};
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type VCType;

      const GWP& gwp;
      const GFS& gfs;
      const IDT& inputdata;
      const SDT& setupdata;
      VCType vc_source;
      VCType vc_source_1;
      VCType vc_source_2;
      VCType vc_source_3;
      VCType vc_source_4;

      const int maxGridLevel; // Important if TPE grid is finer than the GWE grid

      const Passenger::Type passengerType;
      bool well_flag;
      
      int n_RHS_Functions;

      inline void zone_parameters(const REAL& zone_coord, REAL & porosity, REAL & rho_s, REAL & c_s) const{
	
        UINT inside=0;
        // only one zone!
        if(inputdata.yfield_properties.nz<2){
          porosity=inputdata.yfield_properties.zones[0].porosity;
          rho_s=inputdata.yfield_properties.zones[0].rho;
          c_s=inputdata.yfield_properties.zones[0].c_s;
          return; 
        }
	
        for(UINT ii=0; ii<inputdata.yfield_properties.nz-1; ii++){
          if(zone_coord<=inputdata.yfield_properties.zones[ii].bottom ){
            porosity=inputdata.yfield_properties.zones[ii].porosity;
            rho_s=inputdata.yfield_properties.zones[ii].rho;
            c_s=inputdata.yfield_properties.zones[ii].c_s;
            inside=1;
            break;
          }
        }
        //last zone!
        if(zone_coord>inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-2].bottom && inside==0){
          porosity=inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].porosity;
          rho_s=inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].rho;
          c_s=inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].c_s;
        }

      } 
      


      template<typename DFV0>   // Dune::FieldVector of codim 0
      inline bool isPointInsideReachOfWell( const DFV0& elementpoint, 
                                            const REAL& reach,
#ifdef DIMENSION3
                                            const REAL& reach_y,
#endif
                                            const int& iWell ) const
      {

        REAL well_position_x =  setupdata.wdlist.pointdata_vector[iWell].x;
#ifdef DIMENSION3
        REAL well_position_y =  setupdata.wdlist.pointdata_vector[iWell].y;
#endif
        REAL well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
        REAL well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;


        if(
           elementpoint[0] >= well_position_x - reach
           &&
           elementpoint[0] < well_position_x + reach
           &&
#ifdef DIMENSION3
           elementpoint[1] >= well_position_y - reach_y
           &&
           elementpoint[1] < well_position_y + reach_y
           &&
           elementpoint[2] > well_bottom - GEO_EPSILON
           &&
           elementpoint[2] < well_top + GEO_EPSILON
#else
           elementpoint[1] > well_bottom - GEO_EPSILON
           &&
           elementpoint[1] < well_top + GEO_EPSILON           
#endif
           )
          return true;
        else
          return false;
      }

    public:
    
      typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> DGF;
      typedef Dune::PDELab::DiscreteGridFunctionDarcy<GWP,GFS> DARCY_FLUX_DGF;
      typedef REAL RF;
    
      // refined versions evaluating on leaf elements:
      typedef Dune::PDELab::DiscreteRefinedGridFunction<GFS,VCType> DRGF;

      const SourceNature source_nature;
      typedef SourceTermTraits<GV,REAL> Traits;

      // constructor:
      FunctionSource( const GWP& gwp_
                      , const GFS& gfs_
                      , const IDT& inputdata_
                      , const SDT& setupdata_
                      , const int maxGridLevel_=0
                      , const Passenger::Type passengerType_=Passenger::solute
                      , const bool well_flag_=false
                      )
        :
        gwp(gwp_),
        gfs( gfs_ ),
        inputdata( inputdata_ ),
        setupdata( setupdata_ ),
        vc_source( gfs, 0.0 ), 
        vc_source_1( gfs, 0.0 ), 
        vc_source_2( gfs, 0.0 ),
        vc_source_3( gfs, 0.0 ),
        vc_source_4( gfs, 0.0 ),
        maxGridLevel(maxGridLevel_),
        passengerType(passengerType_),
        well_flag(well_flag_),
        source_nature( FUNCTIONAL_SOURCE )
      {
        logger << "FunctionSource constructor with maxGridLevel = " 
               << maxGridLevel
               << std::endl;
      }

      void reset_rhs()
      {
        vc_source=VCType(gfs,0.0);
        vc_source_1=VCType(gfs,0.0);
        vc_source_2=VCType(gfs,0.0);
        vc_source_3=VCType(gfs,0.0);
        vc_source_4=VCType(gfs,0.0);
      }

      void set_rhs( const VCType& vc_solution )
      {
        // Hint: 
        // Copying a VCType is easier than copying a DGF! 
        // That's why the DGF will be created out of the VCType inside the function evaluate(...)
        // However, I don't know yet whether this will have a negative effect on performance.
        vc_source_1=VCType(gfs,0.0);
        vc_source_2=VCType(gfs,0.0);
        vc_source_3=VCType(gfs,0.0);
        vc_source_4=VCType(gfs,0.0);
        vc_source = vc_solution;
        n_RHS_Functions = 1;
      }


      // for h_adj_m0:
      void set_rhs( const VCType& vc_1
                    , const VCType& vc_2
                    )
      {
        vc_source=VCType(gfs,0.0);
        vc_source_1 = vc_1;
        vc_source_2 = vc_2;
        vc_source_3=VCType(gfs,0.0);
        vc_source_4=VCType(gfs,0.0);
        n_RHS_Functions = 2;
      }

      
      void set_rhs( const VCType& vc_1
                    , const VCType& vc_2
                    , const VCType& vc_3
                    , const VCType& vc_4 )
      {
        vc_source=VCType(gfs,0.0);
        vc_source_1 = vc_1;
        vc_source_2 = vc_2;
        vc_source_3 = vc_3;
        vc_source_4 = vc_4;
        n_RHS_Functions = 4;
      }







      template<
        typename EG
        , typename LFSV
        , typename SF
        , typename RESIDUAL_TYPE
        >
      bool evaluate_residual_on_element( const EG& eg
                                         , const LFSV& lfsv
                                         , SF& shapefunctions
                                         , RESIDUAL_TYPE& residual
                                         , const EQ::Mode equationMode = EQ::forward
                                         //, UINT& flag_source  // Ask Ronnie: Wofür?
                                         ) const
      {
        
        RF rho_w = inputdata.transport_parameters.heat_rho_w;
        RF c_w   = inputdata.transport_parameters.heat_c_w;
          
        UINT nSources=setupdata.pdlist.total;  // adding Pointsources only!

        for( UINT i=0; i<nSources; i++ ) {

          RF input,injection_time,rho_s,c_s,porosity;
          RF pumping_rate = setupdata.pdlist.pointdata_vector[i].value;

          if( equationMode == EQ::adjoint )
            pumping_rate = -pumping_rate;

          if( passengerType == Passenger::heat ){
            input = setupdata.pdlist.pointdata_vector[i].temperature;
            injection_time = setupdata.pdlist.pointdata_vector[i].temperature_injection_time;
          }else{
            input = setupdata.pdlist.pointdata_vector[i].concentration;
            injection_time = setupdata.pdlist.pointdata_vector[i].concentration_injection_time;
          }

          if(injection_time<=GEO_EPSILON)
            injection_time=0.0;
          else
            injection_time *= injection_time*0.5;  // Ask Ronnie: Why this?



          if( pumping_rate > GEO_EPSILON*0.5 && std::abs(input) >  GEO_EPSILON*0.5 ) {

            Dune::FieldVector<RF,dim> pointsource_global;

            pointsource_global[0] = setupdata.pdlist.pointdata_vector[i].x;
            pointsource_global[1] = setupdata.pdlist.pointdata_vector[i].y;
#ifdef DIMENSION3
            pointsource_global[2] = setupdata.pdlist.pointdata_vector[i].z;
            if( passengerType == Passenger::heat )
              zone_parameters(pointsource_global[2], porosity, rho_s, c_s);
#else
            if( passengerType == Passenger::heat )
              zone_parameters(pointsource_global[1], porosity, rho_s, c_s);
#endif

            // get the local coordinate of all the source locations:
            Dune::FieldVector<RF,dim> pointsource_local 
              = eg.geometry().local( pointsource_global );
				
            // Check if the point source lies inside the rectangular cell.
            UINT iOutside = 0;
            for( UINT ii=0; ii<dim; ii++ ){
              if( (pointsource_local[ii] < 0.0) || (pointsource_local[ii] >= 1.0) ){
                iOutside++;
              }
            }

				
            if( iOutside==0 ) {
              lfsv.finiteElement().localBasis().evaluateFunction( pointsource_local, shapefunctions );
            
              RF heat_fac=1.0;

              if( passengerType == Passenger::heat )
                heat_fac=(rho_w*c_w)/( (porosity*rho_w*c_w)+((1.0-porosity)*rho_s*c_s) );
            
              RF factor=eg.geometry().volume();
              RF npoints=eg.geometry().corners();
	    
              //UINT flag_source=1; // Ask Ronnie: Wofür?
              for( UINT ii=0; ii<lfsv.size(); ii++ ) {
                residual.accumulate( lfsv, ii, (-pumping_rate/factor*npoints)*input*injection_time*(heat_fac)* shapefunctions[ii]);
              }
            }
          }
        }
      

#ifdef WELL_FRACTURE_MODEL 

        nSources=setupdata.wdlist.total;
   
        for( UINT i=0; i<nSources; i++ ) {
          Dune::FieldVector<RF,dim> global=eg.geometry().center();
          //std::cout<<"i:"<<i<<std::endl;
          std::vector<REAL> meshsize(inputdata.domain_data.yasp_gridsizes);
          std::vector<REAL> meshsize_2(inputdata.domain_data.yasp_gridsizes);
          for(UINT ii=0; ii<meshsize_2.size();ii++)
            meshsize_2[ii]*=0.5;

          RF well_x = setupdata.wdlist.pointdata_vector[i].x;
          RF well_top = setupdata.wdlist.pointdata_vector[i].well_top;
          RF well_bottom = setupdata.wdlist.pointdata_vector[i].well_bottom;
#ifdef DIMENSION3
          RF well_y=setupdata.wdlist.pointdata_vector[i].y;
#endif

          if( global[0] >= well_x-meshsize_2[0] &&  global[0] < well_x+meshsize_2[0] 
#ifdef DIMENSION3             
              && global[1] >= well_y-meshsize_2[1] &&  global[1] < well_y+meshsize_2[1] 
              && global[2] >= well_bottom && global[2] < well_top
#else
              && global[1] >= well_bottom && global[1] < well_top
#endif
              ){
          
            RF pumping_rate, input,injection_time,rho_s, c_s;
            RF porosity(0.3); // will be set in zone_parameters(...);
              
            pumping_rate = setupdata.wdlist.pointdata_vector[i].well_rate;

            if( equationMode == EQ::adjoint )
              pumping_rate = -pumping_rate;

            if( passengerType == PassengerType::heat ){
              input = setupdata.wdlist.pointdata_vector[i].temperature;
              injection_time=setupdata.wdlist.pointdata_vector[i].temperature_injection_time;
#ifdef DIMENSION3                       
              zone_parameters(global[2], porosity, rho_s, c_s);
#else
              zone_parameters(global[1], porosity, rho_s, c_s);
#endif       
            }else{
              input = setupdata.wdlist.pointdata_vector[i].concentration;
              injection_time=setupdata.wdlist.pointdata_vector[i].concentration_injection_time;
            }
         
              
              
            if(injection_time<=GEO_EPSILON)
              injection_time=0.0;
            injection_time*=injection_time*0.5;


            if( pumping_rate > GEO_EPSILON*0.5 && fabs(input) >  GEO_EPSILON*0.5) {
                //UINT flag_source=1; // Ask Ronnie: Wofür?
              RF heat_factor=1.0;
                
              if( passengerType == PassengerType::heat )
                heat_factor=(rho_w*c_w)/( (porosity*rho_w*c_w)+((1.0-porosity)*rho_s*c_s) );
                  
              Dune::GeometryType gt = eg.geometry().type();

                  
              std::vector<Dune::FieldVector<RF,dim>> well_points;
              std::vector<size_t> well_points_cornerindex;
    
              size_t nElementCorners = eg.geometry().corners();
              for( size_t iElementCorner = 0; iElementCorner < nElementCorners; iElementCorner++ ) {
                Dune::FieldVector<RF,dim> elementcorner = eg.geometry().corner( iElementCorner );
#ifdef DIMENSION3
                if( isPointInsideReachOfWell( elementcorner,meshsize_2[0], meshsize_2[1], i) )
#else
                  if( isPointInsideReachOfWell( elementcorner, meshsize_2[0], i ) )
#endif
                    {
                      well_points.push_back( elementcorner );
                      well_points_cornerindex.push_back( iElementCorner );

                    }
              }
                                    
                  
              for(UINT iEdge=0; iEdge<well_points.size(); iEdge++) {
                Dune::FieldVector<RF,dim> local=eg.geometry().local( well_points[iEdge]);
                lfsv.finiteElement().localBasis().evaluateFunction( local, shapefunctions);

                RF factor=eg.geometry().volume();
                RF nPoints=eg.geometry().corners();
                for( UINT ii=0; ii<lfsv.size(); ii++ ) {

                  residual.accumulate( lfsv, ii,  (-pumping_rate/factor*nPoints)*input*injection_time*heat_factor* shapefunctions[ii]);

                }

              }

            }
              
          }
        }

#endif // WELL_FRACTURE_MODEL 
        return true;
      }


      template<typename EG, typename COORDINATES>
      bool evaluate_function( const EG& e
                              , const COORDINATES& xlocal
                              , REAL& fval
                              ) const {
      
        DGF dgf( gfs, vc_source );  
        // For the evaluation of the dgf, a Dune::FieldVector of dim 1 is required:
        Dune::FieldVector<RF, 1> fvec(0.0);
        dgf.evaluate( e, xlocal, fvec );
        fval = fvec[0];
        return true;
      }


      // source of the combi-adjoint:
      // 
      // Unresolved problem: 
      // FEM version of combi-adjoint RHS evaluation
      // works properly only in sequential mode.
      // In parallel mode, the set of leaf elements on the overlap of the coarse
      // grid will not be complete because the refined grid will have only 
      // those child elements on overlap=1.
      // The DG version circumvents this problem by the CCFV discretization of the 
      // RHS sourceterm which requires only values on the cell center :-)
      // These values are being communicated in advance :-))
      // See: void evaluate_residual_on_skeleton(...) used inside "GroundWaterOperatorCCFV.hh"
      template<typename EG, typename COORDINATES>
      bool evaluate_vector( const EG& e
                            , const COORDINATES& xlocal
                            , Dune::FieldVector<RF,dim>& fvec
                            ) const
      {
        const int baselevel = inputdata.domain_data.yasp_baselevel;

        //DGF dgf_1( gfs, vc_source_1 );
        DRGF dgf_1( gfs, vc_source_1 );
        DARCY_FLUX_DGF dgf_2( gwp, gfs, vc_source_2, well_flag );

        // For the evaluation of the dgf_1, a Dune::FieldVector of dim 1 is required:
        Dune::FieldVector<RF, 1> fvec1(0.0);

        //dgf_1.evaluate( e, xlocal, fvec1 );
        dgf_1.evaluate_on_leaf( maxGridLevel, baselevel, 
                                e, xlocal, fvec1 );
      
        // For the evaluation of the dgf_2, a Dune::FieldVector of dim is required:
        Dune::FieldVector<RF, dim> fvec2(0.0);

        //dgf_2.evaluate( e, xlocal, fvec2 );
        dgf_2.evaluate_on_leaf( maxGridLevel, baselevel, 
                                e, xlocal, fvec2 );
      
        fvec = fvec2;
        fvec *= fvec1[0];   // fvec = m0_adj * (-K grad(m0))
        

        if( n_RHS_Functions>2 ){
          //DGF dgf_3( gfs, vc_source_3 );
          DRGF dgf_3( gfs, vc_source_3 );
          DARCY_FLUX_DGF dgf_4( gwp, gfs, vc_source_4, well_flag );
          // For the evaluation of the dgf_3, a Dune::FieldVector of dim 1 is required:
          Dune::FieldVector<RF, 1> fvec3(0.0);
          //dgf_3.evaluate( e, xlocal, fvec3 );
          dgf_3.evaluate_on_leaf( maxGridLevel, baselevel, 
                                  e, xlocal, fvec3 );
      
          // For the evaluation of the dgf_2, a Dune::FieldVector of dim is required:
          Dune::FieldVector<RF, dim> fvec4(0.0);
          //dgf_4.evaluate( e, xlocal, fvec4 );
          dgf_4.evaluate_on_leaf( maxGridLevel, baselevel, 
                                  e, xlocal, fvec4 );
 
          fvec.axpy( fvec3[0], fvec4 );  // fvec = - m0_adj * K grad(m0) - m1_adj * K grad(m1)
        }
        fvec *= -1.0;                 // fvec = + m0_adj * K grad(m0) + m1_adj * K grad(m1)
      
        return true;
      }


      template<typename IG, typename LFSV, typename R>
      void evaluate_on_boundary (const IG& ig, 
                                 const LFSV& lfsv_s, 
                                 R& r_s 
                                 ) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits LocalBasisTraits;
        typedef typename LocalBasisTraits::DomainFieldType DF;
        //typedef typename LocalBasisTraits::RangeFieldType  RF;
        // typedef typename LocalBasisTraits::RangeType RangeType;
        // typedef typename LFSV::Traits::SizeType size_type;
        
        const int dim = IG::dimension;
        Dune::GeometryType gtface = ig.geometryInInside().type();
        
        const Dune::FieldVector<DF,dim-1> 
          face_local = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        typename GWP::Traits::BCType bctype = gwp.bctype(ig.intersection(),face_local);

        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        // evaluate Y-field tensor at cell center, assume it is constant over elements
        Dune::FieldVector<DF,dim> inside_global = ig.inside()->geometry().center();
        const Dune::FieldVector<DF,dim>& localcenter_inside
          = Dune::ReferenceElements<DF,dim>::general( ig.inside()->type() ).position(0,0);

        typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(*(ig.inside()),localcenter_inside));
        RF A_eff = tensor_inside[0][0];

        DRGF dgf_1( gfs, vc_source_1 ); // m0_adj
        DARCY_FLUX_DGF dgf_2( gwp, gfs, vc_source_2, well_flag ); // -K*grad(psi_m0)

        DRGF dgf_3( gfs, vc_source_3 ); // m1_adj
        DARCY_FLUX_DGF dgf_4( gwp, gfs, vc_source_4, well_flag ); // -K*grad(psi_m1)

        const int qorder = 2;
        // select quadrature rule for boundary intergal
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);
        
        if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet){
          // loop over quadrature points and integrate normal flux
          for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it) {
            Dune::FieldVector<DF,dim-1> qPoint_local = it->position(); // local coordinate!

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global( qPoint_local );
            RF factor = it->weight() * ig.geometry().integrationElement( qPoint_local );

            Dune::FieldVector<RF, 1> m0_adj_s(0.0);
            dgf_1.evaluate_on_leaf( maxGridLevel, *(ig.inside()), local, m0_adj_s );
            RF m0_adj = m0_adj_s[0];

            Dune::FieldVector<RF, dim> Kgrad_m0(0.0);
            dgf_2.evaluate_on_leaf( maxGridLevel, *(ig.inside()), local, Kgrad_m0 );
            RF nKgrad_m0 = - (n_F * Kgrad_m0);

            RF contribution_s = - (nKgrad_m0 * m0_adj); // - n * m0_adj * K * gradm0

            if( n_RHS_Functions>2 ){
              Dune::FieldVector<RF, 1> m1_adj_s(0.0);
              dgf_3.evaluate_on_leaf( maxGridLevel, *(ig.inside()), local, m1_adj_s );
              RF m1_adj = m1_adj_s[0];

              Dune::FieldVector<RF, dim> Kgrad_m1(0.0);
              dgf_4.evaluate_on_leaf( maxGridLevel, *(ig.inside()), local, Kgrad_m1 );
              RF nKgrad_m1 = - (n_F * Kgrad_m1);

              contribution_s -= (nKgrad_m1 * m1_adj); // - n * m1_adj * K * gradm1
            }              
            
            contribution_s *= -1.0; // corresponds to the line fvec *= -1.0;
            //std::cout << "contribution_s = " << contribution_s << std::endl;
            r_s.accumulate(lfsv_s,0,contribution_s);

          }

        }

      }



      // See: void evaluate_vector(...) used inside "GroundWaterOperatorGalerkin.hh"
      template<typename IG, typename LFSV, typename R>
      void evaluate_residual_on_skeleton (const IG& ig, 
                                          const LFSV& lfsv_s, const LFSV& lfsv_n, 
                                          R& r_s, R& r_n
                                          ) const {
      
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits LocalBasisTraits;
        typedef typename LocalBasisTraits::DomainFieldType DF;
        //typedef typename LocalBasisTraits::RangeFieldType  RF;
        
        const int dim = IG::dimension;
        
        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,dim> inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,dim> outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();
        
        // evaluate Y-field tensor at cell center, assume it is constant over elements
        const Dune::FieldVector<DF,dim>& localcenter_inside 
          = Dune::ReferenceElements<DF,dim>::general( ig.inside()->type() ).position(0,0);

        // get global coordinate
        //Dune::FieldVector<DF,dim>& globalcenter_inside 
        //  = ig.geometryInInside().global(localcenter_inside);

        const Dune::FieldVector<DF,dim>& localcenter_outside 
          = Dune::ReferenceElements<DF,dim>::general( ig.outside()->type() ).position(0,0);

        typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(*(ig.inside()),localcenter_inside));

        typename GWP::Traits::DiffusionTensorType tensor_outside(gwp.DiffusionTensor(*(ig.outside()),localcenter_outside));


        RF A_s = tensor_inside[0][0];
        RF A_n = tensor_outside[0][0];

        // effective conductivity:
        RF element_volume_s = ig.inside()->geometry().volume();
#ifdef DIMENSION3
        RF element_length_s = std::pow( element_volume_s, 1.0/3.0 );
#else
        RF element_length_s = std::sqrt( element_volume_s );
#endif
        
        RF element_volume_n = ig.outside()->geometry().volume();
#ifdef DIMENSION3
        RF element_length_n = std::pow( element_volume_n, 1.0/3.0 );
#else
        RF element_length_n = std::sqrt( element_volume_n );
#endif
        
        RF A_eff = 0.0;
        if( element_length_s - element_length_n<1E-12 )
          A_eff = Dune::Gesis::General::harmonicAverage(A_s,A_n);
        else
          A_eff = Dune::Gesis::General::harmonicAverageWeightedByDistance( A_s, 
                                                                                  A_n, 
                                                                                  element_length_s, 
                                                                                  element_length_n );
        //std::cout << "A_eff = " << A_eff << std::endl;
        

        // TODO:
        // dgf_1.evaluate_on_leaf() can actually be replaced by
        // dgf_1.evaluate() because we are using the CCFV projected version
        // of the DG solution. That way, the grid dependent 
        // variable 'baselevel' can be removed here.

        DRGF dgf_1( gfs, vc_source_1 ); // m0_adj
        DRGF dgf_2( gfs, vc_source_2 ); // m0
        
        // For the evaluation of a dgf_1, a Dune::FieldVector is required:


        //std::cout << "DEBUG: dgf_1.evaluate_on_leaf(inside), maxGridLevel = " 
        //          << maxGridLevel << std::endl;

#ifdef USE_YASP
        int baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_ALUGRID
        int baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
#ifdef USE_UG
        int baselevel = inputdata.domain_data.ug_baselevel;
#endif

        Dune::FieldVector<RF, 1> m0_adj_s(0.0);
        dgf_1.evaluate_on_leaf( maxGridLevel, baselevel, 
                                *(ig.inside()), localcenter_inside, m0_adj_s );


        //std::cout << "DEBUG: dgf_1.evaluate_on_leaf(outside), maxGridLevel = " 
        // << maxGridLevel << std::endl;

        Dune::FieldVector<RF, 1> m0_adj_n(0.0);
        dgf_1.evaluate_on_leaf( maxGridLevel, baselevel, 
                                *(ig.outside()), localcenter_outside, m0_adj_n );

        RF m0_adj = 0.0;
#ifdef HARMONIC_AVG
        if( element_length_s - element_length_n<1E-12 )
          m0_adj = Dune::Gesis::General::
            harmonicAverage(m0_adj_n[0],m0_adj_s[0]);
        else
          m0_adj = Dune::Gesis::General::
            harmonicAverageWeightedByDistance( m0_adj_s[0], 
                                               m0_adj_n[0], 
                                               element_length_s, 
                                               element_length_n );
#else
        m0_adj = 0.5 * ( m0_adj_n[0] + m0_adj_s[0] );
#endif
        //std::cout << "m0_adj = " << m0_adj << std::endl;

        Dune::FieldVector<RF, 1> m0_s(0.0);
        dgf_2.evaluate_on_leaf( maxGridLevel, baselevel, 
                                *(ig.inside()), localcenter_inside, m0_s );
        Dune::FieldVector<RF, 1> m0_n(0.0);
        dgf_2.evaluate_on_leaf( maxGridLevel, baselevel, 
                                *(ig.outside()), localcenter_outside, m0_n );
        RF Kgrad_m0 = A_eff * ( m0_n[0] - m0_s[0] ) / distance;

        //std::cout << "m0_n[0] - m0_s[0] = " << m0_n[0] - m0_s[0] << std::endl;

        RF face_volume = ig.geometry().volume();
        RF contribution_s = - (Kgrad_m0 * m0_adj) * face_volume;

        if( n_RHS_Functions>2 ) {

          DRGF dgf_3( gfs, vc_source_3 ); // m1_adj
          DRGF dgf_4( gfs, vc_source_4 ); // m1
          Dune::FieldVector<RF, 1> m1_adj_s(0.0);
          dgf_3.evaluate_on_leaf( maxGridLevel, baselevel, 
                                  *(ig.inside()), localcenter_inside, m1_adj_s );
          Dune::FieldVector<RF, 1> m1_adj_n(0.0);
          dgf_3.evaluate_on_leaf( maxGridLevel, baselevel, 
                                  *(ig.outside()), localcenter_outside, m1_adj_n );

          RF m1_adj = 0.0;
#ifdef HARMONIC_AVG
          if( element_length_s - element_length_n<1E-12 )
            m1_adj = Dune::Gesis::General::
              harmonicAverage(m1_adj_n[0],m1_adj_s[0]);
          else
            m1_adj = Dune::Gesis::General::
              harmonicAverageWeightedByDistance( m1_adj_s[0], 
                                                 m1_adj_n[0], 
                                                 element_length_s, 
                                                 element_length_n );
#else
          m1_adj = 0.5 * ( m1_adj_n[0] + m1_adj_s[0] );
#endif
          //std::cout << "m1_adj = " << m1_adj << std::endl;

          Dune::FieldVector<RF, 1> m1_s(0.0);
          dgf_4.evaluate_on_leaf( maxGridLevel, baselevel, 
                                  *(ig.inside()), localcenter_inside, m1_s );
          Dune::FieldVector<RF, 1> m1_n(0.0);
          dgf_4.evaluate_on_leaf( maxGridLevel, baselevel, 
                                  *(ig.outside()), localcenter_outside, m1_n );
          RF Kgrad_m1 = A_eff * ( m1_n[0] - m1_s[0] ) / distance;

          //std::cout << "m1_n[0] - m1_s[0] = " << m1_n[0] - m1_s[0] << std::endl;

          //std::cout << "face_volume = " << face_volume << std::endl;
          contribution_s -= (Kgrad_m1 * m1_adj) * face_volume;
        }
        
        contribution_s *= -1.0; // corresponds to the line fvec *= -1.0;

        //std::cout << "contribution_s = " << contribution_s << std::endl;
        r_s.accumulate(lfsv_s,0,contribution_s);

#ifndef SKELETON_TWOSIDED
        RF contribution_n = - contribution_s;
        r_n.accumulate(lfsv_n,0,contribution_n);
#endif
      } // void evaluate_residual_on_skeleton()

    }; // class FunctionSource


  } // PDELab

} // Dune
#endif
