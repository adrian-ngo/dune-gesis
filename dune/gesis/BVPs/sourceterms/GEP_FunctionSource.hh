#ifndef DUNE_GESIS_GEP_FUNCTION_SOURCE_HH
#define DUNE_GESIS_GEP_FUNCTION_SOURCE_HH

namespace Dune {
  namespace Gesis {

    // Eqn. (66),(147),(153)

    //=========================================================================
    // The template class to define the source term for forward flow simulations
    //========================================================================
    
    template<typename KAPPA_FIELD
             , typename GFS_CG   // eqn (66): m_k
             , typename GFS_GW   // eqn (66): phi0
             , typename IDT
             , typename SDT>
    class GEP_FunctionSource
      : public SourceTermInterface<
      SourceTermTraits<typename GFS_GW::Traits::GridViewType,REAL>,
      GEP_FunctionSource<KAPPA_FIELD,
                         GFS_CG,  // eqn (66): m_k
                         GFS_GW,  // eqn (66): phi0
                         IDT,SDT>
      > 
    {
    private:
      typedef typename GFS_GW::Traits::GridViewType GV_GW;
      typedef typename GFS_CG::Traits::GridViewType GV_TP;
      enum{dim=GV_GW::dimension};

      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;

      const KAPPA_FIELD& kappa_field;
      const GFS_CG& gfs_cg;  // eqn (66): m_k
      const GFS_GW& gfs_gw;  // eqn (66): phi0
      const IDT& inputdata;
      const SDT& setupdata;

      VCType_CG vc_source_1;  // eqn (66): m_k
      VCType_GW vc_source_2;  // eqn (66): phi0

      const int maxGridLevel; // Important if TPE grid is finer than the GWE grid



    public:
    
      typedef REAL RF;
      // refined versions evaluating on leaf elements:
      typedef Dune::Gesis::DiscreteRefinedGridFunction<GFS_CG,VCType_CG> DRGF;
      typedef GradientVectorField<KAPPA_FIELD,GFS_GW> DARCY_FLUX_DGF;

      const SourceNature source_nature;
      //typedef SourceTermTraits<GV,REAL> Traits;

      // constructor:
      GEP_FunctionSource( const KAPPA_FIELD& kappa_field_
                          , const GFS_CG& gfs_cg_  // eqn (66): m_k
                          , const GFS_GW& gfs_gw_  // eqn (66): phi0
                          , const IDT& inputdata_
                          , const SDT& setupdata_
                          , const int maxGridLevel_=0
                          , const bool heat_flag_=false
                          , const bool well_flag_=false
                          )
        :
        kappa_field(kappa_field_),
        gfs_cg( gfs_cg_ ),      // eqn (66): m_k
        gfs_gw( gfs_gw_ ),      // eqn (66): phi0
        inputdata( inputdata_ ),
        setupdata( setupdata_ ),
        vc_source_1( gfs_cg, 0.0 ), // eqn (66): m_k
        vc_source_2( gfs_gw, 0.0 ), // eqn (66): phi0
        maxGridLevel(maxGridLevel_),
        //heat_flag(heat_flag_),
        //well_flag(well_flag_),
        source_nature( GEP_FUNCTIONAL_SOURCE )
      {
        logger << "GEP_FunctionSource constructor with maxGridLevel = " 
               << maxGridLevel
               << std::endl;
      }

      void reset_rhs()
      {
        vc_source_1 = VCType_CG(gfs_cg,0.0);
        vc_source_2 = VCType_GW(gfs_gw,0.0);
      }


      // for equation (66):
      void set_rhs( const VCType_CG& vc_1,
                    const VCType_GW& vc_2 )
      {
        reset_rhs();
        vc_source_1 = vc_1; // eqn (66): m_k
        vc_source_2 = vc_2; // eqn (66): phi0
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
        // For equation (66):
        // This function is used only for the Standard Galerkin Groundwater Operator!
        
        const int baselevel = inputdata.domain_data.yasp_baselevel;

        DRGF dgf_1( gfs_cg, vc_source_1 );
        DARCY_FLUX_DGF dgf_2( kappa_field, gfs_gw, vc_source_2 /*, well_flag*/ );

        // For the evaluation of the dgf_1, a Dune::FieldVector of dim 1 is required:
        Dune::FieldVector<RF,1> fvec1(0.0);

        dgf_1.evaluate_on_leaf( maxGridLevel, baselevel, 
                                e, xlocal, fvec1 );
      
        // For the evaluation of the dgf_2, a Dune::FieldVector of dim is required:
        Dune::FieldVector<RF, dim> fvec2(0.0);

        //dgf_2.evaluate( e, xlocal, fvec2 );
        dgf_2.evaluate_on_leaf( maxGridLevel, baselevel, 
                                e, xlocal, fvec2 );
      
        fvec = fvec2;
        fvec *= fvec1[0];


        fvec *= -1.0;                 // fvec = - m_k * kappa * grad phi0
      
        return true;
      }


      template<typename IG, typename LFSV, typename R>
      void evaluate_on_boundary (const IG& ig, 
                                 const LFSV& lfsv_s, 
                                 R& r_s 
                                 ) const
      {
        // RHS of equation (66), boundary part
        // CCFV version!

        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits LocalBasisTraits;
        typedef typename LocalBasisTraits::DomainFieldType DF;
        //typedef typename LocalBasisTraits::RangeFieldType  RF;
        //typedef typename LocalBasisTraits::RangeType RangeType;
        //typedef typename LFSV::Traits::SizeType size_type;
        
        const int dim = IG::dimension;
        Dune::GeometryType gtface = ig.geometryInInside().type();
        
        const Dune::FieldVector<DF,dim-1> 
          face_local = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        typename KAPPA_FIELD::Traits::BCType bctype = kappa_field.bctype(ig.intersection(),face_local);

        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        // evaluate Y-field tensor at cell center, assume it is constant over elements
        Dune::FieldVector<DF,dim> inside_global = ig.inside()->geometry().center();
        const Dune::FieldVector<DF,dim>& localcenter_inside
          = Dune::ReferenceElements<DF,dim>::general( ig.inside()->type() ).position(0,0);

        typename KAPPA_FIELD::Traits::DiffusionTensorType tensor_inside(kappa_field.DiffusionTensor(*(ig.inside()),localcenter_inside));
        // RF A_eff = tensor_inside[0][0];

        DRGF dgf_1( gfs_cg, vc_source_1 ); // m0_adj
        DARCY_FLUX_DGF dgf_2( kappa_field, gfs_gw, vc_source_2 /*, well_flag*/ ); // -K*grad(psi_m0)


        const int qorder = 2;
        // select quadrature rule for boundary intergal
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);
        
        if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet){
          // loop over quadrature points and integrate normal flux
          for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it) {
            Dune::FieldVector<DF,dim-1> qPoint_local = it->position(); // local coordinate!

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global( qPoint_local );
            // RF factor = it->weight() * ig.geometry().integrationElement( qPoint_local );

            Dune::FieldVector<RF, 1> m0_adj_s(0.0);
            dgf_1.evaluate_on_leaf( maxGridLevel, *(ig.inside()), local, m0_adj_s );
            RF m0_adj = m0_adj_s[0];

            Dune::FieldVector<RF, dim> Kgrad_m0(0.0);
            dgf_2.evaluate_on_leaf( maxGridLevel, *(ig.inside()), local, Kgrad_m0 );
            RF nKgrad_m0 = - (n_F * Kgrad_m0);

            RF contribution_s = - (nKgrad_m0 * m0_adj); // - n * m0_adj * K * gradm0
            
            contribution_s *= -1.0; // corresponds to the line fvec *= -1.0;
            //std::cout << "contribution_s = " << contribution_s << std::endl;
            r_s.accumulate(lfsv_s,0,contribution_s);

          }

        }

      }


      //
      // Compare: void evaluate_vector(...) used inside "GroundWaterOperatorGalerkin.hh"
      //
      template<typename IG, typename LFSV, typename R>
      void evaluate_residual_on_skeleton (const IG& ig, 
                                          const LFSV& lfsv_s, const LFSV& lfsv_n, 
                                          R& r_s, R& r_n
                                          ) const {
        // RHS of equation (66), skeleton part, used by "GroundWaterOperatorCCFV.hh"
        // CCFV version!
      
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

        typename KAPPA_FIELD::Traits::DiffusionTensorType tensor_inside(kappa_field.DiffusionTensor(*(ig.inside()),localcenter_inside));

        typename KAPPA_FIELD::Traits::DiffusionTensorType tensor_outside(kappa_field.DiffusionTensor(*(ig.outside()),localcenter_outside));


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
        
        DRGF dgf_1( gfs_cg, vc_source_1 ); // m_k in equation (66)
        DRGF dgf_2( gfs_gw, vc_source_2 ); // phi0 in equation (66)
        
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


        Dune::FieldVector<RF,1> m0_adj_s(0.0);
        dgf_1.evaluate_on_leaf( maxGridLevel, baselevel, 
                                *(ig.inside()), localcenter_inside, m0_adj_s );


        //std::cout << "DEBUG: dgf_1.evaluate_on_leaf(outside), maxGridLevel = " 
        // << maxGridLevel << std::endl;

        Dune::FieldVector<RF,1> m0_adj_n(0.0);
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
        
        contribution_s *= -1.0; // corresponds to the line fvec *= -1.0;

        //std::cout << "contribution_s = " << contribution_s << std::endl;
        r_s.accumulate(lfsv_s,0,contribution_s);

#ifndef SKELETON_TWOSIDED
        RF contribution_n = - contribution_s;
        r_n.accumulate(lfsv_n,0,contribution_n);
#endif
      } // void evaluate_residual_on_skeleton()




      template<typename EG, typename COORDINATES>
      bool evaluate_function( const EG& e
                              , const COORDINATES& xlocal
                              , REAL& fval
                              ) const {
        // Hint: 
        // The DGF is being created out of the VCType everytime.
        // I don't know yet whether this will have a negative effect on the performance.
        // I hope that this construtor does not do much work here!

#ifdef USE_YASP
        int baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_ALUGRID
        int baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
#ifdef USE_UG
        int baselevel = inputdata.domain_data.ug_baselevel;
#endif

        DARCY_FLUX_DGF dgf_1( kappa_field, // dummy place-holder!
                              gfs_gw, 
                              vc_source_1, 
                              baselevel,
                              true //compute gradient only
                              );
        DARCY_FLUX_DGF dgf_2( kappa_field, gfs_gw, vc_source_2 /*, well_flag*/ ); // -K*grad(psi_m0)

        // For the evaluation of the dgf, a Dune::FieldVector of dim 1 is required:
        Dune::FieldVector<RF, dim> fvec1(0.0);
        dgf_1.evaluate( e, xlocal, fvec1 );
		
        // For the evaluation of the dgf_2, a Dune::FieldVector of dim is required:
        Dune::FieldVector<RF, dim> fvec2(0.0);
        dgf_2.evaluate( e, xlocal, fvec2 );

        //logger<<"flux: "<<fvec2[0]<<","<<fvec2[1]<<std::endl;
		
        //scalar product
        fval = (fvec1*fvec2);	
                
        return true;
      }  // bool evaluate_function()


    }; // class GEP_FunctionSource


  } // PDELab

} // Dune
#endif // DUNE_GESIS_GEP_Function_Source_HH
