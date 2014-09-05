// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef TRANSPORT_OPERATOR_SDFEM_HH
#define TRANSPORT_OPERATOR_SDFEM_HH

//
// A local operator for solving ...
//
//
// ... the stationary transport equation with
// * velocity-field beta
// * Scheidegger Dispersion tensor D and porosity theta
// * and sources/sinks f
//
// - \beta \nabla u + \nabla \cdot ( \theta D \nabla u ) = f in \Omega
//                                                     u = g on \Gamma_D
//                    { -\theta D \nabla u } \cdot \nu   = j on \Gamma_N
//
//
//
//
// special case (bAdjoint==true): 
// ===================================
// beta := -beta 
// ... the adjoint transport problem with reversed velocity field
//

#include<vector>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include"dune/pdelab/common/geometrywrapper.hh"
#include"dune/pdelab/common/function.hh"
#include<dune/pdelab/gridoperator/gridoperator.hh>
//#include"dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh"
#include"dune/pdelab/localoperator/pattern.hh"
#include"dune/pdelab/localoperator/flags.hh"
#include"dune/pdelab/localoperator/idefault.hh"
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include "deltaSD.hh"

#include "dune/gesis/common/eval.hh"
#include "dune/gesis/BVPs/wells/wellposition.hh"


#ifndef CONVECTIVE_FORMULATION
#define DO_ALPHA_BOUNDARY
#endif

namespace Dune {
  namespace GeoInversion {

    template< 
      typename TP,
      typename FEM,
      typename BCType,
      typename SOURCE_TYPE,
      typename IDT,
      typename SDT
      >
	class StreamlineDiffusionOperator 
      :
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags
	{

    private:
      const IDT& inputdata;
      const SDT& setupdata;
      const TP& tp;

      const BCType& bctype;
      SOURCE_TYPE& sourceterm;

      const Passenger::Type passenger;
      const EQ::Mode equationMode;

      template< typename DFV0   // Dune::FieldVector of codim 0
                , typename RF
                >
      bool isPointInsideReachOfWell(
                                    const DFV0& elementpoint, 
                                    const RF& reach,
#ifdef DIMENSION3
                                    const RF& reach_y,
#endif
                                    const int& iWell
                                    ) const
      {

        // extract well-data:
        RF well_position_x =  setupdata.wdlist.pointdata_vector[iWell].x;
#ifdef DIMENSION3
        RF well_position_y =  setupdata.wdlist.pointdata_vector[iWell].y;
#endif
        RF well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
        RF well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;

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



      typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;

      typedef Dune::PDELab::LocalBasisCache<LocalBasisType> Cache;

      // In theory it is possible that one and the same local operator is
      // called first with a finite element of one type and later with a
      // finite element of another type.  Since finite elements of different
      // type will usually produce different results for the same local
      // coordinate they cannot share a cache.  Here we use a vector of caches
      // to allow for different orders of the shape functions, which should be
      // enough to support p-adaptivity.  (Another likely candidate would be
      // differing geometry types, i.e. hybrid meshes.)

      std::vector<Cache> cache;



	public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = false }; // would require public FullSkeletonPattern
      enum { doPatternBoundary = false }; // would require public FullBoundaryPattern which does not exist! Why not?

	  // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doAlphaSkeleton = false };

#ifdef DO_ALPHA_BOUNDARY
      enum { doAlphaBoundary = true };
#else
      enum { doAlphaBoundary = false };
#endif

      enum { doLambdaVolume = true };
      enum { doLambdaSkeleton = false };
      enum { doLambdaBoundary = true };



      template<typename VCType>
      void set_rhs( const VCType& vc_solution )
      {
        sourceterm.set_rhs( vc_solution );
      }
      
      template<typename VType,typename UType>
      void set_rhs( const VType& vc_1,const UType& vc_2 )
      {
        sourceterm.set_rhs( vc_1, vc_2 );
      }

      template<typename COORDINATES>
      void set_PointSourceLocation( const COORDINATES& x )
      {
        sourceterm.set_PointSourceLocation( x );
      }


      // Constructor:
      StreamlineDiffusionOperator(
                                  const IDT& inputdata_,
                                  const SDT& setupdata_,
                                  const TP& tp_,
                                  const BCType& bctype_,
                                  SOURCE_TYPE& sourceterm_,
                                  const Passenger::Type passenger_,
                                  const EQ::Mode equationMode_
                                  )
        : 
        inputdata( inputdata_ ),
        setupdata( setupdata_ ),
        tp(tp_), 
        bctype( bctype_ ),
        sourceterm( sourceterm_ ),
        passenger( passenger_ ),
        equationMode( equationMode_ ),
        cache(20)
      {
        logger << "StreamlineDiffusionOperator constructor" << std::endl;
      }
      
      





     

      template<
        typename EG,
        typename LFSU,
        typename X,
        typename M
        >
      void addPointSourceToJacobianVolumeMatrix(
                                        const int& iWell, 
                                        const EG& eg, 
                                        const LFSU& lfsu,
                                        const X& x, 
                                        M& mat 
                                        ) const
      {

        return; // testwise

		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

		
		// dimensions
        const int dim = EG::Geometry::dimension;
	
        RF pumping_rate = setupdata.pdlist.pointdata_vector[iWell].value;
        if( equationMode == EQ::adjoint )
          pumping_rate=-pumping_rate;

        if(  pumping_rate  > GEO_EPSILON*0.5 ){
          Dune::FieldVector<RF,dim> pointsource_global;
          pointsource_global[0] = setupdata.pdlist.pointdata_vector[iWell].x;
          pointsource_global[1] = setupdata.pdlist.pointdata_vector[iWell].y;
#ifdef DIMENSION3
          pointsource_global[2] = setupdata.pdlist.pointdata_vector[iWell].z;
#endif
          Dune::FieldVector<RF,dim> pointsource_local = eg.geometry().local( pointsource_global );

          UINT iOutside = 0;
          for( UINT i=0; i<dim; i++ ){
            if( (pointsource_local[i] < 0.0) || (pointsource_local[i] >= 1.0) ){
              iOutside++;
            }
          }

          if( iOutside==0 ) {
            std::vector<RangeType> shapefunctions( lfsu.size() );
            evalFunctionBasis( pointsource_local, lfsu, cache, shapefunctions );

            RF factor=eg.geometry().volume();
            RF nPoints=eg.geometry().corners();
            
            for( UINT i=0; i<lfsu.size(); i++ ) {
              for( UINT j=0; j<lfsu.size(); j++ ) {
                mat.accumulate( lfsu, i,lfsu, j, ( pumping_rate/factor*nPoints)* shapefunctions[i]*shapefunctions[j]);
              }
            }
          }
        }
      }
      



      template<typename EG,
               typename LFSU,
               typename X,
               typename M>
      void addWellTo_Jacobian_Volume( const int& iWell, 
                                      const EG& eg, 
                                      const LFSU& lfsu,
                                      const X& x, 
                                      M& mat 
                                      ) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSU::Traits::SizeType size_type;
		
        //Dune::GeometryType gt = eg.geometry().type();

        RF well_rate = setupdata.wdlist.pointdata_vector[iWell].well_rate;
        if( equationMode == EQ::adjoint )
          well_rate *= -1;
          
        if( well_rate > GEO_EPSILON*0.5 ) {
        
          if( w_outside != isElementWithinWellZone( eg, iWell, inputdata, setupdata ) ){

            const int dim = EG::Geometry::dimension;
        
            RF well_top = setupdata.wdlist.pointdata_vector[iWell].well_top;
            RF well_bottom = setupdata.wdlist.pointdata_vector[iWell].well_bottom;

            const int order = lfsu.finiteElement().localBasis().order();
            const int qorder = 2 * order;

#ifdef USE_YASP
            UINT baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_ALUGRID
            UINT baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
#ifdef USE_UG
            UINT baselevel = inputdata.domain_data.ug_baselevel;
#endif
            REAL baselevel_factor = std::pow(2.0,baselevel);
            REAL level_diff = eg.entity().level() - baselevel;
            REAL refinement_factor_dim = std::pow(2.0,level_diff*(dim-1)); // dim is important here: one level =  four cells in 2D, or 8 cells in 3D

            REAL refinement_factor = std::pow(2.0,level_diff);
            const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1] / refinement_factor;

            //Dune::FieldVector<REAL,dim> cellcenter = eg.geometry().center();
            //REAL z=cellcenter[dim-1];
            REAL screening_length = well_top - well_bottom;          
            REAL screening_fraction = dz / screening_length;

            // element integral along well elements:
            Dune::GeometryType gt = eg.geometry().type();
            const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

            // loop over quadrature points
            for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it) {

              // evaluate shape functions 
              std::vector<RangeType> phi(lfsu.size());
              evalFunctionBasis( it->position(), lfsu, cache, phi );

              RF u0 = 1;

              REAL fval = -well_rate * u0 * screening_fraction / eg.geometry().volume(); // good

              // Take refinement into account:
              if( eg.entity().level() > baselevel ){
                fval *= (refinement_factor_dim / baselevel_factor);
              }
            
              RF factor = it->weight() * eg.geometry().integrationElement( it->position() );

              for (size_type j=0; j<lfsu.size(); j++){
                for (size_type i=0; i<lfsu.size(); i++){
                  mat.accumulate( lfsu, i,lfsu, j, -fval*phi[i]*phi[j]*factor );
                  // addwellto_: jacobian_volume (element integral)
                }
              }

            } // end of loop over quadrature points

          } // end if w_outside

        } // end if well_rate

      }


	  // volume integral depending on test and ansatz functions
	  template< typename EG
                , typename LFSU
                , typename X
                , typename LFSV
                , typename R >
	  void alpha_volume( const EG& eg,
                         const LFSU& lfsu,
                         const X& x,
                         const LFSV& lfsv,
                         R& r
                         ) const
	  {

		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        Dune::GeometryType gt = eg.geometry().type();

        const int qorder = 2 * lfsu.finiteElement().localBasis().order();
        // select quadrature rule
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule( gt, qorder );

        // loop over quadrature points
        for( typename Dune::QuadratureRule<DF,dim>::const_iterator it = rule.begin()
               ; it != rule.end()
               ; ++it )
          {            
            // evaluate velocity at it->position(), because this function expects local coordinates!
            typename TP::Traits::RangeType beta( tp.velocity( eg.entity(), it->position(), equationMode ) );

            // the SDFEM delta:
            typename TP::Traits::RangeFieldType deltaSD
              = tp.deltaSDFEM( eg.entity(), it->position() );

            // evaluate shape functions for lfsu
            std::vector<RangeType> phi( lfsu.size() );
            evalFunctionBasis( it->position(), lfsu, cache, phi );

            // get gradient basis: we assume lfsu = lfsv so that we can use gradphi_u = gradphi_v identically
            std::vector<Dune::FieldVector<RF,dim> > gradphi( lfsu.size() );
            evalGradientBasis( eg,
                               it->position(),
                               lfsu,
                               cache,
                               gradphi );

            // compute gradient of u on the global element
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (unsigned int i=0; i<lfsu.size(); i++) {
              gradu.axpy( x(lfsu,i), gradphi[i] );  // w.axpy( a, y ) means: w = w + a*y
            }

            // factor from the integration element:
            RF factor = it->weight() * eg.geometry().integrationElement( it->position() );

            Dune::FieldVector<RF,dim> Dgradu(0.0);

            // get the Scheidegger dispersion tensor
            typename TP::Traits::DispersionTensorType
              DispersionTensor = tp.DispersionTensor( eg.entity(), it->position(), equationMode );
            DispersionTensor.mv( gradu, Dgradu );  // Dgradu = DispersionTensor * gradu

            //============ assemble diffusion term  ============
            //
            // + Dgradu * gradv  (standard Galerkin)
            //
            //
            for( unsigned int i=0; i<lfsv.size(); i++ ) {
              r.accumulate( lfsv, i, ( Dgradu * gradphi[i] ) * factor );
            }

#ifdef CONVECTIVE_FORMULATION
            //============ assemble convection term ============
            //
            // + beta * gradu * v  ( standard Galerkin, convective formulation )
            //
            for( unsigned int i=0; i<lfsv.size(); i++ ) {
              r.accumulate( lfsv, i, ( beta * gradu ) * phi[i] * factor );
            }
#else

            //============ assemble convection term ============
            //
            // - u * beta * gradv ( standard Galerkin, conservative formulation )
            // 
            // compute u
            RF uval( 0.0 );
            for( size_type i=0; i<lfsu.size(); i++)
              uval += x(lfsu,i)*phi[i];

            for( size_type i=0; i<lfsv.size(); i++ ) {
              r.accumulate( lfsv, i, -uval * ( beta * gradphi[i] ) * factor );
            }
#endif

            //============ assemble SDFEM term  ============
            // adding some artificial diffusion in streamline direction
            //
            // + deltaSD * ( beta * gradu ) * ( beta * gradv )
            //
            for( unsigned int i=0; i<lfsv.size(); i++ )
              r.accumulate( lfsv, i, deltaSD * ( beta * gradu ) * ( beta * gradphi[i] ) * factor );


          } // end loop over quadrature points
          

	  } // alpha_volume()






      template<typename EG, typename LFSV, typename R>
      void addWellTo_Lambda_Volume( const UINT iWell, 
                                    const EG& eg, 
                                    const LFSV& lfsv,
                                    R& residual
                                    ) const {

        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;

        const int dim = EG::Geometry::dimension;

        RF well_rate = setupdata.wdlist.pointdata_vector[iWell].well_rate;
        if( well_rate < 1e-12 )
          return; // no inflow contribution

        RF concentration = 0;
        if( passenger == Passenger::heat )
          concentration = setupdata.wdlist.pointdata_vector[iWell].temperature;
        else
          concentration = setupdata.wdlist.pointdata_vector[iWell].concentration;

        if( equationMode == EQ::adjoint )
          concentration = 0;  // adjoint case: extraction well becomes injection well with zero concentration/temperature!

        if(concentration<1e-12)
          return; // no inflow contribution

        RF injectiontime = 0;
        if( passenger == Passenger::heat )
          injectiontime = setupdata.wdlist.pointdata_vector[iWell].temperature_injection_time;
        else
          injectiontime = setupdata.wdlist.pointdata_vector[iWell].concentration_injection_time;

        if(injectiontime<1e-12)
          return; // no inflow contribution
        
        if( w_outside != isElementWithinWellZone( eg, iWell, inputdata, setupdata ) ){


#ifdef USE_YASP
          UINT baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_ALUGRID
            UINT baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
#ifdef USE_UG
            UINT baselevel = inputdata.domain_data.ug_baselevel;
#endif
          REAL baselevel_factor = std::pow(2.0,baselevel);
          REAL level_diff = eg.entity().level() - baselevel;
          REAL refinement_factor_dim = std::pow(2.0,level_diff*(dim-1)); // dim is important here: one level =  four cells in 2D, or 8 cells in 3D

          REAL refinement_factor = std::pow(2.0,level_diff);
          const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1] / refinement_factor;
          
          REAL well_top           = setupdata.wdlist.pointdata_vector[iWell].well_top;
          REAL well_bottom        = setupdata.wdlist.pointdata_vector[iWell].well_bottom;
          REAL screening_length   = well_top - well_bottom;
          REAL screening_fraction = dz / screening_length;

          // element integral along well elements:
          Dune::GeometryType gt = eg.geometry().type();
          const int qorder = 2 * lfsv.finiteElement().localBasis().order();
          const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

          // loop over quadrature points
          for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it) {

            // evaluate shape functions
            std::vector<RangeType> phi(lfsv.size());
            evalFunctionBasis( it->position(), lfsv, cache, phi );

            RF u0 = concentration * injectiontime;
            // scaling with 1./eg.geometry().volume() is important due to refinement
            REAL fval = well_rate * u0 / eg.geometry().volume() * screening_fraction;     // good! 


            // Take refinement into account:
            if( eg.entity().level() > baselevel ){
              fval *= (refinement_factor_dim / baselevel_factor);
            }

            RF factor = it->weight() * eg.geometry().integrationElement( it->position() );

            for (size_type i=0; i<lfsv.size(); i++){
              residual.accumulate(lfsv,i, -fval*phi[i]*factor);
              // addwellto_: lambda_volume (element integral)
            }
          }
        }
      }








 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        UINT flag_source=0;

        // ============= Source term for adjoint flow problem with point source  ==============
        // source term = delta( x - x_j )
        //
        // or:
        //
        // ============= Source term for forward flow problem with pumping source  ==============
        // source term is some pumping source specified in the InputFile-XML (this case is not used here so far)
        //

        // TODO: Adding Point-sources only! Check with DG implementation!
        if( sourceterm.source_nature == POINT_SOURCE
            ||
            sourceterm.source_nature == SOLUTE_PUMPING_SOURCES
            ) {
          std::vector<RangeType> shapefunctions( lfsv.size() );
          sourceterm.evaluate_residual_on_element( eg.entity(),
                                                   lfsv,
                                                   shapefunctions,
                                                   r );
        }
        


        // ============= Source term for adjoint flow problem wrt 
        //
        // theta*m0 
        //
        // or: 
        //
        // m1_adj 
        //
        
        // forward eqn: (63)
        // adjoint eqn: (150)
        if( sourceterm.source_nature == FUNCTIONAL_SOURCE ||
            sourceterm.source_nature == GEP_FUNCTIONAL_SOURCE )
          {
            const int dimw = EG::Geometry::dimensionworld;
            Dune::GeometryType gt = eg.geometry().type();
            
            // get the local meshsize h:
            RF meshsize = sqrt( eg.geometry().volume() );
            
            //
            // TODO: Adding Point-sources only! Check with DG implementation!
            //if( sourceterm.source_nature == FUNCTIONAL_SOURCE ){
            //  std::vector<RangeType> shapefunctions( lfsv.size() );
            //  sourceterm.evaluate_residual_on_element( eg.entity(),
            //                                           lfsv,
            //                                           shapefunctions,
            //                                           r );
            //}

            if(flag_source==0){
              //logger<<"TP lambda_volume: "<<std::endl;
              // select quadrature rule
              const int qorder = 2 * lfsv.finiteElement().localBasis().order();

              const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
              // loop over quadrature points
              for( typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin()
                     ; it!=rule.end()
                     ; ++it )
                {
                  // evaluate velocity at it->position(), because this function expects local coordinates!
                  typename TP::Traits::RangeType beta( tp.velocity( eg.entity(), it->position(), equationMode ) );

                  // the SDFEM delta:
                  typename TP::Traits::RangeFieldType deltaSD
                    = tp.deltaSDFEM( eg.entity(), it->position() );

                
                  // evaluate shape functions 
                  std::vector<RangeType> phi;
                  evalFunctionBasis( it->position(), lfsv, cache, phi );
                
                  // evaluate right hand side parameter function
                  // evaluate source term
                  typename TP::Traits::RangeFieldType fval = 0.0;

                  // forward eqn: (63) for k=1
                  // adjoint eqn: (150) for k=1
                  sourceterm.evaluate_function( eg.entity(),
                                                it->position(),
                                                fval );
                
                  //logger<<"sourceterm value: " << fval << "    at: "<< it->position()<<std::endl;

                  // integrate f
                  RF factor = it->weight() * eg.geometry().integrationElement( it->position() );
                
                  // ============= source term (standard Galerkin) ==============
                  //
                  // - f * v 
                  //
                  for( unsigned int i=0; i<lfsv.size(); i++ )
                    r.accumulate( lfsv, i, - fval * phi[i] * factor );
                
                  // get gradient basis: we assume lfsu = lfsv so that we can use gradphi_u = gradphi_v identically
                  std::vector<Dune::FieldVector<RF,dim> > gradphi( lfsv.size() );
                  evalGradientBasis( eg,
                                     it->position(),
                                     lfsv,
                                     cache,
                                     gradphi );

                  // ============= source term with SDFEM ==============
                  //
                  // - f * deltaSD * beta * gradv
                  //
                  for (unsigned int i=0; i<lfsv.size(); i++)
                    r.accumulate( lfsv, i, - deltaSD * fval * ( beta * gradphi[i] ) * factor );
                
                } // End of Quadrature Loop
              
            }
            
            
          } // end if( sourceterm.source_nature )


        // Add inflow or outlow rates for all the wells
        UINT nWells = setupdata.wdlist.total;
        for( UINT iWell=0; iWell<nWells; iWell++ ) {
          addWellTo_Lambda_Volume( iWell, 
                                   eg, 
                                   lfsv,
                                   r );
        }
        
      }




	  // boundary integral depending on test and ansatz functions
      // Neumann boundary
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, 
                           const X& x_s, 
                           const LFSV& lfsv_s,
                           R& r_s) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;


        //		typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;
    

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        Dune::GeometryType gtface = ig.geometryInInside().type();

        //Dune::FieldVector<DF,dim-1> facecenterlocal = Dune::GenericReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        //Dune::FieldVector<DF,dim> facecenterinelement = ig.geometryInInside().global( facecenterlocal );

        typedef typename IG::EntityPointer CellEntityPointer;
        const CellEntityPointer eg = ig.inside();

        // evaluate velocity
        //typename TP::Traits::RangeType beta( tp.velocity( *eg, facecenterinelement ) );

        const int qorder = 2 * lfsu_s.finiteElement().localBasis().order();
        // select quadrature rule
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // loop over quadrature points and integrate normal flux
        for( typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin()
               ; it!=rule.end()
               ; ++it )
          {
            // skip rest if we are on Dirichlet boundary
            if( bctype.isDirichlet( ig, it->position() ) )
              continue;

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global( it->position() );

            // evaluate velocity
            typename TP::Traits::RangeType beta( tp.velocity( *eg, local, equationMode ) );

            RF beta_n = beta * ig.unitOuterNormal( it->position() );
            
            if( beta_n < 1E-12 )
              continue; // skip if we are not on outflow

            // We are on an outflow boundary now:

            // evaluate test shape functions 
            std::vector<RangeType> phi(lfsu_s.size());
            evalFunctionBasis( local, lfsu_s, cache, phi );
            
            // evaluate u
            RF uval = 0.0;
            for (int i=0; i<lfsu_s.size(); i++)
              uval += x_s[i] * phi[i];
            
            // integrate J
            RF factor = it->weight()*ig.geometry().integrationElement( it->position() );

            // ========== part of the convection term ====================
            // conservative formulation has boundary integral with term
            //
            // + v * u * beta * n
            //
            for (int i=0; i<lfsv_s.size(); i++)
              r_s[i] += phi[i] * uval * beta_n * factor;
          }
      }




      // boundary integral independent of ansatz functions
      // Neumann boundary
 	  template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = IG::dimension;
        //const int dimw = IG::dimensionworld;

        const int qorder = 2 * lfsv.finiteElement().localBasis().order();
        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule( gtface, qorder );

        // loop over quadrature points and integrate normal flux
        for( typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin()
               ; it!=rule.end()
               ; ++it )
          {
            // skip rest if we are on Dirichlet boundary
            if( bctype.isDirichlet( ig, it->position() ) )
              continue;

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global( it->position() );

            // evaluate test shape functions 
            std::vector<RangeType> phi;
            evalFunctionBasis( local, lfsv, cache, phi );
            
            // evaluate flux boundary condition
            typename TP::Traits::RangeFieldType jflux = 0.0; // tp.j( *(ig.inside()), local );
            
            // integrate J
            RF factor = it->weight() * ig.geometry().integrationElement( it->position() );

            //============= Flux term through boundary ===============
            //
            // - j * v
            //
            for (unsigned int i=0; i<lfsv.size(); i++)
              r.accumulate( lfsv, i, - jflux * phi[i] * factor );
          }

      } // lambda_boundary()






      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
	  void jacobian_volume( const EG& eg, 
							const LFSU& lfsu, 
							const X& x, 
							const LFSV& lfsv, 
                            M& mat
							) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const unsigned int dim = EG::Geometry::dimension;
        const unsigned int dimw = EG::Geometry::dimensionworld;

        Dune::GeometryType gt = eg.geometry().type();

        // get the local meshsize h:
        RF meshsize = sqrt( eg.geometry().volume() );

        const int qorder = 2 * lfsu.finiteElement().localBasis().order();
        // select quadrature rule
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for( typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin()
               ; it!=rule.end()
               ; ++it )
          {
            // evaluate velocity at it->position(), because this function expects local coordinates!
            typename TP::Traits::RangeType beta( tp.velocity( eg.entity(), it->position(), equationMode ) );

            // depth is needed to find the porosity of the right zone!
            Dune::FieldVector<DF,dim> qPoint_global = eg.geometry().global( it->position() );
            CTYPE depth = qPoint_global[dim-1];



            // the SDFEM delta:
            typename TP::Traits::RangeFieldType deltaSD
              = tp.deltaSDFEM( eg.entity(), it->position() );

            // evaluate shape functions (we assume Galerkin method lfsu = lfsv)
            std::vector<RangeType> phi( lfsu.size() );
            evalFunctionBasis( it->position(), lfsu, cache, phi );

            // get gradient basis
            std::vector<Dune::FieldVector<RF,dim> > gradphi( lfsu.size() );
            evalGradientBasis( eg,
                               it->position(),
                               lfsu,
                               cache,
                               gradphi );

            // factor from the integration element:
            RF factor = it->weight() * eg.geometry().integrationElement( it->position() );

            // get the Scheidegger dispersion tensor
            typename TP::Traits::DispersionTensorType
              DispersionTensor = tp.DispersionTensor( eg.entity(), it->position(), equationMode );

            // compute D * gradient of shape functions
            std::vector<Dune::FieldVector<RF,dim> > Dgradphi( lfsu.size() );
            for( unsigned int i=0; i<lfsu.size(); i++ ) {
              DispersionTensor.mv( gradphi[i], Dgradphi[i] );  
              // Dgradu = DispersionTensor * gradu
            }


            // Jacobian-term = 
            //   D * gradphi_j * gradphi_i 
            // - phi_j * (beta * gradphi_i)
            // + deltaSD * (beta * gradphi_j) * (beta * gradphi_i)
            // 
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
				{
                  mat.accumulate( lfsu, i, lfsu, j, ( Dgradphi[j] * gradphi[i] ) * factor );
#ifdef CONVECTIVE_FORMULATION
				  mat.accumulate( lfsu, i, lfsu, j, ( beta * gradphi[j] ) * phi[i] * factor );
#else
				  mat.accumulate( lfsu, i, lfsu, j, - phi[j] * ( beta * gradphi[i] ) * factor );
#endif
				  mat.accumulate( lfsu, i, lfsu, j, ( deltaSD * (beta * gradphi[j]) * (beta * gradphi[i]) ) * factor );
				}

          } // loop over quadrature points



        // Ask Ronnie: When will we need this point-source?
        /*
        if( equationMode == EQ::forward && 
            ( sourceterm.source_nature == SOLUTE_PUMPING_SOURCES ||
              sourceterm.source_nature == FUNCTIONAL_SOURCE ) 
            ) {

          UINT nWells = setupdata.wdlist.total;
         
          for( UINT iWell=0; iWell<nWells; iWell++ ) {
            addPointSourceToJacobianVolumeMatrix(
                                                 iWell, 
                                                 eg, 
                                                 lfsu,
                                                 x, 
                                                 mat
                                                 );
          }
        }
        */
          
        
        // Adding well contribution to the stiffness matrix:
        // Note:
        // For well monitoring purposes, the concentration for the TPE should be zero. This is set in the input-XML.
        UINT nWells = setupdata.wdlist.total;
        for( UINT iWell=0; iWell<nWells; iWell++ ) {
          addWellTo_Jacobian_Volume( iWell, 
                                     eg, 
                                     lfsu,
                                     x, 
                                     mat );
        }

      }


      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, 
							  const X& x, 
							  const LFSV& lfsv_s,
                              M& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        
        //typedef typename LFSV::Traits::FiniteElementType::
        //Traits::LocalBasisType::Traits::JacobianType JacobianType;

        // dimensions
        const unsigned int dim = IG::dimension;
        const unsigned int dimw = IG::dimensionworld;

        typedef typename IG::EntityPointer CellEntityPointer;
        const CellEntityPointer eg = ig.inside();

        Dune::GeometryType gtface = ig.geometryInInside().type();

        const int qorder = 2 * lfsu_s.finiteElement().localBasis().order();
        // select quadrature rule
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // loop over quadrature points
        for( typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin()
			   ; it!=rule.end()
			   ; ++it )
          {
            // skip rest if we are on Dirichlet boundary
            if( bctype.isDirichlet( ig, it->position() ) )
              continue;

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global( it->position() );
            
            // evaluate velocity
            typename TP::Traits::RangeType beta( tp.velocity( *eg, local, equationMode ) );

            // evaluate test shape functions 
            std::vector<RangeType> phi_u(lfsu_s.size());
            evalFunctionBasis( local, lfsu_s, cache, phi_u );

            std::vector<RangeType> phi_v(lfsv_s.size());
            evalFunctionBasis( local, lfsv_s, cache, phi_v );

            // integrate J
            RF factor = it->weight() * ig.geometry().integrationElement( it->position() );
            RF beta_n = beta * ig.unitOuterNormal( it->position() );

            for( unsigned int i=0; i<lfsv_s.size(); i++ )
			  for( unsigned int j=0; j<lfsu_s.size(); j++ ){
                mat.accumulate( lfsv_s, i, lfsu_s, j, beta_n * phi_u[j] * phi_v[i] * factor );
			  }
		  }
      }





	}; // class StreamlineDiffusionOperator


    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // TRANSPORT_OPERATOR_SDFEM_HH
