// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_GESIS_GW_LOP_GALERKIN_HH
#define DUNE_GESIS_GW_LOP_GALERKIN_HH



    /** a local operator for solving the stationary groundwater equation
     *
     * \f{align*}{
     * \nabla \cdot\{ - K(x) \nabla u \}  & = & f               \mbox{ in } \Omega,          \\
     *                                  u & = & g               \mbox{ on } \partial\Omega_D \\
     *              ( - K(x) \nabla u ) \cdot \nu & = & j       \mbox{ on } \partial\Omega_N \\
     * \f}
     * with conforming finite elements on all types of grids in any dimension
     */

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

//#include<dune/common/geometrytype.hh>
//#include<dune/grid/common/genericreferenceelements.hh>
//#include<dune/grid/common/quadraturerules.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
//#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

namespace Dune {
  namespace Gesis {

    template< typename GWP
              , typename BCType
              , typename SourceType
              , typename IDT
              , typename SDT
              >
	class GwFlowOperator 
      : 
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags
	{
      
    private:
      enum{dim=GWP::dim};
      const GWP& gwp;
      const BCType& bctype;
      SourceType& sourceterm;
      const IDT& inputdata;
      const SDT& setupdata;

      const EQ::Mode equationMode;

	public:
      
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };     // alpha_volume, jacobian_volume
      enum { doLambdaVolume = true };    // lambda_volume
      enum { doLambdaBoundary = true };  // lambda_boundary
      enum { doLambdaSkeleton = false };  // lambda_skeleton

      
      template<typename VCType>
      void set_rhs( const VCType& vc1 )
      {
        sourceterm.set_rhs( vc1 );
      }

      template<typename VCType>
      void set_rhs( const VCType& vc1, const VCType& vc2 )
      {
        sourceterm.set_rhs( vc1, vc2 );
      }

      template<typename UType,typename VCType>
      void set_rhs( const UType& vc1,   // eqn (66): m_k
                    const VCType& vc2   // eqn (66): phi0
                    )
      {
        sourceterm.set_rhs( vc1, vc2 );
      }

      template<typename VCType>
      void set_rhs( const VCType& vc1, const VCType& vc2, const VCType& vc3, const VCType& vc4 )
      {
        sourceterm.set_rhs( vc1, vc2, vc3, vc4 );
      }

      template<typename COORDINATES>
      void set_PointSourceLocation( const COORDINATES& x )
      {
        sourceterm.set_PointSourceLocation( x );
      }
      
      template<typename COORDINATES>
      void set_PointSourceLocation( const COORDINATES& x1, const COORDINATES& x2)
      {
        sourceterm.set_PointSourceLocations( x1,x2 );
      }

      // The constructor ( with a point source )
      GwFlowOperator( 
                     const GWP& gwp_, 
                     const BCType& bctype_, 
                     SourceType& sourceterm_,
                     const IDT& inputdata_,
                     const SDT& setupdata_,
                     const EQ::Mode& equationMode_
                      ) : 
        gwp(gwp_), 
        bctype(bctype_), 
        sourceterm( sourceterm_ ),
        inputdata( inputdata_ ),
        setupdata( setupdata_ ),
        equationMode( equationMode_ )
      {}




      // Ausgelagerte Funktionen zur Berechnung des Brunnens im Fracture-Model:

      // Checks whether 'elementpoint' lies within epsilon reach of 'well_postion' and above 'well_bottom'
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
        // should not be used in the GPE!!!
        assert( sourceterm.source_nature != GEOELECTRICAL_POTENTIAL_SOURCE && 
                sourceterm.source_nature != GEP_FUNCTIONAL_SOURCE );

        // extract well-data:
        RF well_position_x = setupdata.wdlist.pointdata_vector[iWell].x;
        RF well_top =        setupdata.wdlist.pointdata_vector[iWell].well_top;
        RF well_bottom =     setupdata.wdlist.pointdata_vector[iWell].well_bottom;
#ifdef DIMENSION3
        RF well_position_y = setupdata.wdlist.pointdata_vector[iWell].y;
#endif        
        if( elementpoint[0] >= well_position_x - reach
            && elementpoint[0] < well_position_x + reach
            
#ifdef DIMENSION3
            && elementpoint[1] >= well_position_y - reach_y
            && elementpoint[1] < well_position_y + reach_y
#endif
            && elementpoint[dim-1] > well_bottom - GEO_EPSILON
            && elementpoint[dim-1] < well_top + GEO_EPSILON
            ){
          return true;
        }
        else
          return false;
      }






      template<typename EG, typename LFSV, typename R>
      void addWellTo_Lambda_Volume( const UINT iWell, 
                                    const EG& eg, 
                                    const LFSV& lfsv,
                                    R& residual
                                    ) const {
        //return; // testwise
        
        // should not be used in the GPE!!!
        if( sourceterm.source_nature == GEOELECTRICAL_POTENTIAL_SOURCE ||
            sourceterm.source_nature == GEP_FUNCTIONAL_SOURCE ){
          return;
        }

		// domain and range field type
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		//typedef typename LFSV::Traits::FiniteElementType::
		//  Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;

        // extract well-data:
        RF well_rate =   setupdata.wdlist.pointdata_vector[iWell].well_rate;
        
        if(std::abs(well_rate)<1e-12)
          return; // no contribution


        // dimensions
        //const int dimw = EG::Geometry::dimensionworld;

        Dune::FieldVector<REAL,dim-1> well_position;
        well_position[0] =  setupdata.wdlist.pointdata_vector[iWell].x;
#ifdef DIMENSION3
        well_position[1] =  setupdata.wdlist.pointdata_vector[iWell].y;
#endif

        RF well_top =   setupdata.wdlist.pointdata_vector[iWell].well_top;
        RF well_bottom = setupdata.wdlist.pointdata_vector[iWell].well_bottom;

        if( w_outside != 
            isElementPenetratedByWell(eg,iWell,inputdata,setupdata) ){

          //REAL refinement_factor = std::pow(2.0,eg.entity().level()*dim);

          const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1];
          const REAL screening_fraction = dz / ( well_top - well_bottom );

          const int qorder = 2 * lfsv.finiteElement().localBasis().order();
          // select quadrature rule
          Dune::GeometryType gt = eg.geometry().type();
          const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
        
          // evaluate Y-field tensor at cell center, assume it is constant over elements
          Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
          // typename GWP::Traits::DiffusionTensorType tensor(gwp.DiffusionTensor(eg.entity(),localcenter));

          // loop over quadrature points
          for( typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin()
                 ; it!=rule.end(); ++it ) {

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            //std::vector<JacobianType> js(lfsv.size());
            //lfsv.finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            //const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            //std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsv.size());
            //for (size_type i=0; i<lfsv.size(); i++) {
            //  gradphi[i] = 0.0;
            //  jac.umv(js[i][0],gradphi[i]);
            //}

            // compute gradient of u
            //Dune::FieldVector<RF,dim> gradu(0.0);
            //for (size_type i=0; i<lfsv.size(); i++)
            //  gradu.axpy(x(lfsv,i),gradphi[i]);

            // compute K * gradient of u
            //Dune::FieldVector<RF,dim> Kgradu(0.0);
            //tensor.umv(gradu,Kgradu);

            // evaluate basis functions
            std::vector<RangeType> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            //RF u=0.0;
            //for (size_type i=0; i<lfsv.size(); i++)
            //  u += x(lfsv,i)*phi[i];

            // evaluate reaction term
            //typename GWP::Traits::RangeFieldType c = gwp.c(eg.entity(),it->position());

            // integrate (K grad u)*grad phi_i + a_0*u*phi_i

            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            RF fval = well_rate * screening_fraction / eg.geometry().volume();
            //fval *= refinement_factor;
              
            for (size_type i=0; i<lfsv.size(); i++)
              residual.accumulate( lfsv, i, -fval*phi[i]*factor );
            // r.accumulate( lfsv, i, ( Kgradu*gradphi[i] + c*u*phi[i] )*factor );
          }  // end of loop over quadrature points

        }

      }









      template<
        typename EG,
        typename LFSU, 
        typename X,
        typename M
        >
      void addWellToJacobianVolumeMatrix(
                                         const int& iWell, 
                                         const EG& eg, 
                                         const LFSU& lfsu, 
                                         const X& x, 
                                         M& mat
                                         ) const
      {

#ifdef USE_LINE_INTEGRAL
#else
        return; // testwise
#endif

	// should not be used in the GPE!!!
        assert( sourceterm.source_nature!=GEOELECTRICAL_POTENTIAL_SOURCE &&
                sourceterm.source_nature!=GEP_FUNCTIONAL_SOURCE );

		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
		//typedef typename LFSU::Traits::FiniteElementType::
		//  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // Well...

        // Loop over the corners of the element and find the two corners lying in the vertical well:
        // typename EG::Geometry elementgeometry = eg.geometry();

        
        const RF cellReach = inputdata.domain_data.yasp_gridsizes[0];
        const RF cellReach_y=inputdata.domain_data.yasp_gridsizes[1];
        //std::cout << "cellcenter " << eg.geometry().center() << " has " << cellReach << std::endl;

#ifdef DIMENSION3
	if( !isPointInsideReachOfWell( eg.geometry().center(), cellReach,cellReach_y, iWell) )
#else
        if( !isPointInsideReachOfWell( eg.geometry().center(), cellReach, iWell) )
#endif
          {
            //std::cout << "cellcenter " << eg.geometry().center() << " is not within " << cellReach << " - reach of well no." << iWell << std::endl;
            return;  // Adrian: This is to improve performance of this algorithm! If Well is not within reach of this element center, there is absolutely no need to run over the corners!
          }
        
        std::vector<Dune::FieldVector<RF,dim>> well_points;
        std::vector<size_t> well_points_cornerindex;

        size_t nElementCorners = eg.geometry().corners();
        for( size_t iElementCorner = 0; iElementCorner < nElementCorners; iElementCorner++ )
          {
            const RF epsilonReach = cellReach / 2.0; // GEO_EPSILON;
#ifdef DIMENSION3
	    const RF epsilonReach_y = cellReach_y / 2.0; // GEO_EPSILON;
#endif
            Dune::FieldVector<RF,dim> elementcorner = eg.geometry().corner( iElementCorner );
            
#ifdef DIMENSION3
	    if( isPointInsideReachOfWell( elementcorner, epsilonReach,epsilonReach_y, iWell ) )
#else
            if( isPointInsideReachOfWell( elementcorner, epsilonReach, iWell ) )
#endif
              {
                well_points.push_back( elementcorner );
                well_points_cornerindex.push_back( iElementCorner );
                
                //std::cout << "  iElementCorner=" << iElementCorner
                //          << "  elementcorner=" << elementcorner
                //          << std::endl;
              }
          }

        if( well_points.size() > 2 ){
          std::cout << "WARNING: Something went wrong. There are more than two element corners lying in one well. This cannot be! well_points.size()= "<<well_points.size() << std::endl;
	}

        if( well_points.size() > 1 )
          {
            Dune::FieldVector<RF,dim> edgevector = well_points[1];
            edgevector -= well_points[0];
            //logger << "  edgevector =" << edgevector
            //          << std::endl;

            const int qorder = 2 * lfsu.finiteElement().localBasis().order();

            const Dune::QuadratureRule<DF,1>& rule_on_edge = Dune::QuadratureRules<DF,1>::rule( Dune::GeometryType( Dune::GeometryType::cube, 1), qorder );

            // loop over quadrature points on the edge
            for( typename Dune::QuadratureRule<DF,1>::const_iterator it_qPointOnEdge = rule_on_edge.begin()
                   ; it_qPointOnEdge != rule_on_edge.end()
                   ; ++it_qPointOnEdge )
              {
                Dune::FieldVector<RF,1> qPointOnEdge_edgelocal =  it_qPointOnEdge->position();
                //std::cout << "qPointOnEdge_edgelocal = " << qPointOnEdge_edgelocal << std::endl;
                
                Dune::FieldVector<RF,dim> qPointOnEdge_global =  well_points[0];
                qPointOnEdge_global.axpy( qPointOnEdge_edgelocal, edgevector );
                //std::cout << "qPointOnEdge_global = " << qPointOnEdge_global << std::endl;
                
                Dune::FieldVector<RF,dim> qPointOnEdge_elementlocal =  eg.geometry().local( qPointOnEdge_global );
                //std::cout << "qPointOnEdge_elementlocal = " << qPointOnEdge_elementlocal << std::endl;

                // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
                std::vector<JacobianType> js_edge(lfsu.size());
                lfsu.finiteElement().localBasis().evaluateJacobian( qPointOnEdge_elementlocal, js_edge );

                // transform gradient to real element
                const Dune::FieldMatrix<DF,dimw,dim> jac_edge = eg.geometry().jacobianInverseTransposed( qPointOnEdge_elementlocal );
                std::vector<Dune::FieldVector<RF,dim> > gradphi_qPointOnEdge( lfsu.size() );
                for (size_type i=0; i<lfsu.size(); i++)
                  {
                    gradphi_qPointOnEdge[i] = 0.0;
                    jac_edge.umv( js_edge[i][0], gradphi_qPointOnEdge[i] );
                  }
                
                // compute gradient of u
                //Dune::FieldVector<RF,dim> gradu_qPointOnEdge(0.0);
                //for( size_type i=0; i<lfsu.size(); i++ )
                //  gradu_qPointOnEdge.axpy( x(lfsu,i), gradphi_qPointOnEdge[i] );
                //RF gradu_edge = edgevector * gradu_qPointOnEdge; // = gradient of u projected onto the edge

                RF edge_factor = it_qPointOnEdge->weight() * edgevector.two_norm();
                //std::cout << "edgevector.two_norm() = " << edgevector.two_norm() << std::endl;

#ifdef DIMENSION3
                RF dimensional_factor = 0.25; // each edge is generally surrounded by four elements in 3D
#else
                RF dimensional_factor = 0.5;  // each edge is generally surrounded by two elements in 2D
#endif
                RF well_conductivity=setupdata.wdlist.pointdata_vector[iWell].well_conductivity;     
                
                for (size_type j=0; j<lfsu.size(); j++)
                  for (size_type i=0; i<lfsu.size(); i++)
                    {
                      RF extra_well_contribution = dimensional_factor * well_conductivity;
                      extra_well_contribution *= ( edgevector * gradphi_qPointOnEdge[j] );
                      extra_well_contribution *= ( edgevector * gradphi_qPointOnEdge[i] );
                      extra_well_contribution *= edge_factor;
                      mat.accumulate( lfsu, i, lfsu, j, extra_well_contribution );
                    }
              }
          }
      }







	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume( 
                        const EG& eg
                        , const LFSU& lfsu
                        , const X& x
                        , const LFSV& lfsv
                        , R& r
                         ) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        const int qorder = 2 * lfsu.finiteElement().localBasis().order();
        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
        
        // evaluate Y-field tensor at cell center, assume it is constant over elements
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        typename GWP::Traits::DiffusionTensorType tensor(gwp.DiffusionTensor(eg.entity(),localcenter));

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphi[i]);

            // compute K * gradient of u
            Dune::FieldVector<RF,dim> Kgradu(0.0);
            tensor.umv(gradu,Kgradu);

            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate reaction term
            typename GWP::Traits::RangeFieldType c = gwp.c(eg.entity(),it->position());

            // integrate (K grad u)*grad phi_i + a_0*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate( lfsu, i, ( Kgradu*gradphi[i] + c*u*phi[i] )*factor );
          }  // end of loop over quadrature points
              
	  }




      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
	  void jacobian_volume(
                           const EG& eg
                           , const LFSU& lfsu
                           , const X& x
                           , const LFSV& lfsv
                           , M& mat
                           ) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        const int qorder = 2 * lfsu.finiteElement().localBasis().order();
        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // evaluate Y-field tensor at cell center, assume it is constant over elements
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        typename GWP::Traits::DiffusionTensorType tensor(gwp.DiffusionTensor(eg.entity(),localcenter));

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradient to real element
            const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              {
                gradphi[i] = 0.0;
                jac.umv(js[i][0],gradphi[i]);
              }

            // compute K * gradient of shape functions
            std::vector<Dune::FieldVector<RF,dim> > Kgradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
              tensor.mv(gradphi[i],Kgradphi[i]);
            
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate reaction term
            typename GWP::Traits::RangeFieldType c = gwp.c(eg.entity(),it->position());

            // integrate (K grad phi_j)*grad phi_i + a_0*phi_j*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate( lfsu, i, lfsu, j, ( Kgradphi[j]*gradphi[i] + c*phi[j]*phi[i] )*factor );
          } // end of loop over quadrature points

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
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;

		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        

        // Add to the residual the contributions of the pumping sources or point source in the given locations if needed
        std::vector<RangeType> shapefunctions( lfsv.size() );
        if( sourceterm.source_nature == PUMPING_SOURCES )
          {
            sourceterm.evaluate_residual_on_element( eg.entity(),
                                                     lfsv,
                                                     shapefunctions,
                                                     r );
          }
        
        if( sourceterm.source_nature == POINT_SOURCE || sourceterm.source_nature == GEOELECTRICAL_POINT_SOURCE)
          {
            sourceterm.evaluate_residual_on_element( eg.entity(),
                                                     lfsv,
                                                     shapefunctions,
                                                     r );
          }
        if( equationMode==EQ::forward &&
            sourceterm.source_nature == GEOELECTRICAL_POTENTIAL_SOURCE )
          {
            sourceterm.evaluate_residual_on_element( eg.entity(),
                                                     lfsv,
                                                     shapefunctions,
                                                     r );
          }
          

        

        const int qorder = 2 * lfsv.finiteElement().localBasis().order();
        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
/*	
			if (sourceterm.source_nature == FUNCTIONAL_SOURCE)
logger<<"GWE lambda_volume: "<<std::endl;
else if (sourceterm.source_nature == GP_FUNCTIONAL_SOURCE)
  logger<<"GPE lambda_volume: "<<std::endl;
*/
        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions 
            std::vector<RangeType> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            // evaluate right hand side parameter function
            // evaluate source term
            typename GWP::Traits::RangeFieldType fval = gwp.f(eg.entity(),it->position());

            // ============= Source term for adjoint flow problem wrt m0 or wrt. m0m1:  ==============
            //
            // - ( psi_m0 Kgrad m0 ) * gradv 
            // or:
            // - ( psi_m0 Kgrad m0 + psi_m1 Kgrad m1  ) * gradv 
            //
            if( ( equationMode==EQ::adjoint &&
                  sourceterm.source_nature == FUNCTIONAL_SOURCE) ||
                sourceterm.source_nature == GEP_FUNCTIONAL_SOURCE )
              {

                // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
                std::vector<JacobianType> js(lfsv.size());
                lfsv.finiteElement().localBasis().evaluateJacobian(it->position(),js);

                // transform gradient to real element
                const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
                std::vector<Dune::FieldVector<RF,dim> > gradphi_v(lfsv.size());
                for (size_type i=0; i<lfsv.size(); i++)
                  {
                    gradphi_v[i] = 0.0;
                    jac.umv(js[i][0],gradphi_v[i]);
                  }
                
                Dune::FieldVector<RF,dim> fvec(0.0);
                sourceterm.evaluate_vector( eg.entity(),
                                            it->position(),
                                            fvec );
		//logger<<"sourceterm value: "<<fvec[0]<<", "<<fvec[1]<<", at:"<< it->position()<<std::endl;
                for( unsigned int i=0; i<lfsv.size(); i++ )
                  r.accumulate( lfsv, i, -(fvec * gradphi_v[i]) * factor );
              }


            // integrate f
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate( lfsv, i, -fval*phi[i]*factor );
          }


        // Add inflow or outlow rates for all the wells
        UINT nWells = setupdata.wdlist.total;
        for( UINT iWell=0; iWell<nWells; iWell++ ) {
          addWellTo_Lambda_Volume( iWell, 
                                   eg, 
                                   lfsv,
                                   r
                                   );
        }


      }



      template<typename IG, typename LFSV, typename R>
      void lambda_skeleton (const IG& ig, 
                            const LFSV& lfsv_s, const LFSV& lfsv_n, 
                            R& r_s, R& r_n) const
      {
        if( equationMode==EQ::adjoint &&
            sourceterm.source_nature == FUNCTIONAL_SOURCE )
          sourceterm.evaluate_residual_on_skeleton( ig, 
                                                    lfsv_s, lfsv_n, 
                                                    r_s, r_n );
      }



      // boundary integral independen of ansatz functions
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

        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        
        const int qorder = 2 * lfsv.finiteElement().localBasis().order();
        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {

            Dune::FieldVector<DF,dim-1> qPoint_local = it->position(); // local coordinate!

            // evaluate boundary condition type

            //typename B::Traits::RangeType bctype;
            //b.evaluate( ig, qPoint_local, bctype );
            //if (bctype>0) continue;

            // skip rest if we are on Dirichlet boundary
            if( bctype.isDirichlet( ig, qPoint_local ) )
              continue;

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global( qPoint_local );
            RF factor = it->weight() * ig.geometry().integrationElement( qPoint_local );

            // evaluate test shape functions 
            std::vector<RangeType> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(local,phi);
            
            // evaluate flux boundary condition
            typename GWP::Traits::RangeFieldType j = gwp.j( ig.intersection(), it->position() );
            
            // integrate J
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate( lfsv, i, j*phi[i]*factor );


            // ============= Neumann Boundary term for adjoint flow problem wrt m0 or wrt m0m1:  ==============
            //
            // - ( psi_m0 Kgrad m0 ) * n
            // or:
            // - ( psi_m0 Kgrad m0 + psi_m1 Kgrad m1  ) * n
            //
            if( sourceterm.source_nature == FUNCTIONAL_SOURCE ||
                sourceterm.source_nature == GEP_FUNCTIONAL_SOURCE )
              {
                Dune::FieldVector<RF,dim> fvec(0.0);
                sourceterm.evaluate_vector( *(ig.inside()),
                                            local,
                                            fvec );
            
                for( unsigned int i=0; i<lfsv.size(); i++ )
                  r.accumulate( lfsv, i, (fvec * ig.unitOuterNormal( qPoint_local )) * phi[i] * factor );
              }

          }


      }

	};


  } // namespace PDELab
} // namespace Dune

#endif
