// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_L2PROJECTION_OPERATOR_HH
#define DUNE_PDELAB_L2PROJECTION_OPERATOR_HH

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

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

#include "../common/eval.hh"

namespace Dune {
  namespace GeoInversion {

    template<typename DGF>
	class L2ProjectionOperator 
      : 
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags
	{
      
    private:
      const DGF& dgf;
      const REAL l2_diffusion;

	public:
      
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };     // alpha_volume, jacobian_volume
      enum { doLambdaVolume = true };    // lambda_volume


      // The constructor ( with a point source )
      L2ProjectionOperator(const DGF& dgf_,
                           const REAL l2_diffusion_) 
        : dgf(dgf_),
          l2_diffusion(l2_diffusion_)
      {
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
        //const int dimw = EG::Geometry::dimensionworld;

        //const int qorder = 2 * lfsu.finiteElement().localBasis().order();
        const int qorder = 2 * pMAX;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
        

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];
            
            // integrate u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate( lfsu, i, u*phi[i]*factor );


            // add diffusion:
            // get gradient basis: we assume lfsu = lfsv so that we can use gradphi_u = gradphi_v identically
            std::vector<Dune::FieldVector<RF,dim> > gradphi( lfsu.size() );
            evalGradientBasis( eg,
                               it->position(),
                               lfsu,
                               gradphi );

            // compute gradient of u on the global element
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (unsigned int i=0; i<lfsu.size(); i++) {
              gradu.axpy( x(lfsu,i), gradphi[i] );  // w.axpy( a, y ) means: w = w + a*y
            }


            //============ assemble diffusion term  ============
            //
            // + Dgradu * gradv  (standard Galerkin)
            //
            //
            REAL hMin=0;
            REAL hMax=0;
            General::element_size( eg.geometry(), hMin, hMax );
            REAL L = 0.5*(hMax+hMin);
            RF epsilon = l2_diffusion*0.5*L*L;

            for( unsigned int i=0; i<lfsv.size(); i++ ) {
              r.accumulate( lfsu, i, epsilon * ( gradu * gradphi[i] ) * factor );
            }



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
        //const int dimw = EG::Geometry::dimensionworld;

        //const int qorder = 2 * lfsu.finiteElement().localBasis().order();
        const int qorder = 2 * pMAX;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
           
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // integrate phi_j*phi_i (Mass matrix)
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate( lfsu, i, lfsu, j, phi[j]*phi[i]*factor );


            // add diffusion:
            // get gradient basis: we assume lfsu = lfsv so that we can use gradphi_u = gradphi_v identically
            std::vector<Dune::FieldVector<RF,dim> > gradphi( lfsu.size() );
            evalGradientBasis( eg,
                               it->position(),
                               lfsu,
                               gradphi );


            //============ assemble diffusion term  ============
            //
            // + Dgradu * gradv  (standard Galerkin)
            //
            //
            REAL hMin=0;
            REAL hMax=0;
            General::element_size( eg.geometry(), hMin, hMax );
            REAL L = 0.5*(hMax+hMin);
            RF epsilon = l2_diffusion*0.5*L*L;

            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate( lfsu, i, lfsu, j, epsilon * (gradphi[j] * gradphi[i]) * factor );



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

        //typedef typename LFSV::Traits::FiniteElementType::
		//  Traits::LocalBasisType::Traits::JacobianType JacobianType;

		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        //const int dimw = EG::Geometry::dimensionworld;

        // Add to the residual the contributions of the pumping sources or point source in the given locations if needed

        std::vector<RangeType> shapefunctions( lfsv.size() );

        //const int qorder = 2 * lfsv.finiteElement().localBasis().order();
        const int qorder = 2 * pMAX;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions 
            std::vector<RangeType> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            // evaluate right hand side load vector

            Dune::FieldVector<REAL,1> fval(0.0);
            dgf.evaluate(eg.entity(),it->position(),fval);

            // integrate f
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate( lfsv, i, -fval[0]*phi[i]*factor );
          }
      }

	};

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
