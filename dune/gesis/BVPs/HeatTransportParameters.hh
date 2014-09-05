// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_GESIS_HEAT_TRANSPORT_PARAMETERS_HH
#define DUNE_GESIS_HEAT_TRANSPORT_PARAMETERS_HH

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



#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include "dune/gesis/BVPs/regularize.hh"
#include "dune/gesis/BVPs/FEM/deltaSD.hh"

#include "ParameterTraits.hh"

namespace Dune {
  namespace Gesis {

    //! base class for parameter class
    template<class T, class Imp>
    class HeatTransportSpatialParameterInterface {

    private:
      typedef T Traits;
      typedef typename Traits::RangeFieldType RF;
      enum {
        dim = Traits::dimDomain
      };

      typedef typename Traits::DiscreteGridFunction DGF;
      typedef typename Traits::InputDataType IDT;
      typedef typename Traits::SetupDataType SDT;

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

      Imp& asImp() {
        return static_cast<Imp &> (*this);
      }
      const Imp& asImp() const {
        return static_cast<const Imp &> (*this);
      }

    protected: // visible in the derived classes
      const IDT& inputdata;
      const SDT& setupdata;
      const DGF& dgf;

      BCType convertInt2BCType( const int & i ) const {
        if( i == 1 )
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
        if( i == 0 )
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        if( i == -1 )
          return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;

        // default:
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }




    public:
      // constructor of base_type:
      HeatTransportSpatialParameterInterface<T,Imp>( const IDT& inputdata_,
                                                 const SDT& setupdata_,
                                                 const DGF& dgf_ )
      : inputdata(inputdata_)
        , setupdata(setupdata_)
        , dgf(dgf_)
      {};


      const CTransportParameters& gettransportparams() const {
        return inputdata.transport_parameters;
      }



      //! velocityvector
      typename Traits::RangeType
      velocity( const typename Traits::ElementType& e
                , const typename Traits::DomainType& x
                , const EQ::Mode equationMode = EQ::forward
                ) const {

        typename Traits::RangeType flux(0.0);

        dgf.evaluate_on_root( e, x, flux );

        if( equationMode == EQ::adjoint )
          flux*=-1.0;   // Reverse velocity field for adjoint case!

        return flux;

      }




      typename Traits::RangeFieldType deltaSDFEM( const typename Traits::ElementType& eg,
                                                  const typename Traits::DomainType& xlocal
                                                  ) const
      {
        // SDFEM - Parameter: Source: Wathen et.al.: p. 132, equation (3.44) and the footnote!
        const REAL factor = inputdata.transport_parameters.heat_deltaSD_factor;
        REAL Pe_l = 1;
        REAL Pe_t = 1;
        REAL L = 1; // characteristic length of grid element
        typename Traits::RangeType beta;
        getMeshPecletNumbers( eg, xlocal, beta, L, Pe_l, Pe_t );
        return factor * L / beta.two_norm() * zeta_1( Pe_l );
      }
      




      // Get MeshPlecletNumber of heat transport for an element
      // by evaluating the maximal Pe of all quadrature points
      // Pe = |q| / D_l * L;
      // D_l := alpha_l * |q| + lambda_m/(c_w*rho_w)
      //
      void getMeshPecletNumbers( const typename Traits::ElementType& eg,
                                 const typename Traits::DomainType& xlocal,
                                 typename Traits::RangeType & beta,
                                 REAL & L,
                                 REAL & Pe_l,
                                 REAL & Pe_t ) const {
        
        // characteristic length of grid element        
        REAL hMin=0;
        REAL hMax=0;
        General::element_size( eg.geometry(), hMin, hMax );
        L = hMax;
        beta = velocity( eg, xlocal );

        typename Traits::RangeType xglobal
          = eg.geometry().global( xlocal );

        REAL porosity = 1.0;
        REAL lambda_s = 1.0;
        General::zone_parameters( xglobal[dim-1], inputdata, porosity, lambda_s );

        const RF a_l      = inputdata.transport_parameters.heat_a_longitudinal;
        const RF a_t      = inputdata.transport_parameters.heat_a_transversal;
        const RF rho_w    = inputdata.transport_parameters.heat_rho_w;
        const RF c_w      = inputdata.transport_parameters.heat_c_w;
        const RF lambda_w = inputdata.transport_parameters.heat_lambda_w;
        const RF lambda_m = std::pow( lambda_s, 1.0 - porosity ) * std::pow( lambda_w, porosity );

        Pe_l = beta.two_norm() * L / ( a_l * beta.two_norm() + lambda_m/(c_w*rho_w) ); //(O.Cirpka, Numerik-Skript, Seite 72, (5.24))
        Pe_t = beta.two_norm() * L / ( a_t * beta.two_norm() + lambda_m/(c_w*rho_w) );

      }


      // ! Scheidegger dispersion tensor
      //
      // q = specific discharge = porosity * v, v = see page velocity
      // q = -K * grad h
      //
      // The return value here is not D, but (porosity * D).
      // The transport equation reads : div( q*c - porosity * D * grad c )
      //
      typename Traits::DispersionTensorType
      DispersionTensor( const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        const EQ::Mode equationMode = EQ::forward,
                        const RF levelfactor=1.0  // needed here to satisfy interface, but not needed for computation
                        ) const
      {
        typename Traits::DispersionTensorType thetaD( 0.0 );
        typename Traits::RangeType beta( velocity( e, xlocal, equationMode ) );
        typename Traits::DomainType xglobal = e.geometry().global( xlocal );

        CTYPE depth = xglobal[dim-1];
        REAL porosity = 1.0;
        REAL lambda_s = 1.0;
        General::zone_parameters( depth, inputdata, porosity, lambda_s );

        const RF a_l      = inputdata.transport_parameters.heat_a_longitudinal;
        const RF a_t      = inputdata.transport_parameters.heat_a_transversal;
	    const RF rho_w    = inputdata.transport_parameters.heat_rho_w;
	    const RF c_w      = inputdata.transport_parameters.heat_c_w;
	    const RF lambda_w = inputdata.transport_parameters.heat_lambda_w;

	    const RF lambda_m = std::pow( lambda_s, 1.0 - porosity ) * std::pow( lambda_w, porosity );

        if( beta.two_norm() / porosity < GEO_EPSILON*0.5 ){
          // diffusive heat flow
          for( UINT i=0; i<dim; i++ )
            thetaD[i][i] = lambda_m / ( c_w * rho_w );
        }
        else {
          // Scheidegger:
          for( UINT i=0; i<dim; i++ ) {
            thetaD[i][i] += lambda_m/(c_w*rho_w);
            thetaD[i][i] += a_t * beta.two_norm();
            for( UINT j=0; j<dim; j++ ) {
              thetaD[i][j] += beta[i] * beta[j] * ( a_l - a_t ) / beta.two_norm();
            }
          }
        }

        return thetaD;
      }



      //! source/reaction term on the RHS of the scalar convection diffusion equation:
      typename Traits::RangeFieldType
      q_sourceterm( const typename Traits::ElementType& e
                   , const typename Traits::DomainType& x ) const
      {
        RF value = 0;
        return value;
      }



      //! Neumann boundary condition
      //typename Traits::RangeFieldType
      //j(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
      //  return asImp().j(e, x);
      //}
      //! Neumann boundary condition
      typename Traits::RangeFieldType
      j( const typename Traits::IntersectionType& is
         , const typename Traits::IntersectionDomainType& x
         ) const {

          //typename Traits::RangeType global = is.geometry().global(x);

        return 0.0;
      }

      //! reaction term
      typename Traits::RangeFieldType
      c(
        const typename Traits::ElementType& e
        , const typename Traits::DomainType& x) const {
        return 0.0; // no reaction term!
      }

      //! outflow
      typename Traits::RangeFieldType
      o( const typename Traits::IntersectionType& is
         , const typename Traits::IntersectionDomainType& x
         ) const {
          //typename Traits::RangeType global = is.geometry().global(x);
        return 0.0;
      }

      //! source f
      typename Traits::RangeFieldType
      f(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {

        typename Traits::RangeType global = e.geometry().global(x);
        return 0.0;
      }


      typename Traits::RangeFieldType
      setTime( double t ) const {
        return 0.0;
      }


      //! boundary condition type function
      /**
       * 0 means Neumann
       * 1 means Dirichlet
       * 2 means Outflow (zero diffusive flux, velocity points outward)
       */
      //! boundary condition type function
      // 0 means Neumann
      // 1 means Dirichlet
      // -1 means Outflow
      BCType bctype( const typename Traits::IntersectionType& is
                     , const typename Traits::IntersectionDomainType& x
                     ) const {

        typename Traits::RangeType global = is.geometry().global(x);

        // west boundary(0):
        if (global[0] < GEO_EPSILON)
          return convertInt2BCType( setupdata.heat_transport_equation.boundaries[ 0 ].bctype );

        // east boundary(1):
        if (global[0] > inputdata.domain_data.extensions[0] - GEO_EPSILON)
          return convertInt2BCType( setupdata.heat_transport_equation.boundaries[ 1 ].bctype );

        // north boundary(2):
        if (global[1] > inputdata.domain_data.extensions[1] - GEO_EPSILON)
          return convertInt2BCType( setupdata.heat_transport_equation.boundaries[ 2 ].bctype );

        // south boundary(3):
        if (global[1] < GEO_EPSILON)
          return convertInt2BCType( setupdata.heat_transport_equation.boundaries[ 3 ].bctype );

#ifdef DIMENSION3
        // top boundary(2):
        if (global[2] > inputdata.domain_data.extensions[2] - GEO_EPSILON)
          return convertInt2BCType( setupdata.heat_transport_equation.boundaries[ 4 ].bctype );

        // bottom boundary(3):
        if (global[2] < GEO_EPSILON)
          return convertInt2BCType( setupdata.heat_transport_equation.boundaries[ 5 ].bctype );
#endif

        // default: (required to avoid control missing end of non-void function)
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;

      }

      /*********************************************************************
        The following methods will be different for different
        versions of the derived classes! They are actually the reason
        for deriving.
      *********************************************************************/

      //! Dirichlet boundary condition on inflow
      typename Traits::RangeFieldType
      g_function(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
        return asImp().g_function(e, x);
      }

      typename Traits::RangeFieldType
      g(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
        return asImp().g(e, x);
      }

    }; // class HeatTransportSpatialParameterInterface

    // a helper class from Jö
    // see: http://www.boost.org/doc/libs/1_49_0/libs/utility/base_from_member.html
    /*
    template<class T, std::size_t i = 0>
    struct BaseFromMember {
      T data;
      BaseFromMember() { }
      template<class T1>
      BaseFromMember(const T1 &v1) : data(v1) { }
    };
    */





    /*****************************************************************************
          DERIVED CLASS no.1
     *****************************************************************************/

    //! Transport parameters class for m0
    template<typename GV,
             typename RF,
             typename DGFTYPE,
             typename IDT,
             typename SDT
             >
    class HeatTransportProblemM0
      :
      //private BaseFromMember<DGFTYPE>,
      public HeatTransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
      HeatTransportProblemM0<GV,RF,DGFTYPE,IDT,SDT>
      > {

    private:

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      //typedef BaseFromMember<DGFTYPE> pbase_type;

      typedef HeatTransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
        HeatTransportProblemM0<GV,RF,DGFTYPE,IDT,SDT>
        > base_type;

    public:
      typedef TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT> Traits;

      //! constructor
      HeatTransportProblemM0( const DGFTYPE& dgf_,
                          const IDT& inputdata_,
                          const SDT& setupdata_ )
        :
        //pbase_type(dgf_),
        //base_type( pbase_type::data ),
        base_type( inputdata_, setupdata_, dgf_ ){
      }

      //! convectiondiffusionparameter.hh requires g()
      typename Traits::RangeFieldType
      g(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
        return g_function(e,x);
      }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
      g_function( const typename Traits::ElementType& e,
                  const typename Traits::DomainType& x ) const {

        typename Traits::RangeType global = e.geometry().global(x);
        REAL y = 0;

        // TAKE CARE: Tracer stripes only on the WEST boundary!

        if( global[0] > GEO_EPSILON )  // WEST only
          return y;

        RF u = global[0];
        RF v = global[1];

        RF y1 = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].startpoint[0];
        RF y2 = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].endpoint[0];

#ifdef DIMENSION3
        RF w = global[2];
        RF z1 = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].startpoint[1];
        RF z2 = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].endpoint[1];
#endif // DIMENSION3

        RF regularization_factor =
          base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].regularization_factor;

        bool bfixedwidth =
          base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].bfixedwidth;

#ifdef USE_YASP
        REAL delta_y = base_type::inputdata.domain_data.yasp_gridsizes[1];
#ifdef DIMENSION3
        REAL delta_z = base_type::inputdata.domain_data.yasp_gridsizes[2];
#endif
#endif

#ifdef USE_ALUGRID
        REAL delta_y = base_type::inputdata.domain_data.alugrid_gridsizes[1];
#ifdef DIMENSION3
        REAL delta_z = base_type::inputdata.domain_data.alugrid_gridsizes[2];
#endif
#endif

#ifdef USE_UG
        REAL delta_y = base_type::inputdata.domain_data.ug_gridsizes[1];
#ifdef DIMENSION3
        REAL delta_z = base_type::inputdata.domain_data.ug_gridsizes[2];
#endif
#endif

        REAL concentration = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].value;

        REAL injection_time = base_type::setupdata.heat_transport_equation.injection_time;

        if( u < 1e-12 ) // WEST
          {
            RF gDeltaY = regularization_factor * delta_y;

#ifdef DIMENSION3
            RF gDeltaZ = regularization_factor * delta_z;
            return Dune::Gesis::Regularizer::gRegularYZ( v, w,
                                                                y1, z1,
                                                                y2, z2,
                                                                gDeltaY, gDeltaZ,
                                                                concentration * injection_time,
                                                                bfixedwidth,
                                                                regularization_factor );
#else
            return Dune::Gesis::Regularizer::gRegularY( v,
                                                               y1,
                                                               y2,
                                                               gDeltaY,
                                                               concentration * injection_time,
                                                               bfixedwidth,
                                                               regularization_factor );
#endif // DIMENSION3
          }
        else
          return RF(0.0);

      }

    };




    /*****************************************************************************
          DERIVED CLASS no.2
     *****************************************************************************/

    //! Transport parameters class for m0
    template<
      typename GV
      , typename RF
      , typename DGFTYPE
      , typename IDT
      , typename SDT
      >
    class HeatTransportProblemM1
      :
      //private BaseFromMember<DGFTYPE>,
      public HeatTransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
      HeatTransportProblemM1<GV,RF,DGFTYPE,IDT,SDT>
      > {

    private:

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      //typedef BaseFromMember<DGFTYPE> pbase_type;

      typedef HeatTransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
        HeatTransportProblemM1<GV,RF,DGFTYPE,IDT,SDT>
        > base_type;

    public:
      typedef TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT> Traits;

      HeatTransportProblemM1( const DGFTYPE& dgf_,
                          const IDT& inputdata_,
                          const SDT& setupdata_ )
        :
        base_type( inputdata_, setupdata_, dgf_ )
      {}


      //! convectiondiffusionparameter.hh requires g()
      typename Traits::RangeFieldType
      g(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
        return g_function(e,x);
      }

      //! Dirichlet boundary condition value of m1 =
      // Dirichlet boundary condition value of m0 * 0.5 * injection_time * injection_time
      typename Traits::RangeFieldType
      g_function( const typename Traits::ElementType& e,
                  const typename Traits::DomainType& x) const {

        typename Traits::RangeType global = e.geometry().global(x);
        REAL y = 0;

        // TAKE CARE: Tracer stripes only on the WEST boundary!
        if( global[0] > GEO_EPSILON )
          return y;

        RF u = global[0];
        RF v = global[1];

        RF y1 = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].startpoint[0];
        RF y2 = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].endpoint[0];

#ifdef DIMENSION3
        RF w = global[2];
        RF z1 = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].startpoint[1];
        RF z2 = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].endpoint[1];
#endif // DIMENSION3

        RF regularization_factor =
          base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].regularization_factor;

        bool bfixedwidth =
          base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].bfixedwidth;

#ifdef USE_YASP
        REAL delta_y = base_type::inputdata.domain_data.yasp_gridsizes[1];
#ifdef DIMENSION3
        REAL delta_z = base_type::inputdata.domain_data.yasp_gridsizes[2];
#endif
#endif

#ifdef USE_ALUGRID
        REAL delta_y = base_type::inputdata.domain_data.alugrid_gridsizes[1];
#ifdef DIMENSION3
        REAL delta_z = base_type::inputdata.domain_data.alugrid_gridsizes[2];
#endif
#endif

#ifdef USE_UG
        REAL delta_y = base_type::inputdata.domain_data.ug_gridsizes[1];
#ifdef DIMENSION3
        REAL delta_z = base_type::inputdata.domain_data.ug_gridsizes[2];
#endif
#endif

        REAL concentration = base_type::setupdata.heat_transport_equation.boundaries[0].stripes[0].value;

        REAL injection_time = base_type::setupdata.heat_transport_equation.injection_time;

        // IMPORTANT:
        // m1 = integral over (t * c_in) from 0 to injection_time
        //    = 1/2 * t*t * c_in for t = 0 to injection_time
        if( injection_time>GEO_EPSILON )
          concentration *= 0.5 * injection_time * injection_time;


        if( u < 1e-12 ){
          RF gDeltaY = regularization_factor * delta_y;

#ifdef DIMENSION3
          RF gDeltaZ = regularization_factor * delta_z;
          return Dune::Gesis::Regularizer::gRegularYZ( v, w,
                                                              y1, z1,
                                                              y2, z2,
                                                              gDeltaY, gDeltaZ,
                                                              concentration,
                                                              bfixedwidth,
                                                              regularization_factor );
#else
          return Dune::Gesis::Regularizer::gRegularY( v,
                                                             y1,
                                                             y2,
                                                             gDeltaY,
                                                             concentration,
                                                             bfixedwidth,
                                                             regularization_factor );
#endif // DIMENSION3

        }
        else
          return RF(0.0);

      }

    };


    /*****************************************************************************
          DERIVED CLASS no.3
     *****************************************************************************/

    //! Transport parameters class for adjoint states
    template<
      typename GV
      , typename RF
      , typename DGFTYPE
      , typename IDT
      , typename SDT
      >
    class HeatTransportProblemAdjoint
      :
      public HeatTransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
      HeatTransportProblemAdjoint<GV,RF,DGFTYPE,IDT,SDT>
      > {

    private:
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      typedef HeatTransportSpatialParameterInterface<
        TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
        HeatTransportProblemAdjoint<GV,RF,DGFTYPE,IDT,SDT>
        > base_type;

    public:
      typedef TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT> Traits;

      // constructor:
      HeatTransportProblemAdjoint( const DGFTYPE& dgf_,
                               const IDT& inputdata_,
                               const SDT& setupdata_ )
        :
        base_type( inputdata_, setupdata_, dgf_ ){}



      //! boundary condition type function
      /**
       * 0 means Neumann
       * 1 means Dirichlet
       * 2 means Outflow (zero diffusive flux, velocity points outward)
       */
      //! boundary condition type function
      // 0 means Neumann
      // 1 means Dirichlet
      // -1 means Outflow
      BCType bctype( const typename Traits::IntersectionType& is
                     , const typename Traits::IntersectionDomainType& x
                     ) const {

        typename Traits::RangeType global = is.geometry().global(x);

        // west boundary(0):
        if (global[0] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.heat_transport_equation.boundaries[ 1 ].bctype );

        // east boundary(1):
        if (global[0] > base_type::inputdata.domain_data.extensions[0] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.heat_transport_equation.boundaries[ 0 ].bctype );

        // north boundary(2):
        if (global[1] > base_type::inputdata.domain_data.extensions[1] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.heat_transport_equation.boundaries[ 2 ].bctype );

        // south boundary(3):
        if (global[1] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.heat_transport_equation.boundaries[ 3 ].bctype );

#ifdef DIMENSION3
        // top boundary(2):
        if (global[2] > base_type::inputdata.domain_data.extensions[2] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.heat_transport_equation.boundaries[ 4 ].bctype );

        // bottom boundary(3):
        if (global[2] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.heat_transport_equation.boundaries[ 5 ].bctype );
#endif

        // default: (required to avoid control missing end of non-void function)
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
      }




      //! convectiondiffusionparameter.hh requires g()
      typename Traits::RangeFieldType
      g(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
        return g_function(e,x);
      }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
      g_function( const typename Traits::ElementType& e,
                  const typename Traits::DomainType& x) const {

          //typename Traits::RangeType global = e.geometry().global(x);
        REAL y = 0;
        return y;
      }

    };









  } // namespace PDELab

} // namespace Dune

#endif // DUNE_GESIS_HEAT_TRANSPORT_PARAMETERS_HH
