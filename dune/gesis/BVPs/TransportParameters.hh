// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef TRANSPORT_PARAMETERS_HH
#define TRANSPORT_PARAMETERS_HH

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

#include "../common/ParameterTraits.hh"

namespace Dune {
  namespace GeoInversion {

    //! base class for parameter class
    template<class T, class Imp>
    class TransportSpatialParameterInterface {

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
      TransportSpatialParameterInterface<T,Imp>( const IDT& inputdata_,
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
                , const EQ::Mode equationMode
                ) const {

        typename Traits::RangeType flux(0.0);

        dgf.evaluate_on_root( e, x, flux );

        if( equationMode == EQ::adjoint )
          flux *= -1.0;   // Reverse velocity field for adjoint case!

        return flux;

      }



      typename Traits::RangeFieldType deltaSDFEM( const typename Traits::ElementType& eg,
                                                  const typename Traits::DomainType& xlocal
                                                  ) const
      {
        // SDFEM - Parameter: Source: Wathen et.al.: p. 132, equation (3.44) and the footnote!
        const REAL factor = inputdata.transport_parameters.deltaSD_factor;
        REAL Pe_l = 1;
        REAL Pe_t = 1;
        REAL L = 1; // characteristic length of grid element
        typename Traits::RangeType beta;
        getMeshPecletNumbers( eg, xlocal, beta, L, Pe_l, Pe_t );
        return factor * 0.5 * L / beta.two_norm() * zeta_1( Pe_l );
      }


      // Get MeshPlecletNumber of solute transport for an element
      // by evaluating the maximal Pe of all quadrature points
      // Pe = |q| / D_l * L;
      // D_l := alpha_l * |q| + porosity * D_m
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

        const EQ::Mode equationMode = EQ::forward;// mode does not matter here
        beta = velocity( eg, xlocal, equationMode );

        typename Traits::RangeType xglobal
          = eg.geometry().global( xlocal );

        REAL porosity;
        General::zone_parameters( xglobal[dim-1], inputdata, porosity );

        const RF D_molecular = inputdata.transport_parameters.D_molecular;
        RF extra_diffusion = 1.0;
#if defined USE_UG || defined USE_ALUGRID
        const REAL amplifier = inputdata.transport_parameters.refined_diffusion;
        REAL levelfactor_1 = REAL(eg.level());
        if( levelfactor_1 > 1E-12 && amplifier > 1E-12 )
          extra_diffusion = levelfactor_1 * amplifier;
#endif
        RF D_m = D_molecular * extra_diffusion;

        const REAL a_l = inputdata.transport_parameters.a_longitudinal;
        const REAL a_t = inputdata.transport_parameters.a_transversal;

        Pe_l = 0.5 * beta.two_norm() * L / ( a_l * beta.two_norm() + porosity * D_m); // (Cirpka/Kitanidis 2001)
        Pe_t = 0.5 * beta.two_norm() * L / ( a_t * beta.two_norm() + porosity * D_m );

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
                        const EQ::Mode equationMode
                        ) const 
      {

        typename Traits::DispersionTensorType thetaD( 0.0 );
        typename Traits::RangeType beta( velocity( e, xlocal, equationMode ) );
        typename Traits::DomainType xglobal = e.geometry().global( xlocal );

        RF porosity;
        General::zone_parameters( xglobal[dim-1], inputdata, porosity );

        const RF a_l = inputdata.transport_parameters.a_longitudinal;
        const RF a_t = inputdata.transport_parameters.a_transversal;
        const RF D_molecular = inputdata.transport_parameters.D_molecular;

        RF extra_diffusion = 1.0;

#if defined USE_UG || defined USE_ALUGRID
        const REAL amplifier = inputdata.transport_parameters.refined_diffusion;
        REAL levelfactor_1 = REAL(e.level());
        if( levelfactor_1 > 1E-12 && amplifier > 1E-12 )
          extra_diffusion = levelfactor_1 * amplifier;
#endif

        RF D_m = D_molecular * extra_diffusion;
        
        for( UINT i=0; i<dim; i++ ) {
          thetaD[i][i] += porosity * D_m;
          thetaD[i][i] += a_t * beta.two_norm();
          for( UINT j=0; j<dim; j++ ) {
            thetaD[i][j] += beta[i] * beta[j] * ( a_l - a_t ) / beta.two_norm();
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

#ifdef USE_REFERENCE_SOLUTION_001
        // Reference Solution from "Presentation_SIAM_SEAS_08.pdf" (Hoang)
        // f(x,y) = (x+y)*(1.0 - exp((x-1.0)/epsilon) * exp((y-1.0)/epsilon))
        //        + (x-y)*( exp((y-1.0)/epsilon) - exp((x-1.0)/epsilon) )

        RF D_m = inputdata.transport_parameters.D_molecular;
        RF xfactor = exp((x[0]-1.0)/D_m);
        RF yfactor = exp((x[1]-1.0)/D_m);
        value = ( x[0] + x[1] ) * ( 1.0 - xfactor * yfactor )
          + ( x[0] - x[1] ) * ( yfactor - xfactor );
#endif // USE_REFERENCE_SOLUTION_001

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
          
        //typename Traits::RangeType global = e.geometry().global(x);
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
          return convertInt2BCType( setupdata.transport_equation.boundaries[ 0 ].bctype );
          
        // east boundary(1):
        if (global[0] > inputdata.domain_data.extensions[0] - GEO_EPSILON)
          return convertInt2BCType( setupdata.transport_equation.boundaries[ 1 ].bctype );

        // south boundary(2):
        if (global[1] < GEO_EPSILON)
          return convertInt2BCType( setupdata.transport_equation.boundaries[ 2 ].bctype );
            
        // north boundary(3):
        if (global[1] > inputdata.domain_data.extensions[1] - GEO_EPSILON)
          return convertInt2BCType( setupdata.transport_equation.boundaries[ 3 ].bctype );

#ifdef DIMENSION3
        // bottom boundary(4):
        if (global[2] < GEO_EPSILON)
          return convertInt2BCType( setupdata.transport_equation.boundaries[ 4 ].bctype );

        // top boundary(5):
        if (global[2] > inputdata.domain_data.extensions[2] - GEO_EPSILON)
          return convertInt2BCType( setupdata.transport_equation.boundaries[ 5 ].bctype );            
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

    }; // class TransportSpatialParameterInterface

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
    class TransportProblemM0
      :
      //private BaseFromMember<DGFTYPE>,
      public TransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
      TransportProblemM0<GV,RF,DGFTYPE,IDT,SDT>
      > {

    private:

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      //typedef BaseFromMember<DGFTYPE> pbase_type;

      typedef TransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
        TransportProblemM0<GV,RF,DGFTYPE,IDT,SDT>
        > base_type;

    public:
      typedef TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT> Traits;

      //! constructor
      TransportProblemM0( const DGFTYPE& dgf_,
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

        RF y1 = base_type::setupdata.transport_equation.boundaries[0].stripes[0].startpoint[0];
        RF y2 = base_type::setupdata.transport_equation.boundaries[0].stripes[0].endpoint[0];

#ifdef DIMENSION3
        RF w = global[2];
        RF z1 = base_type::setupdata.transport_equation.boundaries[0].stripes[0].startpoint[1];
        RF z2 = base_type::setupdata.transport_equation.boundaries[0].stripes[0].endpoint[1];
#endif // DIMENSION3

        RF regularization_factor =
          base_type::setupdata.transport_equation.boundaries[0].stripes[0].regularization_factor;

        bool bfixedwidth =
          base_type::setupdata.transport_equation.boundaries[0].stripes[0].bfixedwidth;

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

        REAL concentration = base_type::setupdata.transport_equation.boundaries[0].stripes[0].value;

        REAL injection_time = base_type::setupdata.transport_equation.injection_time;

#ifdef LINKPOINTS_REGULARIZATION
              
        if( u < 1e-12 ) // WEST
        //if( u > 500.0 - 1e-12 ) // EAST
          return Dune::GeoInversion::Regularizer::gRegular2( v, y1, y2, delta_y, concentration * injection_time);
        else
          return RF(0.0);

#else
        if( u < 1e-12 ) // WEST
        //if( u > 500.0 - 1e-12 ) // EAST
          {
            RF gDeltaY = regularization_factor * delta_y;

#ifdef DIMENSION3
            RF gDeltaZ = regularization_factor * delta_z;
            return Dune::GeoInversion::Regularizer::gRegularYZ( v, w, 
                                                                y1, z1,
                                                                y2, z2,
                                                                gDeltaY, gDeltaZ,
                                                                concentration * injection_time,
                                                                bfixedwidth,
                                                                regularization_factor );
#else
            return Dune::GeoInversion::Regularizer::gRegularY( v, 
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
#endif

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
    class TransportProblemM1
      :
      //private BaseFromMember<DGFTYPE>,
      public TransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
      TransportProblemM1<GV,RF,DGFTYPE,IDT,SDT>
      > {

    private:

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      //typedef BaseFromMember<DGFTYPE> pbase_type;

      typedef TransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
        TransportProblemM1<GV,RF,DGFTYPE,IDT,SDT>
        > base_type;

    public:
      typedef TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT> Traits;

      TransportProblemM1( const DGFTYPE& dgf_,
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

        RF y1 = base_type::setupdata.transport_equation.boundaries[0].stripes[0].startpoint[0];
        RF y2 = base_type::setupdata.transport_equation.boundaries[0].stripes[0].endpoint[0];

#ifdef DIMENSION3
        RF w = global[2];
        RF z1 = base_type::setupdata.transport_equation.boundaries[0].stripes[0].startpoint[1];
        RF z2 = base_type::setupdata.transport_equation.boundaries[0].stripes[0].endpoint[1];
#endif // DIMENSION3

        RF regularization_factor =
          base_type::setupdata.transport_equation.boundaries[0].stripes[0].regularization_factor;

        bool bfixedwidth =
          base_type::setupdata.transport_equation.boundaries[0].stripes[0].bfixedwidth;

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

        REAL concentration = base_type::setupdata.transport_equation.boundaries[0].stripes[0].value;

        REAL injection_time = base_type::setupdata.transport_equation.injection_time;

        // IMPORTANT:
        // m1 = integral over (t * c_in) from 0 to injection_time
        //    = 1/2 * t*t * c_in for t = 0 to injection_time
        if( injection_time>GEO_EPSILON )
          concentration *= 0.5 * injection_time * injection_time;


#ifdef LINKPOINTS_REGULARIZATION
              
        if( u < 1e-12 )
          return Dune::GeoInversion::Regularizer::gRegular2( v, y1, y2, delta_y, concentration );
        else
          return RF(0.0);

#else
        if( u < 1e-12 ){
          RF gDeltaY = regularization_factor * delta_y;

#ifdef DIMENSION3
          RF gDeltaZ = regularization_factor * delta_z;
          return Dune::GeoInversion::Regularizer::gRegularYZ( v, w, 
                                                              y1, z1,
                                                              y2, z2,
                                                              gDeltaY, gDeltaZ,
                                                              concentration,
                                                              bfixedwidth,
                                                              regularization_factor );
#else
          return Dune::GeoInversion::Regularizer::gRegularY( v, 
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
#endif
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
    class TransportProblemAdjoint
      :
      public TransportSpatialParameterInterface<
      TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
      TransportProblemAdjoint<GV,RF,DGFTYPE,IDT,SDT>
      > {
    
    private:
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      typedef TransportSpatialParameterInterface<
        TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT>,
        TransportProblemAdjoint<GV,RF,DGFTYPE,IDT,SDT>
        > base_type;

    public:
      typedef TransportParameterTraits<GV,RF,DGFTYPE,IDT,SDT> Traits;

      // constructor:
      TransportProblemAdjoint( const DGFTYPE& dgf_,
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

        // west boundary(0) = east boundary of forward problem:
        if (global[0] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.transport_equation.boundaries[ 1 ].bctype );
          
        // east boundary(1) = west boundary of forward problem:
        if (global[0] > base_type::inputdata.domain_data.extensions[0] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.transport_equation.boundaries[ 0 ].bctype );
          
        // south boundary(2) = north boundary of forward problem:
        if (global[1] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.transport_equation.boundaries[ 3 ].bctype );

        // north boundary(3) = south boundary of forward problem:
        if (global[1] > base_type::inputdata.domain_data.extensions[1] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.transport_equation.boundaries[ 2 ].bctype );
            
#ifdef DIMENSION3
        // bottom boundary(4) = top boundary of forward problem:
        if (global[2] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.transport_equation.boundaries[ 5 ].bctype );
            
        // top boundary(5) = bottom boundary of forward problem:
        if (global[2] > base_type::inputdata.domain_data.extensions[2] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.transport_equation.boundaries[ 4 ].bctype );
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
        /*
        REAL y = 0;
        REAL Lx = base_type::inputdata.domain_data.extensions[0];
        //if( global[0] > GEO_EPSILON ) // WEST only!
        if( global[0] < Lx - GEO_EPSILON ) // EAST only!
          return y;
        if( global[1] > 250 )
          y = 1.0;
        */
        REAL y = 0;
        return y;
      }
      
    };








    
  } // namespace PDELab

} // namespace Dune

#endif //TRANSPORT_PARAMETERS_HH
