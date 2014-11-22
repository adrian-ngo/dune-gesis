// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_GESIS_GEO_ELECTRICAL_POTENTIAL_PARAMETERS_HH
#define DUNE_GESIS_GEO_ELECTRICAL_POTENTIAL_PARAMETERS_HH

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

#include "ParameterTraits.hh"

namespace Dune {
  namespace Gesis {

    //! base class for parameter class
    template<class T, class Imp>
    class GeoElectricalPotentialParameterInterface {

    private:
      typedef T Traits;
      typedef typename Traits::RangeFieldType RF;
      //! velocity field

      typedef typename Traits::InputDataType IDT;
      typedef typename Traits::SetupDataType SDT;
      typedef typename Traits::YFieldGeneratorType YFG;

      //typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      typedef typename Traits::BCType BCType;

      Imp& asImp() {
        return static_cast<Imp &> (*this);
      }
      const Imp& asImp() const {
        return static_cast<const Imp &> (*this);
      }

    protected: // visible in the derived classes
      const IDT& inputdata;
      const SDT& setupdata;
      const YFG& fieldgenerator;

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
      enum {
        dim = Traits::dimDomain
      };
      // constructor of base_type:
      GeoElectricalPotentialParameterInterface<T,Imp>( const IDT& inputdata_,
                                                       const SDT& setupdata_,
                                                       const YFG& fieldgenerator_ )
      : inputdata(inputdata_),
        setupdata(setupdata_),
        fieldgenerator(fieldgenerator_)
      {};


      // ! diffusion tensor K(x), groundwater equation: div(-K(x) grad h(x)) = 0
      //
      typename Traits::DiffusionTensorType
      DiffusionTensor( const typename Traits::ElementType& e
                       , const typename Traits::DomainType& x ) const
      {
        typename Traits::DiffusionTensorType K( 0.0 );
        typename Traits::DomainType xg = e.geometry().global(x);
        typename Traits::RangeFieldType fieldvalue = 0.0;

        typename Traits::RangeFieldType v_Well = 0.0; // dummy variable here!

#ifdef DIMENSION3
        fieldgenerator.evaluate3d(xg, fieldvalue, v_Well);
#else
        fieldgenerator.evaluate2d(xg, fieldvalue, v_Well);
#endif
        for( int i=0; i<dim; i++ )
          K[i][i] = std::exp( fieldvalue );

        return K; // conductivity coefficient
      }


      //! Neumann boundary condition (boundary term)
      typename Traits::RangeFieldType
      j( const typename Traits::IntersectionType& is
         , const typename Traits::IntersectionDomainType& x
         ) const {
        //typename Traits::RangeType xg = is.geometry().global(x);
        return 0.0;
      }

      //! reaction term (volume term)
      typename Traits::RangeFieldType
      c( const typename Traits::ElementType& e
         , const typename Traits::DomainType& x) const {
        //typename Traits::RangeType xg = e.geometry().global(x);
        return 0.0; // no reaction term!
      }

      //! source term f (volume term)
      typename Traits::RangeFieldType
      f( const typename Traits::ElementType& e,
         const typename Traits::DomainType& x) const {
        //typename Traits::RangeType xg = e.geometry().global(x);
        return 0.0;
      }


      //! boundary condition type function
      // 0 means Neumann
      // 1 means Dirichlet
      BCType bctype( const typename Traits::IntersectionType& is
                     , const typename Traits::IntersectionDomainType& x
                     ) const {

        typename Traits::RangeType global = is.geometry().global(x);

        // west boundary(0):
        if (global[0] < GEO_EPSILON)
          return convertInt2BCType( setupdata.geoelectrical_potential_equation.boundaries[ 0 ].bctype );

        // east boundary(1):
        if (global[0] > inputdata.domain_data.extensions[0] - GEO_EPSILON)
          return convertInt2BCType( setupdata.geoelectrical_potential_equation.boundaries[ 1 ].bctype );

        // north boundary(2):
        if (global[1] > inputdata.domain_data.extensions[1] - GEO_EPSILON)
          return convertInt2BCType( setupdata.geoelectrical_potential_equation.boundaries[ 2 ].bctype );

        // south boundary(3):
        if (global[1] < GEO_EPSILON)
          return convertInt2BCType( setupdata.geoelectrical_potential_equation.boundaries[ 3 ].bctype );

#ifdef DIMENSION3
        // top boundary(2):
        if (global[2] > inputdata.domain_data.extensions[2] - GEO_EPSILON)
          return convertInt2BCType( setupdata.geoelectrical_potential_equation.boundaries[ 4 ].bctype );

        // bottom boundary(3):
        if (global[2] < GEO_EPSILON)
          return convertInt2BCType( setupdata.geoelectrical_potential_equation.boundaries[ 5 ].bctype );
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
      typename Traits::RangeFieldType g( const typename Traits::ElementType& e,
                                         const typename Traits::DomainType& x) const {
        return asImp().g(e, x);
      }

    };


    /*****************************************************************************
          DERIVED CLASS no.1
    *****************************************************************************/

    //! Parameters for the geo-electrical potential
    template<
      typename GV
      , typename RF
      , typename IDT
      , typename SDT
      , typename YFG
      >
    class GEPotentialForwardProblem
      :
      public GeoElectricalPotentialParameterInterface<
      GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG>,
      GEPotentialForwardProblem<GV,RF,IDT,SDT,YFG>> {

    private:
      typedef GeoElectricalPotentialParameterInterface< GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG>,
                                                        GEPotentialForwardProblem<GV,RF,IDT,SDT,YFG>
                                                        > base_type;
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      typedef GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG> Traits;

      //! constructor
      GEPotentialForwardProblem( const IDT& inputdata_,
                                 const SDT& setupdata_,
                                 const YFG& fieldgenerator_
                                 )
        :
        //pbase_type(dgf_),
        //base_type( pbase_type::data ),
        base_type( inputdata_, setupdata_, fieldgenerator_ ){}


      //! boundary condition type function
      // 0 means Neumann
      // 1 means Dirichlet
      BCType bctype( const typename Traits::IntersectionType& is
                     , const typename Traits::IntersectionDomainType& x
                     ) const {

        typename Traits::RangeType global = is.geometry().global(x);

        // west boundary(0):
        if (global[0] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 0 ].bctype );

        // east boundary(1):
        if (global[0] > base_type::inputdata.domain_data.extensions[0] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 1 ].bctype );

        // north boundary(2):
        if (global[1] > base_type::inputdata.domain_data.extensions[1] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 2 ].bctype );

        // south boundary(3):
        if (global[1] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 3 ].bctype );

#ifdef DIMENSION3
        // top boundary(2):
        if (global[2] > base_type::inputdata.domain_data.extensions[2] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 4 ].bctype );

        // bottom boundary(3):
        if (global[2] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 5 ].bctype );
#endif
        return base_type::convertInt2BCType(1); // undefined case: default = Dirichlet
      };


      //! convectiondiffusionparameter.hh requires g()
      typename Traits::RangeFieldType
        g( const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
        return g_function(e,x);
      };

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
        g_function(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {

        typename Traits::RangeType global = e.geometry().global(x);

        typename Traits::RangeFieldType y = 0.0;

        // west boundary(0):
        if (global[0] < GEO_EPSILON)
          y = base_type::setupdata.geoelectrical_potential_equation.boundaries[0].stripes[0].value;

        // east boundary(0):
        if( global[0] > base_type::inputdata.domain_data.extensions[0] - GEO_EPSILON )
          y = base_type::setupdata.geoelectrical_potential_equation.boundaries[1].stripes[0].value;


        return y;

      };

    };




    /*****************************************************************************
          DERIVED CLASS no.2
    *****************************************************************************/


    //! Transport parameters class for adjoint states
    template<
      typename GV
      , typename RF
      , typename IDT
      , typename SDT
      , typename YFG
      >
    class GEPotentialAdjointProblem
      :
      public GeoElectricalPotentialParameterInterface<
      GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG>,
      GEPotentialAdjointProblem<GV,RF,IDT,SDT,YFG>
      > {

    private:
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      typedef GeoElectricalPotentialParameterInterface<
        GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG>,
        GEPotentialAdjointProblem<GV,RF,IDT,SDT,YFG>
        > base_type;

    public:
      typedef GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG> Traits;

      // constructor:
      GEPotentialAdjointProblem( const IDT& inputdata_,
                                 const SDT& setupdata_,
                                 const YFG& fieldgenerator_ )
        :
        base_type( inputdata_,
                   setupdata_,
                   fieldgenerator_ ){}

      //! boundary condition type function
      // 0 means Neumann
      // 1 means Dirichlet
      BCType bctype( const typename Traits::IntersectionType& is
                     , const typename Traits::IntersectionDomainType& x
                     ) const {

        typename Traits::RangeType global = is.geometry().global(x);

        // west boundary(0):
        if (global[0] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 0 ].bctype );

        // east boundary(1):
        if (global[0] > base_type::inputdata.domain_data.extensions[0] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 1 ].bctype );

        // north boundary(2):
        if (global[1] > base_type::inputdata.domain_data.extensions[1] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 2 ].bctype );

        // south boundary(3):
        if (global[1] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 3 ].bctype );

#ifdef DIMENSION3
        // top boundary(2):
        if (global[2] > base_type::inputdata.domain_data.extensions[2] - GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 4 ].bctype );

        // bottom boundary(3):
        if (global[2] < GEO_EPSILON)
          return base_type::convertInt2BCType( base_type::setupdata.geoelectrical_potential_equation.boundaries[ 5 ].bctype );
#endif
        return base_type::convertInt2BCType(1); // undefined case: default = Dirichlet

      };


      //! convectiondiffusionparameter.hh requires g()
      typename Traits::RangeFieldType
      g( const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
        return g_function(e,x);
      }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
      g_function( const typename Traits::ElementType& e,
                  const typename Traits::DomainType& x) const {
        return 0.0;
      }

    };

  } // namespace PDELab

} // namespace Dune

#endif // DUNE_GESIS_GEO_ELECTRICAL_POTENTIAL_PARAMETERS_HH
