// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_GESIS_GROUNDWATER_PARAMETERS_HH
#define DUNE_GESIS_GROUNDWATER_PARAMETERS_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

//#include<dune/grid/common/genericreferenceelements.hh>
//#include<dune/grid/common/quadraturerules.hh>

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
//#include"dune/pdelab/localoperator/idefault.hh"  // instationary problems

//#include "inputfile.hh"
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include "ParameterTraits.hh"

namespace Dune {
  namespace Gesis {


    //! base class for parameter class
    template<class T, class Imp>
    class GroundwaterParameterInterface {
      
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
      const YFG& yfieldgenerator;
      
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
      GroundwaterParameterInterface<T,Imp>( const IDT& inputdata_, 
                                            const SDT& setupdata_,
                                            const YFG& yfieldgenerator_ )
      : inputdata(inputdata_),
        setupdata(setupdata_),
        yfieldgenerator(yfieldgenerator_)
      {};


      // ! diffusion tensor K(x), groundwater equation: div(-K(x) grad h(x)) = 0
      //
      typename Traits::DiffusionTensorType
      DiffusionTensor( const typename Traits::ElementType& e
                       , const typename Traits::DomainType& x ) const 
      {
        typename Traits::DiffusionTensorType K( 0.0 );
        typename Traits::DomainType xg = e.geometry().global(x);

        REAL k = 0.0;
        REAL k_Well = 0.0;


#ifdef DIMENSION3
        yfieldgenerator.evaluate3d( xg, k, k_Well );
        K[0][0] = k;
        K[1][1] = k;
        K[2][2] = k_Well;
#else
        yfieldgenerator.evaluate2d( xg, k, k_Well );
        K[0][0] = k;
        K[1][1] = k_Well;
#endif

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
          return convertInt2BCType( setupdata.flow_equation.boundaries[ 0 ].bctype );
          
        // east boundary(1):
        if (global[0] > inputdata.domain_data.extensions[0] - GEO_EPSILON)
          return convertInt2BCType( setupdata.flow_equation.boundaries[ 1 ].bctype );
          
        // south boundary(2):
        if (global[1] < GEO_EPSILON)
          return convertInt2BCType( setupdata.flow_equation.boundaries[ 2 ].bctype );

        // north boundary(3):
        if (global[1] > inputdata.domain_data.extensions[1] - GEO_EPSILON)
          return convertInt2BCType( setupdata.flow_equation.boundaries[ 3 ].bctype );
            
#ifdef DIMENSION3
        // bottom boundary(4):
        if (global[2] < GEO_EPSILON)
          return convertInt2BCType( setupdata.flow_equation.boundaries[ 4 ].bctype );

        // top boundary(5):
        if (global[2] > inputdata.domain_data.extensions[2] - GEO_EPSILON)
          return convertInt2BCType( setupdata.flow_equation.boundaries[ 5 ].bctype );
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

    //! Transport parameters class for m0
    template<
      typename GV
      , typename RF
      , typename IDT
      , typename SDT
      , typename YFG
      >
    class GroundwaterForwardProblem
      :
      public GroundwaterParameterInterface<
      GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG>,
      GroundwaterForwardProblem<GV,RF,IDT,SDT,YFG>> {

    private:
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

      typedef GroundwaterParameterInterface<
      GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG>,
        GroundwaterForwardProblem<GV,RF,IDT,SDT,YFG> > base_type;
    
    public:
      typedef GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG> Traits;

      //! constructor
      GroundwaterForwardProblem( 
                                const IDT& inputdata_,
                                const SDT& setupdata_,
                                const YFG& yfieldgenerator_
                                 )
        :
        //pbase_type(dgf_),
        //base_type( pbase_type::data ),
        base_type( inputdata_, setupdata_, yfieldgenerator_ ){}

      //! convectiondiffusionparameter.hh requires g()
      typename Traits::RangeFieldType
      g( const typename Traits::ElementType& e, const typename Traits::DomainType& x) const {
        return g_function(e,x);
      }

      //! Dirichlet boundary condition value
      typename Traits::RangeFieldType
        g_function( const typename Traits::ElementType& e, const typename Traits::DomainType& x ) const {

        typename Traits::RangeType global = e.geometry().global(x);

        typename Traits::RangeFieldType y = 0.0;
        
        // west boundary(0):
        if (global[0] < GEO_EPSILON){
          y = base_type::setupdata.flow_equation.boundaries[0].stripes[0].value; // default value

          for(UINT ii=0; ii < base_type::setupdata.flow_equation.boundaries[0].stripes.size(); ii++){          
            if( global[1] >= base_type::setupdata.flow_equation.boundaries[0].stripes[ii].startpoint[0] - GEO_EPSILON 
                && 
                global[1] <= base_type::setupdata.flow_equation.boundaries[0].stripes[ii].endpoint[0] + GEO_EPSILON
#ifdef DIMENSION3
                && 
                global[2] >= base_type::setupdata.flow_equation.boundaries[0].stripes[ii].startpoint[1] - GEO_EPSILON 
                && 
                global[2] <= base_type::setupdata.flow_equation.boundaries[0].stripes[ii].endpoint[1] + GEO_EPSILON
#endif	    
                ){
              y = base_type::setupdata.flow_equation.boundaries[0].stripes[ii].value;
            }
          }
          if( std::abs(y)>1E-12 )
            return y; // return before it gets overridden!
        }

        // east boundary(0):

#ifdef OLES_TEST
        if( global[0] > base_type::inputdata.domain_data.extensions[0] - GEO_EPSILON ) {
          REAL yNorth = base_type::setupdata.flow_equation.boundaries[3].stripes[0].value;
          REAL ySouth = base_type::setupdata.flow_equation.boundaries[2].stripes[0].value;
          y = ySouth + global[1] * ( yNorth - ySouth ) / ( base_type::inputdata.domain_data.extensions[1] );
          if( std::abs(y)>1E-12 )
            return y; // return before it gets overridden!
        }
#else
        if( global[0] > base_type::inputdata.domain_data.extensions[0] - GEO_EPSILON ) {
          y = base_type::setupdata.flow_equation.boundaries[1].stripes[0].value;
          if( std::abs(y)>1E-12 )
            return y; // return before it gets overridden!
        }
#endif
        
        // south boundary(0):
        if (global[1] < GEO_EPSILON){
          y = base_type::setupdata.flow_equation.boundaries[2].stripes[0].value;
          if( std::abs(y)>1E-12 )
            return y; // return before it gets overridden!
        }

        // north boundary(0):
        if( global[1] > base_type::inputdata.domain_data.extensions[1] - GEO_EPSILON ){
          y = base_type::setupdata.flow_equation.boundaries[3].stripes[0].value;
          if( std::abs(y)>1E-12 )
            return y; // return before it gets overridden!
        }

        return y;

      }

      const YFG& getYfieldgenerator() const {
        return base_type::yfieldgenerator;
      }
      
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
    class GroundwaterAdjointProblem
      :
      public GroundwaterParameterInterface<
      GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG>,
      GroundwaterAdjointProblem<GV,RF,IDT,SDT,YFG>
      > {
    
    private:
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      typedef GroundwaterParameterInterface<
        GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG>,
        GroundwaterAdjointProblem<GV,RF,IDT,SDT,YFG>
        > base_type;

    public:
      typedef GroundwaterParameterTraits<GV,RF,IDT,SDT,YFG> Traits;

      // constructor:
      GroundwaterAdjointProblem( const IDT& inputdata_,
                                 const SDT& setupdata_,
                                 const YFG& yfieldgenerator_ )
        :
        base_type( inputdata_, 
                   setupdata_,
                   yfieldgenerator_ ){}



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

#endif // DUNE_GESIS_GROUNDWATER_PARAMETERS_HH
