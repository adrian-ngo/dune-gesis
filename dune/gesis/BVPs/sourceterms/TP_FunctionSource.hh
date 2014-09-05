#ifndef NULL_SOURCE_HH
#define NULL_SOURCE_HH

#include "SourceBaseClass.hh"

namespace Dune {

  namespace GeoInversion {


    //===================================================================================
    // The template class to define the functional source term for transport simulations 
    // (needed for adjoint transport for the sensitivity of geoelectrical potential)
    //===================================================================================
    
    template<
      typename GV
      , typename GWP
      , typename GFS
      , typename RF
      , typename VCType
      , typename IDT
      , typename SDT
      >
    class TP_FunctionSource
      : public SourceTermInterface<
      SourceTermTraits<GV,RF>,
      TP_FunctionSource<GV,GWP,GFS,RF,VCType,IDT,SDT>
      >
    {
    private:
      enum{ dim = GV::dimension };

      const GWP& gwp;
      const GFS& gfs;
      const IDT& inputdata;

      VCType vc_source_1;
      VCType vc_source_2;

    public:

      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS,VCType> GRADIENT_DGF;
      typedef Dune::PDELab::DiscreteGridFunctionDarcy<GWP,GFS> DARCY_FLUX_DGF;


      const SourceNature source_nature;
      typedef SourceTermTraits<GV,RF> Traits;

      // constructor:
      TP_FunctionSource( const GWP& gwp_,
                         const GFS& gfs_,
                         const IDT& inputdata_ )
        :
        gwp(gwp_),
        gfs( gfs_ ),
        inputdata( inputdata_ ),
        vc_source_1( gfs, 0.0 ), 
        vc_source_2( gfs, 0.0 ),
        source_nature( TP_FUNCTIONAL_SOURCE )
      {
        if(inputdata.verbosity>8){
          std::cout << "TP_FunctionSource constructor" << std::endl;
        }
      }
      




      void reset_rhs() {
        vc_source_1=VCType(gfs,0.0);
        vc_source_2=VCType(gfs,0.0);
      }
      

      void set_rhs( const VCType& vc_1, const VCType& vc_2 ) {
        vc_source_1 = vc_1;
        vc_source_2 = vc_2;
      }


      template<typename EG, typename COORDINATES>
      bool evaluate_function( const EG& e
                              , const COORDINATES& xlocal
                              , REAL& fval
                              ) const {
		
        // Hint: 
        // The DGF is being created out of the VCType everytime.
        // I don't know yet whether this will have a negative effect on the performance.
        // I hope that this construtor does not do much work here!
        GRADIENT_DGF dgf_1(gfs, vc_source_1);
		
        DARCY_FLUX_DGF dgf_2( gwp, gfs, vc_source_2, inputdata,-1,-1,-1,-1,false );

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

    }; // class TP_FunctionSource


  }
}

#endif
