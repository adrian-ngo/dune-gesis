#ifndef DUNE_GESIS_GEP_POINT_SOURCE_HH
#define DUNE_GESIS_GEP_POINT_SOURCE_HH

#include "SourceBaseClass.hh"

namespace Dune {
  namespace Gesis {

    //=========================================================================
    // The template class to define the source term for 
    // adjoint geoelectrical potential simulations
    //========================================================================

    template<typename GV,
             typename RF,
             typename IDT,
             bool take_cell_average = false
             >
    class GeoelectricalPointSource 
      : public SourceTermInterface<
      SourceTermTraits<GV,RF>,
      GeoelectricalPointSource<GV,RF,IDT,take_cell_average>
      >
    {
    private:

      enum{dim=GV::dimension};

      Dune::FieldVector<RF,dim> pointsource_location1_global;
      Dune::FieldVector<RF,dim> pointsource_location2_global;





    public:
      
      const SourceNature source_nature;
	  
      typedef SourceTermTraits<GV,RF> Traits;
	  
      // The constructor:
      // The 2nd constructor parameter is used to indicate a negative point-source for the adjoint TPE
      GeoelectricalPointSource()
        : 
        pointsource_location1_global( 0.0 ),
        pointsource_location2_global( 0.0 ),
        source_nature( GEOELECTRICAL_POINT_SOURCE )
      {}
      



      template<typename COORDINATES>
      void set_PointSourceLocations( const COORDINATES& pointsource1_global, const COORDINATES& pointsource2_global )
      {
        pointsource_location1_global[ 0 ] = pointsource1_global[0];
        pointsource_location1_global[ 1 ] = pointsource1_global[1];
        pointsource_location2_global[ 0 ] = pointsource2_global[0];
        pointsource_location2_global[ 1 ] = pointsource2_global[1];
#ifdef DIMENSION3
        pointsource_location1_global[ 2 ] = pointsource1_global[2];
        pointsource_location2_global[ 2 ] = pointsource2_global[2];
#endif
        //std::cout << "set pointsource_location_global = " << pointsource_location_global << " ... ";
      }


      template< typename EG
                , typename LFSV
                , typename SF
                , typename R
                >
      bool evaluate_residual_on_element( const EG& eg
                                         , const LFSV& lfsv
                                         , SF& shapefunctions
                                         , R& residual
                                         , const bool bAdjoint=false
                                         ) const {

        Dune::FieldVector<RF,dim> pointsource_local = eg.geometry().local( pointsource_location1_global );

        UINT iOutside = 0;
        for( UINT i=0; i<dim; i++ ){
          if( (pointsource_local[i] < 0.0 ) ||
              (pointsource_local[i] >= 1.0) ){
            iOutside++;
          }
        }
			
        if( iOutside==0 ){
          //std::cout << pointsource_location_global << " ... ";
          // pointsource is inside the rectangular cell!
          lfsv.finiteElement().localBasis().evaluateFunction( pointsource_local, shapefunctions );
          for( UINT i=0; i<lfsv.size(); i++ ){
            //std::cout << " i = " << i;
            //std::cout << " r = " << residual[i];
            residual.accumulate( lfsv, i, -shapefunctions[i] );
            //std::cout << " --> " << residual[i] << std::endl;
          }
        }
	
        pointsource_local = eg.geometry().local( pointsource_location2_global );
			
        // Code by Olaf Cirpka to evaluate shape functions at a single point source!
        // Check if the point source lies inside the rectangular cell.
		
        iOutside = 0;
        for( UINT i=0; i<dim; i++ ){
          if( (pointsource_local[i] < 0.0 ) ||
              (pointsource_local[i] >= 1.0) ){
            iOutside++;
          }
        }
			
        if( iOutside==0 ){
          //std::cout << pointsource_location_global << " ... ";
          // pointsource is inside the rectangular cell!
          lfsv.finiteElement().localBasis().evaluateFunction( pointsource_local, shapefunctions );
          for( UINT i=0; i<lfsv.size(); i++ ){
            //std::cout << " i = " << i;
            //std::cout << " r = " << residual[i];
            residual.accumulate( lfsv, i, shapefunctions[i] );
            //std::cout << " --> " << residual[i] << std::endl;
          }
        }
	
        return true;

      } // bool evaluate_residual_on_element()


    }; // class GeoelectricalPointSource


  }
}

#endif
