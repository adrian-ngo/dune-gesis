#ifndef DUNE_GESIS_PUMPING_SOURCE_HH
#define DUNE_GESIS_PUMPING_SOURCE_HH

#include "SourceBaseClass.hh"

namespace Dune {

  namespace Gesis {


    //=========================================================================
    // The template class to define the source term for forward flow simulations
    //========================================================================

    template< typename GV,
              typename RF,
              typename SDT
              >
    class PumpingSource
      : public SourceTermInterface<
      SourceTermTraits<GV,RF>,
      PumpingSource<GV,RF,SDT>
      >
    {
    private:
      enum{ dim = GV::dimension };
      const SDT& setupdata;

    public:

      const SourceNature source_nature;
      typedef SourceTermTraits<GV,RF> Traits;
      PumpingSource( const SDT& setupdata_ )
        :
        setupdata(setupdata_),
        source_nature( PUMPING_SOURCES )
      {}

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

        UINT nSources=setupdata.pdlist.total;

        for( UINT i=0; i<nSources; i++ ) {

          RF source_value=setupdata.pdlist.pointdata_vector[i].value;

          if( std::abs( source_value ) > GEO_EPSILON*0.5 ) {

            Dune::FieldVector<RF,dim> pointsource_global;

            pointsource_global[0] = setupdata.pdlist.pointdata_vector[i].x;
            pointsource_global[1] = setupdata.pdlist.pointdata_vector[i].y;
#ifdef DIMENSION3
            pointsource_global[2] = setupdata.pdlist.pointdata_vector[i].z;
#endif

            // get the local coordinate of all the source locations:
            Dune::FieldVector<RF,dim> pointsource_local
              = eg.geometry().local( pointsource_global );

            // Check if the point source lies inside the rectangular cell.
            UINT iOutside = 0;
            for( UINT ii=0; ii<dim; ii++ ){
              if( (pointsource_local[ii] < 0.0) || (pointsource_local[ii] >= 1.0) ){
                iOutside++;
              }
            }

            if( iOutside==0 ) {
              //logger << pointsource_global << "FLOW"<<std::endl;
              //logger << "element center: " <<eg.geometry().center()<<std::endl;
              //logger << "pointsource global: " <<pointsource_global <<std::endl;
              // pointsource is inside the rectangular cell!

              lfsv.finiteElement().localBasis().evaluateFunction( pointsource_local, shapefunctions );

              for( UINT ii=0; ii<lfsv.size(); ii++ )
                {
                  //std::cout << " i = " << i;
                  //std::cout << " r = " << residual[i];
                  residual.accumulate( lfsv, ii, - source_value * shapefunctions[ii] );
                  //std::cout << " --> " << residual[i] << std::endl;
                }
            }
          }
        }

        return true;

      } // bool evaluate_residual_on_element()

    }; // class PumpingSource


  }
}

#endif
