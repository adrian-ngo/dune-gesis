#ifndef DUNE_GESIS_GEOELECTRICAL_POTENTIAL_SOURCE_HH
#define DUNE_GESIS_GEOELECTRICAL_POTENTIAL_SOURCE_HH

#include "SourceBaseClass.hh"

namespace Dune {
  namespace Gesis {


    //==========================================================
    // The template class to define the source term
    // for forward flow simulations
    //==========================================================

    template<
      typename GV,
      typename RF,
      typename SDT
      >
    class GeoelectricalPotentialSource
      : public SourceTermInterface<
      SourceTermTraits<GV,RF>,
      GeoelectricalPotentialSource<GV,RF,SDT>
      >
    {
    private:

      enum{ dim = GV::dimension };

      const SDT& setupdata;
      UINT config;


    public:
      const SourceNature source_nature;
      typedef SourceTermTraits<GV,RF> Traits;

      // The constructor:
      // The 2nd constructor parameter is used to indicate a negative point-source for the adjoint TPE
      GeoelectricalPotentialSource( const SDT& setupdata_,
                                    const UINT & config_ )
        :
        setupdata(setupdata_),
        config(config_),
        source_nature( GEOELECTRICAL_POTENTIAL_SOURCE )
      {
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
                                         ) const
      {

        //pumping HERE!
        for(UINT iSource=0; iSource<setupdata.geoelectrical_potential_inversion_data.inoutlist[config].pointdata_vector.size(); iSource++){

          RF source_value=setupdata.geoelectrical_potential_inversion_data.inoutlist[config].pointdata_vector[iSource].value;

          if( std::abs( source_value ) > GEO_EPSILON*0.5 ) {
            Dune::FieldVector<RF,dim> in_global;
            Dune::FieldVector<RF,dim> out_global;

            in_global[0] = setupdata.geoelectrical_potential_inversion_data.inoutlist[config].pointdata_vector[iSource].x;
            in_global[1] = setupdata.geoelectrical_potential_inversion_data.inoutlist[config].pointdata_vector[iSource].y;
            out_global[0] = setupdata.geoelectrical_potential_inversion_data.inoutlist[config].pointdata_vector[iSource].x2;
            out_global[1] = setupdata.geoelectrical_potential_inversion_data.inoutlist[config].pointdata_vector[iSource].y2;
#ifdef DIMENSION3
            in_global[2] = setupdata.geoelectrical_potential_inversion_data.inoutlist[config].pointdata_vector[iSource].z;
            out_global[2] = setupdata.geoelectrical_potential_inversion_data.inoutlist[config].pointdata_vector[iSource].z2;
#endif

            // get the local coordinate of the injection of the current
            Dune::FieldVector<RF,dim> local = eg.geometry().local( in_global );
            // Check if the point source lies inside the rectangular cell.

            UINT iOutside = 0;
            for( UINT ii=0; ii<dim; ii++ ){
              if( (local[ii] < 0.0) || (local[ii] >= 1.0) ){
                iOutside++;
              }
            }

            if( iOutside==0 ) {

              lfsv.finiteElement().localBasis().evaluateFunction( local, shapefunctions );
              for( UINT i=0; i<lfsv.size(); i++ )
                {
                  residual.accumulate( lfsv, i, - source_value * shapefunctions[i] );
                }
            }


            // get the local coordinate of the extraction of the current
            local = eg.geometry().local( out_global );
            iOutside=0;
            for( UINT ii=0; ii<dim; ii++ ){
              if( (local[ii] < 0.0) || (local[ii] >= 1.0) ){
                iOutside++;
              }
            }

            if( iOutside==0 ) {

              lfsv.finiteElement().localBasis().evaluateFunction( local, shapefunctions );
              for( UINT ii=0; ii<lfsv.size(); ii++ )
                {
                  residual.accumulate( lfsv, ii,  source_value * shapefunctions[ii] );
                }
            }
          }
        }

        return true;
      } // bool evaluate_residual_on_element()

    }; // class GeoelectricalPotentialSource


  }
}

#endif
