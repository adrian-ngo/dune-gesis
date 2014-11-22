#ifndef DUNE_GESIS_SOLUTE_PUMPING_SOURCE_HH
#define DUNE_GESIS_SOLUTE_PUMPING_SOURCE_HH

#include "SourceBaseClass.hh"

namespace Dune {
  namespace Gesis {


    //======================================================
    // The template class to define the solute source term
    // for forward solute concentration simulations
    //======================================================

    template<typename GV,
             typename RF,
             typename IDT,
             typename SDT>
    class SolutePumpingSource
      : public SourceTermInterface<
      SourceTermTraits<GV,RF>,
      SolutePumpingSource<GV,RF,IDT,SDT>
      >
    {
    private:
      enum{ dim = GV::dimension };
      const IDT& inputdata;
      const SDT& setupdata;

      inline void zone_porosity(const REAL& zone_coord, REAL & porosity) const{

        UINT inside=0;
        // only one zone!
        if(inputdata.yfield_properties.nz<2){
          porosity=inputdata.yfield_properties.zones[0].porosity;
          return;
        }

        for(UINT ii=0; ii<inputdata.yfield_properties.nz-1; ii++){
          if(zone_coord<=inputdata.yfield_properties.zones[ii].bottom ){
            porosity=inputdata.yfield_properties.zones[ii].porosity;
            inside=1;
            break;
          }
        }
        //last zone!
        if(zone_coord>inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-2].bottom && inside==0){
          porosity=inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].porosity;
        }

      }



      template<typename DFV0>   // Dune::FieldVector of codim 0
      inline bool isPointInsideReachOfWell( const DFV0& elementpoint,
                                            const RF& reach,
#ifdef DIMENSION3
                                            const RF& reach_y,
#endif
                                            const int& iWell
                                            ) const {
        // extract well-data:
        RF well_position_x =  setupdata.wdlist.pointdata_vector[iWell].x;
#ifdef DIMENSION3
        RF well_position_y =  setupdata.wdlist.pointdata_vector[iWell].y;
#endif
        RF well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
        RF well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;

        if(
           elementpoint[0] >= well_position_x - reach
           &&
           elementpoint[0] < well_position_x + reach
           &&
#ifdef DIMENSION3
           elementpoint[1] >= well_position_y - reach_y
           &&
           elementpoint[1] < well_position_y + reach_y
           &&
           elementpoint[2] > well_bottom - GEO_EPSILON
           &&
           elementpoint[2] < well_top + GEO_EPSILON
#else
           elementpoint[1] > well_bottom - GEO_EPSILON
           &&
           elementpoint[1] < well_top + GEO_EPSILON
#endif
           )
          return true;
        else
          return false;
      }

    public:

      const SourceNature source_nature;
      typedef SourceTermTraits<GV,RF> Traits;
      SolutePumpingSource(
                          const IDT& inputdata_,
                          const SDT& setupdata_
                          )
        :
        inputdata(inputdata_),
        setupdata(setupdata_),
        source_nature( SOLUTE_PUMPING_SOURCES )
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
                                         ) const
      {
        UINT nSources=setupdata.pdlist.total;  // adding Point sources only!

        for( UINT i=0; i<nSources; i++ ) {

          // Do something only if the pumping source is non-zero.

          RF pumping_rate = setupdata.pdlist.pointdata_vector[i].value;
          RF in_concentration = setupdata.pdlist.pointdata_vector[i].concentration;
          RF injection_time = setupdata.pdlist.pointdata_vector[i].concentration_injection_time;

          if(injection_time<=GEO_EPSILON)
            injection_time=1.0;


          if(  pumping_rate  > GEO_EPSILON*0.5&& fabs(in_concentration) >  GEO_EPSILON*0.5) {
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

              lfsv.finiteElement().localBasis().evaluateFunction( pointsource_local, shapefunctions );

              RF factor=eg.geometry().volume();
              RF npoints=eg.geometry().corners();

              for( UINT ii=0; ii<lfsv.size(); ii++ ) {

                residual.accumulate( lfsv, ii, (-pumping_rate/factor*npoints)*in_concentration*injection_time* shapefunctions[ii]);

              }
            }

          }

        } // end for

        return true;
      }  // bool evaluate_residual_on_element()


    }; // class SolutePumpingSource



  }
}

#endif
