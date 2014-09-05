#ifndef HEAT_PUMPING_SOURCE_HH
#define HEAT_PUMPING_SOURCE_HH

#include "SourceBaseClass.hh"

namespace Dune {

  namespace GeoInversion {

    //=================================================================================
    // The template class to define the HEAT source term for forward HEAT simulations
    //=================================================================================
	
    template<typename GV,
             typename RF,
             typename IDT,
             typename SDT>
    class HeatPumpingSource
      : public SourceTermInterface<
      SourceTermTraits<GV,RF>,
      HeatPumpingSource<GV,RF,IDT,SDT>
      > 
    {
    private:
      enum{ dim = GV::dimension };
      const IDT& inputdata;
      const SDT& setupdata;
      
      inline void zone_parameters(const REAL& zone_coord, REAL & porosity, REAL & rho_s, REAL & c_s) const{
	
        UINT inside=0;
        // only one zone!
        if(inputdata.yfield_properties.nz<2){
          porosity=inputdata.yfield_properties.zones[0].porosity;
          rho_s=inputdata.yfield_properties.zones[0].rho;
          c_s=inputdata.yfield_properties.zones[0].c_s;
          return; 
        }
	
        for(UINT ii=0; ii<inputdata.yfield_properties.nz-1; ii++){
          if(zone_coord<=inputdata.yfield_properties.zones[ii].bottom ){
            porosity=inputdata.yfield_properties.zones[ii].porosity;
            rho_s=inputdata.yfield_properties.zones[ii].rho;
            c_s=inputdata.yfield_properties.zones[ii].c_s;
            inside=1;
            break;
          }
        }
        //last zone!
        if(zone_coord>inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-2].bottom && inside==0){
          porosity=inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].porosity;
          rho_s=inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].rho;
          c_s=inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].c_s;
        }

      } 
      

      template<typename DFV0>   // DFV0 = Dune::FieldVector of codim 0
      inline bool isPointInsideReachOfWell( const DFV0& elementpoint, 
                                            const RF& reach,
#ifdef DIMENSION3
                                            const RF& reach_y,
#endif
                                            const int& iWell,
                                            const int& Setup
                                            ) const {

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
      HeatPumpingSource( const IDT& inputdata_, const SDT& setupdata_ )
        : 
        inputdata(inputdata_),
        setupdata(setupdata_),
        source_nature( HEAT_PUMPING_SOURCES )
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
        
        UINT nSources=setupdata.pdlist.total;
		  
        RF rho_w = inputdata.transport_parameters.heat_rho_w;
        RF c_w   = inputdata.transport_parameters.heat_c_w;
        
        for( UINT i=0; i<nSources; i++ ) {
          RF pumping_rate, temperature,injection_time, rho_s,c_s,porosity;
          // Do something only if the pumping source is non-zero.

          pumping_rate = setupdata.pdlist.pointdata_vector[i].value;
          temperature = setupdata.pdlist.pointdata_vector[i].temperature;
          injection_time = setupdata.pdlist.pointdata_vector[i].temperature_injection_time;
          if(injection_time<=GEO_EPSILON)
            injection_time=1.0;  // TODO: Ask Ronnie! Why this?

          if(  pumping_rate > GEO_EPSILON*0.5 && fabs(temperature) > GEO_EPSILON*0.5) {
            Dune::FieldVector<RF,dim> pointsource_global;
              
            pointsource_global[0] = setupdata.pdlist.pointdata_vector[i].x;
            pointsource_global[1] = setupdata.pdlist.pointdata_vector[i].y;
#ifdef DIMENSION3
            pointsource_global[2] = setupdata.pdlist.pointdata_vector[i].z;
            zone_parameters(pointsource_global[2], porosity, rho_s, c_s);
#else
            zone_parameters(pointsource_global[1], porosity, rho_s, c_s);
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
              RF rho_m_c_m;
                
              rho_m_c_m=(porosity*rho_w*c_w)+((1.0-porosity)*rho_s*c_s);
                
              RF factor=eg.geometry().volume();
              RF npoints=eg.geometry().corners();
                
              for( UINT ii=0; ii<lfsv.size(); ii++ ) {
					    
                residual.accumulate( lfsv, ii, (-pumping_rate/factor*npoints)*temperature*injection_time*(rho_w*c_w/rho_m_c_m)* shapefunctions[ii]);
                //residual.accumulate( lfsv, i, (- pumping_rate)*temperature*(/*(rho_w*c_w)*/1.0/(rho_w*c_w))* shapefunctions[i] );
                //std::cout << " --> " << residual[i] << std::endl;
              }
            }
          }
        }

#ifdef WELL_FRACTURE_MODEL
        nSources=setupdata.wdlist.total;

        for( UINT i=0; i<nSources; i++ ) {
          Dune::FieldVector<RF,dim> global=eg.geometry().center();
          
          std::vector<REAL> meshsize(inputdata.domain_data.yasp_gridsizes);
          std::vector<REAL> meshsize_2(inputdata.domain_data.yasp_gridsizes);
          for(UINT ii=0; ii<meshsize_2.size();ii++)
            meshsize_2[ii]*=0.5;
          
          RF well_x = setupdata.wdlist.pointdata_vector[i].x;
#ifdef DIMENSION3     
          RF well_y=setupdata.wdlist.pointdata_vector[i].y;
#endif
          RF well_top = setupdata.wdlist.pointdata_vector[i].well_top;
          RF well_bottom = setupdata.wdlist.pointdata_vector[i].well_bottom;

          if( global[0] >= well_x-meshsize_2[0] &&  global[0] < well_x+meshsize_2[0] 
#ifdef DIMENSION3 			  
              && global[1] >= well_y-meshsize_2[1] &&  global[1] < well_y+meshsize_2[1] 
              && global[2] >= well_bottom && global[2] < well_top
#else
              && global[1] >= well_bottom && global[1] < well_top
#endif
              ) {
            
            RF rho_s,c_s,porosity;
            RF pumping_rate = setupdata.wdlist.pointdata_vector[i].well_rate;
            RF temperature = setupdata.wdlist.pointdata_vector[i].temperature;
            RF injection_time = setupdata.wdlist.pointdata_vector[i].temperature_injection_time;
            if(injection_time<=GEO_EPSILON)
              injection_time=1.0;
            
            if( pumping_rate > GEO_EPSILON*0.5 && fabs(temperature) >  GEO_EPSILON*0.5) {
              Dune::GeometryType gt = eg.geometry().type();
#ifdef DIMENSION3                       
              zone_parameters(global[2], porosity, rho_s, c_s);
#else
              zone_parameters(global[1], porosity, rho_s, c_s);
#endif
              RF rho_m_c_m;
              
				
              rho_m_c_m=(porosity*rho_w*c_w)+((1.0-porosity)*rho_s*c_s);

              std::vector<Dune::FieldVector<RF,dim>> well_points;
              std::vector<size_t> well_points_cornerindex;
	
              size_t nElementCorners = eg.geometry().corners();
              for( size_t iElementCorner = 0; iElementCorner < nElementCorners; iElementCorner++ ) {
                Dune::FieldVector<RF,dim> elementcorner = eg.geometry().corner( iElementCorner );
#ifdef DIMENSION3
                if( isPointInsideReachOfWell( elementcorner,meshsize_2[0], meshsize_2[1], i ) )
#else
                  if( isPointInsideReachOfWell( elementcorner, meshsize_2[0], i ) )
#endif
                    {
                      well_points.push_back( elementcorner );
                      well_points_cornerindex.push_back( iElementCorner );
                      
                    }
              }
				  
              for(UINT iEdge=0; iEdge<well_points.size(); iEdge++) {
                Dune::FieldVector<RF,dim> local=eg.geometry().local( well_points[iEdge]);
                lfsv.finiteElement().localBasis().evaluateFunction( local, shapefunctions);

                for( UINT ii=0; ii<lfsv.size(); ii++ ) {
                  residual.accumulate( lfsv, ii,  (-pumping_rate)*temperature*injection_time*(rho_w*c_w/rho_m_c_m)* shapefunctions[ii]);

                }

              }
              
            }

          }

        }
#endif // WELL_FRACTURE_MODEL
        return true;
      } // bool evaluate_residual_on_element()

    }; // class HeatPumpingSource
    
  }
}

#endif // HEAT_PUMPING_SOURCE_HH
