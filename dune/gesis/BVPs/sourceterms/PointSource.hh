#ifndef DUNE_GESIS_POINT_SOURCE_HH
#define DUNE_GESIS_POINT_SOURCE_HH

#include "SourceBaseClass.hh"

namespace Dune {
  namespace Gesis {

    //=========================================================
    // The template class to define the source term 
    // for forward flow simulations
    //=========================================================
    
    template<typename GV,
             typename RF,
             typename IDT,
             bool take_cell_average = false
             >
    class PointSource 
      : public SourceTermInterface<
      SourceTermTraits<GV,RF>,
      PointSource<GV,RF,IDT,take_cell_average>
      >
    {
    private:
      
      enum{dim=GV::dimension};

      const IDT& inputdata;

      Dune::FieldVector<RF,dim> pointsource_location_global;
      REAL sign_factor; // must be == 1.0 for the adjoint GWE and == -1.0 for the adjoint TPE
      REAL SV_y;
#ifdef DIMENSION3 
      REAL SV_z;
#endif

    public:
      
      const SourceNature source_nature;
      
      typedef SourceTermTraits<GV,RF> Traits;
      
      // The constructor:
      // The 2nd constructor parameter is used to indicate a negative point-source for the adjoint TPE
      PointSource( const IDT& inputdata_,
                   const REAL sign_factor_=1.0 )
        : 
        inputdata(inputdata_),
        pointsource_location_global( 0.0 ),
        sign_factor( sign_factor_ ),
        SV_y( 0.0 ),
#ifdef DIMENSION3 
        SV_z( 0.0 ),
#endif
        source_nature( POINT_SOURCE )
      {}
      

      template<typename COORDINATES>
      void set_PointSourceLocation( const COORDINATES& pointsource_global )
      {
        pointsource_location_global[ 0 ] = pointsource_global[0];
        pointsource_location_global[ 1 ] = pointsource_global[1];
#ifdef DIMENSION3
        pointsource_location_global[ 2 ] = pointsource_global[2];
#endif
        //std::cout << "set pointsource_location_global = " << pointsource_location_global << " ... ";
      }
      
      
      /*
       * \param SV radius of smearing
       */
      void set_SV( const std::vector<REAL> SV )
      {
        SV_y = SV[1];   //  <-- TODO: Ask Ronnie: SV[0] und SV[1] ???
#ifdef DIMENSION3
        SV_z = SV[2];
#endif
      }





      template<
        typename EG
        , typename LFSV
        , typename SF
        , typename R
        >
      bool evaluate_residual_on_element( const EG& eg
                                         , const LFSV& lfsv
                                         , SF& shapefunctions
                                         , R& residual
                                         , const bool bAdjoint=false
                                         //, const RF distribution_area = 1.0
                                         ) const {
        
        Dune::FieldVector<RF,dim> current_source_point ( pointsource_location_global );
		
        std::vector<Dune::FieldVector<RF,dim>> source_points_global;
        source_points_global.push_back( current_source_point );
		
        double contribution = 1.0;
        std::vector<double> partialcontributions;
        partialcontributions.push_back( contribution );

        double distribution_area = 1.0;
        if( take_cell_average )
          distribution_area = eg.geometry().volume();

#ifdef DIMENSION3
        // In 3D use x and y SV for smearing! (no smoothing needed, as it is only used in the adjoint (overshots doesn't matter, but maybe we include smoothing one day!))
        if( SV_y > 1e-12 && SV_z > 1e-12 ){
          //distribution_area = SV_y*SV_z;-> not needed. Linear system, at the end I have to multiply the Sensitivity for m0 by this factor -> leave it here larger!
          
          int nStretches_y = std::floor( (SV_y/2.0) / (double)inputdata.domain_data.yasp_gridsizes[1] );
          int nStretches_z = std::floor( (SV_z/2.0) / (double)inputdata.domain_data.yasp_gridsizes[2] );

          for( int iStretch_z = 0; iStretch_z < nStretches_z; ++iStretch_z ){
            current_source_point[ 2 ] += inputdata.domain_data.yasp_gridsizes[ 2 ];
            source_points_global.push_back( current_source_point );
            partialcontributions.push_back( contribution );
            
          }
          current_source_point=pointsource_location_global;
          for( int iStretch_z = 0; iStretch_z < nStretches_z; ++iStretch_z ){
            current_source_point[ 2 ] -= inputdata.domain_data.yasp_gridsizes[ 2 ];
            source_points_global.push_back( current_source_point );
            partialcontributions.push_back( contribution );
				  
          }
			
          current_source_point = pointsource_location_global;
          for( int iStretch_y = 0; iStretch_y < nStretches_y; ++iStretch_y )
            {
              current_source_point[ 1 ] += inputdata.domain_data.yasp_gridsizes[ 1 ];
              Dune::FieldVector<RF,dim> current_source_point_tmp ( current_source_point );
				
              source_points_global.push_back( current_source_point );
              partialcontributions.push_back( contribution );
				
              for( int iStretch_z = 0; iStretch_z < nStretches_z; ++iStretch_z ){
                current_source_point[ 2 ] += inputdata.domain_data.yasp_gridsizes[ 2 ];
                source_points_global.push_back( current_source_point );
                partialcontributions.push_back( contribution );
				  
              }
				
              current_source_point=current_source_point_tmp;
              for( int iStretch_z = 0; iStretch_z < nStretches_z; ++iStretch_z ){
                current_source_point[ 2 ] -= inputdata.domain_data.yasp_gridsizes[ 2 ];
                source_points_global.push_back( current_source_point );
                partialcontributions.push_back( contribution );
				  
              }
              current_source_point=current_source_point_tmp;
				
            }

          current_source_point = pointsource_location_global;

          for( int iStretch_y = 0; iStretch_y < nStretches_y; ++iStretch_y )
            {
              current_source_point[ 1 ] -= inputdata.domain_data.yasp_gridsizes[ 1 ];
              Dune::FieldVector<RF,dim> current_source_point_tmp ( current_source_point );
				
              source_points_global.push_back( current_source_point );
              partialcontributions.push_back( contribution );
				
              for( int iStretch_z = 0; iStretch_z < nStretches_z; ++iStretch_z ){
                current_source_point[ 2 ] += inputdata.domain_data.yasp_gridsizes[ 2 ];
                source_points_global.push_back( current_source_point );
                partialcontributions.push_back( contribution );
				  
              }
              current_source_point=current_source_point_tmp;
              for( int iStretch_z = 0; iStretch_z < nStretches_z; ++iStretch_z ){
                current_source_point[ 2 ] -= inputdata.domain_data.yasp_gridsizes[ 2 ];
                source_points_global.push_back( current_source_point );
                partialcontributions.push_back( contribution );
				  
              }
				
              current_source_point=current_source_point_tmp;
            }
        }
#else
        // in 2D only y coordinate for smearing is used.
        if( SV_y > 1e-12 )
          {
            //distribution_area = SV_y; -> not needed. Linear system, at the end I have to multiply the Sensitivity for m0 by this factor -> leave it here larger!
			
            int nStretches = std::floor( (SV_y/2.0) / (double)inputdata.domain_data.yasp_gridsizes[1] );

            for( int iStretch = 0; iStretch < nStretches; ++iStretch )
              {
                current_source_point[ 1 ] += inputdata.domain_data.yasp_gridsizes[ 1 ];
                source_points_global.push_back( current_source_point );
                partialcontributions.push_back( contribution );
              }

            current_source_point = pointsource_location_global;

            for( int iStretch = 0; iStretch < nStretches; ++iStretch )
              {
                current_source_point[ 1 ] -= inputdata.domain_data.yasp_gridsizes[ 1 ];
                source_points_global.push_back( current_source_point );
                partialcontributions.push_back( contribution );
              }
          }
#endif


        // use delta distribution as it is
        for( UINT iSP = 0; iSP < source_points_global.size(); ++iSP ){
          Dune::FieldVector<RF,dim> pointsource_local = eg.geometry().local( source_points_global[ iSP ] );
          // Code by Olaf Cirpka to evaluate shape functions at a single point source!
          // Check if the point source lies inside the rectangular cell.
          UINT iOutside = 0;
          for( UINT i=0; i<dim; i++ ){
            if( (pointsource_local[i] < 0.0 ) || (pointsource_local[i] >= 1.0) ){
              iOutside++;
            }
          }
          if( iOutside==0 ){
            // pointsource is inside the rectangular cell!
            lfsv.finiteElement().localBasis().evaluateFunction( pointsource_local, shapefunctions );
            for( UINT i=0; i<lfsv.size(); i++ ){
              residual.accumulate( lfsv, i, -partialcontributions[ iSP ] / distribution_area * shapefunctions[i] * sign_factor );
            }
          }
        }

        return true;

      } // bool evaluate_Residual()
      
    }; // class PointSource

  }
}

#endif
