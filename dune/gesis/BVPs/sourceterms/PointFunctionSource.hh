#ifndef DUNE_GESIS_POINT_FUNCTION_SOURCE_HH
#define DUNE_GESIS_POINT_FUNCTION_SOURCE_HH

#include "dune/gesis/common/deltafunction.hh"
#include "SourceBaseClass.hh"


namespace Dune {

  namespace Gesis {

    /*
     * The PointFunctionSource uses the Gaussian to approx. the peak.
     */

    template<typename GFS,typename IDT>
    class PointFunctionSource
      : public SourceTermInterface<
      SourceTermTraits<typename GFS::Traits::GridViewType,REAL>,
      PointFunctionSource<GFS,IDT>
      >
    {
    private:

      typedef typename GFS::Traits::GridViewType GV;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type VCType;
      typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> DGF;

      enum{dim=GV::dimension};
      //const GV& gv;
      const GFS& gfs;
      const IDT& inputdata;
      const REAL pointValue;
      const REAL pointSmearing;
      const bool bfixedwidth;
      VCType vcPointSource;
      Dune::FieldVector<REAL,dim> pointsource_location_global;

    public:
      const SourceNature source_nature;
      typedef SourceTermTraits<GV,REAL> Traits;
      
      PointFunctionSource( //const GV& gv_, 
                           const GFS& gfs_, 
                           const IDT& inputdata_,
                           const REAL pointValue_=REAL(1),
                           const REAL pointSmearing_=REAL(1),
                           const bool bfixedwidth_=false
                           )
        :
        gfs(gfs_),
        //gv(gfs.gridView()),
        inputdata(inputdata_),
        pointValue(pointValue_),
        pointSmearing(pointSmearing_),
        bfixedwidth(bfixedwidth_),
        vcPointSource(gfs,0.0),
        source_nature(POINT_SOURCE)
      {}


      void set_SV(const std::vector<REAL>& setSV) {
        pointSmearing = setSV[0];
      }


      template<typename COORDINATES>
      void set_PointSourceLocation( const COORDINATES& ps_global ) {

        if( pointSmearing > PEAK_THRESHOLD ){
          const GV& gv=gfs.gridView();
          Dune::Gesis::DeltaFunction<GV,REAL,COORDINATES,IDT>
            ps_function( gv, inputdata, ps_global, pointValue, pointSmearing, bfixedwidth );
          Dune::PDELab::interpolate( ps_function, gfs, vcPointSource );
        }
        else{
          // Do not use Gaussian, instead set the peak at:
          pointsource_location_global = ps_global;
        }
      }


      
      void plot2vtu( const std::string& filename ){

        VTKPlot::output2vtu( gfs, 
                             vcPointSource, 
                             filename.c_str(), 
                             "pointsource", 
                             pMAX-1 
                             );

      }




      template<typename EG
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

        if( pointSmearing > PEAK_THRESHOLD ){
          
          DGF dgf( gfs, vcPointSource );
          // For the evaluation of the dgf, a Dune::FieldVector of dim 1 is required:
          Dune::FieldVector<REAL,1> fvec(0.0);

          // dimensions
          const int dim = EG::Geometry::dimension;
          const int order = lfsv.finiteElement().localBasis().order();
          const int intorder = 2*order;
          
          // select quadrature rule
          Dune::GeometryType gt = eg.geometry().type();
          const Dune::QuadratureRule<CTYPE,dim>& rule = Dune::QuadratureRules<CTYPE,dim>::rule(gt,intorder);

          double distribution_area = 1.0;
          distribution_area = eg.geometry().volume();
          
          // loop over quadrature points
          for(typename Dune::QuadratureRule<CTYPE,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it) {
            // evaluate shape functions
            typedef Dune::FieldVector<REAL,1> RangeType;
            std::vector<RangeType> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            dgf.evaluate( eg, it->position(), fvec );
            REAL fval = fvec[0] / distribution_area;

            // integrate fval
            REAL factor = it->weight() * eg.geometry().integrationElement(it->position());
            for(UINT i=0; i<lfsv.size(); i++)
              residual.accumulate(lfsv,i,-fval*phi[i]*factor);
          }

        }

        else {

          // Do not use Gaussian, instead evaluate the peak!
          // Use delta distribution as it is.

          typedef REAL RF;

          Dune::FieldVector<RF,dim> current_source_point ( pointsource_location_global );
		
          std::vector<Dune::FieldVector<RF,dim>> source_points_global;
          source_points_global.push_back( current_source_point );
		
          double contribution = 1.0;
          std::vector<double> partialcontributions;
          partialcontributions.push_back( contribution );

          double distribution_area = 1.0; // Best 2D-results so far using this!
          
          //if( take_cell_average )
          //distribution_area = eg.geometry().volume();

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
                residual.accumulate( lfsv, i, -partialcontributions[ iSP ] / distribution_area * shapefunctions[i] * pointValue );
              }
            }
          }

        }
        return true;
      }
      
    }; // class PointFunctionSource

      
  }
}

#endif
