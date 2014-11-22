#ifndef LinePlot_HH
#define LinePlot_HH

/*
#include "lineplot.hh"

...
      std::stringstream gnuplot_datfile;
      gnuplot_datfile << "testplot_" << step << ".dat";

      LinePlot<GV> lplot(gv);
      lplot.setEndPoints( std::vector<REAL>{0,1},
                          std::vector<REAL>{1,0} );
      lplot.addDGData( udgf );
      lplot.write(gnuplot_datfile.str());

*/


template<typename GV>
class LinePlot {

public:
  enum{dim=GV::dimension};
  typedef Dune::FieldVector<REAL,dim> COORD;
  typedef typename GV::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;

  private:
  const GV& gv;
  COORD startpoint;
  COORD endpoint;
  std::stringstream contents;
  std::vector<std::vector<REAL>> datamatrix;

public:
  LinePlot( const GV& gv_ ) :
    gv(gv_)
  {
  }
  
  void setEndPoints( const std::vector<REAL>& startpoint_,
                     const std::vector<REAL>& endpoint_ ) {
    for(int i=0; i<dim; ++i ){
      startpoint[i] = startpoint_[i];
      endpoint[i] = endpoint_[i];
    }

    COORD tmp(startpoint);
    tmp -= endpoint;
    
    contents << "# StartPoint = " << startpoint << std::endl;
    contents << "# EndPoint = " << endpoint << std::endl;
    contents << "# Length = " << tmp.two_norm() << std::endl;
    contents << "#" << std::endl;

  }



  template<typename GFS,typename VCType>
  void addDGData( const GFS& gfs,
                  const VCType& solution,
                  const std::string& function_name ){
    typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> DGF;
    DGF solution_dgf(gfs,solution);
    addDGData( solution_dgf, function_name );
  }


  template<typename DGF>
  void addDGData( const DGF& dgf, const std::string& function_name ){

    COORD tmp(startpoint);
    tmp -= endpoint;
    REAL length = tmp.two_norm();

    //logger << "DEBUG: addDGData: startpoint = " << startpoint << std::endl;
    //logger << "DEBUG: addDGData: endpoint = " << endpoint << std::endl;
    //logger << "DEBUG: addDGData: startpoint - endpoint = " << tmp << std::endl;

    contents << "# " << function_name << std::endl;
    // Iterate over the grid elements and evaluate the gridfunction along the 
    // line between startpoint and endpoint.

    std::vector<REAL> datavector;
    for( ElementIterator eit=gv.template begin<0,Dune::All_Partition>()
           ; eit!=gv.template end<0,Dune::All_Partition>()
           ; ++eit) {

      Dune::FieldVector<REAL,dim> elementcenter
        = eit->geometry().center();

      COORD loc(startpoint);
      loc -= elementcenter;

      //logger << "DEBUG: addDGData: elementcenter = " << elementcenter << std::endl;
      //logger << "DEBUG: addDGData: startpoint - elementcenter = " << loc << std::endl;

      if( std::abs( tmp[0]/loc[0] - tmp[1]/loc[1] ) < 1e-12 
#ifdef DIMENSION3
          &
          std::abs( tmp[0]/loc[0] - tmp[2]/loc[2] ) < 1e-12
#endif
          ) {

        //logger << "DEBUG: addDGData: elementcenter hit = " << elementcenter << std::endl;
        contents << "# DEBUG elementcenter = " << elementcenter << std::endl;

        int nElementCorners = eit->geometry().corners();
        for( int iElementCorner = 0; iElementCorner < nElementCorners; iElementCorner++ ){

          COORD eck(startpoint);

          Dune::FieldVector<REAL,dim> elementcorner 
            = eit->geometry().corner( iElementCorner );

          Dune::FieldVector<REAL,dim> elementcornerLocal 
            = eit->geometry().local( elementcorner );

          eck -= elementcorner;

          bool bAddToPlot = false;
          
          if( std::abs(elementcorner[0]) < 1e-12 
              && std::abs( elementcorner[1] - 1.0 ) < 1e-12 ){
            bAddToPlot = true;
          } 
          else if( std::abs(elementcorner[1]) < 1e-12 
                   && std::abs( elementcorner[0] - 1.0 ) < 1e-12 ) {
            bAddToPlot = true;
          }
          else if ( std::abs( tmp[0]/eck[0] - tmp[1]/eck[1] ) < 1e-12 
#ifdef DIMENSION3
                    &
                    std::abs( tmp[0]/eck[0] - tmp[2]/eck[2] ) < 1e-12 
#endif
                    ) {
            bAddToPlot = true;
          }

          if( bAddToPlot ){
            contents << "# DEBUG elementcorner = " << elementcorner << std::endl;
            Dune::FieldVector<REAL,1> fval;
            dgf.evaluate( *eit, elementcornerLocal, fval );
            contents << eck.two_norm() - 0.5*length << " " 
                     << fval[0]
                     << std::endl;
          }
          
        }

        contents << std::endl;
      }

    } // loop over cells

    contents << std::endl;

  }

  void write( const std::string& filename ){
    std::ofstream filestr( filename.c_str(), std::ios::out );
    if( filestr.is_open() ) {
      filestr << contents.str();
      filestr.close();
    }
  }

};


#endif //LinePlot_HH
