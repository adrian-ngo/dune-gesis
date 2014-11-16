/* 
 * File:   take_measurements_on_dgf.hh
 * Author: A. Ngo
 * 2010-2014
 *
 * 2012/11/19: evaluate mesuring point only on one process (interior partition)
 * 2014/07/29: add time logger
 */

#ifndef DUNE_GESIS_TAKE_MEASUREMENTS_ON_DGF_HH
#define	DUNE_GESIS_TAKE_MEASUREMENTS_ON_DGF_HH

#include "dune/gesis/common/general.hh"

extern CLogfile logger;

namespace Dune {

  namespace Gesis {

    // Take care:
    // 1.) 
    // This method contains a check which assumes that ALL mesh cells are cuboidal.
    // Points on the lower boundary of the domain will not be detected!
    //
    // 2.) 
    // We take measurements only on interior elements. This is to make sure
    // that each measurement value belongs to a unique process.
    // This will become important later when we use MPI_Allreduce (... MPI_SUM ...)
    // to distribute measuring values to ALL processors.
    //
    template<typename GFS, 
             typename VCType, 
             typename VEC>
    void take_measurements_on_dgf( const GFS& gfs
                                   , const VCType& xSolution
                                   , VEC& measuring_points
                                   , const Dune::MPIHelper& helper) {

      Dune::Timer watch;
  
      typedef typename GFS::Traits::GridViewType GV;
      const GV& gv = gfs.gridView();

      const UINT dim = GV::dimension;
      const UINT nPoints = measuring_points.size();

      typedef typename GV::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator InteriorElementIterator;

      // make discrete function object
      typedef Dune::PDELab::DiscreteGridFunction<GFS, VCType> DGF;
      DGF dgf(gfs, xSolution);

      // Loop over grid cells, check if the current point belongs to that cell.
      // If yes, get the local coordinate and evaluate the dgf on that cell at that point

      // typedef typename GV::Grid GRIDTYPE;
      logger << "Taking measurements... " << std::endl;

      // Loop over list of measuring points
      for (UINT index = 0; index < nPoints; index++) {
    
        //logger << "measuring point index " << index << std::endl;

        Dune::FieldVector<CTYPE, dim> currentpoint(0.0);
        currentpoint[0] = measuring_points[index].x;
        currentpoint[1] = measuring_points[index].y;
#ifdef DIMENSION3
        currentpoint[2] = measuring_points[index].z;
#endif

        //logger << "measuring point global " << currentpoint << std::endl;

        //MPI_Status stat;

        // loop over grid elements, only on interior elements!
        //const typename GV::IndexSet& is = gv.indexSet();
        for (InteriorElementIterator it = gv.template begin<0,Dune::Interior_Partition>()
               ; it != gv.template end<0,Dune::Interior_Partition>()
               ; ++it) {

          //int id = is.index(*it);
          //logger << "grid cell " << id << std::endl;
    

          // Get the local coordinate:
          Dune::FieldVector<CTYPE, dim> currentpoint_local 
            = it->geometry().local( currentpoint );

          //logger << "currentpoint_local " << currentpoint_local << std::endl;

          // Check if the current point belongs to that cell
          // Never ever compare against >=0 or >=-1e-16. Better use >= -1e-12 for double.
          if (currentpoint_local[0] >= -1e-12 && currentpoint_local[0] < 1.0 - 1e-12
              && currentpoint_local[1] >= -1e-12 && currentpoint_local[1] < 1.0 - 1e-12
#ifdef DIMENSION3
              && currentpoint_local[2] >= -1e-12 && currentpoint_local[2] < 1.0 - 1e-12
#endif
              ) 
            {

              //logger << "measuring point local " << currentpoint_local << std::endl;

              Dune::FieldVector<REAL, 1 > measured_value;
              dgf.evaluate( *it, currentpoint_local, measured_value );
              measuring_points[index].value = measured_value[0];

              //logger << "measured value " << measured_value << std::endl;

              break; // Get out of grid loop and to the next measuring point!
            }

        }  // end of loop over leaf elements

        if( General::verbosity >= VERBOSITY_DEBUG_LEVEL ){
          logger << "m" << index
                 << " on process " << helper.rank()
                 << " at " << currentpoint[0] << "," << currentpoint[1] 
#ifdef DIMENSION3
                 << "," << currentpoint[2] 
#endif
                 << "     measured value = " << measuring_points[index].value
                 << std::endl;
        }

      } // end of loop over measurement points

      std::stringstream jobtitle;
      jobtitle << "take_measurements_on_dgf: nPoints = " << nPoints;
      General::log_elapsed_time( watch.elapsed(),
                                 gv.comm(),
                                 General::verbosity,
                                 "EVAL",
                                 jobtitle.str() );
	
    } // end void take_measurements_on_dgf






    // Take care:
    // This method contains a check which assumes that ALL mesh cells are cuboidal.
    // Points on the lower boundary of the domain will not be detected!
    template<typename GV, typename GFS, typename VCType, typename VEC>
    void take_measurements_of_differences_on_dgf( const GV& gv
                                                  , const GFS& gfs
                                                  , const VCType& xSolution
                                                  , VEC& measuring_points
                                                  , const Dune::MPIHelper& helper ) {
      Dune::Timer watch;

      const UINT dim = GV::dimension;
      const UINT nPoints = measuring_points.size();
      
      typedef typename GV::template Codim < 0 > ::Iterator LeafIterator;

      // make discrete function object
      typedef Dune::PDELab::DiscreteGridFunction<GFS, VCType> DGF;
  
      DGF dgf(gfs, xSolution);

      // Loop over grid cells, check if the current point belongs to that cell.
      // If yes, get the local coordinate and evaluate the dgf on that cell at that point

      typedef typename GV::Grid GRIDTYPE;
      // make a mapper for codim 0 entities in the leaf grid
      Dune::LeafMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout> mapper(gv.grid()); // get the underlying hierarchical grid out ot the gridview

      //    logger << "Taking measurements... " << std::endl;
      //    logger << helper.rank() << " mapper-size=" << mapper.size() << std::endl;


      //const typename GV::IndexSet& is = gv.indexSet();


      for (UINT index = 0; index < nPoints; index++){
      
        Dune::FieldVector<REAL, 1 > measured_value1;
        Dune::FieldVector<REAL, 1 > measured_value2;
        measured_value1[0]=-1e300;
        measured_value2[0]=-1e300;
      
        // first iteration!
        for (LeafIterator it = gv.template begin < 0 > ()
               ; it != gv.template end < 0 > ()
               ; ++it) {
          //int id = is.index(*it);
          Dune::FieldVector<CTYPE, dim> currentpoint(0.0);
          currentpoint[0] = measuring_points[index].x;
          currentpoint[1] = measuring_points[index].y;
#ifdef DIMENSION3
          currentpoint[2] = measuring_points[index].z;
#endif
			
          // Get the local coordinate:
          Dune::FieldVector<CTYPE, dim> currentpoint_local = it->geometry().local( currentpoint );
          // Check if the current point belongs to that cell
          if (currentpoint_local[0] >= 0.0 && currentpoint_local[0] < 1.0
              && currentpoint_local[1] >= 0.0 && currentpoint_local[1] < 1.0
#ifdef DIMENSION3
              && currentpoint_local[2] >= 0.0 && currentpoint_local[2] < 1.0
#endif
              ) 
            {
              dgf.evaluate( *it, currentpoint_local, measured_value1 );
            }
        } // end of loop over leaf elements
      
        if( helper.size() > 1 ){
          // send measured value to all processes, important becasue measpoint 1 and 2 might lie on different processors
          MPI_Allreduce( &(measured_value1[0]),
                         &(measured_value1[0]),
                         1,
                         MPI_DOUBLE,
                         MPI_MAX,
                         helper.getCommunicator() );
        }
      
      
        // second iteration!
        for (LeafIterator it = gv.template begin < 0 > ()
               ; it != gv.template end < 0 > ()
               ; ++it) {
          //int id = is.index(*it);
          Dune::FieldVector<CTYPE, dim> currentpoint(0.0);
          currentpoint[0] = measuring_points[index].x2;
          currentpoint[1] = measuring_points[index].y2;
#ifdef DIMENSION3
          currentpoint[2] = measuring_points[index].z2;
#endif
			
          // Get the local coordinate:
          Dune::FieldVector<CTYPE, dim> currentpoint_local = it->geometry().local( currentpoint );
          // Check if the current point belongs to that cell
          if (currentpoint_local[0] >= 0.0 && currentpoint_local[0] < 1.0
              && currentpoint_local[1] >= 0.0 && currentpoint_local[1] < 1.0
#ifdef DIMENSION3
              && currentpoint_local[2] >= 0.0 && currentpoint_local[2] < 1.0
#endif
              ) 
            {
              dgf.evaluate( *it, currentpoint_local, measured_value2 );
	  
	 
	  
            }
        } // end of loop over leaf elements
      
        if( helper.size() > 1 ){
          // send measured value to all processes, important becasue measpoint 1 and 2 might lie on different processors
          MPI_Allreduce( &(measured_value2[0]),
                         &(measured_value2[0]),
                         1,
                         MPI_DOUBLE,
                         MPI_MAX,
                         helper.getCommunicator() );
        }
      
        measuring_points[index].value = (REAL) measured_value1 - (REAL) measured_value2;
      
      }
   
      std::stringstream jobtitle;
      jobtitle << "take_measurements_of_differences_on_dgf: nPoints = " << nPoints;
      General::log_elapsed_time( watch.elapsed(),
                                 gv.comm(),
                                 General::verbosity,
                                 "REDIST",
                                 jobtitle.str() );

    } // end void take_measurements_of_differences_on_dgf

  } // Gesis

} // Dune

#endif	/* DUNE_GESIS_TAKE_MEASUREMENTS_ON_DGF_HH */

