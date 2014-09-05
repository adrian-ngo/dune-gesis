#ifndef CALIBRATE_WELL_CENTERS_HH
#define CALIBRATE_WELL_CENTERS_HH

#include "wellposition.hh"
#include "../common/MPI_Tools.hh"


namespace Dune {
  namespace GeoInversion {

template<typename GV,typename IDT>
void calibrateWellCentersToElementCenters( const GV& gv, 
                                           IDT& inputdata ){
  typedef typename IDT::SDT SDT;
  enum{ dim = GV::dimension };

  UINT nAllWells = 0;
  const UINT nSetups = inputdata.setups.size();
  inputdata.listOfAllWellCenters.resize( nSetups );
  
  //loop over all setups
  for( UINT iSetup=0; iSetup < nSetups; iSetup++ ) {

    SDT& setupdata = inputdata.setups[iSetup];
    UINT nWells = setupdata.wdlist.total;
    inputdata.listOfAllWellCenters[iSetup].resize( nWells );
    nAllWells += nWells;

    typedef typename GV::Grid::template Codim<0>::Entity Entity;
    typedef typename GV::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
    typedef Dune::PDELab::ElementGeometry<Entity> EG;
    //int iElement = 0;

    //loop over all wells in the setup
    for( UINT iWell=0; iWell<nWells; iWell++ ) {

      if( gv.comm().rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION )
        std::cout << "=== Setup " << iSetup << ": Calibrating well # " << iWell << std::endl;

      logger << std::endl;
      logger << "=== Setup " << iSetup << ": Calibrating well # " << iWell << std::endl;

      REAL well_top    =  setupdata.wdlist.pointdata_vector[iWell].well_top;
      REAL well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;
      REAL well_x      =  setupdata.wdlist.pointdata_vector[iWell].x;
      logger << "   well_x = " << well_x << std::endl;
#ifdef DIMENSION3
      REAL well_y      =  setupdata.wdlist.pointdata_vector[iWell].y;
      logger << "   well_y = " << well_y << std::endl;
#endif

      REAL maxcenter = -1e+12;
      REAL mincenter = 1e+12;

      REAL new_center_x = -1e+12;
#ifdef DIMENSION3
      REAL new_center_y = -1e+12;
#endif

      bool bWellDetected = false;

      typedef std::vector< std::pair<REAL,Dune::FieldVector<CTYPE,dim> > > DATAVEC;
      DATAVEC datavector;

      for (ElementIterator it=gv.template begin<0,Dune::All_Partition>();
           it!=gv.template end<0,Dune::All_Partition>();++it) {

        if( w_outside != 
            isElementPenetratedByWell(EG(*it),iWell,inputdata,setupdata) ){

          bWellDetected = true;
          Dune::FieldVector<REAL,dim> cellcenter = it->geometry().center();
          new_center_x = cellcenter[0];
#ifdef DIMENSION3
          new_center_y = cellcenter[1];
#endif
        
          REAL K_well = setupdata.wdlist.pointdata_vector[iWell].well_conductivity;
          
          datavector.push_back( std::make_pair( K_well, cellcenter ) );

          maxcenter = std::max( cellcenter[dim-1], maxcenter );
          mincenter = std::min( cellcenter[dim-1], mincenter );
        
          //std::cout << "DEBUG: cellcenter inside well " << cellcenter << std::endl;

        }

      }
      
      if( inputdata.parallel_fine_tuning.maxNComm > 1 
#ifdef USE_ALUGRID
          ||
          inputdata.domain_data.alugrid_maxsteps > 0
#endif
          ){
        // The grid might get a new partitionining
        // - after adaptive load balancing (ALUGRID 3D)
        // - or during parallel in parallel computation.
        // Then, the well data need to be available on all processors!

        MPI_Tools::redist_well_data( datavector, dim );
        
        if( gv.comm().rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION ){
          std::cout << "Send well data to all processes! Needed later due to possible repartitioning of the grid." << std::endl;
        }

      }

      inputdata.listOfAllWellCenters[iSetup][iWell] = datavector;

      //std::cout << "DEBUG: Send info to all processes! " << std::endl;
      // Send info to all processes!
      if( gv.comm().size() > 1 ){
        maxcenter = gv.comm().max(maxcenter);
        mincenter = gv.comm().min(mincenter);
      }
      //std::cout << "DEBUG: Done for P" << gv.comm().rank() << std::endl;


      // Distribute adjusted values to all processors!
      std::vector<REAL> sendbuffer( 4, -1E+100 );
      std::vector<REAL> recvbuffer( 4, -1E+100 );

      if( bWellDetected ){

        const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1];

        sendbuffer[0] = maxcenter + 0.5 * dz;
        sendbuffer[1] = mincenter - 0.5 * dz;
        sendbuffer[2] = new_center_x;
#ifdef DIMENSION3
        sendbuffer[3] = new_center_y;
#endif

        logger << "Well #" << iWell << std::endl;
        logger << "Original well-top = " 
               << well_top 
               << "    Adjusted well-top = " 
               << sendbuffer[0]
               << std::endl;

        logger << "Original well-bottom = " 
               << well_bottom 
               << "    Adjusted well-bottom = " 
               << sendbuffer[1]
               << std::endl;
      
        logger << "Original well x = " 
               << well_x
               << "    Adjusted well x = " 
               << sendbuffer[2]
               << std::endl;

#ifdef DIMENSION3
        logger << "Original well y = " 
               << well_y
               << "    Adjusted well y = " 
               << sendbuffer[3]
               << std::endl;
#endif
      }


      MPI_Allreduce( &(sendbuffer[0]),
                     &(recvbuffer[0]),
                     4,
                     MPI_DOUBLE,
                     MPI_MAX,
                     MPI_COMM_WORLD );

      setupdata.wdlist.pointdata_vector[iWell].well_top    = recvbuffer[0];
      setupdata.wdlist.pointdata_vector[iWell].well_bottom = recvbuffer[1];
      setupdata.wdlist.pointdata_vector[iWell].x           = recvbuffer[2];
#ifdef DIMENSION3
      setupdata.wdlist.pointdata_vector[iWell].y           = recvbuffer[3];
#endif

    } // end of loop over all wells in the setup

  } // end of loop over all setups

  if( gv.comm().rank()==0 && inputdata.verbosity >= VERBOSITY_EQ_DETAILS ){
    std::cout << "Total number of wells specified in the inputfile = " << nAllWells << std::endl;
  }

  logger << "=== log all wells known to all processors:" << std::endl;
  inputdata.loglistOfAllWellCenters();
}




template<typename GV,typename IDT>
void verifyCalibratedWellsOnRefinedGrid( const GV& gv, 
                                         const IDT& inputdata ){

  Dune::Timer watch;
  logger << "verifyCalibratedWellsOnRefinedGrid:" << std::endl;

  typedef typename IDT::SDT SDT;
  enum{ dim = GV::dimension };

  const UINT nSetups = inputdata.setups.size();
  for( UINT iSetup=0; iSetup < nSetups; iSetup++ ) {

    const SDT& setupdata = inputdata.setups[iSetup];
    UINT nWells = setupdata.wdlist.total;

    typedef typename GV::Grid::template Codim<0>::Entity Entity;
    typedef typename GV::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
    typedef Dune::PDELab::ElementGeometry<Entity> EG;

    //int iElement = 0;

    for( UINT iWell=0; iWell<nWells; iWell++ ) {

      REAL maxcenter = -1e+12;
      REAL mincenter = 1e+12;
    
      bool bWellDetected = false;

      for (ElementIterator it=gv.template begin<0,Dune::All_Partition>();
           it!=gv.template end<0,Dune::All_Partition>();++it) {

        if( w_outside != 
            isElementWithinWellZone(EG(*it),iWell,inputdata,setupdata) ){

          bWellDetected = true;
          Dune::FieldVector<REAL,dim> cellcenter = it->geometry().center();
          maxcenter = std::max(maxcenter,cellcenter[dim-1]);
          mincenter = std::min(mincenter,cellcenter[dim-1]);
        }
      
      }

      if( bWellDetected && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){

        std::cout << "=== A well was calibrated on process " 
                  << gv.comm().rank() << std::endl;
        std::cout << "Well #" << iWell << std::endl;
        std::cout << "Please verify: highest cellcenter inside adjusted well? " << maxcenter << std::endl;
        std::cout << "Please verify: lowest cellcenter inside adjusted well? " << mincenter << std::endl;
      }
    
    }

  }

  General::log_elapsed_time( watch.elapsed(),
                             gv.comm(),
                             inputdata.verbosity,
                             "EVAL",
                             "verifyCalibratedWellsOnRefinedGrid" );

}

  }
}

#endif // CALIBRATE_WELL_CENTERS_HH
