/*
 * This function is to calculate the sensitivity of a lnK measurement wrt to lnK itself.
 *   It calculates also the sensitivity times the covariance matrix (for the current measurement)
 *   and the output.
 * 
 *   "Sensitivity" and "sensitivity times the covariance matrix" are stored on the disk via HDF5
 * 
 * 
 * INPUTS:
 * gv:            grid view
 * inputdata:     the class holding all input information
 * dir:           information about IO Paths and file locations!

 * nCellsExt:     number of all cells of the extended domain_data
 * iteration_number: the current iteration number (only for some VTK output needed)
 * orig_measurements: the measuering points
 * Lambdas:       eigenvalues of the extended covariance matrix (in FFT data distribution format)
 * Y_old:         previous Y field (in format of the gv)
 * helper:        Dune MPI helper (MPI_COMM_WORLD settings)
 * CommunicatorPool : vector containing all information about the different communicators! 
 * 
 * 
 * OUTPUTS
 * JX:            value for the multiplication of the sensitivity times the trend matrix (at the current meas location)
 * J_times_Y_old: value for the multiplication of the sensitivity times old Y field (at the current meas location)
 */

#ifndef DUNE_GESIS_LNK_SENSITIVITIES_HH
#define DUNE_GESIS_LNK_SENSITIVITIES_HH

#include "dune/gesis/common/MyMPIComm.hh"
#include "dune/gesis/common/io/IO_routines.hh"
#include "dune/gesis/QLGA/cross_covariance_JQ.hh"

// Lamndas: needs to fit to the size of the FFT data distribution!!!
// Xzones: zonation matrix -> only on process 0 (process leader of the communicator group) 


namespace Dune{
  namespace Gesis{
    
    template<typename GV,
             typename DIR,
             typename MEASLIST,
             typename IDT
             >
    void lnK_sensitivities( const GV& gv,
                            const IDT& inputdata,
                            const DIR& dir,
                            UINT iSetup,
                            UINT iteration_number,
                            MEASLIST& orig_measurements,
                            const std::vector< Vector<UINT> >& nCellsExt,
                            const std::vector< Vector<REAL> >& Lambdas,
                            const std::vector< Vector<REAL> >& Xzones,
                            Vector<REAL>& Y_old,
                            std::vector<MyMPIComm> CommunicatorPool,
                            //OUTPUT
                            std::vector< Vector<REAL> >& JX,
                            Vector<REAL>& J_times_Y_old
                            )
    { 
      //std::cout << "inside: " << nCellsExt[0][1] << std::endl;
      logger << "lnK_sensitivities(...)" << std::endl;

      // # of lnk measurements
      UINT nPoints_lnk=orig_measurements.nMeasPerSetupPerType(iSetup,0);
  
      // number of zones
      UINT nzones=inputdata.yfield_properties.nz;
  
      //number of available communicators. BE CAREFUL "nComm <= inputdata.maxNComm", because inputdata holds the given input for the pool size and this is a max. value!!!! 
      UINT nComm=CommunicatorPool.size();
  
  
      //some needed variables
      const UINT dim= GV::dimension;
      //char buffer[128]= "";

      //Leaf iterator
      typedef typename GV::template Codim < 0 > ::Iterator LeafIterator;
    
      const typename GV::IndexSet& is = gv.indexSet();
      Dune::FieldVector<double, dim> mp(0.0);
  
      //loop over all lnK measurements
      for(UINT meas_point=0; meas_point<nPoints_lnk; meas_point++ ){
    
        /*
         * distribute the measurements to the available MPI pools!!!! 
         */
        if(CommunicatorPool[meas_point%nComm].I_am_in()){
    
          UINT global_meas_id= orig_measurements.get_global(iSetup,0,meas_point); 
        
          logger<<"working on meas #"<<global_meas_id<<" (parallel part)"<<std::endl;
          
          Vector<REAL> Sensitivity(gv.size(0),0.0);
   
          mp[0]=orig_measurements[global_meas_id].x;
          mp[1]=orig_measurements[global_meas_id].y;
#ifdef DIMENSION3
          mp[2]=orig_measurements[global_meas_id].z;
#endif
    
          //loop over all elements
          for (LeafIterator it = gv.template begin < 0 > ()
                 ; it != gv.template end < 0 > ()
                 ; ++it) {
            int id = is.index(*it);
  
            Dune::FieldVector<double, dim> currentpoint_local = it->geometry().local(mp);
  
            //current measurement point found?
            if (currentpoint_local[0] > 0.0 && currentpoint_local[0] <= 1.0
                && currentpoint_local[1] > 0.0 && currentpoint_local[1] <= 1.0
#ifdef DIMENSION3
                && currentpoint_local[2] > 0.0 && currentpoint_local[2] <= 1.0
#endif
                ) {
              // set sensitivity to 1.0 at this location
              Sensitivity[id]=1.0;
            }
          } // end of loop over leaf elements
  
          // write the sensitivity to the hard disc
          HDF5Tools::
            write_parallel_to_HDF5(
                                   gv
                                   , inputdata
                                   , Sensitivity
                                   , "/Sensitivity"
                                   , inputdata.domain_data.nCells
                                   , dir.Sensitivity_h5file[global_meas_id]
                                   );



       
          //wait in parallel case (means if the used communicator has more than one processor) until all have finished writing -> maybe not needed because the HDF5 might also block!
          if(CommunicatorPool[meas_point%nComm].get_size() >1)
            MPI_Barrier(CommunicatorPool[meas_point%nComm].get_comm());
    
          //calculate cross_covariance_Parallel! (incl. timing for the parallel FFT)
          Dune::Timer watch;
          watch.reset();
          logger<<"calculating JQ PARALLEL ...."<<std::endl;
    
          //calculate JQ!
          cross_covariance_JQ( inputdata, 
                               nCellsExt,
                               Lambdas,
                               dir.Sensitivity_h5file[global_meas_id],
                               dir.JQ_h5file[global_meas_id],
                               CommunicatorPool[meas_point%nComm].get_comm(),
                               global_meas_id,
                               iteration_number,
                               dir
                               );
    
          /*
           * DOES NOT WORK -> I CAN NOT EXCLUDE THE OVERLAPING ELEMENTS, GIVES THEN WRONG RESULTS!!! 
           *                  -> Thus this needs to be sequential in the MPI pool which is responsible for this measurement
           * logger<<"..JQ calculation PARALLEL DONE!!! (took: "<<watch.elapsed()<<" sec)"<<std::endl;
           * for(UINT iCell=0; iCell<Sensitivity.size(); iCell++)
           *     {
           *	 JX[meas_point] += Sensitivity[iCell] * 1.0; //ZonationMatrixX[ iCell ];
           *     J_times_Y_old[meas_point] += Sensitivity[ iCell ] * Y_old[ iCell ];
           *	 }
           *
           */
          logger<<"working on meas #"<<global_meas_id<<" (parallel part Done)"<<std::endl;

          //sequential part: Only P0 (in the responsible MPI pool) works now
          if(CommunicatorPool[meas_point%nComm].get_rank()==0){ // something to do for P0 (in the responsible)
            logger<<"working on meas #"<<global_meas_id<<" (sequential part)"<<std::endl;
     
            // read the sensitivity IF communicator is parallel -> to obtain all data on one processor
      
            if( CommunicatorPool[meas_point%nComm].get_size()>1){
	
              // needed for the HDF5 reading
              Vector<UINT> read_local_count,read_local_offset;
              HDF5Tools::read_sequential_from_HDF5(
                                                   Sensitivity
                                                   , "/Sensitivity"
                                                   , read_local_count
                                                   , read_local_offset
                                                   , dir.Sensitivity_h5file[global_meas_id]
                                                   , inputdata
                                                   );
            }
       
            /*
             * Calculate the other return values (needed for the cokriging matrix)
             */

            const UINT nAllCells = inputdata.domain_data.nCells.componentsproduct();

            for(UINT iCell=0; iCell<nAllCells; iCell++){
              for(UINT ii=0; ii<nzones; ii++)
                JX[global_meas_id][ii] += Sensitivity[iCell] * Xzones[ii][iCell];//X[iCell][ii];//1.0; //ZonationMatrixX[ iCell ];
              J_times_Y_old[global_meas_id] += Sensitivity[ iCell ] * Y_old[ iCell ];
            }
         
            logger<<"working on meas #"<<global_meas_id<<" (sequential part done!)"<<std::endl;
         
	
          } // END: if(CommunicatorPool[meas_point%nComm].get_rank()==0)
      
        }  // END: if(CommunicatorPool[meas_point%nComm].I_am_in())
    
      } // END: for(UINT meas_point=0; meas_point<nPoints_lnk; meas_point++ )
 
    }


  } // Gesis
} // Dune

#endif
