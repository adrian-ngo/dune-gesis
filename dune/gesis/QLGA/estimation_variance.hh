/* 
 * File:   Estimation_Variance.hh
 * Author: ronnie
 */

#ifndef DUNE_GESIS_ESTIMATION_VARIANCE_HH
#define DUNE_GESIS_ESTIMATION_VARIANCE_HH

namespace Dune {
  namespace Gesis {

    template<typename GV,typename IDT,typename DIR>
    void estimation_variance( const GV& gv,
                              const IDT& inputdata,
                              DIR& dir,
                              Vector<REAL>& sigma
                              ){

      int mpi_rank = gv.comm().rank();

      Dune::Gesis::DenseMatrix<REAL> invM(0,0,0.0);

      invM.read_from_HDF5(dir.cokriging_inverse_h5file,"/Minv",inputdata);

      UINT nzones=inputdata.yfield_properties.nz;
      UINT nmeas=invM.n_rows()-nzones;
      dir.set_inversion(nmeas);

      std::vector< Vector<REAL> > JQ(nmeas), X(nzones);
      Vector<UINT> local_count,local_offset;

      for(UINT ii=0; ii<nzones; ii++) {
        HDF5Tools::h5g_pRead<GV,Dune::Interior_Partition>( gv,
                                                           X[ii],
                                                           dir.zonation_matrix[ii],
                                                           "/X",
                                                           local_count,
                                                           local_offset,
                                                           inputdata
                                                           );
      }

      UINT nAllCells1 = X[0].size();
      nAllCells1 = gv.comm().sum( nAllCells1 );

      if( mpi_rank == 0 ){
        std::cout << "DEBUG: Total number of all grid cells read from /X in zone 0 is " << nAllCells1 << std::endl;
        std::cout << "nmeas = " << nmeas << std::endl;
      }

      for(UINT ii=0; ii<nmeas; ii++) {
        HDF5Tools::h5g_pRead<GV,Dune::Interior_Partition>( gv,
                                                           JQ[ii],
                                                           dir.JQ_h5file[ii],
                                                           "/JQ",
                                                           local_offset,
                                                           local_count,
                                                           inputdata
                                                           );
      }

      Dune::Timer watch;

      UINT nAllCells2 = JQ[nmeas-1].size();
      nAllCells2 = gv.comm().sum( nAllCells2 );

      if( mpi_rank == 0 )
        std::cout << "DEBUG: Total number of all grid cells read from last /JQ in zone 0 is " << nAllCells2 << std::endl;

      if( nAllCells2 != nAllCells1 ){
        std::cout << "WARNING: nAllCells2 != nAllCells1" << std::endl;
      }

      //Vector<REAL> sigma( X[0].size(), 0.0 );
      sigma.resize( X[0].size() );
      
      for(UINT ii=0; ii<nzones; ii++)
        for(UINT jj=0; jj<X[ii].size();jj++ )
          sigma[jj]+=X[ii][jj]*inputdata.yfield_properties.zones[ii].variance;

      for(UINT iCell=0; iCell<X[0].size(); iCell++){
        std::vector<REAL> tmp(nmeas+nzones,0.0);
        for(UINT ii=0; ii<nmeas+nzones; ii++){
          for(UINT jj=0; jj<nmeas+nzones; jj++)
            if(jj<nmeas)
              tmp[ii]+=JQ[jj][iCell]*invM(jj,ii);
            else
              tmp[ii]+=X[jj-nmeas][iCell]*invM(jj,ii);
        }
        for(UINT ii=0; ii<nmeas+nzones; ii++){
          if(ii<nmeas)
            sigma[iCell]-=tmp[ii]*JQ[ii][iCell];
          else
            sigma[iCell]-=tmp[ii]*X[ii-nmeas][iCell];
        }

      }// end for iCell

      //std::cout << "V_EST: sigma.size() = " 
      //          << sigma.size()
      //          << std::endl;
      
      General::log_elapsed_time( watch.elapsed(),
                                 gv.comm(),
                                 inputdata.verbosity,
                                 "EVAL",
                                 "estimation_variance: computing sigma2"
                                 );

      HDF5Tools::h5g_pWrite( gv
                             , sigma
                             , dir.estimation_variance_h5file
                             , "/sigma2"
                             , inputdata
                             , inputdata.domain_data.nCells
                             , 1            // blocksize for P0
                             , FEMType::DG  // P0
                             , 0            // grid level 0 for sigma
                             , false        // bWriteOnOverlap
                             );

#ifdef VTK_PLOT_V_EST
      typedef Dune::Gesis::FFTFieldGenerator<IDT,REAL,GV::dimension> YFG;
      YFG yfg_EstVar( inputdata, dir, gv.comm() );
      VTKPlot::output_hdf5data_to_gridview( gv,
                                            inputdata,
                                            dir.estimation_variance_vtu,
                                            dir.estimation_variance_h5file,
                                            "/sigma2",
                                            yfg_EstVar
                                            );
#endif // VTK_PLOT_V_EST

    }


  } // Gesis

} // Dune

#endif // DUNE_GESIS_ESTIMATION_VARIANCE_HH
