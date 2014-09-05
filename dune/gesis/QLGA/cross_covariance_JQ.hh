#ifndef DUNE_GESIS_CROSS_COVARIANCE_JQ_HH
#define DUNE_GESIS_CROSS_COVARIANCE_JQ_HH

#include "covariancematrix_times_vector.hh"
/*
 *
 * Remark:
 * Since Q is symmetric, computing the first row of JQ is the same as computing the
 * first column of QJ' (with J' being the transposed of J)
 *
 * Background:
 * A circulant (sxs)-matrix Q has the nice property that it can be written as Q = 1/s * F' * Lambda * F,
 * where F is the matrix of discrete fourier transform coefficients
 * and Lambda is the diagonal matrix of eigenvalues of Q.
 *
 * Therefore, multiplying Q with a vetor u simplifies to
 *
 * Q u = 1/s (F' Lambda F) u = 1/s FFT_BACKWARD( Lambda * FFT_FORWARD(u) )
 *
 */

namespace Dune{
  namespace Gesis{

    template<typename IDT, typename DIR>
    void cross_covariance_JQ( const IDT& inputdata,                           // inputdata
                              const std::vector< Vector<UINT> >& nCellsExt,  // size of the extended domain per zone
                              const std::vector< Vector<REAL> >& Lambdas,    // eigenvalues of the spatial covariance matrix R_YY per zone, IMPORTANT: size and order need to fit to the parallel data distribution of the FFT
                              const std::string& sensitivity_filename,       // name where the needed sensitivity is stored
                              const std::string& JQ_filename,                // name where the calculated product is stored
                              const MPI_Comm& communicator, // MPI communicator on which the calculation will be performed
                              UINT global_meas_id,
                              UINT iteration_number,
                              const DIR& dir
                              ){

      // If the corresponding JQ.h5 file exists there is 
      // no need to waste time doing the calculation once again!
      if( General::bFileExists( JQ_filename ) ){
        return;
      }



      const UINT dim = inputdata.dim; // domain dimension
      const UINT nzones = nCellsExt.size(); // total number of zones

      std::vector<UINT>nCells_zones(nzones,0); // zone thickness: How many cells is the zone thick?

      // This is the same code as in the constructor of "FFTFieldGenerator.hh":
      for (UINT ii =0; ii<nzones;ii++){
        if(nzones>1){
	  if(ii==0)
	    nCells_zones[ii] = floor( inputdata.yfield_properties.zones[ii].bottom / inputdata.domain_data.virtual_gridsizes[dim-1] );
	  else if (ii==nzones-1)
	    nCells_zones[ii] = inputdata.domain_data.nCells[ dim-1 ] 
              - floor( inputdata.yfield_properties.zones[ii-1].bottom / inputdata.domain_data.virtual_gridsizes[dim-1] );
	  else
	    nCells_zones[ii] = floor( inputdata.yfield_properties.zones[ii].bottom / inputdata.domain_data.virtual_gridsizes[dim-1] ) 
              - floor( inputdata.yfield_properties.zones[ii-1].bottom / inputdata.domain_data.virtual_gridsizes[dim-1] );
	} else
          nCells_zones[ii] = inputdata.domain_data.nCells[ dim-1 ];
      }

      UINT zone_offset = 0; // This number gets incremented by the zone thickness inside the coming loop.
      // loop over zones:
      for( UINT ii=0; ii<nzones; ii++ ){


        Vector<UINT> dims(dim);     // zone in embedded domain
        Vector<UINT> dimsExt(dim);  // zone in embedding domain
        dims[0] = inputdata.domain_data.nCells[0];
        dimsExt[0] = nCellsExt[ii][0];
        dimsExt[1] = nCellsExt[ii][1];
#ifdef DIMENSION3
        dims[1] = inputdata.domain_data.nCells[1];
        dims[2] = nCells_zones[ii];
        dimsExt[2] = nCellsExt[ii][2];
#else
        dims[1] = nCells_zones[ii];
#endif

        CovarianceMatrix<IDT> C_YY( dims, dimsExt, communicator, inputdata );

        Vector<REAL> SensitivityFieldJ;   // variable storing the k-th sensitivity field == k-th column of the sensitivity matrix J^T

        Vector<UINT> local_offset(dim,0);
        Vector<UINT> local_count(dim,0);
        C_YY.prepare( zone_offset, local_offset, local_count );


        //logger << "DEBUG: cross_covariance_JQ: local_offset " << local_offset;
        //logger << "DEBUG: cross_covariance_JQ: local_count " << local_count;

        // read the sensitivities
        HDF5Tools::
          read_parallel_from_HDF5_without_DUNE( inputdata, 
                                                SensitivityFieldJ,
                                                local_count,
                                                local_offset,
                                                communicator,
                                                "/Sensitivity",
                                                sensitivity_filename );
        
        // Important note: 
        // SensitivityFieldJ cannot be plottet on the grid view
        // because local_count and local_offset are based on the 
        // FFT parallelization on the extended field!
        

        Vector<REAL> JQ_needed( SensitivityFieldJ.size() );
        C_YY.mv( Lambdas[ii], SensitivityFieldJ, JQ_needed );

        //write the data to HDF5
        if(ii>0){

          Vector<UINT> local_count;   // Adrian: Ist das nicht zuviel?
          Vector<UINT> local_offset;

          logger << "cross_covariance_JQ: Append JQ of new zone." << std::endl;
          HDF5Tools::
            write_parallel_to_existing_HDF5_without_DUNE( 
                                                         JQ_needed,
                                                         local_count,
                                                         local_offset,
                                                         communicator,
                                                         "/JQ",
                                                         JQ_filename );
        }
        else {

          /*
            Important note: 
            FFT parallelization happens on the extended field:
            The gridview must not be used here. It has a different
            parallel decomposition.
            I.e. write_parallel_to_HDF5( gv ,...) must not be used here.
           */
          logger << "cross_covariance_JQ: write_parallel_to_HDF5_without_DUNE for JQ ..." << std::endl;

          HDF5Tools::
            write_parallel_to_HDF5_without_DUNE( inputdata,
                                                 inputdata.domain_data.nCells,
                                                 JQ_needed,
                                                 local_count,
                                                 local_offset,
                                                 communicator,
                                                 "/JQ",
                                                 JQ_filename );

          logger << "cross_covariance_JQ: write_parallel_to_HDF5_without_DUNE for JQ done." << std::endl;

        }

        //offset correction for the next zone!
        zone_offset += nCells_zones[ii];

      }
    }

  } // namespace Gesis
} // namespace Dune

#endif // DUNE_GESIS_CROSS_COVARIANCE_JQ_HH
