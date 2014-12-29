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

    template<typename FD, typename DIR>
    void cross_covariance_JQ( const FD& fielddata, 
                              const Vector<REAL>& Lambdas,    // eigenvalues of the spatial covariance matrix R_YY per zone, IMPORTANT: size and order need to fit to the parallel data distribution of the FFT
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

      const UINT dim = fielddata.dim; // domain dimension


      CovarianceMatrix<FD> C_YY( fielddata, communicator);

      Vector<REAL> SensitivityFieldJ;   // variable storing the k-th sensitivity field == k-th column of the sensitivity matrix J^T

      Vector<UINT> local_offset(dim,0);
      Vector<UINT> local_count(dim,0);
      C_YY.prepare( local_offset, local_count );


      //logger << "DEBUG: cross_covariance_JQ: local_offset " << local_offset;
      //logger << "DEBUG: cross_covariance_JQ: local_count " << local_count;

      // read the sensitivities
      HDF5Tools::h5_pRead( SensitivityFieldJ,
                           sensitivity_filename,
                           "/Sensitivity",
                           local_offset,
                           local_count,
                           communicator
                           );
        
      // Important note: 
      // SensitivityFieldJ cannot be plottet on the grid view
      // because local_count and local_offset are based on the 
      // FFT parallelization on the extended field!
        

      Vector<REAL> JQ_needed( SensitivityFieldJ.size() );
      C_YY.mv( Lambdas, SensitivityFieldJ, JQ_needed );

      //write the data to HDF5

      /*
        Important note: 
        FFT parallelization happens on the extended field:
        The gridview must not be used here. It has a different
        parallel decomposition.
        I.e. h5g_pWrite( gv ,...) must not be used here.
      */
      logger << "cross_covariance_JQ: h5_pWrite for JQ ..." << std::endl;

      HDF5Tools::h5_pWrite( JQ_needed,
                            JQ_filename,
                            "/JQ",
                            fielddata.nCells,
                            local_offset,
                            local_count,
                            communicator
                            );

      logger << "cross_covariance_JQ: h5_pWrite for JQ done." << std::endl;



    }

  } // namespace Gesis

} // namespace Dune

#endif // DUNE_GESIS_CROSS_COVARIANCE_JQ_HH
