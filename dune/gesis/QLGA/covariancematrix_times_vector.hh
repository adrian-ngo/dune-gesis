#ifndef DUNE_GESIS_COVARIANCE_MATRIX_TIMES_VECTOR_HH
#define DUNE_GESIS_COVARIANCE_MATRIX_TIMES_VECTOR_HH

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

    template<typename IDT>
    class CovarianceMatrix {

    private:
      const Vector<UINT>& dims;
      const Vector<UINT>& dimsExt;

      MPI_Comm communicator;
      const IDT& inputdata;
      UINT nx;
      UINT ny;
      UINT nz;
      UINT Nx;
      UINT Ny;
      UINT Nz;
      UINT dim;
      REAL scalingfactor;

      UINT localn0_Sensitivity;

      // define the FFT variables
      fftw_complex *J_extended;   // The sensitivity field will be extended by zero-padding.
      fftw_complex *FFT_J_extended;
      fftw_complex *Lambdas_FFT_J_extended;

      ptrdiff_t localN0;      // local portion of the 1st dimension (Nz) for the current MPI process
      ptrdiff_t local_0_start; // offset for the current MPI process
      ptrdiff_t alloc_local;   // number of elements to allocate (complex numbers for dft/r2c/c2r plans, real numbers for r2r plans)


    public:
      CovarianceMatrix( const Vector<UINT>& dims_,
                        const Vector<UINT>& dimsExt_,
                        MPI_Comm communicator_,
                        const IDT& inputdata_
                        )
        : dims( dims_ ),
          dimsExt( dimsExt_ ),
          communicator( communicator_ ),
          inputdata( inputdata_ )
      {
        dim = dims.size();
        nx = dims[0];
        ny = dims[1];
        Nx = dimsExt[0];
        Ny = dimsExt[1];
#ifdef DIMENSION3
        nz = dims[2];
        Nz = dimsExt[2];
#endif
        localn0_Sensitivity = 0;
        localN0 = 0;
      };

      void prepare( UINT zone_offset,
                    Vector<UINT>& local_offset,
                    Vector<UINT>& local_count,
                    const bool bBigCount = false )
      {

        Dune::Timer watch;

        local_offset.reset( dim );
        local_count.reset( dim );

        UINT VectorSizeExt = dimsExt.componentsproduct();
        scalingfactor = REAL ( VectorSizeExt ); // normalization factor of the discrete Fourier transform

        /*
         * Definition of the block on each processor for the parallel FFT
         *   Idea: Maybe later another data distribution. "fftw_mpi_local_size_many" is a more general routine
         *         BUT: don't forget to load also the Lambdas in the right format
         *
         * 
         *   --> Documentation: http://www.fftw.org/doc/MPI-Data-Distribution-Functions.html#MPI-Data-Distribution-Functions
         */
        
        ptrdiff_t local_n0;
#ifdef DIMENSION3
        alloc_local = fftw_mpi_local_size_3d( Nz, Ny, Nx,
                                              communicator,
                                              &local_n0,
                                              &local_0_start
                                              );
#else
        alloc_local = fftw_mpi_local_size_2d( Ny, Nx,
                                              communicator,
                                              &local_n0,
                                              &local_0_start
                                              );
#endif


        logger << "C_YY.prepare: alloc_local = " << alloc_local << std::endl;
        logger << "C_YY.prepare: local_n0 = " << local_n0 << std::endl;
        logger << "C_YY.prepare: local_0_start = " << local_0_start << std::endl;

        //allocate memory for the forward transform
        J_extended             = (fftw_complex*) fftw_malloc( sizeof(fftw_complex)*alloc_local);
        FFT_J_extended         = (fftw_complex*) fftw_malloc( sizeof(fftw_complex)*alloc_local);

        /*
         * Here is the loading procedure for the sensitivities.
         * FFTW needs the data in blocks which do not necessarily coincide with processor partitions of the DUNE grid!
         */



        // define the parameters to read the right hyperslap/block from the HDF5 file
        for(UINT jj=0;jj<dim;jj++){
          //local_count[jj]     = inputdata.domain_data.nCells[jj];
          if( bBigCount )
            local_count[jj]     = dimsExt[jj];
          else
            local_count[jj]     = dims[jj];
        }

        UINT nDepth = 0;
#ifdef DIMENSION3
        int iDepth = 2;
        if( bBigCount )
          nDepth = Nz;
        else
          nDepth = nz;
#else
        int iDepth = 1;
        if( bBigCount )
          nDepth = Ny;
        else
          nDepth = ny;
#endif

        UINT local_n0_Sensitivity = 0;
        //correction!! because the local_count only in one direction is different to the # of cells
        //  because of 1-D data distribution
        if( (UINT)local_0_start < nDepth ) {

          local_offset[ iDepth ] = local_0_start + zone_offset;

          if( (UINT)local_n0 + (UINT)local_0_start <= nDepth ){
            local_count[ iDepth ] = local_n0;
            local_n0_Sensitivity  = local_n0;
          }else{
            local_count[ iDepth ] = nDepth - local_0_start;
            local_n0_Sensitivity  = nDepth - local_0_start;
          }
        } else {
          // the sensitivities have dimension of the grid and not the extended grid. So on some processors nothing needes to be read!!!
          local_count[ iDepth ] = 0;
        }

        localN0 = local_n0;
        if( !bBigCount ){
          localn0_Sensitivity = local_n0_Sensitivity;
        }


        General::log_elapsed_time( watch.elapsed(),
                                   communicator,
                                   inputdata.verbosity,
                                   "FFTW",
                                   "covariance matrix prepare" );

      }



      void mv( const Vector<REAL>& Lambdas,
               const Vector<REAL>& SensitivityFieldJ,
               Vector<REAL>& JQ_needed ) {

        Dune::Timer watch;

        //build up the extended sensitivities in 3D
#ifdef DIMENSION3

        for( UINT iz=0; iz < (UINT)localN0; iz++ ){
          for( UINT iy=0; iy < Ny; iy++ ){
            for( UINT ix = 0; ix < Nx; ix++ ){
              UINT L = General::indexconversion_3d_to_1d( iz, iy, ix, Nz, Ny, Nx );
              UINT l = General::indexconversion_3d_to_1d( iz, iy, ix, localn0_Sensitivity, ny, nx );
              // only at these locations the sensitivities needs to filled in
              if( (iz+local_0_start < nz) && (iy < ny) && (ix < nx) ){
                J_extended[ L ][0] = SensitivityFieldJ[ l ]; // real part
                J_extended[ L ][1] = 0.0; // imaginary part
              }
              else{
                // zero-padding:
                J_extended[ L ][0] = 0.0; // real part
                J_extended[ L ][1] = 0.0; // imaginary part
              }
            }
          }
        }

#else //build up the extended sensitivities in 2D


        logger << "C_YY.mv: localN0 = " << localN0 << std::endl;
        logger << "C_YY.mv: localn0_Sensitivity = " << localn0_Sensitivity << std::endl;

        for( UINT iy=0; iy < (UINT)localN0; iy++ ){
          for( UINT ix = 0; ix < Nx; ix++ ){
            UINT L = General::indexconversion_2d_to_1d( iy, ix, Ny, Nx );
            UINT l = General::indexconversion_2d_to_1d( iy, ix, localn0_Sensitivity, nx );
            // only at these locations the sensitivities needs to filled in

            //if( (iy+local_0_start < inputdata.domain_data.nCells[1]) && (ix < inputdata.domain_data.nCells[0]) ){
            //logger<<"zone : "<<ii<<", local_0_start : "<<local_0_start<<std::endl;
            if( (iy+local_0_start < ny) && (ix < nx) ){
              //logger<<"ix:"<<ix<<", iy :"<< iy<<", L : "<<L<<", l : "<<l<<std::endl;
              J_extended[ L ][0] = SensitivityFieldJ[ l ];
              J_extended[ L ][1] = 0.0;
            }
            else{
              J_extended[ L ][0] = 0.0;
              J_extended[ L ][1] = 0.0;
            }
            /*
            logger << " ix: " << ix 
                   <<", iy: " << iy
                   <<", L : " << L
                   <<", l : " << l 
                   << " J_extended = " << J_extended[ L ][0]
                   << std::endl;
            */

          }
        }

#endif


        fftw_plan forward_plan;
        logger << "C_YY.mv: Define forward fftw plan." << std::endl;
        // define the forward FFT plan
#ifdef DIMENSION3
        forward_plan = fftw_mpi_plan_dft_3d( Nz, Ny, Nx,
                                             J_extended,
                                             FFT_J_extended,
                                             communicator,
                                             FFTW_FORWARD,
                                             FFTW_ESTIMATE );
#else
        forward_plan = fftw_mpi_plan_dft_2d( Ny, Nx,
                                             J_extended,
                                             FFT_J_extended,
                                             communicator,
                                             FFTW_FORWARD,
                                             FFTW_ESTIMATE );
#endif

        logger << "C_YY.mv: fftw plan defined." << std::endl;

        // do the forward FFT!
        fftw_execute( forward_plan );
        logger << "C_YY.mv: Forward fftw executed." << std::endl;

        fftw_destroy_plan( forward_plan );
        logger << "C_YY.mv: Forward fftw plan destroyed." << std::endl;

        fftw_free( J_extended );
        logger << "C_YY.mv: J_extended freed." << std::endl;


        //allocate memory for the multiplication
        Lambdas_FFT_J_extended = (fftw_complex*) fftw_malloc( sizeof(fftw_complex)*alloc_local);

        // perform the multiplication important L needs to be smaller than Lambdas.size because FFT_J_extended might be larger than Lambdas size!
        for( UINT L=0; L<Lambdas.size(); L++ ){
          Lambdas_FFT_J_extended[ L ][0]
            = Lambdas[ L ] * FFT_J_extended[ L ][0] / scalingfactor;
          Lambdas_FFT_J_extended[ L ][1]
            = Lambdas[ L ] * FFT_J_extended[ L ][1] / scalingfactor;
        }
        fftw_free( FFT_J_extended );

        //allocate memory for the backward FFT
        fftw_complex *JQ;
        JQ = (fftw_complex*) fftw_malloc( sizeof(fftw_complex)*alloc_local);

        //define the plan for the backward FFT
        fftw_plan backward_plan;
#ifdef DIMENSION3
        backward_plan = fftw_mpi_plan_dft_3d( Nz, Ny, Nx,Lambdas_FFT_J_extended, JQ,
                                              communicator,
                                              FFTW_BACKWARD , FFTW_ESTIMATE);
#else
        backward_plan = fftw_mpi_plan_dft_2d( Ny, Nx, Lambdas_FFT_J_extended, JQ,
                                              communicator,
                                              FFTW_BACKWARD , FFTW_ESTIMATE);
#endif

        logger << "C_YY.mv: Backward fftw plan defined." << std::endl;

        //do the backward FFT!!
        fftw_execute( backward_plan );
        logger << "C_YY.mv: Backward fftw executed." << std::endl;

        fftw_destroy_plan( backward_plan );
        logger << "C_YY.mv: Backward fftw plan destroyed." << std::endl;

        fftw_free( Lambdas_FFT_J_extended );
        logger << "C_YY.mv: Lambdas_FFT_J_extended freed." << std::endl;

        /*
         * Extract the needed data out of the result!
         */


#ifdef DIMENSION3 // extraction in 3D
        for( UINT iz=0; iz < (UINT)localN0; iz++ ){
          for( UINT iy=0; iy < Ny; iy++ ){
            for( UINT ix = 0; ix < Nx; ix++ ){
              UINT L = General::indexconversion_3d_to_1d( iz, iy, ix, Nz, Ny, Nx );
              UINT l = General::indexconversion_3d_to_1d( iz, iy, ix, localn0_Sensitivity, ny, nx );
              if( (iz+local_0_start < nz) && (iy < ny) && (ix < nx) ){
                JQ_needed[ l ] = JQ[ L ][0];
              }
            }
          }
        }
#else  // extraction in 2D
        for( UINT iy=0; iy < (UINT)localN0; iy++ ){
          for( UINT ix = 0; ix < Nx; ix++ ){
            UINT L = General::indexconversion_2d_to_1d( iy, ix, Ny, Nx );
            UINT l = General::indexconversion_2d_to_1d( iy, ix, localn0_Sensitivity, nx );
            if( (iy+local_0_start < ny) && (ix < nx) ){
              JQ_needed[ l ] = JQ[ L ][0];
            }
          }
        }
#endif
        fftw_free( JQ );
        logger << "C_YY.mv: JQ freed." << std::endl;

        General::log_elapsed_time( watch.elapsed(),
                                   communicator,
                                   inputdata.verbosity,
                                   "FFTW",
                                   "covariance matrix times vector" );


      }

    }; // class CovarianceMatrix

  } // Gesis

} // Dune

#endif // DUNE_GESIS_COVARIANCE_MATRIX_TIMES_VECTOR_HH

