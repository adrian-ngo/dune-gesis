/* 
 * File:   FFTFieldGenerator.hh
 * Author: A. Ngo (10/2014)
 *
 * Implementation of the FFT-based random field generator by Dietrich and Newsam
 * 
 * Article: "Fast and exact simulation of stationary Gaussian processes through circulant embedding of the covariance matrix" by Dietrich, C. R. and Newsam, G. N. (1997), SIAM J. Sci. Comp. vol 18(4), http://dx.doi.org/10.1137/S1064827592240555
 *
 *
 */

#ifndef FFT_FIELD_GENERATOR_HH
#define	FFT_FIELD_GENERATOR_HH

#include <random>      // C++11 random number generator
#include <fftw3.h>     // Fastest Fourier Transform in the West (www.fftw.org)
#include <fftw3-mpi.h>

#include <time.h>                      // define time()


#include "FieldData.hh"
#include "dune/gesis/common/io/HDF5Tools.hh"


namespace Dune {
  namespace Gesis {


    template<typename FD,typename DIR,int dim>
    class FFTFieldGenerator {

    private:

      const FD &fielddata;
      const DIR &dir;
  
      // MPI communicator and information about rank and size
      // const Dune::MPIHelper &helper;

      const MPI_Comm &comm;
      int my_rank,comm_size;

      // Note:
      // Always remeber that Y=log(K) where log is the natural logarithm.
      // We store only the Y-field inside the HDF5 files.
      // We compute K=std::exp(Y) each time we read from such a HDF5 file.
      // This is done once in
      // parallel_import_from_local_vector(...)
      //


      // The following variables are the data containers.

      Vector<REAL> YFieldVector;   // used for the sequential case only
      Vector<REAL> KFieldVector;   // used for the sequential case only
      Vector<REAL> YFieldVector_Well;
      Vector<REAL> KFieldVector_Well;
    
      Vector<REAL>  LocalKFieldVector;  // used for the parallel case only
      Vector<REAL>  LocalYFieldVector;  // used for the parallel case only
      Vector<REAL>  LocalKFieldVector_Well;
      Vector<REAL>  LocalYFieldVector_Well;

      Vector<UINT>  LocalCount;  // used for the parallel case only
      Vector<UINT>  LocalOffset; // used for the parallel case only
  
    public:
  
      UINT size(){
        return KFieldVector.size();
      };
        
      UINT localSize(){
        return LocalKFieldVector.size();
      };


      // The constructor:
      FFTFieldGenerator( const FD& fielddata_,
                         const DIR& dir_,
                         const MPI_Comm &comm_ 
                         ) 
        : fielddata(fielddata_),
          dir(dir_),
          comm(comm_) {
        //setting my_rank  (rank within the stored communicator)
        MPI_Comm_rank(comm,&my_rank);
	
        //setting comm_size (size of the stored communicator)
        MPI_Comm_size(comm,&comm_size);

        //initialization for the multi-threaded fft
        //fftw_init_threads(); 
      };



      // This function can be used for dim=3 only!
      void generate_eigenvalues( ) {
    
        Dune::Timer watch;
        // i=2: number of cells in z-direction = #levels
        // i=1: number of cells in y-direction = #rows
        // i=0: number of cells in x-direction = #columns

        Vector<UINT> nExt( fielddata.nCellsExt );

        // indices:
        // 0 <= iz <= (Nz-1)
        // 0 <= iy <= (Ny-1)
        // 0 <= ix <= (Nx-1)

        // Always remember this and you will be on the safe side of life!
        // ix (most inner loop) runs faster than iy (inner loop) which runs faster than iz (most outer loop)
        // iz belongs to Nz
        // iy belongs to Ny
        // ix belongs to Nx

        // Nz - Ny - Nx  or ( iz, iy, ix) is the index order you must always keep in mind

        Vector<REAL> delta(dim,0.0);
        for(int i=0;i<dim;++i)
          delta[i] = fielddata.gridsizes[i];
    
        //bool bWisdomFileExists = false;

        /*
         * Definition of the block on each processor for the parallel FFT
         *   Idea: Maybe later another data distribution. "fftw_mpi_local_size_many" is a more general routine
         *         BUT: don't forget to load also the Lambdas in the right format
         */
        ptrdiff_t alloc_local, local_n0, local_0_start;
  
        // get FFTW data
        if(dim==3) {
          alloc_local = fftw_mpi_local_size_3d( nExt[2],
                                                nExt[1],
                                                nExt[0],
                                                comm,
                                                &local_n0,
                                                &local_0_start );
        } else {
          alloc_local = fftw_mpi_local_size_2d( nExt[1],
                                                nExt[0],
                                                comm,
                                                &local_n0, 
                                                &local_0_start);
        }
        
        std::vector<REAL> R_YY(alloc_local);

        Vector<REAL> X( dim, 0.0 );
        Vector<REAL> Lambda( fielddata.correlations );
        UINT l;

        // build up R_YY
        // writing of the R_YY in HDF5
        Vector<UINT> local_count(dim,0);
        Vector<UINT> local_offset(dim,0);

        if( dim==3 ){

          for (UINT iz = 0; iz < (UINT)local_n0; iz++) {
            for (UINT iy = 0; iy < nExt[1]; iy++) {
              for (UINT ix = 0; ix < nExt[0]; ix++) {

                X[2] = std::min( REAL(iz+local_0_start), REAL(nExt[2]-(iz+local_0_start)) ) * delta[2];
                X[1] = std::min( REAL(iy), REAL(nExt[1]-iy) ) * delta[1];
                X[0] = std::min( REAL(ix), REAL(nExt[0]-ix) ) * delta[0];
                l = General::indexconversion_3d_to_1d( iz, iy, ix, UINT(local_n0), nExt[1], nExt[0] );
                R_YY[ l ] = autocovariancefunction( fielddata.fieldVariance
                                                    , X
                                                    , Lambda
                                                    , fielddata.variogramModel
                                                    );
              }
            }
          }
          local_count[0] = fielddata.nCellsExt[0];
          local_count[1] = fielddata.nCellsExt[1];
          local_count[2] = UINT(local_n0);
          local_offset[2] = UINT(local_0_start);

        } else {

          for (UINT iy = 0; iy < (UINT)local_n0; iy++) {
            for (UINT ix = 0; ix < nExt[0]; ix++) {

              X[1] = std::min( REAL(iy+local_0_start), REAL(nExt[1]-(iy+local_0_start)) ) * delta[1];
              X[0] = std::min( REAL(ix), REAL(nExt[0]-ix) ) * delta[0];
              l = General::indexconversion_2d_to_1d( iy, ix, UINT(local_n0), nExt[0] );

              if(l>R_YY.size()){
                std::cout << "ERROR: P" << my_rank << ": l=" << l 
                          << " iy = " << iy
                          << " ix = " << ix
                          << " local_n0 = " << UINT(local_n0)
                          << " nExt[0] = " << nExt[0]
                          << std::endl;
              }

              R_YY[ l ] = autocovariancefunction( fielddata.fieldVariance
                                                  , X
                                                  , Lambda
                                                  , fielddata.variogramModel
                                                  );
            }
          }

          local_count[0] = fielddata.nCellsExt[0];
          local_count[1] = UINT(local_n0);
          local_offset[1] = UINT(local_0_start);

        }


        HDF5Tools::h5_pWrite( R_YY
                              , dir.R_YY_h5file[0]
                              , "/R_YY"
                              , fielddata.nCellsExt
                              , local_offset
                              , local_count
                              ,  comm
                              );
    

        fftw_complex *Eigenvalues;
        Eigenvalues = (fftw_complex*) fftw_malloc( (alloc_local ) * sizeof (fftw_complex) );

        if(dim==3) {
          for( UINT iz = 0; iz < (UINT)local_n0; iz++ ) {
            for( UINT iy = 0; iy < nExt[1]; iy++ ) {
              for( UINT ix = 0; ix < nExt[0]; ix++ ) {
                UINT l = General::indexconversion_3d_to_1d( iz, iy, ix, UINT(local_n0), nExt[1], nExt[0] );
                Eigenvalues[ l ][0] = R_YY[ l ];                  // real part
                Eigenvalues[ l ][1] = 0.0;                        // imaginary part
              }
            }
          }
        }
        else {
          for( UINT iy = 0; iy < (UINT)local_n0; iy++ ) {
            for( UINT ix = 0; ix < nExt[0]; ix++ ) {
              UINT l = General::indexconversion_2d_to_1d( iy, ix, UINT(local_n0), nExt[0] );
              Eigenvalues[ l ][0] = R_YY[ l ];                  // real part
              Eigenvalues[ l ][1] = 0.0;                        // imaginary part
            }
          }
        }

        General::log_elapsed_time( watch.elapsed(),
                                   comm,
                                   General::verbosity,
                                   "EVAL",
                                   "generate_eigenvalues()" );
        watch.reset();

        //fftw_plan_with_nthreads(FFT_THREADS); // tell the FFT routine how many threads should be used
    
        fftw_plan plan_forward;
        //logger << "FFTW_FORWARD(3d): Create 3d forward plan using FFTW_ESTIMATE flag..." << std::endl;
        if(dim ==  3){
          plan_forward = fftw_mpi_plan_dft_3d( nExt[2],
                                               nExt[1],
                                               nExt[0],
                                               Eigenvalues,
                                               Eigenvalues,
                                               comm,
                                               FFTW_FORWARD,
                                               FFTW_ESTIMATE );
        } else {
          plan_forward = fftw_mpi_plan_dft_2d( nExt[1],
                                               nExt[0],
                                               Eigenvalues,
                                               Eigenvalues,
                                               comm,
                                               FFTW_FORWARD,
                                               FFTW_ESTIMATE );
        }
        //logger << "FFTW_FORWARD(3d): executing plan..." << std::endl;

        fftw_execute(plan_forward);
        fftw_destroy_plan(plan_forward);

        General::log_elapsed_time( watch.elapsed(),
                                   comm,
                                   General::verbosity,
                                   "FFTW",
                                   "generate_eigenvalues()" );
        watch.reset();

        //save the FFT of the R_YY
        Vector<REAL> tmp_vec;
	
        // Count positive eigenvalues.
        int nPositive=0; // Eigenvalue > 1E-6
        int nNegative=0; // Eigenvalue < -1E-6
        int nNearZero=0; // -1E-6 <= Eigenvalue <= 1E6

        //extract the FFT result in the vector which should be saved!
        tmp_vec.resize(alloc_local);
        for(UINT i=0; i<(UINT)alloc_local;i++){
          tmp_vec[i] = (double)Eigenvalues[i][0]; // real part
          if( tmp_vec[i] > 1E-6 )
            nPositive++;
          else if( tmp_vec[i] < -1E-6 )
            nNegative++;
          else
            nNearZero++;
        }

        if( fielddata.showEV ){
          std::cout << "Properties of Syy: " << std::endl;
          std::cout << alloc_local << " eigenvalues" << std::endl;          
          std::cout << nPositive << " eigenvalues > 1E-6" << std::endl;
          std::cout << nNegative << " eigenvalues < -1E-6" << std::endl;
          std::cout << nNearZero << " eigenvalues near zero" << std::endl;
        }

        
        General::log_elapsed_time( watch.elapsed(),
                                   comm,
                                   General::verbosity,
                                   "EVAL",
                                   "generate_eigenvalues(): count negative Eigenvalues" );

        // save FFT of R_YY!
        HDF5Tools::h5_pWrite( tmp_vec
                              , dir.EV_h5file[0]
                              , "/FFT_R_YY"
                              , fielddata.nCellsExt
                              , local_offset
                              , local_count
                              , comm
                              );

        fftw_free(Eigenvalues);
	
	
        if(my_rank==0)
          save_yfield_properties( dir.yfield_properties_file );
	
      }; // End of generate_eigenvalues()

    

      void generateY( const bool bNewEigenvalues,
                      const bool bNewField,
                      const bool CR=false )
      {
        Dune::Timer watch;
        
        if( bNewEigenvalues 
            ||
            !load_yfield_properties( dir.yfield_properties_file, my_rank )
            ) {
          std::cout << "New geostatistical data: compute eigenvalues of the circulant matrix." << std::endl;
          generate_eigenvalues();        
        }
        else
          std::cout << "Given geostatistical data match existing field." << std::endl;


        if( my_rank == 0 ){
          if( !bNewField ){
            std::cout << "New field will be generated due to changed geostatistical data." << std::endl;
          } else {
            std::cout << "New field will be generated (user's choice)." << std::endl;
          }
        }
        
        Vector<UINT> nExt( fielddata.nCellsExt );
          
        REAL scalingfactor = REAL( nExt.componentsproduct() );
	
        /*
         * Definition of the block on each processor for the parallel FFT
         *   Idea: Maybe later another data distribution. "fftw_mpi_local_size_many" is a more general routine
         *         BUT: don't forget to load also the Lambdas in the right format
         */
        ptrdiff_t alloc_local, local_n0, local_0_start;
  
        // get the FFT data

        Vector<UINT> local_count(dim,0);
        Vector<UINT> local_offset(dim,0);
        if(dim==3){
          alloc_local = fftw_mpi_local_size_3d( nExt[2],
                                                nExt[1],
                                                nExt[0],
                                                comm,
                                                &local_n0,
                                                &local_0_start );
          local_count[2] = UINT(local_n0);
          local_count[1] = fielddata.nCellsExt[1];
          local_count[0] = fielddata.nCellsExt[0];
          local_offset[2] = UINT(local_0_start);
        }
        else {
          alloc_local = fftw_mpi_local_size_2d( nExt[1],
                                                nExt[0],
                                                comm,
                                                &local_n0,
                                                &local_0_start );
          local_count[1] = UINT(local_n0);
          local_count[0] = fielddata.nCellsExt[0];
          local_offset[1] = UINT(local_0_start);
        }


        //loading in HDF5
        Vector<REAL> tmp_vec; // temporary vector for the loaded eigenvalues
        HDF5Tools::h5_pRead( tmp_vec
                             , dir.EV_h5file[0]
                             , "/FFT_R_YY"
                             , local_offset
                             , local_count
                             , comm
                             );

        watch.reset();

        // initialize pseudo-random generator
        std::random_device rd;
        std::mt19937_64 gen; // 64-bit Mersenne Twister
        std::normal_distribution<> ndist(0,1); // mean = 0, sigma = 1

        int seed = fielddata.seed;
        if( seed == 0 ){
          seed = rd();
        }
        // different seed for different processors -> very IMPORTANT to obtain the right result!
        seed += my_rank;

        gen.seed( seed );
            
        if(my_rank==0)
          std::cout << "seed = " << seed << std::endl;

	
        fftw_complex *KField;
        KField = (fftw_complex*) fftw_malloc( ( alloc_local ) * sizeof (fftw_complex) );
	
        fftw_complex random_epsilon;
	
        // add randomness into the variation of each field element
        //  important to loop over tmp_vec.size() not over all alloc_local, because tmp_vec.size() <= alloc_local
        for (UINT jj = 0; jj < tmp_vec.size() ; jj++) {

          // normal distribution with mean 0 and standard deviation 1
          random_epsilon[0] = ndist(gen);
          random_epsilon[1] = ndist(gen);

          REAL lambda = tmp_vec[ jj ] / scalingfactor;

          if (lambda > 0.0) {
            KField[ jj ][0] = sqrt( lambda ) * random_epsilon[0];
            KField[ jj ][1] = sqrt( lambda ) * random_epsilon[1];
          } 
          else {
            KField[ jj ][0] = 0.0; // sqrt( -lambda ) * random_epsilon[0];
            KField[ jj ][1] = 0.0; // sqrt( -lambda ) * random_epsilon[1];
          }
        }
        
        tmp_vec.resize(0);
	
        General::log_elapsed_time( watch.elapsed(),
                                   comm,
                                   General::verbosity,
                                   "EVAL",
                                   "generateY() - part 1" );
        watch.reset();

        // fftw_plan_with_nthreads(1); // tell the FFT routine how many threads should be used

        fftw_plan plan_forward;
            
        if(dim==3)
          plan_forward = fftw_mpi_plan_dft_3d( nExt[2], 
                                               nExt[1], 
                                               nExt[0], 
                                               KField,
                                               KField,
                                               comm,
                                               FFTW_FORWARD,
                                               FFTW_ESTIMATE );
        else
          plan_forward = fftw_mpi_plan_dft_2d( nExt[1], 
                                               nExt[0], 
                                               KField,
                                               KField,
                                               comm,
                                               FFTW_FORWARD,
                                               FFTW_ESTIMATE );

        //std::cout << "FFTW_FORWARD(3d): execute forward plan..." << std::endl;
        fftw_execute( plan_forward );  
        fftw_destroy_plan( plan_forward );

        General::log_elapsed_time( watch.elapsed(),
                                   comm,
                                   General::verbosity,
                                   "FFTW3",
                                   "generateY() - part 2" );
        watch.reset();

        REAL beta_mean=0.0;  

        if(my_rank==0){
          beta_mean = fielddata.beta;
        }
            
        MPI_Bcast(&(beta_mean),1,MPI_DOUBLE,0,comm);

        
        if(my_rank==0)
          std::cout << "FFTFieldGenerator(3d): beta_mean = " << beta_mean << std::endl;

        UINT nz = fielddata.nCells[ 2 ];
        UINT ny = fielddata.nCells[ 1 ];
        UINT nx = fielddata.nCells[ 0 ];

        local_count[0] = nx;
        local_count[1] = ny;
        local_count[2] = nz;
	
        //local_count correction!
        UINT local_nField = 0;

        if( dim==3 ){

          if( (UINT)local_0_start < nz ){
            local_offset[2] = local_0_start;
            if( (UINT)local_n0 + (UINT)local_0_start <= nz ){ 
              local_count[2] = UINT(local_n0);
              local_nField = UINT(local_n0);
            } else {
              local_count[2] = nz - UINT(local_0_start);
              local_nField = nz - UINT(local_0_start);
            } 
          } else {
            // Actually, this case should not occur!
            local_count[2] = 0;
          }

        }
        else{

          if( (UINT)local_0_start < ny ){
            local_offset[1] = local_0_start;
            if( (UINT)local_n0 + (UINT)local_0_start <= ny ){ 
              local_count[1] = UINT(local_n0);
              local_nField = UINT(local_n0);
            } else {
              local_count[1] = ny - UINT(local_0_start);
              local_nField = ny - UINT(local_0_start);
            } 
          } else {
            // Actually, this case should not occur!
            local_count[1] = 0;
          }

        }
	
        UINT L,l;
        Vector<REAL> tmp( local_count.componentsproduct() );

        if(dim==3){

          for (UINT iz = 0; iz < UINT(local_n0); iz++) {
            for (UINT iy = 0; iy < nExt[1]; iy++) {
              for (UINT ix = 0; ix < nExt[0]; ix++) {
                L = General::indexconversion_3d_to_1d( iz, iy, ix, UINT(local_n0), nExt[1], nExt[0] );
                l = General::indexconversion_3d_to_1d( iz, iy, ix, UINT(local_nField), ny, nx );
                if( (iz+local_0_start < nz) && (iy < ny) && (ix < nx) ){
                  tmp[ l ] = KField[ L ][0] + beta_mean;
                }
              }
            }
          }

        }
        else{
          for (UINT iy = 0; iy < (UINT) local_n0; iy++) {
            for (UINT ix = 0; ix < nExt[0]; ix++) {
              L = General::indexconversion_2d_to_1d( iy, ix, UINT(local_n0), nExt[0] );
              l = General::indexconversion_2d_to_1d( iy, ix, local_nField, nx );
              if( (iy+local_0_start < ny) && (ix < nx) ){
                tmp[ l ] = KField[ L ][0] + beta_mean;
              }
            }
          }
        }
            
        fftw_free( KField );

        General::log_elapsed_time( watch.elapsed(),
                                   comm,
                                   General::verbosity,
                                   "EVAL",
                                   "generateY() - part 3" );

        // save the field to disc
        HDF5Tools::h5_pWrite( tmp,
                              dir.yfield_h5file,
                              "/YField",
                              fielddata.nCells,
                              local_offset,
                              local_count,
                              comm );

        // set zonation matrix right!
        l = tmp.size();
        tmp.resize(0);
        tmp.resize(l,1.0);
        HDF5Tools::h5_pWrite( tmp,
                              dir.zonation_matrix[0],
                              "/X",
                              fielddata.nCells,
                              local_offset,
                              local_count,
                              comm );

        if(comm_size==1){
          HDF5Tools::h5_Read( YFieldVector,
                              dir.yfield_h5file,
                              "/YField" );
          
          Dune::Timer watch;
          UINT N = YFieldVector.size();
          YFieldVector_Well = YFieldVector;
          
          KFieldVector.resize( N );
          for( UINT i=0; i<N; ++i ){
            KFieldVector[i] = std::exp( YFieldVector[i] );
          }
          KFieldVector_Well = KFieldVector;
          General::log_elapsed_time( watch.elapsed(),
                                     General::verbosity,
                                     "EVAL",
                                     "Evaluate K-Field once at the beginning." );

        }

      }; // End of void generateY()







      // If fielddata.newEV==false and fielddata.newField==false,
      // we try to load an existing field.
      bool load_from_file() {
        if( !load_yfield_properties( dir.yfield_properties_file, my_rank )
            )
          {
            return false;
          }
        else
          {          
	  
            if( !General::bFileExists( fielddata.location + "/YField.h5" ) )
              {
                std::cout << "Warning: File is missing: YField.h5" 
                          << std::endl;
                return false;
              }
        
            std::cout << "Found suitable Y-Field characteristics! Load the field from YField.dat" 
                      << std::endl;
	       
            if(comm_size==1){
              HDF5Tools::h5_Read( YFieldVector,
                                  dir.yfield_h5file,
                                  "/YField"
                                  );
              UINT N = YFieldVector.size();
              YFieldVector_Well = YFieldVector;
          
              KFieldVector.resize( N );
              for( UINT i=0; i<N; ++i ){
                KFieldVector[i] = std::exp( YFieldVector[i] );
              }
              KFieldVector_Well = KFieldVector;
            }
            return true;
          }

      }



      //
      // Read Y-field from HDF5 file.
      //
      void h5_Read( const std::string& filename, const std::string& groupname ){
        Vector<REAL> yfield_vector;
        HDF5Tools::h5_Read( yfield_vector,
                            filename,
                            groupname );
        import_from_vector( yfield_vector );
      }


      void import_from_vector(const Vector<REAL>& yfield_vector) {

        Dune::Timer watch;

        UINT VectorSize = 0;
        if (dim == 3) {
          const UINT L = fielddata.nCells[2];
          const UINT M = fielddata.nCells[1];
          const UINT N = fielddata.nCells[0];
          VectorSize = L * M * N;
        } 
        else {
          const UINT M = fielddata.nCells[1];
          const UINT N = fielddata.nCells[0];
          VectorSize = M * N;
        }
        if( yfield_vector.size() != VectorSize ){
          std::cout << "Warning: mismatch at import_from_vector: " << std::endl;
          std::cout << "yfield_vector.size() = " << yfield_vector.size() << std::endl;
          std::cout << "VectorSize = " << VectorSize << std::endl;
        }

        YFieldVector      = yfield_vector;
        YFieldVector_Well = yfield_vector;

        KFieldVector.resize( VectorSize );
        for (UINT l=0; l<VectorSize; l++)
          KFieldVector[l] = std::exp( yfield_vector[l] );

        KFieldVector_Well = KFieldVector;

        General::log_elapsed_time( watch.elapsed(),
                                   General::verbosity,
                                   "EVAL",
                                   "import_from_vector: Evaluate KFieldVector once at the beginning." );

      };





      //
      // Read Y-field from HDF5 file in parallel mode onto the grid.
      //
      template<typename GV>
      void h5g_pRead( GV gv,
                      const std::string& filename,
                      const std::string& groupname                      
                      ){
        Vector<REAL> local_Yfield_vector;
        Vector<UINT> local_count;
        Vector<UINT> local_offset;
        HDF5Tools::h5g_pRead( local_Yfield_vector,
                              filename,
                              groupname,
                              gv,
                              fielddata.nCells,
                              fielddata.gridsizes,
                              local_offset,
                              local_count
                              );
        if( my_rank==0 )
          std::cout << "Parallel import of Yfield data from " << filename << std::endl;

        parallel_import_from_local_vector( local_Yfield_vector,
                                           local_offset,
                                           local_count
                                           );
      }


      
      void parallel_import_from_local_vector( const Vector<REAL>& local_Yfield_vector,
                                              const Vector<UINT>& local_offset,
                                              const Vector<UINT>& local_count
                                              )
      {
        Dune::Timer watch;

        UINT Nlocal = local_Yfield_vector.size();

        LocalYFieldVector.clear();
        LocalYFieldVector_Well.clear();

        LocalYFieldVector = local_Yfield_vector; // this vector contains the Y-field data from the current hyperslab        
        LocalYFieldVector_Well = local_Yfield_vector;
        
        LocalKFieldVector.clear();
        LocalKFieldVector_Well.clear();

        LocalKFieldVector.resize( Nlocal );
        for( UINT i=0; i<Nlocal; ++i )
          LocalKFieldVector[i] = std::exp( local_Yfield_vector[i] );

        LocalKFieldVector_Well = LocalKFieldVector;

        LocalCount.resize(dim);
        LocalOffset.resize(dim);

        LocalCount = local_count; // this vector contains the dimensions (the number of cells per dimension) of the current hyperslab
        LocalOffset = local_offset; // this vector describes the distance of the current hyperslab from the origin

        General::log_elapsed_time( watch.elapsed(),
                                   General::verbosity,
                                   "EVAL",
                                   "parallel_import_from_local_vector" );

      };



      void export_to_vector(Vector<REAL>& yfield_vector ) 
      {
        UINT VectorSize = 0;
        if (dim == 3) 
          {
            const UINT L = fielddata.nCells[2];
            const UINT M = fielddata.nCells[1];
            const UINT N = fielddata.nCells[0];
            VectorSize = L * M * N;
          }
        else 
          {
            const UINT M = fielddata.nCells[1];
            const UINT N = fielddata.nCells[0];
            VectorSize = M * N;
          }
	
        yfield_vector.resize( VectorSize );
	
        for (UINT l = 0; l < VectorSize; l++)
          yfield_vector[l] = YFieldVector[l];
	
        return;
      };



      // Given a location in space, return the appropriate
      // index of the multi-dim. field array.
      std::size_t getFieldIndex( const Dune::FieldVector<REAL,dim>& xglobal ) const {

        std::size_t myIndex = 0;
        Vector<UINT> global_index;
        global_index.resize(dim);

        /*
          Translate the global coordinate 'xglobal' into its corresponding multi-dim. Yfield array index 'global_index'.
          It requires fielddata.gridsizes.
        */

        for (UINT i=0; i<dim; i++) {
          REAL fraction = xglobal[i] / fielddata.gridsizes[i];
          global_index[i] = static_cast<UINT>( fraction );
          if( std::abs( global_index[i] - fraction ) < 1E-12
              &&
              global_index[i] > 0 )
            global_index[i] -= 1;
        }

        if( comm_size == 1 ) {
          // single-processore case:
          myIndex = General::indexconversion_to_1d( global_index, fielddata.nCells );
        }
        else { // multi-processore case:

          // convert global index pair to processor-local-index pair
          Vector<UINT> local_index( dim, 0 );

          for( UINT i=0; i<dim; i++ ) {
            int localindex = global_index[i] - LocalOffset[ i ];
            if(localindex < 0)
              localindex = 0;

            local_index[i] = static_cast<UINT>( localindex );

            if( local_index[i] > 100000 ){
              std::cout << "DEBUG: local_index[" << i << "] = "
                        << local_index[i] << std::endl;
              std::cout << "DEBUG: global_index[" << i << "] = "
                        << global_index[i] << std::endl;
              std::cout << "DEBUG: LocalOffset[" << i << "] = "
                        << LocalOffset[i] << std::endl;
            }
          }

          myIndex = General::indexconversion_to_1d( local_index, LocalCount );
        }

        return myIndex;
      }




      template<typename DFV>
      void evaluateField( const DFV& xglobal, 
                          REAL& output,
                          bool lnK=false ) const 
      {
        /* 
           Translate the global coordinate 'xglobal' into its multi-dim. Yfield array index 'l'. 
        */
        UINT l = getFieldIndex( xglobal );

        if( comm_size == 1 ) {

          if( lnK )
            output = YFieldVector[l];
          else
            output = KFieldVector[l];

        }
        else {

          if( lnK )
            output = LocalYFieldVector[l];
          else
            output = LocalKFieldVector[l];

        }

      }




      template<typename DFV>
      void evaluateFieldAndWell( const DFV& xglobal, 
                                 REAL& output,
                                 REAL& output_Well,
                                 bool lnK=false ) const 
      {
        /* 
           Translate the global coordinate 'xglobal' into its multi-dim. Yfield array index 'l'. 
        */
        UINT l = getFieldIndex( xglobal );

        if( comm_size == 1 ) {
          if( lnK ) {
            output      = YFieldVector[l];
            output_Well = YFieldVector_Well[l];
          }
          else{
            output      = KFieldVector[l];
            output_Well = KFieldVector_Well[l];
          }
        }
        else {
          if( lnK ){
            output      = LocalYFieldVector[l];
            output_Well = LocalYFieldVector_Well[l];
          }
          else{
            output      = LocalKFieldVector[l];
            output_Well = LocalKFieldVector_Well[l];
          }
        }

      }


      template<typename DFV>
      void evaluateY( const DFV& xglobal, REAL& output ) const 
      {
        evaluateField( xglobal, output, true );
      }

      template<typename DFV>
      void evaluateYW( const DFV& xglobal, REAL& output, REAL& output_well ) const 
      {
        evaluateFieldAndWell( xglobal, output, output_well, true );
      }

      template<typename DFV>
      void evaluateK( const DFV& xglobal, REAL& output ) const 
      {
        evaluateField( xglobal, output, false );    
      }

      template<typename DFV>
      void evaluateKW( const DFV& xglobal, REAL& output, REAL& output_well ) const 
      {
        evaluateFieldAndWell( xglobal, output, output_well, false );
      }



      void init()
      {
        if(my_rank==0)
          std::cout << "------------------------------------\ngenerate_Yfield()" << std::endl;
    
        // Be aware of this: 
        // The virtual YField grid is a structured grid (FFT requires this!)
        // no matter what kind of grid we are using
        // in DUNE for the solution of the PDEs.

        UINT nAllCells = fielddata.nCells[0] * fielddata.nCells[1];
        if (dim == 3)
          nAllCells *= fielddata.nCells[2];

        int generate_new=0;

        if( my_rank == 0 )
          {
            General::createDirectory( fielddata.location );

            //loading is in HDF5 or *.dat, (see FFTFieldGenerator.hh) 
            if( !fielddata.newField 
                && 
                !fielddata.newEV 
                &&
                this->load_from_file( ) ) // loads data only if sequential program, in parallel only consistency check
              {
                std::cout << "Reading existing Y-field."
                          << std::endl;		
              }
            else 
              {
                std::cout << "Generating new data ..."
                          << std::endl;
                generate_new = 1;
              }
	  
          } // endif helper.rank() == 0
	
        MPI_Bcast(&(generate_new), 1, MPI_INT, 0, comm );


        if(generate_new){
      
          this->generateY( fielddata.newEV,
                           fielddata.newField
                           );
      
          if( my_rank == 0 )
            std::cout << "Generating new Y-field done." << std::endl;
      
        }
    
      }; // end of void init()





      //
      // This method maps the multi-dim. field array to a 1d vector, 
      // following the  order of the grid cells.
      //
      template<typename GV>
      void export_field_to_vector_on_grid( const GV& gv,
                                           Vector<REAL>& log_conductivity_field,
                                           Vector<REAL>& well_conductivity_field
                                           ) const {

        Dune::Timer watch;

        // Get the conductivity field on the grid for the vtkwriter
        typedef typename GV::Grid GRIDTYPE; // get the grid-type out of the gridview-type

        typedef typename GV::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementLeafIterator;


        // make a mapper for codim 0 entities in the leaf grid
        Dune::LeafMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout>
          mapper(gv.grid()); // get the underlying hierarchical grid out ot the gridview


        const int nGridCells = mapper.size();

        //std::cout << my_rank << ": mapper.size() = " << nGridCells << std::endl;

        log_conductivity_field.resize( nGridCells );
        well_conductivity_field.resize( nGridCells );

        if( my_rank == 0 )
          std::cout << "Plot Yfield: Loop over leaf elements..." << std::endl;
    
        const typename GV::IndexSet& indexset = gv.indexSet();

        for (ElementLeafIterator it = gv.template begin<0,Dune::All_Partition> ()
               ; it != gv.template end<0,Dune::All_Partition> (); ++it) {
          int index = indexset.index(*it);

          //std::cout << my_rank << ": index = " << index << std::endl;
      
          Dune::FieldVector<REAL,dim> xglobal = it->geometry().center();
          REAL logK = 0.0;
          REAL logK_Well = 0.0;

          this->evaluateYW(xglobal,logK,logK_Well);

          log_conductivity_field[index] = logK;
          well_conductivity_field[index] = logK_Well;

        }


        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   General::verbosity,
                                   "EVAL",
                                   "export_field_to_vector_on_grid" );

      } // end of export_field_to_vector_on_grid()




      //
      // This method can be used to create a VTK plot of the field.
      //
      template<typename GV>
      void plot2vtu( const GV& gv
                     , const std::string filename
                     , const std::string title
                     , int nSubSampling = 0
                     , bool bWells = false
                     ) const {

        Vector<REAL> log_conductivity_field;
        Vector<REAL> well_conductivity_field;

        export_field_to_vector_on_grid( gv, 
                                        log_conductivity_field,
                                        well_conductivity_field );

        VTKPlot::output_vector_to_vtu( gv, 
                                       log_conductivity_field,
                                       filename,
                                       title,
                                       nSubSampling
                                       );

        if( bWells )
          VTKPlot::output_vector_to_vtu( gv, 
                                         well_conductivity_field,
                                         filename + "_Well",
                                         title + "_Well",
                                         nSubSampling
                                         );

      }





      //
      // This method implements various variogram models.
      //
      inline REAL autocovariancefunction( const REAL variance, 
                                        const std::vector<REAL> x, 
                                        const std::vector<REAL> lambda, 
                                        const std::string model  ) {

        //UINT dim = x.size();
        if( variance == 0.0 )
          return 0.0;

        REAL C1 = 0.0;
        REAL sum = 0.0;
        for(UINT i=0; i<dim; i++) {
          sum += (x[i] * x[i]) / (lambda[i] * lambda[i]);
        }
        if( !strcmp( model.c_str(), "spherical" ) ) {
          REAL h_eff = sqrt( sum );
          //std::cout << "...spherical variogram..." << std::endl;
          REAL alpha=2.2;
          if (h_eff > alpha)
            C1 = 0.0;
          else
            C1 = variance * (1.0 - 1.5 * h_eff/alpha + 0.5 * std::pow(h_eff/alpha, 3.0));
        }
        else if( !strcmp( model.c_str(), "power" ) ){
          C1 = variance * (1.0 - sum) * exp(-sum);
        }
        else if( !strcmp( model.c_str(), "expomodulus" ) ){
          REAL sum_modulus = 0.0;
          for(UINT i=0; i<dim; i++) {
            sum_modulus += std::abs(x[i]) / lambda[i];
          }
          //std::cout << "...exponential exponential with modulus..." << std::endl;
          C1 = variance * exp(- sum_modulus);
        }
        else if( !strcmp( model.c_str(), "exponential" ) ){
          REAL h_eff = sqrt( sum );
          //std::cout << "...exponential variogram..." << std::endl;
          C1 = variance * exp(-h_eff);
        }
        else{ // gaussian
              //std::cout << "...gaussian variogram..." << std::endl;
          C1 = variance * exp(-sum);
        }
        return C1;
      }



      //
      // This function stores the properties of the generated random field to a file called "YField.dat".
      //
      inline bool save_yfield_properties( const std::string filename )
      {
        std::ofstream filestr( filename.c_str(), std::ios::out );
        if( filestr.is_open() ) {
        
          std::cout << "Saving Y-Field properties to "
                    << filename
                    << std::endl;

          // DO NOT insert a blank into "dim="
          filestr << "dim= " << dim << std::endl;

          filestr << "extensions= ";
          for(UINT i=0; i<dim; i++)
            filestr << " " << fielddata.extensions[i];
          filestr << std::endl;

          filestr << "nCells= ";
          for(UINT i=0; i<dim; i++)
            filestr << " " << fielddata.nCells[i];
          filestr << std::endl;

        
          filestr << "VM= " << fielddata.variogramModel;
          filestr << std::endl;

          filestr << "beta= " << fielddata.beta;
          filestr << std::endl;

          filestr << "variance= " << fielddata.fieldVariance;
          filestr << std::endl;

          filestr << "EF= " << fielddata.embeddingFactor;
          filestr << std::endl;

          filestr << "correlations= ";
          for(UINT i=0; i<dim; i++)
            filestr << " " << fielddata.correlations[i];
          filestr << std::endl;

          filestr << "seed= " << fielddata.seed;
          filestr << std::endl;

          filestr.close();
        }
        else
          {
            std::cout << "Error: Could not save Y-Field properties to file "
                      << filename.c_str()
                      << std::endl;
            return false;
          }

        return 0;
      }



      //
      // This function reads the properties of a random field generated in a previous run.
      // The field properties are read from a file called "YField.dat".
      //
      inline bool load_yfield_properties( const std::string filename, const int my_rank )
      {    
        std::ifstream filestr( filename.c_str(), std::ios::in );
        if( filestr.is_open() ) {

          if( my_rank==0 )
            std::cout << "Loading " << filename << std::endl;
      
          std::string str;

          UINT read_dim;
          filestr >> str >> read_dim;
          if( dim != read_dim )
            {
              if(my_rank==0)
                std::cout << "Drop obsolete Y-field with : dim = " << read_dim << std::endl;
              return false;
            }

          CTYPE read_extensions[3];
          if(dim == 3){
            filestr >> str >> read_extensions[0] >> read_extensions[1] >> read_extensions[2];
            if( read_extensions[0] != fielddata.extensions[0] 
                || 
                read_extensions[1] != fielddata.extensions[1] 
                || 
                read_extensions[2] != fielddata.extensions[2] )
              {
                if( my_rank==0 )
                  std::cout << "New 3D random field with Lx = " << fielddata.extensions[0]
                            << " Ly = " << fielddata.extensions[1]
                            << " Lz = " << fielddata.extensions[2]
                            << std::endl;
                return false;
              }
          }
          else
            {
              filestr >> str >> read_extensions[0] >> read_extensions[1];
              if( read_extensions[0] != fielddata.extensions[0] || 
                  read_extensions[1] != fielddata.extensions[1] )
                {
                  if( my_rank==0 )
                    std::cout << "New 2D random field with Lx = " << fielddata.extensions[0]
                              << " Ly = " << fielddata.extensions[1]
                              << std::endl;
                  return false;
                }
            }

          UINT read_nCells[3];
          if(dim == 3){
            filestr >> str >> read_nCells[0] >> read_nCells[1] >> read_nCells[2];
            if( read_nCells[0] != fielddata.nCells[0] || 
                read_nCells[1] != fielddata.nCells[1] || 
                read_nCells[2] != fielddata.nCells[2] )
              {
                if( my_rank==0 )
                  std::cout << "New 3D random field with Nx = " << fielddata.nCells[0]
                            << " Ny = " << fielddata.nCells[1]
                            << " Nz = " << fielddata.nCells[2]
                            << std::endl;
                return false;
              }
          }
          else
            {
              filestr >> str >> read_nCells[0] >> read_nCells[1];
              if( read_nCells[0] != fielddata.nCells[0] || 
                  read_nCells[1] != fielddata.nCells[1] )
                {
                  if( my_rank==0 )
                    std::cout << "New 2D random field with Nx = " << fielddata.nCells[0]
                              << " Ny = " << fielddata.nCells[1]
                              << std::endl;
                  return false;
                }
            }

        


          std::string read_model;
          filestr >> str >> read_model;
          if( fielddata.variogramModel != read_model ) {
            if(my_rank==0)
              std::cout << "Drop obsolete Y-field with : model = " << read_model.c_str() << std::endl;
            return false;
          }

          REAL read_beta;
          filestr >> str >> read_beta;
          if( fielddata.beta != read_beta ) {
            if(my_rank==0)
              std::cout << "Drop obsolete Y-field with : beta = " << read_beta << std::endl;
            return false;
          }

          REAL read_variance;
          filestr >> str >> read_variance;
          if( fielddata.fieldVariance != read_variance ) {
            if(my_rank==0)
              std::cout << "Drop obsolete Y-field with : variance = " << read_variance << std::endl;
            return false;
          }

          REAL read_embedding_factor;
          filestr >> str >> read_embedding_factor;
          if( fielddata.embeddingFactor != read_embedding_factor ) {
            if(my_rank==0)
              std::cout << "Drop obsolete Y-field with : embedding_factor = " << read_embedding_factor << std::endl;
            return false;
          }


          REAL read_correlation_lambda[3];
          if(dim == 3){
            filestr >> str >> read_correlation_lambda[0] >> read_correlation_lambda[1] >> read_correlation_lambda[2];
            if( fielddata.correlations[0] != read_correlation_lambda[0]
                || 
                fielddata.correlations[1] != read_correlation_lambda[1]
                ||
                fielddata.correlations[2] != read_correlation_lambda[2]
                )
              {
                if( my_rank==0 )
                  std::cout << "New 3D random field with lambda_x = " << fielddata.correlations[0]
                            << " lambda_y = " << fielddata.correlations[1]
                            << " lambda_z = " << fielddata.correlations[2]
                            << std::endl;
                return false;
              }
          }
          else
            {
              filestr >> str >> read_correlation_lambda[0] >> read_correlation_lambda[1];
              if( fielddata.correlations[0] != read_correlation_lambda[0]
                  || 
                  fielddata.correlations[1] != read_correlation_lambda[1] )
                {
                  if( my_rank==0 )
                    std::cout << "New 2D random field with lambda_x = " 
                              << fielddata.correlations[0]
                              << " lambda_y = " << fielddata.correlations[1]
                              << std::endl;
                  return false;
                }
            }

          REAL read_seed;
          filestr >> str >> read_seed;
          if( fielddata.seed != read_seed ) {
            if(my_rank==0)
              std::cout << "Drop obsolete Y-field with : seed = " << read_seed << std::endl;
            return false;
          }

          filestr.close();
        }
        else
          {
            if( my_rank == 0 ){
              std::cout << "No existing field properties found for the given inputs! Generate new field." << std::endl;
            }
            return false;
          }

        return true;
      }




      template<typename GV_GW,typename IDT>
      void setWellConductivities( const GV_GW& gv_gw, const IDT& inputdata ){

        Dune::Timer watch;
        logger << "setWellConductivities()" << std::endl;

        const UINT nSetups = inputdata.setups.size();
        for( UINT iSetup=0; iSetup < nSetups; iSetup++ ) {
          UINT nWells = inputdata.listOfAllWellCenters[iSetup].size();
          for( UINT iWell=0; iWell<nWells; iWell++ ) {
            for( UINT iCenter=0;
                 iCenter<inputdata.listOfAllWellCenters[iSetup][iWell].size();
                 iCenter++ ){

              Dune::FieldVector<REAL,dim> wellCenter = inputdata.listOfAllWellCenters[iSetup][iWell][iCenter].second;

              bool bWellIsInsidePartition = false;
              if( inputdata.parallel_fine_tuning.maxNComm > 1
#ifdef USE_ALUGRID
                  ||
                  inputdata.domain_data.alugrid_maxsteps > 0
#endif
                  ){
                // Parallel-in-parallel needs this extra check:
                if( General::doesPointBelongToGridPartition( wellCenter, gv_gw ) )
                  bWellIsInsidePartition = true;
                else{

                  //std::cout << "FTTFieldGenerator::setWellConductivities(), Parallel-in-parallel extra check: Well-center " <<  wellCenter
                  //<< " may not belong to rank " << my_rank << " anymore inside new pool."
                  //<< std::endl;

                }
              }
              else {
                // The existing wells are always inside the current partition as long as the partitions have not change since
                // the calibration at the beginning!
                bWellIsInsidePartition = true;
              }


              if( bWellIsInsidePartition ) {

#ifdef DEBUG_LOG
                std::cout << "P" << my_rank
                          << ": setWellConductivities() at: " <<  wellCenter
                          << std::endl;
#endif

                UINT l = getFieldIndex( inputdata.listOfAllWellCenters[iSetup][iWell][iCenter].second );
                REAL K_well = inputdata.listOfAllWellCenters[iSetup][iWell][iCenter].first;
                if( comm_size == 1 ){
                  KFieldVector_Well[l] = K_well;
                  YFieldVector_Well[l] = std::log( K_well );
                }
                else{
                  LocalKFieldVector_Well[l] = K_well;
                  LocalYFieldVector_Well[l] = std::log( K_well );
                }

              }

            }
          }
        }

        General::log_elapsed_time( watch.elapsed(),
                                   gv_gw.comm(),
                                   General::verbosity,
                                   "EVAL",
                                   "setWellConductivities()" );

      } // void setWellConductivities(...)


      std::string getFilename() const {
        return dir.yfield_h5file;
      }


    }; // class FFTFieldGenerator


  }
}

#endif	/* FFT_FIELD_GENERATOR_HH */

