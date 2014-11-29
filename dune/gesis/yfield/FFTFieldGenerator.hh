/*
 * File:   FFTFieldGenerator.hh
 * Author: ngo
 *
 * Created on June 25, 2010, 3:43 PM
 * Last Modified on Friday July 28, 2014
 */

#ifndef DUNE_GESIS_FFT_FIELD_GENERATOR_HH
#define	DUNE_GESIS_FFT_FIELD_GENERATOR_HH


#include <fftw3.h>
#include <fftw3-mpi.h>

#include <time.h>                      // define time()

#include <random>      // C++11 random number generator (This works with g++-4.8.)


#include <dune/gesis/common/io/HDF5Tools.hh>

extern CLogfile logger;


namespace Dune {
  namespace Gesis {

    template<typename IDT,typename REAL,int dim>
    class FFTFieldGenerator {

    private:

      //UINT dim;
      UINT nzones;
      const IDT &inputdata;
      const IO_Locations &dir;

      // MPI communicator and information about rank and size
      // const Dune::MPIHelper &helper;

      const MPI_Comm &comm;
      int my_rank,comm_size;

      std::vector< Vector<UINT> > nCells_ExtendedDomain;
      //  std::vector<std::string> R_YY_filename;
      //  std::vector<std::string> ev_filename;
      //  std::vector<std::string> X_filename;

      // Note:
      // Always remeber that Y=log(K) where log is the natural logarithm.
      // We store the Y field inside the HDF5 files, because the Y field is
      // being computed by the random field generator and by the inversion
      // scheme.
      // This means that we need to compute K=std::exp(Y) each time we read from
      // such a HDF5 file! This is done once during
      // parallel_import_from_local_vector(...)
      // and
      // each time after h5g_Read(...)
      //
      //
      Vector<REAL> YFieldVector;
      Vector<REAL> YFieldVector_Well;
      Vector<REAL> KFieldVector;
      Vector<REAL> KFieldVector_Well;

      Vector<UINT> nCells_zones;

      Vector<REAL>  LocalKFieldVector;
      Vector<REAL>  LocalKFieldVector_Well;
      Vector<REAL>  LocalYFieldVector;
      Vector<REAL>  LocalYFieldVector_Well;

      Vector<UINT>  LocalCount;
      Vector<UINT>  LocalOffset;

    public:

      UINT size(){
        return KFieldVector.size();
      };
      UINT localSize(){
        return LocalKFieldVector.size();
      };

      // first constructor:
      FFTFieldGenerator( const IDT &inputdata_
                         , const IO_Locations &dir_
                         )
        :
        inputdata(inputdata_),dir(dir_)
      {
        nzones=inputdata.yfield_properties.nz;
        //dim=inputdata.dim;

        nCells_ExtendedDomain.resize(nzones);
        for( UINT i=0; i<nzones; ++i ){
          nCells_ExtendedDomain[i].resize(dim);
        }

        nCells_zones.resize(nzones);

        for (UINT ii =0; ii<nzones;ii++){
          nCells_ExtendedDomain[ii].resize(dim);
          if(nzones>1){
            if(ii==0)
              nCells_zones[ii]=floor(inputdata.yfield_properties.zones[ii].bottom/inputdata.domain_data.gridsizes[dim-1]);
            else if (ii==nzones-1)
              nCells_zones[ii]=inputdata.domain_data.nCells[ dim-1 ]-floor(inputdata.yfield_properties.zones[ii-1].bottom/inputdata.domain_data.gridsizes[dim-1]);
            else
              nCells_zones[ii]=floor(inputdata.yfield_properties.zones[ii].bottom/inputdata.domain_data.gridsizes[dim-1])-floor(inputdata.yfield_properties.zones[ii-1].bottom/inputdata.domain_data.gridsizes[dim-1]);
          }else
            nCells_zones[ii] = inputdata.domain_data.nCells[ dim-1 ];
        }


        if( my_rank==0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
          std::cout << "Zonation:" << std::endl;
          for (UINT ii =0; ii<nzones;ii++){
            std::cout << "nCells_zones[" << ii <<"] = " << nCells_zones[ii] << std::endl;
          }
        }

#ifdef FFT_THREADS
        //initialization for the multi-threaded fft
        fftw_init_threads();
#endif
      };


      // second constructor:
      FFTFieldGenerator( const IDT &inputdata_
                         , const IO_Locations &dir_
                         , const MPI_Comm &comm_
                         )
        :
        inputdata( inputdata_ )
        , dir(dir_)
        //, helper( helper_ )
        //, comm( helper_.getCommunicator() )
        , comm(comm_)
      {
        //setting my_rank  (rank within the stored communicator)
        MPI_Comm_rank(comm,&my_rank);

        //setting comm_size (size of the stored communicator)
        MPI_Comm_size(comm,&comm_size);

        nzones=inputdata.yfield_properties.nz;
        //dim=inputdata.dim;

        nCells_ExtendedDomain.resize(nzones);
        for( UINT i=0; i<nzones; ++i ){
          nCells_ExtendedDomain[i].resize(dim);
        }

        nCells_zones.resize(nzones);
	for (UINT ii =0; ii<nzones;ii++){
          if(nzones>1){
            if(ii==0)
              nCells_zones[ii]=floor(inputdata.yfield_properties.zones[ii].bottom/inputdata.domain_data.virtual_gridsizes[dim-1]);
            else if (ii==nzones-1)
              nCells_zones[ii]=inputdata.domain_data.nCells[ dim-1 ]-floor(inputdata.yfield_properties.zones[ii-1].bottom/inputdata.domain_data.virtual_gridsizes[dim-1]);
            else
              nCells_zones[ii]=floor(inputdata.yfield_properties.zones[ii].bottom/inputdata.domain_data.virtual_gridsizes[dim-1])-floor(inputdata.yfield_properties.zones[ii-1].bottom/inputdata.domain_data.virtual_gridsizes[dim-1]);
          }else
            nCells_zones[ii] = inputdata.domain_data.nCells[ dim-1 ];
        }


        if( my_rank==0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
          std::cout << "Zonation:" << std::endl;
          for (UINT ii =0; ii<nzones;ii++){
            std::cout << "nCells_zones[" << ii <<"] = " << nCells_zones[ii] << std::endl;
          }
        }

        /* -> not needed!!!! makes it for later use only complitcated! -> no special case for a homogenous field!!!
         * REASONS: - seldom used, mostly for test cases
         *          - makes the inversion to complicated. you need always this two cases! homo-or hetero!
         if( inputdata.yfield_properties.variance == 0 )
         {
         KFieldVector.resize( 1 );
         int seed = (int) time(0); // create seed out ot the current time
         logger << std::endl;
         logger << "Seed for homogeneous field with new random mean = " << seed << std::endl;
         logger << std::endl;
         StochasticLib1 stochastic( seed ); // make instance of random library
         #ifdef SEED_OFF
         REAL beta_mean = inputdata.yfield_properties.beta;
         #else
         REAL beta_mean = inputdata.yfield_properties.beta + stochastic.Normal(0.0, inputdata.yfield_properties.qbb_y);
         #endif
         KFieldVector[ 0 ] = beta_mean;
         }*/
#ifdef FFT_THREADS
        //initialization for the multi-threaded fft
        fftw_init_threads();
#endif
      };

      // This function can be used for dim=2 or dim=3


      // compute the size of the extended domain required for circulant embedding
      void extend_domain()
      {
        for(UINT ii=0;ii<nzones;ii++){
          for (UINT i = 0; i < dim; i++) {
            REAL b = 0;
            if(i==dim-1 && nzones>1){
              b = nCells_zones[ii] * inputdata.domain_data.virtual_gridsizes[i];
            } else {
              b = inputdata.domain_data.extensions[i];
            }
            REAL a = 2.0 * b;
            REAL f = inputdata.yfield_properties.zones[ii].embedding_factor;
            REAL c = inputdata.yfield_properties.zones[ii].correlation_lambda[i];
            REAL d = inputdata.domain_data.virtual_gridsizes[i];

            nCells_ExtendedDomain[ii][i] = std::ceil( std::max( a , b + f * c ) / d );

            logger << "Extended domain for zone #"<<ii+1<<" has " << nCells_ExtendedDomain[ii][i] << " elements in direction " << i << std::endl;
            //std::cout << "Extended domain for zone #"<<ii+1<<" has " << nCells_ExtendedDomain[ii][i] << " elements in direction " << i << std::endl;

          }
        }
      };


      void export_nCellsExt( std::vector< Vector<UINT> >& nCellsExt ) const
      {

        nCellsExt.resize( nzones );
        for(UINT ii=0;ii<nzones;ii++){
          nCellsExt[ii].resize(dim);
          for( UINT i=0; i<dim; i++ )
            nCellsExt[ii][i] = nCells_ExtendedDomain[ii][i];
        }
      }



      // This function can be used for dim=3 only!
      void generate_eigenvalues_3d( ) {

        Dune::Timer watch;

        for(UINT ii=0;ii<nzones;ii++){

          const UINT Nz = nCells_ExtendedDomain[ii][2]; // number of cells in z-direction = #levels
          const UINT Ny = nCells_ExtendedDomain[ii][1]; // number of cells in y-direction = #rows
          const UINT Nx = nCells_ExtendedDomain[ii][0]; // number of cells in x-direction = #columns

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

          const REAL dx = inputdata.domain_data.virtual_gridsizes[0];
          const REAL dy = inputdata.domain_data.virtual_gridsizes[1];
          const REAL dz = inputdata.domain_data.virtual_gridsizes[2];

          //bool bWisdomFileExists = false;

          /*
           * Definition of the block on each processor for the parallel FFT
           *   Idea: Maybe later another data distribution. "fftw_mpi_local_size_many" is a more general routine
           *         BUT: don't forget to load also the Lambdas in the right format
           */
          ptrdiff_t alloc_local, local_n0, local_0_start;

          // get FFTW data
          alloc_local = fftw_mpi_local_size_3d(Nz,Ny,Nx, comm,
                                               &local_n0, &local_0_start);


          /*
          // Only for the special tensors 'z_mirror', 'y_mirror' and 'x_mirror' remember this cycle :
          // Nz - Ny - Nx - Nz - Ny - ...
          // iz - iy - ix - iz - iy - ...

          std::vector<std::vector< std::vector<REAL> > > z_mirror; // Nz x Ny x Nx
          z_mirror.resize( Nz );
          for (UINT iz = 0; iz < Nz; iz++) {
          z_mirror[ iz ].resize( Ny );
          for (UINT iy = 0; iy < Ny; iy++) {
          z_mirror[ iz ][ iy ].resize( Nx );
          for( UINT ix = 0; ix < Nx; ix++ ) {
          z_mirror[ iz ][ iy ][ ix ] = std::min((REAL) iz, (REAL) (Nz-iz)) * dz;
          }
          }
          }

          std::vector<std::vector< std::vector<REAL> > > y_mirror; // Ny x Nx x Nz
          y_mirror.resize( Ny );
          for( UINT iy = 0; iy < Ny; iy++ ) {
          y_mirror[ iy ].resize( Nx );
          for( UINT ix = 0; ix < Nx; ix++ ) {
          y_mirror[ iy ][ ix ].resize( Nz );
          for( UINT iz = 0; iz < Nz; iz++) {
          y_mirror[ iy ][ ix ][ iz ] = std::min( (REAL) iy, (REAL) (Ny-iy) ) * dy;
          }
          }
          }

          std::vector<std::vector< std::vector<REAL> > > x_mirror; // Nx x Nz x Ny
          x_mirror.resize( Nx );
          for (UINT ix = 0; ix < Nx; ix++) {
          x_mirror[ ix ].resize( Nz );
          for (UINT iz = 0; iz < Nz; iz++) {
          x_mirror[ ix ][ iz ].resize( Ny );
          for( UINT iy = 0; iy < Ny; iy++ ) {
          x_mirror[ ix ][ iz ][ iy ] = std::min((REAL) ix, (REAL) (Nx-ix)) * dx;
          }
          }
          }



          std::vector< std::vector< std::vector<REAL> > > R_YY; // Nz x Ny x Nx  ( iz, iy, ix)
          */
          std::vector<REAL> R_YY(alloc_local);

          Vector<REAL> X( dim, 0.0 );
          Vector<REAL> Lambda( inputdata.yfield_properties.zones[ii].correlation_lambda );
          UINT l;

          // build up R_YY
          for (UINT iz = 0; iz < (UINT)local_n0; iz++) {
            for (UINT iy = 0; iy < Ny; iy++) {
              for (UINT ix = 0; ix < Nx; ix++) {

                X[0] = std::min( (REAL) ix, (REAL) (Nx-ix) ) * dx;
                X[1] = std::min( (REAL) iy, (REAL) (Ny-iy) ) * dy;
                X[2] = std::min( (REAL) (iz+local_0_start), (REAL) (Nz-(iz+local_0_start)) ) * dz;

                l = General::indexconversion_3d_to_1d( iz, iy, ix, local_n0, Ny, Nx );
                R_YY[ l ] = General::autocovariancefunction(
                                                            inputdata.yfield_properties.zones[ii].variance
                                                            , X
                                                            , Lambda
                                                            , inputdata.yfield_properties.zones[ii].model
                                                            );
              }
            }
          }
          /*
            x_mirror.clear();
            y_mirror.clear();
            z_mirror.clear();
          */

          // writing of the R_YY in HDF5
          Vector<UINT> local_count(3,0), local_offset(3,0);
          local_count[0]=nCells_ExtendedDomain[ii][0];
          local_count[1]=nCells_ExtendedDomain[ii][1];
          local_count[2]=local_n0;

          local_offset[2]=local_0_start;

          General::log_elapsed_time( watch.elapsed(),
                                     comm,
                                     inputdata.verbosity,
                                     "EVAL",
                                     "generate_eigenvalues_3d() - part 1" );

          HDF5Tools::h5_pWrite( R_YY
                                , dir.R_YY_h5file[ii]
                                , "/R_YY"
                                , inputdata
                                , nCells_ExtendedDomain[ii]
                                , local_offset
                                , local_count
                                , comm );

          watch.reset();

          //logger << "Use FFTW_FORWARD( 3d ) to generate eigenvalues... allocate fftw_complex vector of size "
          //	   << Nz * Ny * Nx << std::endl;

          fftw_complex *Eigenvalues;
          Eigenvalues = (fftw_complex*) fftw_malloc( (alloc_local ) * sizeof (fftw_complex) );

          for( UINT iz = 0; iz < (UINT)local_n0; iz++ ) {
            for( UINT iy = 0; iy < Ny; iy++ ) {
              for( UINT ix = 0; ix < Nx; ix++ ) {

                UINT l = General::indexconversion_3d_to_1d( iz, iy, ix, local_n0, Ny, Nx );
                Eigenvalues[ l ][0] = R_YY[ l ];                  // real part
                Eigenvalues[ l ][1] = 0.0;                        // imaginary part

              }
            }
          }


#ifdef FFT_THREADS
          fftw_plan_with_nthreads(FFT_THREADS); // tell the FFT routine how many threads should be used
#endif
          fftw_plan plan_forward;
          //logger << "FFTW_FORWARD(3d): Create 3d forward plan using FFTW_ESTIMATE flag..." << std::endl;
          plan_forward = fftw_mpi_plan_dft_3d( Nz, Ny, Nx, Eigenvalues, Eigenvalues, comm, FFTW_FORWARD, FFTW_ESTIMATE );
          //logger << "FFTW_FORWARD(3d): executing plan..." << std::endl;

          fftw_execute(plan_forward);
          fftw_destroy_plan(plan_forward);


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

          if( inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
            std::cout << "Properties of Syy: " << std::endl;
            std::cout << alloc_local << " eigenvalues" << std::endl;
            std::cout << nPositive << " eigenvalues > 1E-6" << std::endl;
            std::cout << nNegative << " eigenvalues < -1E-6" << std::endl;
            std::cout << nNearZero << " eigenvalues near zero" << std::endl;
          }

          General::log_elapsed_time( watch.elapsed(),
                                     comm,
                                     inputdata.verbosity,
                                     "FFTW",
                                     "generate_eigenvalues_3d() - part 2" );

          // save FFT of R_YY!
          HDF5Tools::h5_pWrite( tmp_vec
                                , dir.EV_h5file[ii]
                                , "/FFT_R_YY"
                                , inputdata
                                , nCells_ExtendedDomain[ii]
                                , local_offset
                                , local_count
                                , comm );

          fftw_free(Eigenvalues);

        } // END: for-loop over the zones

        if(my_rank==0)
          General::save_kfield_properties(
                                          dir.kfield_properties_file
                                          , dim
                                          , inputdata.domain_data.extensions
                                          , inputdata.domain_data.nCells
                                          , inputdata.yfield_properties );

      }; // End of generate_eigenvalues_3d()



      void generate3d( const bool bNewEigenvalues,
                       const bool bNewField,
                       const bool CR=false )
      {

        if( bNewEigenvalues
            ||
            ( !bNewEigenvalues
              &&
              !General::load_kfield_properties(
                                               dir.kfield_properties_file
                                               , dim
                                               , inputdata.domain_data.extensions
                                               , inputdata.domain_data.nCells
                                               , inputdata.yfield_properties
                                               , my_rank
                                               )
              )
            )
          {
            generate_eigenvalues_3d();
          }
        else
          logger << "Found suitable Y-Field characteristics!." << std::endl;


        if( !bNewField ){
          if( my_rank == 0 ){
            std::cout << "Warning: New Y-field will be generated due to changed Y-field characteristics." << std::endl;
          }
        }

        for(UINT ii=0;ii<nzones;ii++){

          // size of the extended domain
          const UINT Nz = nCells_ExtendedDomain[ii][2];
          const UINT Ny = nCells_ExtendedDomain[ii][1];
          const UINT Nx = nCells_ExtendedDomain[ii][0];

          REAL scalingfactor = REAL( Nz * Ny * Nx );

          /*
           * Definition of the block on each processor for the parallel FFT
           *   Idea: Maybe later another data distribution. "fftw_mpi_local_size_many" is a more general routine
           *         BUT: don't forget to load also the Lambdas in the right format
           */
          ptrdiff_t alloc_local, local_n0, local_0_start;

          // get the FFT data
          alloc_local = fftw_mpi_local_size_3d(Nz ,Ny,Nx, comm,
                                               &local_n0, &local_0_start);


          //loading in HDF5
          Vector<REAL> tmp_vec; // temporary vector for the loaded eigenvalues

          Vector<UINT> local_count(3,0), local_offset(3,0);
          local_count[0]=nCells_ExtendedDomain[ii][0];
          local_count[1]=nCells_ExtendedDomain[ii][1];
          local_count[2]=local_n0;

          local_offset[2]=local_0_start;

          HDF5Tools::h5_pRead( tmp_vec
                               , dir.EV_h5file[ii]
                               , "/FFT_R_YY"
                               , inputdata
                               , local_offset
                               , local_count
                               , comm
                               );

          Dune::Timer watch;
          watch.reset();


          // initialize pseudo-random generator
          std::random_device rd;
          std::mt19937_64 gen; // 64-bit Mersenne Twister
          std::normal_distribution<> ndist(0,1); // mean = 0, sigma = 1

          int seed = inputdata.yfield_properties.random_seed;
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

#ifdef FFT_THREADS
          fftw_plan_with_nthreads(FFT_THREADS); // tell the FFT routine how many threads should be used
#endif
          fftw_plan plan_backward;

          //logger << "FFTW_BACKWARD(3d): create backward plan using FFTW_ESTIMATE..." << std::endl;
          plan_backward = fftw_mpi_plan_dft_3d( Nz, Ny, Nx, KField, KField, comm, FFTW_BACKWARD, FFTW_ESTIMATE );

          //logger << "FFTW_BACKWARD(3d): execute backward plan..." << std::endl;
          fftw_execute( plan_backward );
          fftw_destroy_plan( plan_backward );


          REAL beta_mean=0.0;

          if(!CR){
            if(my_rank==0){
              if(inputdata.yfield_properties.zones[ii].qbb_y>0.0){
                REAL rv = ndist(gen);
                beta_mean = inputdata.yfield_properties.zones[ii].beta + rv * std::sqrt( inputdata.yfield_properties.zones[ii].qbb_y );
              }
              else
                beta_mean= inputdata.yfield_properties.zones[ii].beta;
            }
          }



          MPI_Bcast(&(beta_mean),1,MPI_DOUBLE,0,comm);

          logger << "FFTFieldGenerator(3d): beta_mean = " << beta_mean << std::endl;

          //correction du to previous zones
          UINT zone_offset=0;
          for(UINT jj=0; jj<ii; jj++)
            zone_offset+=nCells_zones[jj];

          UINT nz = nCells_zones[ii];
          UINT ny = inputdata.domain_data.nCells[ 1 ];
          UINT nx = inputdata.domain_data.nCells[ 0 ];

          // extract the zone with the original domain size:

          local_count[0]=nx;
          local_count[1]=ny;
          local_count[2]=nz;

          //local_count correction!
          UINT local_nField=0;
          if((UINT)local_0_start<nz){
            local_offset[2]=local_0_start+zone_offset;
            if((UINT)local_n0+(UINT)local_0_start<=nz){
              local_count[2]=local_n0;
              local_nField=local_n0;
            }else{
              local_count[2]=nz-local_0_start;
              local_nField=nz-local_0_start;
            }
          }else{
            // the sensitivities have dimension of the grid and not the extended grid. So one some processors nothing needes to be read!!!
            local_count[2]=0;
          }

          UINT L,l;
          Vector<REAL> tmp(local_count[0]*local_count[1]*local_count[2]);
          for (UINT iz = 0; iz < (UINT) local_n0; iz++) {
            for (UINT iy = 0; iy < Ny; iy++) {
              for (UINT ix = 0; ix < Nx; ix++) {
                L = General::indexconversion_3d_to_1d( iz, iy, ix, local_n0, Ny, Nx );
                l = General::indexconversion_3d_to_1d( iz, iy, ix, local_nField, ny, nx );
                if( (iz+local_0_start < nz) && (iy < ny) && (ix < nx) ){
                  tmp[ l ] = KField[ L ][0] + beta_mean;
                }
              }
            }
          }

          fftw_free( KField );

          General::log_elapsed_time( watch.elapsed(),
                                     comm,
                                     inputdata.verbosity,
                                     "FFTW",
                                     "generate_3d()" );

          // save the field to disc
          if(CR){
            if(ii>0)
              HDF5Tools::h5_pAppend( tmp,
                                     dir.unconditionalField_h5file,
                                     "/YField",
                                     local_count,
                                     local_offset,
                                     comm );
            else
              HDF5Tools::h5_pWrite( tmp,
                                    dir.unconditionalField_h5file,
                                    "/YField",
                                    inputdata,
                                    inputdata.domain_data.nCells,
                                    local_offset,
                                    local_count,
                                    comm );
          }else{
            if(ii>0)
              HDF5Tools::h5_pAppend( tmp,
                                     dir.kfield_h5file,
                                     "/YField",
                                     local_count,
                                     local_offset,
                                     comm );
            else
              HDF5Tools::h5_pWrite( tmp,
                                    dir.kfield_h5file,
                                    "/YField",
                                    inputdata,
                                    inputdata.domain_data.nCells,
                                    local_offset,
                                    local_count,
                                    comm );
          }

          // set zonation matrix right!
          //if( !General::bFileExists( dir.zonation_matrix[ii] ) ){
          l=tmp.size();
          tmp.resize(0);
          tmp.resize(l,1.0);
          HDF5Tools::h5_pWrite( tmp,
                                dir.zonation_matrix[ii],
                                "/X",
                                inputdata,
                                inputdata.domain_data.nCells,
                                local_offset,
                                local_count,
                                comm );
          //}

        }// END loop over zones

        if(comm_size==1){
          Vector<UINT> local_count,local_offset;
          HDF5Tools::h5g_Read( YFieldVector,
                               dir.kfield_h5file,
                               "/YField",
                               local_offset,
                               local_count,
                               inputdata
                               );

          Dune::Timer watch;

          UINT N = YFieldVector.size();
          //YFieldVector_Well.resize( N );
          YFieldVector_Well = YFieldVector;

          KFieldVector.resize( N );
          for( UINT i=0; i<N; ++i ){
            KFieldVector[i] = std::exp( YFieldVector[i] );
          }
          KFieldVector_Well = KFieldVector;

          General::log_elapsed_time( watch.elapsed(),
                                     inputdata.verbosity,
                                     "EVAL",
                                     "Evaluate KFieldVector once at the beginning." );

        }

      }; // End of void generate_3d()



      // This function can be used for dim=2 only!

      void generate_eigenvalues_2d()
      {

        Dune::Timer watch;
        watch.reset();

        for(UINT ii=0;ii<nzones;ii++){
          //size of the extended domain
          const UINT Ny = nCells_ExtendedDomain[ii][1]; // number of cells in y-direction = #rows
          const UINT Nx = nCells_ExtendedDomain[ii][0]; // number of cells in x-direction = #columns

          // Very important:
          // The indices iy and ix follow this rule:
          // 0 <= iy <= (Ny-1)
          // 0 <= ix <= (Nx-1)
          //
          // Always remember this and you will be on the safe side of life!
          // ix (inner loop) runs faster than iy (outer loop)
          // iy belongs to Ny
          // ix belongs to Nx

          // get the grid size
          const REAL dx = inputdata.domain_data.virtual_gridsizes[0];
          const REAL dy = inputdata.domain_data.virtual_gridsizes[1];


          /*
           * Definition of the block on each processor for the parallel FFT
           *   Idea: Maybe later another data distribution. "fftw_mpi_local_size_many" is a more general routine
           *         BUT: don't forget to load also the Lambdas in the right format
           */
          ptrdiff_t alloc_local, local_n0, local_0_start;

          // get FFTW data
          alloc_local = fftw_mpi_local_size_2d(Ny,Nx, comm,
                                               &local_n0, &local_0_start);

          bool bWisdomFileExists = false;

          Vector<REAL> R_YY(alloc_local);

          Vector<REAL> X(dim, 0.0);
          Vector<REAL> Lambda(inputdata.yfield_properties.zones[ii].correlation_lambda);



          // build up R_YY
          UINT l;
          for( UINT iy = 0; iy < (UINT)local_n0; iy++ ){
            for( UINT ix = 0; ix < Nx; ix++ ) {
              X[0]=std::min( (REAL) ix, (REAL) (Nx-ix) ) * dx;
              X[1]=std::min( (REAL) iy+local_0_start, (REAL) (Ny-(iy+local_0_start)) ) * dy;
              l = General::indexconversion_2d_to_1d( iy, ix, local_n0, Nx );
              R_YY[l]=General::autocovariancefunction(
                                                      inputdata.yfield_properties.zones[ii].variance
                                                      , X
                                                      , Lambda
                                                      , inputdata.yfield_properties.zones[ii].model
                                                      );
            }
          }

          // save the R_YY in HDF5
          Vector<REAL> tmp_vec;

          Vector<UINT> local_count(2,0), local_offset(2,0);
          local_count[0]=nCells_ExtendedDomain[ii][0];
          local_count[1]=local_n0;

          local_offset[1]=local_0_start;

          General::log_elapsed_time( watch.elapsed(),
                                     comm,
                                     inputdata.verbosity,
                                     "EVAL",
                                     "generate_eigenvalues_2d() - part 1" );

          HDF5Tools::h5_pWrite( R_YY
                                , dir.R_YY_h5file[ii]
                                , "/R_YY"
                                , inputdata
                                , nCells_ExtendedDomain[ii]
                                , local_offset
                                , local_count
                                , comm );

          watch.reset();

          // plot this file in octave using
          // load "RYY.dat"REAL
          // pcolor(log(R_YY)); shading flat; colorbar


          //logger << "Use FFTW_FORWARD(2d) to generate eigenvalues..." << std::endl;


          /*
           *
           * calculate the FFT of R_YY in the extended domain!!!
           *
           */
          /*
            char sWisdomFileName[100];
            sprintf( sWisdomFileName, "wisdom_%ux%u.dat", Ny, Nx );

            FILE * rFile;
            rFile = fopen(sWisdomFileName, "r");
            if (rFile != NULL) {
            //logger << "Importing wisdom file "<< sWisdomFileName << std::endl;
            bWisdomFileExists = true;
            fftw_import_wisdom_from_file(rFile);
            fclose(rFile);
            }else{
            logger << "FFTW(2d): wisdom file " << sWisdomFileName << " not found." << std::endl;
            }
          */
          //fftw_complex *R_YY_complex; // Ny x Nx;
          //R_YY_complex = (fftw_complex*) fftw_malloc((/*Ny * Nx*/alloc_local) * sizeof (fftw_complex));

          fftw_complex *Eigenvalues;
          Eigenvalues = (fftw_complex*) fftw_malloc((/*Ny * Nx*/alloc_local) * sizeof (fftw_complex));

          for (UINT iy = 0; iy < /*Ny*/(UINT)local_n0; iy++) {
            for (UINT ix = 0; ix < Nx; ix++) {
              UINT l = General::indexconversion_2d_to_1d( iy, ix, local_n0, Nx );
              Eigenvalues[l][0] = R_YY[ l ];
              Eigenvalues[l][1] = 0.0;
            }
          }

#ifdef FFT_THREADS
          fftw_plan_with_nthreads(FFT_THREADS); // tell the FFT routine how many threads should be used
#endif
          fftw_plan plan_forward_2d;
          if (!bWisdomFileExists) {
            plan_forward_2d = fftw_mpi_plan_dft_2d( Ny, Nx, Eigenvalues, Eigenvalues, comm, FFTW_FORWARD, FFTW_ESTIMATE);
            //logger << "FFTW_FORWARD(2d): create forward plan using FFTW_ESTIMATE" << std::endl;
          } else {
            plan_forward_2d = fftw_mpi_plan_dft_2d( Ny, Nx, Eigenvalues, Eigenvalues, comm, FFTW_FORWARD, FFTW_MEASURE);
            //logger << "FFTW_FORWARD(2d): create forward plan using FFTW_MEASURE" << std::endl;
          }

          fftw_execute(plan_forward_2d);
          fftw_destroy_plan(plan_forward_2d);


          // Count positive eigenvalues.
          int nPositive=0; // Eigenvalue > 1E-6
          int nNegative=0; // Eigenvalue < -1E-6
          int nNearZero=0; // -1E-6 <= Eigenvalue <= 1E6
          tmp_vec.resize(alloc_local);
          for(UINT i=0;i<(UINT)alloc_local;i++){
            tmp_vec[i] = (double)Eigenvalues[i][0];
            if( tmp_vec[i] > 1E-6 )
              nPositive++;
            else if( tmp_vec[i] < -1E-6 )
              nNegative++;
            else
              nNearZero++;
          }

          if( inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
            std::cout << "Properties of Syy: " << std::endl;
            std::cout << alloc_local << " eigenvalues" << std::endl;
            std::cout << nPositive << " eigenvalues > 1E-6" << std::endl;
            std::cout << nNegative << " eigenvalues < -1E-6" << std::endl;
            std::cout << nNearZero << " eigenvalues near zero" << std::endl;
          }

          General::log_elapsed_time( watch.elapsed(),
                                     comm,
                                     inputdata.verbosity,
                                     "FFTW",
                                     "generate_eigenvalues_2d() - part 2" );

          // Save eigenvalues of Syy (= FFT of its first row) in HDF5.
          HDF5Tools::h5_pWrite( tmp_vec
                                , dir.EV_h5file[ii]
                                , "/FFT_R_YY"
                                , inputdata
                                , nCells_ExtendedDomain[ii]
                                , local_offset
                                , local_count
                                , comm );


          fftw_free(Eigenvalues);
          //fftw_free(R_YY_complex);
        } //end for-loop of zones


        // save the Y-Field properties in text format
        if(my_rank==0)
          General::save_kfield_properties(
                                          dir.kfield_properties_file
                                          , dim
                                          , inputdata.domain_data.extensions
                                          , inputdata.domain_data.nCells
                                          , inputdata.yfield_properties);

      }; // End of void generate_eigenvalues_2d()



      void generate2d( const bool bNewEigenvalues,
                       const bool bNewField,
                       const bool CR=false
                       )
      {
        if( bNewEigenvalues
            ||
            ( !bNewEigenvalues
              &&
              !General::load_kfield_properties(
                                               dir.kfield_properties_file
                                               , dim
                                               , inputdata.domain_data.extensions
                                               , inputdata.domain_data.nCells
                                               , inputdata.yfield_properties
                                               , my_rank
                                               )
              )
            ) {
          // generate new eigenvalues
          generate_eigenvalues_2d();
        }
        else {
          logger << "Found suitable Y-field characteristics!" << std::endl;
        }

        if( !bNewField ){
          if( my_rank == 0 ){
            std::cout << "Warning: New Y-field will be generated due changed Y-field characteristics." << std::endl;
          }
        }

        for(UINT ii=0;ii<nzones;ii++){

          // size of the extended domain
          const UINT Ny = nCells_ExtendedDomain[ii][1];
          const UINT Nx = nCells_ExtendedDomain[ii][0];

          //logger<<"Ny : "<<Ny<<std::endl;
          //logger<<"Nx : "<<Nx<<std::endl;

          REAL scalingfactor = REAL( Ny * Nx );

          bool bWisdomFileExists = false;

          /*
           * Definition of the block on each processor for the parallel FFT
           *   Idea: Maybe later another data distribution. "fftw_mpi_local_size_many" is a more general routine
           *         BUT: don't forget to load also the Lambdas in the right format
           */
          ptrdiff_t alloc_local, local_n0, local_0_start;

          // get the FFT data
          alloc_local = fftw_mpi_local_size_2d(Ny,Nx, comm,
                                               &local_n0, &local_0_start);


          // load the eigenvalues from HDF5 file
          Vector<REAL> tmp_vec;

          Vector<UINT> local_count(2,0), local_offset(2,0);
          local_count[0]=nCells_ExtendedDomain[ii][0];
          local_count[1]=local_n0;

          local_offset[1]=local_0_start;

          HDF5Tools::h5_pRead( tmp_vec
                               , dir.EV_h5file[ii]
                               , "/FFT_R_YY"
                               , inputdata
                               , local_offset
                               , local_count
                               , comm
                               );

          Dune::Timer watch;
          watch.reset();

          // initialize pseudo-random generator
          std::random_device rd;
          std::mt19937_64 gen; // 64-bit Mersenne Twister
          std::normal_distribution<> ndist(0,1); // mean = 0, sigma = 1

          int seed = inputdata.yfield_properties.random_seed;
          if( seed == 0 ){
            seed = rd();
          }
          // different seed for different processors -> very IMPORTANT to obtain the right result!
          seed += my_rank +ii;

          gen.seed( seed );
            
          if(my_rank==0)
            std::cout << "seed = " << seed << std::endl;


          fftw_complex *KField;
          KField = (fftw_complex*) fftw_malloc(( /*Ny * Nx*/alloc_local ) * sizeof (fftw_complex));

          fftw_complex random_epsilon;


          REAL lambda=0.0;
          // add some rendomness
          //  important to loop over tmp_vec.size() not over all alloc_local, because tmp_vec.size() <= alloc_local

          int countNegativeValues = 0;
          REAL lambda_min = 1e+12;

          for (UINT jj = 0; jj < tmp_vec.size(); jj++) {

            random_epsilon[0] = ndist(gen);
            random_epsilon[1] = ndist(gen);
      
            lambda = tmp_vec[ jj ] / scalingfactor;

            lambda_min = std::min( lambda_min, lambda );

            if (lambda > 1e-12) {

              KField[ jj ][0] = sqrt(lambda) * random_epsilon[0];
              KField[ jj ][1] = sqrt(lambda) * random_epsilon[1];

            }
            else {

              countNegativeValues++;
              //std::cout << "2D-case ---------> WARNING: negative eigenvalue lambda = "
              //          << lambda
              //          << " Choose 0."
              //          << std::endl;
              KField[ jj ][0] = 0.0; // sqrt(-lambda) * random_epsilon[0];
              KField[ jj ][1] = 0.0; // sqrt(-lambda) * random_epsilon[1];

            }
          }


          if( inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
            std::cout << "2D-case ---------> countNegativeValues = "
                      << countNegativeValues
                      << std::endl;
            std::cout << "2D-case ---------> lambda_min = "
                      << lambda_min
                      << std::endl;
          }


          tmp_vec.resize(0);

#ifdef FFT_THREADS
          fftw_plan_with_nthreads(FFT_THREADS); // tell the FFT routine how many threads should be used
#endif
          fftw_plan plan_backward_2d;

          if (!bWisdomFileExists) {
            //logger << "FFTW(2d): create backward plan using FFTW_ESTIMATE" << std::endl;
            plan_backward_2d = fftw_mpi_plan_dft_2d( Ny, Nx, KField, KField, comm,
                                                     FFTW_FORWARD, // FFTW_BACKWARD,
                                                     FFTW_ESTIMATE );
          } else {
            //logger << "FFTW(2d): create backward plan using FFTW_MEASURE" << std::endl;
            plan_backward_2d = fftw_mpi_plan_dft_2d( Ny, Nx, KField, KField, comm,
                                                     FFTW_FORWARD, // FFTW_BACKWARD,
                                                     FFTW_MEASURE );
          }

          //fftw_execute_dft(plan_backward_2d, EigenvaluesRandomized, KField);
          fftw_execute(plan_backward_2d);
          fftw_destroy_plan(plan_backward_2d);


          /*
           * VERY IMPORTANT: beta_mean needs to be the same on all Prozesses. So a MPI_Bcast is needed!
           */
          REAL beta_mean=0.0;
          if(!CR){
            if(my_rank==0){
              if(inputdata.yfield_properties.zones[ii].qbb_y>0.0)
                {
                  REAL rv = ndist(gen);
                  beta_mean = inputdata.yfield_properties.zones[ii].beta + rv * std::sqrt( inputdata.yfield_properties.zones[ii].qbb_y );
                }
              else
                beta_mean = inputdata.yfield_properties.zones[ii].beta;
            }
            MPI_Bcast(&(beta_mean),1,MPI_DOUBLE,0,comm);
          }
          logger << "FFTFieldGenerator(2d): beta_mean = " << beta_mean << std::endl;


          //correction due to previous zones
          UINT zone_offset=0;
          for(UINT jj=0; jj<ii; jj++)
            zone_offset+=nCells_zones[jj];

          UINT nx = inputdata.domain_data.nCells[0];

          UINT ny = nCells_zones[ii];

          // extract the zone with the original domain size:

          local_count[0]=nx;
          local_count[1]=ny;

          //local_count correction!
          UINT local_nField=0;
          if((UINT)local_0_start<ny){
            local_offset[1]=local_0_start+zone_offset;
            if((UINT)local_n0+(UINT)local_0_start<=ny){
              local_count[1]=local_n0;
              local_nField=local_n0;
            }else{
              local_count[1]=ny-local_0_start;
              local_nField=ny-local_0_start;
            }
          }else{
            // the sensitivities have dimension of the grid and not the extended grid. So one some processors nothing needes to be read!!!
            local_count[1]=0;
          }


          Vector<REAL> tmp(local_count[0]*(local_count[1]));
          Vector<REAL> extended_tmp(Nx*(UINT)local_n0);
          UINT L,l;
          //logger<<" l for zone #"<<ii<<", local_n0 = "<<local_n0<<std::endl;
          for( UINT iy=0; iy < (UINT)local_n0; iy++ ){
            for( UINT ix = 0; ix < Nx; ix++ ){
              L = General::indexconversion_2d_to_1d( iy, ix, local_n0, Nx );
              l = General::indexconversion_2d_to_1d( iy, ix, local_nField , nx );
              if( (iy+local_0_start < ny) && (ix < nx) ){
                //if(zone_offset)
                //    logger<<"l : "<<l<<beta_mean<<std::endl;

                tmp[ l ] = KField[ L ][0] + beta_mean;
              }

              extended_tmp[L] = KField[ L ][0] + beta_mean;

            }
          }

          fftw_free( KField );

          General::log_elapsed_time( watch.elapsed(),
                                     comm,
                                     inputdata.verbosity,
                                     "FFTW",
                                     "generate_2d()" );

          /*
            logger<<"Zone #"<<ii<<":"<<std::endl;
            logger<<"nCells :"<<inputdata.domain_data.nCells[0]<<", "<<inputdata.domain_data.nCells[1]<<std::endl;
            logger<<"KFieldVector[ii].size() :"<<KFieldVector[ii].size()<<std::endl;
            logger<<"local_count :"<<local_count[0]<<", "<<local_count[1]<<std::endl;
            logger<<"local_offset :"<<local_offset[0]<<", "<<local_offset[1]<<std::endl;
            logger<<"kfield_filename : "<<kfield_filename<<std::endl;
            logger<<"zone_offset : "<<zone_offset<<std::endl;
          */
          // save the field to disc
          if(CR){
            if(ii>0)
              HDF5Tools::h5_pAppend( tmp,
                                     dir.unconditionalField_h5file,
                                     "/YField",
                                     local_count,
                                     local_offset,
                                     comm );
            else
              HDF5Tools::h5_pWrite( tmp,
                                    dir.unconditionalField_h5file,
                                    "/YField",
                                    inputdata,
                                    inputdata.domain_data.nCells,
                                    local_offset,
                                    local_count,
                                    comm );
          }else{
            if(ii>0)
              HDF5Tools::h5_pAppend( tmp,
                                     dir.kfield_h5file,
                                     "/YField",
                                     local_count,
                                     local_offset,
                                     comm );
            else {
              HDF5Tools::h5_pWrite( tmp,
                                    dir.kfield_h5file,
                                    "/YField",
                                    inputdata,
                                    inputdata.domain_data.nCells,
                                    local_offset,
                                    local_count,
                                    comm );

              Vector<UINT> extended_local_count(2,0), extended_local_offset(2,0);
              extended_local_count[0] = nCells_ExtendedDomain[ii][0];
              extended_local_count[1] = local_n0;
              extended_local_offset[1] = local_0_start;
              HDF5Tools::h5_pWrite( extended_tmp,
                                    dir.extended_yfield_h5file,
                                    "/YField",
                                    inputdata,
                                    nCells_ExtendedDomain[ii],
                                    extended_local_offset,
                                    extended_local_count,
                                    comm );
            }
          }

          // set zonation matrix right!
          // if( !General::bFileExists( dir.zonation_matrix[ii] ) ) {
          l=tmp.size();
          tmp.resize(0);
          tmp.resize(l,1.0);
          HDF5Tools::h5_pWrite( tmp
                                , dir.zonation_matrix[ii]
                                , "/X"
                                , inputdata
                                , inputdata.domain_data.nCells
                                , local_offset
                                , local_count
                                ,  comm );
          //}

        }// END for-loop over zones!

        if(comm_size==1){
          Vector<UINT> local_count,local_offset;
          HDF5Tools::h5g_Read( YFieldVector,
                               dir.kfield_h5file,
                               "/YField",
                               local_offset,
                               local_count,
                               inputdata
                               );
          UINT N = YFieldVector.size();
          //YFieldVector_Well.resize( N );
          YFieldVector_Well = YFieldVector;

          KFieldVector.resize( N );
          for( UINT i=0; i<N; ++i ){
            KFieldVector[i] = std::exp( YFieldVector[i] );
          }
          KFieldVector_Well = KFieldVector;
        }

      }; // End of void generate2d()



      bool load_from_file()
      {
        if( !General::load_kfield_properties(
                                             dir.kfield_properties_file
                                             , dim
                                             , inputdata.domain_data.extensions
                                             , inputdata.domain_data.nCells
                                             , inputdata.yfield_properties
                                             , my_rank
                                             )
            )
          {
            return false;
          }
        else
          {

            if( !General::bFileExists( dir.kfield_h5file ) )
              {
                std::cout << "Warning: File is missing: "
                          <<  dir.kfield_h5file
                          << std::endl;
                return false;
              }

            logger << "Found suitable Y-Field characteristics! Load the kfield from "
                   <<  dir.kfield_h5file
                   << std::endl;

            if(comm_size==1){
              Vector<UINT> local_count,local_offset;
              HDF5Tools::h5g_Read( YFieldVector,
                                   dir.kfield_h5file,
                                   "/YField",
                                   local_offset,
                                   local_count,
                                   inputdata
                                   );
              UINT N = YFieldVector.size();
              //YFieldVector_Well.resize( N );
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



      template<typename GV_GW>
      void setWellConductivities( const GV_GW& gv_gw ){

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
                                   inputdata.verbosity,
                                   "EVAL",
                                   "setWellConductivities()" );

      }



      void import_from_vector(const Vector<REAL>& yfield_vector) {

        Dune::Timer watch;

        UINT VectorSize = 0;
        if (dim == 3) {
          const UINT L = inputdata.domain_data.nCells[2];
          const UINT M = inputdata.domain_data.nCells[1];
          const UINT N = inputdata.domain_data.nCells[0];
          VectorSize = L * M * N;
        }
        else {
          const UINT M = inputdata.domain_data.nCells[1];
          const UINT N = inputdata.domain_data.nCells[0];
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
                                   inputdata.verbosity,
                                   "EVAL",
                                   "import_from_vector: Evaluate KFieldVector once at the beginning." );

      };



      void parallel_import_from_local_vector( const Vector<REAL>& local_Yfield_vector,
                                              const Vector<UINT>& local_offset,
                                              const Vector<UINT>& local_count )
      {
        Dune::Timer watch;

        UINT Nlocal = local_Yfield_vector.size();
        logger << "Nlocal = " << Nlocal << std::endl;

        // TODO: Maybe we do not need this LocalYFieldVector anymore!
        LocalYFieldVector.clear();
        LocalYFieldVector_Well.clear();

        LocalYFieldVector = local_Yfield_vector; // this vector contains the Y-field data from the current hyperslab
        LocalYFieldVector_Well = local_Yfield_vector;
        logger << "LocalYFieldVector.size() = " 
               << LocalYFieldVector.size() << std::endl;
        logger << "LocalYFieldVector_Well.size() = " 
               << LocalYFieldVector_Well.size() << std::endl;

        LocalKFieldVector.clear();
        LocalKFieldVector_Well.clear();

        LocalKFieldVector.resize( Nlocal );
        for( UINT i=0; i<Nlocal; ++i )
          LocalKFieldVector[i] = std::exp( local_Yfield_vector[i] );

        LocalKFieldVector_Well = LocalKFieldVector;
        logger << "LocalKFieldVector_Well.size() = " 
               << LocalKFieldVector_Well.size() << std::endl;

        LocalCount = local_count; // this vector contains the dimensions (the number of cells per dimension) of the current hyperslab
        LocalOffset = local_offset; // this vector describes the distance of the current hyperslab from the origin

        logger << "LocalCount = " << LocalCount << std::endl;
        logger << "LocalOffset = " << LocalOffset << std::endl;

        /*
        General::log_elapsed_time( watch.elapsed(),
                                   comm,
                                   inputdata.verbosity,
                                   "EVAL",
                                   "parallel_import_from_local_vector: Evaluate KFieldVector once at the beginning." );
        */
      };


      // not used now, just test-wise:
      // overloaded with extra parameter for the well:
      /*
        void parallel_import_from_local_vector(
        const Vector<REAL>& local_Yfield_vector,
        const Vector<REAL>& local_Yfield_Well_vector,
        const Vector<UINT>& local_count,
        const Vector<UINT>& local_offset
        )
        {

        LocalYFieldVector.clear();

        LocalYFieldVector.resize(local_Yfield_vector.size());
        LocalYFieldVector = local_Yfield_vector; // this vector contains the Y-field data from the current hyperslab

        LocalYFieldVector_Well.resize(local_Yfield_Well_vector.size());
        LocalYFieldVector_Well = local_Yfield_Well_vector;

        LocalCount = local_count; // this vector contains the dimensions (the number of cells per dimension) of the current hyperslab
        LocalOffset = local_offset; // this vector describes the distance of the current hyperslab from the origin
        };

      */



      // not used now, just test-wise:
      template<typename GV>
      void parallel_add_wells_to_hdf5( const GV& gv,
                                       const std::string& data_name,
                                       const Vector<UINT> &nCells,   // nKnotsPerDim
                                       const Vector<CTYPE> &gridsizes,
                                       const std::string& data_filename ) {

        const int blocksizePerDimension = 1;
        HDF5Tools::h5g_pWrite( gv,
                               LocalYFieldVector_Well,
                               data_filename,
                               data_name,
                               inputdata,
                               nCells,
                               blocksizePerDimension
                               , FEMType::DG  // P0
                               );

      }


      void export_to_vector(Vector<REAL>& kfield_vector )
      {
        /*
         * ->deleted!
         if( inputdata.yfield_properties.variance == 0 )
         {
         kfield_vector[ 0 ] = KFieldVector[ 0 ];
         return;
         }
        */
        UINT VectorSize = 0;
        if (dim == 3)
          {
            const UINT L = inputdata.domain_data.nCells[2];
            const UINT M = inputdata.domain_data.nCells[1];
            const UINT N = inputdata.domain_data.nCells[0];
            VectorSize = L * M * N;
          }
        else
          {
            const UINT M = inputdata.domain_data.nCells[1];
            const UINT N = inputdata.domain_data.nCells[0];
            VectorSize = M * N;
          }

        kfield_vector.resize( VectorSize );

        for (UINT l = 0; l < VectorSize; l++)
          kfield_vector[l] = KFieldVector[l];

        return;
      };



      // Given a location in space, return the appropriate
      // index of the parameter field.
      std::size_t getFieldIndex( const Dune::FieldVector<REAL,dim>& xglobal ) const {

        std::size_t myIndex = 0;
        Vector<UINT> global_index;
        global_index.resize(dim);

        /*
          Translate the global coordinate 'x' into its corresponding Yfield-tensor-index 'global_index'.
          This requires only the virtual gridsize of the virtual Yfield grid.
          So 'global_index' is just a virtual index.
          In the same way, 'local_index' is a virtual index that is used to access virtual grid cells values
          stored in 'KFieldVector'.
        */

        for (UINT i=0; i<dim; i++) {
          REAL fraction = xglobal[i] / inputdata.domain_data.virtual_gridsizes[i];
          global_index[i] = static_cast<UINT>( fraction );
          if( std::abs( global_index[i] - fraction ) < GEO_EPSILON
              &&
              global_index[i] > 0 )
            global_index[i] -= 1;
        }

        if( comm_size == 1 ) {
          // single-processore case:
          myIndex = General::indexconversion_to_1d( global_index,
                                                    inputdata.domain_data.nCells );
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


      void evaluate2d( const Dune::FieldVector<REAL,2>& xglobal,
                       REAL& output,
                       REAL& output_Well,
                       bool lnK=false ) const
      {

        // Special Test Case 2d:
        //output = -6.0 + x[0]/500.0;
        //output = -6.0 + 3.0 * ((x[1]-250.0)+(x[0]-250.0))/500.0;

#ifdef HOMOGEN
        {
          output = -6.0;
          return;
        }
#endif

#ifdef CENTER_SQUARE
        {
          // Center Square:
          REAL a = std::abs( x[0]-5.0 );
          REAL b = std::abs( x[1]-5.0 );
          //if( a*a + b*b < 2500.0 )
          if( std::max(a,b) < 3.0 )
            output = -6.0;
          else
            output = -2.0;
          return;
        }
#endif

#ifdef CHESS_PATTERN
        {
          // chess board pattern:
          int mx = (static_cast<int>( 5.0 * x[0]/500.0 ))%2;
          int my = (static_cast<int>( 5.0 * x[1]/500.0 ))%2;
          if ( my==0 && mx==0 ) output = -2.0; // 1.0 * 1E+9;
          if ( my==0 && mx==1 ) output = -4.0; // 1.0 * 1E-6;
          if ( my==1 && mx==0 ) output = -4.0; // 1.0 * 1E-2;
          if ( my==1 && mx==1 ) output = -2.0; // 1.0 * 1E+5;
          return;
        }
#endif

        output_Well = output;

        UINT l = getFieldIndex( xglobal );


        if( comm_size == 1 ) // single-processore case:
          {
            if( lnK )
              output_Well = YFieldVector_Well[ l ];
            else
              output_Well = KFieldVector_Well[ l ];

            if( inputdata.yfield_properties.well_type == "cell" ){
              output = output_Well;
            }
            else{
              if( lnK )
                output = YFieldVector[ l ];
              else
                output = KFieldVector[ l ];
            }

          }

        else // multi-processore case:

          {


            if( lnK )
              output_Well = LocalYFieldVector_Well[l];
            else
              output_Well = LocalKFieldVector_Well[l];


            if( inputdata.yfield_properties.well_type == "cell" ){
              output = output_Well;
            }
            else{
              if( lnK )
                output = LocalYFieldVector[l];
              else
                output = LocalKFieldVector[l];
            }
          }

      }


      void evaluate3d( Dune::FieldVector<REAL,3>& xglobal,
                       REAL& output,
                       REAL& output_Well,
                       bool lnK=false ) const
      {
#ifdef HOMOGEN
        {
          output = -6.0;
          return;
        }
#endif

        /*
          Translate the global coordinate 'x' into its corresponding Yfield-tensor-index 'global_index'.
          This requires only the virtual gridsize of the virtual Yfield grid.
          So 'global_index' is just a virtual index.
          In the same way, 'local_index' is a virtual index that is used to access virtual grid cells values
          stored in 'KFieldVector'.
        */
        UINT l = getFieldIndex( xglobal );

        if( comm_size == 1 ) {

          if( lnK )
            output_Well = YFieldVector_Well[ l ];
          else
            output_Well = KFieldVector_Well[ l ];

          if( inputdata.yfield_properties.well_type == "cell" ){
            output = output_Well;
          }
          else {
            if( lnK )
              output = YFieldVector[l];
            else
              output = KFieldVector[l];
          }

        }
        else {

          if( lnK )
            output_Well = LocalYFieldVector_Well[l];
          else
            output_Well = LocalKFieldVector_Well[l];

          if( inputdata.yfield_properties.well_type == "cell" ){
            output = output_Well;
          }
          else {
            if( lnK )
              output = LocalYFieldVector[l];
            else
              output = LocalKFieldVector[l];
          }
        }

      }





      void init()
      {
        logger << "------------------------------------\ngenerate_Yfield()" << std::endl;

        //const int dim = inputdata.dim;

        // Be aware of this:
        // The virtual YField grid is a structured grid (FFT requires this!)
        // no matter what kind of grid we are using
        // in DUNE for the solution of the PDEs.

        UINT nAllCells = inputdata.domain_data.nCells[0] * inputdata.domain_data.nCells[1];
        if (dim == 3)
          nAllCells *= inputdata.domain_data.nCells[2];

        this->extend_domain();
        //this->export_nCellsExt( nCellsExt );

        Dune::Timer watch,watch2;
        int generate_new=0;
        // do something on P0 (sequential) only

        //if( helper.rank() == 0 )
        if( my_rank == 0 )
          {
            /*
              if( inputdata.yfield_properties.variance != 0 )
              {
            */
            watch.reset();

            //loading is in HDF5 or *.dat, (see FFTFieldGenerator.hh)
            if( !inputdata.problem_types.new_YField
                &&
                !inputdata.problem_types.new_Eigenvalues
                &&
                this->load_from_file( ) ) // loads data only if sequential program, in parallel only consistency check
              {
                logger << "Reading existing Y-field took: "
                       << watch.elapsed() << " sec."
                       << std::endl;
                std::cout << "Reading existing Y-field took: "
                          << watch.elapsed() << " sec."
                          << std::endl;
              }
            else
              {

                logger << "Generating new data ..." << std::endl;
                std::cout << "Generating new data ..."
                          << std::endl;

                generate_new = 1;

              }

          } // endif helper.rank() == 0

        MPI_Bcast(&(generate_new), 1, MPI_INT, 0, comm );


        if(generate_new){

          if( dim == 3 )
            this->generate3d( inputdata.problem_types.new_Eigenvalues,
                              inputdata.problem_types.new_YField
                              );
          else
            this->generate2d( inputdata.problem_types.new_Eigenvalues,
                              inputdata.problem_types.new_YField
                              );

          logger << "Generating new Y-field took: "
                 << watch.elapsed() << " sec. (this duration includes generating/reading the eigenvalues)"
                 << std::endl;

          //if( helper.rank() == 0 )
          if( my_rank == 0 )
            std::cout
              << "Generating new Y-field took:  "
              << watch.elapsed() << " sec. (this duration includes generating/reading the eigenvalues)"
              << std::endl;

        }

      }; // end of void init()





      template<typename GV>
      void export_field_to_vector_on_grid( const GV& gv,
                                           Vector<REAL>& log_conductivity_field,
                                           Vector<REAL>& well_conductivity_field,
                                           const int gv_level=0
                                           ) const {

        Dune::Timer watch;

        // Get the conductivity field on the grid for the vtkwriter
        typedef typename GV::Grid GRIDTYPE; // get the grid-type out of the gridview-type
        //typedef typename GV::Grid::template Codim<0>::Entity Entity;
        //typedef typename GV::Grid::template Codim<0>::EntityPointer EntityPointer;

        //typedef typename GV::Traits::template Codim < 0 > ::Iterator ElementLeafIterator;
        typedef typename GV::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementLeafIterator;

        //typedef typename ElementLeafIterator::Entity::Geometry LeafGeometry;

        // make a mapper for codim 0 entities in the leaf grid
        //Dune::LeafMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout>
        //  mapper(gv.grid()); // get the underlying hierarchical grid out ot the gridview

        // make a mapper for codim 0 entities in the level(0) grid
        Dune::LevelMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout>
          mapper(gv.grid(),gv_level); // get the underlying hierarchical grid out ot the gridview


        const int nGridCells = mapper.size();

#ifdef DEBUG_LOG
        logger << "DEBUG: mapper.size() = " << nGridCells << std::endl;
#endif

        log_conductivity_field.resize( nGridCells );
        well_conductivity_field.resize( nGridCells );

        if( my_rank == 0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL )
          std::cout << "Plot Yfield: Loop over leaf elements..." << std::endl;

        const typename GV::IndexSet& indexset = gv.indexSet();

        for (ElementLeafIterator it = gv.template begin<0,Dune::All_Partition> ()
               ; it != gv.template end<0,Dune::All_Partition> (); ++it) {
          int index = indexset.index(*it);

          Dune::FieldVector<REAL,dim> xglobal = it->geometry().center();
          REAL logK = 0.0;
          REAL logK_Well = 0.0;

          bool lnK = true;
#ifdef DIMENSION3
          this->evaluate3d(xglobal,logK,logK_Well,lnK);
#else
          this->evaluate2d(xglobal,logK,logK_Well,lnK);
#endif
          log_conductivity_field[index] = logK;
          well_conductivity_field[index] = logK_Well;
        }


        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   inputdata.verbosity,
                                   "EVAL",
                                   "export_field_to_vector_on_grid" );

      } // end of export_field_to_vector_on_grid()


      template<typename GV>
      void plot2vtu( const GV& gv,
                     const std::string filename,
                     const std::string title,
                     const int gv_level=0
                     ) const {

        Vector<REAL> log_conductivity_field;
        Vector<REAL> well_conductivity_field;

        export_field_to_vector_on_grid( gv,
                                        log_conductivity_field,
                                        well_conductivity_field,
                                        gv_level
                                        );

        VTKPlot::output_vector_to_vtu( gv,
                                       log_conductivity_field,
                                       filename,
                                       title,
                                       inputdata.verbosity,
                                       true,
                                       0
                                       );

        if( inputdata.plot_options.vtk_plot_wells )
          VTKPlot::output_vector_to_vtu( gv,
                                         well_conductivity_field,
                                         filename + "_Wells",
                                         title + "_Wells",
                                         inputdata.verbosity,
                                         true,
                                         0
                                         );

      }


      std::string getFilename() const {
        return dir.kfield_h5file;
      }


    }; // class FFTFieldGenerator

  } // namespace Gesis
} // namespace Dune



#endif	/* DUNE_GESIS_FFT_FIELD_GENERATOR_HH */
