/*
 * File    : MPI_Tools.hh 
 * Authors : Ronnie L. Schwede & Adrian Ngo (2009-2013)
 *
 * 2014/07/29: add time logger
 *
 *
 */
#ifndef DUNE_GEOINVERSION_MPI_ROUTINES_HH
#define DUNE_GEOINVERSION_MPI_ROUTINES_HH

namespace Dune {

  namespace GeoInversion {

    class MPI_Tools {

    private:
      MPI_Tools(){};
      
    public:

      template<typename MeasurementElement>
      static void redist_measurements( std::vector<MeasurementElement>& mp )
      {
        Dune::Timer watch;
        
        int mpi_size;
        MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );

        if( mpi_size < 2 )
          return; // Do nothing, otherwise redistribute to ALL.

        int mpi_rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
        //logger << "redist_measurements: mpi_size = " << mpi_size
        //       << " mpi_rank = " << mpi_rank
        //       << std::endl;

        std::vector<REAL> out(mp.size(),0.0);
        std::vector<REAL> in(mp.size(),0.0);

        for (UINT ii=0;ii<mp.size();ii++){
          out[ii] = mp[ii].value;
        }

        //for (UINT i=0;i<out.size();i++){
        //  logger << "DEBUG: before send: out[" << i << "] = " << out[i] << std::endl;
        //}

        // Note:
        // Negative value must be redistributed as well. Therefore, we must use MPI_SUM instead of MPI_MAX.
        // And: A measuring value must not belong to more than one processor.
        // And: mp[*].value must be initialized with 0.
        MPI_Allreduce( &(out[0]),
                       &(in[0]),
                       mp.size(),
                       MPI_DOUBLE,
                       MPI_SUM,
                       MPI_COMM_WORLD );

        //for (UINT i=0;i<in.size();i++){
        //  logger << "DEBUG: after send: in[" << i << "] = " << in[i] << std::endl;
        //}

        for (UINT ii=0;ii<mp.size();ii++){
          mp[ii].value = in[ii];
        }


        General::log_elapsed_time( watch.elapsed(),
                                   MPI_COMM_WORLD,
                                   General::verbosity,
                                   "REDIST",
                                   "redist_measurements" );

      } // redist_measurements






      template<typename VEC>
      static void redist_well_data( VEC& mp, const int dim )
      {

        Dune::Timer watch;

        int mpi_size;
        int mpi_rank;
        MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
        MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );


        if( mpi_size < 2 )
          return; // Do nothing, otherwise redistribute well locations to ALL.


        //for(UINT ii=0;ii<mp.size();ii++){
        //  logger << "DEBUG: well = " << mp[ii].second << " k=" << mp[ii].first << std::endl;
        //}

        int datasize = (dim+1)*mp.size();

        MPI_Allreduce( &datasize,
                       &datasize,
                       1,
                       MPI_INT,
                       MPI_MAX,
                       MPI_COMM_WORLD );
        

        std::vector<REAL> out( datasize );
        std::vector<REAL> in( datasize );
  
        //logger << "MPI_Tools: DEBUG: mp.size() = " << mp.size() << std::endl;
        //logger << "MPI_Tools: DEBUG: datasize = " << datasize << std::endl;

        for( UINT ii=0; ii<mp.size(); ii++ ){
          out[ii] = mp[ii].first;
        }

        for( int jj=0; jj<dim; jj++ ){
          for( UINT ii=0; ii<mp.size(); ii++ ){
            out[ii+mp.size()*(jj+1)] = mp[ii].second[jj];
          }
        }


        //logger << "MPI_Tools: DEBUG: out.size() = " << out.size() << std::endl;


        //for(UINT i=0;i<datasize;i++){
        //  logger << "DEBUG: out[" << i << "] = " << out[i] << std::endl;
        //}

        MPI_Allreduce( &(out[0]),
                       &(in[0]),
                       datasize,
                       MPI_DOUBLE,
                       MPI_MAX,
                       MPI_COMM_WORLD );

        //for(UINT i=0;i<datasize;i++){
        //  logger << "DEBUG: in[" << i << "] = " << in[i] << std::endl;
        //}

        int mpsize = datasize/(dim+1);
        mp.resize( mpsize );
          
        for(int ii=0;ii<mpsize;ii++){
          mp[ii].first = in[ii];
        }

        for (int jj=0;jj<dim;jj++){
          for(int ii=0;ii<mpsize;ii++){
            mp[ii].second[jj] = in[ii+mpsize*(jj+1)];
          }
        }

        //for(UINT ii=0;ii<mpsize;ii++){
        //  logger << "DEBUG: well = " << mp[ii].second << " k=" << mp[ii].first << std::endl;
        //}

        General::log_elapsed_time( watch.elapsed(),
                                   MPI_COMM_WORLD,
                                   General::verbosity,
                                   "REDIST",
                                   "redist_well_data" );

      } // redist_well_data


    }; // class MPI_Tools

  } // GeoInversion

} // Dune

#endif //DUNE_GEOINVERSION_MPI_ROUTINES_HH






