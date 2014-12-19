/*
 * File:   HDF5Tools.hh
 * Author: Adrian Ngo / Ronnie Schwede, 2010-2014
 *
 *  This requires the MPI version of HDF5.
 *
 *  Overview of HDF5 methods:
 *  h5_Write   : sequential write of multi-dim array to HDF5 file
 *  h5_Read    : sequential read of multi-dim array from HDF5 file
 *  h5_pWrite  : parallel write of multi-dim array to HDF5 file
 *  h5_pRead   : parallel read of multi-dim array from HDF5 file
 *  h5g_Write  : sequential write of grid data to HDF5 file
 *  h5g_Read   : sequential read of grid data from HDF5 file
 *  h5g_pWrite : parallel write of grid data to HDF5 file
 *  h5g_pRead  : parallel read of grid data from HDF5 file
 *  h5_pAppend : parallel append of multi-dim array to exisiting HDF5 file
 *  h5_pWrite1d: parallel write of local 1d arrays to a global 1d array on one continuous disjointly hyperslabbed HDF5 file
 *  h5_pRead1d : parallel read of local 1d arrays from a global 1d array on one continuous disjointly hyperslabbed HDF5 file
 *
 */

#ifndef DUNE_GESIS_HDF5_TOOLS_HH
#define	DUNE_GESIS_HDF5_TOOLS_HH

#include <dune/pdelab/common/geometrywrapper.hh>
#include <assert.h>
#include <sstream>
#include <dune/gesis/common/MyMPIComm.hh>

/***
    http://www.hdfgroup.org/HDF5/release/obtain5.html
    HDF5 (Hierarchical Data Format 5) Software Library and Utilities
    Copyright 2006-2014 by The HDF Group.
 ***/
#if HAVE_HDF5
#include "hdf5.h"
#endif


extern CLogfile logger;


namespace Dune {
  namespace Gesis {


    class HDF5Tools{

    private:
      HDF5Tools(){};

    public:

      /*
       * Note that in hdf5, all array indices are ordered the other way round!
       *
       * So far, we had this:
       *
       *   nCells[0] = number of cells in x-direction
       *   nCells[1] = number of cells in y-direction
       *   nCells[2] = number of cells in z-direction
       *
       * But storing data in hdf5 requires to define a vector 'dims' (in 3d) with
       *
       *  dims[2] = (hsize_t) ( nCells[0] ); // x-direction
       *  dims[1] = (hsize_t) ( nCells[1] ); // y-direction
       *  dims[0] = (hsize_t) ( nCells[2] ); // z-direction
       *
       * and in 2d:
       *
       * dims[1] = (hsize_t) ( nCells[0] ); // x-direction
       * dims[0] = (hsize_t) ( nCells[1] ); // y-direction
       *
       */






      /** Function to create hdf5 file. It stores the data of the virtual Y-field vector in spatial format such that hyperslabs of data can be cut out later ( see function h5g_pRead() )
       *
       \tparam data = vector that is been stored. (one processor only!!!)
       \tparam data_name = name/location where the data a stored in the hdf5-file
       \tparam nCells = numbers of virtual cells in each dimension
       \tparam nChunks = numbers of chunks (needed for hdf5) in each dimension
       \tparam gridsizes = virtual grid-sizes, needed only for extra data
       \tparam filename = name of the hdf5 file where all the data will be stored
       *
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */
      template<typename IDT>
      static void h5g_Write(
                            const Vector<REAL> & data
                            , const std::string & filename
                            , const std::string & groupname
                            , const IDT & inputdata
                            )
      {
        Dune::Timer watch;
        logger << "h5g_Write: " << filename << std::endl;

        UINT dim = inputdata.domain_data.nCells.size();

        //std::cout << std::endl << "h5g_Write: dim = " << dim << std::endl;
        logger << "h5g_Write: groupname = " << groupname << std::endl;

        /* Create a new file using default properties. */
        hid_t file_id = H5Fcreate(
                                  filename.c_str()    // name of the file to be created
                                  , H5F_ACC_TRUNC          // if you are trying to create a file that exists already, the existing file will be truncated, i.e., all data stored on the original file will be erased
                                  , H5P_DEFAULT            // default file creation property list
                                  , H5P_DEFAULT            // default file access property list
                                  );

        assert( file_id > -1 );
        //logger << "h5g_Write: h5-file created: " << filename.c_str() << std::endl;


        hsize_t mdims[1];
        mdims[0] = data.size();
        hid_t memspace_id = H5Screate_simple(           // H5Screate_simple creates a new simple dataspace and opens it for access, returning a dataset identifier.
                                             1          // rank=1 is the number of dimensions used in the dataspace
                                             , mdims    // mdims is a one-dimensional array of size rank specifying the size of each dimension of the dataset.
                                             , NULL     // maxdims is an array of the same size specifying the upper limit on the size of each dimension. maxdims may be NULL, in which case the maxdims == mdims.
                                                        ); // must be released with H5Sclose or resource leaks will occur

        assert( memspace_id > -1 );
        //logger << "h5g_Write: memspace_id created." << std::endl;


        /* Create the dataspace for the dataset.
         * The dataspace describes the dimensions of the dataset array. */

        hsize_t mindims = 50000;

        hsize_t dims[ dim ];

        for(UINT i=0;i<dim;i++){
          dims[dim-i-1] = (hsize_t) ( inputdata.domain_data.nCells[i] );
          mindims = std::min( mindims, dims[dim-i-1] );
        }


        //logger << "dims [0] = "<< dims[0]<<std::endl;
        //logger << "dims [1] = "<< dims[1]<<std::endl;
        hid_t dataspace_id = H5Screate_simple(
                                              dim     // number of dimensions = rank
                                              , dims  // vector containing sizes per dimension
                                              , NULL  // maxdims == dims
                                              );
        assert( dataspace_id > -1 );
        //logger << "h5g_Write: dataspace_id created." << std::endl;



        hid_t plist_id = H5Pcreate(                    // The new property list is initialized with default values for the specified class.
                                   H5P_DATASET_CREATE  // Properties for dataset creation
                                                       );
        assert( plist_id > -1 );
        //logger << "h5g_Write: plist_id created." << std::endl;


        herr_t status;

        hsize_t maxnchunks = 1;
        for(UINT i=0;i<dim;i++){
          maxnchunks = std::max( maxnchunks, (hsize_t) inputdata.domain_data.nChunks[i] );
        }


        hsize_t chunk_dims[ dim ];
        if( maxnchunks >= mindims ){
          logger << "Warning: maxnchunks > mindims, set chunk_dims[i]=0" << std::endl;
          for(UINT i=0;i<dim;i++)
            chunk_dims[i]=(hsize_t)1;
        }
        else{
          for(UINT i=0;i<dim;i++){
            chunk_dims[dim-i-1] = (hsize_t) inputdata.domain_data.nChunks[i];
          }
        }


        status = H5Pset_chunk( plist_id
                               , dim             // must be == rank of the dataset
                               , chunk_dims      // The values of the check_dims array define the size of the chunks to store the dataset's raw data. The unit of measure for check_dims values is dataset elements.
                               );

        if( status < 0 )
          logger << "Warning : H5Pset_chunk() failed." << std::endl;

        assert( status > -1 );
        //logger << "h5g_Write: H5Pset_chunk() o.k." << std::endl;


        status = H5Pset_shuffle( plist_id ); // Sets the shuffle filter, H5Z_FILTER_SHUFFLE, in the dataset creation property list. This re-orders data to simplify compression.
        assert( status > -1 );
        //logger << "h5g_Write: H5Pset_shuffle() o.k." << std::endl;


        status = H5Pset_deflate( plist_id, 1 ); // Sets deflate (GNU gzip) compression method and compression level. ( 0 < level < 9, lower = faster, but less compression )
        assert( status > -1 );
        //logger << "h5g_Write: H5Pset_deflate() o.k." << std::endl;



        /* Create the dataset. */
        hid_t dataset_id = H5Dcreate( file_id,            // Location identifier: id of the file or the group within which to create the dataset
                                      groupname.c_str(),  // Dataset name: may be either an absolute path in the file or a relative path from file_id naming the dataset.
                                      HDF5_DATA_TYPE,     // Datatype identifier, here: IEEE floating point 32-bit or 64-bit, see my_macros.hh
                                      dataspace_id ,      // Dataspace identifier
                                      plist_id            // Dataset creation property list identifier
                                      );
        assert( dataset_id > -1 );
        //logger << "h5g_Write: dataset_id created." << std::endl;


        /* Write the dataset. */
        status = H5Dwrite(                         // Writes raw data from a buffer to a dataset.
                          dataset_id               // dataset identifier
                          , H5T_NATIVE_DOUBLE      // memory datatype id
                          , memspace_id            // specifies the memory dataspace and the selection within it
                          , H5S_ALL                // specifies the selection within the file dataset's dataspace. H5S_ALL indicates that the entire file dataspace, as defined by the current dimensions of the dataset, is to be selected
                          , H5P_DEFAULT            // Identifier of a transfer property list for this I/O operation. H5P_DEFAULT: The default data transfer properties are used.
                          , &(data[0])    // application memory buffer
                                                   );
        assert( status > -1 );
        //logger << "h5g_Write: H5Dwrite( " <<groupname.c_str()<<" ) o.k." << std::endl;

        status = H5Sclose( dataspace_id );
        assert( status > -1 );
        //logger << "h5g_Write: dataspace_id closed." << std::endl;

        status = H5Sclose( memspace_id );
        assert( status > -1 );
        //logger << "h5g_Write: memspace_id closed." << std::endl;

        status = H5Dclose( dataset_id );
        assert( status > -1 );
        //logger << "h5g_Write: dataset_id closed." << std::endl;


        /* Close the property list. */
        status = H5Pclose( plist_id );
        assert( status > -1 );
        //logger << "h5g_Write: H5Pclose(plist_id) done." << std::endl;

        /* Close the file. */
        status = H5Fclose( file_id );
        assert( status > -1 );
        //logger << "h5g_Write: H5Fclose( file_id ) done." << std::endl;

        std::stringstream jobtitle;
        jobtitle << "h5g_Write: writing " << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   inputdata.verbosity,
                                   "IO",
                                   jobtitle.str() );

      }







      /** function to write data (in parallel) to the hdf5 file
       *
       * \tparam gv is the gridview
       * \tparam inputdata is the inputdata object
       * \tparam data are the data which are written into the file (current hyperslab)
       * \tparam groupname is the name/path where the data are stored in the HDF5 file
       * \tparam nCells  the number of virtual(global) cells for each dimension
       * \tparam nChunks = numbers of chunks (needed for hdf5) in each dimension
       * \tparam gridsizes = virtual grid-sizes, needed only for extra data (not used so far)
       * \tparam filename = name of the hdf5 file where all the data will be stored
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */
      template<typename GV>
      static void h5g_pWrite( const GV& gv,
                              const Vector<REAL>& datavector,
                              const std::string& filename,
                              const std::string& groupname,
                              const CInputData& inputdata,
                              const Vector<UINT>& nCells,
                              const int blocksizePerDimension = 1,
                              const FEMType::Type femtype=FEMType::DG,
                              const int current_grid_level=0,
                              const bool preserve_structure=true
                              )
      {

        Dune::Timer watch;

        logger << "h5g_pWrite: ... "
               << filename << std::endl;

        //get the communicator of the dune grid
        //const Dune::MPIHelper::MPICommunicator &comm = gv.comm();

        //Info variable need for the HDF5
        MPI_Info mpiInfo = MPI_INFO_NULL;

        //some needed variables and typdefs
        const UINT dim = GV::dimension;  // dimensionality of the problem (2-d or 3-d)

        // DO NOT CHANGE THIS EVER: Leave it to be "Dune::All_Partition".
        // Since datavector contains ALL elements of the vector backend,
        // we must consider over ALL cells. Neglecting overlap cells would
        // shift the solution on the HDF file, leading to wrong inversion results!
        const Dune::PartitionIteratorType partitiontype
          = Dune::All_Partition;


        typedef typename GV::template Codim<0>::template Partition<partitiontype>::Iterator LeafIterator;

        //typedef typename GV::Grid GRIDTYPE;

        //variables for defining the offset and the count in the different dimensions
        Vector<UINT> max_index( dim, 0 );
        Vector<UINT> min_index( dim, 100000 );

        /*
          1.) Loop over the current hyperslab (hyperslab = all the leaf grid-cells that belong to the current processor) and get the global coordinate 'x' of each grid-cell center.
          Remark: From now on, we are working solemnly with the virtual Y-field grid, using virtual indices. Later, in the evaluate2d() and evaluate3d() member functions of the class kfieldgenerator, we are doing the very same conversion from a global coordinate to the virtual global index. That's why it can work independently of the grid type!
          2.) Translate this global coordinate 'x' into its corresponding Yfield-tensor-index 'global_index'. This requires only the virtual gridsize of the virtual Yfield grid. So 'global_index' is just a virtual index.
          3.) 'min_index' corresponds to the virtual cell with the lowest coordinate values in the current hyperslab. It is the offset (= translational displacement) from the very first virtual cell (0,0).
          4.) 'local_count[i]' is the number of virutal cells in the dimension 'i' of the current hyperslab.
          Remark: The volume of the current hyperslab is 'nAllLocalCells'. 'mappersize' == 'nAllLocalCells' only if the resolution of the yaspgrid == resolution of the virtual grid !

        */


        Vector<REAL> gridsizes (inputdata.domain_data.virtual_gridsizes);

        int gridLevelFactor = 1;
        if(current_grid_level>0)
          gridLevelFactor = std::pow(2,current_grid_level);

        if( current_grid_level > 0 )
          gridsizes /= REAL( gridLevelFactor );

        logger << "h5g_pWrite: gridsizes = "
               << gridsizes << std::endl;

        for (LeafIterator it = gv.template begin<0,partitiontype> ()
               ; it != gv.template end<0,partitiontype> ()
               ; ++it)
          {

            if( !preserve_structure ){
              // Do not consider cells on overlap!
              if( it->partitionType() != Dune::InteriorEntity ){
                continue;
              }
            }

            // Get the local coordinate:
            Dune::FieldVector<CTYPE,dim> x = it->geometry().center();

            std::vector<UINT> global_index(dim,0);
            for(UINT i=0; i<dim; i++ ) {
              global_index[i] = static_cast<UINT>( x[i] / gridsizes[i] );
              max_index[i] = std::max( global_index[i], max_index[i] );
              min_index[i] = std::min( global_index[i], min_index[i] );
            }

            //logger << "DEBUG: global-index = ("  << global_index[0] << "," << global_index[1]
            //<< ")" << std::endl;

            // Idee für UG später:
            // lese direkt punktweise aus der h5-Datei, über den offset[] = global_index[], mit count = (1 1 1) ?
          } // end of loop over leaf elements


#ifdef DEBUG_LOG
        logger << "DEBUG: h5g_pWrite: nCells = "
               << nCells
               << std::endl;

        logger << "DEBUG: h5g_pWrite: min_index = "
               << min_index
               << std::endl;

        logger << "DEBUG: h5g_pWrite: max_index = "
               << max_index
               << std::endl;

        logger << "DEBUG: h5g_pWrite: hdf5 functions ... "
               << std::endl;
#endif

        MPI_Barrier( gv.comm() );

        //int mpi_rank = gv.comm().rank();


        // Set up file access property list with parallel I/O access
        hid_t plist_id = H5Pcreate( H5P_FILE_ACCESS );
        //logger << "h5g_pWrite: H5Pcreate() done. plist_id = "
        //       << plist_id << std::endl;
        assert( plist_id > -1 );


        herr_t status = H5Pset_fapl_mpio( plist_id, gv.comm(), mpiInfo );
        //logger << "h5g_pWrite: H5Pset_fapl_mpio: status = "
        // << status << std::endl;
        if( status < 0 )
          logger << "Warning : H5Pset_fapl_mpio() failed." << std::endl;
        assert( status > -1 );

        //  Create a new file collectively and release property list identifier.
        //logger << "h5g_pWrite: "
        //       << " create file "
        //       << filename.c_str()
        //       << std::endl;

        hid_t file_id = H5Fcreate( filename.c_str(),
                                   H5F_ACC_TRUNC,
                                   H5P_DEFAULT,
                                   plist_id );
        H5Pclose( plist_id );
        //logger << "h5g_pWrite: "
        //       << filename.c_str()
        //       << " file created!" << std::endl;
        assert( file_id > -1 );

        //==========================================================
        // "dims_global" defines the filespace of the HDF5 file:
        // It is process-independent and contains all DOFs.
        //==========================================================
        hsize_t dims_global[ dim];

        for(UINT i=0; i<dim; i++) {
          dims_global[dim -i-1] = nCells[i] * blocksizePerDimension * gridLevelFactor;
          if(femtype==FEMType::CG){
            dims_global[dim -i-1] += 1;
          }
        }

        /*
          std::cout << "DEBUG:  p=" << mpi_rank
          << "  dims_global = " << dims_global[0] << "," << dims_global[1]
          #ifdef DIMENSION3
          << "," << dims_global[2]
          #endif
          << std::endl;
        */

        //  Create the dataspace for the dataset.
        hid_t filespace = H5Screate_simple( dim, dims_global, NULL );
        assert(filespace>-1);
        //logger << "h5g_pWrite:  fileSPACE created!" << std::endl;


        // Create the dataset with default properties and close filespace.
        hid_t dset_id = H5Dcreate( file_id,
                                   groupname.c_str(),
                                   HDF5_DATA_TYPE,    // Datatype identifier, here: IEEE floating point 32-bit or 64-bit, see my_macros.hh
                                   filespace,
                                   H5P_DEFAULT );
        H5Sclose(filespace);
        assert(dset_id>-1);
        //logger<< "h5g_pWrite:  dataset created!" << std::endl;


        // set the count in the different dimensions (determine the size of the hyperslab)
        hsize_t count[ dim ];


        for(UINT i=0; i<dim; i++ ) {

          if(femtype==FEMType::CG){
            // be careful!!! needed when vertex-based data instead of cell based data
            count[dim-i-1] = max_index[i] - min_index[i] + 2; // Q1-FEM: in each dimension, #nodes = #elements+1
          }
          else{
            count[dim-i-1] = max_index[i] - min_index[i] + 1;
            count[dim-i-1] *= blocksizePerDimension;
          }
#ifdef DEBUG_LOG
          logger << "DEBUG: h5g_pWrite: count[" << dim-i-1 << "] = " << count[dim-i-1] << std::endl;
#endif
        }


        //define the total size of the local data

        hsize_t nAllLocalCells = count[0] * count[1];
#ifdef DIMENSION3
        nAllLocalCells *= count[2];
#endif

        //set the offset of the data!
        hsize_t offset[ dim ];
        for(UINT i=0; i<dim; i++ ) {

          if(femtype==FEMType::CG)
            offset[dim-i-1] = min_index[i];
          else
            offset[dim-i-1] = min_index[i] * blocksizePerDimension;

        }

        /*
          std::cout << "DEBUG:  p=" << mpi_rank
          << "  count =" << count[0] << "," << count[1]
          #ifdef DIMENSION3
          << "," << count[2]
          #endif
          << std::endl;
        */

        // Each process defines dataset in memory and writes it to the hyperslab in the file.

        hid_t memspace_id = -1;
        if(nAllLocalCells!=0){  // -> otherwise HDF5 warning, because of writing nothing!
          memspace_id = H5Screate_simple(dim, count, NULL);
          assert(memspace_id>-1);
        }


        //logger<< "h5g_pWrite:  memspace created!" << std::endl;

        // Select hyperslab in the file.
        filespace = H5Dget_space( dset_id );

        /*
          REAL d = datavector.one_norm();
          std::cout << "DEBUG:  p= " << mpi_rank
          << "  |datavector| = " << d
          << "  min_index = " << min_index[0] << "," << min_index[1]
          #ifdef DIMENSION3
          << "," << min_index[2]
          #endif
          << "  max_index= " << max_index[0] << "," << max_index[1]
          #ifdef DIMENSION3
          << "," << max_index[2]
          #endif
          << std::endl;

          std::cout << "DEBUG:  p=" << mpi_rank
          << "  offset = "  << offset[0] << "," << offset[1]
          #ifdef DIMENSION3
          << "," << offset[2]
          #endif
          << std::endl;
        */

        H5Sselect_hyperslab( filespace,
                             H5S_SELECT_SET,
                             offset,
                             NULL,   // <==> stride={1,1,1}
                             count,
                             NULL    // <==> block={1,1,1}
                             );
        //logger<< "h5g_pWrite:  hyperslab selected!" << std::endl;

        // Create property list for collective dataset write.
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_COLLECTIVE );
        //H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );

        //logger<< "h5g_pWrite:  properties set!" << std::endl;

        // finally write the data to the disk
        //logger<< "h5g_pWrite:  writing ... " << std::endl;

        if( nAllLocalCells != 0 ){
          status = H5Dwrite( dset_id,
                             H5T_NATIVE_DOUBLE,
                             memspace_id,
                             filespace,
                             plist_id,
                             &( datavector[0] )
                             );
          assert(status>-1);
        }
        else{ // avoid hdf5 blocking
          status = H5Dwrite( dset_id,
                             H5T_NATIVE_DOUBLE,
                             H5S_ALL,
                             filespace,
                             plist_id,
                             &( datavector[0] )
                             );
          assert(status>-1);
        }

        logger << "h5g_pWrite:  .... done (writing) "
               << filename
               << std::endl;

        // close the used identifyers
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace_id);
        H5Pclose(plist_id);
        H5Fclose(file_id);

        std::stringstream jobtitle;
        jobtitle << "h5g_pWrite: writing " << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   inputdata.verbosity,
                                   "IO",
                                   jobtitle.str() );


      }

      /** function to retrieve data (sequential) from the hdf5 file
       *
       * \tparam local_data is a 'return value' which gets the data that belongs to the current processor (current hyperslab)
       * \tparam groupname is the name/path where the data are stored in the HDF5 file
       * \tparam local_count is a 'return value' which gets the number of virtual cells for each dimension
       * \tparam local_offset is a 'return value' which gets the distance of the current hyperslab from the origin.
       * \tparam filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */



      template<typename IDT>
      static void h5g_Read( Vector<REAL>& local_data
                            , const std::string& filename
                            , const std::string& groupname
                            , Vector<UINT>& local_offset
                            , Vector<UINT>& local_count
                            , const IDT& inputdata
                            , const int iPart=1
                            , const int nParts=1
                            , const bool bTimerOn=true
                            )
      {
        Dune::Timer watch;

        //logger << General::getDateAndTime() << "h5g_Read: ... " << filename << std::endl;

        // open the file for reading
        hid_t file_id=  H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        assert( file_id > -1 );
        //logger << "h5g_Read(): h5-file open: " << filename.c_str() << std::endl;

        // open the dataset
        hid_t dataset_id = H5Dopen(file_id, groupname.c_str());
        assert( dataset_id > -1 );
        //logger << "h5g_Read(): H5Dopen dataset done." << std::endl;

        // gwt the dataspace
        hid_t dataspace_id=H5Dget_space (dataset_id);
        assert(dataspace_id>-1);
        //logger << "h5g_Read(): H5Dget_space  done." << std::endl;

        // some needed variables
        herr_t status;
        hsize_t dim,*dims;

        // get the dimension (2-d or 3-d)
        //logger << "Getting parameters for reading from data (read all)"<<std::endl;
        dim = H5Sget_simple_extent_ndims( dataspace_id );
        //logger << "Got dim = " << dim << std::endl;
        //logger << "h5g_Read()! ( dim = "<<dim<< " )" << std::endl;

        // get the size of the problem
        dims=(hsize_t*)malloc(dim * sizeof (hsize_t));

        status = H5Sget_simple_extent_dims( dataspace_id , dims, 0);
        if( status < 0 )
          logger << "Warning : H5Sget_simple_extent_dims() failed." << std::endl;

        assert( status > -1 );

        // file the return value. the local_count
        local_count.resize(dim);
        //logger << "Got local_count from dataspace:" << std::endl;
        for( UINT i=0; i<dim; i++ ) {
          local_count[dim-1-i] = dims[i];
          //logger << "local_count[" << dim-i-1 << "] = " << local_count[dim-1-i] << std::endl;
        }

        // set offset(for reading) and offset_local(as return value) to 0 (read all the data, maybe later also used to read only parts, maybe we can kick it out first)
        // like this it is not nice but works, maybe we change it later!!!
        // local_size is the size of the data
        hsize_t offset[dim];
        local_offset.resize(dim);
        hsize_t  local_size=1;
        for( UINT i=0; i<dim; i++ ){
          local_offset[i]=0;
          offset[i]=0;
        }

        if(nParts>1){
          local_count[1]/=nParts;
          local_offset[1]=local_count[1]*iPart;
          offset[1]=local_offset[1];
          if(iPart+1==nParts){
            local_count[1]=dims[1]-local_offset[1];
          }
          for( UINT i=0; i<dim; i++ ){
            dims[i]=local_count[dim-1-i];
          }
        }

        //logger << "local_count[" << 1 << "] = " << local_count[1] << std::endl;
        for( UINT i=0; i<dim; i++ ){
          local_size *= local_count[i];
        }
        // create the memory space
        hid_t memspace_id = H5Screate_simple (dim, dims, NULL);
        //hid_t memspace_id = H5Screate_simple (1, &local_size, NULL);

        //select the hyperslab-
        status = H5Sselect_hyperslab( dataspace_id,
                                      H5S_SELECT_SET,
                                      offset,
                                      NULL, // <==> stride={1,1,1}
                                      dims,
                                      NULL  // <==> block={1,1,1}
                                      );
        assert(status>-1);
        logger << General::getDateAndTime() << "h5g_Read(): hyperslab selected." << std::endl;

        //resize the return data
        if( local_data.size() != local_size )
          local_data.resize( local_size /*mappersize*/ );

        /* set up the collective transfer properties list */
        hid_t xferPropList = H5Pcreate (H5P_DATASET_XFER);
        assert( xferPropList > -1 );
        logger << General::getDateAndTime() << "xferPropList created." << std::endl;

        // finally the reading from the file
        status = H5Dread( dataset_id
                          , H5T_NATIVE_DOUBLE
                          , memspace_id
                          , dataspace_id // H5S_ALL //
                          , xferPropList //H5P_DEFAULT //
                          , &( local_data[0] )
                          );
        assert(status>-1);
        logger << General::getDateAndTime() << "h5g_Read(): reading data done." << std::endl;

        // close the identifyers
        H5Dclose (dataset_id);
        H5Sclose (dataspace_id);
        H5Sclose (memspace_id);
        free(dims);
        status = H5Fclose( file_id );
        assert( status > -1 );
        logger << General::getDateAndTime() << "h5g_Read: H5Fclose( file_id ) done." << std::endl;


        std::stringstream jobtitle;
        jobtitle << "h5g_Read: reading " << filename;
        if( bTimerOn ){
          // log and add to elapsed_IO
          General::log_elapsed_time( watch.elapsed(),
                                     inputdata.verbosity,
                                     "IO",
                                     jobtitle.str() );
        }
        else {
          // log only to logger
          if( inputdata.verbosity >= VERBOSITY_TIMER_DETAILS ){
            logger << "Duration:["
                   << std::setw(6)
                   << "IO" << "]"
                   << std::setw(15)
                   << std::setprecision(4)
                   << watch.elapsed() << " sec."
                   << " for "
                   << jobtitle
                   << std::endl;
          }

        }

      }










#ifdef PARALLEL
      /** function to retrieve data (in parallel) from the hdf5 file
       *
       * \tparam gv is the gridview
       * \tparam inputdata is the inputdata object
       * \tparam local_data is a 'return value' which gets the data that belongs to the current processor (current hyperslab)
       * \tparam groupname is the name/path where the data are stored in the HDF5 file
       * \tparam local_count is a 'return value' which gets the number of virtual cells for each dimension
       * \tparam local_offset is a 'return value' which gets the distance of the current hyperslab from the origin.
       * \tparam filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */

      template<typename GV,
               Dune::PartitionIteratorType partitiontype = Dune::All_Partition>
      static void h5g_pRead( const GV& gv,
                             Vector<REAL>& local_data,
                             const std::string& filename,
                             const std::string& groupname,
                             Vector<UINT>& local_offset,
                             Vector<UINT>& local_count,
                             const CInputData& inputdata,
                             const int blocksizePerDimension=1,
                             const FEMType::Type femtype=FEMType::DG,
                             const int current_grid_level=0
                             //, const bool bAllPartition=true
                             )
      {

        // Dune::PartitionIteratorType partitiontype = Dune::All_Partition;
        // if( !bAllPartition )
        //   partitiontype = Dune::Interior_Partition;

        Dune::Timer watch;

        logger << General::getDateAndTime() << "h5g_pRead: " << filename << std::endl;

        const Dune::MPIHelper::MPICommunicator &comm = gv.comm();

        MPI_Info mpiInfo = MPI_INFO_NULL;

        const UINT dim = GV::dimension;
        typedef typename GV::template Codim<0>::template Partition<partitiontype>::Iterator LeafIterator;
        //typedef typename GV::Grid GRIDTYPE;

        std::vector<UINT> max_index( dim, 0 );
        std::vector<UINT> min_index( dim, 100000 );

        /*
          1.) Loop over the current hyperslab (hyperslab = all the leaf grid-cells that belong to the current processor) and get the global coordinate 'x' of each grid-cell center.
          Remark: Form now on, we are working solemnly with the virtual Y-field grid, using virtual indices. Later, in the evaluate2d() and evaluate3d() member functions of the class kfieldgenerator, we are doing the very same conversion from a global coordinate to the virtual global index. That's why it can work independently of the grid type!
          2.) Translate this global coordinate 'x' into its corresponding Yfield-tensor-index 'global_index'. This requires only the virtual gridsize of the virtual Yfield grid. So 'global_index' is just a virtual index.
          3.) 'min_index' corresponds to the virtual cell with the lowest coordinate values in the current hyperslab. It is the offset (= translational displacement) from the very first virtual cell (0,0).
          4.) 'local_count[i]' is the number of virutal cells in the dimension 'i' of the current hyperslab.
          Remark: The volume of the current hyperslab is 'nAllLocalCells'. 'mappersize' == 'nAllLocalCells' only if the resolution of the yaspgrid == resolution of the virtual grid !

        */

        Vector<REAL> gridsizes (inputdata.domain_data.virtual_gridsizes);
        if( current_grid_level > 0 )
          gridsizes /= REAL( std::pow(2,current_grid_level) );


        for (LeafIterator it = gv.template begin<0,partitiontype> ()
               ; it != gv.template end<0,partitiontype> ()
               ; ++it) {
          // Get the local coordinate:
          Dune::FieldVector<CTYPE, dim> x = it->geometry().center();

          std::vector<UINT> global_index(dim,0);
          for(UINT i=0; i<dim; i++ ) {
            global_index[i] = static_cast<UINT>( x[i] / gridsizes[i] );
            max_index[i] = std::max( global_index[i], max_index[i] );
            min_index[i] = std::min( global_index[i], min_index[i] );
          }
          // logger << "DEBUG: global-index = "  << global_index << std::endl;
        } // end of loop over leaf elements

        logger << General::getDateAndTime() << "min/max_index loop done" << std::endl;



        /* Do some preparations for accessing the hdf5 file...
         * This hyperslab selection works only for a structured grid, I guess!
         * For an unstructured grid, a pointwise hyperslab selection should still be possible.
         */

        hsize_t count[ dim ];
        local_count.resize( dim );
        local_offset.resize( dim );

        for(UINT i=0; i<dim; i++ )
          {
            local_count[i] = max_index[i] - min_index[i] + 1;  // So far, this is the number of elements per dimension.
#ifdef DEBUG_LOG
            logger << "DEBUG: h5g_pRead: local_count[" << i << "] = " << local_count[i] << std::endl;
#endif
          }

        for(UINT i=0; i<dim; i++ )
          {
            count[i] = (hsize_t) local_count[ dim - i - 1 ];
#ifdef DEBUG_LOG
            logger << "DEBUG: h5g_pRead: count[" << i <<"] = " << count[i] << std::endl;
#endif
          }

        hsize_t offset[ dim ];
        for(UINT i=0; i<dim; i++ ) {
          if(femtype==FEMType::CG)
            local_offset[i] = min_index[i];
          else
            local_offset[i] = min_index[i] * blocksizePerDimension;
#ifdef DEBUG_LOG
          logger << "DEBUG: h5g_pRead: local_offset[" << i << "] = " << local_offset[i] << std::endl;
#endif
        }

        for(UINT i=0; i<dim; i++ )
          {
            offset[i] = (hsize_t) local_offset[ dim - i - 1 ];
#ifdef DEBUG_LOG
            logger << "DEBUG: h5g_pRead: offset[" << i <<"] = " << offset[i] << std::endl;
#endif
          }


        /* setup file access template with parallel IO access. */
        hid_t access_pList = H5Pcreate( H5P_FILE_ACCESS );
        assert( access_pList > -1 );
        //logger << "access_pList created." << std::endl;

        herr_t status;             /* Generic return value */
        status = H5Pset_fapl_mpio( access_pList, comm, mpiInfo ); // Take care, this requires the MPI version of hdf5!
        if( status < 0 )
          logger << "Warning : H5Pset_fapl_mpio() failed." << std::endl;
        assert( status > -1 );
        //logger << "access_pList set" << std::endl;



        logger << General::getDateAndTime() << "H5Fopen ..." << filename << std::endl;
        /* open the file collectively */
        hid_t infileId = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, access_pList );
        assert( infileId > -1 );
        //logger << "infileId opened." << std::endl;





        /* Release file-access template */
        status = H5Pclose( access_pList );
        assert( status > -1 );
        //logger << "access_pList closed." << std::endl;

        hid_t indatasetId = H5Dopen( infileId, groupname.c_str() );
        assert( indatasetId > -1 );
        //logger << "indatasetId opened." << std::endl;

        hid_t indataspaceId = H5Dget_space( indatasetId );
        assert( indataspaceId > -1 );
        //logger << "indataspaceId opened." << std::endl;

        //hsize_t rank = H5Sget_simple_extent_ndims( indataspaceId );
        //logger << "Got rank = " << rank << std::endl;

        hsize_t dims[ dim ];
        status = H5Sget_simple_extent_dims( indataspaceId, dims, 0);
        assert( status > -1 );
        // logger << "Got dims from dataspace:" << std::endl;

        /*
         * Be careful!!!
         * correction if vertex based data a read!
         */
        if(femtype==FEMType::CG) {
          for(UINT ii=0;ii<dim;ii++){

            //logger << "dims["<<dim-1-ii<<"] = " << dims[dim-1-ii] << std::endl;
            //logger << "local_count["<<ii<<"] = " << local_count[ii] << std::endl;
            //logger << "count["<<dim-1-ii<<"] = " << count[dim-1-ii] << std::endl;

            if(dims[dim-1-ii]-inputdata.domain_data.nCells[ii]){
              local_count[ii]+=1;    // Q1-FEM: #nodes = #elements + 1 per dimension
              count[dim-1-ii]+=1;
            }

            //logger << "local_count["<<ii<<"] = " << local_count[ii] << std::endl;
            //logger << "count["<<dim-1-ii<<"] = " << count[dim-1-ii] << std::endl;
          }
        }
        else{
          for(UINT ii=0;ii<dim;ii++){


            //if(dims[dim-1-ii]-inputdata.domain_data.nCells[ii]){
            local_count[ii] = blocksizePerDimension * local_count[ii];
            count[dim-1-ii] = blocksizePerDimension * count[dim-1-ii];
            //}

#ifdef DEBUG_LOG
            logger << "DEBUG: h5g_pRead: dims["<<dim-1-ii<<"] = " << dims[dim-1-ii] << std::endl;
            logger << "DEBUG: h5g_pRead: local_count["<<ii<<"] = " << local_count[ii] << std::endl;
            logger << "DEBUG: h5g_pRead: count["<<dim-1-ii<<"] = " << count[dim-1-ii] << std::endl;
#endif
          }
        }


        hsize_t nAllLocalCells = local_count[0] * local_count[1];
#ifdef DIMENSION3
        nAllLocalCells *= local_count[2];
#endif

        hid_t memspaceId = H5Screate_simple(1, &nAllLocalCells, NULL);
        assert(memspaceId > -1 );

        logger << General::getDateAndTime() << "H5Sselect_hyperslab() ..." << std::endl;

        status = H5Sselect_hyperslab( indataspaceId
                                      , H5S_SELECT_SET
                                      , offset
                                      , NULL // <==> stride={1,1,1}
                                      , count
                                      , NULL // <==> block={1,1,1}
                                      );
        assert( status > -1 );

        logger << General::getDateAndTime() << "H5Sselect_hyperslab() done.  H5Dread() ..." << std::endl;

        /* set up the collective transfer properties list */
        hid_t xferPropList = H5Pcreate (H5P_DATASET_XFER);
        assert( xferPropList > -1 );
        //logger << "DEBUG: xferPropList created." << std::endl;

        //Vector<double> local_Yfield_vector;

        if( local_data.size() != nAllLocalCells )
          local_data.resize( nAllLocalCells );

        status = H5Dread( indatasetId
                          , H5T_NATIVE_DOUBLE //image.DataType
                          , memspaceId
                          , indataspaceId // H5S_ALL //
                          , xferPropList // H5P_DEFAULT
                          , &( local_data[0] )
                          );
        assert( status > -1 );
        logger << General::getDateAndTime() << "H5Dread() done." << std::endl;


        status = H5Sclose(indataspaceId);
        assert( status > -1 );
        //logger << "indataspaceId closed." << std::endl;

        status = H5Sclose(memspaceId);
        assert( status > -1 );
        //logger << "memspaceId closed." << std::endl;

        status = H5Pclose( xferPropList );
        assert( status > -1 );
        //logger << "xferPropList closed." << std::endl;

        /* Close the dataset. */
        status = H5Dclose( indatasetId );
        assert( status > -1 );
        //logger << "indatasetId closed." << std::endl;

        /* Close the file. */
        status = H5Fclose( infileId );
        assert( status > -1 );
        logger << General::getDateAndTime() << "infileId closed." << std::endl;


        /*

          The following code is just for debugging purposes!

        */



        std::stringstream jobtitle;
        if( partitiontype != Dune::Interior_Partition )
          jobtitle << "h5g_pRead: overlapped reading " << filename;
        else
          jobtitle << "h5g_pRead: nonoverlapped reading " << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   comm,
                                   inputdata.verbosity,
                                   "IO",
                                   jobtitle.str() );


        return;
      }
#endif // ifdef PARALLEL

      /** function to write a vector (sequential) to a hdf5 file
       *
       * \tparam dimensions the dimensions of the written data vector
       * \tparam data this data will be written to the file
       * \tparam groupname is the name/path where the data are stored in the HDF5 file
       * \tparam filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */
      static void h5_Write( const Vector<REAL> & data
                            , const std::string & filename
                            , const std::string & groupname
                            , const Vector<UINT> dimensions
                            )
      {
        // ge the dimensionality of the data
        UINT dim=dimensions.size();

        //number of elements the data vector should have
        UINT n=1;
        for(UINT i=0;i<dim;i++)
          n*=dimensions[i];

        //the data make no sense!
        assert(n==data.size());

        /* Create a new file using default properties. */
        hid_t file_id = H5Fcreate(
                                  filename.c_str()    // name of the file to be created
                                  , H5F_ACC_TRUNC          // if you are trying to create a file that exists already, the existing file will be truncated, i.e., all data stored on the original file will be erased
                                  , H5P_DEFAULT            // default file creation property list
                                  , H5P_DEFAULT            // default file access property list
                                  );

        assert( file_id > -1 );

        // create the memspace
        hsize_t mdims[1];
        mdims[0] = n;
        hid_t memspace_id = H5Screate_simple(           // H5Screate_simple creates a new simple dataspace and opens it for access, returning a dataset identifier.
                                             1          // rank=1 is the number of dimensions used in the dataspace
                                             , mdims    // mdims is a one-dimensional array of size rank specifying the size of each dimension of the dataset.
                                             , NULL     // maxdims is an array of the same size specifying the upper limit on the size of each dimension. maxdims may be NULL, in which case the maxdims == mdims.
                                                        ); // must be released with H5Sclose or resource leaks will occur

        assert( memspace_id > -1 );

        /* Create the dataspace for the dataset.
         * The dataspace describes the dimensions of the dataset array.
         */
        hsize_t dims[ dim ];
        for(UINT i=0;i<dim;i++)
          dims[i]=dimensions[i];

        hid_t dataspace_id = H5Screate_simple(
                                              dim     // number of dimensions = rank
                                              , dims  // vector containing sizes per dimension
                                              , NULL  // maxdims == dims
                                              );
        assert( dataspace_id > -1 );


        hid_t plist_id = H5Pcreate(                    // The new property list is initialized with default values for the specified class.
                                   H5P_DATASET_CREATE  // Properties for dataset creation
                                                       );
        assert( plist_id > -1 );

        //needed variable for HDF5 operations
        herr_t status;


        /*
         * chunks
         * !!! not optimized. Might be better if larger! -> use the inputdata
         */
        hsize_t chunk_dims[ dim ];
        for(UINT i=0;i<dim;i++)
          chunk_dims[i]=1;  //maybe to small for fast saving!

        // set the chunk size!
        status = H5Pset_chunk(
                              plist_id
                              , dim             // must be == rank of the dataset
                              , chunk_dims      // The values of the check_dims array define the size of the chunks to store the dataset's raw data. The unit of measure for check_dims values is dataset elements.
                              );
        assert( status > -1 );

        status = H5Pset_shuffle( plist_id ); // Sets the shuffle filter, H5Z_FILTER_SHUFFLE, in the dataset creation property list. This re-orders data to simplify compression.
        assert( status > -1 );

        status = H5Pset_deflate( plist_id, 1 ); // Sets deflate (GNU gzip) compression method and compression level. ( 0 < level < 9, lower = faster, but less compression )
        assert( status > -1 );

        /* Create the dataset. */
        hid_t dataset_id = H5Dcreate(  file_id,            // Location identifier: id of the file or the group within which to create the dataset
                                       groupname.c_str(),  // Dataset name: may be either an absolute path in the file or a relative path from file_id naming the dataset.
                                       HDF5_DATA_TYPE,    // Datatype identifier, here: IEEE floating point 32-bit or 64-bit, see my_macros.hh
                                       dataspace_id,      // Dataspace identifier
                                       plist_id           // Dataset creation property list identifier
                                       );
        assert( dataset_id > -1 );

        /* Write the dataset. */
        status = H5Dwrite(                         // Writes raw data from a buffer to a dataset.
                          dataset_id               // dataset identifier
                          , H5T_NATIVE_DOUBLE      // memory datatype id
                          , memspace_id            // specifies the memory dataspace and the selection within it
                          , H5S_ALL                // specifies the selection within the file dataset's dataspace. H5S_ALL indicates that the entire file dataspace, as defined by the current dimensions of the dataset, is to be selected
                          , H5P_DEFAULT            // Identifier of a transfer property list for this I/O operation. H5P_DEFAULT: The default data transfer properties are used.
                          , &(data[0])    // application memory buffer
                                                   );
        assert( status > -1 );

        /*
         * close everything
         */
        status = H5Sclose( dataspace_id );
        if( status < 0 )
          logger << "Warning : H5Sclose() failed." << std::endl;
        assert( status > -1 );

        status = H5Sclose( memspace_id );
        if( status < 0 )
          logger << "Warning : H5Sclose() failed." << std::endl;
        assert( status > -1 );

        status = H5Dclose( dataset_id );
        if( status < 0 )
          logger << "Warning : H5Dclose() failed." << std::endl;
        assert( status > -1 );

        status = H5Pclose( plist_id );
        if( status < 0 )
          logger << "Warning : H5Pclose() failed." << std::endl;
        assert( status > -1 );

        status = H5Fclose( file_id );
        if( status < 0 )
          logger << "Warning : H5Fclose() failed." << std::endl;
        assert( status > -1 );
      };



      /** function to retrieve a vector (sequential) from the hdf5 file
       *
       * \tparam local_data is a 'return value' which gets the data that belongs to the current processor (current hyperslab)
       * \tparam groupname is the name/path where the data are stored in the HDF5 file
       * \tparam filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */

      static void h5_Read(
                          Vector<REAL>& local_data
                          , const std::string & filename
                          , const std::string & groupname
                          )
      {


        // open the file for reading
        hid_t file_id=  H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        assert( file_id > -1 );

        // open the dataset
        hid_t dataset_id = H5Dopen(file_id, groupname.c_str());
        assert( dataset_id > -1 );

        // get the dataspace
        hid_t dataspace_id = H5Dget_space (dataset_id);
        assert(dataspace_id>-1);

        // some needed variables
        herr_t status;
        hsize_t dim,*dims;

        // get the dimension (2-d or 3-d)
        dim = H5Sget_simple_extent_ndims( dataspace_id );

        // get the size of the problem
        dims=(hsize_t*)malloc(dim * sizeof (hsize_t));
        status = H5Sget_simple_extent_dims( dataspace_id , dims, 0);
        if( status < 0 )
          logger << "Warning : H5Sget_simple_extent_dims() failed." << std::endl;
        assert( status > -1 );

        UINT local_size=1;
        hsize_t offset[dim];
        for( UINT i=0; i<dim; i++ )
          {
            local_size*=dims[i];
            offset[i]=0;
          }

        // create the memory space
        hid_t memspace_id = H5Screate_simple (dim, dims, NULL);

        //select the hyperslab-
        status = H5Sselect_hyperslab (memspace_id, H5S_SELECT_SET, offset, NULL,
                                      dims, NULL);
        assert(status>-1);

        //resize the return data
        if( local_data.size() != local_size )
          local_data.resize( local_size );

        /* set up the collective transfer properties list */
        hid_t xferPropList = H5Pcreate (H5P_DATASET_XFER);
        assert( xferPropList > -1 );

        // finally the reading from the file
        status = H5Dread(                     dataset_id
                                              , H5T_NATIVE_DOUBLE //image.DataType
                                              , memspace_id
                                              , dataspace_id // H5S_ALL //
                                              , xferPropList //H5P_DEFAULT //
                                              , &( local_data[0] )
                                              );
        assert(status>-1);

        // close the identifyers
        H5Dclose (dataset_id);
        H5Sclose (dataspace_id);
        H5Sclose (memspace_id);
        free(dims);
        status = H5Fclose( file_id );
        assert( status > -1 );
      };

      /** function to write a vector (parallel) to a hdf5 file
       *
       * \tparam gobal_dim the global dimension of the stored data (total size)
       * \tparam data data which will be written to the file
       * \tparam local_count give the size of the local data
       * \tparam local_offset the offset of the data (in each direction)
       * \tparam helper the DUNE MPIHelper for the MPI communication
       * \tparam groupname is the name/path where the data are stored in the HDF5 file
       * \tparam filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */
      template<typename IDT>
      static void h5_pWrite( const Vector<REAL> &data
                             , const std::string& filename
                             , const std::string& groupname
                             , const IDT& inputdata
                             , const Vector<UINT>& global_dim
                             , const Vector<UINT>& local_offset
                             , const Vector<UINT>& local_count
                             , MPI_Comm communicator
                             )
      {
        Dune::Timer watch;

        logger << "h5_pWrite: " << filename << std::endl;

        //Info varibale need for the HDF5
        MPI_Info mpiInfo = MPI_INFO_NULL;
        herr_t status;             /* Generic return value */

        //get the dimension of the problem -> no cross checking ig the given data make any sense!
        UINT dim = local_count.size();

        // Set up file access property list with parallel I/O access
        hid_t plist_id= H5Pcreate( H5P_FILE_ACCESS );


        // Set up file access property list with parallel I/O access
        status=H5Pset_fapl_mpio(plist_id, communicator, mpiInfo);  //collcetive MPI!!! needs to be called on each processor of the communicator
        if( status < 0 )
          logger << "Warning : H5Pset_fapl_mpio() failed." << std::endl;
        assert(plist_id>-1);

        // Create a new file using default properties.
        hid_t file_id= H5Fcreate(
                                 filename.c_str()    // name of the file to be created
                                 , H5F_ACC_TRUNC          // if you are trying to create a file that exists already, the existing file will be truncated, i.e., all data stored on the original file will be erased
                                 , H5P_DEFAULT            // default file creation property list
                                 , plist_id
                                 );
        assert( file_id > -1 );
        H5Pclose(plist_id);

        // set the global size of the grid into a vector of type hsize_t (needed for HDF5 routines)
        hsize_t global_dim_HDF5[ dim];
        for(UINT i=0; i<dim;i++)
          global_dim_HDF5[dim -i-1]=global_dim[i];

        // set the count and offset in the different dimensions (determine the size of the hyperslab) (in hsize_t format, needed for HDF5 routines)
        hsize_t count[ dim ],offset[ dim ];
        for(UINT i=0;i<dim;i++){
          count[dim-i-1]=local_count[i];
          offset[dim-i-1]=local_offset[i];
        }

        //define the total size of the local data
        hsize_t nAllLocalCells = count[0] * count[1];
#ifdef DIMENSION3
        nAllLocalCells *= count[2];
#endif

        //  Create the dataspace for the dataset.
        hid_t filespace = H5Screate_simple(dim, global_dim_HDF5, NULL);
        assert(filespace>-1);

        // Create the dataset with default properties and close filespace.
        hid_t dset_id = H5Dcreate( file_id,
                                   groupname.c_str(),
                                   HDF5_DATA_TYPE,    // Datatype identifier, here: IEEE floating point 32-bit or 64-bit, see my_macros.hh
                                   filespace,
                                   H5P_DEFAULT );
        H5Sclose(filespace);
        assert(dset_id>-1);

        //get the memoryspace (but only if something needs to be written on this processor!)
        hid_t memspace_id;
        if(nAllLocalCells!=0){  // -> otherwise HDF5 warning, because of writing nothing!
          memspace_id = H5Screate_simple(dim, count, NULL);
          assert(memspace_id>-1);
        }

        // Select hyperslab in the file.
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
        //logger<< "h5_pWrite:  hyperslab selected!" << std::endl;

        // Create property list for collective dataset write.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // all MPI processors of this communicator need to call this (if they write or not)

        // fineally write the data to the disk
        // even if nothing should be written H5Dwrite needs to be called!!
        if(nAllLocalCells!=0){ // -> otherwise HDF5 warning, because of writing nothing!
          status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace,
                            plist_id, &( data[0] )  );
          assert(status>-1);
        }else{ // IMPORTANT. otherwise the H5Dwrite() blocks!!!
          status = H5Dwrite(dset_id,
                            H5T_NATIVE_DOUBLE,
                            H5S_ALL,
                            filespace,
                            plist_id,
                            &( data[0] )  );
          assert(status>-1);
        }

        // Close the property list;
        status=H5Pclose(plist_id);
        assert( status > -1 );

        // Close the filespace;
        status=H5Sclose(filespace);
        assert( status > -1 );

        //if something writen close the memspace
        if(nAllLocalCells!=0){
          // Close the mem space;
          status=H5Sclose(memspace_id);
          assert( status > -1 );
        }

        // Close the dataset;
        status=H5Dclose(dset_id);
        assert( status > -1 );

        // Close the file.
        status = H5Fclose( file_id );
        assert( status > -1 );

        //propably not needed. because the H5Dwrite blocks anyway!!
        MPI_Barrier( communicator );


        std::stringstream jobtitle;
        jobtitle << "h5_pWrite: writing "
                 << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   communicator,
                                   inputdata.verbosity,
                                   "IO",
                                   jobtitle.str() );

      };



      /** function to retrieve a vector (sequential) from the hdf5 file
       *
       * \tparam local_data is a 'return value' which gets the data that belongs to the current processor (current hyperslab)
       * \tparam local_count give the size of the local data
       * \tparam local_offset the offset of the data (in each direction)
       *  \tparam helper the DUNE MPIHelper for the MPI communication
       * \tparam groupname is the name/path where the data are stored in the HDF5 file
       * \tparam filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */
      template<typename IDT>
      static void h5_pRead( Vector<REAL>& local_data
                            , const std::string& filename
                            , const std::string& groupname
                            , const IDT& inputdata
                            , const Vector<UINT>& local_offset
                            , const Vector<UINT>& local_count
                            , MPI_Comm communicator
                            )
      {
        Dune::Timer watch;

        logger << "h5_pRead: " << filename << std::endl;

        /* setup file access template with parallel IO access. */
        hid_t access_pList = H5Pcreate( H5P_FILE_ACCESS );
        assert( access_pList > -1 );

        herr_t status;             /* Generic return value */
        MPI_Info mpiInfo = MPI_INFO_NULL;
        status = H5Pset_fapl_mpio( access_pList, communicator, mpiInfo ); // Take care, this requires the MPI version of hdf5!
        if( status < 0 )
          logger << "Warning : H5Pset_fapl_mpio() failed." << std::endl;
        assert( status > -1 );

        // open the file for reading
        hid_t file_id=  H5Fopen (filename.c_str(), H5F_ACC_RDONLY, access_pList);
        assert( file_id > -1 );

        /* Release file-access template */
        status = H5Pclose( access_pList );
        assert( status > -1 );

        // open the dataset
        hid_t dataset_id = H5Dopen(file_id, groupname.c_str());
        assert( dataset_id > -1 );

        // get the dataspace
        hid_t dataspace_id=H5Dget_space (dataset_id);
        assert(dataspace_id>-1);

        // some needed variables
        hsize_t dim,*dims;

        // get the dimension (2-d or 3-d)
        dim = H5Sget_simple_extent_ndims( dataspace_id );

        // get the size of the data structure
        dims=(hsize_t*)malloc(dim * sizeof (hsize_t));
        status = H5Sget_simple_extent_dims( dataspace_id , dims, 0);
        assert( status > -1 );

        //set the local, offset, and count as hsize_t, which is needed by the HDF5 routines
        hsize_t local_size=1;
        hsize_t offset[dim],count[dim];
        for( UINT i=0; i<dim; i++ )
          {
            local_size*=local_count[i];
            offset[dim-i-1]=local_offset[i];
            count[dim-i-1]=local_count[i];
          }

        // create the memory space, if something needes to be read on this processor
        hid_t memspace_id = -1; // just an initialization
        if(local_size!=0)
          memspace_id = H5Screate_simple (1, &local_size, NULL);

        //select the hyperslab
        status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL,
                                      count, NULL);
        assert(status>-1);

        //resize the return data
        if( local_data.size() != local_size )
          local_data.resize( local_size );

        // set up the collective transfer properties list
        hid_t xferPropList = H5Pcreate (H5P_DATASET_XFER);
        assert( xferPropList > -1 );
        //logger << "xferPropList created." << std::endl;

        // finally the reading from the file, only if something needes to be written
        if(local_size!=0){
          status = H5Dread(
                           dataset_id
                           , H5T_NATIVE_DOUBLE //image.DataType
                           , memspace_id
                           , dataspace_id // H5S_ALL //
                           , xferPropList //H5P_DEFAULT //
                           , &( local_data[0] )
                           );
          assert(status>-1);
        }

        // close the identifyers
        H5Dclose (dataset_id);
        H5Sclose (dataspace_id);
        if(local_size!=0) //this identifier only exists if somethings needs to be read
          H5Sclose (memspace_id);
        free(dims);
        status = H5Fclose( file_id );
        assert( status > -1 );


        std::stringstream jobtitle;
        jobtitle << "h5_pRead(): reading "
                 << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   communicator,
                                   inputdata.verbosity,
                                   "IO",
                                   jobtitle.str() );

      };

      /** function to write a vector (parallel) to an EXISTING hdf5 file
       *
       * \tparam gobal_dim the global dimension of the stored data (total size)
       * \tparam data data which will be written to the file
       * \tparam local_count give the size of the local data
       * \tparam local_offset the offset of the data (in each direction)
       * \tparam helper the DUNE MPIHelper for the MPI communication
       * \tparam groupname is the name/path where the data are stored in the HDF5 file
       * \tparam filename is the filename of the data file
       *
       *
       See
       http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html
       for a documentation of the HDF5 API.
       *
       */

      static void h5_pAppend(
                             const Vector<REAL> &data
                             , const std::string& filename
                             , const std::string& groupname
                             , const Vector<UINT>& local_offset
                             , const Vector<UINT>& local_count
                             , MPI_Comm communicator
                             )
      {

        // setup file access template with parallel IO access.
        hid_t access_pList = H5Pcreate( H5P_FILE_ACCESS );
        assert( access_pList > -1 );

        herr_t status;             // Generic return value
        MPI_Info mpiInfo = MPI_INFO_NULL;
        status = H5Pset_fapl_mpio( access_pList, communicator, mpiInfo ); // Take care, this requires the MPI version of hdf5!
        if( status < 0 )
          logger << "Warning : H5Pset_fapl_mpio() failed." << std::endl;
        assert( status > -1 );

        // open the file for reading
        hid_t file_id=  H5Fopen (filename.c_str(), H5F_ACC_RDWR, access_pList);
        assert( file_id > -1 );
        H5Pclose(access_pList);

        // open the dataset
        hid_t dataset_id = H5Dopen(file_id, groupname.c_str());
        assert( dataset_id > -1 );

        // get the dataspace
        hid_t dataspace_id=H5Dget_space (dataset_id);
        assert(dataspace_id>-1);


        // some needed variables
        hsize_t dim,*dims;

        // get the dimension (2-d or 3-d)
        dim = H5Sget_simple_extent_ndims( dataspace_id );

        // get the size of the data structure
        dims=(hsize_t*)malloc(dim * sizeof (hsize_t));
        status = H5Sget_simple_extent_dims( dataspace_id , dims, 0);
        assert( status > -1 );


        // set the count and offset in the different dimensions (determine the size of the hyperslab) (in hsize_t format, needed for HDF5 routines)
        hsize_t count[ dim ],offset[ dim ];
        for(UINT i=0;i<dim;i++){
          count[dim-i-1]=local_count[i];
          offset[dim-i-1]=local_offset[i];
        }


        //define the total size of the local data
        hsize_t nAllLocalCells = count[0] * count[1];
#ifdef DIMENSION3
        nAllLocalCells *= count[2];
#endif

        //get the memoryspace (but only if something needs to be written on this processor!)
        hid_t memspace_id;
        if(nAllLocalCells!=0){  // -> otherwise HDF5 warning, because of writing nothing!
          memspace_id = H5Screate_simple(dim, count, NULL);
          assert(memspace_id>-1);
        }

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

        access_pList = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(access_pList, H5FD_MPIO_COLLECTIVE); // all MPI processors of this communicator need to call this (if they write or not)

        // fineally write the data to the disk
        // even if nothing should be written H5Dwrite needs to be called!!
        if(nAllLocalCells!=0){ // -> otherwise HDF5 warning, because of writing nothing!
          status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
                             access_pList, &( data[0] )  );
          assert(status>-1);
        }else{ // IMPORTANT. otherwise the H5Dwrite() blocks!!!
          status = H5Dwrite( dataset_id,
                             H5T_NATIVE_DOUBLE,
                             H5S_ALL,
                             dataspace_id,
                             access_pList,
                             &( data[0] )  );
          assert(status>-1);
        }

        // close the identifyers
        H5Dclose (dataset_id);
        H5Sclose (dataspace_id);
        // Release file-access template
        status = H5Pclose( access_pList );
        assert( status > -1 );
        if(nAllLocalCells!=0) //this identifier only exists if somethings needs to be read
          H5Sclose (memspace_id);
        free(dims);
        status = H5Fclose( file_id );
        assert( status > -1 );

        return;
      };





      // ==========
      // Take care:
      // ==========
      // write_vector_to_HDF5( ) is used to store a DOF vector to a one-dimensional HDF5 file.
      // For each process, it reserves a section which has the size of the local number of degrees of freedom.
      // This implementation has the advantage that it works indepently of the used FEM basis.
      // Disadvantage:
      // Since a one-dimensional HDF5 file structure contains no information about the 3D structure of the
      // solution vector, you need to make sure, that the same number of processes and the same grid partitioning
      // will be available when you retrieve the DOF vector later using read_vector_from_HDF5( ).
      //
      // counterpart:
      // void read_vector_from_HDF5( )
      //
      template<typename NUM,typename GFS>
      static void write_vector_to_HDF5( std::vector<NUM>& vData,
                                        const GFS &gfs,
                                        const std::string & filename,
                                        const std::string & datasetname
                                        ){

        Dune::Timer watch;
        logger << "write_vector_to_HDF5: " << filename << std::endl;

        const typename GFS::Traits::GridViewType& gv = gfs.gridView();

        UINT local_dofs = gfs.size();
        UINT total_dofs = gfs.size();
        total_dofs = gv.comm().sum( total_dofs );

        UINT global_dofs = gfs.size();
        global_dofs = gv.comm().max( global_dofs );

#ifdef DEBUG_LOG
        logger << "DEBUG: write_vector_to_HDF5: local_dofs = " << local_dofs << std::endl;
        logger << "DEBUG: write_vector_to_HDF5: total_dofs = " << total_dofs << std::endl;
#endif
        int mpi_size = gv.comm().size();
        int mpi_rank = gv.comm().rank();
        std::vector<int> v1(mpi_size,0);
        v1[ mpi_rank ] = local_dofs;          // This means: v1[i] = 0 for i!=mpi_rank.

        std::vector<int> v2(mpi_size,0);

        MPI_Allreduce( &(v1[0]),  // v1 contains local_dofs only of the current process.
                       &(v2[0]),  // v2 will contain info about local_dofs of ALL processes.
                       mpi_size,
                       MPI_INT,
                       MPI_SUM,
                       gv.comm() );


        // The whole DOF vector vData will be stored as a one-dimensional HDF5 file.
        // Each process will occupy a section of the hdf5 file.
        // Each section starts at the position "my_offset" and has the length "local_dofs".

        int my_offset = 0;
        for(int i=0;i<mpi_rank;i++){
          my_offset += v2[i];
        }

#ifdef DEBUG_LOG
        logger << "DEBUG: write_vector_to_HDF5: my_offset = " << my_offset << std::endl;
#endif
        const int h5rank = 2;
        hid_t plist_id;          /* property list identifier */
        hid_t file_id;
        hid_t dset_id;           /* file and dataset identifiers */
        hid_t filespace_id;
        hid_t memspace_id;       /* file and memory dataspace identifiers */
        hsize_t dimsf[h5rank];   /* dataset dimensions */
        hsize_t	count[h5rank];	 /* hyperslab selection parameters */
        hsize_t	offset[h5rank];
        //double *data;            /* pointer to data buffer to write */
        herr_t status;


        //Info varibale need for the HDF5
        MPI_Info mpiInfo = MPI_INFO_NULL;

        //some needed variables and typdefs
        //const UINT dim = GV::dimension;  // dimensionality of the problem (2-d or 3-d)

        /*
         * Set up file access property list with parallel I/O access
         */
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        status = H5Pset_fapl_mpio(plist_id,gv.comm(),mpiInfo);
        if( status < 0 )
          logger << "Warning : H5Pset_fapl_mpio() failed." << std::endl;


        /*
         * Create a new file collectively and release property list identifier.
         */
        file_id = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id );
        H5Pclose(plist_id);

        /*
         * Create the dataspace for the dataset.
         */
        dimsf[0] = total_dofs;
        dimsf[1] = 1;
        filespace_id = H5Screate_simple(h5rank,dimsf,NULL);
        /*
         * Create the dataset with default properties and close filespace.
         */
        dset_id = H5Dcreate(file_id,
                            datasetname.c_str(),
                            H5T_NATIVE_DOUBLE,
                            filespace_id,
                            H5P_DEFAULT);
        H5Sclose(filespace_id);
        /*
         * Each process defines dataset in memory and writes it to the hyperslab
         * in the file.
         */
        count[0] = local_dofs;
        count[1] = 1;
        offset[0] = my_offset;
        offset[1] = 0;
        memspace_id = H5Screate_simple(h5rank,count,NULL);
        /*
         * Select hyperslab in the file.
         */
        filespace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace_id,H5S_SELECT_SET,offset,NULL,count,NULL);
        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite( dset_id,
                           H5T_NATIVE_DOUBLE,
                           memspace_id,
                           filespace_id,
                           plist_id,
                           &(vData[0]) );
        assert(status>-1);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace_id);
        H5Sclose(memspace_id);
        H5Pclose(plist_id);
        H5Fclose(file_id);

        std::stringstream jobtitle;
        jobtitle << "write_vector_to_HDF5: writing "
                 << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   General::verbosity,
                                   "IO",
                                   jobtitle.str() );

      }





      // Save BackendVector to hdf5 (This implementation works not only for Q1 finite elements!)
      // counterparts:
      // void read_BackendVector_from_HDF5()
      template<typename GFS>
      static void write_BackendVector_to_HDF5( const GFS& gfs,
                                               const CInputData& inputdata,
                                               const std::string & filename,
                                               const std::string & groupname,
                                               const typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type& backend_vector,
                                               bool preserve_structure=true,
                                               int current_grid_level=0
                                               ){
        Dune::Timer watch;
        logger << "write_BackendVector_to_HDF5: " << filename << std::endl;


        Vector<REAL> standard_vector;
        General::copy_to_std( backend_vector, standard_vector );
        //backend_vector.std_copy_to( standard_vector );
        //for( auto it = backend_vector.begin(); it!=backend_vector.end(); ++it )
        //  standard_vector.push_back(*it);

        std::stringstream jobtitle;
        jobtitle << "write_BackendVector_to_HDF5: copying data";
        General::log_elapsed_time( watch.elapsed(),
                                   gfs.gridView().comm(),
                                   inputdata.verbosity,
                                   "IO",
                                   jobtitle.str() );

        if( !preserve_structure ) {
          write_vector_to_HDF5( standard_vector,
                                gfs,
                                filename,
                                groupname
                                );
        }
        else {

          typedef typename GFS::Traits::GridViewType GV;
          enum{ dim = GV::dimension };
          //typedef typename GFS::Traits::BackendType VBE;

          //enum{blocksize = VBE::BlockSize};
          int blocksize = int( backend_vector.N() );

          // P0: blocksize=1 for dim=2 or dim=3
          // Q1: blocksize=4 for dim=2
          // Q1: blocksize=8 for dim=3

          int blocksizePerDimension = 1;       // P0 or CG
          FEMType::Type femtype=FEMType::CG;

          if( gfs.size() == (UINT)gfs.gridView().size(0) ){
            // This is the CCFV case.
            femtype=FEMType::DG;
          }

          if( blocksize == 4 || blocksize == 8 ){
            // DG(1):
            blocksizePerDimension = 2;
            femtype=FEMType::DG;
          }

          if( blocksize == 9 || blocksize == 27 ){
            // DG(2):
            blocksizePerDimension = 3;
            femtype=FEMType::DG;
          }

          if( blocksize == 16 || blocksize == 64 ){
            // DG(3):
            blocksizePerDimension = 4;
            femtype=FEMType::DG;
          }

          //for(int i=0;i<standard_vector.size();i++){
          //  logger << "DEBUG: standard_vector["<<i<<"] = " << standard_vector[i] << std::endl;
          //}

          Vector<REAL> gridsizes (inputdata.domain_data.virtual_gridsizes);
          if( current_grid_level > 0 )
            gridsizes /= REAL( std::pow(2,current_grid_level) );

          h5g_pWrite( gfs.gridView(),
                      standard_vector,
                      filename,
                      groupname,
                      inputdata,
                      inputdata.domain_data.nCells,
                      blocksizePerDimension,
                      femtype,
                      current_grid_level,
                      preserve_structure
                      );
        }

      }
      // end of void write_BackendVector_to_HDF5








      // ==========
      // Take care:
      // ==========
      // read_vector_from_HDF5( ) is used to read a DOF vector from a one-dimensional HDF5 file that was previously stored using the function
      // write_vector_to_HDF5( ). If the number of processes or the grid partioning has changed since then, you will
      // get a mess!
      //
      // void write_vector_to_HDF5( )
      template<typename NUM,typename GFS>
      static void read_vector_from_HDF5( std::vector<NUM>& vData,
                                         const GFS& gfs,
                                         const std::string & filename,
                                         const std::string & datasetname
                                         ){

        Dune::Timer watch;
        logger << "read_vector_from_HDF5: " << filename << std::endl;

        const typename GFS::Traits::GridViewType& gv = gfs.gridView();

        UINT local_dofs = gfs.size();
        UINT total_dofs = gfs.size();
        total_dofs = gv.comm().sum( total_dofs );

        UINT global_dofs = gfs.size();
        global_dofs = gv.comm().max( global_dofs );

#ifdef DEBUG_LOG
        logger << "DEBUG: read_vector_from_HDF5: local_dofs = " << local_dofs << std::endl;
        logger << "DEBUG: read_vector_from_HDF5: total_dofs = " << total_dofs << std::endl;
#endif

        int mpi_size = gv.comm().size();
        int mpi_rank = gv.comm().rank();
        std::vector<int> v1(mpi_size,0);
        v1[ mpi_rank ] = local_dofs;

        std::vector<int> v2(mpi_size,0);

        MPI_Allreduce( &(v1[0]),
                       &(v2[0]),
                       mpi_size,
                       MPI_INT,
                       MPI_SUM,
                       gv.comm() );

        // The whole DOF vector vData was stored as a one-dimensional HDF5 file.
        // There were np sections in the hdf5 file.
        // Each section started at the position "my_offset" and had the length "local_dofs".

        int my_offset = 0;
        for(int i=0;i<mpi_rank;i++){
          my_offset += v2[i];
        }

#ifdef DEBUG_LOG
        logger << "DEBUG: my_offset = " << my_offset << std::endl;
#endif

        const int h5rank = 2;
        hid_t plist_id;          /* property list identifier */
        hid_t file_id;
        hid_t dset_id;           /* file and dataset identifiers */
        //hid_t filespace_id;
        hid_t memspace_id;       /* file and memory dataspace identifiers */
        //hsize_t dimsf[h5rank];   /* dataset dimensions */
        hsize_t	count[h5rank];	 /* hyperslab selection parameters */
        hsize_t	offset[h5rank];
        //double *data;            /* pointer to data buffer to write */
        herr_t status;


        //Info varibale need for the HDF5
        MPI_Info mpiInfo = MPI_INFO_NULL;

        /*
         * Set up file access property list with parallel I/O access
         */
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        status = H5Pset_fapl_mpio(plist_id,gv.comm(),mpiInfo);
        if( status < 0 )
          logger << "Warning : H5Pset_fapl_mpio() failed." << std::endl;

        /* open the file collectively */
        file_id = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, plist_id );

        H5Pclose(plist_id);

        dset_id = H5Dopen( file_id, datasetname.c_str() );

        hid_t dspace_id = H5Dget_space( dset_id );

        hsize_t dims[ h5rank ];
        status = H5Sget_simple_extent_dims( dspace_id, dims, 0 );


        /*
         * Each process defines dataset in memory and reads it from the hyperslab
         */
        count[0] = local_dofs;
        count[1] = 1;
        offset[0] = my_offset;
        offset[1] = 0;

        memspace_id = H5Screate_simple( h5rank, count, NULL );

        status = H5Sselect_hyperslab( dspace_id
                                      , H5S_SELECT_SET
                                      , offset
                                      , NULL
                                      , count
                                      , NULL
                                      );

        if( vData.size() != local_dofs )
          vData.resize( local_dofs );
        status = H5Dread( dset_id
                          , H5T_NATIVE_DOUBLE
                          , memspace_id
                          , dspace_id
                          , H5P_DEFAULT
                          , &( vData[0] )
                          );

#ifdef DEBUG_LOG
        logger << "DEBUG: read_Vector_From_HDF5: H5read status = " << status << std::endl;
#endif
        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(dspace_id);
        H5Sclose(memspace_id);
        H5Fclose(file_id);

        //vData.resize( gfs.size() ); // cut vector back to original size
        std::stringstream jobtitle;
        jobtitle << "read_vector_from_HDF5: reading "
                 << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   General::verbosity,
                                   "IO",
                                   jobtitle.str() );

      }





      // Read BackendVector from hdf5 (This implementation works only for Q1 finite elements!)
      // <---> write_BackendVector_to_HDF5()
      template<typename GFS>
      static void read_BackendVector_from_HDF5( const GFS& gfs,
                                                const CInputData& inputdata,
                                                const std::string & filename,
                                                const std::string & groupname,
                                                typename Dune::PDELab::BackendVectorSelector<GFS,REAL>::Type& backend_vector,
                                                bool preserve_structure=true,
                                                int current_grid_level=0
                                                ) {

        Dune::Timer watch;
        logger << "read_BackendVector_from_HDF5: " << filename << std::endl;

        Vector<REAL> read_local_data;

        if( !preserve_structure ) {
          read_vector_from_HDF5( read_local_data, gfs, filename, groupname );
        }
        else {

          typedef typename GFS::Traits::GridViewType GV;
          enum{ dim = GV::dimension };
          //typedef typename GFS::Traits::BackendType VBE;

          int blocksize = int( backend_vector.N() );
          // enum{blocksize = VBE::BlockSize};
          // P0: blocksize=1 for dim=2 or dim=3
          // Q1: blocksize=4 for dim=2
          // Q1: blocksize=8 for dim=3

          int blocksizePerDimension = 1;       // P0 or CG
          FEMType::Type femtype=FEMType::CG;

          if( gfs.size() == (UINT)gfs.gridView().size(0) ){
            femtype=FEMType::DG;
          }

          if( blocksize == 4 || blocksize == 8 ){
            // DG(1):
            blocksizePerDimension = 2;
            femtype=FEMType::DG;
          }

          if( blocksize == 9 || blocksize == 27 ){
            // DG(2):
            blocksizePerDimension = 3;
            femtype=FEMType::DG;
          }

          if( blocksize == 16 || blocksize == 64 ){
            // DG(3):
            blocksizePerDimension = 4;
            femtype=FEMType::DG;
          }


          Vector<UINT> local_count;
          Vector<UINT> local_offset;

          h5g_pRead( gfs.gridView(),
                     read_local_data,
                     filename,
                     groupname,
                     local_count,
                     local_offset,
                     inputdata,
                     blocksizePerDimension,
                     femtype,
                     current_grid_level
                     );
        }

        watch.reset();

        UINT localdatasize = read_local_data.size();
        if( backend_vector.base().size() != localdatasize )
          backend_vector.base().resize( localdatasize );

        for( int i=0; i<localdatasize; ++i )
          backend_vector.base()[i] = read_local_data[i];

        //logger << "Read backend-vector from " << filename.c_str() << std::endl;
        //logger << "standard vector size = " << read_local_data.size() << std::endl;

        //for(int i=0;i<read_local_data.size();i++){
        //  logger << "DEBUG: read_local_data["<<i<<"] = " << read_local_data[i] << std::endl;
        //}

        std::stringstream jobtitle;
        jobtitle << "read_BackendVector_from_HDF5(): copying data";
        General::log_elapsed_time( watch.elapsed(),
                                   gfs.gridView().comm(),
                                   inputdata.verbosity,
                                   "IO",
                                   jobtitle.str() );

      } // end of void read_BackendVector_from_HDF5


      template< typename GV,
                typename PGV,
                typename IDT,
                typename V1,
                typename V2
                >
      static void write_to_HDF5_on_root( const GV& gv,
                                         const PGV pRootGridView,
                                         const IDT& inputdata,
                                         const V1& iData, // datavector on current process
                                         const V2& iGIdx, // global grid index on current process
                                         const std::string& filename,
                                         const std::string& groupname
                                         ){

        Dune::Timer watch;
        logger << "write_to_HDF5_on_root: " << filename << std::endl;

        int mpi_size = gv.comm().size();
        int mpi_rank = gv.comm().rank();
        UINT local_dofs = iData.size();
        UINT global_dofs = gv.comm().sum( local_dofs );

        /*
          std::cout << "QQ" << mpi_rank
          << " local_dofs = " << local_dofs
          << " global_dofs = " << global_dofs
          << std::endl;
        */

        V1 gData( global_dofs, 0 ); // datavector on all processors
        V2 gGIdx( global_dofs, 0 ); // global grid index on all processors

        V1 gData2( global_dofs, 0 ); // datavector on all processors
        V2 gGIdx2( global_dofs, 0 ); // global grid index on all processors

        std::vector<int> v1(mpi_size,0);
        v1[ mpi_rank ] = local_dofs;          // This means: v1[i] = 0 for i!=mpi_rank.

        std::vector<int> v2(mpi_size,0);

        MPI_Allreduce( &(v1[0]),  // v1 contains local_dofs only of the current process.
                       &(v2[0]),  // v2 will contain info about local_dofs of ALL processes.
                       mpi_size,
                       MPI_INT,
                       MPI_SUM,
                       gv.comm() );


        // The whole DOF vector vData will be stored as a one-dimensional HDF5 file.
        // Each process will occupy a section of the hdf5 file.
        // Each section starts at the position "my_offset" and has the length "local_dofs".

        int my_offset = 0;
        for(int i=0;i<mpi_rank;i++){
          my_offset += v2[i];
        }

        /*
          std::cout << "QQ" << mpi_rank
          << " my_offset: " << my_offset
          << std::endl;
        */

        for(int j=0;j<local_dofs;j++){
          gData[my_offset+j] = iData[j];
          /*
            std::cout << "QQ" << mpi_rank
            << " j: " << j
            << " iData[j]: " << iData[j]
            << std::endl;
          */
          gGIdx[my_offset+j] = iGIdx[j];

          /*
            std::cout << "QQ" << mpi_rank
            << " >> iGIdx[" << j << "]: " << iGIdx[j]
            << std::endl;
          */
        }

        MPI_Reduce( &(gData[0]),
                    &(gData2[0]),
                    global_dofs,
                    MPI_DOUBLE,
                    MPI_SUM,
                    0, // root is the target
                    gv.comm() );

        const int dim = GV::dimension;
        MPI_Reduce( &(gGIdx[0]),
                    &(gGIdx2[0]),
                    dim*global_dofs,
                    MPI_INT,
                    MPI_SUM,
                    0, // root is the target
                    gv.comm() );

        if( mpi_rank == 0 ){

          V1 s0;

          typedef typename GV::Grid::GlobalIdSet::IdType GlobalIdType;
          typedef std::map<GlobalIdType,REAL> ContainerType;
          ContainerType data_container;
          for( int i=0; i<global_dofs; ++i ){
            /*
              std::cout << "QQ" << mpi_rank
              << " >> gGIdx2[" << i << "]: " << gGIdx2[i]
              << std::endl;
            */
            data_container[ gGIdx2[i] ] = gData2[i];
          }
          for( auto eit=pRootGridView->template begin<0>()
                 ; eit!=pRootGridView->template end<0>()
                 ; ++eit) {

            auto idx = pRootGridView->grid().globalIdSet().id(*eit);

            typename ContainerType::const_iterator data_it =
              data_container.find( idx );
            if(data_it != data_container.end())
              s0.push_back( data_it->second );
            else
              logger << "WARNING: data_it not found for element idx = "
                     << idx << std::endl;

            Dune::FieldVector<REAL,dim> center = eit->geometry().center();
            /*
              std::cout << "QQ" << mpi_rank
              << " idx: " << idx
              << " center: " << center
              << std::endl;
            */
          }
          /*
            std::cout << "QQ" << mpi_rank
            << " s0.size(): " << s0.size()
            << std::endl;
          */

          HDF5Tools::h5g_Write( s0
                                , groupname
                                , inputdata
                                , filename
                                );
        }

        std::stringstream jobtitle;
        jobtitle << "write_to_HDF5_on_root(): writing "
                 << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   inputdata.verbosity,
                                   "IO",
                                   jobtitle.str() );

      } // end of void write_to_HDF5_on_root()


    }; // class HDF5Tools


  } // Gesis
} // Dune








#endif	/* DUNE_GESIS_HDF5_TOOLS_HH */
