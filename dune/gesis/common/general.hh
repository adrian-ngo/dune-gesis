/* 
 * File:   general.hh
 * Author: ngo
 *
 * Created on June 25, 2010, 3:46 PM
 * Last Modified on July 29 2014
 */

#ifndef DUNE_GESIS_GENERAL_HH
#define	DUNE_GESIS_GENERAL_HH

#include <iostream>
#include <string>
#include <fstream>

#include <string.h> // needed for strcmp()
#include <stdio.h>
#include <stdlib.h>

#include <sys/stat.h> // for mkdir()

#include <vector>

//#include "boost/filesystem.hpp"

extern CLogfile logger; // declared + initalized in the main function!



namespace Dune {
  namespace Gesis {


    class General{
      
    private:
      General(){};
      
    public:
      
      static int verbosity; // must be initialized down below

      static double elapsed_MARK;
      static double elapsed_CHECK;

      static double elapsed_IO;
      static double elapsed_FFTW;
      static double elapsed_PDE;
      static double elapsed_EVAL;
      static double elapsed_REDIST;


      
      static bool bDirectoryExists( const char *dirname )
      {
        struct stat st;
        if( stat( dirname, &st ) == 0 )
          return true;
        else
          return false;
      };


      static void makedir( const std::string foldername ){

        int stat = 0;
        stat = mkdir( foldername.c_str(), S_IRWXU | S_IRWXG | S_IRWXO );

        if( stat == 0 )
          std::cout << "Folder created: " << foldername.c_str() << std::endl;
        else if( stat == -1 )
          std::cout << "Folder exists: " << foldername.c_str() << std::endl;
        else
          std::cout << "Error: Cannot create folder "
                    << foldername
                    << ". Please check your write permissions in the current working directory!"
                    << std::endl;

      };

      static void createDirectory( const std::string dirname ){
        if( bDirectoryExists( dirname.c_str() ) )
          std::cout << "Directory already existing: " << dirname.c_str() << std::endl;
        else{
          makedir( dirname.c_str() );
          std::cout << "Directory created: " << dirname.c_str() << std::endl;
        }
      };


      static void deleteFile( const std::string filename ){
        if( 0 != remove( filename.c_str() ) ){
          std::cout << "Error deleting file: " << filename << std::endl;
        }
        else{
          std::cout << "File deleted: " << filename << std::endl;
        }
      };



      template<typename BVEC, typename VEC>
      inline static void copy_to_std( const BVEC& backend_vector, VEC& standard_vector ){
        standard_vector.resize(0);
        for( auto it = backend_vector.begin(); it!=backend_vector.end(); ++it )
          standard_vector.push_back(*it);
      }

      template<typename BVEC, typename VEC>
      inline static void copy_from_std( const VEC& standard_vector, BVEC& backend_vector ){
        for( int i=0; i<standard_vector.size(); ++i )
          backend_vector.base()[i] = standard_vector[i];
      }
      

      inline static UINT indexconversion_3d_to_1d( 
                                                  const UINT iz
                                                  , const UINT iy
                                                  , const UINT ix
                                                  , const UINT Nz
                                                  , const UINT Ny
                                                  , const UINT Nx
                                                   )
      {
        // Adrian: I got this formula from "http://www.fftw.org/doc/Dynamic-Arrays-in-C.html#Dynamic-Arrays-in-C"!
        UINT l = ix + Nx * ( iy + Ny*iz );
        return l;
      };


      inline static UINT indexconversion_to_1d( const Vector<UINT>& index,
                                                const Vector<UINT>& nCells )
      {
        // Adrian: I got this formula from "http://www.fftw.org/doc/Dynamic-Arrays-in-C.html#Dynamic-Arrays-in-C"!
#ifdef DIMENSION3
        UINT l = index[0] + nCells[0] * ( index[1] + nCells[1] * index[2] ); // 3d to 1d
#else
        UINT l = index[0] + nCells[0] * index[1]; // 2d to 1d
#endif
        return l;
      };



      // index conversion from 1d-index to 3d-index, opposite of "indexconversion_3d_to_1d"
      inline static void indexconversion_1d_to_3d( 
                                                  UINT i 
                                                  ,UINT& iz
                                                  ,UINT& iy
                                                  ,UINT& ix
                                                  , const UINT Nz
                                                  , const UINT Ny
                                                  , const UINT Nx
                                                   )
      {
        // Adrian: I got this formula from "http://www.fftw.org/doc/Dynamic-Arrays-in-C.html#Dynamic-Arrays-in-C"!
        UINT tmp;
  
        tmp=i%(Nx*Ny);
        iz=(i-tmp)/(Nx*Ny);
  
        i=tmp;
        ix=i%Nx;
        iy=(i-ix)/Nx;
      };



      // index conversion from 2d-index to 1d-index
      inline static UINT indexconversion_2d_to_1d( 
                                                  const UINT iy
                                                  , const UINT ix
                                                  , const UINT Ny
                                                  , const UINT Nx
                                                   )
      {
        UINT l = iy * Nx + ix;
        return l;
      };

      // index conversion from 1d-index to 2d-index, opposite of indexconversion_2d_to_1d
      inline static void indexconversion_1d_to_2d( 
                                                  UINT i
                                                  , UINT& iy
                                                  , UINT& ix
                                                  , const UINT Ny
                                                  , const UINT Nx
                                                   )
      {
        ix=i%Nx;
        iy=(i-ix)/Nx;
      };






      template<typename REAL>
      inline static REAL autocovariancefunction( 
                                         REAL variance, 
                                         std::vector<REAL> x, 
                                         std::vector<REAL> lambda, 
                                         std::string model  ) {

        // Ref. Olaf Cirpka's "Stochastic Methods" Script, p. 25, equation (2.9) and (2.15)-(2.7)

        UINT dim = x.size();
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


      template<typename PROP>
      static bool save_kfield_properties(
                                         const std::string filename
                                         , const UINT dim
                                         , const std::vector<CTYPE>& extensions
                                         , const std::vector<UINT>& nCells
                                         , const PROP& kfield_data )
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
            filestr << " " << extensions[i];
          filestr << std::endl;

          filestr << "nCells= ";
          for(UINT i=0; i<dim; i++)
            filestr << " " << nCells[i];
          filestr << std::endl;
	
          filestr << "vertical_zones= " << kfield_data.nz << std::endl;
	
          for(UINT ii=0;ii<kfield_data.nz;ii++){
        
            filestr << "model= " << kfield_data.zones[ii].model;
            filestr << std::endl;

            filestr << "beta= " << kfield_data.zones[ii].beta;
            filestr << std::endl;

            filestr << "qbb_y= " << kfield_data.zones[ii].qbb_y;
            filestr << std::endl;

            filestr << "variance= " << kfield_data.zones[ii].variance;
            filestr << std::endl;

            filestr << "embedding_factor= " << kfield_data.zones[ii].embedding_factor;
            filestr << std::endl;

            filestr << "correlation_lambda= ";
            for(UINT i=0; i<dim; i++)
              filestr << " " << kfield_data.zones[ii].correlation_lambda[i];
            filestr << std::endl;
          }
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

      template<typename PROP>
      static bool load_kfield_properties(
                                         const std::string filename, 
                                         const UINT dim, 
                                         const std::vector<CTYPE>& extensions, 
                                         const std::vector<UINT>& nCells, 
                                         const PROP& kfield_data,
                                         const int my_rank
                                         )
      {    
        std::ifstream filestr( filename.c_str(), std::ios::in );
        if( filestr.is_open() ) {

          if( verbosity>=VERBOSITY_INVERSION && my_rank==0 )
            std::cout << "Loading " << filename << std::endl;
      
          std::string str;

          UINT read_dim;
          filestr >> str >> read_dim;
          if( dim != read_dim )
            {
              if(my_rank==0)
                logger << "Drop obsolete Y-field with : dim = " << read_dim << std::endl;
              return false;
            }

          CTYPE read_extensions[3];
          if(dim == 3){
            filestr >> str >> read_extensions[0] >> read_extensions[1] >> read_extensions[2];
            if( read_extensions[0]!=extensions[0] || read_extensions[1] != extensions[1] || read_extensions[2] != extensions[2] )
              {
                if( verbosity>=VERBOSITY_INVERSION && my_rank==0 )
                  std::cout << "New 3D random field with Lx = " << extensions[0]
                            << " Ly = " << extensions[1]
                            << " Lz = " << extensions[2]
                            << std::endl;
                return false;
              }
          }
          else
            {
              filestr >> str >> read_extensions[0] >> read_extensions[1];
              if( read_extensions[0]!=extensions[0] || read_extensions[1] != extensions[1] )
                {
                  if( verbosity>=VERBOSITY_INVERSION && my_rank==0 )
                    std::cout << "New 2D random field with Lx = " << extensions[0]
                              << " Ly = " << extensions[1]
                              << std::endl;
                  return false;
                }
            }

          UINT read_nCells[3];
          if(dim == 3){
            filestr >> str >> read_nCells[0] >> read_nCells[1] >> read_nCells[2];
            if( read_nCells[0]!=nCells[0] || read_nCells[1] != nCells[1] || read_nCells[2] != nCells[2] )
              {
                if( verbosity>=VERBOSITY_INVERSION && my_rank==0 )
                  std::cout << "New 3D random field with Nx = " << nCells[0]
                            << " Ny = " << nCells[1]
                            << " Nz = " << nCells[2]
                            << std::endl;
                return false;
              }
          }
          else
            {
              filestr >> str >> read_nCells[0] >> read_nCells[1];
              if( read_nCells[0]!=nCells[0] || read_nCells[1] != nCells[1] )
                {
                  if( verbosity>=VERBOSITY_INVERSION && my_rank==0 )
                    std::cout << "New 2D random field with Nx = " << nCells[0]
                              << " Ny = " << nCells[1]
                              << std::endl;
                  return false;
                }
            }

          UINT vertical_zones;
          filestr >> str >> vertical_zones;
          if(vertical_zones != kfield_data.nz ){
            logger << "Drop obsolete Y-field with : vertical_zones = " << vertical_zones << std::endl;
            return false;
          }
        
          for(UINT ii=0; ii<kfield_data.nz;ii++){

            std::string read_model;
            filestr >> str >> read_model;
            if( kfield_data.zones[ii].model != read_model )
              {
                logger << "Drop obsolete Y-field with : model = " << read_model.c_str() << std::endl;
                return false;
              }

            REAL read_beta;
            filestr >> str >> read_beta;
            if( kfield_data.zones[ii].beta != read_beta )
              {
                logger << "Drop obsolete Y-field with : beta = " << read_beta << std::endl;
                return false;
              }

            REAL read_qbb_y;
            filestr >> str >> read_qbb_y;
            if( kfield_data.zones[ii].qbb_y != read_qbb_y )
              {
                logger << "Drop obsolete Y-field with : qbb_y = " << read_qbb_y << std::endl;
                return false;
              }

            REAL read_variance;
            filestr >> str >> read_variance;
            if( kfield_data.zones[ii].variance != read_variance )
              {
                logger << "Drop obsolete Y-field with : variance = " << read_variance << std::endl;
                return false;
              }

            REAL read_embedding_factor;
            filestr >> str >> read_embedding_factor;
            if( kfield_data.zones[ii].embedding_factor != read_embedding_factor )
              {
                logger << "Drop obsolete Y-field with : embedding_factor = " << read_embedding_factor << std::endl;
                return false;
              }


            REAL read_correlation_lambda[3];
            if(dim == 3){
              filestr >> str >> read_correlation_lambda[0] >> read_correlation_lambda[1] >> read_correlation_lambda[2];
              if(        kfield_data.zones[ii].correlation_lambda[0] != read_correlation_lambda[0]
                         || kfield_data.zones[ii].correlation_lambda[1] != read_correlation_lambda[1]
                         || kfield_data.zones[ii].correlation_lambda[2] != read_correlation_lambda[2] )
                {
                  if( verbosity>=VERBOSITY_INVERSION && my_rank==0 )
                    std::cout << "New 3D random field with lambda_x = " << kfield_data.zones[ii].correlation_lambda[0]
                              << " lambda_y = " << kfield_data.zones[ii].correlation_lambda[1]
                              << " lambda_z = " << kfield_data.zones[ii].correlation_lambda[2]
                              << " inside zone " << ii
                              << std::endl;
                  return false;
                }
            }
            else
              {
                filestr >> str >> read_correlation_lambda[0] >> read_correlation_lambda[1];
                if(        kfield_data.zones[ii].correlation_lambda[0] !=read_correlation_lambda[0]
                           || kfield_data.zones[ii].correlation_lambda[1] != read_correlation_lambda[1] )
                  {
                  if( verbosity>=VERBOSITY_INVERSION && my_rank==0 )
                      std::cout << "New 2D random field with lambda_x = " << kfield_data.zones[ii].correlation_lambda[0]
                                << " lambda_y = " << kfield_data.zones[ii].correlation_lambda[1]
                                << " inside zone " << ii
                                << std::endl;
                    return false;
                  }
              }
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

      static REAL time_start( const std::string displaytext ){

        logger << displaytext << std::endl;
        //std::cout << displaytext << std::endl;
        return clock();

      }


      static REAL time_stop( const REAL t_start, const std::string displaytext="" )
      {

        REAL t_end = clock() - t_start;
        t_end /= CLOCKS_PER_SEC;
        std::cout << displaytext << " took " << t_end << " sec." << std::endl;
        logger << displaytext << " took " << t_end << " sec." << std::endl;

        return t_end;
      }



      // addToTimeChecker() is used to set a time marker
      // General::elapsed_MARK is used to measure the time interval between two markers
      // General::elapsed_CHECK is used to sum up all time intervals
      static void addToTimeChecker( const int rank, Dune::Timer& check, const std::string marker ){
        if(rank == 0){
          double check_interval = check.elapsed();

          std::cout << "Duration:|" 
                    << std::setw(6)
                    << "CHECK" << "|"
                    << std::setw(15)
                    << std::setprecision(4)
                    << check_interval << " sec." 
                    << " for "
                    << marker
                    << std::endl;

          std::cout << "Duration:|" 
                    << std::setw(6)
                    << "MARK" << "|"
                    << std::setw(15)
                    << std::setprecision(4)
                    << General::elapsed_MARK << " sec." 
                    << " for "
                    << marker
                    << std::endl;

          General::elapsed_CHECK += check_interval;
        }
        General::elapsed_MARK = 0; // reset after each time marker
        check.reset();
      }


      // parallel version:
      static void log_elapsed_time( REAL elapsed_time,
                                    MPI_Comm communicator,
                                    int verbosity,
                                    std::string jobtype,
                                    std::string jobtitle 
                                    ){

        if( verbosity >= VERBOSITY_TIMER_SUMMARY ){
          int rank=1;
          int size=1;
          MPI_Comm_rank( communicator, &rank );
          MPI_Comm_size( communicator, &size );

          double res;
          if( size>1 )
            MPI_Allreduce( &elapsed_time, &res, 1, MPI_DOUBLE, MPI_MAX, communicator );

          if( rank==0 && verbosity >= VERBOSITY_TIMER_DETAILS ){
            std::cout << "Duration:|"
                      << std::setw(6)
                      << jobtype << "|"
                      << std::setw(15)
                      << std::setprecision(4)
                      << elapsed_time << " sec." 
                      << " for "
                      << jobtitle
                      << std::endl;
          }

          if( rank==0 ){
            if( jobtype == "IO" )
              elapsed_IO += elapsed_time;

            if( jobtype == "FFTW" )
              elapsed_FFTW += elapsed_time;

            if( jobtype == "GWE" )
              elapsed_PDE += elapsed_time;
            if( jobtype == "TPE" )
              elapsed_PDE += elapsed_time;
            if( jobtype == "L2Proj" )
              elapsed_PDE += elapsed_time;
            if( jobtype == "RGV" )
              elapsed_PDE += elapsed_time;

            if( jobtype == "EVAL" )
              elapsed_EVAL += elapsed_time;
            if( jobtype == "REDIST" )
              elapsed_REDIST += elapsed_time;

            elapsed_MARK += elapsed_time;

          }
        }
      }
      
      // sequential version:
      static void log_elapsed_time( REAL elapsed_time,
                                    int verbosity,
                                    std::string jobtype,
                                    std::string jobtitle ){

        if( verbosity >= VERBOSITY_TIMER_DETAILS ){
          std::cout << "Duration:[" 
                    << std::setw(6)
                    << jobtype << "]"
                    << std::setw(15)
                    << std::setprecision(4)
                    << elapsed_time << " sec." 
                    << " for "
                    << jobtitle
                    << std::endl;
        }
        if( verbosity >= VERBOSITY_TIMER_SUMMARY ){
          if( jobtype == "IO" )
            elapsed_IO += elapsed_time;

          if( jobtype == "FFTW" )
            elapsed_FFTW += elapsed_time;

          if( jobtype == "GWE" )
            elapsed_PDE += elapsed_time;
          if( jobtype == "TPE" )
            elapsed_PDE += elapsed_time;
          if( jobtype == "L2Proj" )
            elapsed_PDE += elapsed_time;
          if( jobtype == "RGV" )
            elapsed_PDE += elapsed_time;

          if( jobtype == "EVAL" )
            elapsed_EVAL += elapsed_time;
          if( jobtype == "REDIST" )
            elapsed_REDIST += elapsed_time;

          elapsed_MARK += elapsed_time;
        }
      }



      static REAL reportTotalTime( const int my_rank ){

        REAL totalCountedTime 
          = General::elapsed_PDE 
          + General::elapsed_IO 
          + General::elapsed_FFTW 
          + General::elapsed_EVAL
          + General::elapsed_REDIST
          ;

        if( my_rank==0
            && General::verbosity >= VERBOSITY_TIMER_SUMMARY ){

          std::cout << std::fixed << std::endl;
          std::cout << "TIMER: ========= OVERVIEW OF TIME CONSUMPTION ====================" << std::endl;
          std::cout << "TIMER:" 
                    << std::setw(10) << "IO:" 
                    << std::setw(15) << General::elapsed_IO << " sec."
                    << std::setw(15) << 100.0 * General::elapsed_IO / totalCountedTime << " %"
                    << std::endl;
          std::cout << "TIMER:" 
                    << std::setw(10) << "PDEs:" 
                    << std::setw(15) << General::elapsed_PDE << " sec."
                    << std::setw(15) << 100.0 * General::elapsed_PDE / totalCountedTime << " %"
                    << std::endl;
          std::cout << "TIMER:" 
                    << std::setw(10) << "FFTW:" 
                    << std::setw(15) << General::elapsed_FFTW << " sec."
                    << std::setw(15) << 100.0 * General::elapsed_FFTW / totalCountedTime << " %"
                    << std::endl;
          std::cout << "TIMER:" 
                    << std::setw(10) << "REDIST:" 
                    << std::setw(15) << General::elapsed_REDIST << " sec."
                    << std::setw(15) << 100.0 * General::elapsed_REDIST / totalCountedTime << " %"
                    << std::endl;
          std::cout << "TIMER:" 
                    << std::setw(10) << "EVAL:" 
                    << std::setw(15) << General::elapsed_EVAL << " sec."
                    << std::setw(15) << 100.0 * General::elapsed_EVAL / totalCountedTime << " %"
                    << std::endl;
          std::cout << "TIMER:" 
                    << "===========================================================" << std::endl;
          std::cout << "TIMER:" 
                    << std::setw(10) << "TOTAL:" 
                    << std::setw(15) << totalCountedTime << " sec."
                    << std::endl;
        }


        if( my_rank==0 
            && General::verbosity >= VERBOSITY_TIMER_SUMMARY ){
          std::cout << std::endl;
          std::cout << "Duration CHECK (total) = " << General::elapsed_CHECK << std::endl;
        }

        return totalCountedTime;

      }




      static void createPVDfromVtuList(
                                const std::string& filename
                                , const std::vector< std::string >& vtu_list 
                                )
      {
        logger << "Create PVD out of list of VTUs..." << std::endl;

        std::ofstream file( filename.c_str(), std::ios::out );
        if( file.is_open() )
          {
            file << "<?xml version=\"1.0\"?>" << std::endl;
            file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
            file << "<Collection>" << std::endl;

            //for( std::vector<std::string>::iterator it = vtu_list.begin()
            //	 ; it != vtu_list.end() 
            //	 ; it++ )
            logger << "Collecting VTU-files:  ";
            for( UINT i=0; i<vtu_list.size(); i++ )
              {
                logger << vtu_list[i].c_str() << ".vtu\t";
                file << "<DataSet timestep=\"" << i << "\" part=\"" << i << "\" file=\"" << vtu_list[i].c_str() <<  ".vtu\"/>" << std::endl;
              }
            logger << "\t";
            logger << " -----> " << filename.c_str() << std::endl;


            file << "</Collection>" << std::endl;
            file << "</VTKFile>" << std::endl;
            file << "\n";
            file.close();
          }

      }

      /*
       * This simple function checks if a file exists.
       */
      static bool bFileExists( const std::string& filename )
      {
        struct stat buf;
        if( stat( filename.c_str(), &buf ) != -1 ){

          UINT filesize = (UINT) buf.st_size;
          logger << "File found: " << filename 
                 << ", file size = " << filesize
                 << std::endl;
          if( filesize > 0 ){
            return true;
          }
          else{
            deleteFile( filename );
            return false;
          }
        }
        else{
          return false;
        }
      }



      template<typename T>
      static std::string number2string( const T& x){
        std::stringstream ss;
        ss << x;
        return ss.str();
      }



      static std::string textfile2string( const std::string filename ){
        
        std::ifstream is;
        is.open( filename.c_str(), std::ios::binary );
        
        // get length of file:
        is.seekg (0, std::ios::end);
        long length = is.tellg();
        is.seekg (0, std::ios::beg);
        
        // allocate memory:
        char *buffer = new char [length];
        
        // read data as a block:
        is.read (buffer,length);
        
        std::string contents( buffer );
        delete[] buffer;
        is.close();

        return contents;
      }
      

      static std::string readTextFile2String( const std::string filename, 
                                              const Dune::MPIHelper& helper ) 
      {
        std::istringstream instream; // Declare an input string stream
        // ------------------------------------------------------------------
        // Read text file as a stringstream only on one processor.
        // Then, send this stream to all other processes using 
        // a dynamic character buffer.
        // ------------------------------------------------------------------
        std::string contents;
        int filesize = 0;

        if( helper.rank() == 0 ) {
          logger << "readTextFile2String: Read text file as string stream on P0..." << std::endl;
          
          char *textbuffer;
          std::filebuf *pFileBuffer;
          std::ifstream filestream ( filename.c_str() );
          if( filestream.is_open() ){
            pFileBuffer = filestream.rdbuf();
            // get file size using buffer's members
            filesize=pFileBuffer->pubseekoff (0,std::ios::end,std::ios::in);
            pFileBuffer->pubseekpos (0,std::ios::in);
            
            // allocate memory to contain file data, reserve one more char for ending null
            textbuffer = new char[filesize+1];
            
            // get file data  
            pFileBuffer->sgetn (textbuffer,filesize);

            // Important: append null character to string makes sure that the string ends here!
            textbuffer[filesize] = '\0';
        
            contents = std::string( textbuffer );
            
            filestream.close();
            
          }
          else {
            std::cout << "Warning: Unable to open textfile " << filename.c_str() << std::endl;
            return std::string("");
          }
          
          if( helper.size() > 1 ) {
            // Sender:
            // prepare to broadcast data to other processes
            // Send filesize
            MPI_Bcast( &filesize, 1, MPI_INT, 0, helper.getCommunicator());
            // Send text and do not forget the trailing '\0'
            MPI_Bcast( textbuffer, filesize+1, MPI_CHAR, 0, helper.getCommunicator());
            delete[] textbuffer;
          }

        }

        else {

          // Receive length of character buffer:
          MPI_Bcast( &filesize, 1, MPI_INT, 0, helper.getCommunicator());

          char textbuffer[filesize]; // Do not use char* with new and delete here!
          
          // Receive character buffer with trailing '\0':
          MPI_Bcast( textbuffer, filesize+1, MPI_CHAR, 0, helper.getCommunicator());
          //std::cout.write (textbuffer,filesize);
          
          contents = std::string( textbuffer );
        }
        return contents;
      }







      /*
       * calculate the greatest common divisor (GCD) of integers!
       */
      template<typename A>
      static A GCD(A a, A b){
        if( b==0 )
          return a;
        else
          return GCD( b, a%b );    // recursive
      }






      // This is required only for the DispersionTensor calculation.
      // The porosity is supposed to be dependent on the depth only.
      // for solute TPE:
      template<typename IDT>
      static void zone_parameters( const REAL& zone_coord,
                                   const IDT& inputdata,
                                   REAL& porosity ) {

        // only one zone!
        if(inputdata.yfield_properties.nz<2)
          porosity = inputdata.yfield_properties.zones[0].porosity;

        // at least 2 zones!
        UINT inside=0;

        for(UINT ii=0; ii<inputdata.yfield_properties.nz-1; ii++){
          if(zone_coord<=inputdata.yfield_properties.zones[ii].bottom ){
            porosity = inputdata.yfield_properties.zones[ii].porosity;
            inside=1;
            break;
          }
        }
        //last zone!
        if(zone_coord>inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].bottom && inside==0)
          porosity = inputdata.yfield_properties.zones[inputdata.yfield_properties.nz-1].porosity;

      }

      // for HeatTPE:
      template<typename IDT>
      static void zone_parameters( const REAL& zone_coord,
                                   const IDT& inputdata,
                                   REAL & porosity,
                                   REAL & lambda_s ) {

        UINT inside=0;
        UINT nz = inputdata.yfield_properties.nz;
        // only one zone!
        if(nz<2){
          porosity = inputdata.yfield_properties.zones[0].porosity;
          lambda_s = inputdata.yfield_properties.zones[0].lambda_s;
          return;
        }

        for(UINT ii=0; ii < nz-1; ii++){
          if(zone_coord <= inputdata.yfield_properties.zones[ii].bottom ){
            porosity = inputdata.yfield_properties.zones[ii].porosity;
            lambda_s = inputdata.yfield_properties.zones[ii].lambda_s;
            inside=1;
            break;
          }
        }
        //last zone!
        if(zone_coord > inputdata.yfield_properties.zones[ nz-2 ].bottom && inside==0){
          porosity = inputdata.yfield_properties.zones[ nz-1 ].porosity;
          lambda_s = inputdata.yfield_properties.zones[ nz-1 ].lambda_s;
        }

      }






      // some tools for adaptive mesh refinement

      template< class PowerGFS,
                class UPowerType,
                class GFS,
                class UType
                >
      static void extract_two_indicators( const PowerGFS& powergfs, 
                                   const UPowerType& xpower, 
                                   const GFS& gfs, 
                                   UType & u1, 
                                   UType & u2 
                                   ){
        /*
        for(std::size_t j=0; j<gfs.globalSize(); j++){
          u1.base()[j] = xpower.base()[j];
        }
  
        for(std::size_t j=0; j<gfs.globalSize(); j++){
          u2.base()[j] = xpower.base()[gfs.size()+j];
        }
        */

        for(std::size_t j=0; j<gfs.globalSize(); j++){
          u1.base()[j] = xpower.base()[powergfs.ordering().blockOffset(0)+j];
          u2.base()[j] = xpower.base()[powergfs.ordering().blockOffset(1)+j];
        }

      }


      template<typename IntegerType>
      static void factorize_two_optimally( const IntegerType product, 
                                           IntegerType& factor1, 
                                           IntegerType& factor2 )
      {
        REAL squareroot = sqrt( (REAL) product );
        IntegerType n1 = static_cast<IntegerType>( squareroot );
        IntegerType n2 = 1;
        bool bFound = false;
        for( IntegerType i=n1; i > 1; i-- )
          {
            //std::cout << "i = " << i << std::endl;
            if( (product % i) == 0 )
              {
                n1 = product / i;
                n2 = i;
                bFound = true;

                //std::cout << "n1 = " << n1 << std::endl;
                //std::cout << "n2 = " << n2 << std::endl;

                break;
              }
          }

        if( !bFound )
          {
            factor1 = product;
            factor2 = 1;
          }
        else
          {
            factor1 = n1;
            factor2 = n2;
          }
      }


      template<typename IntegerType>
      static void factorize_three_optimally( const IntegerType product, 
                                             IntegerType& factor1, 
                                             IntegerType& factor2, 
                                             IntegerType& factor3 )
      {
        REAL cuberoot = pow( (REAL) product, 1.0/3.0 );
        IntegerType n1 = static_cast<IntegerType>( cuberoot );
        IntegerType n2 = 1;
        IntegerType n3 = 1;
        bool bFound = false;

        IntegerType nTemp = 1;

        for( IntegerType i=n1; i > 1; i-- )
          {
            if( (product % i) == 0 )
              {
                nTemp = product / i;
                n1 = i;
                bFound = true;
                break;
              }
          }

        if( !bFound )
          {
            factorize_two_optimally( product, n1, n2 );
            factor1 = n1;
            factor2 = n2;
            factor3 = 1;
          }
        else
          {
            factorize_two_optimally( nTemp, n2, n3 );
            factor1 = n1;
            factor2 = n2;
            factor3 = n3;
          }
      }

      /* Check if p is a factor of q */
      static bool isFactorOf( const UINT p, const int q )
      {
        if( (q % p) == 0 )
          {
            return true;
          }
        else
          {
            return false;
          }
      };




      // harmonic averaging of two K-Field values between two cells of the same size:
      template<typename T>
      static T harmonicAverage (T a, T b)
      {
        T eps = 1E-20;
        //return 2.0/(1.0/(a+eps) + 1.0/(b+eps));
        return T(2.0)*a*b / ( a + b + eps );
      }

      //
      // Harmonic averaging of two K-Field values between two rectangular cells 
      // of different sizes:
      // Required for CCFV with locally adaptive mesh refinement 
      // on hanging nodes refinement with cubical cells:
      //
      // Let x1 be the center of C1.
      // Let x2 be the center of C2.
      // Let x1-x2 will cut the big face between C1 and C2 in some point z0.
      // 
      // k1 = K-value of the cell C1
      // k2 = K-value of the neighbouring cell C2
      // l1 = |x1 - z0|
      // l2 = |x2 - z0|
      //
      // Note: 
      // Since we are only interested in the proportion l1/(l1+l2) and l2/(l1+l2),
      // we can as well consider l1 = |x1-f1| and l2 = |x1-f2| where 
      // f1 is the face-center of the face of C1 facing the direction of C2
      // and
      // f2 is the face-center of the face of C2 facing the direction of C1.
      // Or, even easier: We consider the face-lengths instead.
      // 
      template<typename T>
      static T harmonicAverageWeightedByDistance(T k1, T k2, T l1, T l2)
      {
        T eps = T(1E-20);
        T average = (l1+l2)*k1*k2 / ( l1*k2 + l2*k1 + eps );
        /*
          std::cout << " k1=" << k1
          << " k2=" << k2
          << " l1=" << l1
          << " l2=" << l2
          << " average=" << average 
          << std::endl;
        */
        return average;
      }



      /*
        This class is used to get the minimum and the maximum of a 
        solution whose FEM is represented by the Lagrange Basis!
      */
      template<typename VCType, typename COMM>
      static void logMinAndMax( const VCType& vc, 
                                const COMM& communicator,
                                REAL& minimum,
                                REAL& maximum
                                ){
        std::vector<REAL> vc_flat;
        copy_to_std(vc,vc_flat);
        for( UINT i=0; i<vc_flat.size(); i++ ) {
          minimum = std::min( minimum, vc_flat[i] );
          maximum = std::max( maximum, vc_flat[i] );
        }
        logger << "=====> process local maximum/minimum = " 
               << maximum << " / "
               << minimum << std::endl;
        maximum = communicator.max( maximum );
        minimum = communicator.min( minimum );
        if(communicator.rank()==0 && verbosity>=VERBOSITY_EQ_SUMMARY){
          std::cout << "=====> Global maximum/minimum = " 
                    << maximum << " / "
                    << minimum << std::endl;
        }
      }


      /*
        This class is used to get the minimum and the maximum of a 
        solution whose FEM is represented by the Lagrange Basis!
      */
      template<typename VCType, typename COMM>
      static void logMinAndMax( const VCType& vc, 
                                const COMM& communicator
                                ){
        REAL minimum = 1E+100;
        REAL maximum = -1E+100;
        std::vector<REAL> vc_flat;
        copy_to_std(vc,vc_flat);
        for( UINT i=0; i<vc_flat.size(); i++ ) {
          minimum = std::min( minimum, vc_flat[i] );
          maximum = std::max( maximum, vc_flat[i] );
        }
        logger << "=====> process local maximum/minimum = " 
               << maximum << " / "
               << minimum << std::endl;
        maximum = communicator.max( maximum );
        minimum = communicator.min( minimum );
        if(communicator.rank()==0 && verbosity>=VERBOSITY_EQ_SUMMARY){
          std::cout << "=====> Global maximum/minimum = " 
                    << maximum << " / "
                    << minimum << std::endl;
        }
      }


      //
      // Multiply a solution vector with the zone-wise porosity
      //
      template<typename VC_Type,typename GV,typename IDT>
      static void vc_times_theta( VC_Type& out, 
                                  const GV& gv,
                                  const IDT& inputdata )
      {
        
        // types and constants
        int id = 0;
        const int dim = GV::dimension;
        const int dim_minus1=dim-1;
        typedef typename GV::template Codim < dim > ::Iterator ElementLeafIterator;

        const typename GV::IndexSet& is = gv.indexSet();
        
        Dune::FieldVector<REAL,dim> global(0.0);
        
        UINT nzones = inputdata.yfield_properties.nz;
        UINT zone;
        
        if(nzones>1){
          // loop over the grid
          for (ElementLeafIterator it = gv.template begin < dim > (); it != gv.template end < dim > (); ++it) {
    
            id = is.index(*it);
    
            global=it->geometry().center();
    
            zone=0;
            while(zone<nzones-1){
              if(global[dim_minus1]<inputdata.yfield_properties.zones[zone].bottom){
                break;
              }else{
                zone+=1; 
              }
            }
            out.base()[id]*=inputdata.yfield_properties.zones[zone].porosity;
    
          }
        }else{
          out*= inputdata.yfield_properties.zones[0].porosity;
        }

      }






      /*
       * calculate vector container times theta*RT(retardation factor of heat)
       */
      template<typename VC_Type,typename GV,typename IDT>
      static void vc_times_thetaRT( VC_Type& out, 
                                    const GV& gv,
                                    const IDT& inputdata )
      {
        // types and constants
        int id = 0;
        const int dim = GV::dimension;
        const int dim_minus1 = dim-1;
        typedef typename GV::template Codim<dim>::Iterator ElementLeafIterator;
        
        const typename GV::IndexSet& is = gv.indexSet();
        
        Dune::FieldVector<REAL,dim> global(0.0);
        
        UINT nzones = inputdata.yfield_properties.nz;
        UINT zone;

        REAL porosity;
        REAL rho_s;
        REAL c_s;
        REAL rho_m_c_m;

        REAL rho_w = inputdata.transport_parameters.heat_rho_w;
        REAL c_w   = inputdata.transport_parameters.heat_c_w;

        if(nzones>1){
          // loop over the grid
          for (ElementLeafIterator it = gv.template begin<dim> ()
                 ; it != gv.template end<dim> (); ++it) {
            
            id = is.index(*it);
    
            global=it->geometry().center();
    
            zone=0;
            while(zone<nzones-1){
              if(global[dim_minus1] < inputdata.yfield_properties.zones[zone].bottom){
                break;
              }else{
                zone+=1; 
              }
            }
            porosity  = inputdata.yfield_properties.zones[zone].porosity;
            rho_s     = inputdata.yfield_properties.zones[zone].rho;
            c_s       = inputdata.yfield_properties.zones[zone].c_s;
            rho_m_c_m = porosity*rho_w*c_w + (1.0-porosity)*rho_s*c_s;
      
            out.base()[id] *= rho_m_c_m / ( rho_w * c_w );
    
          }
        }else{
          porosity  = inputdata.yfield_properties.zones[0].porosity;
          rho_s     = inputdata.yfield_properties.zones[0].rho;
          c_s       = inputdata.yfield_properties.zones[0].c_s;
          rho_m_c_m = porosity * rho_w * c_w + ( 1.0 - porosity ) * rho_s * c_s;
    
          out *= rho_m_c_m / ( rho_w * c_w );
        }

      }





      // for logging purposes, returns a string with current date and time
      inline static std::string getDateAndTime(){

        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );

        std::stringstream current_time;
        current_time << (now->tm_year + 1900) << '-' 
                     << std::setw(2) << std::setfill('0') << (now->tm_mon + 1) << '-'
                     << std::setw(2) << std::setfill('0') <<  now->tm_mday
                     << " "
                     << std::setw(2) << std::setfill('0') << now->tm_hour << ":"
                     << std::setw(2) << std::setfill('0') << now->tm_min << ":"
                     << std::setw(2) << std::setfill('0') << now->tm_sec
                     << " ";

        return current_time.str();
      }




      // This function is taken from "convectiondiffusiondg.hh"
      template<class GEO>
      static void element_size( const GEO& geo, 
                                typename GEO::ctype& hmin, 
                                typename GEO::ctype& hmax )
      {
        typedef typename GEO::ctype DF;
        hmin = 1.0E100;
        hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        if (dim==1)
          {
            Dune::FieldVector<DF,dim> x = geo.corner(0);
            x -= geo.corner(1);
            hmin = hmax = x.two_norm();
            return;
          }
        else
          {
            Dune::GeometryType gt = geo.type();
            for (int i=0; i<Dune::ReferenceElements<DF,dim>::general(gt).size(dim-1); i++)
              {
                Dune::FieldVector<DF,dim> x = geo.corner(Dune::ReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,0,dim));
                x -= geo.corner(Dune::ReferenceElements<DF,dim>::general(gt).subEntity(i,dim-1,1,dim));
                hmin = std::min(hmin,x.two_norm());
                hmax = std::max(hmax,x.two_norm());
              }
            return;
          }
      }



      template<typename IndexType>
      static IndexType getIndex( REAL current, REAL delta, bool bCloseUp=true ){

        REAL fraction = current / delta;
        IndexType index = static_cast<IndexType>( fraction );

        if( std::fabs( (REAL)index - fraction ) < GEO_EPSILON  && index && bCloseUp )
          index -= 1;

        if(index<0){
          std::cout << "Error: getIndex: index = " << index << " < 0" << std::endl;
          std::cout << "Error: getIndex: Set index = 0" << std::endl;
          index = 0;
        }
        return index;
        
      }

      template<typename IDT>
      static void add_wells_to_field( const IDT& inputdata,
                                      Vector<REAL> & YField ){

        const int dim = inputdata.domain_data.dim;

        for(UINT iSetup=0; iSetup<inputdata.setups.size(); iSetup++){
          for(UINT iWell=0; iWell<inputdata.setups[iSetup].wdlist.pointdata_vector.size(); iWell++){

            std::size_t myIndex = 0;

            /* 
               Translate the global coordinate '(x,y,z)' into its corresponding Yfield-tensor-index 'global_index'. 
               This requires only the virtual gridsize of the virtual Yfield grid. 
               So 'global_index' is just a virtual index.
            */
            Vector<UINT> global_index;
            global_index.resize(dim);

            global_index[0] = getIndex<UINT>( inputdata.setups[iSetup].wdlist.pointdata_vector[iWell].x,
                                              inputdata.domain_data.virtual_gridsizes[0] );

#ifdef DIMENSION3
            global_index[1] = getIndex<UINT>( inputdata.setups[iSetup].wdlist.pointdata_vector[iWell].y,
                                              inputdata.domain_data.virtual_gridsizes[1] );
#endif

            UINT l_int = getIndex<UINT>( inputdata.setups[iSetup].wdlist.pointdata_vector[iWell].well_bottom,
                                         inputdata.domain_data.virtual_gridsizes[inputdata.dim-1],
                                         false );

            UINT u_int = getIndex<UINT>( inputdata.setups[iSetup].wdlist.pointdata_vector[iWell].well_top,
                                         inputdata.domain_data.virtual_gridsizes[inputdata.dim-1] );

            double well_conductivity = inputdata.setups[iSetup].wdlist.pointdata_vector[iWell].well_conductivity;

            for(UINT ii=l_int; ii<=u_int; ii++){

              global_index[dim-1] = ii;
              int myIndex = General::indexconversion_to_1d( global_index,
                                                            inputdata.domain_data.nCells );

              if( myIndex >= YField.size() )
                std::cout << "ERROR: add_wells_to_field: myIndex = " << myIndex 
                          << " > YField.size() = " << YField.size()
                          << std::endl;
              
              YField[ myIndex ] = well_conductivity;
              
              std::cout << " myIndex = " << myIndex 
                        << " <=> globalindex = (" << global_index[0]
#ifdef DIMENSION3
                        << "," << global_index[1]
#endif
                        << "," << global_index[dim-1]
                        << ") to " << well_conductivity 
                        << " for Setup " << iSetup
                        << " for iWell " << iWell 
                        << std::endl;
            }
            
            /*            
            if(l_int>=0){
#ifdef DIMENSION3
              int myindex = indexconversion_3d_to_1d(l_int-1,y_int,x_int,inputdata.domain_data.nCells[2],inputdata.domain_data.nCells[1],inputdata.domain_data.nCells[0]);
#else
              int myindex = indexconversion_2d_to_1d(l_int-1,x_int,inputdata.domain_data.nCells[1],inputdata.domain_data.nCells[0]);
#endif    
              YField[ myindex ]=-100.0;
              std::cout << " l_int = " << l_int
                        << " myindex = " << myindex << std::endl;
            }
            
            if(u_int<inputdata.domain_data.nCells[inputdata.dim-1]-1){
#ifdef DIMENSION3
              int myindex = indexconversion_3d_to_1d(u_int+1,y_int,x_int,inputdata.domain_data.nCells[2],inputdata.domain_data.nCells[1],inputdata.domain_data.nCells[0]);
#else
              int myindex = indexconversion_2d_to_1d(u_int+1,x_int,inputdata.domain_data.nCells[1],inputdata.domain_data.nCells[0]);
#endif    
              YField[ myindex ]=-100.0;
              std::cout << " u_int = " << u_int
                        << " myindex = " << myindex << std::endl;
            }
            */
           
          }
        
        }

      }


      
      template<typename Coord,typename GV>
      static bool doesPointBelongToGridPartition( const Coord& xglobal, const GV& gv ){

        enum{ dim = GV::dimension };
        typedef typename GV::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementLeafIterator;
        for (ElementLeafIterator it = gv.template begin<0,Dune::All_Partition> ()
               ; it != gv.template end<0,Dune::All_Partition> (); ++it) {
          
          Dune::FieldVector<REAL,dim> xlocal = it->geometry().local(xglobal);
          Dune::GeometryType gt = it->geometry().type();
          
          if( Dune::ReferenceElements<CTYPE,dim>::general(gt).checkInside( xlocal ) ){
            return true;
          }
          
        }
        return false;
        
      }




      template<typename GV,typename REAL>
      static void vector2histogram( const GV& gv,
                                    const Vector<REAL>& v, 
                                    const int nBins,
                                    const std::string filename
                                    )
      {
        REAL maximum = -1e12;
        REAL minimum = 1e12;
        for( int i=0; i<v.size(); ++i ){
          maximum = std::max( maximum, v[i] );
          minimum = std::min( minimum, v[i] );
        }

        REAL interval = ( maximum - minimum ) / REAL( nBins ) + 1E-6;

        Vector<REAL> bins(nBins,0.0);
        int N = v.size();
        for( int i=0; i<v.size(); ++i ){
          for( int j=0; j<nBins; ++j ){
            if( v[i] >= minimum + j*interval 
                && 
                v[i] < minimum + (j+1)*interval ){
              bins[j] += 1.0;
              //break;
            }
          }
        }

        Vector<REAL> sbins(nBins,0.0);
        if( gv.comm().size() > 1 ) {
          // parallel case: sum over all processes!
          MPI_Reduce( &bins[0],
                      &sbins[0],
                      nBins,
                      MPI_DOUBLE,
                      MPI_SUM,
                      0,
                      gv.comm() );
          N = gv.comm().sum( N );
        }
        else {
          sbins = bins;
        }


        if( gv.comm().rank() == 0 ) {

          std::cout << "Histogram for " 
                    << N << " values on " 
                    << nBins << " bins written to data file " 
                    << filename.c_str()
                    << " on process " << gv.comm().rank()
                    << std::endl;

          std::ofstream outfile( filename.c_str(), std::ios::out );
          if( outfile.is_open() ) {

            REAL total = 0.0;
            for( int j=0; j<nBins; ++j ){
              outfile << std::fixed << std::setw(10) << std::setprecision(2) << minimum + (j+0.5) * interval
                      << "  "
                      << std::scientific << std::setw(10) << std::setprecision(2)
                      << sbins[j] / REAL(N) / interval
                      << std::scientific << std::setw(10) << std::setprecision(2)
                      << sbins[j]
                      << std::scientific << std::setw(10) << std::setprecision(2)
                      << sbins[j] / REAL(N)
                      << std::endl;
              total += sbins[j];
            }
            std::cout << "Total number of values = " << total << std::endl;
            std::cout << "Normalized integral = " << total / REAL(N) / interval << std::endl;
          
            outfile.close();

          }


          std::ofstream outfile2( filename + ".vec.dat", std::ios::out );
          if( outfile2.is_open() ) {
            for( int i=0; i<v.size(); ++i ){
              outfile2
                << std::scientific << std::setw(10) << std::setprecision(2)
                << v[i]
                << std::endl;
            }
            outfile2.close();
          }

        }


      }


      static void systemCall( const std::string cmdline ){
        Dune::Timer watch;
        int tmp = system( cmdline.c_str() );
        if(tmp!=0)
          std::cout << "Warning: systemCall to " << cmdline 
                    << " returns " << tmp << std::endl;
        assert(tmp==0);
        std::stringstream jobtitle;
        jobtitle << "systemCall: " << cmdline;
        General::log_elapsed_time( watch.elapsed(),
                                   General::verbosity,
                                   "IO",
                                   jobtitle.str() );
      }



    }; // class General

    int General::verbosity=0;         // must be initialized here

    double General::elapsed_MARK=0;   // must be initialized here
    double General::elapsed_CHECK=0;  // must be initialized here
    
    double General::elapsed_IO=0;     // must be initialized here
    double General::elapsed_FFTW=0;   // must be initialized here
    double General::elapsed_PDE=0;    // must be initialized here
    double General::elapsed_EVAL=0;   // must be initialized here
    double General::elapsed_REDIST=0;   // must be initialized here



  } // Gesis
} // Dune


#endif	/* DUNE_GESIS_GENERAL_HH */

