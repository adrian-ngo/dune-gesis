// Note: "config.h" should be included first!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// Set the world and grid dimension here!
#define DIMENSION3
#include "my_macros.hh"

#ifndef GEO_EPSILON
#define GEO_EPSILON 1e-6
#endif

#ifdef FLOATING_POINT_EXEPTION_CHECK
#include "fenv.h"
#endif

#include<iostream>
#include<fstream>
#include<vector>
#include<random>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
//#include<dune/common/static_assert.hh>

#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/rt0cube3dfem.hh>
#include<dune/pdelab/constraints/raviartthomas0.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/finiteelementmap/opbfem.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/p0.hh>

#ifdef USE_ALUGRID
#include<dune/pdelab/finiteelementmap/p0ghostconstraints.hh>
#endif

#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>

#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>


#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>

#include <fftw3.h>
#include <fftw3-mpi.h>

#include <unistd.h>

#include <dune/gesis/common/datatypes.hh>
#include <dune/gesis/common/io/logfile.hh>
#include <dune/gesis/common/io/inputfile.hh>
#include <dune/gesis/common/general.hh>
#include <dune/gesis/common/io/IO_locations.hh>
#include <dune/gesis/common/io/IO_routines.hh>
#include <dune/gesis/common/io/VTKPlot.hh>
#include <dune/gesis/yfield/FFTFieldGenerator.hh>
#include <dune/gesis/common/YaspPartition.hh>
#include <dune/gesis/common/initialize_DuneGrid.hh>


CLogfile logger;


int main(int argc, char** argv) {

  Dune::Timer main_timer;
  main_timer.reset();

#ifdef DIMENSION3
  const UINT dim = 3;
#else
  const UINT dim = 2;
#endif
  
  //define the outputpath 
  std::string outputpath="./";
  std::string inputfile="./InputFile3d.xml";

#ifdef FLOATING_POINT_EXEPTION_CHECK
    // activate floating point exception handler
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    try {
      //Maybe initialize Mpi
      Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

      fftw_mpi_init();

      if(helper.rank()==0){
        std::cout << Dune::GeoInversion::General::getDateAndTime() 
                  << "Program starts." 
                  << std::endl
                  << "DUNE Geostatistical Inversion - stationary (dune-gsis)"
                  << std::endl
                  << "3D-Version (Build number: " 
                  << MAJOR_VERSION << "."
                  << MINOR_VERSION << "."
                  << PDELAB_VERSION << "."
                  << BUILD_VERSION << ")"
                  << std::endl
                  << std::endl
          ;
      }

      bool bContinue = false;
      
      // work on options!
        for(INT ii=0; ii<argc; ii++){
            
          // outputpath  
          if( strcmp(argv[ii],"-o")==0 && argc>ii+1){
               if(helper.rank()==0)
                std::cout<<"Output path specified!!"<<std::endl;
               outputpath=argv[ii+1];
          }
          // inputfile 
          if( strcmp(argv[ii],"-i")==0 && argc>ii+1){
               if(helper.rank()==0)
                std::cout<<"Input file specified!!"<<std::endl;
               inputfile=argv[ii+1];
          }

          // --continue
          if( strcmp(argv[ii],"-c")==0 ){
            if(helper.rank()==0)
              std::cout << "Continuing computations with existing results ..." << std::endl;
            bContinue = true;
          }

          if( strcmp(argv[ii],"--continue")==0 ){
            if(helper.rank()==0)
              std::cout << "Continuing computations with existing results ..." << std::endl;
            bContinue = true;
          }

        }
      
        if(helper.rank()==0)
          std::cout<<"The inputfile is: "<<inputfile<<std::endl; 
        
        //check if writing is guranteed in the outputpath 
        if( access(outputpath.c_str(), W_OK)!=0){
           if(helper.rank()==0){
            std::cout<<" The output directory >>"<<outputpath<<"<< is NOT writeable!!! --> QUIT!"<<std::endl;
            std::cerr<<" The output directory >>"<<outputpath<<"<< is NOT writeable!!! --> QUIT!"<<std::endl;
           }

           return -1;
        }else{
           if(outputpath[outputpath.size()-1]!='/')
               outputpath+="/";
           if(helper.rank()==0)
             std::cout<<"The output directory is: "<<outputpath<<std::endl; 
        }
        
        if(access(outputpath.c_str(), W_OK)!=0){
            if(helper.rank()==0){
                std::cout<<" The output directory >>"<<outputpath<<"<< does NOT exist!!! --> QUIT! ( type: mkdir "<<outputpath<<" )"<<std::endl;
                std::cerr<<" The output directory >>"<<outputpath<<"<< does NOT exist!!! --> QUIT! ( type: mkdir "<<outputpath<<" )"<<std::endl;
            }
            return -1;
        }

        typedef Dune::GeoInversion::IO_Locations DIR;
        DIR dir(dim,outputpath,inputfile);
      
      if (CLEAN_OUTPUT && helper.rank()==0){

        std::cout << "Cleaning old LOG and VTU files ..." << std::endl;
        Dune::GeoInversion::General::systemCall( "rm -rvf " + dir.logdir + "/*" );
        Dune::GeoInversion::General::systemCall( "rm -rvf " + dir.vtudir + "/*" );
        Dune::GeoInversion::General::systemCall( "rm -rvf " + dir.hdf5dir + "/*" );

      }

      // All processes please wait until the cleaning on root is done!
      if(helper.size()>1)
        MPI_Barrier(helper.getCommunicator());

      if( !bContinue && helper.rank()==0){

        int tmp;
        std::string cmdline;

        /*
        //delete buffer BUT keep the needed stuff if using_existing_Yold
        cmdline="test -f " + dir.Y_estimated_h5file;
        if(! system(cmdline.c_str())){
          cmdline="mv "+dir.Y_estimated_h5file+" "+outputpath+".3DY_estimated.h5";
          tmp=system(cmdline.c_str());
        }

        cmdline="test -f " + dir.ksi_old_h5file;
        if(! system(cmdline.c_str())){
          cmdline="mv "+dir.ksi_old_h5file+" "+outputpath+".3Dksi_old.h5";
          tmp=system(cmdline.c_str());
        }
    
        cmdline="test -f " + dir.beta_old_h5file;
        if(! system(cmdline.c_str())){
          cmdline="mv "+dir.beta_old_h5file+" "+outputpath+".3Dbeta_old.h5";
          tmp=system(cmdline.c_str());
        }
      
        cmdline="test -f " + dir.L_prior_h5file;
        if(! system(cmdline.c_str() )){
          cmdline="mv "+dir.L_prior_h5file+" "+outputpath+".3DL_prior.h5";
          tmp=system(cmdline.c_str());
        }   
      
        cmdline="test -f " + dir.estimation_variance_h5file;
        if(! system(cmdline.c_str())){
          cmdline="mv "+ dir.estimation_variance_h5file +" "+outputpath+".3DEstimation_Variance.h5";
          tmp=system(cmdline.c_str());
        }
        */

        std::cout << "Cleaning old BUFFER files ..." << std::endl;
        Dune::GeoInversion::General::systemCall( "rm -rvf " + dir.bufferdimdir + "/*" );

        /*
        cmdline="test -f "+outputpath+".3DL_prior.h5";
        if(! system(cmdline.c_str())){
          cmdline="mv "+outputpath+".3DL_prior.h5 "+dir.L_prior_h5file;
          tmp=system(cmdline.c_str());
        }
      
        cmdline="test -f "+outputpath+".3Dbeta_old.h5";
        if(! system(cmdline.c_str())){
          cmdline="mv "+outputpath+".3Dbeta_old.h5 "+dir.beta_old_h5file; 
          tmp=system(cmdline.c_str());
        } 
      
        cmdline="test -f "+outputpath+".3Dksi_old.h5";
        if(! system(cmdline.c_str())){
          cmdline="mv "+outputpath+".3Dksi_old.h5 "+dir.ksi_old_h5file;
          tmp=system(cmdline.c_str());
        } 
      
        cmdline="test -f "+outputpath+".3DY_estimated.h5";
        if(! system(cmdline.c_str())){
          cmdline="mv "+outputpath+".3DY_estimated.h5 "+dir.Y_estimated_h5file;
          tmp=system(cmdline.c_str());
        } 

        cmdline="test -f "+outputpath+".3DEstimation_Variance.h5";
        if(! system(cmdline.c_str())){
          cmdline="mv "+outputpath+".3DEstimation_Variance.h5 " + dir.estimation_variance_h5file;
          tmp=system(cmdline.c_str());
        }

        */

      }

      // All processes please wait until the cleaning on root is done!
      if(helper.size()>1)
        MPI_Barrier(helper.getCommunicator());

        std::stringstream slogfilename;
        slogfilename << dir.logdir << "/Logfile_P" << helper.rank() << ".log";
        
        if(helper.rank()==0){

          Dune::GeoInversion::General::createDirectory( dir.datadir );
          Dune::GeoInversion::General::createDirectory( dir.bufferdir );
          Dune::GeoInversion::General::createDirectory( dir.outputdir );
          
          Dune::GeoInversion::General::createDirectory( dir.datadimdir );
          Dune::GeoInversion::General::createDirectory( dir.bufferdimdir );
          Dune::GeoInversion::General::createDirectory( dir.outputdimdir );
          
          Dune::GeoInversion::General::createDirectory( dir.yfielddir );
          Dune::GeoInversion::General::createDirectory( dir.CRdir );
          Dune::GeoInversion::General::createDirectory( dir.logdir );
          Dune::GeoInversion::General::createDirectory( dir.vtudir );
          Dune::GeoInversion::General::createDirectory( dir.hdf5dir );
          
        }

        if(helper.size()>1)
          MPI_Barrier(helper.getCommunicator());

	logger.init( slogfilename.str() );
	logger << "Logfile " << slogfilename.str() << " for process " << helper.rank() << std::endl;

        // All processes please wait until all logfiles are opened!
        if(helper.size()>1)
          MPI_Barrier(helper.getCommunicator());
		
        if (Dune::MPIHelper::isFake)
		  logger << "This is a sequential program." << std::endl;
        else {
		  if (helper.rank() == 0)
			{
			  std::cout <<  "--------------------------------------------------------------" << std::endl;
			  std::cout << "Parallel version run with " << helper.size() << " process(es)" << std::endl;
			  std::cout << "dim = " << dim << std::endl;
			  std::cout << "--------------------------------------------------------------" << std::endl;
			}
			logger << "--------------------------------------------------------------" << std::endl;
			logger << "Parallel version run with " << helper.size() << " process(es)" << std::endl;
			logger << "dim = " << dim << std::endl;
			logger << "--------------------------------------------------------------" << std::endl;
        }

    CInputData inputdata( helper );
    if( !inputdata.readInputFileXml(dir.inputfile) )
      return 2; // missing input-file

    inputdata.writeToScreen();

    Dune::GeoInversion::General::verbosity = inputdata.verbosity;
    inputdata.problem_types.using_existing_Yold = bContinue;


    // set the directories for zonation!
    dir.set_zonation(inputdata.yfield_properties.nz);
    
    // set the directories for inversion!
    dir.set_setups(inputdata.setups.size());
       
    //dir.set_heat_tomography(inputdata.heat_tomography_inversion_data.size());
    for(UINT ii=0; ii<inputdata.setups.size(); ii++ )
        dir.set_geoelectrical_potential_tomography(ii , inputdata.setups[ii].geoelectrical_potential_inversion_data.nconfig);
/*	
	// if a command line argument is given, then it is the number of MPI Pools which should be used for the sensitivity calculations
	if(argc>1){
	  //set the maxNComm value in the inputdata structure
	  // which defines the maximum number of comunicators in each MPI pool.
	  inputdata.maxNComm=atoi(argv[1]);
	  //check if as many pools as requested can be created
	  if(inputdata.maxNComm>(UINT)helper.size()){
	    // more pools requested than processor available -> Quitting
	    if(helper.rank()==0)
	      std::cout<<"can't use "<<inputdata.maxNComm<<" when only "<<helper.size()<<" processors available! -> QUITTING!"<<std::endl;
	    logger <<"can't use "<<inputdata.maxNComm<<" when only "<<helper.size()<<" processors available! -> QUITTING!"<<std::endl;
	    return -1;
	  }else{
        if(helper.rank()==0)
          std::cout<<"using max. "<<inputdata.maxNComm<<" comunicators for parallel in parallel calculations."<<std::endl;
        logger<<"using max. "<<inputdata.maxNComm<<" comunicators for parallel in parallel calculations."<<std::endl;
      }
	}
	*/

	if(helper.size()>1)
	  MPI_Barrier(helper.getCommunicator());
	   
        if( helper.rank() == 0 ) {
          if( inputdata.problem_types.generate_measurement_data ) {
            std::cout << "Cleaning *meas* files inside datadimdir ..." << std::endl;
            Dune::GeoInversion::General::systemCall( "rm -rvf " + dir.datadimdir + "/*meas*" );
          }
        }

        //definitions for the Y-Field
        typedef Dune::GeoInversion::FFTFieldGenerator<CInputData,REAL,dim> YFG;
        YFG Yfieldgenerator( inputdata,dir,helper.getCommunicator() );
        //generate the Y_field
        Yfieldgenerator.init();


        // start the main-work-flow
        double timeCounted = Dune::GeoInversion::initialize_DuneGrid<CInputData,YFG,DIR,dim>( inputdata, Yfieldgenerator, dir, helper );

        if( helper.rank()==0 && inputdata.verbosity > 0 ){
          double elapsedTime = main_timer.elapsed();

          std::cout << "TIMER: The program took " 
                    << elapsedTime
                    << " sec. in total. Time measured in detail covers "
                    << 100.0 * timeCounted / elapsedTime
                    << " percent."
                    << std::endl;
        }
        //everything is done!

        if(helper.rank()==0){
          std::cout << Dune::GeoInversion::General::getDateAndTime() 
                    << "3D-Program ends."
                    << std::endl;
        }
	
        return 0;

    } catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}
