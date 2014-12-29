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

#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

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
#include<dune/pdelab/finiteelementmap/rt0cube2dfem.hh>
#include<dune/pdelab/constraints/raviartthomas0.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/finiteelementmap/opbfem.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/p0.hh>


#include<dune/pdelab/constraints/p0ghost.hh>



#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
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
#include <dune/gesis/yfield/FFTFieldGenerator.hh>

#if HAVE_DUNE_ALUGRID || HAVE_ALUGRID
#include <dune/gesis/driver/driverALU.hh>
#endif

CLogfile logger;


int main(int argc, char** argv) {

  Dune::Gesis::General::verbosity = 3;

  Dune::Timer main_timer;

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
      std::cout << Dune::Gesis::General::getDateAndTime() 
                << "Program starts."
                << std::endl;
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


    typedef Dune::Gesis::IO_Locations DIR;
    DIR dir(dim, outputpath,inputfile);

    if (CLEAN_OUTPUT && helper.rank()==0){

      std::cout << "Cleaning old LOG and VTU files ..." << std::endl;
      Dune::Gesis::General::systemCall( "rm -rvf " + dir.logdir + "/*" );
      Dune::Gesis::General::systemCall( "rm -rvf " + dir.vtudir + "/*" );
      Dune::Gesis::General::systemCall( "rm -rvf " + dir.hdf5dir + "/*" );

    }
    // All processes please wait until the cleaning on root is done!
    if(helper.size()>1)
      MPI_Barrier(helper.getCommunicator());

    if( !bContinue && helper.rank()==0){

      //int tmp;
      std::string cmdline;

      /*
      //delete buffer BUT keep the needed stuff if using_existing_Yold
      cmdline="test -f " + dir.Y_estimated_h5file;
      if(! system(cmdline.c_str())){
        cmdline="mv "+dir.Y_estimated_h5file+" "+outputpath+".2DY_estimated.h5";
        tmp=system(cmdline.c_str());
      }

      cmdline="test -f " + dir.ksi_old_h5file;
      if(! system(cmdline.c_str())){
        cmdline="mv "+dir.ksi_old_h5file+" "+outputpath+".2Dksi_old.h5";
        tmp=system(cmdline.c_str());
      }
    
      cmdline="test -f " + dir.beta_old_h5file;
      if(! system(cmdline.c_str())){
        cmdline="mv "+dir.beta_old_h5file+" "+outputpath+".2Dbeta_old.h5";
        tmp=system(cmdline.c_str());
      }
 	  
      cmdline="test -f " + dir.L_prior_h5file;
      if(! system(cmdline.c_str() )){
        cmdline="mv "+dir.L_prior_h5file+" "+outputpath+".2DL_prior.h5";
        tmp=system(cmdline.c_str());
      }	  
	  
      cmdline="test -f " + dir.estimation_variance_h5file;
      if(! system(cmdline.c_str())){
        cmdline="mv "+ dir.estimation_variance_h5file +" "+outputpath+".2DEstimation_Variance.h5";
        tmp=system(cmdline.c_str());
      }
      */

      std::cout << "Cleaning old BUFFER files ..." << std::endl;
      Dune::Gesis::General::systemCall( "rm -rvf " + dir.bufferdimdir + "/*" );

      /*
      cmdline="test -f "+outputpath+".2DL_prior.h5";
      if(! system(cmdline.c_str())){
        cmdline="mv "+outputpath+".2DL_prior.h5 "+dir.L_prior_h5file;
        tmp=system(cmdline.c_str());
      }
 	  
      cmdline="test -f "+outputpath+".2Dbeta_old.h5";
      if(! system(cmdline.c_str())){
        cmdline="mv "+outputpath+".2Dbeta_old.h5 "+dir.beta_old_h5file; 
        tmp=system(cmdline.c_str());
      } 
  	  
      cmdline="test -f "+outputpath+".2Dksi_old.h5";
      if(! system(cmdline.c_str())){
        cmdline="mv "+outputpath+".2Dksi_old.h5 "+dir.ksi_old_h5file;
        tmp=system(cmdline.c_str());
      } 
 	  
      cmdline="test -f "+outputpath+".2DY_estimated.h5";
      if(! system(cmdline.c_str())){
        cmdline="mv "+outputpath+".2DY_estimated.h5 "+dir.Y_estimated_h5file;
        tmp=system(cmdline.c_str());
      } 

      cmdline="test -f "+outputpath+".2DEstimation_Variance.h5";
      if(! system(cmdline.c_str())){
        cmdline="mv "+outputpath+".2DEstimation_Variance.h5 " + dir.estimation_variance_h5file;
        tmp=system(cmdline.c_str());
      }
      */

    }
    // All processes please wait until the cleaning on root is done!
    if(helper.size()>1)
      MPI_Barrier(helper.getCommunicator());

        
    if(helper.rank()==0){

      Dune::Gesis::General::createDirectory( dir.datadir );
      Dune::Gesis::General::createDirectory( dir.bufferdir );
      Dune::Gesis::General::createDirectory( dir.outputdir );
          
      Dune::Gesis::General::createDirectory( dir.datadimdir );
      Dune::Gesis::General::createDirectory( dir.bufferdimdir );
      Dune::Gesis::General::createDirectory( dir.outputdimdir );

      Dune::Gesis::General::createDirectory( dir.yfielddir );
      Dune::Gesis::General::createDirectory( dir.CRdir );
      Dune::Gesis::General::createDirectory( dir.logdir );
      Dune::Gesis::General::createDirectory( dir.vtudir );
      Dune::Gesis::General::createDirectory( dir.hdf5dir );
          
    }

    std::stringstream slogfilename;
    slogfilename << dir.logdir << "/Logfile_P" << helper.rank() << ".log";
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
          std::cout << "--------------------------------------------------------------" << std::endl;
          std::cout << "Parallel version run with " << helper.size() << " process(es)" << std::endl;
          std::cout << "dim = " << dim << std::endl;
          std::cout << "--------------------------------------------------------------" << std::endl;
        }
      logger << "--------------------------------------------------------------" << std::endl;
      logger << "Parallel version run with " << helper.size() << " process(es)" << std::endl;
      logger << "dim = " << dim << std::endl;
      logger << "--------------------------------------------------------------" << std::endl;
			
    }


    typedef Dune::Gesis::CInputData IDT;
    IDT inputdata( helper );
    if( !inputdata.readInputFileXml(dir.inputfile) )
      exit(2); // missing input-file
      
    inputdata.writeToScreen();

    Dune::Gesis::General::verbosity = inputdata.verbosity;
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
        Dune::Gesis::General::systemCall( "rm -rvf " + dir.datadimdir + "/*meas*" );
      }
    }

    // All processes please wait until the cleaning on root is done!
    if(helper.size()>1)
      MPI_Barrier(helper.getCommunicator());

      
    typedef Dune::Gesis::FieldData FD;
    FD fielddata(inputdata);

    //definitions for the Y-Field
    typedef Dune::Gesis::FFTFieldGenerator<FD,DIR,dim> YFG;
    YFG Yfieldgenerator( fielddata, dir, helper.getCommunicator() );
    //generate the Y_field
    Yfieldgenerator.init();

    double timeCounted = 0;
#if HAVE_DUNE_ALUGRID || HAVE_ALUGRID
    // start the main-work-flow
    timeCounted = Dune::Gesis::driverALU<IDT,YFG,DIR,dim>( inputdata, Yfieldgenerator, dir, helper );
#else
    std::cout << "The DUNE module dune-alugrid is required!" << std::endl;
#endif

    if( helper.rank()==0 && inputdata.verbosity > 0 ){
      double elapsedTime = main_timer.elapsed();

      std::cout << "TIMER: The program took " 
                << elapsedTime
                << " sec. in total. Time measured in detail covers "
                << std::fixed
                << 100.0 * timeCounted / elapsedTime
                << " percent."
                << std::endl;
    }
    //everything is done!


    if(helper.rank()==0){
      std::cout << Dune::Gesis::General::getDateAndTime() 
                << "Program ends."
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

