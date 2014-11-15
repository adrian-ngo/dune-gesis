/* 
 * File:   VTKPlot.hh
 * Author: A. Ngo
 *
 * 2010-2014
 */

#ifndef DUNE_GESIS_VTU_ROUTINES_HH
#define	DUNE_GESIS_VTU_ROUTINES_HH

#include <dune/pdelab/common/geometrywrapper.hh>
#include <assert.h>
#include <sstream>

#include "IO_routines.hh"

extern CLogfile logger;


namespace Dune {
  namespace Gesis {


    class VTKPlot{

    private:
      VTKPlot(){};

    public:

      // There are three functions with different signatures serving 
      // different input types 
      // 
      // 1.) output2vtu:                    VCType    ---> vertex data
      // 2.) output_vector_to_vtu:          VECTOR    ---> cell-data / vertex data
      // 3.) output_dgf_to_vtu:             DGF       ---> vertex-data
      // 4.) output_hdf5data_to_gridview    HDF5 file ---> cell-data / vertex data
      // 5.) plotDataToRefinedGrid          HDF5 file ---> HDF5 file
      // 


      // *******************************************************************
      // 1.) output2vtu:                    VCType    ---> vertex data
      // *******************************************************************
      template<typename GFS,typename VCType>
      static void output2vtu( const GFS& gfs
                              , const VCType& xSolution
                              , const std::string filename
                              , const std::string title 
                              , const int verbosity=0
                              , bool bSubSampling = false
                              , int nSubSampling = 0
                              , bool bCellData = false
                              )
      {
        Dune::Timer watch;
        logger << "output2vtu: " << filename << std::endl;

        //Dune::VTK::OutputType vtkOutputType = Dune::VTK::ascii;
        Dune::VTK::OutputType vtkOutputType = Dune::VTK::appendedraw;


        typedef typename GFS::Traits::GridViewType GV;
        GV gv = gfs.gridView();
        typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> DGF;
        DGF dgf(gfs,xSolution);

        if( gv.comm().rank()==0 
            &&
            verbosity >= VERBOSITY_DEBUG_LEVEL )
          std::cout << "VTK output of vertex-data " << title.c_str() << " to file '" << filename.c_str() << ".vtu'" <<  std::endl;

        if(bSubSampling){
    
#ifdef PLOT_GHOST
          Dune::SubsamplingVTKWriter<GV,Dune::All_Partition> vtkwriter(gv,nSubSampling);
#else
          Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,nSubSampling);
#endif
          if( bCellData )
            vtkwriter.addCellData( new Dune::PDELab::VTKGridFunctionAdapter<DGF > ( dgf, "celldata" ) );
          else
            vtkwriter.addVertexData( new Dune::PDELab::VTKGridFunctionAdapter<DGF> ( dgf, title.c_str() ) );
          
          //find '/' or '\' in the filename
          size_t found=filename.find_last_of("/\\");
          
          //if parallel and complete path is given, use the pwrite() function
          if(found<filename.size()&& gv.comm().size() > 1) {
            vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", vtkOutputType );
          }
          else {
            vtkwriter.write( filename.c_str(), vtkOutputType );
          }
        }
        else{

          Dune::VTKWriter<GV> vtkwriter(gv);
          if( bCellData )
            vtkwriter.addCellData( new Dune::PDELab::VTKGridFunctionAdapter<DGF > ( dgf, "celldata" ) );
          else
            vtkwriter.addVertexData( new Dune::PDELab::VTKGridFunctionAdapter<DGF> ( dgf, title.c_str() ) );

          //find '/' or '\' in the filename
          size_t found=filename.find_last_of("/\\");
  
          //if parallel and complete path is given, use the pwrite() function
          if(found<filename.size()&& gv.comm().size() > 1) {
            vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", vtkOutputType );
          }
          else {
            vtkwriter.write( filename.c_str(), vtkOutputType );
          }

        }

        std::stringstream jobtitle;
        jobtitle << "VTK output to "
                 << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   verbosity,
                                   "IO",
                                   jobtitle.str() );


        return;

      }




      // *************************************************************************
      // 2.) output_vector_to_vtu:          VECTOR    ---> cell-data / vertex data
      // *************************************************************************
      template<typename GV, typename VECTOR >
      static void output_vector_to_vtu( 
                                       const GV& gv
                                       , const VECTOR& vCellData
                                       , const std::string filename
                                       , const std::string title 
                                       , const int verbosity=0
                                       , bool bSubSampling = false
                                       , int nSubSampling = 0
                                       , bool bCellData = false // special plot
                                        )
      {
        Dune::Timer watch;
        logger << "output_vector_to_vtu: " << filename << std::endl;

        //Dune::VTK::OutputType vtkOutputType = Dune::VTK::ascii;
        Dune::VTK::OutputType vtkOutputType = Dune::VTK::appendedraw;

        if( gv.comm().rank()==0 
            &&
            verbosity >= VERBOSITY_DEBUG_LEVEL )
          std::cout << "VTK output of cell-data " << title.c_str() << " to file '" << filename.c_str() << ".vtu'" <<  std::endl;
  

        if(bSubSampling){

          // This will plot the cell data as they are, 
          // namely in a non-conforming way.

          const int dim = GV::dimension;
          typedef Dune::PDELab::P0LocalFiniteElementMap<CTYPE,REAL,dim> P0FEM;
          typedef Dune::PDELab::NoConstraints NOCONS;
          typedef Dune::PDELab::ISTLVectorBackend<> VBE1;
    
#ifdef USE_SIMPLEX
          P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));
#else
          P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
#endif

          typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,NOCONS,VBE1> P0GFS;
          P0GFS p0gfs(gv,p0fem);

          typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,REAL>::Type P0VCType;
          P0VCType p0Celldata(p0gfs,0.0);
          //p0Celldata.std_copy_from( vCellData );
          General::copy_from_std( vCellData, p0Celldata );


          std::stringstream jobtitle;
          jobtitle << "output_vector_to_vtu: copying data";
          General::log_elapsed_time( watch.elapsed(),
                                     gv.comm(),
                                     verbosity,
                                     "IO",
                                     jobtitle.str() );

          output2vtu( p0gfs,
                      p0Celldata,
                      filename.c_str(),
                      title.c_str(),
                      verbosity,
                      bSubSampling,
                      nSubSampling,
                      bCellData
                      );
        }
        else{

#ifdef PLOT_GHOST
          Dune::VTKWriter<GV,Dune::All_Partition> vtkwriter( gv );
#else
          Dune::VTKWriter<GV> vtkwriter(gv);
#endif
          vtkwriter.addCellData( vCellData, title.c_str() );
  
          //find '/' or '\' in the filename
          size_t found=filename.find_last_of("/\\");
  
          //if parallel and complete path is given, use the pwrite() function
          if(found<filename.size()&& gv.comm().size() > 1) {
            vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", vtkOutputType );
          }
          else {
            vtkwriter.write( filename.c_str(), vtkOutputType );
          }

          std::stringstream jobtitle;
          jobtitle << "output_vector_to_vtu: writing " << filename;
          General::log_elapsed_time( watch.elapsed(),
                                     gv.comm(),
                                     verbosity,
                                     "IO",
                                     jobtitle.str() );

        }

      }



      // *******************************************************************
      // 3.) output_dgf_to_vtu:             DGF       ---> vertex-data
      // *******************************************************************
      template<typename GV, typename GFS, typename DGF >
      static void output_dgf_to_vtu( 
                                    const GV& gv
                                    , const GFS& gfs
                                    , const DGF& dgf
                                    , const std::string filename
                                    , const std::string title 
                                    , const int verbosity=0
                                    , bool bSubSampling = false
                                    , int nSubSampling = 0
                                     )
      {

        Dune::Timer watch;
        logger << "output_dgf_to_vtu: " << filename << std::endl;

        //Dune::VTK::OutputType vtkOutputType = Dune::VTK::ascii;
        Dune::VTK::OutputType vtkOutputType = Dune::VTK::appendedraw;

        if( gv.comm().rank()==0 &&
            verbosity >= VERBOSITY_DEBUG_LEVEL )
          std::cout << "VTK output of vertex-data '" << title.c_str() << "' to file '" << filename.c_str() << ".vtu'" <<  std::endl;
  
        if( bSubSampling )
          {
#ifdef PLOT_GHOST
            Dune::SubsamplingVTKWriter<GV,Dune::All_Partition> ss_vtkwriter( gv, nSubSampling );
#else
            Dune::SubsamplingVTKWriter<GV> ss_vtkwriter( gv, nSubSampling );
#endif
            ss_vtkwriter.addVertexData( new Dune::PDELab::VTKGridFunctionAdapter<DGF > ( dgf, title.c_str() ) );
            ss_vtkwriter.addCellData( new Dune::PDELab::VTKGridFunctionAdapter<DGF > ( dgf, "celldata" ) );
          
            //find '/' or '\' in the filename
            size_t found=filename.find_last_of("/\\");
	  
            //if parallel and complete path is given, use the pwrite() function
            if(found<filename.size()&& gv.comm().size() > 1) {
              ss_vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", vtkOutputType );
            }
            else {
              ss_vtkwriter.write( filename.c_str(), vtkOutputType );
            }
          }
        else
          {

            Dune::VTKWriter<GV> vtkwriter( gv ); // Default is conforming!
            //Dune::VTKWriter<GV> vtkwriter( gv, Dune::VTK::nonconforming );

            vtkwriter.addVertexData( new Dune::PDELab::VTKGridFunctionAdapter<DGF> ( dgf, title.c_str() ) );
	  
            //find '/' or '\' in the filename
            size_t found=filename.find_last_of("/\\");
	  
            //if parallel and complete path is given, use the pwrite() function
            if(found<filename.size()&& gv.comm().size() > 1) {
              vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", vtkOutputType );
            }
            else {
              vtkwriter.write( filename.c_str(), vtkOutputType );
            }
          }

        std::stringstream jobtitle;
        jobtitle << "DGF output to "
                 << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   verbosity,
                                   "IO",
                                   jobtitle.str() );



      }



      // *************************************************************************
      // 4.) output_hdf5data_to_gridview    HDF5 file ---> cell-data / vertex data
      // *************************************************************************
      template<typename GV,typename YFG>
      inline static void output_hdf5data_to_gridview( const GV &gv_0,
                                                      const CInputData& inputdata,
                                                      const std::string & vtu_filename,
                                                      const std::string & Y_old_filename,
                                                      const std::string & hd5_data_name,
                                                      YFG& yfg_Y_old,
                                                      const int verbosity=0
                                                      ){

        Dune::Timer watch;
        logger << "output_hdf5data_to_gridview: " << vtu_filename << std::endl;


        Vector<REAL> local_Y_old;
        Vector<UINT> local_count,local_offset;

        HDF5Tools::blob();

        HDF5Tools::read_parallel_from_HDF5( gv_0
                                            , inputdata
                                            , local_Y_old
                                            , hd5_data_name
                                            , local_count
                                            , local_offset
                                            , Y_old_filename
                                            );

        if( gv_0.comm().size() > 1 ){
          yfg_Y_old.parallel_import_from_local_vector( local_Y_old,
                                                       local_count,
                                                       local_offset );
        }
        else {
          yfg_Y_old.import_from_vector( local_Y_old );
        }

        yfg_Y_old.plot2vtu( gv_0, vtu_filename, hd5_data_name );

        std::stringstream jobtitle;
        jobtitle << "VTK output of HDF5 data to "
                 << vtu_filename;
        General::log_elapsed_time( watch.elapsed(),
                                   gv_0.comm(),
                                   verbosity,
                                   "IO",
                                   jobtitle.str() );

      }






      // *******************************************************************
      // 5.) plotDataToRefinedGrid          HDF5 file ---> HDF5 file
      // *******************************************************************
      template< typename GRID,
                typename GV_GW,
                typename IDT,
                typename DIR,
                typename YFG
                >
      static void plotDataToRefinedGrid( GRID& theGrid,
                                         const GV_GW& gv_0,
                                         const std::string filename1,
                                         const std::string filename2,
                                         const std::string groupname,
                                         const IDT& inputdata,
                                         const DIR& dir
                                         ) {
      
        //const int dim = GV_GW::dimension;
        // What is this good for? This is for step 3 of ...
        // Idea:
        // 1.) Generate measurements on grid level L1.
        // 2.) Copy the measurement data from BUFFER into DATA sub-dir.
        // 3.) Run a first inversion on a coarse grid level L0 with refine_estimate = "1". This will produce "Y_estimated_2.h5" with "/Y_est".
        // 4a) mv Y_estimated_2.h5 Y_estimated.h5
        // 4b) rm L_prior.h5 (very important: Its value from L0 will be too small for L1.)
        // 5.) Run a second inversion on level L1 with using_existing_Yold="yes" and start geo_inversion with -c.
        // 
        // Background: Step 3) works as a speeded up the generation of an initial estimation for step 5)
        // For example, one can start with head inversion in step 3) and carry on with m0m1 inversion in step 5)
        //
          
        Vector<UINT> local_count,local_offset;
        Vector<REAL> Y_est_parallel;

        HDF5Tools::
          read_parallel_from_HDF5( gv_0
                                 , inputdata
                                 , Y_est_parallel
                                 , groupname
                                 , local_count
                                 , local_offset
                                 , filename1
                                 , 1 // P0 blocksize
                                 , FEMType::DG // P0
                                 , 0 // structure is on grid level 0
                                 );
        YFG yfg_Y_est( inputdata, dir, gv_0.comm() );

        if( gv_0.comm().size() > 1 )
          yfg_Y_est.parallel_import_from_local_vector( Y_est_parallel,
                                                       local_count,
                                                       local_offset );
        else
          yfg_Y_est.import_from_vector( Y_est_parallel );

        theGrid.globalRefine( 1 );
        const GV_GW& gv_1 = theGrid.levelGridView(1);
        Vector<REAL> Y_est2;
        Vector<REAL> Y_est2_well; // dummy, not used

        yfg_Y_est.export_field_to_vector_on_grid( gv_1,
                                                  Y_est2,
                                                  Y_est2_well, // dummy
                                                  1 // grid level for mapper
                                                  );


        VTKPlot::output_vector_to_vtu( gv_1, 
                                       Y_est2,
                                       dir.vtudir + "/Y_est2",
                                       "Y_est2",
                                       inputdata.verbosity,
                                       true, // subsampling
                                       0 // subsampling degree
                                       );


        HDF5Tools::
          write_parallel_to_HDF5(
                                 gv_1
                               , inputdata
                               , Y_est2
                               , groupname
                               , inputdata.domain_data.nCells
                               , filename2
                               , 1
                               , FEMType::DG
                               , 1
                               , true
                               );
        theGrid.globalRefine(-1);

      } // end of void plotDataToRefinedGrid



      /**
         Function to plot the grid ordering to a cell-wise output.
         \tparam GV The gridview type.
         \param gv The gridview.
         \param filename The name of the VTK file.
       */

      template<typename GV>
      static void outputGridviewIndexToDGF( const GV& gv,
                                            const std::string filename,
                                            const std::string meshtype="cube"
                                            ){
        Dune::Timer watch;

        // (1) Construct finite element space
        const int dim = GV::dimension;
        Dune::GeometryType gt;

        if( meshtype == "cube" )
          gt.makeCube(dim);
        else
          gt.makeSimplex(dim);
  
        typedef Dune::PDELab::P0LocalFiniteElementMap<double,double,dim> P0FEM;
        P0FEM p0fem( gt );

        // (2) Set up grid function space
        typedef Dune::PDELab::ISTLVectorBackend<> VBE;

        typedef Dune::PDELab::NoConstraints NOCON;
        NOCON nocon;

        typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,NOCON,VBE> P0GFS;
        P0GFS p0gfs(gv, p0fem, nocon );

        typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,double>::Type UType;    
        UType u0(p0gfs, 0.0);


        for (unsigned int i=0;i<u0.flatsize();i++) {
          u0.base()[i] = REAL(i)/(u0.flatsize()-1);
          // plot the processor rank for debugging purposes:
          // u0[i] = gv.comm().rank();
        }

        /*
        // Alternative loop implementation:
        typedef typename GV::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
        int iElement = 0;
        for (ElementIterator it=gv.template begin<0,Dune::All_Partition>();
        it!=gv.template end<0,Dune::All_Partition>();++it) {

        // plotting element level
        //u0[gv.indexSet().index(*it)] = (*it).level(); 
        // plotting element index
        u0[gv.indexSet().index(*it)] = gv.indexSet().index(*it);
        iElement++;
        }
        */
    
        typedef Dune::PDELab::DiscreteGridFunction<P0GFS,UType> P0DGF;
        P0DGF u0dgf(p0gfs,u0);

        std::stringstream jobtitle;
        jobtitle << "outputGridviewIndexToDGF: copying ";
        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   General::verbosity,
                                   "IO",
                                   jobtitle.str() );

        output2vtu( p0gfs, 
                    u0, 
                    filename, 
                    "elementorder", 
                    General::verbosity, 
                    true, 
                    0
                    );

      }


    }; // class VTKPlot

  } // Gesis

} // Dune

#endif
