/* 
 * File:   output2vtu.hh
 * Author: ngo
 *
 * Created on July 8, 2010, 1:32 PM
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

      //
      // Take care!
      // There are three functions with different signatures serving 
      // different input types 
      // 
      // 1.) output2vtu:            VCType ---> vertex data
      // 2.) output_vector_to_vtu:  VECTOR ---> cell-data
      // 3.) output_dgf_to_vtu:     DGF ---> vertex-data
      //

      // Version 1:
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

        REAL elapsed_time;

        typedef typename GFS::Traits::GridViewType GV;
        GV gv = gfs.gridView();
        typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> DGF;
        DGF dgf(gfs,xSolution);

        //std::cout << "DEBUG: xSolution.flatsize() = " << xSolution.flatsize() << std::endl;
        //std::cout << "DEBUG: gfs.size() = " << gfs.size() << std::endl;
        //std::cout << "DEBUG: gridview.size(0) = " << gv.size(0) << std::endl;

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
          if(found<filename.size()&& gv.comm().size() > 1)
            {
#ifdef VTK_PLOT_ASCII
              vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::ascii );
#else
              vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::appendedraw );
#endif
            }
          else
            {
#ifdef VTK_PLOT_ASCII
              vtkwriter.write( filename.c_str(), Dune::VTK::ascii );
#else
              vtkwriter.write( filename.c_str(), Dune::VTK::appendedraw );
#endif
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
          if(found<filename.size()&& gv.comm().size() > 1)
            {
#ifdef VTK_PLOT_ASCII
              vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::ascii );
#else
              vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::appendedraw );
#endif
            }
          else
            {
#ifdef VTK_PLOT_ASCII
              vtkwriter.write( filename.c_str(), Dune::VTK::ascii );
#else
              vtkwriter.write( filename.c_str(), Dune::VTK::appendedraw );
#endif
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




      // Take care: This output will not work, if there is a mismatch between the resolution of the used grid and the virtual grid!
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

        if( gv.comm().rank()==0 
            &&
            verbosity >= VERBOSITY_DEBUG_LEVEL )
          std::cout << "VTK output of cell-data " << title.c_str() << " to file '" << filename.c_str() << ".vtu'" <<  std::endl;
  

        if(bSubSampling){

          // Adrian (2012-05-02): 
          // This will plot the cell data as they are, 
          // namely in a non-conforming way.

          const int dim = GV::dimension;
          typedef Dune::PDELab::P0LocalFiniteElementMap<CTYPE,REAL,dim> P0FEM;
          typedef Dune::PDELab::NoConstraints NOCONS;
          typedef Dune::PDELab::ISTLVectorBackend<> VBE1;
    
#ifdef USE_CUBE
          P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
#else
          P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));
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
          if(found<filename.size()&& gv.comm().size() > 1)
            {
#ifdef VTK_PLOT_ASCII
              vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::ascii );
#else
              vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::appendedraw );
#endif
            }
          else
            {
#ifdef VTK_PLOT_ASCII
              vtkwriter.write( filename.c_str(), Dune::VTK::ascii );
#else
              vtkwriter.write( filename.c_str(), Dune::VTK::appendedraw );
#endif
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


      // Version 3:
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
            if(found<filename.size()&& gv.comm().size() > 1)
              {
#ifdef VTK_PLOT_ASCII
                ss_vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::ascii );
#else
                ss_vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::appendedraw );
#endif
              }
            else
              {
#ifdef VTK_PLOT_ASCII
                ss_vtkwriter.write( filename.c_str(), Dune::VTK::ascii );
#else
                ss_vtkwriter.write( filename.c_str(), Dune::VTK::appendedraw );
#endif
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
            if(found<filename.size()&& gv.comm().size() > 1)
              {
#ifdef VTK_PLOT_ASCII
                vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::ascii );
#else
                vtkwriter.pwrite( filename.substr(found+1),filename.substr(0,found), "", Dune::VTK::appendedraw );
#endif
              }
            else
              {
#ifdef VTK_PLOT_ASCII
                vtkwriter.write( filename.c_str(), Dune::VTK::ascii );
#else
                vtkwriter.write( filename.c_str(), Dune::VTK::appendedraw );
#endif
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

      /*
       * function for writing the Y_old to vtu on Level-0 grid
       */
      template<typename GV,typename YFG>
      inline static void output_hdf5data_to_gridview
      ( const GV &gv_0,
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
        HDF5Tools::
          read_parallel_from_HDF5( gv_0
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

    }; // class VTKPlot

  }
}

#endif
