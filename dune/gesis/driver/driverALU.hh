#ifndef DUNE_GESIS_DRIVER_ALU_HH
#define DUNE_GESIS_DRIVER_ALU_HH

#ifdef HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>   // <--- This requires the dune-alugrid module.
#endif
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>   // <--- This requires the external library ALUGrid.
#endif

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/gesis/common/io/HDF5Tools.hh>

#ifdef LINE_PLOT_M0
#include "lineplot.hh"
#endif

#ifndef USE_YASP
#include <dune/grid/io/file/dgfparser.hh>
#endif

#include <dune/gesis/common/Vector.hh>
#include <dune/gesis/common/DenseMatrix.hh>

#include <dune/gesis/BVPs/wells/calibrateWellCentersToElementCenters.hh>

#include <dune/gesis/BVPs/DG/RT0flux.hh>
#include <dune/gesis/BVPs/DarcyVelocityCache.hh>
#include <dune/gesis/common/gfs_utilities.hh>
#include <dune/gesis/BVPs/GroundWaterParameters.hh>
#include <dune/gesis/BVPs/GeoElectricalPotentialParameters.hh>


#include <dune/gesis/BVPs/DG/dgf_pressurefield.hh>
#include <dune/gesis/BVPs/DG/rt0_pressurefield.hh>
#include <dune/gesis/BVPs/DG/reorderedgridview.hh>

#include <dune/gesis/BVPs/obs/MeasurementList.hh>

#include <dune/gesis/BVPs/adaptive/driver_h_adaptive_M0.hh>

#include <dune/gesis/BVPs/totalMass.hh>


extern CLogfile logger;



namespace Dune {
  namespace Gesis {

    template<typename IDT,typename YFG,typename DIR,int dim>
    REAL driverALU( IDT& inputdata,   // must be non-const because of well-calibration
                    YFG& yfg_orig,
                    DIR& dir,
                    const Dune::MPIHelper& helper )
    {
      if( helper.rank()==0 && inputdata.verbosity>0 )
        std::cout << "driverALU(...)" << std::endl;

      if( helper.rank()==0 && inputdata.verbosity>0 ){
        std::cout << std::endl;
        std::cout << "*************************" << std::endl;
        std::cout << "    CCFV / DG(" << pMAX << ") version" << std::endl;
        std::cout << "*************************" << std::endl;
        std::cout << std::endl;
      }

      Dune::Timer watch;

      typedef Dune::Gesis::FieldData FD;
      const FD fielddata(inputdata);

      // Total number of all cells required to resolve the parameter field
      UINT nAllCells = inputdata.domain_data.nCells[0] * inputdata.domain_data.nCells[1];
      if(dim==3)
        nAllCells *= inputdata.domain_data.nCells[2];

      Dune::FieldVector<CTYPE,dim> extensions(1.0);
      for(UINT i=0;i<dim;i++)
        extensions[i] = inputdata.domain_data.extensions[i];


      /*
        Wait until the new DGF file is fully written before every processor starts reading data from it!
        Otherwise, we might run into trouble because some fast process is reading old info!
      */
      logger << "---MPI_Barrier---" << std::endl;
      MPI_Barrier( helper.getCommunicator() );

      UINT baselevel = inputdata.domain_data.alugrid_baselevel;

      // make grid
      Dune::FieldVector<REAL,dim> LowerLeft(0);
      Dune::FieldVector<REAL,dim> UpperRight;
      UpperRight[0] = inputdata.domain_data.extensions[0];
      UpperRight[1] = inputdata.domain_data.extensions[1];

      Dune::array<unsigned int,dim> elements;
      elements[0] = inputdata.domain_data.nCells[0];
      elements[1] = inputdata.domain_data.nCells[1];

      UpperRight[2] = inputdata.domain_data.extensions[2];
      elements[2] = inputdata.domain_data.nCells[2];

      typedef Dune::ALUGrid<dim,dim,Dune::cube,Dune::nonconforming> GRID;

      watch.reset();
      Dune::StructuredGridFactory<GRID> structuredGridFactory;

      std::cout << "=== Dune::StructuredGridFactory.createCubeGrid() building 3D ALUGRID ..." << std::endl;

      Dune::shared_ptr<GRID> gridptr =
        structuredGridFactory.createCubeGrid(LowerLeft, UpperRight, elements);

      REAL elapsed_time;
      elapsed_time = watch.elapsed();
      std::cout << "=== Dune::StructuredGridFactory building 3D ALUGRID took "
                << elapsed_time << " sec." << std::endl;

      watch.reset();
      GRID& grid = *gridptr;

      if(grid.loadBalance()){
        elapsed_time = watch.elapsed();
        std::cout << "=== Initial load balancing for ALUGRID took "
                  << elapsed_time << " sec." << std::endl;
      }

      std::cout << "ALUGRID grid baselevel = " << grid.maxLevel() << std::endl;


      // 1.) GV and GFS for the Groundwater Equation:
      logger << "Get level gridview for the flow equation." << std::endl;
      typedef typename GRID::LevelGridView GV_GW;

      //const GV_GW &gv_gw = grid.levelView(baselevel);
      GV_GW gv_gw = grid.levelGridView(baselevel); // non-const because load-balancing may change the grid for rt0_pressurefield


      //well-calibration:
      calibrateWellCentersToElementCenters( gv_gw,
                                            inputdata // must be non-const because well data is being adjusted here!
                                            );

      inputdata.loglistOfAllWellCenters();

      //**************************************************************
      // CONSTRAINTS section:
      // Define appropriate constraints for the
      // different numerical schemes
      //
      //**************************************************************


#if USE_YASP
      typedef Dune::PDELab::P0ParallelConstraints CONSTRAINTS;
#else
      // This must be used for the parallel ALUGRID.
      // We will use an overlapping linear solver on a non-overlapping grid
      // with ghost-cells. Therefore, the ghost-cells are treated as
      // if they were overlap cells of an overlapping grid with overlap=1.
      typedef Dune::PDELab::P0ParallelGhostConstraints CONSTRAINTS;
#endif



      typedef Dune::PDELab::ISTLVectorBackend<> VBE_GW; // blocksize=1 for GWE
      typedef Dune::PDELab::P0LocalFiniteElementMap<CTYPE,REAL,dim> FEM_ELLIP;

#ifdef USE_Pk
      // 2.) Transport Equation using DG with Pk elements on cubes
      typedef Dune::PDELab::OPBLocalFiniteElementMap<CTYPE,REAL,pMAX,dim,bt> FEM_HYPER;
      const int blocksize = Dune::PB::PkSize<pMAX,dim>::value;
      logger << "DG: Using PkDG for the hyperbolic PDE with blocksize "
             << blocksize << std::endl;
#else
      // 2.) Transport Equation using DG with Qk elements on cubes
      typedef Dune::PDELab::QkDGLocalFiniteElementMap<CTYPE,REAL,pMAX,dim> FEM_HYPER;
      const int blocksize = Dune::QkStuff::QkSize<pMAX,dim>::value;
      logger << "DG: Using QkDG for the hyperbolic PDE with blocksize "
             << blocksize << std::endl;
#endif



      //**************************************************************
      // YField section:
      //
      // Sequential mode: yfg_orig contains random field data
      // Parallel mode: Retrieve random field data from HDF5 file
      //
      //**************************************************************

      // Get conductivity field: parallel fetching of random field data
      Vector<REAL> local_Yfield_vector;
      Vector<UINT> local_offset;
      Vector<UINT> local_count;

      if( ( helper.size() > 1 )
          /*&& ( inputdata.yfield_properties.variance != 0 )*/
          ) {
        //std::string YField_h5_filename = dir.yfielddir + "/Yfield.h5";


        if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_EQ_DETAILS )
          std::cout << "Start parallel fetching of Yfield data..." << std::endl;

        HDF5Tools::h5g_pRead(
                             gv_gw
                             , local_Yfield_vector
                             , dir.yfield_h5file
                             , "/YField"
                             , local_offset
                             , local_count
                             , fielddata
                             , 1 // P0 blocksize
                             , FEMType::DG // P0
                             , 0 // YField is on grid level 0.
                             );

        yfg_orig.parallel_import_from_local_vector( local_Yfield_vector, local_offset, local_count );

      }


      if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ) {
        std::cout << "yfg_orig.size() = " << yfg_orig.size() << std::endl;
        std::cout << "yfg_orig.localSize() = " << yfg_orig.localSize() << std::endl;
      }

      yfg_orig.setWellConductivities( gv_gw, inputdata );



      if( inputdata.plot_options.vtk_plot_yfield ){
        yfg_orig.plot2vtu( gv_gw, dir.Y_orig_vtu, "Y_orig", 0 );
      }


      // 1.) GFS for the Flow Equation:
      typedef Dune::PDELab::GridFunctionSpace<GV_GW,FEM_ELLIP,CONSTRAINTS,VBE_GW> GFS_GW;


      typedef MeasurementList<CPointData> MEASLIST;
      /*
       * define number of measurements but only if this type of inversion is performed (inputdata.problem_types.***_inversion) for the geo inversion.
       */
      MEASLIST orig_measurements(inputdata,inputdata.problem_types.synthetic);


      //**************************************************************
      // FEM section:
      // Define appropriate finite element maps for the
      // different numerical schemes
      //
      //**************************************************************


#ifdef USE_SIMPLEX
      const Dune::GeometryType::BasicType bt = Dune::GeometryType::simplex;
#else
      const Dune::GeometryType::BasicType bt = Dune::GeometryType::cube;
#endif
      // 1.) Groundwater Equation using CCFV
      Dune::GeometryType gt = Dune::GeometryType(bt,dim);
      FEM_ELLIP fem_ellip(gt);
      logger << "FEM: Using P0FEM for the elliptic PDE" << std::endl;
      FEM_HYPER fem_hyper;



      CONSTRAINTS con_gw;
      GFS_GW gfs_gw(gv_gw,fem_ellip,con_gw);


      Dune::Gesis::hAdaptiveLoopForM0<GRID,
                                      GFS_GW,
                                      IDT,
                                      YFG,
                                      DIR>( grid,
                                            gfs_gw,
                                            inputdata,
                                            yfg_orig,
                                            baselevel,
                                            dir,
                                            helper
                                            );

      // log the original measurements
      //orig_measurements.write_to_logger("original measurements for inversion");


      return General::reportTotalTime( helper.rank() );

    } // end of driverALU()

  } // Gesis

} // Dune
#endif // DRIVER_ALU_HH
