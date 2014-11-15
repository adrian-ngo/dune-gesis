#ifndef DUNE_GESIS_DRIVER_UG_HH
#define DUNE_GESIS_DRIVER_UG_HH


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
#include <dune/gesis/common/my_gfs_utilities.hh>
#include <dune/gesis/BVPs/GroundWaterParameters.hh>
#include <dune/gesis/BVPs/GeoElectricalPotentialParameters.hh>


#include <dune/gesis/BVPs/DG/dgf_pressurefield.hh>
#include <dune/gesis/BVPs/DG/rt0_pressurefield.hh>
#include <dune/gesis/BVPs/DG/reorderedgridview.hh>


#include <dune/gesis/BVPs/obs/MeasurementList.hh>

//#include <dune/gesis/BVPs/projectors/L2SubspaceProjector.hh>
//#include <dune/gesis/BVPs/projectors/CoarseGridP0Projector.hh>
//#include <dune/gesis/BVPs/ForwardSimulator.hh>

#if HAVE_UG
#include <dune/gesis/BVPs/adaptive/driver_h_adaptive.hh>
#endif

#include <dune/gesis/BVPs/totalMass.hh>


extern CLogfile logger;

namespace Dune {
  namespace Gesis {

    template<
      typename IDT,
      typename YFG,
      typename DIR,
      int dim
      >
    REAL driverUG( IDT& inputdata,   // must be non-const because of well-calibration
                   YFG& yfg_orig,
                   DIR& dir,
                   const Dune::MPIHelper& helper )
    {
      if( helper.rank()==0 && inputdata.verbosity>0 )
        std::cout << "driverUG(...)" << std::endl;

      if( helper.rank()==0 && inputdata.verbosity>0 ){
        std::cout << std::endl;
        std::cout << "*************************" << std::endl;
#ifdef USE_FEM
        std::cout << "    FEM / SDFEM version" << std::endl;
#else
        std::cout << "    CCFV / DG(" << pMAX << ") version" << std::endl;
#endif
        std::cout << "*************************" << std::endl;
        std::cout << std::endl;
      }

      Dune::Timer watch;

      // Total number of all cells required to resolve the parameter field
      UINT nAllCells = inputdata.domain_data.nCells[0] * inputdata.domain_data.nCells[1];
      if(dim==3)
        nAllCells *= inputdata.domain_data.nCells[2];

      Dune::FieldVector<CTYPE,dim> extensions(1.0);
      for(UINT i=0;i<dim;i++)
        extensions[i] = inputdata.domain_data.extensions[i];



#if HAVE_UG
      // If we use the UG grid, we need to first setup an appropriate base grid,
      // which has at least enough cells
      // to be distributed among the number of processes before we can start refining.
      // This is done by creating a DGF (Dune Grid File) which is then being
      // read by the Dune Grid Pointer.

      std::string sDGF_filename( "mygrid.dgf" );

      if( helper.rank() == 0 ) {
        DGFTools::
          output_Domain_DGF( sDGF_filename.c_str()
                             , inputdata.domain_data.extensions
                             , inputdata.domain_data.nCells
                             , "cube"
                             , inputdata.domain_data.ug_heapsize
                             , helper.size()
                             );
      }


      /*
        Wait until the new DGF file is fully written before every processor starts reading data from it!
        Otherwise, we might run into trouble because some fast process is reading old info!
      */
      MPI_Barrier( helper.getCommunicator() );
      typedef Dune::UGGrid<dim> GRID;
      UINT baselevel = inputdata.domain_data.ug_baselevel;
      std::cout << "Create UG grid with " << sDGF_filename.c_str() << std::endl;
      Dune::GridPtr<GRID> gridptr( sDGF_filename.c_str(), helper.getCommunicator() );
      std::cout << "Get UG grid pointer" << std::endl;
      GRID& grid = *gridptr;

      if(grid.loadBalance()){
        std::cout << "load balance successful" << std::endl;
      }

      if( baselevel==0 ){
        std::cout << "No global refinement for UG grid" << std::endl;
      }
      else{
        std::cout << "Globally refine UG grid with level " << baselevel << std::endl;
        grid.globalRefine( baselevel );
      }

      grid.setClosureType( Dune::UGGrid<dim>::NONE );  // we work with hanging nodes!
      grid.setRefinementType( Dune::UGGrid<dim>::COPY );  // Gets full gridview on all levels, BUT this works only for conforming refinement. (Oliver Sander)
      std::cout << "ug grid baselevel = " << grid.maxLevel() << std::endl;

#endif // HAVE_UG


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

#ifdef USE_FEM

#ifdef PARALLEL

#ifdef USE_NOVLP_MODE
      typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CONSTRAINTS;
#else
      typedef Dune::PDELab::OverlappingConformingDirichletConstraints CONSTRAINTS;
#endif

#else // sequential

      typedef Dune::PDELab::ConformingDirichletConstraints CONSTRAINTS;

#endif // parallel/seqential

#else // using DG

#if HAVE_UG
      // This must be used for the parallel ALUGRID (Rebecca).
      // We will use an overlapping linear solver on a non-overlapping grid
      // with ghost-cells. Therefore, the ghost-cells are treated as
      // if they were overlap cells of an overlapping grid with overlap=1.
      typedef Dune::PDELab::P0ParallelGhostConstraints CONSTRAINTS;
#else
      typedef Dune::PDELab::P0ParallelConstraints CONSTRAINTS;
#endif

#endif // USE_FEM

#ifndef NEW_GV
      //====================================================================
      // Use default (lexicographical) sorting of grid elements for SDFEM!
      // Resorting would not reduce the number of iterations for the
      // AMG based linear solver.
      //====================================================================
      typedef typename GRID::LeafGridView GV_TP;
#endif

      typedef Dune::PDELab::ISTLVectorBackend<> VBE_GW; // blocksize=1 for GWE
#ifdef USE_FEM
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV_GW,CTYPE,REAL,1> FEM_ELLIP;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV_TP,CTYPE,REAL,1> FEM_HYPER;
#else
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
      Vector<UINT> local_count;
      Vector<UINT> local_offset;


#ifdef PARALLEL

      if( ( helper.size() > 1 )
          /*&& ( inputdata.yfield_properties.variance != 0 )*/
          ) {
        //std::string YField_h5_filename = dir.yfielddir + "/Yfield.h5";


        if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_EQ_DETAILS )
          std::cout << "Start parallel fetching of Yfield data..." << std::endl;

        HDF5Tools::
          read_parallel_from_HDF5(
                                  gv_gw
                                  , inputdata
                                  , local_Yfield_vector
                                  , "/YField"
                                  , local_count
                                  , local_offset
                                  , dir.kfield_h5file
                                  , 1 // P0 blocksize
                                  , FEMType::DG // P0
                                  , 0 // YField is on grid level 0.
                                  );

        yfg_orig.parallel_import_from_local_vector( local_Yfield_vector, local_count, local_offset );

      }

#endif // #ifdef PARALLEL

      if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ) {
        std::cout << "yfg_orig.size() = " << yfg_orig.size() << std::endl;
        std::cout << "yfg_orig.localSize() = " << yfg_orig.localSize() << std::endl;
      }

      yfg_orig.setWellConductivities( gv_gw );

      if( inputdata.plot_options.vtk_plot_yfield ){
        yfg_orig.plot2vtu( gv_gw, dir.Y_orig_vtu, "Y_orig", baselevel );
      }


      // Dimensions of extended field per zone:
      std::vector< Vector<UINT> > nCellsExt;
      yfg_orig.export_nCellsExt( nCellsExt );

      typedef typename IDT::SDT SDT;

      typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD;
      // 1.) GFS for the Flow Equation:
      typedef Dune::PDELab::GridFunctionSpace<GV_GW,FEM_ELLIP,CONSTRAINTS,VBE_GW> GFS_GW;

      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;

#ifdef NEW_GV
#ifdef USE_DGF_PressureField
      logger << "Using DGF_PressureField ... " << std::endl;
      typedef DGF_PressureField<GFS_GW,VCType_GW> RT0_PF;
#else
      logger << "Using RT0_PressureField ... " << std::endl;
      typedef RT0_PressureField<GV_GW,GFS_GW,GWP_FWD> RT0_PF;
#endif

      typedef typename GRID::LeafGridView OLD_GV_TP;
      typedef ReorderedGridView
        < OLD_GV_TP
          , RT0_PF
          , PressureLikeOrdering
          > GV_TP;
      // Grid elements will be sorted according to the appropriate hydraulic head.
#endif


      // 2.) GFS for the Transport Equation:
#ifdef USE_FEM
      typedef Dune::PDELab::ISTLVectorBackend<> VBE_TP; // blocksize=1 for TPE with SDFEM
#else
      typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking,blocksize> VBE_TP;
#endif
      typedef Dune::PDELab::GridFunctionSpace<GV_TP,FEM_HYPER,CONSTRAINTS,VBE_TP> GFS_TP;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;

      // WARNING: Here, CG == DG, L2-Projection is switched OFF!
      typedef FEM_HYPER FEMCG;
      typedef VBE_TP VBE_CG;
      typedef CONSTRAINTS CGCONSTRAINTS;
      typedef Dune::PDELab::GridFunctionSpace<GV_TP,FEMCG,CGCONSTRAINTS,VBE_CG> GFS_CG;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;


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

#ifdef USE_FEM


#ifndef NEW_GV
      logger << "Get leaf gridview for the transport equation." << std::endl;
      const GV_TP &gv_tp = grid.leafGridView();
#endif

      // 1.) Groundwater Equation using Standard Galerkin
      FEM_ELLIP fem_ellip(gv_gw);
      logger << "FEM: Using Q1FEM for the elliptic PDE" << std::endl;

      // 2.) Transport Equation using Streamline Diffusion
      FEM_HYPER fem_hyper(gv_tp);
      logger << "FEM: Using Q1FEM for the hyperbolic PDE" << std::endl;

#else

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

#endif

      CONSTRAINTS con_gw;
      GFS_GW gfs_gw(gv_gw,fem_ellip,con_gw);


      // CONSTRAINTS con_tp;

#if defined USE_FEM && defined PARALLEL && defined USE_NOVLP_MODE
      con_gw.compute_ghosts( gfs_gw );
      con_tp.compute_ghosts( gfs_tp );
#endif

#if HAVE_UG
      Dune::Gesis::adaptiveDrive_h_rgv
        <GRID,
         GFS_GW,
         IDT,
         YFG,
         DIR>
        ( grid,
          gfs_gw,
          inputdata,
          yfg_orig,
          baselevel,
          dir,
          helper
          );

#endif // USE_UG

      // log the original measurements
      //orig_measurements.write_to_logger("original measurements for inversion");


      return General::reportTotalTime( helper.rank() );

    } // end of driverUG()

  } // Gesis

} // Dune
#endif // DRIVER_UG_HH
