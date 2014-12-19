#ifndef DUNE_GESIS_DRIVER_HH
#define DUNE_GESIS_DRIVER_HH

#ifdef LINE_PLOT_M0
#include <dune/gesis/common/io/lineplot.hh>
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

#include <dune/gesis/BVPs/projectors/L2SubspaceProjector.hh>
#include <dune/gesis/BVPs/projectors/CoarseGridP0Projector.hh>

#include <dune/gesis/BVPs/ForwardSimulator.hh>

#include <dune/gesis/QLGA/inversionkernel.hh>

#include <dune/gesis/yfield/generate_CR.hh>

#include <dune/gesis/BVPs/totalMass.hh>


extern CLogfile logger;



namespace Dune {
  namespace Gesis {

    template<typename IDT,typename YFG,typename DIR,int dim>
    REAL driver( IDT& inputdata,   // must be non-const because of well-calibration
                 YFG& yfg_orig,
                 DIR& dir,
                 const Dune::MPIHelper& helper )
    {
      if( helper.rank()==0 && inputdata.verbosity>0 )
        std::cout << "driver(...)" << std::endl;

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

#ifdef USE_YASP
      Dune::array<int,dim> nYaspGridCells( Dune::fill_array<int,dim>(1) );
      for (UINT i = 0; i < dim; i++)
        nYaspGridCells[i] = inputdata.domain_data.yasp_nCells[i];

      std::bitset<dim> periodic(false);

#ifdef PARALLEL
      // parallel code: define overlap for the YASP grid
#ifdef USE_NOVLP_MODE
      int overlap = 0;
#else // overlapping:
      int overlap = inputdata.domain_data.yasp_overlap;
#endif

      if( helper.rank()==0 && inputdata.verbosity>0 )
        std::cout << "Define parallel Yaspgrid with overlap " << overlap << std::endl;


#ifdef DIMENSION3
      Vector<UINT> yasp_partitions(inputdata.domain_data.nPartitions);

      typedef YaspPartition<Vector<UINT>,dim> YP;

      YP* yp = (YP*) Dune::YaspGrid<dim>::defaultLoadbalancer();

      if( yasp_partitions.componentsproduct()==0 ){
        if( helper.rank()==0 && inputdata.verbosity>0 )
          std::cout << "Using default partitioning of YASP." << std::endl;
      }
      else {

        if( (int)yasp_partitions.componentsproduct() != helper.size() ){
          // calculate the optimal partitioning
          UINT n1,n2, larger, smaller;
          General::factorize_two_optimally(UINT(helper.size()),n1,n2);

          larger=(n1<n2)?n2:n1;
          smaller=(n1<n2)?n1:n2;

          // Note: We should have as few partitions in x direction as possible if the flow goes in x direction.
          if( inputdata.domain_data.yasp_nCells[1] > larger ){
            // enough cells in y direction to take up this larger number of partitions
            yasp_partitions[0] = smaller;
            yasp_partitions[1] = larger;
            yasp_partitions[2] = 1;
          } else {
            // NOT enough cells in y direction to take up this larger number of partitions
            yasp_partitions[0] = larger;
            yasp_partitions[1] = smaller;
            yasp_partitions[2] = 1;
          }
          yp = new YP( yasp_partitions );
        }
        else {
          yp = new YP( yasp_partitions );
        }

        logger << "Partitioning of YASP: " << yasp_partitions << std::endl;
        if( helper.rank()==0 && inputdata.verbosity>0 )
          std::cout << "Partitioning of YASP: " << yasp_partitions << std::endl;

      }
#endif


      typedef Dune::YaspGrid<dim> GRID;
      GRID grid( helper.getCommunicator()
                 , extensions
                 , nYaspGridCells
                 , periodic
                 , overlap
#ifdef DIMENSION3
                 , yp
#endif
                 );


      // It is very important to keep the number of overlaps constant during refinement!
      // Otherwise, the solution of the transport problem will go wrong.
      bool keepPhysicalOverlap = false;
      grid.refineOptions( keepPhysicalOverlap );

#else // sequential:

      int overlap = 0;
      if( inputdata.verbosity>0 )
        std::cout << "Define sequential Yaspgrid..." << std::endl;
      typedef Dune::YaspGrid<dim> GRID;
      GRID grid( extensions, nYaspGridCells, periodic, overlap );

#endif // PARALLEL

      UINT baselevel = inputdata.domain_data.yasp_baselevel;
      if( baselevel > 0 ){
        std::cout << "Globally refine YASP grid with level "
                  << baselevel
                  << std::endl;
        grid.globalRefine( baselevel );
      }


      // Define root grid:
      Dune::shared_ptr<GRID> pRootGrid;
      if( helper.rank()==0 )
        pRootGrid = Dune::make_shared<GRID>( extensions, nYaspGridCells, periodic, 0 );

#endif // USE_YASP



      // 1.) GV and GFS for the Groundwater Equation:
      logger << "Get level gridview for the flow equation." << std::endl;
      typedef typename GRID::LevelGridView GV_GW;

      //const GV_GW &gv_gw = grid.levelView(baselevel);
      GV_GW gv_gw = grid.levelGridView(baselevel); // non-const because load-balancing may change the grid for rt0_pressurefield

      typedef Dune::shared_ptr<GV_GW> PGV;
      PGV pRootGridView;
      if( helper.rank()==0 )
        pRootGridView = Dune::make_shared<GV_GW>(pRootGrid->levelGridView(0));

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

      typedef Dune::PDELab::P0ParallelConstraints CONSTRAINTS;

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
        //std::string YField_h5_filename = dir.yfielddir + "/YField.h5";


        if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_EQ_DETAILS )
          std::cout << "Start parallel fetching of YField data..." << std::endl;

        HDF5Tools::h5g_pRead( gv_gw
                              , local_Yfield_vector
                              , dir.kfield_h5file
                              , "/YField"
                              , local_offset
                              , local_count
                              , inputdata
                              , 1 // P0 blocksize
                              , FEMType::DG // P0
                              , 0 // YField is on grid level 0.
                              );

        yfg_orig.parallel_import_from_local_vector( local_Yfield_vector, local_offset, local_count );

      }

#endif // #ifdef PARALLEL

      if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ) {
        std::cout << "yfg_orig.size() = " << yfg_orig.size() << std::endl;
        std::cout << "yfg_orig.localSize() = " << yfg_orig.localSize() << std::endl;
      }

      yfg_orig.setWellConductivities( gv_gw );


      /*
        if( helper.size()>1 ){
        yfg_orig.parallel_add_wells_to_hdf5( gv_gw,
        "/YField",
        inputdata.domain_data.nCells,
        inputdata.domain_data.virtual_gridsizes,
        dir.kfield_Well_h5file );
        }
      */

      if( inputdata.plot_options.vtk_plot_yfield ){
        yfg_orig.plot2vtu( gv_gw, dir.Y_orig_vtu, "YField", baselevel );
      }




      // Dimensions of extended field per zone:
      std::vector< Vector<UINT> > nCellsExt;
      yfg_orig.export_nCellsExt( nCellsExt );

      // smooth plot of original field:
      //{

      // WARNING: assume we have only one zone here!
      CovarianceMatrix<IDT> C_YY( inputdata.domain_data.nCells,
                                  nCellsExt[0],
                                  helper.getCommunicator(),
                                  inputdata );

      local_offset.resize(dim);
      local_count.resize(dim);

      C_YY.prepare( 0, local_offset, local_count, true );

      logger << "DEBUG: Read extended field " << std::endl;
      logger << "DEBUG: nCellsExt[0] = " << nCellsExt[0];
      logger << "DEBUG: local_offset = " << local_offset;
      logger << "DEBUG: local_count  = " << local_count;

      Vector<REAL> Lambdas;
      HDF5Tools::h5_pRead( Lambdas,
                           dir.EV_h5file[0],
                           "/FFT_R_YY",
                           inputdata,
                           local_offset,
                           local_count,
                           helper.getCommunicator()
                           );


      C_YY.prepare( 0, local_offset, local_count, false );

      logger << "DEBUG: Read embedded field " << std::endl;
      logger << "DEBUG: local_offset = " << local_offset;
      logger << "DEBUG: local_count  = " << local_count;

      Vector<REAL> orig_YField(1,0);
      HDF5Tools::h5_pRead( orig_YField,
                           dir.kfield_h5file,
                           "/YField",
                           inputdata,
                           local_offset,
                           local_count,
                           helper.getCommunicator()
                           );

      orig_YField -= inputdata.yfield_properties.zones[0].beta;


      Vector<REAL> smoothed1_orig_YField( orig_YField.size() );
      Vector<REAL> smoothed_orig_YField( orig_YField.size() );
      C_YY.mv( Lambdas, orig_YField, smoothed_orig_YField );

      REAL a1 = orig_YField * smoothed_orig_YField;
      a1 = gv_gw.comm().sum( a1 );
      REAL a2 = smoothed_orig_YField.two_norm_sqr();
      a2 = gv_gw.comm().sum( a2 );
      a1 /= a2;

      if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION )
        std::cout << "smooth scaling factor = " << a1 << std::endl;

      smoothed_orig_YField *= a1;
      smoothed_orig_YField += inputdata.yfield_properties.zones[0].beta;

      logger << "Parallel write to hdf5: smoothed_orig_YField." << std::endl;

      HDF5Tools::h5_pWrite( smoothed_orig_YField,
                            dir.kfield_h5file + "_smoothed",
                            "/YFieldSmoothed",
                            inputdata,
                            inputdata.domain_data.nCells,
                            local_offset,
                            local_count,
                            helper.getCommunicator()
                            );

      logger << "Parallel read of smoothed_orig_YField from hdf5 to gv_gw." << std::endl;



      HDF5Tools::h5g_pRead( gv_gw
                            , smoothed_orig_YField
                            , dir.kfield_h5file + "_smoothed"
                            , "/YFieldSmoothed"
                            , local_offset
                            , local_count
                            , inputdata
                            , 1 // P0 blocksize
                            , FEMType::DG // P0
                            , 0 // structure is on grid level 0
                            );
      //}


      if( inputdata.plot_options.vtk_plot_y_smooth ) {
        YFG yfg_smoothed_Yorig(inputdata,dir,helper.getCommunicator());
        if( helper.size() > 1 )
          yfg_smoothed_Yorig.parallel_import_from_local_vector( smoothed_orig_YField, local_offset, local_count );
        else
          yfg_smoothed_Yorig.import_from_vector( smoothed_orig_YField );

        yfg_smoothed_Yorig.plot2vtu( gv_gw,
                                     dir.Y_orig_vtu + "_smoothed",
                                     "YField",
                                     baselevel );
      }


      typedef typename IDT::SDT SDT; // Setup Data Type

      typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD;
      // 1.) GFS for the Flow Equation:
      typedef Dune::PDELab::GridFunctionSpace<GV_GW,FEM_ELLIP,CONSTRAINTS,VBE_GW> GFS_GW;

      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;

#ifdef NEW_GV
      // Grid elements will be sorted according to the appropriate hydraulic head.
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
#endif


      // 2.) GFS for the Transport Equation:
#ifdef USE_FEM
      typedef Dune::PDELab::ISTLVectorBackend<> VBE_TP; // blocksize=1 for TPE with SDFEM
#else
      typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking,blocksize> VBE_TP;
#endif
      typedef Dune::PDELab::GridFunctionSpace<GV_TP,FEM_HYPER,CONSTRAINTS,VBE_TP> GFS_TP;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;


      // L2-Projection of DG Solution to the CG space
#if defined L2ProjectionOfM0
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV_TP,CTYPE,REAL,1> FEMCG;
      typedef Dune::PDELab::ISTLVectorBackend<> VBE_CG;
      typedef Dune::PDELab::OverlappingConformingDirichletConstraints CGCONSTRAINTS;
      typedef Dune::PDELab::GridFunctionSpace<GV_TP,FEMCG,CGCONSTRAINTS,VBE_CG> GFS_CG;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;
#else
      // WARNING: Here, CG == DG, L2-Projection is switched OFF!
      typedef FEM_HYPER FEMCG;
      typedef VBE_TP VBE_CG;
      typedef CONSTRAINTS CGCONSTRAINTS;
      typedef Dune::PDELab::GridFunctionSpace<GV_TP,FEMCG,CGCONSTRAINTS,VBE_CG> GFS_CG;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;
#endif


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


      CONSTRAINTS con_tp;

#if defined USE_FEM && defined PARALLEL && defined USE_NOVLP_MODE
      con_gw.compute_ghosts( gfs_gw );
      con_tp.compute_ghosts( gfs_tp );
#endif


#ifdef USE_YASP

      if( inputdata.problem_types.estimation_variance_only ){

        // Read JQ files and the cokriging matrix and compute the estimation variance directly!
        GV_GW gv_0 = grid.levelGridView(0);
        watch.reset();

        Vector<REAL> EstimatedVariances;

        estimation_variance( gv_0, inputdata, dir,
                             EstimatedVariances );

        General::log_elapsed_time( watch.elapsed(),
                                   gv_0.comm(),
                                   inputdata.verbosity,
                                   "UQ",
                                   "estimation_variance: total time computing sigma2"
                                   );

        // Report to standard output.
        if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
          std::cout << "inversionkernel: Estimation variance stored in '"
                    << dir.estimation_variance_h5file
                    << "'. All JQ h5-files of the very last iteration are removed from the buffer directory." << std::endl;
        }

        return General::reportTotalTime( helper.rank() );

      }


#endif // USE_YASP




      //*******************************************************************
      // GFS section:
      //
      // Make a grid function space for each of the equations to be solved
      //
      //*******************************************************************


      typedef GradientVectorField<GWP_FWD,GFS_GW> DARCY_FLUX_BASE;
      typedef DarcyVelocityCache<GWP_FWD,GFS_GW,DARCY_FLUX_BASE> DARCY_FLUX_DGF;

      /*
       * take measurements of lnK (if needed)
       */
      if( inputdata.problem_types.synthetic && inputdata.problem_types.lnK_inversion){
        orig_measurements.take_measurements_lnK( yfg_orig,
                                                 gv_gw,
                                                 helper );
      }




      /*
       * needed precalculation for geoelectrical potential calculation!
       */
      if( inputdata.problem_types.moments_geoeletric_potential_inversion ||
          inputdata.problem_types.moments_geoeletric_potential_forward ) {
        /*
         * IMPORTANT:
         *
         * Here is something to be done if really the hdf5 files for each setup are given!
         *
         *
         */

        // generate a field on one processor
        if(helper.rank()==0){
          Vector<REAL> logsigma0(nAllCells,0.0);
          Vector<REAL> logkappa(nAllCells,0.0);
          std::vector<Vector<REAL>>  X(inputdata.yfield_properties.nz); // zonation matrix
          // read the needed data from the HDF5 file -> and again for the different MPI Pools
          for(UINT jj=0; jj<inputdata.yfield_properties.nz; jj++)
            HDF5Tools::h5_Read(X[jj],dir.zonation_matrix[jj],"/X");

          for(UINT jj=0; jj<inputdata.yfield_properties.nz; jj++)
            for(UINT ii=0; ii<nAllCells; ii++){
              logsigma0[ii] += X[jj][ii]*inputdata.yfield_properties.zones[jj].sigma0;
              logkappa[ii]  += X[jj][ii]*inputdata.yfield_properties.zones[jj].kappa;
            }

          for(UINT ii=0; ii<nAllCells; ii++){
            logsigma0[ii]=log(logsigma0[ii]);
            logkappa[ii]=log(logkappa[ii]);
          }

          HDF5Tools::h5g_Write( logsigma0
                                , dir.logsigma0_h5file
                                , "/logsigma0"
                                , inputdata
                                );

          HDF5Tools::h5g_Write(
                               logkappa
                               , dir.logkappa_h5file
                               , "/logkappa"
                               , inputdata
                               );
        }
      }



      //orig_measurements.write_to_logger("measurements for inversion");

      YFG log_electricalConductivity(inputdata,dir,helper.getCommunicator());
      YFG log_kappafield(inputdata,dir,helper.getCommunicator());

      if( inputdata.problem_types.synthetic
          && ( inputdata.problem_types.head_forward ||
               inputdata.problem_types.head_inversion ||
               inputdata.problem_types.transport_forward ||
               inputdata.problem_types.transport_inversion_m0 ||
               inputdata.problem_types.transport_inversion_m1 ||
               inputdata.problem_types.heat_forward ||
               inputdata.problem_types.heat_mean_arrival_time_inversion ||
               inputdata.problem_types.moments_geoeletric_potential_forward ||
               inputdata.problem_types.moments_geoeletric_potential_inversion ) ) {


        // Generate the log electrical Conductivity field and the kappa field.
        // Read the field from hdf5!
        if(inputdata.problem_types.moments_geoeletric_potential_forward ||
           inputdata.problem_types.moments_geoeletric_potential_inversion ){

          Vector<REAL> logsigma0;

          HDF5Tools::h5g_pRead( gv_gw
                                , logsigma0
                                , dir.logsigma0_h5file
                                , "/logsigma0"
                                , local_offset
                                , local_count
                                , inputdata
                                , 1 // P0 blocksize
                                , FEMType::DG // P0
                                , 0 // structure is on grid level 0
                                );




          // export the data to the sigma0 Generator
          if(helper.size()>1)
            log_electricalConductivity.parallel_import_from_local_vector(logsigma0, local_offset, local_count );
          else
            log_electricalConductivity.import_from_vector( logsigma0 );


          Vector<REAL> logkappa;

          HDF5Tools::h5g_pRead( gv_gw
                                , logkappa
                                , dir.logkappa_h5file
                                , "/logkappa"
                                , local_offset
                                , local_count
                                , inputdata
                                , 1 // P0 blocksize
                                , FEMType::DG // P0
                                , 0 // structure is on grid level 0
                                );

          // export the data to the sigma0 Generator
          if(helper.size()>1)
            log_kappafield.parallel_import_from_local_vector(logkappa, local_offset, local_count );
          else
            log_kappafield.import_from_vector( logkappa );

        }


        typedef GEPotentialForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GEP_FWD;

        // ***********************************************************************
        // BIG LOOP over setups:
        //
        // ***********************************************************************

        for( UINT iSetup=0; iSetup<inputdata.setups.size(); iSetup++ ) {


          const SDT& setupdata = inputdata.setups[iSetup];


          GWP_FWD gwp_fwd( inputdata, setupdata,
                           yfg_orig
                           );
          GEP_FWD gep_fwd( inputdata, setupdata, log_electricalConductivity );

          GEP_FWD log_kappawrapper_gw( inputdata, setupdata, log_kappafield ); // equation (66), (147): kappa

          typedef ForwardSimulator<GWP_FWD,
                                   GEP_FWD,
                                   GFS_GW,
                                   GFS_TP,
                                   GFS_CG,
                                   IDT,
                                   SDT,
                                   DIR> FWD_SIM;


          FWD_SIM fwd_sim( gwp_fwd,
                           gep_fwd,
                           log_kappawrapper_gw,
                           gfs_gw,
                           inputdata,
                           setupdata,
                           dir);

          if(helper.rank()==0){
            std::cout << General::getDateAndTime()
                      << "Forward simulation of original measurements for setup #"
                      << iSetup+1 << std::endl;
          }
          // Allocate vectors to hold head, m0 and m1:

          VCType_GW vchead_orig( gfs_gw, 0.0 );

          if( !General::bFileExists( dir.vchead_orig_h5file[iSetup] ) ){

            //solve for forward head
            fwd_sim.head_simulation( vchead_orig );

            // store result in hdf5
            HDF5Tools::write_BackendVector_to_HDF5<GFS_GW>( gfs_gw
                                                            , inputdata
                                                            , dir.vchead_orig_h5file[iSetup]
                                                            , "/vchead_orig"
                                                            , vchead_orig
                                                            , true // PRESERVE_STRUCTURE
                                                            );


          }
          else{

            // read result from hdf5
            HDF5Tools
              ::read_BackendVector_from_HDF5<GFS_GW>( gfs_gw
                                                      , inputdata
                                                      , dir.vchead_orig_h5file[iSetup]
                                                      , "/vchead_orig"
                                                      , vchead_orig
                                                      , true // preserve_structure!
                                                      );
          }

          if( inputdata.problem_types.head_forward ||
              inputdata.problem_types.head_inversion ){

            std::stringstream datafile_head_measurements;
            datafile_head_measurements << dir.datadimdir << "/head_orig.meas.iSetup" << iSetup;
            bool bReadheadMeasurementDone = false;
            if( !inputdata.problem_types.new_YField ){
              if( General::bFileExists( datafile_head_measurements.str() ) ){
                // Read the measurements from another file!
                bReadheadMeasurementDone = orig_measurements.read_measurements( 1, iSetup, datafile_head_measurements.str(), helper );
              }
            } else {
              std::cout << "Generate new head measurement data for new parameter field." << std::endl;
            }

            

            if( bReadheadMeasurementDone ){
              if(helper.rank()==0){
                std::cout << "Reading head measurements from "
                          << datafile_head_measurements.str()
                          << " was successful." << std::endl;
              }
            }
            else{
              //Take the head measurements!
              if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION )
                std::cout << "Take new measurements of head." << std::endl;
              orig_measurements.take_measurements( 1, vchead_orig, gfs_gw, helper, iSetup);
            }


            std::stringstream bufferfile_head_measurements;
            bufferfile_head_measurements << dir.bufferdimdir
                                         << "/head_orig.meas.iSetup" << iSetup
                                         << ".rank" << helper.rank();

            orig_measurements.store_measurements( 1, iSetup, bufferfile_head_measurements.str() );

            if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
              std::cout << "Head measurements stored to "
                        << bufferfile_head_measurements.str() << " etc." << std::endl;
              // copy file to DATA:
              std::stringstream datafile_per_Setup;
              datafile_per_Setup << dir.datadimdir
                                 << "/head_orig.meas.iSetup" << iSetup;

              if( inputdata.problem_types.generate_measurement_data ||
                  inputdata.problem_types.new_YField ){
                General::systemCall( "cp -v " + bufferfile_head_measurements.str() + " " + datafile_per_Setup.str() );
              }
            }

          }

          if(inputdata.plot_options.vtk_plot_head ){
            VTKPlot::output2vtu( gfs_gw,
                                 vchead_orig,
                                 dir.head_orig_vtu[iSetup],
                                 "h_orig",
                                 inputdata.verbosity,
                                 true,
                                 0
                                 );
          }


          // Calculate the Darcy flux out of the head "vchead_orig"
          DARCY_FLUX_DGF darcyflux_dgf( gwp_fwd,
                                        gfs_gw,
                                        vchead_orig,
                                        baselevel );

          if( inputdata.plot_options.vtk_plot_q ){
            VTKPlot::output_dgf_to_vtu( gv_gw,
                                        gfs_gw,
                                        darcyflux_dgf.exportDGF(),
                                        dir.q_orig_vtu[iSetup],
                                        "q_orig",
                                        inputdata.verbosity,
                                        true,
                                        0 );
          }


#ifdef USE_YASP
          UINT maxlevel=inputdata.domain_data.yasp_maxlevel;
#endif

          if(maxlevel>baselevel)
            grid.globalRefine( maxlevel - baselevel );

          if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_DEBUG_LEVEL ){
            std::cout << "YASP grid baselevel = " << baselevel << std::endl;
            std::cout << "YASP grid leaflevel = " << maxlevel << std::endl;
          }




#ifdef NEW_GV

#ifdef USE_DGF_PressureField
          Dune::shared_ptr<RT0_PF> rt0_pressurefield
            = Dune::make_shared<RT0_PF>(gfs_gw,vchead_orig,baselevel);
#else

          Dune::shared_ptr<RT0_PF> rt0_pressurefield
            = Dune::make_shared<RT0_PF>( gv_gw,
                                         gwp_fwd,
                                         baselevel );

          rt0_pressurefield->assemble(gfs_gw,vchead_orig);
#endif

          //ReversePressureLikeOrdering comparefunctor;
          PressureLikeOrdering comparefunctor;

          if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_EQ_SUMMARY ){
            std::cout << "driver: Reorder gridview" << std::endl;
          }

          GV_TP gv_tp( grid.leafGridView(),
                       rt0_pressurefield,
                       comparefunctor );

          verifyCalibratedWellsOnRefinedGrid(gv_tp,inputdata);


          if( inputdata.plot_options.vtk_plot_element_ordering ){
            std::stringstream vtu_gv_tp;
            vtu_gv_tp << dir.vtudir << "/gv_tp_orig_" << iSetup;
            VTKPlot::outputGridviewIndexToDGF( gv_tp, vtu_gv_tp.str() );
          }


#endif


          // L2-Projection of DG Solution to the CG space
          FEMCG femcg(gv_tp);
          CGCONSTRAINTS cgcons;

          GFS_TP gfs_tp(gv_tp,fem_hyper,con_tp);
          GFS_CG gfs_cg(gv_tp,femcg,cgcons );


          // L2SubspaceProjector<GFS_TP,GFS_CG,IDT,dim> L2sp(gfs_tp,gfs_cg,inputdata);
          //logger << "gfs_cg.maxLocalSize() = " << gfs_cg.maxLocalSize() << std::endl;
          //logger << "gfs_cg.size() = " << gfs_cg.size() << std::endl;

          if( inputdata.problem_types.transport_forward ||
              inputdata.problem_types.transport_inversion_m0 ||
              inputdata.problem_types.transport_inversion_m1 ||
              inputdata.problem_types.moments_geoeletric_potential_forward ||
              inputdata.problem_types.moments_geoeletric_potential_inversion ){


            std::stringstream datafile_M0_measurements;
            datafile_M0_measurements << dir.datadimdir << "/m0_orig.meas.iSetup" << iSetup;
            bool bReadM0MeasurementDone = false;
            if( !inputdata.problem_types.new_YField ){
              if( General::bFileExists( datafile_M0_measurements.str() ) ){
                // Read the measurements from another file!
                bReadM0MeasurementDone = orig_measurements.read_measurements( 2, iSetup, datafile_M0_measurements.str(), helper );
              }
            } else {
              std::cout << "Generate new m0 measurement data for new parameter field." << std::endl;
            }

            std::stringstream datafile_M1_measurements;
            datafile_M1_measurements << dir.datadimdir << "/m1_orig.meas.iSetup" << iSetup;
            bool bReadM1MeasurementDone = false;
            if( !inputdata.problem_types.new_YField ){
              if( General::bFileExists( datafile_M1_measurements.str() ) ){
                // Read the measurements from another file!
                bReadM1MeasurementDone = orig_measurements.read_measurements( 3, iSetup, datafile_M1_measurements.str(), helper );
              }
            } else {
              std::cout << "Generate new m1 measurement data for new parameter field." << std::endl;
            }

            VCType_TP vcM0_orig( gfs_tp, 0.0 );
            VCType_CG vcM0_cg( gfs_cg, 0.0 );

            VCType_TP vcM1_orig( gfs_tp, 0.0 );
            VCType_CG vcM1_cg( gfs_cg, 0.0 );


            if( !bReadM0MeasurementDone || !bReadM1MeasurementDone ){

              if(!General::bFileExists( dir.vcM0_cg_orig_h5file[iSetup] ) ){

                fwd_sim.soluteM0_simulation( darcyflux_dgf,
                                             gfs_tp,
                                             gfs_cg,
                                             vcM0_orig,
                                             vcM0_cg );

                HDF5Tools
                  ::write_BackendVector_to_HDF5<GFS_CG>( gfs_cg
                                                         , inputdata
                                                         , dir.vcM0_cg_orig_h5file[iSetup]
                                                         , "/vcM0_cg"
                                                         , vcM0_cg
                                                         , true // preserve_structure!
                                                         );

              }else{
                HDF5Tools
                  ::read_BackendVector_from_HDF5<GFS_CG>( gfs_cg
                                                          , inputdata
                                                          , dir.vcM0_cg_orig_h5file[iSetup]
                                                          , "/vcM0_cg"
                                                          , vcM0_cg
                                                          , true // preserve_structure!
                                                          );
              }

              if(inputdata.plot_options.vtk_plot_m0 ){
                VTKPlot::output2vtu( gfs_cg,
                                     vcM0_cg,
                                     dir.m0_orig_vtu[iSetup] + "_cg",
                                     "m0_orig_cg",
                                     inputdata.verbosity,
                                     true,
                                     0 );

                VTKPlot::output2vtu( gfs_tp,
                                     vcM0_orig,
                                     dir.m0_orig_vtu[iSetup],
                                     "m0_orig",
                                     inputdata.verbosity,
                                     true,
                                     std::max(0,pMAX-1) );
              }


              REAL m0dg_negMass(0), m0dg_posMass(0), m0dg_totMass(0);
              GridFunctionTools::totalMass( gfs_tp, vcM0_orig, m0dg_negMass, m0dg_posMass, m0dg_totMass );
              if( helper.rank() == 0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
                std::cout << "pos./neg. mass (m0) = "
                          << std::fixed << std::setprecision(5)
                          << m0dg_posMass << " " << m0dg_negMass  << std::endl;
                std::cout << "total mass (m0) = "
                          << std::fixed << std::setprecision(5)
                          << m0dg_totMass << std::endl;
              }

              REAL m0cg_negMass(0), m0cg_posMass(0), m0cg_totMass(0);
              GridFunctionTools::totalMass( gfs_cg, vcM0_cg, m0cg_negMass, m0cg_posMass, m0cg_totMass );
              if( helper.rank() == 0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
                std::cout << "pos./neg. mass (m0-L2) = "
                          << std::fixed << std::setprecision(5)
                          << m0cg_posMass << " " << m0cg_negMass << std::endl;
                std::cout << "total mass (m0-L2) = "
                          << std::fixed << std::setprecision(5)
                          << m0cg_totMass << std::endl;
              }


#ifdef LINE_PLOT_M0
              // *********************************
              // Plot over line for gnuplot:
              // *********************************
              std::stringstream gnuplot_datfile;
              gnuplot_datfile << dir.m0_orig_vtu[iSetup] + "_gnu.dat";

              LinePlot<GV_TP> lplot(gv_tp);
              lplot.setEndPoints( inputdata.domain_data.lineplot.startpoint,
                                  inputdata.domain_data.lineplot.endpoint );
              lplot.addDGData( gfs_tp,
                               vcM0_orig,
                               "m0_orig" );
#ifdef L2ProjectionOfM0
              lplot.addDGData( gfs_cg,
                               vcM0_cg,
                               "m0_orig_cg" );
#endif
              lplot.write(gnuplot_datfile.str());
#endif // LINE_PLOT_M0


            }

            if( bReadM0MeasurementDone ){
              if(helper.rank()==0){
                std::cout << "Reading measurements of m0 from DATA was successful." << std::endl;
              }
            }
            else{
              //Take the measurements!
              if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION )
                std::cout << "Take new measurements of m0." << std::endl;
              orig_measurements.take_measurements( 2, vcM0_cg, gfs_cg, helper, iSetup);
            }


            std::stringstream bufferfile_M0_measurements;
            bufferfile_M0_measurements << dir.bufferdimdir << "/m0_orig.meas.iSetup" << iSetup << ".rank" << helper.rank();

            orig_measurements.store_measurements( 2, iSetup, bufferfile_M0_measurements.str() );

            if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
              std::cout << "m0 measurements stored to "
                        << bufferfile_M0_measurements.str() << " etc." << std::endl;
              // copy file to DATA:
              std::stringstream datafile_per_Setup;
              datafile_per_Setup << dir.datadimdir
                                 << "/m0_orig.meas.iSetup" << iSetup;
              if( inputdata.problem_types.generate_measurement_data ||
                  inputdata.problem_types.new_YField ){
                General::systemCall( "cp -v " + bufferfile_M0_measurements.str() + " " + datafile_per_Setup.str() );
              }
            }


            if( inputdata.problem_types.transport_forward ||
                inputdata.problem_types.transport_inversion_m1 ||
                inputdata.problem_types.moments_geoeletric_potential_forward ||
                inputdata.problem_types.moments_geoeletric_potential_inversion ){

              if( bReadM1MeasurementDone ){

                if(helper.rank()==0){
                  std::cout << "Reading measurements of m1 from DATA was successful." << std::endl;
                }

              }
              else{

                if(!General::bFileExists( dir.vcM1_cg_orig_h5file[iSetup] ) ){

                  //solve solute M1 forward
                  fwd_sim.soluteM1_simulation( darcyflux_dgf,
                                               gfs_tp,
                                               gfs_cg,
                                               vcM0_cg,   // input
                                               vcM1_orig, // output
                                               vcM1_cg    // output
                                               );

                  HDF5Tools
                    ::write_BackendVector_to_HDF5<GFS_CG>( gfs_cg
                                                           , inputdata
                                                           , dir.vcM1_cg_orig_h5file[iSetup]
                                                           , "/vcM1_cg"
                                                           , vcM1_cg
                                                           , true // preserve_structure!
                                                           );
                }else{
                  HDF5Tools
                    ::read_BackendVector_from_HDF5<GFS_CG>( gfs_cg
                                                            , inputdata
                                                            , dir.vcM1_cg_orig_h5file[iSetup]
                                                            , "/vcM1_cg"
                                                            , vcM1_cg
                                                            , true // preserve_structure!
                                                            );
                }


                if(inputdata.plot_options.vtk_plot_m1 ){
                  VTKPlot::output2vtu( gfs_cg,
                                       vcM1_cg,
                                       dir.m1_orig_vtu[iSetup] + "_cg",
                                       "m1_orig_cg",
                                       inputdata.verbosity,
                                       true,
                                       0 );

                  VTKPlot::output2vtu( gfs_tp,
                                       vcM1_orig,
                                       dir.m1_orig_vtu[iSetup],
                                       "m1_orig",
                                       inputdata.verbosity,
                                       true,
                                       std::max(0,pMAX-1) );
                }

                REAL m1dg_negMass(0), m1dg_posMass(0), m1dg_totMass(0);
                GridFunctionTools::totalMass( gfs_tp, vcM1_orig, m1dg_negMass, m1dg_posMass, m1dg_totMass );
                if( helper.rank() == 0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
                  std::cout << "pos./neg. mass (m1) = "
                            << std::scientific << std::setprecision(5)
                            << m1dg_posMass << " " << m1dg_negMass  << std::endl;
                  std::cout << "total mass (m1) = "
                            << std::scientific << std::setprecision(5)
                            << m1dg_totMass << std::endl;
                }

                REAL m1cg_negMass(0), m1cg_posMass(0), m1cg_totMass(0);
                GridFunctionTools::totalMass( gfs_cg, vcM1_cg, m1cg_negMass, m1cg_posMass, m1cg_totMass );
                if( helper.rank() == 0 && inputdata.verbosity>=VERBOSITY_EQ_SUMMARY ){
                  std::cout << "pos./neg. mass (m1-L2) = "
                            << std::scientific << std::setprecision(5)
                            << m1cg_posMass << " " << m1cg_negMass << std::endl;
                  std::cout << "total mass (m1-L2) = "
                            << std::scientific << std::setprecision(5)
                            << m1cg_totMass << std::endl;
                }

                //Take the measurements!
                if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION )
                  std::cout << "Take new measurements of m1." << std::endl;
                orig_measurements.take_measurements( 3, vcM1_cg, gfs_cg, helper, iSetup);

              }

              std::stringstream bufferfile_M1_measurements;
              bufferfile_M1_measurements << dir.bufferdimdir << "/m1_orig.meas.iSetup" << iSetup << ".rank" << helper.rank();

              if( helper.rank() == 0 )
                std::cout << "Absolute m1 error for Setup #" << iSetup << " is "
                          << orig_measurements.calibrateAbsErrorForM1( iSetup )
                          << std::endl;

              orig_measurements.store_measurements( 3, iSetup, bufferfile_M1_measurements.str() );


              if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
                std::cout << "m1 measurements stored to "
                          << bufferfile_M1_measurements.str() << " etc." << std::endl;
                // copy file to DATA:
                std::stringstream datafile_per_Setup;
                datafile_per_Setup << dir.datadimdir
                                   << "/m1_orig.meas.iSetup" << iSetup;
                if( inputdata.problem_types.generate_measurement_data ||
                    inputdata.problem_types.new_YField ){
                  General::systemCall( "cp -v " + bufferfile_M1_measurements.str() + " " + datafile_per_Setup.str() );
                }
              }


              // geoelectrical potential calculations!
              if( inputdata.problem_types.moments_geoeletric_potential_forward ||
                  inputdata.problem_types.moments_geoeletric_potential_inversion ){

                UINT n_GP_config=inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig;

                for(UINT iconfig=0; iconfig<n_GP_config; iconfig++){

                  VCType_GW vcphi0_orig(gfs_gw,0.0);
                  VCType_GW vcM0phi_orig(gfs_gw,0.0);
                  VCType_GW vcM1phi_orig(gfs_gw,0.0);
                  VCType_GW vcATphi_orig(gfs_gw,0.0);

                  fwd_sim.GE_phi0_simulation( iconfig,
                                              vcphi0_orig );

                  fwd_sim.GE_simulation( gfs_cg,
                                         vcM0_cg,
                                         vcM1_cg,
                                         vcphi0_orig,
                                         vcM0phi_orig,  // output
                                         vcM1phi_orig   // output
                                         );

                  //calculate arrival time!
                  for(UINT ii=0; ii<gfs_gw.size();ii++){
                    //vcATphi_orig[ii]=vcM1phi_orig[ii]/vcM0phi_orig[ii];
                    vcATphi_orig.base()[ii] = vcM1phi_orig.base()[ii] / vcM0phi_orig.base()[ii];
                  }

                  HDF5Tools::write_BackendVector_to_HDF5<GFS_GW>( gfs_gw
                                                                  , inputdata
                                                                  , dir.vcphi0_orig_h5file[iSetup][iconfig]
                                                                  , "/vcphi0_orig"
                                                                  , vcphi0_orig
                                                                  );
                  HDF5Tools::write_BackendVector_to_HDF5<GFS_GW>( gfs_gw,
                                                                  inputdata,
                                                                  dir.vcM0phi_orig_h5file[iSetup][iconfig],
                                                                  "/vcM0phi_orig",
                                                                  vcM0phi_orig
                                                                  );
                  HDF5Tools::write_BackendVector_to_HDF5<GFS_GW>( gfs_gw,
                                                                  inputdata,
                                                                  dir.vcM1phi_orig_h5file[iSetup][iconfig],
                                                                  "/vcM1phi_orig",
                                                                  vcM1phi_orig
                                                                  );
                  HDF5Tools::write_BackendVector_to_HDF5<GFS_GW>( gfs_gw,
                                                                  inputdata,
                                                                  dir.vcATphi_orig_h5file[iSetup][iconfig],
                                                                  "/vcATphi_orig",
                                                                  vcATphi_orig
                                                                  );


                  if( inputdata.plot_options.vtk_plot_el_potential_field ){
                    VTKPlot::output2vtu( gfs_gw,
                                         vcphi0_orig,
                                         dir.phi0_orig_vtu[iSetup][iconfig],
                                         "phi0_orig",
                                         inputdata.verbosity,
                                         true,
                                         0 );
                    VTKPlot::output2vtu( gfs_gw,
                                         vcM0phi_orig,
                                         dir.M0phi_orig_vtu[iSetup][iconfig],
                                         "M0phi_orig",
                                         inputdata.verbosity,
                                         true,
                                         0 );
                    VTKPlot::output2vtu( gfs_gw,
                                         vcM1phi_orig,
                                         dir.M1phi_orig_vtu[iSetup][iconfig],
                                         "M1phi_orig",
                                         inputdata.verbosity,
                                         true,
                                         0 );
                    VTKPlot::output2vtu( gfs_gw,
                                         vcATphi_orig,
                                         dir.ATphi_orig_vtu[iSetup][iconfig],
                                         "ATphi_orig",
                                         inputdata.verbosity,
                                         true,
                                         0 );
                  }


                  if(inputdata.problem_types.moments_geoeletric_potential_inversion){
                    //Take measurements !!!
                    logger << "GE: take_measurements_AT of original field" << std::endl;
                    orig_measurements.take_measurements_AT( 5,
                                                            vcM0phi_orig,
                                                            vcM1phi_orig,
                                                            gfs_gw,
                                                            helper,
                                                            iSetup,
                                                            iconfig );
                  }
                }
              }  // END if(inputdata.problem_types.moments_geoeletric_potential_forward || inputdata.problem_types.moments_geoeletric_potential_inversion)

            }  // END if (inputdata.problem_types.transport_forward || inputdata.problem_types.transport_inversion_m1|| inputdata.problem_types.moments_geoeletric_potential_forward || inputdata.problem_types.moments_geoeletric_potential_inversion)

          }
          // END if (inputdata.problem_types.transport_forward
          //         || inputdata.problem_types.transport_inversion_m0
          //         || inputdata.problem_types.transport_inversion_m1
          //         || inputdata.problem_types.moments_geoeletric_potential_forward
          //         || inputdata.problem_types.moments_geoeletric_potential_inversion )


          if( inputdata.problem_types.heat_forward || inputdata.problem_types.heat_mean_arrival_time_inversion ){

            VCType_TP vcheatM0_orig( gfs_tp, 0.0 );
            VCType_CG vcheatM0_cg( gfs_cg, 0.0 );
            VCType_TP vcheatM1_orig( gfs_tp, 0.0 );
            VCType_CG vcheatM1_cg( gfs_cg, 0.0 );

            VCType_CG vcheatAT_orig( gfs_cg, 0.0 );


            if(!General::bFileExists( dir.vcheatM0_cg_orig_h5file[iSetup] ) ||
               !General::bFileExists( dir.vcheatM1_cg_orig_h5file[iSetup] ) ||
               !General::bFileExists( dir.vcheatArrival_Time_orig_h5file[iSetup] ) ){

              fwd_sim.heat_simulation( darcyflux_dgf,
                                       gfs_tp,
                                       gfs_cg,
                                       vcheatM0_orig,  // unused output
                                       vcheatM1_orig,  // unused output
                                       vcheatM0_cg,    // output
                                       vcheatM1_cg     // output
                                       );


#if defined L2ProjectionOfM0 || defined USE_FEM
              //calculate the original arrival time!
              for(UINT ii=0; ii < gfs_cg.size();ii++){
                if( vcheatM1_cg.base()[ii] < GEO_EPSILON*0.5)
                  vcheatAT_orig.base()[ii] = 0.0;
                else if (vcheatM0_cg.base()[ii] < GEO_EPSILON*0.5)
                  vcheatAT_orig.base()[ii] = 0.0;
                else if(vcheatM0_cg.base()[ii] < GEO_EPSILON*0.5)
                  vcheatAT_orig.base()[ii] = vcheatM1_cg.base()[ii]/(GEO_EPSILON*0.5);
                else
                  vcheatAT_orig.base()[ii] = vcheatM1_cg.base()[ii]/vcheatM0_cg.base()[ii];
              }
#else
              // Careful! CG == DG
              for(UINT ii=0; ii < vcheatM1_cg.flatsize();ii++){
                if(vcheatM1_cg.base()[ii/blocksize].base()[ii%blocksize] < GEO_EPSILON*0.5)
                  vcheatAT_orig.base()[ii/blocksize].base()[ii%blocksize] = 0.0;
                else if (vcheatM0_cg.base()[ii/blocksize].base()[ii%blocksize] < GEO_EPSILON*0.5)
                  vcheatAT_orig.base()[ii/blocksize].base()[ii%blocksize] = 0.0;
                else if(vcheatM0_cg.base()[ii/blocksize].base()[ii%blocksize] < GEO_EPSILON*0.5)
                  vcheatAT_orig.base()[ii/blocksize].base()[ii%blocksize] = vcheatM1_cg.base()[ii/blocksize].base()[ii%blocksize]/(GEO_EPSILON*0.5);
                else
                  vcheatAT_orig.base()[ii/blocksize].base()[ii%blocksize] = vcheatM1_cg.base()[ii/blocksize].base()[ii%blocksize]/vcheatM0_cg.base()[ii/blocksize].base()[ii%blocksize];
              }
#endif

              HDF5Tools
                ::write_BackendVector_to_HDF5<GFS_CG>( gfs_cg
                                                       , inputdata
                                                       , dir.vcheatM0_cg_orig_h5file[iSetup]
                                                       , "/vcheatM0_cg"
                                                       , vcheatM0_cg
                                                       );

              HDF5Tools
                ::write_BackendVector_to_HDF5<GFS_CG>( gfs_cg
                                                       , inputdata
                                                       , dir.vcheatM1_cg_orig_h5file[iSetup]
                                                       , "/vcheatM1_cg"
                                                       , vcheatM1_cg
                                                       );


              HDF5Tools
                ::write_BackendVector_to_HDF5<GFS_CG>( gfs_cg
                                                       , inputdata
                                                       , dir.vcheatArrival_Time_orig_h5file[iSetup]
                                                       , "/vcheatArrival_Time_orig"
                                                       , vcheatAT_orig
                                                       );

            }else{
              HDF5Tools
                ::read_BackendVector_from_HDF5<GFS_CG>( gfs_cg
                                                        , inputdata
                                                        , dir.vcheatM0_cg_orig_h5file[iSetup]
                                                        , "/vcheatM0_cg"
                                                        , vcheatM0_cg
                                                        );

              HDF5Tools
                ::read_BackendVector_from_HDF5<GFS_CG>( gfs_cg
                                                        , inputdata
                                                        , dir.vcheatM1_cg_orig_h5file[iSetup]
                                                        , "/vcheatM1_cg"
                                                        , vcheatM1_cg
                                                        );
              HDF5Tools
                ::read_BackendVector_from_HDF5<GFS_CG>( gfs_cg
                                                        , inputdata
                                                        , dir.vcheatArrival_Time_orig_h5file[iSetup]
                                                        , "/vcheatArrival_Time_orig"
                                                        , vcheatAT_orig
                                                        );
            }


            //set original measurements
            /*
              logger << "HEAT: take_measurements of heat M0 and M1 of original field to inputdata" << std::endl;
              take_measurements_on_dgf( gfs_cg,
              vcheatM0_cg,
              inputdata.setups[iSetup].heat_inversion_data.mpM0list.pointdata_vector,
              helper );

              take_measurements_on_dgf( gfs_cg,
              vcheatM1_cg,
              inputdata.setups[iSetup].heat_inversion_data.mpM1list.pointdata_vector,
              helper );
            */


            if(inputdata.problem_types.heat_mean_arrival_time_inversion){
              //Take measurements !!!
              logger << "HEAT: take_measurements_AT of original field" << std::endl;
              orig_measurements.take_measurements_AT( 4,
                                                      vcheatM0_cg,
                                                      vcheatM1_cg,
                                                      gfs_cg,
                                                      helper,
                                                      iSetup );
            }

            if( inputdata.plot_options.vtk_plot_heat ){
              VTKPlot::output2vtu( gfs_tp, vcheatM0_orig, dir.heatM0_orig_vtu[iSetup],
                                   "heatM0",
                                   inputdata.verbosity, true, std::max(0,pMAX-1) );

              VTKPlot::output2vtu( gfs_cg, vcheatM0_cg, dir.heatM0_orig_vtu[iSetup] + "_cg",
                                   "heatM0_cg",
                                   inputdata.verbosity, true, 0 );

              VTKPlot::output2vtu( gfs_cg, vcheatM1_cg, dir.heatM1_orig_vtu[iSetup] + "_cg",
                                   "heatM1_cg",
                                   inputdata.verbosity, true, 0 );

              VTKPlot::output2vtu( gfs_cg, vcheatAT_orig, dir.heatArrival_Time_orig_vtu[iSetup] + "_cg",
                                   "heat_arrival_time_orig",
                                   inputdata.verbosity, true, 0 );
            }

          } // END: if(inputdata.problem_types.heat_forward || inputdata.problem_types.heat_mean_arrival_time_inversion)

          //} //END: for all setups

        } // end of loop over iSetup

      } // if( inputdata.problem_types.synthetic &&( inputdata.problem_types.head_forward || inputdata.problem_types.head_inversion ||inputdata.problem_types.transport_forward || inputdata.problem_types.transport_inversion_m0 || inputdata.problem_types.transport_inversion_m1|| inputdata.problem_types.heat_forward||inputdata.problem_types.head_HT_inversion || inputdata.problem_types.heat_mean_arrival_time_inversion || inputdata.problem_types.moments_geoeletric_potential_forward || inputdata.problem_types.moments_geoeletric_potential_inversion))
      else /*ONLY the synthetic case */if(
                                          inputdata.problem_types.head_inversion ||
                                          inputdata.problem_types.transport_inversion_m0 ||
                                          inputdata.problem_types.transport_inversion_m1 ||
                                          inputdata.problem_types.heat_mean_arrival_time_inversion ||
                                          inputdata.problem_types.moments_geoeletric_potential_inversion ) {

          // here might be something coming!
          // in our implementation a Yfield is generated if syntetic yes or no. This is some overhead. The Yfield is only needed in the synthetic case. BUT: Then here is the generation of the eigenvalues needed and the writing of the zonation matrix. In the ERT_inversion code it was already changed!

          //phi 0 needs to be generated!
          if(inputdata.problem_types.moments_geoeletric_potential_inversion){
            Vector<REAL> logsigma0;

            HDF5Tools::h5g_pRead( gv_gw
                                  , logsigma0
                                  , dir.logsigma0_h5file
                                  , "/logsigma0"
                                  , local_offset
                                  , local_count
                                  , inputdata
                                  , 1 // P0 blocksize
                                  , FEMType::DG // P0
                                  , 0 // structure is on grid level 0
                                  );




            // export the data to the sigma0 Generator
            if(helper.size()>1)
              log_electricalConductivity.parallel_import_from_local_vector(logsigma0, local_offset, local_count );
            else
              log_electricalConductivity.import_from_vector( logsigma0 );

            Vector<REAL> logkappa;

            HDF5Tools::h5g_pRead( gv_gw
                                  , logkappa
                                  , dir.logkappa_h5file
                                  , "/logkappa"
                                  , local_offset
                                  , local_count
                                  , inputdata
                                  , 1 // P0 blocksize
                                  , FEMType::DG // P0
                                  , 0 // structure is on grid level 0
                                  );

            // export the data to the sigma0 Generator
            if(helper.size()>1)
              log_kappafield.parallel_import_from_local_vector(logkappa, local_offset, local_count );
            else
              log_kappafield.import_from_vector( logkappa );


            typedef GEPotentialForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GEP_FWD;
            for( UINT iSetup=0; iSetup<inputdata.setups.size(); iSetup++ ) {
              const SDT& setupdata = inputdata.setups[iSetup];


              GWP_FWD gwp_fwd( inputdata, setupdata, yfg_orig );
              GEP_FWD gep_fwd( inputdata, setupdata, log_electricalConductivity );

              //kappa is not needed here. I don't like the changes in the Forward simulator here. This works better in the old code. SEE trunk version!
              GEP_FWD log_kappawrapper_gw( inputdata, setupdata, log_kappafield ); // equation (66), (147): kappa


              typedef ForwardSimulator<GWP_FWD,
                                       GEP_FWD,
                                       GFS_GW,
                                       GFS_TP,
                                       GFS_CG,
                                       IDT,
                                       SDT,
                                       DIR> FWD_SIM;

              FWD_SIM fwd_sim( gwp_fwd,
                               gep_fwd,
                               log_kappawrapper_gw,
                               gfs_gw,
                               inputdata,
                               setupdata,
                               dir);
              for(UINT iconfig=0; iconfig<inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig; iconfig++){

                VCType_GW vcphi0_orig(gfs_gw,0.0);


                fwd_sim.GE_phi0_simulation( iconfig,
                                            vcphi0_orig );
                HDF5Tools::write_BackendVector_to_HDF5<GFS_GW>( gfs_gw
                                                                , inputdata
                                                                , dir.vcphi0_orig_h5file[iSetup][iconfig]
                                                                , "/vcphi0_orig"
                                                                , vcphi0_orig
                                                                );
                if( inputdata.plot_options.vtk_plot_el_potential_field ){
                  VTKPlot::output2vtu( gfs_gw,
                                       vcphi0_orig,
                                       dir.phi0_orig_vtu[iSetup][iconfig],
                                       "phi0_orig",
                                       inputdata.verbosity,
                                       true,
                                       0 );
                }
              }



            }

          }

        }

      // log the original measurements
      orig_measurements.write_to_logger("original measurements for inversion");


      if (inputdata.problem_types.CR){
        //start generation of conditional realizations
        // start the inversion kernel
        if( inputdata.CR_parameters.total
            && (
                inputdata.problem_types.lnK_inversion
                ||
                inputdata.problem_types.head_inversion
                ||
                inputdata.problem_types.transport_inversion_m0
                ||
                inputdata.problem_types.transport_inversion_m1
                ||
                inputdata.problem_types.heat_mean_arrival_time_inversion
                ||
                inputdata.problem_types.moments_geoeletric_potential_inversion
                )
            )
          {

            generate_CR<GRID,
                        PGV,
                        GFS_GW,
                        GFS_TP,
                        GFS_CG,
                        MEASLIST,
                        DIR,
                        IDT,
                        //SDT,
                        YFG
                        >
              ( grid,
                pRootGridView,
                gfs_gw,
                inputdata,
                yfg_orig,
                orig_measurements,
                //ACHTUNG ÄNDERUNG ZONES
                nCellsExt,
                dir,
                log_electricalConductivity,
                log_kappafield,
                helper
                );

          }else
          logger<<"Nothing to do for generating "<<inputdata.CR_parameters.total<<" conditional realizations! (Maybe no inversion data!)"<<std::endl;
      }else {
        // start the inversion kernel
        if(
           inputdata.problem_types.lnK_inversion
           ||
           inputdata.problem_types.head_inversion
           ||
           inputdata.problem_types.transport_inversion_m0
           ||
           inputdata.problem_types.transport_inversion_m1
           ||
           inputdata.problem_types.heat_mean_arrival_time_inversion
           ||
           inputdata.problem_types.moments_geoeletric_potential_inversion
           )
          {


            if( inputdata.problem_types.generate_measurement_data ){
              if( helper.rank()==0
                  &&
                  inputdata.verbosity>=VERBOSITY_INVERSION
                  ){
                std::cout << "This run was just to generate measurement data. Exit here." << std::endl;
              }
            }
            else {

              REAL L = inversionkernel<GRID,
                                       PGV,
                                       GFS_GW,
                                       GFS_TP,
                                       GFS_CG,
                                       MEASLIST,
                                       DIR,
                                       IDT,
                                       YFG
                                       >
                ( grid,
                  pRootGridView,
                  gfs_gw,
                  inputdata,
                  orig_measurements,
                  //ACHTUNG ÄNDERUNG ZONES
                  nCellsExt,
                  dir,
                  log_electricalConductivity,
                  log_kappafield,
                  helper
                  );

              if( helper.rank()==0
                  &&
                  inputdata.verbosity>=VERBOSITY_INVERSION )
                std::cout << "Return value of the inversion kernel L = " << L << std::endl;

              logger << "Return value of the inversion kernel L = " << L << std::endl;

            }
          }

      }

      return General::reportTotalTime( helper.rank() );

    }


  }
}
#endif // DRIVER_HH
