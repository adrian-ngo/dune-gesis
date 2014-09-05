/* 
 * File:   SensitivityClass.hh
 * Author: Ronnie L. Schwede and Adrian Ngo, 2010-2014
 */

#ifndef SENSITIVITYCLASS_HH
#define SENSITIVITYCLASS_HH


#include "lnK_sensitivities.hh"
#include "head_sensitivities.hh"
#include "m0_sensitivities.hh"
#include "m1_sensitivities.hh"
#include "heat_sensitivities.hh"
#include "geoelectrical_potential_sensitivities.hh"

#include "../common/MyMPIComm.hh"

#include "../common/io/IO_routines.hh"

namespace Dune {
  namespace GeoInversion {

    template<
      typename POOL_GRID
      , typename PGV
      , typename GFS_GW
      , typename GFS_TP
      , typename GFS_CG
      , typename MEASLIST
      , typename YFG
      , typename DIR
      , typename IDT
      //, typename SDT
      >
    class SensitivityClass{
      
      // Extract types from Gridfunction Spaces
      typedef typename GFS_GW::Traits::GridViewType GV_GW;
      typedef typename GFS_GW::Traits::FiniteElementMapType FEM_ELLIP;
      typedef typename GFS_GW::Traits::BackendType VBE_GW;
      typedef typename GFS_GW::Traits::ConstraintsType CON_GW;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;

      typedef typename GFS_TP::Traits::GridViewType GV_TP;
      typedef typename GFS_TP::Traits::FiniteElementMapType FEM_HYPER;
      typedef typename GFS_TP::Traits::BackendType VBE_TP;
      typedef typename GFS_TP::Traits::ConstraintsType CON_TP;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;

      typedef typename GFS_CG::Traits::FiniteElementMapType FEM_CG;
      typedef typename GFS_CG::Traits::BackendType VBE_CG;
      typedef typename GFS_CG::Traits::ConstraintsType CON_CG;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;

      //typedef typename POOL_GRID::LevelGridView POOL_GV_GW;

      const PGV pRootGridView;
      const IDT& inputdata;


      const std::vector< Vector<UINT> >& nCellsExt;
      DIR& dir;


      UINT nSetups;
      UINT nzones;

      MEASLIST& orig_measurements;
      UINT nMeas;

      DenseMatrix<REAL> Gyy_OLD;
      DenseMatrix<REAL> Gyy_NEW;

      const Dune::MPIHelper& helper;

      //data of the MPI pools
      UINT number_of_MPIPools;
      std::vector< std::vector<UINT> > pool_lookup;
      std::vector< std::vector<MyMPIComm> > CommunicatorPool;
      std::vector<UINT> final_pool_size;
     
      
      // Finite Element Maps:
      //Dune::shared_ptr<FEM_CG> pfem_cg;

      // POOL_GRID, POOL_GV, and GFS as pointer. They need to be known on all processes but with different values!
      std::vector<Dune::shared_ptr<POOL_GRID>> pool_grid;

      std::vector<GV_GW> pool_gv_0;
      std::vector<GV_GW> pool_gv_gw;
      std::vector<Dune::shared_ptr<CON_GW>> pool_con_gw;
      std::vector<Dune::shared_ptr<GFS_GW>> pool_gfs_gw;

      /*
      std::vector<Dune::shared_ptr<GV_TP>> pool_gv_tp;
      std::vector<Dune::shared_ptr<CON_TP>> pool_con_tp;
      std::vector<Dune::shared_ptr<GFS_TP>> pool_gfs_tp;
      std::vector<Dune::shared_ptr<CON_CG>> pool_con_cg;
      std::vector<Dune::shared_ptr<GFS_CG>> pool_gfs_cg;
      */

      std::vector<std::vector< Vector<REAL> > > pool_Lambdas;
      std::vector<std::vector< Vector<REAL> > > pool_X; // zonation matrix

      std::vector<YFG*> YfieldGenerator_old;

      //initialize the sigm0 and kappa field for the pools! only for geoelectrical potential measurements!
      std::vector<Dune::shared_ptr<YFG>> log_electricalConductivity_pool;
      std::vector<Dune::shared_ptr<YFG>> log_kappafield_pool;
      
      std::vector< Vector<REAL> > JX;
      std::vector< Vector<REAL> > JQJ;
      Vector<REAL> J_times_Y_old;

    public:
      // The Destructor:
      ~ SensitivityClass
      <POOL_GRID,
       PGV,
       GFS_GW,
       GFS_TP,
       GFS_CG,
       MEASLIST, 
       YFG,
       DIR,
       IDT
       //, SDT
       >(){
        // free the memory space
        for(UINT ii=0; ii<number_of_MPIPools;ii++){
          if( YfieldGenerator_old[ii] != NULL )
            delete( YfieldGenerator_old[ii] );
        }
        // Note: Dune::shared_ptr<...> get deleted automatically the class gets deleted!
      }



      // The Constructor:
      SensitivityClass
      <POOL_GRID,
       PGV,
       GFS_GW,
       GFS_TP,
       GFS_CG,
       MEASLIST, 
       YFG,
       DIR,
       IDT
       //, SDT
       >
      (
       POOL_GRID & theGrid,
       const PGV pRootGridView_,
       const IDT& inputdata_,
       const std::vector< Vector<UINT> >& nCellsExt_, // related to Lambdas
       DIR& dir_,
       MEASLIST& orig_measurements_,
       const Dune::MPIHelper& helper_
       )
      : pRootGridView(pRootGridView_),
        inputdata(inputdata_),
        nCellsExt(nCellsExt_),
        dir(dir_),
        orig_measurements(orig_measurements_),
        nMeas(orig_measurements_.nMeasurements()),
        Gyy_OLD(nMeas,nMeas,0.0),
        Gyy_NEW(nMeas,nMeas,0.0),
        helper(helper_) 
      {
        
        typedef typename GFS_GW::Traits::GridViewType GV_GW;
        enum{dim=GV_GW::dimension};

        // Total number of all cells required to resolve the parameter field
        // const UINT nAllCells = inputdata.domain_data.nCells.componentsproduct();
        // number of zones

        nzones = inputdata.yfield_properties.nz;
        //std::cout << "DEBUG: nzones = " << nzones << std::endl;

        nSetups=inputdata.setups.size();


        typedef Dune::FieldVector<CTYPE,dim> COORDINATES;
        //small stuff for making a grid
        COORDINATES extensions(1.0);
        for (UINT i = 0; i < dim; i++)
          extensions[i] = inputdata.domain_data.extensions[i];

#ifdef USE_YASP

        Dune::array<int,dim> nYaspGridCells( Dune::fill_array<int,dim>(1) );
        for (UINT i = 0; i < dim; i++)
          nYaspGridCells[i] = inputdata.domain_data.yasp_nCells[i];

        std::bitset<dim> periodic(false);

#endif        
        // set the number of implemented measurements --> needed for the MPIPool 
        /*
         * 0 = lnK
         * 1 = head
         * 2 = m0
         * 3 = m1
         * 4 = heat (AT)
         * 5 = geoelectrical potential (mean arrival time (AT))
         */
        UINT number_of_implemented_measurements=6;
        // this is the number of different MPIPools
        number_of_MPIPools=0;
  
        // look up Matrix: each row of the matrix belongs to a measurement setup.
        // each entry in each row corresponds to the measurement type. The value is the used MPI pool
        pool_lookup.resize(nSetups);
        
        //loop over all setups
        for(UINT iSetup=0; iSetup<nSetups; iSetup++){       
            
          pool_lookup[iSetup].resize(number_of_implemented_measurements,0);

            
          /*
           * SET THE MPI POOL
           */
          //logger<<"MPI Pool(s): "<<inputdata.parallel_fine_tuning.maxNComm<<" pools (max. possible number) out of "<<helper.size()<<" processes!"<<std::endl;
  
          // vector of vectors holding all the communicator information
          // first index is for different sizes of communicator pools, 
          //  e.g if you have 4 lnK and 36 head measurements and you want 12 pools then for the lnK measurements only a 4 pools (with all processors) are generated
          // second vector refers to the size of the communicator pool, in the example CommunicatorPool[0].size()==4 and CommunicatorPool[1].size()=12
          // MyMPIComm class holds information about the communicator -> see MyMPIComm.hh
          //CommunicatorPool
        
          //Communicator pool size for the measurement type 
          //  (if needed this vector will grow, if you only want 2 pools and you have for each type more than 2 measurements this vector will have one element)
          //std::vector<UINT> final_pool_size;

          // if more than one MPI Pool is requested
          if(inputdata.parallel_fine_tuning.maxNComm>1){

            //size of the global MPI communicator
            INT world_size=helper.size();

            //loop over all implemented measurements -> IMPORTANT! needes to be extended for each new measurment type
            // last measurement type = GE measurements has different configurations. here the configurations are set to be he measurements
            //  if the configurations have very uneaven distributed # of measurements this leads to unbalanced work load!!!
            for (UINT meas_type=0; meas_type<number_of_implemented_measurements; meas_type++){
                
              int n=0;
                    
              if(meas_type==5){
                if(inputdata.problem_types.moments_geoeletric_potential_inversion)
                  n=orig_measurements.nGE_config(iSetup);
              }else{ 
                n=orig_measurements.nMeasPerSetupPerType(iSetup,meas_type);
              }
              if(n){
                // set the final pool size. It is the max. of requested pools and numbner of measutrements
                  UINT tmp_pool_size;
                if (n<inputdata.parallel_fine_tuning.maxNComm)
                  tmp_pool_size=n;
                else
                  tmp_pool_size=General::GCD(n,inputdata.parallel_fine_tuning.maxNComm);
            
                bool newpool=true;
                // check if a new communicator pool is needed.
                //   it is needed only if not a communicator pool of this size will be generated! 
                for (UINT ii=0; ii<number_of_MPIPools;ii++){
                  // adequate pool found?
                  if(tmp_pool_size==final_pool_size[ii]){
                    // an adequate pool is found:
                    // set the pool_lookup to the right value and set that no new pool is needed
                    pool_lookup[iSetup][meas_type]=ii;
                    newpool=false;
                    break;
                  }
                }
                if(newpool){
                  // set the values for a new communicator pool
                  pool_lookup[iSetup][meas_type]=number_of_MPIPools;
                  number_of_MPIPools+=1;
                  final_pool_size.push_back(tmp_pool_size);
                }
              }

            }


            if( helper.rank() == 0 && inputdata.verbosity >= VERBOSITY_INVERSION ){
              std::cout << "SensitivityClass: Parallel-in-parallel active!" << std::endl;
              std::cout << "SensitivityClass: world_size = " 
                        << world_size
                        << std::endl;
              std::cout << "SensitivityClass: maxNComm = " 
                        << inputdata.parallel_fine_tuning.maxNComm
                        << std::endl;
              std::cout << "SensitivityClass: We will work with " 
                        << number_of_MPIPools << " MPI-Pool(s):"
                        << std::endl;
              for (UINT ii=0; ii<number_of_MPIPools;ii++){
                std::cout << "SensitivityClass: final_pool_size["<<ii<<"] = " 
                          << final_pool_size[ii] << " Communicators"
                          << std::endl;
              }
            }


            //number_of_MPIPools has the right number of needed communicator POOLS!
            CommunicatorPool.resize(number_of_MPIPools);
        
            /*
             * generate the needed communicator pool(s)
             */
            for(UINT iiPool=0; iiPool<number_of_MPIPools; iiPool++){
        
              //needed size of the vector(for all communicators) for this communicator pool
              CommunicatorPool[iiPool].resize(final_pool_size[iiPool]);
        
              // vector of vector!! holds in each row "i" the process ranks (global numbers) which should belong to the group of "i" 
              std::vector<int*> processor_ranks(final_pool_size[iiPool]);
        
              // vector containing the  sizes of this MPI pools
              std::vector<int> PoolSizes(final_pool_size[iiPool]);
        
              // calculate the pool size
              for(UINT ii=0;ii<final_pool_size[iiPool];ii++){
        
                // this is the poolsize
                PoolSizes[ii]=(int)(world_size/final_pool_size[iiPool]);
        
                //distribute to the last pools the left over processors
                //if(ii>=final_pool_size[iPool]-(world_size%final_pool_size[iPool]))  // larger last groups at the end
                if( ii<(world_size%final_pool_size[iiPool]) ) // larger first groups
                  PoolSizes[ii]++;


                if( helper.rank() == 0 && inputdata.verbosity >= VERBOSITY_INVERSION ){
                  std::cout << "SensitivityClass: Pool " << iiPool
                            << ", Communicator " << ii 
                            << " uses " << PoolSizes[ii] << " processor(s)."
                            << std::endl;
                }

                //allocate memory for the prozessor ranks
                processor_ranks[ii] = new int[PoolSizes[ii]];
              }
        
              /*
               * get the ranks (global ranks) for each MPI pool
               */
              //needed counter
              int iPool=0;
              int counter=0;
        
              //loop over all available ranks(global ranks) 
              for(int ii=0;ii<world_size;ii++){
        
                //normal case the counter is smaller than the considered MPI pool size
                if(counter<PoolSizes[iPool]){
        
                  //add the processor with rank "ii" to the list
                  processor_ranks[iPool][counter]=ii;
        
                  //if this processor is the one we are working one, then ...
                  if(helper.rank()==ii) 
                    // ... set in the MyMPIComm the flag, that this processor belongs to this group.
                    // This is very important, because it is not done automatically by any function provided by MyMPIComm class!
                    // Needed when ever collactive calls are done on the communicator, because this code should only be called on the processors of the group!!!!
                    CommunicatorPool[iiPool][iPool].set_I_am_in();
        
                  //finally add 1 to the counter
                  counter++;
        
                }else{ // if counter == PoolSizes[iPool], then the pool is full and we should start with the next one!
            
                  // add 1 to iPool
                  iPool++;
        
                  //set processor with rank "ii" at first position (will become P0 in this communicator)
                  processor_ranks[iPool][0]=ii;
        
                  //if this processor is the one we are working one, then ...
                  if(helper.rank()==ii)
                    // ... set in the MyMPIComm the flag, that this processor belongs to this group.
                    CommunicatorPool[iiPool][iPool].set_I_am_in();
        
                  // set the counter to 1
                  counter=1;
                }
              } //END: for(int ii=0;ii<world_size;ii++){
        
              /*
               * generating the communicators (MPI routines!!)
               */
              MPI_Group world_group;
              MPI_Comm_group(helper.getCommunicator(),&world_group);

              // loop over the needed MPI pools
              for(UINT iPool=0;iPool<final_pool_size[iiPool];iPool++){ 
        
                // generate a new MPI group out of the list of processors and store it in the right location in CommunicatorPool
                MPI_Group_incl(world_group,PoolSizes[iPool],processor_ranks[iPool],CommunicatorPool[iiPool][iPool].get_group_adress());
        
                // generate a new MPI communicator out of the new generated group and store it in the right location in CommunicatorPool
                MPI_Comm_create(helper.getCommunicator(),CommunicatorPool[iiPool][iPool].get_group(),CommunicatorPool[iiPool][iPool].get_comm_adress());
        
                //update the MyMPIComm:
                // important so that the rank and the size (local in the stored group/communicator) are set! (the flag if a process is in a group or not needs to be set before by hand (see above) )
                CommunicatorPool[iiPool][iPool].update();
              }
        
              //delete the allocated processor rank list!
              for(UINT ii=0;ii<final_pool_size[iiPool];ii++){    
                delete [] (processor_ranks[ii]);
              }
                    
        
            } // END: for(UINT iiPool=0; iiPool<number_of_MPIPools; iiPool++)    

          }else{  // only ONE MPI pool per measurement is needed. This pool/communicator is already stored in the "MPIHelper helper"
            number_of_MPIPools=1;
            CommunicatorPool.resize(1);
            final_pool_size.resize(1,1);
            MPI_Group world_group;
            MPI_Comm_group(helper.getCommunicator(),&world_group);
            CommunicatorPool[0].resize(1); //final_pool_size[ii]==1  !!!!
            //CommunicatorPool[0][0].set(helper.getCommunicator(),world_group, true);
            CommunicatorPool[0][0].set(MPI_COMM_WORLD,world_group, true);
          }  
          /*
           * SET THE MPI POOL: DONE!!!
           */
            
        }// END         for(UINT iSetup=0; iSetup<nSetups; iSetup++){       

            
#ifdef DEBUG_LOG
        //plot pool information!!
        logger<<std::endl<<"HERE COMES THE POOL INFO ( of "<<number_of_MPIPools<<" Pools) : "<<std::endl<<std::endl;
        for (UINT iPool=0; iPool<number_of_MPIPools; iPool++){
          logger<<std::endl<<"Communicator Pool #"<<iPool<<" has "<<final_pool_size[iPool]<<" Communicators : "<<std::endl;
          for(UINT ii=0; ii<final_pool_size[iPool];ii++){
            logger<<"Communicator # "<<ii<<std::endl;
            CommunicatorPool[iPool][ii].print_info();
          }
        }
        logger<<std::endl<<"POOL INFO : DONE!!! "<<std::endl<<std::endl; 
        logger<<"pool_lookup :"<<std::endl;
        for(UINT iSetup=0; iSetup<nSetups; iSetup++){ 
          logger<<"For Setup #"<<iSetup+1<<" :"<<std::endl;
          logger<< "pool_lookup : lnK pool#"<<pool_lookup[iSetup][0]<<", head pool#"<<pool_lookup[iSetup][1]<<", m0 pool#"<<pool_lookup[iSetup][2]<<", m1 pool#"<<pool_lookup[iSetup][3]<<", heat pool#"<<pool_lookup[iSetup][4]<<", GE pool#"<<pool_lookup[iSetup][5];
          logger<<std::endl; 
        }
        logger<<std::endl;
#endif



#ifdef USE_CUBE
        const Dune::GeometryType::BasicType bt = Dune::GeometryType::cube;
#else
        const Dune::GeometryType::BasicType bt = Dune::GeometryType::simplex;
#endif

        //pfem_cg    = Dune::make_shared<FEM_CG>(gv_tp);

        pool_Lambdas.resize(number_of_MPIPools);
        pool_X.resize(number_of_MPIPools); // zonation matrix
        Vector<UINT> local_count,local_offset;
        std::vector<ptrdiff_t> alloc_local(nzones),local_n0(nzones), local_0_start(nzones);
        std::vector<UINT> local_FFT_size(nzones),local_1D_size(nzones),local_1D_start(nzones);
        
        //loop over the needed communicator pool(s)
        for (UINT iPool=0; iPool<number_of_MPIPools; iPool++){
          // loop over all MPI pools to set the values of pool_grid, pool_gv, and,pool_gfs
          for(UINT ii=0;ii<final_pool_size[iPool];ii++){
            // do only something if the processor belongs to the considered group,
            if(CommunicatorPool[iPool][ii].I_am_in()){  

#ifdef USE_YASP
#if defined USE_NOVLP_MODE
              int overlap = 0;
#else 
              // overlapping-case: get the needed overlap
              int overlap = inputdata.domain_data.yasp_overlap;
#endif
#endif
              // define the grid for this MPI pool
#ifdef DIMENSION3

#ifdef USE_YASP
              Vector<UINT> yasp_partitions(inputdata.domain_data.nPartitions);

              // How many processes are there in the current pool?
              const UINT np = CommunicatorPool[iPool][ii].get_size();
          
              if( yasp_partitions.componentsproduct() != np ){
                // calculate the optimal partitioning
                UINT n1,n2, larger, smaller;
                General::factorize_two_optimally(np,n1,n2);

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
                
                logger << "Warning: Partitioning of YASP for pool " << iPool << " might have changed!" << std::endl;
                if( helper.rank()==0 && inputdata.verbosity>0 )
                  std::cout << "Warning in SensitivityClass: Partitioning of YASP for pool " << iPool << " might have changed!" << std::endl;
                
              }

              logger << "Partitioning of YASP for pool " << iPool << ": " << yasp_partitions << std::endl;
              if( helper.rank()==0 && inputdata.verbosity>0 )
                std::cout << "SensitivityClass: Partitioning of YASP for pool " << iPool << ": " << yasp_partitions << std::endl;
              
              typedef YaspPartition<Vector<UINT>,dim> YP;
              //Dune::shared_ptr<YP> yp = Dune::make_shared<YP>( yasp_partitions );
              YP* yp = new YP( yasp_partitions );
          
              pool_grid.push_back( Dune::make_shared<POOL_GRID>( CommunicatorPool[iPool][ii].get_comm(), 
                                                                 extensions, 
                                                                 nYaspGridCells, 
                                                                 periodic, 
                                                                 overlap, 
                                                                 yp
                                                                 ) );
              if( yp != NULL )
                delete yp;      // TODO: Check if this really works with many pools!

#endif // USE_YASP

#else // 2D:

#ifdef USE_YASP
      pool_grid.push_back( Dune::make_shared<POOL_GRID>( CommunicatorPool[iPool][ii].get_comm(), 
                                                         extensions, 
                                                         nYaspGridCells, 
                                                         periodic, 
                                                         overlap 
                                                         ) );
#endif

#endif // DIMENSION



#ifdef USE_ALUGRID
      // make grid
      Dune::FieldVector<REAL,dim> LowerLeft(0);
      Dune::FieldVector<REAL,dim> UpperRight;
      UpperRight[0] = inputdata.domain_data.extensions[0];
      UpperRight[1] = inputdata.domain_data.extensions[1];
      Dune::array<unsigned int,dim> elements;
      elements[0] = inputdata.domain_data.nCells[0];
      elements[1] = inputdata.domain_data.nCells[1];

#ifdef DIMENSION3
      UpperRight[2] = inputdata.domain_data.extensions[2];
      elements[2] = inputdata.domain_data.nCells[2];
#endif

      Dune::Timer watch;
      Dune::StructuredGridFactory<POOL_GRID> structuredGridFactory;
      Dune::shared_ptr<POOL_GRID> gridptr =
        structuredGridFactory.createCubeGrid(LowerLeft, UpperRight, elements);
      REAL elapsed_time = watch.elapsed();
      std::cout << "=== Dune::StructuredGridFactory building ALUGRID took "
                << elapsed_time << " sec." << std::endl;

      watch.reset();
      if( gridptr->loadBalance() ){
        elapsed_time = watch.elapsed();
        std::cout << "=== Initial load balancing for POOL ALUGRID took "
                  << elapsed_time << " sec." << std::endl;
      }

      pool_grid.push_back( gridptr );
#endif // USE_ALUGRID


#ifdef USE_YASP
      int baselevel = inputdata.domain_data.yasp_baselevel;
      // It is very important to keep the number of overlaps constant during refinement!
      // Otherwise, the solution of the transport problem will go wrong.
      bool keepPhysicalOverlap = false;
      pool_grid[iPool]->refineOptions( keepPhysicalOverlap );

#endif

#ifdef USE_UG
      int baselevel = inputdata.domain_data.ug_baselevel;
#endif

#ifdef USE_ALUGRID
      int baselevel = inputdata.domain_data.alugrid_baselevel;
#endif

      if( baselevel>0 )
        pool_grid[iPool]->globalRefine(baselevel);

#ifdef USE_YASP
              int maxlevel=inputdata.domain_data.yasp_maxlevel;
#endif
#ifdef USE_ALUGRID
              int maxlevel = inputdata.domain_data.alugrid_baselevel + inputdata.domain_data.alugrid_maxsteps;
#endif
#ifdef USE_UG
              int maxlevel = inputdata.domain_data.ug_baselevel + inputdata.domain_data.ug_maxsteps;
#endif

              if(maxlevel>0)
                pool_grid[iPool]->globalRefine(maxlevel);
    
              if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
                std::cout << "SensitivityClass: Pool grid baselevel = " << baselevel << std::endl;
                std::cout << "SensitivityClass: Pool grid leaflevel = " << maxlevel << std::endl;
              }


              // gfs on all pools (not needed for the lnk-pool, but maybe this pool is also used for other measurement types)

              // define the grid view for this MPI pool -> needed on all measurements
              GV_GW gv_0 = pool_grid[iPool]->levelGridView(0);
              GV_GW pgv_gw = pool_grid[iPool]->levelGridView(baselevel);
              Dune::shared_ptr<CON_GW> pcon_gw = Dune::make_shared<CON_GW>();

              Dune::shared_ptr<FEM_ELLIP> pfem_ellip;
#ifdef USE_FEM
              pfem_ellip = Dune::make_shared<FEM_ELLIP>(gv_0);
#else
              Dune::GeometryType gt = Dune::GeometryType(bt,dim);
              pfem_ellip = Dune::make_shared<FEM_ELLIP>(gt);
#endif
              Dune::shared_ptr<GFS_GW> pgfs_gw = Dune::make_shared<GFS_GW>( pgv_gw, 
                                                                            pfem_ellip, 
                                                                            pcon_gw );

              /*
              Dune::shared_ptr<GV_TP> pgv_tp = Dune::make_shared<GV_TP>( pool_grid[iPool]->leafView());
              Dune::shared_ptr<CON_TP> pcon_tp = Dune::make_shared<CON_TP>();
              Dune::shared_ptr<GFS_TP> pgfs_tp = Dune::make_shared<GFS_TP>( *pgv_tp, 
                                                                            *pfem_hyper, 
                                                                            *pcon_tp );

              Dune::shared_ptr<CON_CG> pcon_cg = Dune::make_shared<CON_CG>();
              Dune::shared_ptr<GFS_CG> pgfs_cg = Dune::make_shared<GFS_CG>( *pgv_tp, 
                                                                            *pfem_cg, 
                                                                            *pcon_cg );
                                                                            */

#ifdef USE_NOVLP_MODE
              pcon_gw->compute_ghosts( *pgfs_gw );
              //pcon_tp->compute_ghosts( *pgfs_tp );
#endif     

              pool_gv_0.push_back( gv_0 );
              pool_gv_gw.push_back( pgv_gw );
              pool_con_gw.push_back( pcon_gw );

              pool_gfs_gw.push_back( pgfs_gw );

              /*
                pool_gv_tp.push_back( pgv_tp );
              pool_con_tp.push_back( pcon_tp );
              pool_gfs_tp.push_back( pgfs_tp );              
              pool_con_cg.push_back( pcon_cg );
              pool_gfs_cg.push_back( pgfs_cg );
              */           

              if(inputdata.problem_types.moments_geoeletric_potential_inversion){
                for(UINT iSetup=0; iSetup<nSetups; iSetup++){
                  if(iPool==pool_lookup[iSetup][5]){
                    log_electricalConductivity_pool.push_back( Dune::make_shared<YFG>( inputdata,dir,CommunicatorPool[iPool][ii].get_comm() ) );
                    log_kappafield_pool.push_back( Dune::make_shared<YFG>( inputdata,dir,CommunicatorPool[iPool][ii].get_comm() ) );
                    Vector<REAL> logsigma0;
                    HDF5Tools::
                      read_parallel_from_HDF5( pool_gv_gw[iPool]
                                               , inputdata
                                               , logsigma0
                                               , "/logsigma0"
                                               , local_count
                                               , local_offset
                                               , dir.logsigma0_h5file
                                               );
                    if(CommunicatorPool[iPool][ii].get_size()>1)
                      log_electricalConductivity_pool[iPool]->parallel_import_from_local_vector(logsigma0, local_count, local_offset );
                    else
                      log_electricalConductivity_pool[iPool]->import_from_vector( logsigma0 );
                
                    Vector<REAL> logkappa;
                    HDF5Tools::
                      read_parallel_from_HDF5(
                                              pool_gv_gw[iPool]
                                              , inputdata
                                              , logkappa
                                              , "/logkappa"
                                              , local_count
                                              , local_offset
                                              , dir.logkappa_h5file
                                              );
                
        
                    if(CommunicatorPool[iPool][ii].get_size()>1)
                      log_kappafield_pool[iPool]->parallel_import_from_local_vector(logkappa, local_count, local_offset );
                    else
                      log_kappafield_pool[iPool]->import_from_vector( logkappa );
                  }
                }
              }
          
              //get the parameters how the FFT distributes the data
              for(UINT jj=0; jj<nzones; jj++){
#ifdef DIMENSION3
                alloc_local[jj] = fftw_mpi_local_size_3d( nCellsExt[jj][2],
                                                          nCellsExt[jj][1],
                                                          nCellsExt[jj][0], 
                                                          CommunicatorPool[iPool][ii].get_comm(),
                                                          &(local_n0[jj]), 
                                                          &(local_0_start[jj]) );
#else
                alloc_local[jj] = fftw_mpi_local_size_2d( nCellsExt[jj][1],
                                                          nCellsExt[jj][0], 
                                                          CommunicatorPool[iPool][ii].get_comm(),
                                                          &(local_n0[jj]), 
                                                          &(local_0_start[jj]) );
#endif
    
              }

            } // END: if(CommunicatorPool[iPool][ii].I_am_in()){

          }// END:  for(UINT ii=0;ii<final_pool_size[iPool];ii++)
            
          // save this value into UINT variables!
          for(UINT ii=0; ii<nzones; ii++){
            local_FFT_size[ii]=(UINT)alloc_local[ii];
            local_1D_size[ii]=(UINT)local_n0[ii];
            local_1D_start[ii]=(UINT)local_0_start[ii];
          }
    
          for(UINT jj=0; jj<nzones; jj++){
            //setting parameter for the HDF5 read of the eigenvalues "Lambdas"
            local_count.resize(0);
            local_offset.resize(0);
            local_count.resize(dim,0);
            local_offset.resize(dim,0);
   
#ifdef DIMENSION3
            local_count[0]=nCellsExt[jj][0];
            local_count[1]=nCellsExt[jj][1];
            local_count[2]=local_1D_size[jj];
            local_offset[2]=local_1D_start[jj];
#else
            local_count[0]=nCellsExt[jj][0];
            local_count[1]=local_1D_size[jj];
            local_offset[1]=local_1D_start[jj];
#endif
            pool_Lambdas[iPool].resize(nzones);
            // read the needed data from the HDF5 file -> and again for the different MPI Pools
            for(UINT ii=0; ii<final_pool_size[iPool];ii++){
              if(CommunicatorPool[iPool][ii].I_am_in()){
                HDF5Tools::
                  read_parallel_from_HDF5_without_DUNE( inputdata, 
                                                        pool_Lambdas[iPool][jj],
                                                        local_count,
                                                        local_offset,
                                                        CommunicatorPool[iPool][ii].get_comm(),
                                                        "/FFT_R_YY", 
                                                        dir.EV_h5file[jj] );
     
                // read zones on pool leader
                if(CommunicatorPool[iPool][ii].get_rank()==0){
                  pool_X[iPool].resize(nzones);
                  for(UINT jj=0; jj<nzones; jj++)
                    HDF5Tools::
                      read_sequential_from_HDF5_without_DUNE(pool_X[iPool][jj],"/X",dir.zonation_matrix[jj]);
                }
              }
            }
        
          }       
        } // END: for(UINT iPool=0; iPool<number_of_MPIPools; iPool++)
        
        // be sure that all processors have all information about the groups
        if(helper.size()>1)
          MPI_Barrier(helper.getCommunicator()); // Barrier on MPI_COMM_WORLD
    
        YfieldGenerator_old.resize(number_of_MPIPools);
        
      }// END constructor!
    










      /*
       * set the data from the previous iteration step!
       */
      void set_OLD( int it_counter ){

        logger << "void set_OLD()" << std::endl;

        Vector<REAL> local_Y_old;
        //Vector<REAL> local_Y_old_Well;
        Vector<UINT> local_count,local_offset;
        /*
         * load the YfieldGenerator_old for all communicator pool(s) (also loaded for the lnK pool (little overhead))
         * idea for a speed up: load here also the old head, m0, m1 and so on!
         */
        
        for(UINT iPool=0; iPool<number_of_MPIPools;iPool++){
          // the YfieldGenerator_old is only needed for head (and maybe others in future) measurments
          for(UINT ii=0;ii<final_pool_size[iPool];ii++){
            if(CommunicatorPool[iPool][ii].I_am_in()){

              delete(YfieldGenerator_old[iPool]);
              
              YfieldGenerator_old[iPool] = new YFG( inputdata,
                                                    dir, 
                                                    CommunicatorPool[iPool][ii].get_comm() );
              
              // read the data of Y_old ->parallel on each MPI pool
              //logger << "read the old Y-field " << dir.Y_old_h5file << std::endl;
              HDF5Tools::
                read_parallel_from_HDF5(
                                        pool_gv_gw[iPool]
                                        , inputdata
                                        , local_Y_old
                                        , "/Y_old"
                                        , local_count
                                        , local_offset
                                        , dir.Y_old_h5file
                                        );
              /*
              HDF5Tools::
                read_parallel_from_HDF5(
                *(pool_gv_gw[iPool])
                                        , inputdata
                                        , local_Y_old_Well
                                        , "/Y_old_Well"
                                        , local_count
                                        , local_offset
                                        , dir.Y_old_Well_h5file
                                        );
              */


              if(CommunicatorPool[iPool][ii].get_size()>1)
                YfieldGenerator_old[iPool]->parallel_import_from_local_vector( local_Y_old,
                                                                               //local_Y_old_Well,
                                                                               local_count,
                                                                               local_offset );
              else {
                YfieldGenerator_old[iPool]->import_from_vector( local_Y_old);
                //inputdata.loglistOfAllWellCenters();
                //YfieldGenerator_old[iPool]->setWellConductivities();
              }

              logger << "SensitivityClass: " << std::endl;
              inputdata.loglistOfAllWellCenters();
              YfieldGenerator_old[iPool]->setWellConductivities( pool_gv_gw[iPool] );
              
#ifdef VTK_PLOT_YTRY

#ifdef USE_YASP
              int baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_UG
      int baselevel = inputdata.domain_data.ug_baselevel;
#endif

#ifdef USE_ALUGRID
      int baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
              std::stringstream file_YoldWell;
              file_YoldWell << dir.vtudir 
                            << "/Y_old_" << it_counter
                            << "_Well_pool" << iPool;
              YfieldGenerator_old[iPool]->plot2vtu( pool_gv_gw[iPool],
                                                    file_YoldWell.str(),
                                                    "Y_try",
                                                    baselevel );
#endif

            } // END: if(CommunicatorPool[iPool][ii].I_am_in())
          } // END: for(UINT ii=0;ii<final_pool_size[iPool];ii++){

        }// END: for(UINT iPool=0; iPool<number_of_MPIPools;iPool++)   
        
      }// END set_OLD
    





      /*
       * calculate the sensitivities
       */
      void calculate( UINT it_counter, 
                      Vector<REAL>& Y_old,
                      const UINT& nAllCells, 
                      MEASLIST &measurements, 
                      UINT & SV_increased, 
                      const INT& CR=-1 ) {

        logger << "void calculate(...)" << std::endl;


#ifdef USE_YASP
        int baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_UG
      int baselevel = inputdata.domain_data.ug_baselevel;
#endif

#ifdef USE_ALUGRID
      int baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
        Vector<REAL> Y_u(0);
        /*
         * important update Y_old for all process group leaders exept P0(global).
         */
        // update the new Y_old
        //read Y_old for all processor pool leaders(LOCAL P0)!
        Vector<UINT> local_count,local_offset;
        //P0 (global) has the data already!
        if(helper.rank()!=0){
          for (UINT iPool=0; iPool<number_of_MPIPools; iPool++){
            for(UINT ii=0;ii<final_pool_size[iPool];ii++){
              if(CommunicatorPool[iPool][ii].get_rank()==0){
                //logger<<std::endl<<"!!! UPDATE Y_OLD!!!"<<std::endl<<std::endl;
                if(Y_old.size()==0){
                  HDF5Tools::
                    read_sequential_from_HDF5(  Y_old
                                                , "/Y_old"
                                                , local_count
                                                , local_offset
                                                , dir.Y_old_h5file
                                                , inputdata
                                                );
                }
                if(Y_u.size()==0 && CR>-1){
                  HDF5Tools::
                    read_sequential_from_HDF5(Y_u
                                              , "/YField"
                                              , local_count
                                              , local_offset
                                              , dir.unconditionalField_h5file
                                              , inputdata
                                              );
                }
              } // END: if(CommunicatorPool[iPool][ii].get_rank()==0)
            } // END: for (INT ii=0;ii<final_pool_size[iPool];ii++)
          } // END: for (UINT iPool=0; iPool<number_of_MPIPools; iPool++)
        }else{
          if(CR>-1){
            HDF5Tools::
              read_sequential_from_HDF5( Y_u
                                         , "/YField"
                                         , local_count
                                         , local_offset
                                         , dir.unconditionalField_h5file
                                         , inputdata
                                         );
          }
        }// END: if(helper.rank()!=0)
        for(UINT ii=0; ii<Y_u.size(); ii++){
          Y_old[ii]-=Y_u[ii];
        }
        
        JX.resize(0);
        J_times_Y_old.resize(0);


        std::vector< Dune::shared_ptr<GV_TP> > pool_gv_tp(number_of_MPIPools);
        std::vector< Dune::shared_ptr<CON_TP> > pool_con_tp(number_of_MPIPools);
        std::vector< Dune::shared_ptr<GFS_TP> > pool_gfs_tp(number_of_MPIPools);
        std::vector< Dune::shared_ptr<CON_CG> > pool_con_cg(number_of_MPIPools);
        std::vector< Dune::shared_ptr<GFS_CG> > pool_gfs_cg(number_of_MPIPools);

        /*
        //loop over the needed communicator pool(s)
        for (UINT iPool=0; iPool<number_of_MPIPools; iPool++){
          // loop over all MPI pools to set the values of pool_grid, pool_gv, and,pool_gfs
          for(UINT ii=0;ii<final_pool_size[iPool];ii++){
            // do only something if the processor belongs to the considered group,
            if(CommunicatorPool[iPool][ii].I_am_in()){  
              
              Dune::shared_ptr<GV_TP> pgv_tp = Dune::make_shared<GV_TP>( pool_grid[iPool]->leafView());
              Dune::shared_ptr<CON_TP> pcon_tp = Dune::make_shared<CON_TP>();
              Dune::shared_ptr<GFS_TP> pgfs_tp = Dune::make_shared<GFS_TP>( *pgv_tp, 
                                                                            *pfem_hyper, 
                                                                            *pcon_tp );
              
              Dune::shared_ptr<CON_CG> pcon_cg = Dune::make_shared<CON_CG>();
              Dune::shared_ptr<GFS_CG> pgfs_cg = Dune::make_shared<GFS_CG>( *pgv_tp, 
                                                                            *pfem_cg, 
                                                                            *pcon_cg );


              
            }
          }
        }

        */

        //std::cout << "DEBUG: nzones = " << nzones << std::endl;
        
        JX.resize(nMeas);
        for(UINT ii=0; ii<nMeas; ii++)
          JX[ii].resize(nzones,0.0);
        J_times_Y_old.resize(nMeas);
        
        enum{ dim=GV_GW::dimension };
        Dune::Timer watch;
        REAL time=0.0,max_time=0.0;

        std::vector<VCType_GW*> head_old(number_of_MPIPools);
        std::vector<VCType_CG*> soluteM0_old(number_of_MPIPools);
        std::vector<VCType_CG*> soluteM1_old(number_of_MPIPools);

        
        typedef typename IDT::SDT SDT;
        //loop over all setups
        for(UINT iSetup=0; iSetup<nSetups; iSetup++){   

          const SDT& setupdata = inputdata.setups[iSetup];
      
          // preloaded multiple needed data (suffix _old):
          // head_old(needed everywhere exept lnK), solute M0(needed in M0,M1, and GP), solute M1(needed in M1 and GP)
                
          //load needed old data
          //for(UINT iPool=0; iPool<number_of_MPIPools; iPool++){
          //loop over measurement types. Start from 1 because lnK does not need the old data  
          for(UINT t_meas=1; t_meas<pool_lookup[iSetup].size(); t_meas++){
            //if inversion of the type is needed
            if(orig_measurements.nMeasPerSetupPerType(iSetup,t_meas)){      
              /*
               * head
               */
              if(head_old[pool_lookup[iSetup][t_meas]]==NULL){                            
                // load the old head data
                head_old[pool_lookup[iSetup][t_meas]] = new VCType_GW( *(pool_gfs_gw[pool_lookup[iSetup][t_meas]]) , 0.0 );
                
                HDF5Tools::
                  read_BackendVector_from_HDF5(
                                               *(pool_gfs_gw[pool_lookup[iSetup][t_meas]])
                                               , inputdata
                                               , dir.vchead_old_h5file[iSetup] // important: read from the right file. To get the head for the solute transport setup!
                                               , "/vchead_old"
                                               , *(head_old[pool_lookup[iSetup][t_meas]])
                                               , PRESERVE_STRUCTURE
                                               , baselevel
                                                );

#ifdef DEBUG_PLOT_H_OLD
                std::stringstream vtu_h_old;
                vtu_h_old << dir.vtudir
                          << "/h_old"
                          << "_s" << iSetup
                          << "_i" << it_counter;
                
                VTKPlot::output2vtu( *(pool_gfs_gw[pool_lookup[iSetup][t_meas]]),
                                     *(head_old[pool_lookup[iSetup][t_meas]]),
                                     vtu_h_old.str(),
                                     "h_old",
                                     inputdata.verbosity,
                                     true,
                                     0
                                     );
#endif //DEBUG_PLOT_H_OLD



              }

              // The reordered gridview is needed only for meas.type=2,3,4,5, alles ausser 0 und 1
              if( t_meas>1 ){
#ifdef NEW_GV

#ifdef USE_DGF_PressureField
                logger << "Using DGF_PressureField ... " << std::endl;
                typedef DGF_PressureField<GFS_GW,VCType_GW> RT0_PF;
                Dune::shared_ptr<RT0_PF> rt0_pressurefield 
                  = Dune::make_shared<RT0_PF>( *(pool_gfs_gw[pool_lookup[iSetup][t_meas]]), *(head_old[pool_lookup[iSetup][t_meas]]), baselevel );
#else

                logger << "Using RT0_PressureField ... " << std::endl;
                typedef GroundwaterForwardProblem<GV_GW,REAL,IDT,SDT,YFG> GWP_FWD_OLD;
                typedef RT0_PressureField<GV_GW,GFS_GW,GWP_FWD_OLD> RT0_PF;
                GWP_FWD_OLD gwp_fwd_old( inputdata, 
                                         setupdata,
                                         *(YfieldGenerator_old[ pool_lookup[iSetup][t_meas] ]) );
                Dune::shared_ptr<RT0_PF> rt0_pressurefield
                  = Dune::make_shared<RT0_PF>( pool_gv_gw[pool_lookup[iSetup][t_meas]], gwp_fwd_old, baselevel );
                rt0_pressurefield->assemble( *(pool_gfs_gw[pool_lookup[iSetup][t_meas]]), *(head_old[pool_lookup[iSetup][t_meas]]) );

#endif
                
                PressureLikeOrdering comparefunctor;
                //ReversePressureLikeOrdering comparefunctor;
                
                if( CommunicatorPool[pool_lookup[iSetup][t_meas]][0].get_rank()==0 
                    && inputdata.verbosity >= VERBOSITY_EQ_SUMMARY ){
                  std::cout << "SensitivityClass: Reorder gridview for MPI-pool " 
                            << pool_lookup[iSetup][t_meas] 
                            << " , setup " << iSetup
                            << " , measurement type " << t_meas
                            << std::endl;
                }
                
                Dune::shared_ptr<GV_TP> pgv_tp 
                  = Dune::make_shared<GV_TP>( pool_grid[ pool_lookup[iSetup][t_meas] ]->leafGridView(),
                                              rt0_pressurefield,
                                              comparefunctor );

#else
                // default gridview:
                Dune::shared_ptr<GV_TP> pgv_tp = Dune::make_shared<GV_TP>( pool_grid[ pool_lookup[iSetup][t_meas] ]->leafGridView() );
#endif
                

                verifyCalibratedWellsOnRefinedGrid(*pgv_tp,inputdata);

#ifdef DEBUG_PLOT
                std::stringstream vtu_elementorder;
                vtu_elementorder << dir.vtudir
                                 << "/pgv_tp"
                                 << "_s" << iSetup
                                 << "_i" << it_counter;
                outputGridviewIndexToDGF( *pgv_tp, inputdata, vtu_elementorder.str() );
#endif

                pool_gv_tp[ pool_lookup[iSetup][t_meas] ] = pgv_tp;

                Dune::shared_ptr<CON_TP> pcon_tp = Dune::make_shared<CON_TP>();
                pool_con_tp[ pool_lookup[iSetup][t_meas] ] = pcon_tp;

                Dune::shared_ptr<FEM_HYPER> pfem_hyper;
#ifdef USE_FEM
                pfem_hyper = Dune::make_shared<FEM_HYPER>(*pgv_tp);
#else
                pfem_hyper = Dune::make_shared<FEM_HYPER>();
#endif
                Dune::shared_ptr<GFS_TP> pgfs_tp = Dune::make_shared<GFS_TP>( *pgv_tp, 
                                                                              pfem_hyper, 
                                                                              pcon_tp );
                pool_gfs_tp[ pool_lookup[iSetup][t_meas] ] = pgfs_tp;

                Dune::shared_ptr<CON_CG> pcon_cg = Dune::make_shared<CON_CG>();
                pool_con_cg[ pool_lookup[iSetup][t_meas] ] = pcon_cg;
                Dune::shared_ptr<FEM_CG> pfem_cg
                  = Dune::make_shared<FEM_CG>(*pgv_tp);
                Dune::shared_ptr<GFS_CG> pgfs_cg = Dune::make_shared<GFS_CG>( *pgv_tp, 
                                                                              pfem_cg, 
                                                                              pcon_cg );
                pool_gfs_cg[ pool_lookup[iSetup][t_meas] ] = pgfs_cg;

              }


              /*
               * solute M0
               */
              if((t_meas==2 || t_meas==3 || t_meas==5) && soluteM0_old[pool_lookup[iSetup][t_meas]]==NULL){

                /*
                Dune::shared_ptr<GV_TP> pgv_tp = Dune::make_shared<GV_TP>( pool_grid[ pool_lookup[iSetup][t_meas] ]->leafGridView());
                Dune::shared_ptr<CON_CG> pcon_cg = Dune::make_shared<CON_CG>();
                Dune::shared_ptr<GFS_CG> pgfs_cg = Dune::make_shared<GFS_CG>( *pgv_tp, 
                                                                              *pfem_cg, 
                                                                              *pcon_cg );
                                                                              */
                
                // load the old forward solution for m0
                soluteM0_old[pool_lookup[iSetup][t_meas]] = new VCType_CG( *(pool_gfs_cg[pool_lookup[iSetup][t_meas]]) 
                                                                           , 0.0 );

                logger << "DEBUG: read /vcM0_old" << std::endl;
                HDF5Tools::
                  read_BackendVector_from_HDF5( 
                                               *(pool_gfs_cg[pool_lookup[iSetup][t_meas]])
                                               , inputdata
                                               , dir.vcM0_old_h5file[iSetup] 
                                               , "/vcM0_old"
                                               , *(soluteM0_old[pool_lookup[iSetup][t_meas]])
                                               , PRESERVE_STRUCTURE
                                               , baselevel
                                                );


#ifdef DEBUG_PLOT_M0_OLD
                std::stringstream vtu_m0_old;
                vtu_m0_old << dir.vtudir
                           << "/m0_old"
                           << "_s" << iSetup
                           << "_i" << it_counter;
                VTKPlot::output2vtu( *(pool_gfs_cg[pool_lookup[iSetup][t_meas]]),
                                     *(soluteM0_old[pool_lookup[iSetup][t_meas]]),
                                     vtu_m0_old.str(),
                                     "M0_old",
                                     inputdata.verbosity,
                                     true,
                                     0
                                     );
#endif // DEBUG_PLOT_M0_OLD

              }
              /*
               * solute M1
               */
              if(( t_meas==3 || t_meas==5) && soluteM1_old[pool_lookup[iSetup][t_meas]]==NULL){
                // load the old head data

                /*
                Dune::shared_ptr<GV_TP> pgv_tp = Dune::make_shared<GV_TP>( pool_grid[ pool_lookup[iSetup][t_meas] ]->leafGridView());                
                Dune::shared_ptr<CON_CG> pcon_cg = Dune::make_shared<CON_CG>();
                Dune::shared_ptr<GFS_CG> pgfs_cg = Dune::make_shared<GFS_CG>( *pgv_tp, 
                                                                              *pfem_cg, 
                                                                              *pcon_cg );
                                                                              */
                
                soluteM1_old[pool_lookup[iSetup][t_meas]]= new VCType_CG( *(pool_gfs_cg[pool_lookup[iSetup][t_meas]])
                                                                          , 0.0 );
                HDF5Tools::
                  read_BackendVector_from_HDF5( 
                                               *(pool_gfs_cg[pool_lookup[iSetup][t_meas]])
                                               , inputdata
                                               , dir.vcM1_old_h5file[iSetup] 
                                               , "/vcM1_old"
                                               , *(soluteM1_old[pool_lookup[iSetup][t_meas]])
                                               , PRESERVE_STRUCTURE
                                               , baselevel
                                                );


              }


              /*
               * solute heat M0 and heat M1
               *

              if(( t_meas==4 || t_meas==5) && 
                 heatM0_old[pool_lookup[iSetup][t_meas]]==NULL &&
                 heatM1_old[pool_lookup[iSetup][t_meas]]==NULL
                 ){
                
                heatM0_old[pool_lookup[iSetup][t_meas]]= new VCType_CG( *(pool_gfs_cg[pool_lookup[iSetup][t_meas]])
                                                                        , 0.0 );
                HDF5Tools::
                  read_BackendVector_from_HDF5( 
                                               *(pool_gfs_cg[pool_lookup[iSetup][t_meas]])
                                               , inputdata
                                               , dir.vcheatM0_old_h5file[iSetup]
                                               , "/vcheatM0_old"
                                               , *(heatM0_old[pool_lookup[iSetup][t_meas]])
                                                );
                
                heatM1_old[pool_lookup[iSetup][t_meas]]= new VCType_CG( *(pool_gfs_cg[pool_lookup[iSetup][t_meas]])
                                                                        , 0.0 );
                HDF5Tools::
                  read_BackendVector_from_HDF5( 
                                               *(pool_gfs_cg[pool_lookup[iSetup][t_meas]])
                                               , inputdata
                                               , dir.vcheatM1_old_h5file[iSetup]
                                               , "/vcheatM1_old"
                                               , *(heatM1_old[pool_lookup[iSetup][t_meas]])
                                                );
              }
              */

            }
          }
          //}


               
          // lnK SENSITIVITIES!
          if(orig_measurements.nMeasPerSetupPerType(iSetup,0)){

            //std::cout << "outside: " << nCellsExt[0][1] << std::endl;
            
            // important to use the pool_lookup vector to get the right pool_gv
            lnK_sensitivities<GV_GW,
                              DIR,
                              MEASLIST,
                              IDT
                              >( pool_gv_0[ pool_lookup[iSetup][0] ],
                                 inputdata,
                                 dir,
                                 iSetup,
                                 it_counter,
                                 orig_measurements,
                                 nCellsExt,
                                 pool_Lambdas[ pool_lookup[iSetup][0] ],
                                 pool_X[ pool_lookup[iSetup][0] ],
                                 Y_old,
                                 CommunicatorPool[ pool_lookup[iSetup][0] ],
                                 //OUTPUT
                                 JX,  // only on P0 
                                 J_times_Y_old // only on P0
                                 );
        
          }
          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());

                
          // head sensitivities
          if(orig_measurements.nMeasPerSetupPerType(iSetup,1)){

            head_sensitivities<GV_GW,
                               GFS_GW,
                               VCType_GW,
                               PGV,
                               DIR,
                               MEASLIST,
                               IDT,
                               SDT,
                               YFG>
              ( pool_gv_0[ pool_lookup[iSetup][0] ],
                *(pool_gfs_gw[pool_lookup[iSetup][1]]),
                *(head_old[pool_lookup[iSetup][1]]),
                pRootGridView,
                dir,
                orig_measurements,
                inputdata,
                setupdata,
                nCellsExt,
                pool_Lambdas[pool_lookup[iSetup][1]],
                pool_X[pool_lookup[iSetup][1]],
                Y_old,
                *(YfieldGenerator_old[pool_lookup[iSetup][1]]),
                it_counter,
                helper,
                CommunicatorPool[pool_lookup[iSetup][1]],
                // OUTPUT
                JX, // only on P0
                J_times_Y_old // only on P0
                );

          }

          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());


          // M0 sensitivities
      
          if(orig_measurements.nMeasPerSetupPerType(iSetup,2)){
        
            // important to use the pool_lookup vector to get the right pool gridfunction space
            m0_sensitivities<POOL_GRID,
                             GFS_GW,
                             GFS_TP,
                             GFS_CG,
                             VCType_GW,
                             VCType_CG,
                             MEASLIST,
                             YFG,
                             DIR,
                             IDT,
                             SDT
                             >
              (
               *(pool_grid[ pool_lookup[iSetup][2] ]),
               *(pool_gfs_gw[ pool_lookup[iSetup][2] ]),
               *(pool_gfs_tp[ pool_lookup[iSetup][2] ]),
               *(pool_gfs_cg[ pool_lookup[iSetup][2] ]),

               inputdata,
               setupdata,
               dir,

               orig_measurements,
               measurements,

               nCellsExt,
               pool_Lambdas[ pool_lookup[iSetup][2] ], // important to use the pool_lookup vector to get the right Lambdas
               pool_X[ pool_lookup[iSetup][2] ],

               Y_old,
               *(YfieldGenerator_old[ pool_lookup[iSetup][2] ]), // important to use the pool_lookup vector to get the right YfieldGenerator_old

               *(head_old[ pool_lookup[iSetup][2] ]),
               *(soluteM0_old[ pool_lookup[iSetup][2] ]),

               it_counter,
               helper,
               CommunicatorPool[pool_lookup[iSetup][2]], // important to use the pool_lookup vector to get the right communicator pool
               // OUTPUT
               JX, // only on P0
               J_times_Y_old, // only on P0
               SV_increased
               );   
          }              
          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());

          // M1 sensitivities

          if(orig_measurements.nMeasPerSetupPerType(iSetup,3)){

            m1_sensitivities<POOL_GRID,
                             GFS_GW,
                             GFS_TP,
                             GFS_CG,
                             VCType_GW,
                             VCType_CG,
                             MEASLIST,
                             YFG,
                             DIR,
                             IDT,
                             SDT
                             >
              (
               *(pool_grid[ pool_lookup[iSetup][3] ]),
               *(pool_gfs_gw[pool_lookup[iSetup][3]]),
               *(pool_gfs_tp[pool_lookup[iSetup][3]]),
               *(pool_gfs_cg[pool_lookup[iSetup][3]]),

               inputdata,
               setupdata,
               dir,

               orig_measurements,
               measurements,
               
               nCellsExt,
               pool_Lambdas[pool_lookup[iSetup][3]], // important to use the pool_lookup vector to get the right Lambdas
               pool_X[pool_lookup[iSetup][3]],

               Y_old,
               *(YfieldGenerator_old[pool_lookup[iSetup][3]]), // important to use the pool_lookup vector to get the right YfieldGenerator_old

               *(head_old[pool_lookup[iSetup][3]]),
               *(soluteM0_old[pool_lookup[iSetup][3]]),
               *(soluteM1_old[pool_lookup[iSetup][3]]),
               
               it_counter,
               helper,
               CommunicatorPool[pool_lookup[iSetup][3]], // important to use the pool_lookup vector to get the right communicator pool
               // OUTPUT
               JX, // only on P0
               J_times_Y_old // only on P0
               );
          }
          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());


           
          // heat sensitivities
          if(orig_measurements.nMeasPerSetupPerType(iSetup,4)){

            heat_sensitivities<POOL_GRID,
                               GFS_GW,
                               GFS_TP,
                               GFS_CG,
                               VCType_GW,
                               VCType_CG,
                               MEASLIST,
                               YFG,
                               DIR,
                               IDT,
                               SDT
                               >
              (
               *(pool_grid[ pool_lookup[iSetup][4] ]),
               *(pool_gfs_gw[pool_lookup[iSetup][4]]),
               *(pool_gfs_tp[pool_lookup[iSetup][4]]),
               *(pool_gfs_cg[pool_lookup[iSetup][4]]),

               inputdata,
               setupdata,
               dir,

               orig_measurements,
               measurements,
               
               nCellsExt,
               pool_Lambdas[pool_lookup[iSetup][4]], // important to use the pool_lookup vector to get the right Lambdas
               pool_X[pool_lookup[iSetup][4]],

               Y_old,
               *(YfieldGenerator_old[pool_lookup[iSetup][4]]), // important to use the pool_lookup vector to get the right YfieldGenerator_old

               *(head_old[pool_lookup[iSetup][4]]),

               it_counter,
               helper,
               CommunicatorPool[pool_lookup[iSetup][4]], // important to use the pool_lookup vector to get the right communicator pool
               // OUTPUT
               JX, // only on P0
               J_times_Y_old // only on P0
               );

          }
          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());




                
          // GP sensitivities
          if(orig_measurements.nMeasPerSetupPerType(iSetup,5)){
            geoelectrical_potential_sensitivities<
              POOL_GRID,
              GFS_GW,
              GFS_TP,
              GFS_CG,
              VCType_GW,
              VCType_CG,
              MEASLIST,
              YFG,
              DIR,
              IDT,
              SDT
              >(
                *(pool_grid[ pool_lookup[iSetup][5] ]),
                *(pool_gfs_gw[pool_lookup[iSetup][5]]),
                *(pool_gfs_tp[pool_lookup[iSetup][5]]),
                *(pool_gfs_cg[pool_lookup[iSetup][5]]),

                inputdata,
                setupdata,
                dir,

                orig_measurements,
                measurements,

                nCellsExt,
                pool_Lambdas[pool_lookup[iSetup][5]], // important to use the pool_lookup vector to get the right Lambdas
                pool_X[pool_lookup[iSetup][5]],
                
                Y_old,
                *(YfieldGenerator_old[ pool_lookup[iSetup][5] ]),
                *(log_electricalConductivity_pool[ pool_lookup[iSetup][5] ]),
                *(log_kappafield_pool[ pool_lookup[iSetup][5] ]),

                *(head_old[ pool_lookup[iSetup][5] ]),
                *(soluteM0_old[ pool_lookup[iSetup][5] ]),
                *(soluteM1_old[ pool_lookup[iSetup][5] ]),

                it_counter,
                helper,
                CommunicatorPool[ pool_lookup[iSetup][5] ], // important to use the pool_lookup vector to get the right communicator pool
                // OUTPUT
                JX, // only on P0
                J_times_Y_old // only on P0
                );
          }
                
          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());
                

     
          //free the old vectors
          for(UINT iPool=0; iPool<number_of_MPIPools; iPool++){
            if(head_old[iPool]!=NULL){
              delete(head_old[iPool]);
              head_old[iPool]=NULL;
            }
            if(soluteM0_old[iPool]!=NULL){
              delete(soluteM0_old[iPool]);
              soluteM0_old[iPool]=NULL;
            }
            if(soluteM1_old[iPool]!=NULL){
              delete(soluteM1_old[iPool]);
              soluteM1_old[iPool]=NULL;
            }
          }
          if(helper.rank()!=0){
            Y_old.resize(0);   
          } else if(CR>-1){

            HDF5Tools::
              read_sequential_from_HDF5( Y_old
                                         , "/Y_old"
                                         , local_count
                                         , local_offset
                                         , dir.Y_old_h5file
                                         , inputdata
                                         );
          }
                
        }// END: for(UINT iSetup=0; iSetup<nSetups; iSetup++)
        
        //collect JX and J_times_Y_old on P0
        REAL tmp;
        for(UINT ii=0;ii<nMeas;ii++){
          for(UINT jj=0; jj<nzones; jj++){
            tmp=0.0;
            MPI_Reduce(&(JX[ii][jj]),&(tmp),1,MPI_DOUBLE,MPI_SUM,0,helper.getCommunicator());
            if(helper.rank()==0)
              JX[ii][jj]=tmp;
          }
          tmp=0.0;
          MPI_Reduce(&(J_times_Y_old[ii]),&(tmp),1,MPI_DOUBLE,MPI_SUM,0,helper.getCommunicator());
          if(helper.rank()==0)
            J_times_Y_old[ii]=tmp;
        }
 

        logger << "DEBUG:" << std::endl;
        for(UINT ii=0;ii<nMeas;ii++){
          logger << "J_times_Y_old[" << ii << "] = " << J_times_Y_old[ii] << std::endl;
        }
        for(UINT ii=0;ii<nMeas;ii++){
          logger << "JX[" << ii << "][0] = " << JX[ii][0] << std::endl;
        }

        //wait until all sensitivities are calculated! VERY VERY VERY IMPORTANT
        if(helper.size()>1)
          MPI_Barrier(helper.getCommunicator());


      }// END: calculate





   
      void calculate_JQJ(){
        //parallel calculation of JQJ (parallel on MPI_COMM_WORLD)
        Dune::Timer watch;

        Vector<UINT> local_count,local_offset;
        Vector<REAL> vecSensitivity;
        Vector<REAL> vecJQ;

        UINT nJQJ;
        UINT shift;
        if(inputdata.parallel_fine_tuning.JQJ_max >= (UINT)helper.size()){
          nJQJ=helper.size();
          shift=1;
        }else{
          nJQJ=inputdata.parallel_fine_tuning.JQJ_max;
          shift=std::floor(helper.size()/nJQJ);
        }
            
        
        logger<<"JQJ ( "<<nJQJ<<" processors working)...(of "<<nMeas<<" )"<<std::endl;
        JQJ.resize(0);
        JQJ.resize(nMeas);
        
        //parallel temporary matrix!
        std::vector< std::vector< REAL > > JQJ_p(nMeas);
       
       
        if( helper.rank()==0 && General::verbosity >= VERBOSITY_INVERSION )
          std::cout << "calculate_JQJ: Upper triangular matrix loop ... " << std::endl;

        for( UINT iPoint=0; iPoint < nMeas; iPoint++ ){
          //resize the matrix!
          JQJ_p[ iPoint ].resize( nMeas,0.0 );
          JQJ[ iPoint ].resize( nMeas, 0.0 );
          
          if((int)((iPoint%nJQJ)*shift)==helper.rank()){ // distribute the work to the available processors!
            for(UINT iParts=0; iParts<inputdata.parallel_fine_tuning.JQJ_slices; iParts++){

              //std::cout << "calculate_JQJ(): Process " << helper.rank() 
              //          << " reading JQ_" << iPoint
              //          << std::endl;
              HDF5Tools::read_sequential_from_HDF5(
                                          vecJQ
                                          , "/JQ"
                                          , local_count
                                          , local_offset
                                          , dir.JQ_h5file[iPoint]
                                          , inputdata
                                          , iParts
                                          , inputdata.parallel_fine_tuning.JQJ_slices
                                          );


              // Count sequential reading time of sensitivity files only for the root process
              // because it has the highest number of files to read and the timer for the sequential reads 
              // on all other processes cannot be synchronized
              bool bSwitchTimerON = true;
              if( helper.rank() > 0 )
                bSwitchTimerON = false; 

              // inner loop forming the upper triangle of the matrix JQJ:
              for( UINT jPoint=iPoint; jPoint < nMeas; jPoint++ ){
                //std::cout << "calculate_JQJ(): Process " << helper.rank() 
                //          << " reading sens_" << iPoint
                //          << std::endl;

                //load Sensitivity!  
                HDF5Tools::read_sequential_from_HDF5(
                                            vecSensitivity
                                            , "/Sensitivity"
                                            , local_count
                                            , local_offset
                                            , dir.Sensitivity_h5file[jPoint]
                                            , inputdata
                                            , iParts
                                            , inputdata.parallel_fine_tuning.JQJ_slices
                                            , bSwitchTimerON
                                            );

                // inner product of two vectors:
                JQJ_p[ iPoint ][ jPoint ] = vecJQ * vecSensitivity;

              } // END: for( UINT jPoint=iPoint; jPoint < nMeas; jPoint++ )

            } // END: for(UINT iParts=0; iParts<inputdata.parallel_fine_tuning.JQJ_slices; iParts++)

          } // END:  if((int)((iPoint%nJQJ)*shift)==helper.rank()) 

        }

        // All processes please wait until all processes are done with their JQJ_p
        MPI_Barrier(helper.getCommunicator());


        watch.reset();
        // MPI REDUCE (for each element of the matrix)-> for the whole matrix does not work (maybe row-wise will work)
        for(UINT ii=0;ii<nMeas;ii++)
          for(UINT jj=ii;jj<nMeas;jj++)
            MPI_Reduce( &(JQJ_p[ii][jj]),
                        &(JQJ[ii][jj]),
                        1,
                        MPI_DOUBLE,
                        MPI_SUM,
                        0,
                        helper.getCommunicator() );

        //mirror upper half to the lower half!
        for(UINT ii=0;ii<nMeas-1;ii++)
          for(UINT jj=ii+1;jj<nMeas;jj++)
            JQJ[jj][ii] = JQJ[ii][jj];

        General::log_elapsed_time( watch.elapsed(),
                                   helper.getCommunicator(),
                                   General::verbosity,
                                   "REDIST",
                                   "JQJ reduce all to root and mirroring upper half to lower half" );

        if(helper.size()>1)
          MPI_Barrier(helper.getCommunicator());
        /*        
                  logger<<"JQJ : "<<std::endl;
                  for( UINT iPoint=0; iPoint < nMeas; iPoint++ ){
                  for( UINT jPoint=0; jPoint < nMeas; jPoint++ ){
                  logger<<JQJ[iPoint][jPoint]<<" ";
                  }
                  logger<<std::endl;
                  }
        */

      } // END calculate_JQJ



      void cleanupFiles(){
        for(UINT ii=0;ii<nMeas;ii++){
          General::deleteFile( dir.Sensitivity_h5file[ ii ] );
        }
      }
        

    
      template<typename DenseMatrix>
      void get_JQJ(DenseMatrix & M){
        for(UINT ii=0; ii<nMeas; ii++)
          for(UINT jj=0; jj<nMeas; jj++)
            M(ii,jj) = JQJ[ii][jj];
      }



    
      REAL get_JX(UINT ii, UINT jj){
        return JX[ii][jj];    
      }




    
      REAL get_J_times_Y_old(UINT ii){
        return J_times_Y_old[ii];
      }
    
    
      /*
       * evaluate the prior term of the objective function!
       *
       * Remark: call only on P0(so far ksi_try is only known there)
       *
       */
      REAL calculate_prior_term( const Vector<REAL>& ksi_try,
                                 const int wCounter
                                 ) {

        Dune::Timer watch;

        REAL L_prior = 0.0;

        if( wCounter > 1 ){
          // weighting is active: Interpolate Gyy component-wise.
          REAL ww = std::pow(0.5,REAL(wCounter-1));

          for( UINT iPoint=0; iPoint < nMeas; iPoint++ ){
            REAL sum1 = 0.0;

            for( UINT jPoint=0; jPoint < nMeas; jPoint++ ){

              REAL localGyy = (1.0 - ww) * Gyy_OLD( jPoint , iPoint )
                + ww * Gyy_NEW( jPoint, iPoint );

              sum1 += ksi_try[ jPoint ] * localGyy;
            }

            REAL contribution = sum1 * ksi_try[ iPoint ];

            if( inputdata.verbosity >= VERBOSITY_L_CONTRIB ){
              std::cout << "Lp: iMeas = " << iPoint
                        << " ksi*Gyy*ksi = " << contribution
                        << std::endl;
            }

            L_prior += contribution;      // (6.30) von Nowak
          }

        }
        else {

          for( UINT iPoint=0; iPoint < nMeas; iPoint++ ){
            REAL sum1 = 0.0;

            for( UINT jPoint=0; jPoint < nMeas; jPoint++ ){
              Gyy_NEW(jPoint,iPoint) = JQJ[ jPoint ][ iPoint ];

              for(UINT ii=0; ii<nzones; ii++) {
                REAL R_bb=inputdata.yfield_properties.zones[ii].qbb_y;
                Gyy_NEW( jPoint , iPoint ) += R_bb * JX[ iPoint ][ii] * JX[ jPoint ][ii];
              }
              sum1 += ksi_try[ jPoint ] * Gyy_NEW( jPoint , iPoint );
            }

            REAL contribution = sum1 * ksi_try[ iPoint ];

            if( inputdata.verbosity >= VERBOSITY_L_CONTRIB ){
              std::cout << "Lp: iMeas = " << iPoint
                        << " ksi*Gyy*ksi = " << contribution
                        << std::endl;
            }

            L_prior += contribution;      // (6.30) von Nowak
          }

        }

        General::log_elapsed_time( watch.elapsed(),
                                   General::verbosity,
                                   "EVAL",
                                   "Compute L_p" );
        return L_prior;  
      }
    

      void store_Gyy(){
        Gyy_OLD = Gyy_NEW;
      }



      template<typename GV>
      void vtu_sensitivity(const GV& gv_0, UINT it_counter){
        Vector<UINT> local_count,local_offset;
        for(UINT ii=0; ii<nMeas; ii++){
          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());

          std::stringstream vtu_file;
          vtu_file << dir.Sensitivity_vtu_prefix << "_m" << ii << "_i" << it_counter;

          YFG yfg_Sensitivity( inputdata, dir, helper.getCommunicator() );
          VTKPlot::output_hdf5data_to_gridview( gv_0,
                                                inputdata,
                                                vtu_file.str(),
                                                dir.Sensitivity_h5file[ii],
                                                "/Sensitivity",
                                                yfg_Sensitivity
                                                );



        }
      }


      template<typename GV>
      void vtu_JQ(const GV& gv_0,UINT it_counter){
        Vector<UINT> local_count,local_offset;
        for(UINT ii=0; ii<nMeas; ii++){

          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());

          std::stringstream vtu_file;
          vtu_file << dir.JQ_vtu_prefix << "_m" << ii << "_i" << it_counter;

          YFG yfg_JQ( inputdata, dir, helper.getCommunicator() );
          VTKPlot::output_hdf5data_to_gridview( gv_0,
                                                inputdata,
                                                vtu_file.str(),
                                                dir.JQ_h5file[ii],
                                                "/JQ",
                                                yfg_JQ
                                                );
        }
      }
    };




  } // namespace GeoInversion  

} // namespace Dune

#endif // #ifndef SENSITIVITYCLASS_HH
