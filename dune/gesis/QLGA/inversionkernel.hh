#ifndef DUNE_GESIS_INVERSION_KERNEL_HH
#define DUNE_GESIS_INVERSION_KERNEL_HH

#include <boost/math/special_functions/gamma.hpp>

//          Copyright Joe Coder 2004 - 2006.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <dune/gesis/QLGA/SensitivityClass.hh>
#include <dune/gesis/QLGA/objectivefunction.hh>

#ifdef ESTIMATION_VARIANCE
#include <dune/gesis/QLGA/estimation_variance.hh>
#endif


namespace Dune {
  namespace Gesis {

    // grid view on the Dune grid with all processors (global grid)
    // helper: for the global MPI communication
    template<typename GRID,
             typename PGV,
             typename GFS_GW,
             typename GFS_TP,
             typename GFS_CG,
             typename MEASLIST,
             typename DIR,
             typename IDT,
             typename YFG
             >
    double inversionkernel( GRID& theGrid,
                            const PGV pRootGridView,
                            const GFS_GW& gfs_gw,
                            const IDT& inputdata,
                            MEASLIST& orig_measurements,
                            const std::vector< Vector<UINT> >& nCellsExt,
                            DIR& dir,
                            const YFG& log_electricalConductivity,
                            const YFG& log_kappafield,
                            const Dune::MPIHelper& helper,
                            const INT& CR=-1
                            ){

      typedef typename GFS_GW::Traits::GridViewType GV_GW;

      const UINT nAllCells = inputdata.domain_data.nCells.componentsproduct();

      /*
       * NO GLOBALREFINE and only yasp grid!!
       */
#ifndef USE_YASP
      if(helper.rank()==0)
        std::cerr<<std::endl<<"Inversion works ONLY with yasp grid! --> QUIT"<<std::endl<<std::endl;
      logger<<std::endl<<"Inversion works ONLY with yasp grid! --> QUIT"<<std::endl<<std::endl;
      return -1.0;
#endif


#ifdef USE_YASP
      UINT baselevel = inputdata.domain_data.yasp_baselevel;
      if(baselevel>1){
        if(helper.rank()==0)
          std::cerr<<std::endl<<"Inversion is tested ONLY with globalrefine<2 (set the virtual_nCells to the needed values)! --> QUIT"<<std::endl<<std::endl;
        logger<<std::endl<<"Inversion is tested ONLY with globalrefine<2 (set the virtual_nCells to the needed values)! --> QUIT"<<std::endl<<std::endl;
        return -1.0;
      }
#else
      UINT baselevel = 0;
      std::cout << "Inversion currently works only with YASP !!!" << std::endl;
      return -1.0;
#endif

      /*
       * END: NO GLOBALREFINE and only yasp grid
       */


      Dune::Timer watch;
      watch.reset();

      std::ofstream likelihood_ostream;
      std::ofstream measurement_ostream;
      if(helper.rank()==0){
        if(CR>-1){
          likelihood_ostream.open (dir.likelihood_cond_file[CR].c_str());
          measurement_ostream.open (dir.measurement_cond_file[CR].c_str());
        }else{
          likelihood_ostream.open (dir.likelihood_file.c_str());
          measurement_ostream.open (dir.measurement_file.c_str());
        }
        if (likelihood_ostream.is_open()){
          likelihood_ostream<<"Values of the objective function!"<<std::endl<<"FORMAT: Iteration number:value"<<std::endl;
        }else{
          std::cerr<<std::endl<<"UNABLE TO OPEN: "<<dir.likelihood_file<<", --> QUIT"<<std::endl<<std::endl;
          return -1.0;
        }
        if (measurement_ostream.is_open()){
          measurement_ostream<<"Measurement output!"<<std::endl<<"FORMAT: measurement #:original value:estimated value: abs error: rel error"<<std::endl;
        }else{
          std::cerr<<std::endl<<"UNABLE TO OPEN: "<<dir.likelihood_file<<", --> QUIT"<<std::endl<<std::endl;
          return -1.0;
        }

      }


      if(helper.rank()==0 && inputdata.verbosity>2)
        std::cout << std::endl << "inversionkernel()" << std::endl;
      logger << std::endl << "inversionkernel()" << std::endl;

      /*
       * variables and definitions
       */
      enum{dim=GV_GW::dimension};
      Vector<UINT> local_count,local_offset;
      Vector<REAL> Y_u(0);
      if(CR>-1){
        logger<<"Calculate conditional realization #"<<CR<<std::endl;
        if(helper.rank()==0){
          HDF5Tools::h5g_Read( Y_u
                               , dir.unconditionalField_h5file
                               , "/YField"
                               , local_offset
                               , local_count
                               , inputdata
                               ); // TODO: Ask Ronnie: Why not parallel read?
          // Answer: Becasue this field is only needed on P0. The inversion kernel is only done on P0.
          //         Solving the kriging system and then updating the YField!
        }
      }

      UINT nMeas=orig_measurements.nMeasurements();
      UINT nzones=inputdata.yfield_properties.nz;
      Vector<REAL> beta_star(nzones); // mean of the trend coefficients
      Vector<REAL> R_bb(nzones);     // uncertainty about the mean

      // set the values from the inputdata
      for(UINT ii=0; ii<nzones; ii++){
        beta_star[ii] = inputdata.yfield_properties.zones[ii].beta;
        R_bb[ii] = inputdata.yfield_properties.zones[ii].qbb_y;
      }

      //initialize the measurement list for the inversion
      // it is always a synthetic measurement -> true
      MEASLIST sim_measurements(inputdata,true);

      typedef SensitivityClass
        <GRID,
         PGV,
         GFS_GW,
         GFS_TP,
         GFS_CG,
         MEASLIST,
         YFG,
         DIR,
         IDT
         > SensCl;

      //initalize the sensitivity class
      SensCl sensitivity( theGrid,
                          pRootGridView,
                          inputdata,
                          nCellsExt,
                          dir,
                          orig_measurements,
                          helper );


      GV_GW gv_0 = theGrid.levelGridView(0);

      //const GV_GW& gv_gw = gfs_gw.gridView();


      // Frag Ronnie:
      // Ist es nicht genug, "log_electricalConductivity"
      // und "log_kappafield"
      // einmal in "initialize_DuneGrid.hh" zu lesen?

      /*
      // define the sigma0_wrapper and kappa_wrapper! works on all processors!
      YFG log_electricalConductivity_world(inputdata,dir,helper.getCommunicator());
      YFG log_kappafield_world(inputdata,dir,helper.getCommunicator());

      if(inputdata.problem_types.moments_geoeletric_potential_inversion){
      Vector<REAL> logsigma0;
      HDF5Tools::h5g_pRead( gv_gw
      , inputdata
      , logsigma0
      , "/logsigma0"
      , local_count
      , local_offset
      , dir.logsigma0_h5file
      );
      // export the data to the sigma0 Generator
      if(helper.size()>1)
      log_electricalConductivity_world.parallel_import_from_local_vector(logsigma0, local_count, local_offset );
      else
      log_electricalConductivity_world.import_from_vector( logsigma0 );


      Vector<REAL> logkappa;
      HDF5Tools::
      h5g_pRead( gv_gw
      , inputdata
      , logkappa
      , "/logkappa"
      , local_count
      , local_offset
      , dir.logkappa_h5file
      );
      // export the data to the sigma0 Generator
      if(helper.size()>1)
      log_kappafield_world.parallel_import_from_local_vector(logkappa, local_count, local_offset );
      else
      log_kappafield_world.import_from_vector( logkappa );
      }
      */


      //const FieldWrapperType log_sigma0_wrapper_world( inputdata, setupdata, log_electricalConductivity_world );
      //const FieldWrapperType log_kappa_wrapper_world( inputdata, setupdata, log_kappafield_world );

      // chi 2 test for acceptance
      REAL L_accept_value = 2.0 * boost::math::gamma_p_inv( 0.5*nMeas, inputdata.inversion_parameters.L_accept_confidence_interval );

      // set the filenames for the needed files (now when the number of measurements is known!)
      dir.set_inversion(nMeas);

      //vector for the measurement error
      std::vector< REAL > MeasurementErrorsV;

      // only on P0 (global rank)
      if(helper.rank()==0){
        //set the vector of the measurement error
        MeasurementErrorsV.resize( nMeas, 0.0 );
        for(UINT ii=0;ii<nMeas;ii++){
          REAL epsilon = orig_measurements[ii].value * orig_measurements[ ii ].rel_error  + orig_measurements[ ii ].abs_error;
          MeasurementErrorsV[ ii ] = epsilon*epsilon;
        }
      }

      Vector<REAL> Y_old(0);
      Vector<REAL> ksi_old( nMeas , 0.0 );
      REAL L_prior_ini = 1e+12;  // Let the initial value for L_prior better be very large!
      Vector<REAL> beta_old(nzones);
      for(UINT ii=0; ii<nzones; ii++)
        beta_old[ii] = inputdata.yfield_properties.zones[ii].beta;
      std::vector< Vector<REAL> >  Xzones; // zonation matrix
      std::vector<REAL> L_objective;

      /*
       * INITIAL GUESS
       * load an existing initial guess or generate a homogeneous guess
       */
      // and do something on the P0
      if(helper.rank()==0){ // PO do some work
        //read the zonation matrix on P0
        Xzones.resize(nzones);
        for(UINT jj=0; jj<nzones; jj++)
          HDF5Tools::h5_ReadDirect( Xzones[jj],
                                    dir.zonation_matrix[jj],
                                    "/X" );

        /*
         * get the initial log K field
         *   -> either homogeneous or read an existing one!!!
         */

        if ( inputdata.problem_types.using_existing_Yold
             // && General::bFileExists( dir.ksi_old_h5file )
             // && General::bFileExists( dir.L_prior_h5file )
             && General::bFileExists( dir.beta_old_h5file )
             && ( General::bFileExists( dir.Y_estimated_h5file )
                  || General::bFileExists( dir.Y_old_h5file )
                  ) ) {  // read an existing Y field!

          if( General::bFileExists( dir.Y_estimated_h5file ) ){

            if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION )
              std::cout << "For the initial guess: Read /Y_est from "
                        << dir.Y_estimated_h5file << std::endl;
            logger << "For the initial guess: Read /Y_est from "
                   << dir.Y_estimated_h5file << std::endl;

            HDF5Tools::h5g_Read( Y_old
                                 , dir.Y_estimated_h5file
                                 , "/Y_est"
                                 , local_offset
                                 , local_count
                                 , inputdata
                                 );
          }
          else {

            if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION )
              std::cout << "For the initial guess: Read /Y_old from "
                        << dir.Y_old_h5file << std::endl;
            logger << "For the initial guess: Read /Y_old from "
                   << dir.Y_old_h5file << std::endl;

            HDF5Tools::h5g_Read( Y_old
                                 , dir.Y_old_h5file
                                 , "/Y_old"
                                 , local_offset
                                 , local_count
                                 , inputdata
                                 );
          }


          if( Y_old.size() != nAllCells ){ // stop if the size of the Y field has changed!
            std::cout<<"YOU CANNOT USE THE OLD Y_OLD! (because dimension has changed!)"<<std::endl;
            logger<<"YOU CANNOT USE THE OLD Y_OLD! (because dimension has changed!)"<<std::endl;
            return -1.0;
          }


          if( General::bFileExists( dir.ksi_old_h5file ) ){
            HDF5Tools::h5_ReadDirect( ksi_old,
                                      dir.ksi_old_h5file,
                                      "/ksi_old" );
            std::cout << "ksi_old read from file " << dir.ksi_old_h5file
                      << std::endl;
          }
          else {
            std::cout << "ksi_old initialized to 0 vector." << std::endl;
          }


          Vector<REAL> tmp;
          if( General::bFileExists( dir.L_prior_h5file ) ){
            HDF5Tools::h5_ReadDirect( tmp,
                                      dir.L_prior_h5file,
                                      "/L_prior" );
            L_prior_ini = tmp[0];
          }
          else {
            L_prior_ini = 1e+32; // Let the initial value for L_prior better be very large!
          }

          //std::cout << "DEBUG: L_prior_ini read from file = " << std::endl;
          //std::cout << L_prior_ini << std::endl;

          HDF5Tools::h5_ReadDirect( tmp,
                                    dir.beta_old_h5file,
                                    "/beta_old" );
          for(UINT ii=0; ii<nzones; ii++)
            beta_old[ii]=tmp[ii];
          //std::cout << "DEBUG: beta_old read from file = " << std::endl;
          //std::cout << beta_old << std::endl;

          if(CR>-1){
            for(UINT ii=0; ii<Y_u.size(); ii++)
              Y_old[ii]+=Y_u[ii];
          }

          // write the inital guess to hdf5 (as Y_old)
          HDF5Tools::h5g_Write(
                               Y_old
                               , dir.Y_old_h5file
                               , "/Y_old"
                               , inputdata
                               );

          /*
            Vector<REAL> Y_old_Well( Y_old );
            General::add_wells_to_field( inputdata, Y_old_Well );

            HDF5Tools::
            h5g_Write(
            Y_old_Well
            , "/Y_old_Well"
            , inputdata
            , dir.Y_old_Well_h5file
            , "Y_old_Well"
            , dir.Y_old_xmffile
            );
          */

        } else { // generate a homogeneous intial guess!

          if(inputdata.problem_types.using_existing_Yold){
            logger<<"MISSING FILE for using existing Y_old as initial guess!"<<std::endl<<"FALL BACK to homogeneous initial guess!"<<std::endl;
            std::cout<<"MISSING FILE for using existing Y_old as initial guess!"<<std::endl<<"FALL BACK to homogeneous initial guess!"<<std::endl;
          }


          // Initial guess for K:
          Y_old.resize(0);
          Y_old.resize( nAllCells, 0.0 );

          for(UINT jj=0; jj<nzones; jj++)
            for(UINT ii=0; ii<nAllCells; ii++)
              Y_old[ii]+=Xzones[jj][ii]*beta_old[jj];

          // write the inital guess to hdf5 (as Y_old)
          HDF5Tools::h5g_Write(
                               Y_old
                               , dir.Y_old_h5file
                               , "/Y_old"
                               , inputdata
                               );
          /*
            Vector<REAL> Y_old_Well (nAllCells,0.0);
            for(int ii=0;ii<nAllCells;ii++)
            Y_old_Well[ii] = Y_old[ii];

            General::add_wells_to_field( inputdata, Y_old_Well );

            HDF5Tools::
            h5g_Write(
            Y_old_Well
            , "/Y_old_Well"
            , inputdata
            , dir.Y_old_Well_h5file
            , "Y_old_Well"
            , dir.Y_old_xmffile
            );
          */

        } // end: else: if ( inputdata.problem_types.using_existing_Yold ...
      }// end:  if(helper.rank()==0)

      // be sure Y_old (as initial guess) is written.
      if(helper.size()>1)
        MPI_Barrier(helper.getCommunicator());


      if( inputdata.problem_types.refine_estimate ){

        if( General::bFileExists( dir.Y_old_h5file ) ){
          VTKPlot::
            plotDataToRefinedGrid< GRID,
                                   GV_GW,
                                   IDT,
                                   DIR,
                                   YFG
                                   >( theGrid,
                                      gv_0,
                                      dir.Y_old_h5file,
                                      dir.Y_old2_h5file,
                                      "/Y_old",
                                      inputdata,
                                      dir );
        }

        if( General::bFileExists( dir.estimation_variance_h5file ) ){
          VTKPlot::
            plotDataToRefinedGrid< GRID,
                                   GV_GW,
                                   IDT,
                                   DIR,
                                   YFG
                                   >( theGrid,
                                      gv_0,
                                      dir.estimation_variance_h5file,
                                      dir.V_est2_h5file,
                                      "/sigma2",
                                      inputdata,
                                      dir );
        }

        if( inputdata.inversion_parameters.max_iter == 0 ){
          if( helper.rank()==0 )
            std::cout << "Y_estimate and Estimation-Variance stored on refined grid. Stop here." << std::endl;
          return -1;
        }

      }

      // A container listing the history of the objective function value
      typedef std::pair<REAL,REAL> REAL_PAIR;
      std::vector< std::vector<REAL_PAIR> > L_history;

      /*
       * calculate the value of the objective function for the initial guess.
       */

      if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
        std::cout << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++ "<< std::endl;
        std::cout << General::getDateAndTime()
                  << "   INVERSION starts" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++ "<< std::endl;
        std::cout << std::endl;
        std::cout << "Calculate the value of the objective function for the initial guess." << std::endl;
        std::cout << std::endl;
      }

      REAL L_m = 0.0;
      REAL L_p = 0.0;
      objectivefunction
        <GRID,
         GFS_GW,
         GFS_TP,
         GFS_CG,
         DIR,
         IDT,
         YFG,
         MEASLIST,
         SensCl
         >
        ( theGrid,
          gfs_gw,
          dir,
          inputdata,
          orig_measurements,
          sim_measurements,
          Vector<REAL>(0),
          1, // 1=Y_old
          log_electricalConductivity,
          log_kappafield,
          sensitivity,
          helper,
          0, // it_counter
          0,  // wCounter
          L_m,
          L_p
          );

      REAL L_old = L_m;

      sim_measurements.write_to_logger("initial guess: simulated measurements");

      // L_old plus the prior term! Only for the initial guess. In all other cases it will be done in the objective function!
      if(helper.rank()==0){ // on P0
        //add the initial prior term (non zero only if use_existing_Y_old is used)!
        L_old += L_prior_ini;
        logger<<"objective function for the initial guess = "<<L_old<<" ( L_prior_ini = "<<L_prior_ini<<" )"<<std::endl;
        //printf( "\n\nobjective function for the initial guess = %10.5f \n\n", L_old );
        L_objective.push_back(L_old);
        L_history.push_back( std::vector<REAL_PAIR>(1,std::make_pair(0,L_old)) );

        logger << "L_objective (init) = " << L_old << std::endl;
        std::cout << "L_objective (init) = " << L_old << std::endl;

      }// end:  if(helper.rank()==0)



      if( inputdata.inversion_parameters.max_iter == 0 ){
        return L_old;
      }

      //----------------------------------------------------------
      // iterative updating
      //----------------------------------------------------------
      unsigned int it_counter = 0;   // a counter for iteration number
      int loopflag = 1;     // flag for performing updates
      int acceptflag = 0;    // flag for acceptable estimates, which are close to optimal
      REAL delta_L = 1E+8;
      std::vector< std::string > vtu_list_Y_old;

      logger << " ..initilization done! (in " << watch.elapsed() << " sec )" << std::endl;

      // L_p corresponding old iteration:
      Vector<REAL> L_p_OLD;
      L_p_OLD.push_back( 0.0 );

      //DenseMatrix<REAL> localGyy(nMeas,nMeas,0.0);
      //Vector<DenseMatrix<REAL> > Gyy;
      //Gyy.push_back( localGyy );



      Vector<REAL> smoothed_orig_YField;
      HDF5Tools::h5g_pRead<GV_GW,Dune::Interior_Partition>( gv_0
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
      /*
      // DEBUG CODE:
      YFG yfg_smoothed_Yorig(inputdata,dir,helper.getCommunicator());
      if( helper.size() > 1 )
      yfg_smoothed_Yorig.parallel_import_from_local_vector( smoothed_orig_YField, local_count, local_offset );
      else
      yfg_smoothed_Yorig.import_from_vector( smoothed_orig_YField );

      yfg_smoothed_Yorig.plot2vtu( gv_0,
      dir.Y_orig_vtu + "_smoothed_inv",
      "Y_smoothed",
      baselevel );
      */

      /*==================
       * Inversion scheme
       *==================
       */
      // <BIG LOOP>: START
      while( loopflag ){


        bool bAllowRemovingJQ_h5file = true;
        //increase iteration counter
        it_counter++;

        //flag which indicates if the SV is increased for a calculation (only needed in M0 inversion)
        UINT SV_increased=0;

        if(helper.rank()==0){
          std::cout << std::endl;
          std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++ "<< std::endl;
          std::cout << General::getDateAndTime()
                    << "ITERATION NUMBER:  "<< it_counter<< std::endl;
          std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++ "<< std::endl;
        }

        logger << std::endl;
        logger << "++++++++++++++++++++++++++++++++++++++++++++++++++++ "<< std::endl;
        logger << General::getDateAndTime()
               << "ITERATION NUMBER:  "<< it_counter<< std::endl;
        logger << "++++++++++++++++++++++++++++++++++++++++++++++++++++ "<< std::endl;

        /*
         * update the sensitivity class
         * important for the calculation of the sensitivities
         */
        sensitivity.set_OLD( it_counter );

        /*
         * calculate sensitivities!!!
         */
        logger<<"Calculate sensitivities ..."<<std::endl;
        if( helper.rank()==0 && inputdata.verbosity >= VERBOSITY_INVERSION ){
          std::cout<<"Calculate sensitivities ..."<<std::endl;
        }
        watch.reset();


#ifdef DEBUG_PLOT
        std::stringstream vtu_Y_old;
        vtu_Y_old << dir.vtudir
                  << "/Y_old"
                  << "_i" << it_counter;



        YFG yfg_Y_old( inputdata, dir, helper.getCommunicator() );
        VTKPlot::output_hdf5data_to_gridview( gv_0,
                                              inputdata,
                                              vtu_Y_old.str(),
                                              dir.Y_old_h5file,
                                              "/Y_old",
                                              yfg_Y_old
                                              );
#endif

        sensitivity.calculate( it_counter,
                               Y_old,
                               nAllCells,
                               sim_measurements,
                               SV_increased,
                               CR );
        logger<<" ...calculation of sensitivities done! (in "<<watch.elapsed()<<" sec )"<<std::endl;


        if( inputdata.plot_options.vtk_plot_sensitivities ){
          sensitivity.vtu_sensitivity( gv_0, it_counter );
        }

        if( inputdata.plot_options.vtk_plot_cross_covariances ){
          sensitivity.vtu_JQ( gv_0, it_counter );
        }

        /*
         * calculate JQJ
         */
        sensitivity.calculate_JQJ();


        /*
         * Clean up hdf5 files that are not used anymore:
         */
#ifdef CLEAN_S
        // Do the file deletion on process 0.
        if( helper.rank()==0 ){
          sensitivity.cleanupFiles();
        }
        // Wait until files are deleted by process 0.
        if( helper.size()>1 )
          MPI_Barrier(helper.getCommunicator());
        // Report to standard output.
        if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
          std::cout << "inversionkernel: All sensitivity h5-files of the current iteration removed from the buffer directory. Next iteration step will produce new files!" << std::endl;
          std::cout << std::endl;
        }
#endif

        // some variables for the weighting loop!
        REAL weighting(1.0);
        Vector<REAL> Y_new;
        Vector<REAL> ksi_new(nMeas);
        Vector<REAL> beta_new(nzones);
        REAL stepsize;

        /*
         * sequential part!
         * for getting the new guess!
         */
        // 'Stick together' the blocks to form the cokriging matrix M:
        DenseMatrix<REAL> CokrigingMatrixM( nMeas + nzones, nMeas + nzones, 0.0 );
        Vector<REAL> xNew( nMeas + nzones, 0.0 );

        // seqential code:
        if(helper.rank()==0){

          watch.reset();

          sensitivity.get_JQJ( CokrigingMatrixM );

          for(UINT iPoint=0; iPoint < nMeas; iPoint++){
            CokrigingMatrixM( iPoint, iPoint ) += MeasurementErrorsV[iPoint];
          }

          // Last lines and last columns of M
          for(UINT jPoint=0; jPoint < nMeas; jPoint++){
            for(UINT ii=0;ii<nzones; ii++){
              CokrigingMatrixM( nMeas+ii , jPoint ) = sensitivity.get_JX(jPoint,ii);
              CokrigingMatrixM( jPoint, nMeas+ii )  = sensitivity.get_JX(jPoint,ii);
            }
          }

          // lower-left block of M:
          for(UINT ii=0;ii<nzones; ii++){
            CokrigingMatrixM( nMeas+ii, nMeas+ii ) = -1.0 / R_bb[ii];
          }

          Vector<REAL> RHS( nMeas + nzones, 0.0 );
          for( UINT iPoint=0; iPoint<nMeas; iPoint++){
            RHS[ iPoint ] = orig_measurements[iPoint].value - sim_measurements[ iPoint ].value;
            RHS[ iPoint ] += sensitivity.get_J_times_Y_old(iPoint);
          }

          // Last elements of the RHS:
          for(UINT ii=0;ii<nzones; ii++){
            RHS[ nMeas+ii ] = -1.0 / R_bb[ii] * beta_star[ii];
          }

          // Now it's time to solve this equation to get the next iteration
          logger << "=== Solve cokriging system" << std::endl;
          CokrigingMatrixM.row_equilibration(); // This function must not be called twice!
          CokrigingMatrixM.gauss_solver( xNew, RHS );
          std::cout << "Solving Cokriging System on root in iteration #" << it_counter
                    << " DONE! ( build up and solving took " << watch.elapsed() << " sec )."
                    << std::endl;

          // Check if the result is good enough:
          Vector<REAL> residual( nMeas + nzones, 0.0 );
          residual = CokrigingMatrixM * xNew;
          residual -= RHS;
          std::cout << "||Ax-b|| = " << residual.two_norm()
                    << std::endl;

          if(inputdata.verbosity>18){
            // debug log for the inversion kernel!
            logger<<"Cokriging Matrix: "<<std::endl;
            if( nMeas+nzones > 20 ){
              logger<<"Larger than 20x20. Do not print."<<std::endl;
            }
            else {
              for(UINT ii=0; ii<nMeas+nzones; ii++){
                for(UINT jj=0; jj<nMeas+nzones; jj++){
                  char buffer[128];
                  sprintf(buffer, "%2.4e ", CokrigingMatrixM(ii, jj));
                  logger<<buffer;
                }
                logger<<std::endl;
              }
            }


            logger<<" RHS: ";
            if( RHS.size() > 20 ) {
              logger<<"Larger than 20x20. Do not print."<<std::endl;
            }
            else{
              for(UINT ii=0; ii<RHS.size(); ii++)
                logger<<RHS[ii]<<" ";
              logger<<std::endl;
            }

            logger<<" Solution Cokriging: ";
            if( xNew.size() > 20 ) {
              logger<<"Larger than 20x20. Do not print."<<std::endl;
            }
            else {
              for(UINT ii=0; ii<xNew.size(); ii++)
                logger<<xNew[ii]<<" ";
              logger<<std::endl;
            }
          }

          watch.reset();
          // extract ksi out of RHS:
          for( UINT iPoint = 0; iPoint<nMeas; iPoint++ ){
            ksi_new[ iPoint ] = xNew[ iPoint ];
          }
        }

#ifdef PARALLEL_KSI_JQ
        // Broadcast ksi_new from 0 to all other processors using MPI_BCast
        MPI_Bcast( &(ksi_new[0]),
                   nMeas,
                   MPI_DOUBLE,
                   0,
                   gv_0.comm() );
#endif

        // sequential code:
        if(helper.rank()==0){

          // extract beta_new out of RHS:
          for( UINT iZone=0; iZone < nzones; ++iZone )
            beta_new[ iZone ] = xNew[ nMeas + iZone ];

          Y_new.resize(0);
          Y_new.resize( nAllCells, 0.0 );
          for(UINT iZone=0; iZone < nzones; ++iZone )
            for(UINT iCell=0; iCell<nAllCells; ++iCell)
              Y_new[ iCell ] = Xzones[ iZone ][ iCell ] * beta_new[ iZone ];

        }

#ifdef PARALLEL_KSI_JQ
        // parallel code:
        Vector<REAL> ksi_JQ(nAllCells);
        Vector<REAL> ksi_JQ_sum(nAllCells);
        for( UINT iPoint=0; iPoint<nMeas; ++iPoint ){

          Vector<REAL> JQ;
          UINT mGroup = iPoint % helper.size();
          if( mGroup == helper.rank() ){
            //std::cout << "Process " << helper.rank()
            //          << " reading JQ " << iPoint
            //          << std::endl;
            HDF5Tools::h5g_Read( JQ
                                 , dir.JQ_h5file[iPoint]
                                 , "/JQ"
                                 , local_offset
                                 , local_count
                                 , inputdata
                                 );
            JQ *= ksi_new[ iPoint ]; // multipliation by a scalar
            ksi_JQ += JQ;
          }
        }

        MPI_Reduce( &(ksi_JQ[0]),
                    &(ksi_JQ_sum[0]),
                    nAllCells,
                    MPI_DOUBLE,
                    MPI_SUM,
                    0,
                    gv_0.comm() );

#endif

        // sequential code:
        if(helper.rank()==0){

#ifdef PARALLEL_KSI_JQ
          for( UINT iCell=0; iCell<nAllCells; iCell++) {
            Y_new[ iCell ] += ksi_JQ_sum[ iCell ];
          }
#else
          // Calculate Y_new = transpose(JQ) * ksi // could be done in parallel! but it seems also sequential very fast!
          for(UINT iPoint=0; iPoint<nMeas; iPoint++  ){
            Vector<REAL> JQ;
            HDF5Tools::h5g_Read( JQ
                                 , dir.JQ_h5file[iPoint]
                                 , "/JQ"
                                 , local_offset
                                 , local_count
                                 , inputdata
                                 );
            for( UINT iCell=0; iCell<nAllCells; iCell++) {
              Y_new[ iCell ] += JQ[ iCell ] * ksi_new[ iPoint ];
            }

            // Clean up JQ h5 files if this is not the very last iteration step!?
          }

#endif


          if(CR>-1){
            for( UINT iCell=0; iCell<nAllCells; iCell++) {
              Y_new[ iCell ] += Y_u[ iCell ];
            }
          }


          logger<<"generating new solution took "<<watch.elapsed()<<" sec."<<std::endl;

          // Compute the current step size
          stepsize = 0.0;
          for( unsigned int iCell=0; iCell<nAllCells; iCell++){
            stepsize = std::max( stepsize, std::abs( Y_new[ iCell ] - Y_old[ iCell ] ) );
          }

          // If the step size is too big, do a line search
          REAL stepsize_limit = inputdata.inversion_parameters.lim_stepsize;
          logger<<"stepsize( "<<it_counter<<" ) = "<<stepsize<<std::endl;
          logger<<"stepsize_limit = "<<stepsize_limit<<std::endl;

          if( stepsize > stepsize_limit ) {
            weighting = stepsize_limit / stepsize;
          }


          //send around the weighting variable -> needed to do the loop in parallel
          for(int ii=1;ii<helper.size();ii++)
            MPI_Send(&weighting,1, MPI_DOUBLE,ii,ii,helper.getCommunicator());

        } // END: if(helper.rank()==0)

        // receive the weighting on all processes except P0
        if(helper.rank()>0){
          MPI_Recv(&weighting,1, MPI_DOUBLE,0,helper.rank(),helper.getCommunicator(),MPI_STATUS_IGNORE);
        }



        //some needed variables in the weighting loop
        UINT wCounter = 0;
        const REAL weighting_limit = inputdata.inversion_parameters.weighting_lim;
        REAL L_try;
        Vector<REAL> ksi_try;
        Vector<REAL> Y_try;
        Vector<REAL> beta_try(nzones);
        int weighting_loop=1;
        logger<<"Weighting loop ..."<<std::endl;
        watch.reset();


        REAL L_p_NEW = 0.0; // variable to store L_p from latest iteration

        /*
         * weighting loop
         */
        std::vector<REAL_PAIR> L_weighting_history;
        while( weighting > weighting_limit && weighting_loop){

          wCounter++;

          if(helper.rank()==0){ // something to do on P0

            //std::cout << "DEBUG: nAllCells = " << nAllCells << std::endl;

            // intermediate solution of the parameters
            Y_try.resize( nAllCells, 0.0 );
            for( UINT iCell=0; iCell<nAllCells; iCell++){
              Y_try[ iCell ] = weighting * Y_new[ iCell ] + (1.0 - weighting) * Y_old[ iCell ];
            }

            ksi_try.resize( nMeas, 0.0 ) ;
            for( UINT iPoint=0; iPoint < nMeas; iPoint++ )
              ksi_try[ iPoint ] = weighting * ksi_new[ iPoint ] + (1.0 - weighting) * ksi_old[ iPoint ];

            //std::cout << "DEBUG: ksi_try = " << std::endl;
            //std::cout << ksi_try << std::endl;

            for( UINT iZone=0; iZone < nzones; ++iZone )
              beta_try[ iZone ] = weighting * beta_new[ iZone ] + (1.0 - weighting) * beta_old[ iZone ];

            //std::cout << "DEBUG: beta_try = " << std::endl;
            //std::cout << beta_try << std::endl;

            //write SEQUENTIAL the new trial solution
            HDF5Tools::h5g_Write(
                                 Y_try
                                 , dir.Y_try_h5file
                                 , "/Y_try"
                                 , inputdata
                                 );

            logger<<"Evaluate the value of the objective function based on the intermediate solution Y_try"<<std::endl;


          } // END: if(helper.rank()==0)

          // be sure Y_try is written.
          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());
          /*
           * calculate the value of the objective function for the initial guess.
           * parallel!!!
           */
          MEASLIST measurements_try(inputdata,true);

          if(helper.rank()==0){
            std::cout << std::endl
                      << General::getDateAndTime()
                      << "Compute objective function for the trial field..." << std::endl;
          }


          REAL L_m_try = 0.0;
          REAL L_p_try = 0.0;
          objectivefunction
            <GRID,
             GFS_GW,
             GFS_TP,
             GFS_CG,
             DIR,
             IDT,
             YFG,
             MEASLIST,
             SensCl
             >
            ( theGrid,
              gfs_gw,
              dir,
              inputdata,
              orig_measurements,
              measurements_try,
              ksi_try,
              0,  // <--- 0 means: we use Y_try
              log_electricalConductivity,
              log_kappafield,
              sensitivity,
              helper,
              it_counter,
              wCounter,
              L_m_try,
              L_p_NEW );


          //REAL ww = std::pow(0.5,REAL(wCounter-1));
          //L_p_try = L_p_OLD[it_counter-1]
          //  + ww*ww * (L_p_NEW - L_p_OLD[it_counter-1]);
          L_p_try = L_p_NEW;

          if( inputdata.inversion_parameters.L_prior == "on" ) {

            L_try = L_m_try + L_p_try;

          }
          else {

            if(helper.rank()==0) {
              std::cout << "L_prior switched off!" << std::endl;
            }

            L_try = L_m_try;

          }

          // sequential again
          if(helper.rank()==0){  // something to do for P0
            logger<<"objective function for the trial (global, over all processes) = "<<L_try<<std::endl;


            std::cout << General::getDateAndTime()
                      << "L_objective: weighting = " << weighting
                      << " it_counter = " << it_counter
                      << " wCounter = " << wCounter
                      << std::endl;
            std::cout << General::getDateAndTime()
                      << "   L_p_OLD = " << L_p_OLD[it_counter-1]
                      << "   L_old   = " << L_old
                      << std::endl;
            std::cout << General::getDateAndTime()
                      << "   L_m_try = " << L_m_try
                      << "   L_p_try = " << L_p_try
                      << "   L_try   = " << L_try
                      << std::endl;

            logger << General::getDateAndTime()
                   << "L_objective: weighting = " << weighting
                   << std::endl;
            logger << General::getDateAndTime()
                   << "   L_p_OLD = " << L_p_OLD[it_counter-1]
                   << "   L_old   = " << L_old
                   << std::endl;
            logger << General::getDateAndTime()
                   << "   L_m_try = " << L_m_try
                   << "   L_p_try = " << L_p_try
                   << "   L_try   = " << L_try
                   << std::endl;


            delta_L = std::abs( L_old - L_try );

            L_weighting_history.push_back(std::make_pair(weighting,L_try));
            /*
             * updating old values for the next iteration:
             */
            if( L_try < L_old || SV_increased){

              if(L_try<L_old){
                std::cout << "New iteration is an improvement!!!" << std::endl;
                std::cout << "Updating ..." << std::endl;
                logger<<"New iteration is an improvement!!!"<<std::endl;
                logger<<"Updating ..."<<std::endl;
              }


              L_old = L_try;

              L_p_OLD.push_back( L_p_try );
              sensitivity.store_Gyy();

              L_objective.push_back(L_old);
              L_history.push_back( L_weighting_history );

              for( UINT iPoint=0; iPoint < nMeas; ++iPoint ){
                ksi_old[ iPoint ] = ksi_try[ iPoint ];
              }

              for( UINT iCell=0; iCell < nAllCells; ++iCell )
                Y_old[ iCell ] = Y_try[ iCell ];

              for( UINT iZone=0; iZone < nzones; ++iZone )
                beta_old[ iZone ] = beta_try[ iZone ];

              HDF5Tools::h5g_Write(
                                   Y_old
                                   , dir.Y_old_h5file
                                   , "/Y_old"
                                   , inputdata
                                   );

              std::ostringstream cmdline;
              cmdline << "cp " << dir.Y_old_h5file
                      << " " << dir.bufferdimdir
                      << "/Y_old_" << it_counter
                      << ".h5";

              General::systemCall( cmdline.str() );

              std::cout << General::getDateAndTime()
                        << "Backup copy of estimated field: Y_old_"
                        << it_counter << ".h5"
                        << std::endl;


              /*
                Vector<REAL> Y_old_Well( Y_old );
                General::add_wells_to_field( inputdata, Y_old_Well );

                HDF5Tools::
                h5g_Write(
                Y_old_Well
                , "/Y_old_Well"
                , inputdata
                , dir.Y_old_Well_h5file
                , "Y_old_Well"
                , dir.Y_old_xmffile
                );
              */

              Vector<UINT> dimensions(1,ksi_old.size());

              //std::cout << "DEBUG: ksi_old write to file = " << std::endl;
              //std::cout << ksi_old << std::endl;

              HDF5Tools::h5_Write( ksi_old,
                                   dir.ksi_old_h5file,
                                   "/ksi_old",
                                   dimensions
                                   );

              dimensions[0]=nzones;

              //std::cout << "DEBUG: beta_old write file = " << std::endl;
              //std::cout << beta_old << std::endl;
              HDF5Tools::h5_Write( beta_old,
                                   dir.beta_old_h5file,
                                   "/beta_old",
                                   dimensions );


              if( L_try < L_accept_value ){
                std::cout << "Trial solution is close to optimal and it is acceptable!!!" << std::endl;
                logger << "Trial solution is close to optimal and it is acceptable!!!" << std::endl;
                acceptflag = 1;
              }

              //break doesn't work in parallel -> because break is then only done on P0
              weighting_loop=0;
            } else {
              weighting = weighting / 2.0;
              logger<<"No improvement of L, therefore try it again with a reduced weighting = "<<weighting<<std::endl;
            }

            // send the needed information of the weighting loop to P1,...
            for(int ii=1;ii<helper.size();ii++){
              MPI_Send(&weighting,1, MPI_DOUBLE,ii,ii,helper.getCommunicator());
              MPI_Send(&weighting_loop,1, MPI_INT,ii,ii,helper.getCommunicator());

            }
          } // END: if(helper.rank()==0)

          // be sure Y_old is there
          if(helper.size()>1)
            MPI_Barrier(helper.getCommunicator());


          if( inputdata.plot_options.vtk_plot_y_old ){
            if(CR==-1){
              std::stringstream vtu_file;
              vtu_file << dir.Y_old_vtu_prefix << it_counter;
              YFG yfg_Y_old( inputdata, dir, helper.getCommunicator() );
              VTKPlot::output_hdf5data_to_gridview( gv_0,
                                                    inputdata,
                                                    vtu_file.str(),
                                                    dir.Y_old_h5file,
                                                    "/Y_old",
                                                    yfg_Y_old
                                                    );  // function defined in this file
              vtu_list_Y_old.push_back( vtu_file.str() );
            }
          }

          // receive the needed loop-data (weighting loop)
          if(helper.rank()>0){
            MPI_Recv(&weighting,1, MPI_DOUBLE,0,helper.rank(),helper.getCommunicator(),MPI_STATUS_IGNORE);
            MPI_Recv(&weighting_loop,1, MPI_INT,0,helper.rank(),helper.getCommunicator(),MPI_STATUS_IGNORE);
          }
          //set the measurements correct
          if(weighting_loop==0){
            if(helper.rank()==0)
              std::cout << "Re-assign measurements..." << std::endl;
            sim_measurements = measurements_try;
            sim_measurements.write_to_logger("Y_try: simulated measurements");
          }


          //MPI_Bcast(&(weighting_loop[0]),1,MPI_INT,0,helper.getCommunicator()); DOES NOT WORK WHY!???
        } // end of while( weighting > weighting_limit )
        /*
         * DONE: weighting loop
         */

        logger<<"...weighting loop DONE! ( took "<<watch.elapsed()<<" sec )"<<std::endl;

        // now the weighting loop is done and it will be calculated sequential again
        if(helper.rank()==0){ // something to do for P0


          likelihood_ostream<<it_counter<<":"<<L_try<<std::endl;


          logger <<std::endl<< "Iteration "<< it_counter<<" done: "<<std::endl;
          logger <<"L_try = "<< L_try<<std::endl;
          logger<< "delta_L = "<< delta_L<<std::endl;

          // Stopping criterions:
          // decision making about the iterative estimation procedure
          if( !SV_increased
              &&
              (
               stepsize < inputdata.inversion_parameters.dY_lim
               ||
               delta_L < inputdata.inversion_parameters.dL_lim_fac * nMeas
               ||
               (weighting <= weighting_limit && acceptflag == 1)
               ) ){

            std::cout << std::endl << "Iteration #"<< it_counter <<" done: "
                      << std::endl;
            std::cout << "CONVERGED !!! " << std::endl;

            // && L_try<1.05*nmeas
            if( stepsize < inputdata.inversion_parameters.dY_lim ){
              std::cout << "Convergence Criterion 1 reached: "
                        << "Max. change in the parameter field is "
                        << stepsize
                        << ". This is below the predefined absolute limit "
                        << inputdata.inversion_parameters.dY_lim
                        << std::endl;
            }

            if( delta_L < inputdata.inversion_parameters.dL_lim_fac * nMeas ){
              std::cout << "Convergence Criterion 2 reached: "
                        << "Last change in L_objective is "
                        << delta_L
                        << ". This is below the predefined limit "
                        << inputdata.inversion_parameters.dL_lim_fac * 100.0
                        << " percent of " << nMeas << " measurement values."
                        << std::endl;
            }

            std::cout << "L_try = " << L_try << std::endl;

            if( acceptflag == 1 ) {
              std::cout << "The solution is within the confidence interval."
                        << std::endl;
            }
            else {
              std::cout << "The solution is not within the confidence interval."
                        << std::endl;
            }

            logger<<"CONVERGED !!!"<<std::endl;
            loopflag = 0;

          }else if( it_counter == inputdata.inversion_parameters.max_iter ){

            std::cout << "L_try = " << L_try << std::endl;
            std::cout<<"QUITTING ALGORITHM: MAXIMUM ITERATION NUMBER REACHED"<<std::endl;
            logger<<"QUITTING ALGORITHM: MAXIMUM ITERATION NUMBER REACHED"<<std::endl;
            if(SV_increased){
              std::cout<<"SV is increasedon m0 measurements"<<std::endl;
              logger<<"SV is increasedon m0 measurements"<<std::endl;
            }
            loopflag = 0;

          }else if( weighting <= weighting_limit && acceptflag == 0 && !SV_increased){

            std::cout << "L_try = " << L_try << std::endl;
            std::cout << "weighting = " << weighting << std::endl;
            std::cout<<"QUITTING ALGORITHM: GOT STUCK"<<std::endl;

            logger << "L_try = " << L_try << std::endl;
            logger << "weighting = " << weighting << std::endl;
            logger<<"QUITTING ALGORITHM: GOT STUCK"<<std::endl;
            loopflag = 0;
          }

          // send the loop flag to all processors
          for(int ii=1;ii<helper.size();ii++){
            MPI_Send(&loopflag,1, MPI_INT,ii,ii,helper.getCommunicator());
          }

          bAllowRemovingJQ_h5file = loopflag;

#ifdef ESTIMATION_VARIANCE
          if(!loopflag && CR==-1){
            DenseMatrix<REAL> inv_CokrigingMatrixM( CokrigingMatrixM.n_rows(), CokrigingMatrixM.n_cols(), 0.0 );

            watch.reset();
            CokrigingMatrixM.inverse( inv_CokrigingMatrixM ); // Calculation of inverse might take very long!
            std::cout << "Computing inverse of cokriging matrix on root took: "
                      << watch.elapsed() << " sec !"
                      << std::endl;

            inv_CokrigingMatrixM.write_to_HDF5( dir.cokriging_inverse_h5file, "/Minv" );
          }
#endif

        } // END: if(helper.rank()==0)

        // be sure inverse cokriging is written.
        if(helper.size()>1)
          MPI_Barrier(helper.getCommunicator());

        //receive the loopflag
        if(helper.rank()>0){
          MPI_Recv(&loopflag,1, MPI_INT,0,helper.rank(),helper.getCommunicator(),MPI_STATUS_IGNORE);
        }
        //logger<<"loopflag = "<<loopflag<<std::endl;


#ifdef CLEAN_JQ
        // Clean up JQ h5 files here unless we are finished with the inversion scheme!
        if( bAllowRemovingJQ_h5file ){
          // Do the file deletion on the root process.
          if( helper.rank()==0 ){
            for( UINT iPoint=0; iPoint < nMeas; iPoint++ ){
              General::deleteFile( dir.JQ_h5file[iPoint] );
            }
          }
          // Report to standard output.
          if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
            std::cout << "inversionkernel: All JQ h5-files of the current iteration are removed from the buffer directory. Next iteration step will produce new files!" << std::endl;
          }
        }
        else{
          std::cout << "inversionkernel: Keeping very last set of JQ h5-files for calculation of estimation variance." << std::endl;
        }

        // Wait until files are deleted by process 0.
        if( helper.size()>1 )
          MPI_Barrier(helper.getCommunicator());

#endif


        // Compute quality of the solution by comparing the
        // structure of the estimated field Y_old
        // with a smoothed version of the original field.
        // The smoothing is achieved by a scaled multiplication
        // alpha * C_YY * Y_orig.
        // The quality measure is then defined as
        // || Y_old - alpha * C_YY * Y_orig || in the 2-norm.
        //
        // Read /Y_old in parallel on gv_0 and store it to Y_current.

        Vector<REAL> Y_current;
        HDF5Tools::h5g_pRead<GV_GW,Dune::Interior_Partition>( gv_0
                                                              , Y_current
                                                              , dir.Y_old_h5file
                                                              , "/Y_old"
                                                              , local_offset
                                                              , local_count
                                                              , inputdata
                                                              , 1 // P0 blocksize
                                                              , FEMType::DG // P0
                                                              , 0 // structure is on grid level 0
                                                              );

        if( smoothed_orig_YField.size() == Y_current.size() ) {

          REAL Quali = ( smoothed_orig_YField - Y_current ).two_norm_sqr();
          Quali = gv_0.comm().sum( Quali );

          if( helper.rank() == 0 ) {
            std::cout << "Quality(Y_est) = "
                      << Quali
                      << std::endl;
          }

        }
        else {
          std::cout << "WARNING: Cannot compute Quality(Y_est) due to different size: "
                    << "smoothed_orig_YField.size() = " << smoothed_orig_YField.size()
                    << "Y_current.size() = " << Y_current.size()
                    << std::endl;
        }


      }// END:  while( loopflag )
      // <BIG LOOP>: END


      if( inputdata.plot_options.vtk_plot_y_old ){
        if(CR==-1)
          General::createPVDfromVtuList( dir.Y_old_pvd, vtu_list_Y_old );
      }


      /*
       * save the estimates as vtu files
       */

      //save the estimated Y field
      if(CR>-1){
        if(helper.rank()==0){
          HDF5Tools::h5g_Write( Y_old
                                , dir.Y_cond_h5file[CR]
                                , "/Y_cond"
                                , inputdata
                                );
        }

        YFG yfg_Y_old( inputdata, dir, helper.getCommunicator() );
        VTKPlot::output_hdf5data_to_gridview( gv_0,
                                              inputdata,
                                              dir.Y_cond_vtu[CR],
                                              dir.Y_old_h5file,
                                              "/Y_old",
                                              yfg_Y_old
                                              );

      }else{

        if(helper.rank()==0){
          HDF5Tools::h5g_Write( Y_old
                                , dir.Y_estimated_h5file
                                , "/Y_est"
                                , inputdata
                                );
        }

        if( inputdata.problem_types.refine_estimate ){
          VTKPlot::
            plotDataToRefinedGrid< GRID,
                                   GV_GW,
                                   IDT,
                                   DIR,
                                   YFG
                                   >( theGrid,
                                      gv_0,
                                      dir.Y_estimated_h5file,
                                      dir.Y_estimated2_h5file,
                                      "/Y_est",
                                      inputdata,
                                      dir );
        }

      }



      if(helper.rank()==0){

        std::cout << "total measuring points : "<<nMeas<<std::endl;
        for( UINT iPoint=0; iPoint < nMeas; iPoint++ ){

          REAL discrepancy = std::abs( orig_measurements[iPoint].value - sim_measurements[iPoint].value );
          REAL relativerfehler = std::abs( ( orig_measurements[iPoint].value - sim_measurements[iPoint].value ) / orig_measurements[iPoint].value );

          measurement_ostream << iPoint << ":"
                              << orig_measurements[iPoint].value<<":"
                              << sim_measurements[iPoint].value<<":"
                              << discrepancy << ":"
                              << relativerfehler
                              << std::endl;
        }

        likelihood_ostream.close();
        measurement_ostream.close();
      }


#ifdef ESTIMATION_VARIANCE
      if(CR==-1){
        logger<<"calculate estimation variance (parallel)..."<<std::endl;
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



        // Clean up JQ h5 files finally.
        // Do the file deletion on the root process.
        if( helper.rank()==0 ){
          for( UINT iPoint=0; iPoint < nMeas; iPoint++ ){
            General::deleteFile( dir.JQ_h5file[iPoint] );
          }
        }
        // Report to standard output.
        if( helper.rank()==0 && inputdata.verbosity>=VERBOSITY_INVERSION ){
          std::cout << "inversionkernel: Estimation variance stored in '"
                    << dir.estimation_variance_h5file
                    << "'. All JQ h5-files of the very last iteration are removed from the buffer directory." << std::endl;
        }


#ifdef USE_HISTO_PLOT

        // Compute NormalizedError[i] := ( Y_orig[i] - Y_est[i] ) / sigma[i]
        Vector<REAL> NormalizedError;
        HDF5Tools::h5g_pRead<GV_GW,Dune::Interior_Partition>( gv_0
                                                              , dir.kfield_h5file
                                                              , "/YField"
                                                              , NormalizedError
                                                              , local_offset
                                                              , local_count
                                                              , inputdata
                                                              , 1 // P0 blocksize
                                                              , FEMType::DG // P0
                                                              , 0 // structure is on grid level 0
                                                              );

        UINT originalFieldSize = NormalizedError.size();

        Vector<REAL> Y_current;
        HDF5Tools::h5g_pRead<GV_GW,Dune::Interior_Partition>( gv_0
                                                              , Y_current
                                                              , dir.Y_old_h5file
                                                              , "/Y_old"
                                                              , local_offset
                                                              , local_count
                                                              , inputdata
                                                              , 1 // P0 blocksize
                                                              , FEMType::DG // P0
                                                              , 0 // structure is on grid level 0
                                                              );

        UINT estimatedFieldSize = Y_current.size();

        if( estimatedFieldSize != originalFieldSize ){
          std::cout << "WARNING: originalFieldSize(" << originalFieldSize << ") and estimatedFieldSize(" << estimatedFieldSize << ") are different. Cannot compute weighted error for histogram." << std::endl;
        }
        else {

          REAL shift = NormalizedError.mean() - Y_current.mean();

          for(UINT i=0;i<NormalizedError.size();++i){
            NormalizedError[i] -= (Y_current[i] - shift); // This shift is necessary to hit the center of the standary deviation!
            NormalizedError[i] /= std::sqrt( EstimatedVariances[i] );
          }

          // plot NormalizedError:
          YFG yfg_NormalizedError( inputdata, dir, helper.getCommunicator() );

          if( helper.size() > 1 )
            yfg_NormalizedError.parallel_import_from_local_vector( NormalizedError, local_offset, local_count );
          else
            yfg_NormalizedError.import_from_vector( NormalizedError );

          yfg_NormalizedError.plot2vtu( gv_0,
                                        dir.vtudir + "/NormalizedError",
                                        "NormalizedError",
                                        baselevel );

          const int nBins = inputdata.inversion_parameters.histogram_bins;
          General::vector2histogram( gv_0, NormalizedError, nBins, dir.vtudir + "/histogram.dat" );
          General::vector2histogram( gv_0, NormalizedError, 2*nBins, dir.vtudir + "/histogram2.dat" );

          std::cout << "Write histogram to " << dir.vtudir << std::endl;

          // to plot the histogram using gnuplot and compare it against the
          // standard deviation, do this:
          //
          // Norm(x,m,s) = 1./(sqrt(2*pi)*s) * exp( -(x-m)**2 / (2*s*s) )
          // plot [-5:5] Norm(x,0,1), "histogram.dat" using 1:2 with steps
          //

        }

#endif //USE_HISTO_PLOT

        // Wait until files are deleted by process 0.
        if( helper.size()>1 )
          MPI_Barrier(helper.getCommunicator());

      }
#endif


      double L_return=0;


      /************************************************************************
       Output history of objective function to the screen and to rank0-logfile
      *************************************************************************/

      if(helper.rank()==0){

        logger<<std::endl<<"Development of objective function:"<<std::endl;
        std::cout<<std::endl<<"Development of objective function:"<<std::endl;

        REAL LStart = L_history[0][0].second;

        UINT iCounter=0;
        for( std::vector<std::vector<REAL_PAIR>>::const_iterator L_it=L_history.begin()
               ; L_it!=L_history.end() ; L_it++ ){

          if(iCounter==0){
            std::cout<<"Initial guess: " << std::endl;
            logger<<"Initial guess: " << std::endl;
          }
          else{
            std::cout<<"it# "<< iCounter <<" : " << std::endl;
            logger<<"it# "<< iCounter <<" : " << std::endl;
          }
          ++iCounter;
          UINT wCounter=0;
          for( std::vector<REAL_PAIR>::const_iterator L_w_it=(*L_it).begin()
                 ; L_w_it!=(*L_it).end()
                 ; L_w_it++ ){

            REAL w_Current = (*L_w_it).first;
            REAL L_current = (*L_w_it).second;
            REAL percentage = L_current/LStart;

            if(iCounter!=0){
              std::cout << "weighting = "
                        << std::fixed << std::setw(6) << std::setprecision(4)
                        << w_Current << ", ";
            }
            std::cout << "L = "
                      << std::scientific << std::setw(6) << std::setprecision(2)
                      << L_current
                      << " ( "
                      << std::scientific << std::setw(6) << std::setprecision(2)
                      << percentage*100 << "% )"
                      << std::endl;
            ++wCounter;

            // Be careful:
            // L_Try might exceed LStart, but the
            // fraction for the plot must never exceed 1
            REAL fraction = std::min( percentage, 1.0 );
            UINT tmp_n = std::floor( std::floor( fraction*50 + 0.5 ) + 0.5 );

            //std::cout << "tmp_n = " << tmp_n << std::endl;

            logger << "w = "
                   << std::fixed << std::setw(6) << std::setprecision(4)
                   << w_Current << ": |";
            for(UINT jj=0; jj<tmp_n; jj++ )
              logger<<"=";
            for(UINT jj=0; jj<50-tmp_n; jj++)
              logger<<" ";
            logger<<"|; L = " << L_current << ", ( "<< percentage*100 <<"% )"
                  << std::endl;
          }

        };

        L_return = L_objective[ L_objective.size() - 1 ];

      }

      orig_measurements.write_to_logger("original measurements");
      sim_measurements.write_to_logger("estimated measurements");

      if(helper.rank()==0){
        logger<<std::endl
              <<"Number of measurements : "<<nMeas
              <<", chi2 confidence interval at : "<< inputdata.inversion_parameters.L_accept_confidence_interval
              <<", chi2inv( "<< inputdata.inversion_parameters.L_accept_confidence_interval
              <<" , "<<nMeas
              <<" ) = "<<L_accept_value
              <<std::endl;
        std::cout<<std::endl
                 <<"Number of measurements : "<<nMeas
                 <<", chi2 confidence interval at : "<< inputdata.inversion_parameters.L_accept_confidence_interval
                 <<", chi2inv( "
                 << inputdata.inversion_parameters.L_accept_confidence_interval
                 <<" , "<<nMeas<<" ) = "<<L_accept_value
                 <<std::endl;
      }

      if( helper.size()>1 ){
        // Broadcast L_return to other processes:
        MPI_Bcast(&(L_return),1,MPI_DOUBLE,0,helper.getCommunicator());
      }

      return L_return;
    }


  }
}

#endif // DUNE_GESIS_INVERSION_KERNEL_HH
