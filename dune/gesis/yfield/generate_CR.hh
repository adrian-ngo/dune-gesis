#ifndef _GENERATE_CR_HH
#define _GENERATE_CR_HH


#include <time.h>                      // define time()

namespace Dune {
  namespace GeoInversion {

      
      
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
void generate_CR( GRID& theGrid,
                  const PGV pRootGridView,
                  const GFS_GW& gfs_gw,
                  const IDT& inputdata,
                  const YFG& yfg_orig,
                  MEASLIST& orig_measurements,
                  const std::vector< Vector<UINT> >& nCellsExt,
                  DIR& dir,
                  const YFG& log_electricalConductivity,
                  const YFG& log_kappafield,
                  const Dune::MPIHelper& helper
                  ){     

    logger<<"Generate "<<inputdata.CR_parameters.total<<" conditional realizations!!!"<<std::endl;
    Dune::Timer watch_total;
    watch_total.reset();
    
    // counter if some realizations need to be done more often! 
    UINT redo=0;

    UINT nmeas=orig_measurements.nMeasurements();
    std::vector< REAL > MeasurementErrorsV(nmeas, 0.0 );

    for(UINT ii=0;ii<nmeas;ii++)
        MeasurementErrorsV[ ii ] = pow(orig_measurements[ii].value * orig_measurements[ ii ].rel_error 
                         + orig_measurements[ ii ].abs_error,2);
        
#ifdef SEED_OFF
    int seed = 100;
    // different seed for different processors -> very IMPORTANT to obtain the right result!
    seed+=helper.rank();
#else
    int seed = (int) (time(0)); // create seed out ot the current time
    // different seed for different processors -> very IMPORTANT to obtain the right result!
    seed+=helper.rank();
#endif
    StochasticLib1 stochastic( seed ); // make instance of random library
    
    YFG YField_unconditional( yfg_orig );
   
    dir.set_CR(inputdata.CR_parameters.total,inputdata.CR_parameters.start);
    
    std::vector<UINT> bad_realizations(0);
    std::vector<UINT> realization_info(inputdata.CR_parameters.total,0);
    std::vector<REAL> L(inputdata.CR_parameters.total,0.0);
    
    // chi 2 test for acceptance
    REAL L_accept_value=boost::math::gamma_p_inv(((double)nmeas )*0.5, inputdata.inversion_parameters.L_accept_confidence_interval )*2.0;
    
    
    for(UINT iCR=inputdata.CR_parameters.start; iCR<inputdata.CR_parameters.start+inputdata.CR_parameters.total; iCR++ ){
        
        Dune::Timer watch;
        watch.reset();
        
        
        /*
         * some Output for the information of the user!
         */
        if(redo){
            logger<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION ( for the " <<redo+1<<" time) #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ... "<<std::endl<<std::endl;
            if(helper.rank()==0)
                std::cout<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION ( for the " <<redo+1<<" time) #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ..."<<std::endl<<std::endl;       
        }else{
            logger<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ... "<<std::endl<<std::endl;
            if(helper.rank()==0)
                std::cout<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ..."<<std::endl<<std::endl;
        }
        
        
        
        /*
         * add measurement error to orig_measurements!!
         * 
         */
        MEASLIST measurements_rand(orig_measurements);
        //measurements_rand.write_to_logger("measurements rand, before randomizing");
        for(UINT ii=0;ii<nmeas;ii++){
            measurements_rand[ii].value+=stochastic.Normal( 0.0, std::sqrt(MeasurementErrorsV[ ii ]) );
        }
        //measurements_rand.write_to_logger("measurements rand, after randomizing");
        
        
        /*
         * generate unconditional realization with zero mean and R_YY
         * 
         */

#ifdef DIMENSION3
        YField_unconditional.generate3d( false, true );
#else
        YField_unconditional.generate2d( false, true  );
#endif
        
        
        
        /*
         * call the inverse kernel for doing the real work!!
         * 
         */
        
        L[iCR-inputdata.CR_parameters.start]
          =inversionkernel<GRID,
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
            ( theGrid,
              pRootGridView,
              gfs_gw,
              inputdata,
              orig_measurements,
              //ACHTUNG ÄNDERUNG ZONES
              nCellsExt,
              dir,
              log_electricalConductivity,
              log_kappafield,
              helper,
              iCR-inputdata.CR_parameters.start
              );

        /*
         * check if the realization is good
         *   (if not set iCR -=1; // means re-do the calculation for this number!)
         *   set redo flag and check
         */
        if(L[iCR-inputdata.CR_parameters.start]>L_accept_value){
            redo+=1;
        }else{
            realization_info[iCR-inputdata.CR_parameters.start]=redo;
            redo=0;
        }
        if(redo>realization_info[iCR-inputdata.CR_parameters.start])
            realization_info[iCR-inputdata.CR_parameters.start]=redo;
           
         /*
         * some Output for the information of the user!
         */      
        if(redo && redo<inputdata.CR_parameters.max_redo){
            logger<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ...NOT SUCESSFUL (will try again)! (took "<<watch.elapsed()<<" sec.)"<<std::endl<<std::endl;
            if(helper.rank()==0)
                std::cout<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ...NOT SUCESSFUL (will try again)! (took "<<watch.elapsed()<<" sec.)"<<std::endl<<std::endl;
        }else if(redo>=inputdata.CR_parameters.max_redo){
            logger<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ...NOT SUCESSFUL. MAX NUMBER  OF TRIES ("<<inputdata.CR_parameters.max_redo<<") IS REACHED --> Next realization! (took "<<watch.elapsed()<<" sec.)"<<std::endl<<std::endl;
            if(helper.rank()==0)
                std::cout<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ...NOT SUCESSFUL. MAX NUMBER  OF TRIES ("<<inputdata.CR_parameters.max_redo<<") IS REACHED --> Next realization! (took "<<watch.elapsed()<<" sec.)"<<std::endl<<std::endl;
            redo=0;
            continue;
        }else{
            logger<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ...DONE! (took "<<watch.elapsed()<<" sec.)"<<std::endl<<std::endl;
            if(helper.rank()==0)
                std::cout<<std::endl<<std::endl<<"GENERATE CONDITIONAL REALIZATION #"<<iCR<<" (this is the "<<iCR-inputdata.CR_parameters.start+1<<" of "<<inputdata.CR_parameters.total<<" realizations) ...DONE! (took "<<watch.elapsed()<<" sec.)"<<std::endl<<std::endl;
        
        }
        if(redo)
            iCR-=1;
        
    } // END : for(UINT iCR=inputdata.CR_parameters.start; iCR<inputdata.CR_parameters.start+inputdata.CR_parameters.total; iCR++ )
    
    logger<<std::endl<<"Conditional realization generation finished. It took "<<watch_total.elapsed()<<" sec. to generate "<<inputdata.CR_parameters.total<<" (numbered from "<<inputdata.CR_parameters.start<<" to "<<inputdata.CR_parameters.total+inputdata.CR_parameters.start-1<<")!"<<std::endl;
    logger<<"Information about the realizations ( L acceptance value = "<<L_accept_value<<" )  :"<<std::endl;
    for(UINT ii=0; ii<inputdata.CR_parameters.total; ii++){
        logger<<inputdata.CR_parameters.start+ii<<": L = "<<L[ii];
        if(realization_info[ii]){
            if(realization_info[ii]>=inputdata.CR_parameters.max_redo){
                logger<<" ( "<<realization_info[ii]<<" tries were performed. Conditional realization is still NOT satisfactory good!)"<<std::endl;
            }else{
                logger<<" ( "<<realization_info[ii]+1<<" tries were necessary to obtain a conditional realization!)"<<std::endl;
            }
        }else{
            logger<<" (First try accepted!)"<<std::endl;
        }
    }
    
    return;
}

  }
}


#endif // _GENERATE_CR_HH
