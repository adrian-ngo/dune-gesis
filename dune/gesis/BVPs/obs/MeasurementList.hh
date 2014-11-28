/*
 * MeasurementList.hh
 * Author: R.L.Schwede and A. Ngo
 */
#ifndef DUNE_GESIS_MEASUREMENTLIST_HH
#define DUNE_GESIS_MEASUREMENTLIST_HH

#include "take_measurements_on_dgf.hh"
#include "take_measurements_on_dgf_lnK.hh"
#include "dune/gesis/common/MPI_Tools.hh"


namespace Dune {

  namespace Gesis {


    template<typename MeasurementElement>
    class MeasurementList{

      const CInputData& inputdata;
      bool synthetic;
      UINT nMeas;
      UINT nSetups;
      std::vector<UINT> meas_per_setup;
      std::vector<std::vector<UINT> > meas_per_setup_per_type;
      std::vector<std::vector<UINT> > GEmeas_per_setup_per_config;

      std::vector<std::vector<std::vector<std::vector<UINT> > > > local_to_global;

      std::vector<std::vector<MeasurementElement> > lnK_meas;
      std::vector<std::vector<MeasurementElement> > head_meas;
      std::vector<std::vector<MeasurementElement> > M0_meas;
      std::vector<std::vector<MeasurementElement> > M1_meas;

      std::vector<std::vector<MeasurementElement> >  heatM0_meas;
      std::vector<std::vector<MeasurementElement> >  heatM1_meas;
      std::vector<std::vector<MeasurementElement> > heat_meas;

      std::vector<std::vector<std::vector<MeasurementElement> > > GEM0_meas;
      std::vector<std::vector<std::vector<MeasurementElement> > > GEM1_meas;
      std::vector<std::vector<std::vector<MeasurementElement> > > GE_meas;

      std::vector<MeasurementElement*> all_meas;



    public:

      //constructor
      MeasurementList<MeasurementElement>( const CInputData& inputdata_,
                                           const bool synthetic_ )
      : inputdata(inputdata_),
        synthetic(synthetic_),
        nMeas(0),
        nSetups( inputdata.setups.size() )
      {
        meas_per_setup.resize(nSetups,0);
        meas_per_setup_per_type.resize(nSetups);
        GEmeas_per_setup_per_config.resize(nSetups);

        local_to_global.resize(nSetups);

        lnK_meas.resize(nSetups);
        head_meas.resize(nSetups);
        M0_meas.resize(nSetups);
        M1_meas.resize(nSetups);

        if(synthetic){
          heatM0_meas.resize(nSetups);
          heatM1_meas.resize(nSetups);
        }
        heat_meas.resize(nSetups);

        if(synthetic){
          GEM0_meas.resize(nSetups);
          GEM1_meas.resize(nSetups);
        }
        GE_meas.resize(nSetups);

        all_meas.resize(0);

        /*
         * Setup the vectors
         */
        for(UINT iSetup=0; iSetup<nSetups; iSetup++){
          meas_per_setup_per_type[iSetup].resize(6,0);
          local_to_global[iSetup].resize(6);

          //lnK meas!
          local_to_global[iSetup][0].resize(1);
          if(inputdata.problem_types.lnK_inversion){
            meas_per_setup[iSetup]+= inputdata.setups[iSetup].lnk_inversion_data.mplist.pointdata_vector.size();
            meas_per_setup_per_type[iSetup][0]=inputdata.setups[iSetup].lnk_inversion_data.mplist.pointdata_vector.size();
            for(UINT ii=0; ii< meas_per_setup_per_type[iSetup][0]; ii++){
              lnK_meas[iSetup].push_back(inputdata.setups[iSetup].lnk_inversion_data.mplist.pointdata_vector[ii]);
            }
          }

          //head meas!
          local_to_global[iSetup][1].resize(1);
          if(inputdata.problem_types.head_inversion){
            meas_per_setup[iSetup]+= inputdata.setups[iSetup].head_inversion_data.mplist.pointdata_vector.size();
            meas_per_setup_per_type[iSetup][1]=inputdata.setups[iSetup].head_inversion_data.mplist.pointdata_vector.size();
            for(UINT ii=0; ii< meas_per_setup_per_type[iSetup][1]; ii++){
              head_meas[iSetup].push_back( inputdata.setups[iSetup].head_inversion_data.mplist.pointdata_vector[ii] );
              //    all_meas.push_back( &(head_meas[iSetup][head_meas[iSetup].size()-1]) );
            }
          }

          //transport M0 meas!
          local_to_global[iSetup][2].resize(1);
          if(inputdata.problem_types.transport_inversion_m0){
            meas_per_setup[iSetup]+= inputdata.setups[iSetup].solute_concentration_inversion_data.m0_inversion_data.mplist.pointdata_vector.size();
            meas_per_setup_per_type[iSetup][2]=inputdata.setups[iSetup].solute_concentration_inversion_data.m0_inversion_data.mplist.pointdata_vector.size();
            for(UINT ii=0; ii< meas_per_setup_per_type[iSetup][2]; ii++){
              M0_meas[iSetup].push_back( inputdata.setups[iSetup].solute_concentration_inversion_data.m0_inversion_data.mplist.pointdata_vector[ii] );
              //   all_meas.push_back( &(M0_meas[iSetup][M0_meas[iSetup].size()-1]) );
            }
          }

          //transport M1 meas!
          local_to_global[iSetup][3].resize(1);
          if(inputdata.problem_types.transport_inversion_m1){
            meas_per_setup[iSetup]+= inputdata.setups[iSetup].solute_concentration_inversion_data.m1_inversion_data.mplist.pointdata_vector.size();
            meas_per_setup_per_type[iSetup][3]=inputdata.setups[iSetup].solute_concentration_inversion_data.m1_inversion_data.mplist.pointdata_vector.size();
            for(UINT ii=0; ii< meas_per_setup_per_type[iSetup][3]; ii++){
              M1_meas[iSetup].push_back( inputdata.setups[iSetup].solute_concentration_inversion_data.m1_inversion_data.mplist.pointdata_vector[ii] );
              // all_meas.push_back( &(M1_meas[iSetup][M1_meas[iSetup].size()-1]) );
            }
          }

          //temp. meas!
          local_to_global[iSetup][4].resize(1);
          if(inputdata.problem_types.heat_mean_arrival_time_inversion){
            meas_per_setup[iSetup]+= inputdata.setups[iSetup].heat_inversion_data.mplist.pointdata_vector.size();
            meas_per_setup_per_type[iSetup][4]=inputdata.setups[iSetup].heat_inversion_data.mplist.pointdata_vector.size();
            for(UINT ii=0; ii< meas_per_setup_per_type[iSetup][4]; ii++){
              heat_meas[iSetup].push_back( inputdata.setups[iSetup].heat_inversion_data.mplist.pointdata_vector[ii] );
              if(synthetic){
                heatM0_meas[iSetup].push_back( inputdata.setups[iSetup].heat_inversion_data.mpM0list.pointdata_vector[ii] );
                heatM1_meas[iSetup].push_back( inputdata.setups[iSetup].heat_inversion_data.mpM1list.pointdata_vector[ii] );
              }
              // all_meas.push_back( &(heat_meas[iSetup][heat_meas[iSetup].size()-1]) );
            }
          }

          //GE meas!
          GEmeas_per_setup_per_config[iSetup].resize(inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig,0);
          local_to_global[iSetup][5].resize(inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig);
          if(inputdata.problem_types.moments_geoeletric_potential_inversion){
            GE_meas[iSetup].resize(inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig);
            if(synthetic){
              GEM0_meas[iSetup].resize(inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig);
              GEM1_meas[iSetup].resize(inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig);
            }
            for(UINT iGE_config=0; iGE_config< inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig; iGE_config++){
              GEmeas_per_setup_per_config[iSetup][iGE_config] = inputdata.setups[iSetup].geoelectrical_potential_inversion_data.mplist[iGE_config].pointdata_vector.size();
              meas_per_setup[iSetup]+= inputdata.setups[iSetup].geoelectrical_potential_inversion_data.mplist[iGE_config].pointdata_vector.size();
              meas_per_setup_per_type[iSetup][5]+=inputdata.setups[iSetup].geoelectrical_potential_inversion_data.mplist[iGE_config].pointdata_vector.size();
              for(UINT ii=0; ii< GEmeas_per_setup_per_config[iSetup][iGE_config]; ii++){
                GE_meas[iSetup][iGE_config].push_back( inputdata.setups[iSetup].geoelectrical_potential_inversion_data.mplist[iGE_config].pointdata_vector[ii] );
                if(synthetic){
                  GEM0_meas[iSetup][iGE_config].push_back( inputdata.setups[iSetup].geoelectrical_potential_inversion_data.mpM0list[iGE_config].pointdata_vector[ii] );
                  GEM1_meas[iSetup][iGE_config].push_back( inputdata.setups[iSetup].geoelectrical_potential_inversion_data.mpM1list[iGE_config].pointdata_vector[ii] );
                }
                //all_meas.push_back( &(GE_meas[iSetup][iGE_config][GE_meas[iSetup][iGE_config].size()-1]) );

              }
            }
          }

          nMeas+=meas_per_setup[iSetup];
          // set the all_meas vector!
          local_to_global[iSetup][0][0].resize(meas_per_setup_per_type[iSetup][0],0);
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][0]; ii++){
            local_to_global[iSetup][0][0][ii]=all_meas.size();
            all_meas.push_back( &(lnK_meas[iSetup][ii]) );
          }

          local_to_global[iSetup][1][0].resize(meas_per_setup_per_type[iSetup][1],0);
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][1]; ii++){
            local_to_global[iSetup][1][0][ii]=all_meas.size();
            all_meas.push_back( &(head_meas[iSetup][ii]) );
          }

          local_to_global[iSetup][2][0].resize(meas_per_setup_per_type[iSetup][2],0);
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][2]; ii++){
            local_to_global[iSetup][2][0][ii]=all_meas.size();
            all_meas.push_back( &(M0_meas[iSetup][ii]) );
          }
          local_to_global[iSetup][3][0].resize(meas_per_setup_per_type[iSetup][3],0);
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][3]; ii++){
            local_to_global[iSetup][3][0][ii]=all_meas.size();
            all_meas.push_back( &(M1_meas[iSetup][ii]) );
          }
          local_to_global[iSetup][4][0].resize(meas_per_setup_per_type[iSetup][4],0);
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][4]; ii++){
            local_to_global[iSetup][4][0][ii]=all_meas.size();
            all_meas.push_back( &(heat_meas[iSetup][ii]) );
          }
          for(UINT ii =0; ii< inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig; ii++){
            local_to_global[iSetup][5][ii].resize(GEmeas_per_setup_per_config[iSetup][ii],0);
            for(UINT jj=0; jj< GEmeas_per_setup_per_config[iSetup][ii]; jj++  ){
              local_to_global[iSetup][5][ii][jj]=all_meas.size();
              all_meas.push_back( &(GE_meas[iSetup][ii][jj]) );
            }
          }
        }// END: for loop of setups

        if(synthetic){
          set_value_zero();
        }

      } // END: constructor


      //copy-constructor
      MeasurementList<MeasurementElement>(const MeasurementList<MeasurementElement>& ml)
      : inputdata(ml.inputdata)
        , synthetic(ml.synthetic)
        , nMeas(ml.nMeas)
        , nSetups(ml.nSetups)
        , meas_per_setup(ml.meas_per_setup)
        , meas_per_setup_per_type(ml.meas_per_setup_per_type)
        , GEmeas_per_setup_per_config(ml.GEmeas_per_setup_per_config)
        , local_to_global(ml.local_to_global)
        , lnK_meas(ml.lnK_meas)
        , head_meas(ml.head_meas)
        , M0_meas(ml.M0_meas)
        , M1_meas(ml.M1_meas)
        , heatM0_meas(ml.heatM0_meas)
        , heatM1_meas(ml.heatM1_meas)
        , heat_meas(ml.heat_meas)
        , GEM0_meas(ml.GEM0_meas)
        , GEM1_meas(ml.GEM1_meas)
        , GE_meas(ml.GE_meas)
      {
        all_meas.resize(0);

        for(UINT iSetup=0; iSetup<nSetups; iSetup++){

          // set the all_meas vector!

          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][0]; ii++){
            all_meas.push_back(&(lnK_meas[iSetup][ii]));
          }
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][1]; ii++){
            all_meas.push_back( &(head_meas[iSetup][ii]) );
          }
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][2]; ii++){
            all_meas.push_back( &(M0_meas[iSetup][ii]) );
          }
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][3]; ii++){
            all_meas.push_back( &(M1_meas[iSetup][ii]) );
          }
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][4]; ii++){
            all_meas.push_back( &(heat_meas[iSetup][ii]) );
          }

          //for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][5]; ii++){
          if(meas_per_setup_per_type[iSetup][5] > 0 ){

            for(UINT ii =0; ii<GEmeas_per_setup_per_config[iSetup].size(); ii++){

              //std::cout << "GEmeas_per_setup_per_config[iSetup]["<<ii<<"] = "
              //          << GEmeas_per_setup_per_config[iSetup][ii]
              //          << std::endl;

              for(UINT jj=0; jj< GEmeas_per_setup_per_config[iSetup][ii]; jj++  ){
                all_meas.push_back( &(GE_meas[iSetup][ii][jj]) );
              }
            }

          }

        }// END: for loop of setups

      } // END copy-constructor



      //
      // assignment operator
      //
      MeasurementList<MeasurementElement> & operator=(const MeasurementList<MeasurementElement> & rhs){
        if(this !=&rhs){
          this->~MeasurementList<MeasurementElement>();//destroy
          new (this) MeasurementList<MeasurementElement>(rhs);//reconstruct self using placement new and copy ctor
        }
        return *this;
      }


#ifdef DEBUG
      template<typename VecVecVecVec>
      void clear_t4( VecVecVecVec& vec ){
        if( vec.size()>0 ){
          for(int i=0;i<vec.size();i++){
            for(int j=0;j<vec[i].size();j++){
              for(int k=0;k<vec[i][j].size();k++){
                vec[i][j][k].clear();
              }
              vec[i][j].clear();
            }
            vec[i].clear();
          }
          vec.clear();
        }
      }

      template<typename VecVecVec>
      void clear_t3( VecVecVec& vec ){
        if( vec.size()>0 ){
          for(int i=0;i<vec.size();i++){
            if( vec[i].size()> 0 ){
              for(int j=0;j<vec[i].size();j++){
                if( vec[i][j].size()> 0 ){
                  vec[i][j].clear();
                }
              }
              vec[i].clear();
            }
          }
          vec.clear();
        }
      }

      template<typename VecVec>
      void clear_t2( VecVec& vec ){
        if( vec.size()>0 ){
          for(int i=0;i<vec.size();i++){
            vec[i].clear();
          }
          vec.clear();
        }
      }


      // Destructor:
      ~MeasurementList<MeasurementElement>(){
        all_meas.clear();

        clear_t3( GEM1_meas );
        clear_t3( GEM0_meas );
        clear_t3( GE_meas );

        clear_t2( heatM1_meas );
        clear_t2( heatM0_meas );
        clear_t2( heat_meas );
        clear_t2( M1_meas );
        clear_t2( M0_meas );
        clear_t2( head_meas );
        clear_t2( lnK_meas );

        clear_t4( local_to_global );

        clear_t2( GEmeas_per_setup_per_config );
        clear_t2( meas_per_setup_per_type );

        meas_per_setup.clear();

      }
#endif //#ifdef DEBUG



      void write_to_logger( std::string info="" ){
        logger << "Writing list >>" << info << "<< : " << std::endl;
        logger << nMeas << " Measurements in total!" << std::endl;
        logger<<nSetups<<" total Setups!"<<std::endl;
        for (UINT iSetup=0; iSetup<nSetups; iSetup++){
          logger<<"Setup #"<<iSetup+1<<" : "<<std::endl;
          logger<<meas_per_setup_per_type[iSetup][0]<<" ln K measurements:"<<std::endl;
          for(UINT ii=0; ii<meas_per_setup_per_type[iSetup][0]; ii ++)
            logger<<ii+1<<". : "<<lnK_meas[iSetup][ii].x<<" , "<<lnK_meas[iSetup][ii].y
#ifdef DIMENSION3
                  <<" , "<<lnK_meas[iSetup][ii].z
#endif
                  <<" ; value : "<<lnK_meas[iSetup][ii].value<<std::endl;
          logger<<meas_per_setup_per_type[iSetup][1]<<" head measurements:"<<std::endl;
          for(UINT ii=0; ii<meas_per_setup_per_type[iSetup][1]; ii ++)
            logger<<ii+1<<". : "<<head_meas[iSetup][ii].x<<" , "<<head_meas[iSetup][ii].y
#ifdef DIMENSION3
                  <<" , "<<head_meas[iSetup][ii].z
#endif
                  <<" ; value : "<<head_meas[iSetup][ii].value<<std::endl;
          logger<<meas_per_setup_per_type[iSetup][2]<<" solute concentration M0 measurements:"<<std::endl;
          for(UINT ii=0; ii<meas_per_setup_per_type[iSetup][2]; ii ++)
            logger<<ii+1<<". : "<<M0_meas[iSetup][ii].x<<" , "<<M0_meas[iSetup][ii].y
#ifdef DIMENSION3
                  <<" , "<<M0_meas[iSetup][ii].z
#endif
                  <<" ; value : "<<M0_meas[iSetup][ii].value<<std::endl;
          logger<<meas_per_setup_per_type[iSetup][3]<<" solute concentration M1 measurements:"<<std::endl;
          for(UINT ii=0; ii<meas_per_setup_per_type[iSetup][3]; ii ++)
            logger<<ii+1<<". : "<<M1_meas[iSetup][ii].x<<" , "<<M1_meas[iSetup][ii].y
#ifdef DIMENSION3
                  <<" , "<<M1_meas[iSetup][ii].z
#endif
                  <<" ; value : "<<M1_meas[iSetup][ii].value<<std::endl;
          logger<<meas_per_setup_per_type[iSetup][4]<<" temp. measurements:"<<std::endl;
          for(UINT ii=0; ii<meas_per_setup_per_type[iSetup][4]; ii ++)
            logger<<ii+1<<". : "<<heat_meas[iSetup][ii].x<<" , "<<heat_meas[iSetup][ii].y
#ifdef DIMENSION3
                  <<" , "<<heat_meas[iSetup][ii].z
#endif
                  <<" ; value : "<<heat_meas[iSetup][ii].value<<std::endl;

          logger<<meas_per_setup_per_type[iSetup][5]<<" GEmeasurements in "<<GEmeas_per_setup_per_config[iSetup].size()<<" configurations:"<<std::endl;
          if(meas_per_setup_per_type[iSetup][5]){
            for(UINT jj=0; jj<GEmeas_per_setup_per_config[iSetup].size(); jj++){
              logger<<"Configuration #"<<jj+1<<" has "<< GEmeas_per_setup_per_config[iSetup][jj]<<" measurements:"<<std::endl;
              for(UINT kk=0; kk<GEmeas_per_setup_per_config[iSetup][jj]; kk++){
                logger<<kk+1<<". : (x1)"<<GE_meas[iSetup][jj][kk].x<<" , (y1)"<<GE_meas[iSetup][jj][kk].y
#ifdef DIMENSION3
                      <<" , (z1)"<<GE_meas[iSetup][jj][kk].z
#endif
                      <<" ; (x2)"<<GE_meas[iSetup][jj][kk].x2<<" , (y2)"<<GE_meas[iSetup][jj][kk].y2
#ifdef DIMENSION3
                      <<" , (z2)"<<GE_meas[iSetup][jj][kk].z2
#endif
                      <<" ; value : "<<GE_meas[iSetup][jj][kk].value<<std::endl;
              }
            }
          }

        }
        logger<<"Writing list >>"<<info<<"<< : DONE!"<<std::endl;
      }



      void store_measurements( const UINT type, const UINT iSetup, const std::string filename ){

        Dune::Timer watch;

        switch( type ){
        case 1:
          {
            std::fstream file(filename,std::ios::out);
            if( file.is_open() ){
              file << "# \t x \t y "
#ifdef DIMENSION3
                   << "\t z "
#endif
                   << "\t h" << std::endl;
              for( int i=0; i<(int)head_meas[iSetup].size(); i++ ){
                file << i << " \t "
                     << head_meas[iSetup][i].x << " \t "
                     << head_meas[iSetup][i].y << " \t "
#ifdef DIMENSION3
                     << head_meas[iSetup][i].z << " \t "
#endif
                     << head_meas[iSetup][i].value << std::endl;
              }
            }
          }
          break;
        case 2:
          {
            std::fstream file(filename,std::ios::out);
            if( file.is_open() ){
              file << "# \t x \t y "
#ifdef DIMENSION3
                   << "\t z "
#endif
                   << "\t m0" << std::endl;
              for( int i=0; i<(int)M0_meas[iSetup].size(); i++ ){
                file << i << " \t "
                     << M0_meas[iSetup][i].x << " \t "
                     << M0_meas[iSetup][i].y << " \t "
#ifdef DIMENSION3
                     << M0_meas[iSetup][i].z << " \t "
#endif
                     << M0_meas[iSetup][i].value << std::endl;
              }
            }
          }
          break;
        case 3:
          {
            std::fstream file(filename,std::ios::out);
            if( file.is_open() ){
              file << "# \t x \t y "
#ifdef DIMENSION3
                   << "\t z "
#endif
                   << "\t m1" << std::endl;
              for( int i=0; i<(int)M1_meas[iSetup].size(); i++ ){
                file << i << " \t "
                     << M1_meas[iSetup][i].x << " \t "
                     << M1_meas[iSetup][i].y << " \t "
#ifdef DIMENSION3
                     << M1_meas[iSetup][i].z << " \t "
#endif
                     << M1_meas[iSetup][i].value << std::endl;
              }
            }
          }
          break;
        default:
          {
            std::cerr << ">>>>> ERROR: store_measurements: not implemented for the type "
                      << type << std::endl;
          }
        }

        std::stringstream jobtitle;
        jobtitle << "store_measurements: writing " << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   General::verbosity,
                                   "IO",
                                   jobtitle.str() );
      }



      bool read_measurements( const UINT type,
                              const UINT iSetup,
                              const std::string filename,
                              const Dune::MPIHelper& helper ){
        Dune::Timer watch;

        REAL x;
        REAL y;
        REAL current_x;
        REAL current_y;
#ifdef DIMENSION3
        REAL z;
        REAL current_z;
#endif

        std::string contents( Dune::Gesis::General::readTextFile2String( filename, helper ) );


        logger << "read_measurements: DEBUG LOG: contents = " << std::endl;
        logger << contents << std::endl;

        logger << "head_meas[iSetup].size() = " << head_meas[iSetup].size() << std::endl;


        std::stringstream instream( contents );
        REAL value;

        int iMeas = 0;
        int iLine = -1;
        bool bExitWhileLoop = false;

        std::random_device rd;
        std::mt19937_64 gen; // 64-bit Mersenne Twister

        int seed = rd();
        gen.seed( seed );
        std::normal_distribution<> distr(0,1);
        //std::uniform_real_distribution<REAL> distr(-1.0, 1.0);


        int nMeasPointsPerTypePerSetup = 0;
        while( !instream.eof() && !bExitWhileLoop ) {
          if(iLine > -1){


            switch( type ){
            case 1:
              {
                nMeasPointsPerTypePerSetup = int(head_meas[iSetup].size());

                logger << "iLine = " << iLine << std::endl;
                logger << "head_meas[iSetup].size() = " << nMeasPointsPerTypePerSetup << std::endl;
                if( iLine >= nMeasPointsPerTypePerSetup ){
                  logger << "read_measurements: WARNING: number of meas. points in input XML is less than in head data file." << std::endl;
                  bExitWhileLoop = true;
                  break;
                }

                current_x = head_meas[iSetup][iLine].x;
                current_y = head_meas[iSetup][iLine].y;
#ifdef DIMENSION3
                current_z = head_meas[iSetup][iLine].z;
#endif
                break;
              }
            case 2:
              {

                nMeasPointsPerTypePerSetup = int(M0_meas[iSetup].size());
                if( iLine >= nMeasPointsPerTypePerSetup ){
                  logger << "read_measurements: WARNING: number of meas. points in input XML is less than in m0 data file." << std::endl;
                  bExitWhileLoop = true;
                  break;
                }

                current_x = M0_meas[iSetup][iLine].x;
                current_y = M0_meas[iSetup][iLine].y;
#ifdef DIMENSION3
                current_z = M0_meas[iSetup][iLine].z;
#endif
                break;
              }
            case 3:
              {

                nMeasPointsPerTypePerSetup = int(M1_meas[iSetup].size());
                if( iLine >= nMeasPointsPerTypePerSetup ){
                  logger << "read_measurements: WARNING: number of meas. points in input XML is less than in m1 data file." << std::endl;
                  bExitWhileLoop = true;
                  break;
                }

                current_x = M1_meas[iSetup][iLine].x;
                current_y = M1_meas[iSetup][iLine].y;
#ifdef DIMENSION3
                current_z = M1_meas[iSetup][iLine].z;
#endif
                break;
              }
            default:
              std::cerr << "read_measurements: not yet implemented for the type " << type << std::endl;
            }

            if( bExitWhileLoop ){
              break; // exit loop over data file if the number of selected measurements in the inputXML is less than the number of lines in the data file
            }

            instream >> iMeas;

            instream >> x;
            if( std::abs(x - current_x) > 0.5 ){
              logger << "read_measurements:  faulty x = " << x << std::endl;
              logger << "read_measurements: current_x = " << current_x << std::endl;
              return false;
            }
            instream >> y;
            if( std::abs(y - current_y) > 0.5 ){
              logger << "read_measurements:  faulty y = " << y << std::endl;
              logger << "read_measurements: current_y = " << current_y << std::endl;
              return false;
            }
#ifdef DIMENSION3
            instream >> z;
            if( std::abs(z - current_z) > 0.5 ){
              logger << "read_measurements:  faulty z = " << z << std::endl;
              logger << "read_measurements: current_z = " << current_z << std::endl;
              return false;
            }
#endif
            instream >> value;



            REAL stochasticDisturbance = inputdata.inversion_parameters.disturbance;
            if( stochasticDisturbance > 1E-12 )
              stochasticDisturbance *= distr(gen);

            switch( type ){
            case 1:
              {
                value += 0.05 * stochasticDisturbance; // Allow maximally 1% disturbance for head measurement (100.0 - 99.90)x0.01
                head_meas[iSetup][iLine].value = value;
                break;
              }
            case 2:
              {
                value += value * stochasticDisturbance;
                M0_meas[iSetup][iLine].value = value;
                break;
              }
            case 3:
              {
                value += value * stochasticDisturbance;
                M1_meas[iSetup][iLine].value = value;
                break;
              }
            default:
              std::cerr << "read_measurements: not yet implemented for the type " << type << std::endl;
            }


            char c=instream.peek();
            if(c=='\n'){
              instream.readsome(&c,1);
              c=instream.peek();
              while (c==' '){
                instream.readsome(&c,1);
                c=instream.peek();
              }
            }
          }
          else{
            char dummyline[256];
            instream.getline(dummyline,256);
          }


          iLine++;
        }


        if( iLine < nMeasPointsPerTypePerSetup ){
          std::cout << "WARNING: mismatch in number of type#"
                    << type
                    << " - measuring points for setup \#"
                    << iSetup << ": "
                    << " data-file: " << iLine
                    << " input-XML: " << nMeasPointsPerTypePerSetup
                    << std::endl;
        }


        std::stringstream jobtitle;
        jobtitle << "read_measurements: reading " << filename;
        General::log_elapsed_time( watch.elapsed(),
                                   General::verbosity,
                                   "IO",
                                   jobtitle.str() );

        return true;
      }




      template<typename YFieldGenerator, typename GV>
      void take_measurements_lnK( const YFieldGenerator& data,
                                  const GV& gv,
                                  const Dune::MPIHelper& helper){

        if( !synthetic )
          return;  // No need to take measurements if the data in the inputfile are to be used.

        for(UINT iSetup=0; iSetup<nSetups; iSetup++ ){
          take_measurements_on_dgf_lnK( gv, data , lnK_meas[iSetup], helper );
          MPI_Tools::redist_measurements( lnK_meas[iSetup] );
        }
      }


      template<typename DATA, typename GFS>
      void take_measurements( UINT type,
                              const DATA& data,
                              const GFS& gfs,
                              const Dune::MPIHelper& helper,
                              UINT iSetup ){

        if( !synthetic )
          return;  // No need to take measurements if the data in the inputfile are to be used.

        // typedef typename GFS::Traits::GridViewType GV;
        //const GV& gv = gfs.gridView();

        switch (type){
        case 1:
          {

            //for( int i=0; i<head_meas[iSetup].size(); i++ ){
            // logger << "DEBUG0: head_meas[" << iSetup << "][" << i << "] = "
            //         << head_meas[iSetup][i].x << " \t "
            //         << head_meas[iSetup][i].y << " \t "
            //         << head_meas[iSetup][i].value << std::endl;
            //}

            take_measurements_on_dgf( gfs, data, head_meas[iSetup], helper );

            //for( int i=0; i<head_meas[iSetup].size(); i++ ){
            //  logger << "DEBUG1: head_meas[" << iSetup << "][" << i << "] = "
            //         << head_meas[iSetup][i].x << " \t "
            //         << head_meas[iSetup][i].y << " \t "
            //         << head_meas[iSetup][i].value << std::endl;
            //}

            MPI_Tools::redist_measurements( head_meas[iSetup] );

            //for( int i=0; i<head_meas[iSetup].size(); i++ ){
            //  logger << "DEBUG2: head_meas[" << iSetup << "][" << i << "] = "
            //         << head_meas[iSetup][i].x << " \t "
            //         << head_meas[iSetup][i].y << " \t "
            //         << head_meas[iSetup][i].value << std::endl;
            //}

          }
          break;
        case 2:
          {
            take_measurements_on_dgf( gfs, data, M0_meas[iSetup], helper );
            MPI_Tools::redist_measurements( M0_meas[iSetup] );
          }
          break;
        case 3:
          {
            take_measurements_on_dgf( gfs, data, M1_meas[iSetup], helper );
            MPI_Tools::redist_measurements( M1_meas[iSetup] );
          }
          break;
        default:
          logger<<"ERROR in take measurements. Meas. type not known or type function not used for this type of measurement(check for special functions)!!!"<<std::endl;
        }

      }



      template<typename DATA, // == VCType
               typename GFS>
      void take_measurements_AT( UINT type,
                                 const DATA& data1,
                                 const DATA& data2,
                                 const GFS& gfs,
                                 const Dune::MPIHelper& helper,
                                 UINT Setup,
                                 UINT config=0 ){

        if( !synthetic )
          return;  // No need to take measurements if the data in the inputfile are to be used.

        typedef typename GFS::Traits::GridViewType GV;
        const GV& gv = gfs.gridView();

        switch (type){
        case 4:
          {
            logger << "Take measurements of heatM1" << std::endl;
            take_measurements_on_dgf( gfs, data2, heatM1_meas[Setup], helper);
            logger << "Take measurements of heatM0" << std::endl;
            take_measurements_on_dgf( gfs, data1, heatM0_meas[Setup], helper);

            for(UINT ii=0; ii<heatM1_meas[Setup].size(); ii++){
              if(heatM1_meas[Setup][ii].value<0.0)
                heat_meas[Setup][ii].value = 0.0;
              else if (heatM0_meas[Setup][ii].value<0.0)
                heat_meas[Setup][ii].value = 0.0;
              else if(heatM0_meas[Setup][ii].value<GEO_EPSILON*0.5)
                heat_meas[Setup][ii].value = heatM1_meas[Setup][ii].value/(GEO_EPSILON*0.5);
              else
                heat_meas[Setup][ii].value = heatM1_meas[Setup][ii].value/heatM0_meas[Setup][ii].value;
              /*
                std::cout << "heat_meas = " << heat_meas[Setup][ii].value
                << " at location " << heat_meas[Setup][ii].x << "," << heat_meas[Setup][ii].y
                << std::endl;
              */
            }

            MPI_Tools::redist_measurements( heat_meas[Setup] );
            MPI_Tools::redist_measurements( heatM0_meas[Setup] );
            MPI_Tools::redist_measurements( heatM1_meas[Setup] );
          }
          break;

        case 5:
          {
            logger << "Take measurements of GEP M1 and GEP M0" << std::endl;
            logger<<"Setup "<<Setup<<std::endl;
            logger<<"config "<<config<<std::endl;
            logger<<"GEM1_meas.size() "<<GEM1_meas.size()<<std::endl;
            logger<<"GEM1_meas[Setup].size() "<<GEM1_meas[Setup].size()<<std::endl;
            logger<<"GEM1_meas[Setup][config].size() "<<GEM1_meas[Setup][config].size()<<std::endl;

            take_measurements_of_differences_on_dgf( gv, gfs, data2, GEM1_meas[Setup][config] , helper );
            take_measurements_of_differences_on_dgf( gv, gfs, data1, GEM0_meas[Setup][config] , helper );


            for(UINT ii=0; ii<GEM1_meas[Setup][config].size(); ii++){
              GE_meas[Setup][config][ii].value
                = GEM1_meas[Setup][config][ii].value / GEM0_meas[Setup][config][ii].value;

              logger <<"GEM0_meas[Setup][config][ii].value = "
                     <<GEM0_meas[Setup][config][ii].value
                     <<std::endl;

              logger <<"GEM1_meas[Setup][config][ii].value = "
                     <<GEM1_meas[Setup][config][ii].value
                     <<std::endl;

              logger <<"GE_meas[Setup][config][ii].value = "
                     <<GE_meas[Setup][config][ii].value
                     <<std::endl;
            }

            MPI_Tools::redist_measurements( GEM0_meas[Setup][config] );
            MPI_Tools::redist_measurements( GEM1_meas[Setup][config] );
            MPI_Tools::redist_measurements( GE_meas[Setup][config] );
          }

          break;
        default:
          logger<<"ERROR in take measurements. Meas. type not known or type function not used for this type of measurement(check for special functions)!!!"<<std::endl;
        }

      }

      //returns how many lnK measurements are there per setup
      UINT nMeasPerSetupPerType(UINT setup, UINT t) const {
        //if(t!=5)
        return meas_per_setup_per_type[setup][t];
        ////number of GE configurations will be returned!
        //return GEmeas_per_setup_per_config[setup].size()
      }

      //returns how many GE configurations are there per setup
      UINT nGE_config(UINT setup) const {
        return GEmeas_per_setup_per_config[setup].size();
      }

      //returns how many GE configurations are there per setup
      UINT nGE_Meas_per_setup_per_config(UINT setup,UINT config) const {
        return GEmeas_per_setup_per_config[setup][config];
      }

      //return the total number of measurements
      UINT nMeasurements() const {
        return nMeas;
      }

      MeasurementElement& operator[](const UINT& id) const {
        return *all_meas[id];
      }

      MeasurementElement& operator[](const UINT& id) {
        return *all_meas[id];
      }

      UINT get_global(UINT setup, UINT t, UINT i, UINT config=0) const {
        return local_to_global[setup][t][config][i];
      }

      REAL get_value_heatM0(UINT Setup, UINT local_ii) const {
        return heatM0_meas[Setup][local_ii].value;
      }

      REAL get_value_heatM1(UINT Setup, UINT local_ii) const {
        return heatM1_meas[Setup][local_ii].value;
      }

      REAL get_value_GEM0(UINT Setup, UINT GE_config, UINT local_ii) const {
        return GEM0_meas[Setup][GE_config][local_ii].value;
      }

      REAL get_value_GEM1(UINT Setup, UINT GE_config, UINT local_ii) const {
        return GEM1_meas[Setup][GE_config][local_ii].value;
      }


      REAL calibrateAbsErrorForM1( const UINT iSetup ) {
        UINT nPoints = M1_meas[iSetup].size();
        REAL maxValue = -1E+12;
        for(UINT iPoint=0; iPoint<nPoints; iPoint++){
          maxValue = std::max( maxValue,
                               M1_meas[iSetup][iPoint].value );
        }
        REAL absErrorForM1 = maxValue * inputdata.inversion_parameters.m1relerror;
        for(UINT iPoint=0; iPoint<nPoints; iPoint++){
          M1_meas[iSetup][iPoint].abs_error = absErrorForM1;
        }
        return absErrorForM1;
      }


      void set_value_zero(){
        for(UINT iSetup=0; iSetup<nSetups; iSetup++){
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][0]; ii++){
            lnK_meas[iSetup][ii].value=0.0;
          }
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][1]; ii++){
            head_meas[iSetup][ii].value=0.0;
          }
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][2]; ii++){
            M0_meas[iSetup][ii].value=0.0;
          }
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][3]; ii++){
            M1_meas[iSetup][ii].value=0.0;
          }
          for(UINT ii =0; ii< meas_per_setup_per_type[iSetup][4]; ii++){
            heat_meas[iSetup][ii].value=0.0;
            if(synthetic){
              heatM0_meas[iSetup][ii].value=0.0;
              heatM1_meas[iSetup][ii].value=0.0;
            }
          }

          for(UINT ii =0; ii< inputdata.setups[iSetup].geoelectrical_potential_inversion_data.nconfig; ii++)
            for(UINT jj=0; jj< GEmeas_per_setup_per_config[iSetup][ii]; jj++  ){
              GE_meas[iSetup][ii][jj].value=0.0;
              GEM0_meas[iSetup][ii][jj].value=0.0;
              GEM1_meas[iSetup][ii][jj].value=0.0;
            }
        }
      }

    };

  } // Gesis

} // Dune

#endif  /* DUNE_GESIS_MEASUREMENTLIST_HH */
