#ifndef DUNE_GEO_INVERSION_IO_LOCATIONS_HH
#define	DUNE_GEO_INVERSION_IO_LOCATIONS_HH

namespace Dune {
  namespace GeoInversion {

    //This class holds the most important paths and filenames that are being used repeatedly!
    class IO_Locations{
    public:

      std::string outputpath;
      std::string inputfile;
    
      std::string datadir;
      std::string bufferdir;
      std::string outputdir;

      std::string datadimdir;
      std::string bufferdimdir;
      std::string outputdimdir;

      std::string yfielddir;
      std::string CRdir;
      std::string logdir;
      std::string vtudir;
      std::string hdf5dir;
  
      std::string likelihood_file;
      std::string measurement_file;
  
      std::vector<std::string> R_YY_h5file;
      std::vector<std::string> EV_h5file;
      std::vector<std::string> zonation_matrix;
      std::string kfield_h5file;
      std::string extended_yfield_h5file;
      //std::string kfield_Well_h5file;
#ifndef WELL_FRACTURE_MODEL
      std::string kfield_h5file_without_wells;
#endif
  
  
      std::string unconditionalField_h5file;
      std::string kfield_properties_file;
      std::string Y_old_h5file;
      std::string Y_old2_h5file;  // same as Y_old, but on a finer grid
      //std::string Y_old_Well_h5file;
      std::string Y_old_xmffile;
      std::string Y_try_h5file;
      std::string Y_try_xmffile;
      std::string Y_estimated_h5file;
      std::string Y_estimated2_h5file;
  
      std::string ksi_old_h5file;
      std::string beta_old_h5file;
      std::string L_prior_h5file;
  
      std::vector<std::string> likelihood_cond_file;
      std::vector<std::string> measurement_cond_file;
      std::vector<std::string> Y_cond_h5file;
#ifdef PLOT_VTU_CR
      std::vector<std::string> Y_cond_vtu;
#endif
  
      std::string Y_est_h5file;

      std::vector<std::string> vchead_old_h5file;
      std::vector<std::string> vchead_est_h5file;

      std::vector<std::string> vcM0_old_h5file;
      std::vector<std::string> vcM0_est_h5file;

      std::vector<std::string> vcM1_old_h5file;
      std::vector<std::string> vcM1_est_h5file;
  
      std::vector<std::string> vchead_orig_h5file;
      std::vector<std::string> vcM0_cg_orig_h5file;
      std::vector<std::string> vcM1_cg_orig_h5file;
  
      std::vector<std::string> vcheatM0_cg_orig_h5file;
      std::vector<std::string> vcheatM1_cg_orig_h5file;
      std::vector<std::string> vcheatArrival_Time_orig_h5file;
  
      std::vector<std::string> vcheatM0_old_h5file;
      std::vector<std::string> vcheatM1_old_h5file;

      std::vector<std::vector<std::string> > vcM0phi_old_h5file;
      std::vector<std::vector<std::string> > vcM1phi_old_h5file;
  
      std::string logsigma0_h5file; 
      std::string logkappa_h5file;
      std::vector<std::vector<std::string> > vcphi0_orig_h5file;
      std::vector<std::vector<std::string> > vcM0phi_orig_h5file;
      std::vector<std::vector<std::string> > vcM1phi_orig_h5file;
      std::vector<std::vector<std::string> > vcATphi_orig_h5file;

  
      std::vector<std::string> JQ_h5file;
      std::vector<std::string> Sensitivity_h5file;
      std::vector<std::string> iSensitivity_h5file;
      std::vector<std::string> aSensitivity_h5file;
  
      std::string vtu_suffix;
      std::string Y_estimated_vtu;
  
#ifdef ESTIMATION_VARIANCE
      std::string estimation_variance_h5file;
      std::string V_est2_h5file;
      std::string estimation_variance_vtu;
      std::string cokriging_inverse_h5file;
#endif
  
#ifdef VTK_PLOT_ESTIMATES
      std::vector<std::string> vchead_estimated_vtu;
      std::vector<std::string> vcM0_estimated_vtu;
      std::vector<std::string> vcM1_estimated_vtu;
      std::vector<std::string> vcheatM0_estimated_vtu;
      std::vector<std::string> vcheatM1_estimated_vtu;
      std::vector<std::string> vcheatArrival_Time_estimated_vtu;
      std::vector<std::vector<std::string> > vcGPM0_estimated_vtu;
      std::vector<std::vector<std::string> > vcGPM1_estimated_vtu;
      std::vector<std::vector<std::string> > vcGPAT_estimated_vtu;
  
#endif
  
  
#ifdef VTK_PLOT_Y_FIELD
      std::string Y_orig_vtu;
#endif
#ifdef VTK_PLOT_P_FIELD
      std::vector<std::string> head_orig_vtu;
#endif
#ifdef VTK_PLOT_Q_FIELD
      std::vector<std::string> q_orig_vtu;
#endif
#ifdef VTK_PLOT_C_FIELD    
      std::vector<std::string> m0_orig_vtu;
      std::vector<std::string> m1_orig_vtu;
#endif 
    
#ifdef VTK_PLOT_HEAT_FIELD
      std::vector<std::string> heatM0_orig_vtu;
      std::vector<std::string> heatM1_orig_vtu;
      std::vector<std::string> heatArrival_Time_orig_vtu;
#endif
#ifdef VTK_PLOT_EL_POTENTIAL_FIELD
      std::vector<std::vector<std::string> > phi0_orig_vtu;
      std::vector<std::vector<std::string> > M0phi_orig_vtu;
      std::vector<std::vector<std::string> > M1phi_orig_vtu;
      std::vector<std::vector<std::string> > ATphi_orig_vtu;
#endif
    
#ifdef VTK_PLOT_PSI_HEAD
      std::string vchead_adj_prefix;
      std::string vchead_heat_tomography_adj_prefix;
      std::string vchead_GP_adj_prefix;
#endif
#ifdef VTK_PLOT_PSI_transport
      std::string vcM0_adj_prefix;
      std::string vchead_adj_given_M0_prefix;
      std::string vcM1_adj_prefix;
      std::string vcM0_adj_given_M1_prefix;
      std::string vchead_adj_given_M0M1_prefix;
#endif
#ifdef VTK_PLOT_PSI_heat_transport
      std::string vcheatM0_adj_prefix;
      std::string vchead_adj_given_heatM0_prefix;
      std::string vcheatM0_adj_given_heatM1_prefix;
      std::string vchead_adj_given_heatM0heatM1_prefix;
#endif
#ifdef VTK_PLOT_PSI_GP
      std::string vcphi_adj_prefix;
      std::string vcM0gp_adj_phi_adj_prefix;
      std::string vcheadgp_adj_M0_prefix;
      std::string vcM0gp_adj_M1_prefix;
      std::string vcheadgp_adj_M0M1_prefix;
#endif
#ifdef VTK_PLOT_Y_OLD 
      std::string Y_old_vtu_prefix;
      std::string Y_old_pvd;
#endif
#ifdef VTK_PLOT_S
      std::string Sensitivity_vtu_prefix;
#endif
#ifdef VTK_PLOT_JQ
      std::string JQ_vtu_prefix;
#endif
#ifdef VTK_PLOT_YTRY
      std::string Y_try_vtu;
#endif
   
  
      IO_Locations( const int dim,
                    const std::string _outputpath,
                    const std::string _inputfile ) 
        : outputpath(_outputpath),
          inputfile(_inputfile) {

        std::stringstream sDim;
        sDim << dim;
    
        datadir   = outputpath+"DATA";

#if defined TEST_GLOBAL_REFINE

#ifdef USE_FEM
        bufferdir    = outputpath+"TEST_BUFFER_FEM";
        outputdir    = outputpath+"TEST_OUTPUT_FEM";
#else
        bufferdir    = outputpath+"TEST_BUFFER_DG";
        outputdir    = outputpath+"TEST_OUTPUT_DG";
#endif

#else

#ifdef USE_FEM
        bufferdir    = outputpath+"BUFFER_FEM";
        outputdir    = outputpath+"OUTPUT_FEM";
#else
        bufferdir    = outputpath+"BUFFER_DG";
        outputdir    = outputpath+"OUTPUT_DG";
#endif

#endif

        bufferdimdir = bufferdir+"/" + sDim.str() + "D";
        outputdimdir = outputdir+"/" + sDim.str() + "D";
        logdir       = outputdir+"/" + sDim.str() + "D/LOG";
        vtudir       = outputdir+"/" + sDim.str() + "D/VTU";
        hdf5dir      = outputdir+"/" + sDim.str() + "D/HDF5";

        datadimdir   = outputpath+"DATA/"   + sDim.str() + "D";

        yfielddir = outputpath+"DATA/"   + sDim.str() + "D/YField";
        CRdir     = outputpath+"DATA/"   + sDim.str() + "D/CR";

        likelihood_file=logdir+"/Likelihood.txt";
        measurement_file=logdir+"/measurements.txt";
    
    
        kfield_properties_file=yfielddir+"/yfield_properties.dat";
        kfield_h5file = yfielddir+"/Yfield.h5";
        extended_yfield_h5file = yfielddir+"/ext_Yfield.h5";
        //kfield_Well_h5file=yfielddir+"/Yfield_Well.h5";
#ifndef WELL_FRACTURE_MODEL
        kfield_h5file_without_wells=yfielddir+"/Yfield_without_wells.h5";
#endif
        unconditionalField_h5file=bufferdimdir+"/Y_unconditional.h5";

        Y_old_h5file  = bufferdimdir + "/Y_old.h5";
        Y_old2_h5file  = bufferdimdir + "/Y_old2.h5";
        //Y_old_Well_h5file  = bufferdimdir + "/Y_old_Well.h5";
        Y_old_xmffile = bufferdimdir + "/Y_old.xmf";

        Y_try_h5file  = bufferdimdir + "/Y_try.h5";
        Y_try_xmffile = bufferdimdir + "/Y_try.xmf";
    
        Y_estimated_h5file = bufferdimdir + "/Y_estimated.h5";
        Y_estimated2_h5file = bufferdimdir + "/Y_estimated2.h5"; // same as Y_estimated, but on a finer grid
    
        ksi_old_h5file    =bufferdimdir + "/ksi_old.h5";
        beta_old_h5file   =bufferdimdir + "/beta_old.h5";
        L_prior_h5file    =bufferdimdir + "/L_prior.h5";

    
        logsigma0_h5file=yfielddir+"/logsigma0_field.h5";
        logkappa_h5file=yfielddir+"/logkappa_field.h5";

        vtu_suffix=".vtu";
        Y_estimated_vtu = vtudir + "/Y_estimated";
    
#ifdef ESTIMATION_VARIANCE
        estimation_variance_h5file = bufferdimdir + "/Estimation_Variance.h5";
        V_est2_h5file = bufferdimdir + "/V_est2.h5";
        estimation_variance_vtu = vtudir + "/V_est";
        cokriging_inverse_h5file = bufferdimdir + "/cokriging_inverse.h5";
#endif

#ifdef VTK_PLOT_Y_FIELD    
        Y_orig_vtu = vtudir + "/Y_orig";
#endif

#ifdef VTK_PLOT_YTRY
        Y_try_vtu = vtudir + "/Y_try";
#endif

#ifdef VTK_PLOT_PSI_HEAD
        vchead_adj_prefix=vtudir+"/vcHead_Adj";
        vchead_heat_tomography_adj_prefix=vtudir+"/vcHead_Heat_Adj";
        vchead_GP_adj_prefix=vtudir+"/vcHead_GP_Adj";
#endif
#ifdef VTK_PLOT_PSI_transport
        vcM0_adj_prefix=vtudir+"/vcM0_Adj";
        vchead_adj_given_M0_prefix=vtudir+"/vchead_Adj_given_M0";
        vcM1_adj_prefix=vtudir+"/vcM1_Adj";
        vcM0_adj_given_M1_prefix=vtudir+"/vcM0_Adj_given_M1";
        vchead_adj_given_M0M1_prefix=vtudir+"/vchead_Adj_given_M0M1";
#endif
#ifdef VTK_PLOT_PSI_heat_transport
        vcheatM0_adj_prefix=vtudir+"/vcheatM0_Adj";
        vchead_adj_given_heatM0_prefix=vtudir+"/vchead_heatM0_Adj";
        vcheatM0_adj_given_heatM1_prefix=vtudir+"/vcheatM0_Adj_given_heatM1";
        vchead_adj_given_heatM0heatM1_prefix=vtudir+"/vchead_Adj_given_heatM0heatM1";
#endif
#ifdef VTK_PLOT_PSI_GP
        vcphi_adj_prefix=vtudir+"/vcphi_Adj";
        vcM0gp_adj_phi_adj_prefix=vtudir+"/vcM0GP_adj_phi_Adj";
        vcheadgp_adj_M0_prefix=vtudir+"/vcheadGP_adj_M0_Adj";
        vcM0gp_adj_M1_prefix=vtudir+"/vcM0GP_adj_M1_Adj";
        vcheadgp_adj_M0M1_prefix=vtudir+"/vcheadGP_adj_M0M1_Adj";
#endif
#ifdef VTK_PLOT_Y_OLD 
        Y_old_vtu_prefix=vtudir + "/Y_old_";
        Y_old_pvd=vtudir + "/Y_old.pvd";
#endif
#ifdef VTK_PLOT_S
        Sensitivity_vtu_prefix=vtudir + "/Sensitivity";
#endif
#ifdef VTK_PLOT_JQ
        JQ_vtu_prefix=vtudir + "/JQ";
#endif
      }
  
      void set_zonation(const UINT& nzones){
        EV_h5file.resize(nzones);
        zonation_matrix.resize(nzones);
        R_YY_h5file.resize(nzones);
        for(UINT ii=0; ii<nzones; ii++){
          std::stringstream sii;
          sii << ii;
          EV_h5file[ii] = yfielddir + "/EV_Zone_" + sii.str() + ".h5";
          zonation_matrix[ii] = yfielddir + "/X_Zone_" + sii.str() + ".h5";
          R_YY_h5file[ii] = yfielddir + "/R_YY_Zone_" + sii.str() + ".h5";
        }
      }
  
  
      void set_setups(const UINT& setups){
        vchead_orig_h5file.resize(setups);
        vcM0_cg_orig_h5file.resize(setups);
        vcM1_cg_orig_h5file.resize(setups);
        vcheatM0_cg_orig_h5file.resize(setups);
        vcheatM1_cg_orig_h5file.resize(setups);
        vcheatArrival_Time_orig_h5file.resize(setups);
      
        vcphi0_orig_h5file.resize(setups);
        vcM0phi_orig_h5file.resize(setups);
        vcM1phi_orig_h5file.resize(setups);
        vcATphi_orig_h5file.resize(setups);
      
        vchead_old_h5file.resize(setups);
        vchead_est_h5file.resize(setups);
      
        vcM0_old_h5file.resize(setups);
        vcM0_est_h5file.resize(setups);
        vcM1_old_h5file.resize(setups);
        vcM1_est_h5file.resize(setups);
      
        vcheatM0_old_h5file.resize(setups);
        vcheatM1_old_h5file.resize(setups);
      
        vcM0phi_old_h5file.resize(setups);
        vcM1phi_old_h5file.resize(setups);
      
      
#ifdef VTK_PLOT_P_FIELD 
        head_orig_vtu.resize(setups);
#endif
#ifdef VTK_PLOT_Q_FIELD
        q_orig_vtu.resize(setups);
#endif   
#ifdef VTK_PLOT_C_FIELD    
        m0_orig_vtu.resize(setups);
        m1_orig_vtu.resize(setups);
#endif   
#ifdef VTK_PLOT_HEAT_FIELD
        heatM0_orig_vtu.resize(setups);
        heatM1_orig_vtu.resize(setups);
        heatArrival_Time_orig_vtu.resize(setups);
#endif
    
#ifdef VTK_PLOT_EL_POTENTIAL_FIELD
        phi0_orig_vtu.resize(setups);
        M0phi_orig_vtu.resize(setups);
        M1phi_orig_vtu.resize(setups);
        ATphi_orig_vtu.resize(setups);
#endif
    
#ifdef VTK_PLOT_ESTIMATES
        vchead_estimated_vtu.resize(setups);
        vcM0_estimated_vtu.resize(setups);
        vcM1_estimated_vtu.resize(setups);
        vcheatM0_estimated_vtu.resize(setups);
        vcheatM1_estimated_vtu.resize(setups);
        vcheatArrival_Time_estimated_vtu.resize(setups);
        vcGPM0_estimated_vtu.resize(setups);
        vcGPM1_estimated_vtu.resize(setups);
        vcGPAT_estimated_vtu.resize(setups);
#endif

        Y_est_h5file = hdf5dir + "/Y_est.h5";
      
        for(UINT ii=0; ii<setups; ii++){
          std::stringstream sii;
          sii << ii;
          vchead_orig_h5file[ii] = bufferdimdir + "/vchead_s" +sii.str() + "_orig.h5";
          vcM0_cg_orig_h5file[ii] = bufferdimdir + "/vcM0_cg_s"+sii.str() +"_orig.h5";
          vcM1_cg_orig_h5file[ii] = bufferdimdir + "/vcM1_cg_s"+sii.str() +"_orig.h5";
        
          vcheatM0_cg_orig_h5file[ii] = bufferdimdir + "/vcheatM0_cg_s" + sii.str() + "_orig.h5";
          vcheatM1_cg_orig_h5file[ii] = bufferdimdir + "/vcheatM1_cg_s" + sii.str() + "_orig.h5";
          vcheatArrival_Time_orig_h5file[ii] = bufferdimdir + "/vcheatArrival_Time_s" + sii.str() + "_orig.h5";
        
          vchead_old_h5file[ii] = bufferdimdir + "/vchead_old_s" + sii.str() + ".h5";
          vchead_est_h5file[ii] = hdf5dir + "/vchead_est_s" + sii.str() + ".h5";

          vcM0_old_h5file[ii] = bufferdimdir + "/vcM0_old_s" + sii.str() + ".h5";
          vcM0_est_h5file[ii] = hdf5dir + "/vcM0_est_s" + sii.str() + ".h5";

          vcM1_old_h5file[ii] = bufferdimdir + "/vcM1_old_s" + sii.str() + ".h5";
          vcM1_est_h5file[ii] = hdf5dir + "/vcM1_est_s" + sii.str() + ".h5";
        
          vcheatM0_old_h5file[ii] = bufferdimdir + "/vcheatM0_s" + sii.str() + "_old.h5";
          vcheatM1_old_h5file[ii] = bufferdimdir + "/vcheatM1_s" + sii.str() + "_old.h5";
        
#ifdef VTK_PLOT_P_FIELD
          head_orig_vtu[ii] = vtudir + "/head_orig_"+sii.str();
#endif
#ifdef VTK_PLOT_Q_FIELD
          q_orig_vtu[ii] = vtudir + "/q_orig_"+sii.str();
#endif
#ifdef VTK_PLOT_C_FIELD    
          m0_orig_vtu[ii] = vtudir + "/m0_orig_"+sii.str();
          m1_orig_vtu[ii] = vtudir + "/m1_orig_"+sii.str();
#endif 
#ifdef VTK_PLOT_HEAT_FIELD
          heatM0_orig_vtu[ii] = vtudir + "/heatM0_orig_"+sii.str();
          heatM1_orig_vtu[ii] = vtudir + "/heatM1_orig_"+sii.str();
          heatArrival_Time_orig_vtu[ii]=vtudir + "/heatArrival_Time_orig_"+sii.str();
#endif  
#ifdef VTK_PLOT_ESTIMATES
          vchead_estimated_vtu[ii]=vtudir + "/head_estimated_"+sii.str();
          vcM0_estimated_vtu[ii]=vtudir + "/SoluteM0_estimated_"+sii.str();
          vcM1_estimated_vtu[ii]=vtudir + "/SoluteM1_estimated_"+sii.str();
          vcheatM0_estimated_vtu[ii]=vtudir + "/HeatM0_estimated_"+sii.str();
          vcheatM1_estimated_vtu[ii]=vtudir + "/HeatM1_estimated_"+sii.str();
          vcheatArrival_Time_estimated_vtu[ii]=vtudir + "/HeatArrival_Time_estimated_"+sii.str();
#endif
        }
      }

      void set_geoelectrical_potential_tomography(const UINT& iSetup, const UINT& nconfigs){
        vcphi0_orig_h5file[iSetup].resize(nconfigs);
        vcM0phi_orig_h5file[iSetup].resize(nconfigs);
        vcM1phi_orig_h5file[iSetup].resize(nconfigs);
        vcATphi_orig_h5file[iSetup].resize(nconfigs);
    
        vcM0phi_old_h5file[iSetup].resize(nconfigs);
        vcM1phi_old_h5file[iSetup].resize(nconfigs);
    

#ifdef VTK_PLOT_EL_POTENTIAL_FIELD
        phi0_orig_vtu[iSetup].resize(nconfigs);
        M0phi_orig_vtu[iSetup].resize(nconfigs);
        M1phi_orig_vtu[iSetup].resize(nconfigs);
        ATphi_orig_vtu[iSetup].resize(nconfigs);
#endif
#ifdef VTK_PLOT_ESTIMATES  
        vcGPM0_estimated_vtu[iSetup].resize(nconfigs);
        vcGPM1_estimated_vtu[iSetup].resize(nconfigs);
        vcGPAT_estimated_vtu[iSetup].resize(nconfigs);
#endif
    
        for(UINT ii=0; ii<nconfigs; ii++){
          std::stringstream sjj;
          sjj << iSetup;
          std::stringstream sii;
          sii << ii;
        
          vcphi0_orig_h5file[iSetup][ii] = bufferdimdir + "/vcphi0_" + sjj.str() + "_config_" + sii.str() + "_orig.h5";
          vcM0phi_orig_h5file[iSetup][ii] = bufferdimdir + "/vcM0phi_" + sjj.str() + "_config_"  + sii.str() + "_orig.h5";
          vcM1phi_orig_h5file[iSetup][ii] = bufferdimdir + "/vcM1phi_" + sjj.str() + "_config_"  + sii.str() + "_orig.h5";
          vcATphi_orig_h5file[iSetup][ii] = bufferdimdir + "/vcATphi_" + sjj.str() + "_config_" + sii.str() + "_orig.h5";
        
          vcM0phi_old_h5file[iSetup][ii] = bufferdimdir + "/vcGPM0_old_" + sjj.str() + "_config_" + sii.str() + "_orig.h5";
          vcM1phi_old_h5file[iSetup][ii] = bufferdimdir + "/vcGPM1_old_" + sjj.str() + "_config_" + sii.str() + "_orig.h5";
        
#ifdef VTK_PLOT_EL_POTENTIAL_FIELD
          phi0_orig_vtu[iSetup][ii] = vtudir + "/phi0_" + sjj.str() + "_config_" + sii.str() + "_orig";
          M0phi_orig_vtu[iSetup][ii] = vtudir + "/M0phi_" + sjj.str() + "_config_" + sii.str() + "_orig";
          M1phi_orig_vtu[iSetup][ii] = vtudir + "/M1phi_" + sjj.str() + "_config_" + sii.str() + "_orig";
          ATphi_orig_vtu[iSetup][ii] = vtudir + "/ATphi_" + sjj.str() + "_config_" + sii.str() + "_orig";
#endif
#ifdef VTK_PLOT_ESTIMATES       
          vcGPM0_estimated_vtu[iSetup][ii]=vtudir + "/GP_M0_estimated_"+ sjj.str() + "_config_" + sii.str();
          vcGPM1_estimated_vtu[iSetup][ii]=vtudir + "/GP_M1_estimated_"+ sjj.str() + "_config_" + sii.str();
          vcGPAT_estimated_vtu[iSetup][ii]=vtudir + "/GP_AT_estimated_"+ sjj.str() + "_config_" + sii.str();
#endif
        }
      };


      void set_inversion(const UINT nMeas){
        JQ_h5file.resize(nMeas);
        Sensitivity_h5file.resize(nMeas);
        iSensitivity_h5file.resize(nMeas);
        aSensitivity_h5file.resize(nMeas);
        for(UINT ii=0; ii<nMeas; ii++){
          std::stringstream sii;
          sii << ii;
          JQ_h5file[ii] = bufferdimdir + "/JQ_" + sii.str() + ".h5";
          Sensitivity_h5file[ii] = bufferdimdir + "/sens_" + sii.str() + ".h5";
          iSensitivity_h5file[ii] = bufferdimdir + "/iSens_" + sii.str() + ".h5";
          aSensitivity_h5file[ii] = bufferdimdir + "/aSens_" + sii.str() + ".h5";
        }
      };
      

      void set_CR(const UINT& nCR,const UINT& start){
        likelihood_cond_file.resize(nCR);
        measurement_cond_file.resize(nCR);
        Y_cond_h5file.resize(nCR);
#ifdef PLOT_VTU_CR
        Y_cond_vtu.resize(nCR);
#endif  
        for(UINT ii=0; ii<nCR; ii++){
          std::stringstream sii;
          sii << ii+start;

          likelihood_cond_file[ii]=logdir+"/Likelihood_cond_#"+sii.str()+".txt";
          measurement_cond_file[ii]=logdir+"/measurements_cond_#"+sii.str()+".txt";
          Y_cond_h5file[ii]= CRdir + "/Y_cond_#"+sii.str()+".h5";
#ifdef PLOT_VTU_CR
          Y_cond_vtu[ii]= vtudir + "/Y_cond_#"+sii.str();
#endif
        }
      };
      
    }; // class IO_locations

  }

}

#endif // DUNE_GEO_INVERSION_IO_LOCATIONS_HH

