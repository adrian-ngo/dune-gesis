/*
 * File:   inputfile.hh
 * Author: A. Ngo
 *
 * Created on June 23, 2010
 * Ported to new XML Parser (boost property tree) on October 28, 2012
 *
 * 2014/08/01: Add timer
 *
 */
#ifndef DUNE_GESIS_INPUTFILE_HH
#define	DUNE_GESIS_INPUTFILE_HH

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


/***
    Boost.PropertyTree (http://www.boost.org/doc/libs/1_57_0/doc/html/property_tree.html)
    Copyright © 2008 Marcin Kalicinski
    Distributed under the Boost Software License, Version 1.0. 
    (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 ***/
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/detail/xml_parser_error.hpp>

/***
    Boost.Foreach (http://www.boost.org/doc/libs/1_57_0/doc/html/foreach.html)
    Copyright © 2004 Eric Niebler
    Distributed under the Boost Software License, Version 1.0. 
    (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 ***/
#include <boost/foreach.hpp>


#include "dune/gesis/common/Vector.hh"
#include "dune/gesis/common/general.hh"


/* Get global logfile object */
#include "logfile.hh"
extern CLogfile logger;


namespace Dune{
  namespace Gesis{


    struct CLinePlot {
      std::vector<REAL> startpoint;
      std::vector<REAL> endpoint;
    };


    class CDomainData {

    public:
      UINT dim;
      bool structured;
      std::vector<CTYPE> extensions;
      Vector<UINT> nCells;       // = virtual_nCells = number of cells per dimension in the virtual grid!
      std::vector<UINT> nPartitions;
      std::vector<UINT> nChunks;
      CLinePlot lineplot;

      Vector<CTYPE> virtual_gridsizes;   // resolution of our conductivity field
      Vector<CTYPE> baselevel_gridsizes; // resolution of Groundwater Grid

#ifdef USE_YASP
      std::vector<UINT> yasp_nCells; // = number of cells per dimension in the YaspGrid!
      std::vector<CTYPE> yasp_gridsizes; // = the length of a YaspGrid element per dimension, after global refinement
      UINT yasp_baselevel;
      UINT yasp_maxlevel;
      UINT yasp_overlap;
#endif

#ifdef USE_ALUGRID
      std::vector<CTYPE> alugrid_gridsizes; // = the length of a AluGrid element along coordinate axis per dimension, after global refinement
      UINT alugrid_globalrefine;
      UINT alugrid_baselevel;
      UINT alugrid_maxsteps;
      REAL alugrid_refinementfraction;
      REAL alugrid_coarseningfraction;
      REAL alugrid_tolerance;
      UINT alugrid_strategy;
#endif

#ifdef USE_UG
      std::vector<CTYPE> ug_gridsizes; // = the length of a UG element along coordinate axis per dimension, after global refinement
      UINT ug_heapsize;
      UINT ug_globalrefine;
      UINT ug_baselevel;
      UINT ug_maxsteps;
      std::string ug_refinementtype;
      REAL ug_refinementfraction;
      REAL ug_coarseningfraction;
      REAL ug_tolerance;
      UINT ug_strategy;
#endif

      CDomainData() {
      };

      ~CDomainData() {
      };

      // Important: This function must be invoked diretly after object instantiation in order to reserve enough memory for the std-vectors!

      void initialize(UINT dim_) {
        dim = dim_;
        extensions.resize(dim);
        nCells.resize(dim);

        virtual_gridsizes.resize(dim);
        baselevel_gridsizes.resize(dim);
        nChunks.resize(dim);
        lineplot.startpoint.resize(dim);
        lineplot.endpoint.resize(dim);
        nPartitions.resize(dim);
#ifdef USE_YASP
        yasp_nCells.resize(dim);
        yasp_gridsizes.resize(dim);
#endif
#ifdef USE_ALUGRID
        alugrid_gridsizes.resize(dim);
#endif
#ifdef USE_UG
        ug_gridsizes.resize(dim);
#endif
      }
    };



    /**********************************************************************************
     *
     *  CProblemTypes is an object map of the XML-node <problem_types>
     *
     **********************************************************************************
     */
    struct CProblemTypes {
      bool new_YField;           // inside XML: new_Yfield="yes" => Generate a new random Y-field.
      bool new_Eigenvalues;      // optional flag: inside XML: new_Eigenvalues="yes" => Generate new Eigenvalues out of the domain data.
      bool synthetic;
      bool using_existing_Yold;  // set by commandline argument --continue or -c
      int refine_estimate;       // Y_est gets replicated on a finer level
      bool head_forward;         //
      bool lnK_inversion;        //
      bool head_inversion;       //
      bool transport_forward;    //
      bool heat_forward;    //
      bool heat_mean_arrival_time_inversion;
      bool moments_geoeletric_potential_forward;
      bool moments_geoeletric_potential_inversion;
      bool transport_inversion_m0;  //
      bool transport_inversion_m1;  //
      bool CR;
      bool estimation_variance_only;
      bool generate_measurement_data;
    };

    /**********************************************************************************
     *
     *  CPlotOptions is an object map of the XML-node <plot_options>
     *
     **********************************************************************************
     */
    struct CPlotOptions {
      bool vtk_plot_yfield;
      bool vtk_plot_head;
      bool vtk_plot_q;
      bool vtk_plot_m0;
      bool vtk_plot_m1;
      bool vtk_plot_element_ordering;
      bool vtk_plot_el_potential_field;
      bool vtk_plot_heat;
      bool vtk_plot_wells;
      bool vtk_plot_adjoint_head;
      bool vtk_plot_adjoint_m0;
      bool vtk_plot_adjoint_m1;
      bool vtk_plot_heat_adjoint;
      bool vtk_plot_gp_adjoint;
      bool vtk_plot_sensitivities;
      bool vtk_plot_cross_covariances;
      bool vtk_plot_y_smooth;
      bool vtk_plot_y_try;
      bool vtk_plot_y_old;
      bool vtk_plot_h_old;
      bool vtk_plot_m0_old;
      bool vtk_plot_trialfields;
    };


    /**********************************************************************************
     *
     *  CYFieldZone is an object map of the XML-node <Yfield_zone>
     *
     **********************************************************************************
     */
    struct CYFieldZone{
    public:
      std::string model;
      REAL bottom;
      REAL beta;
      REAL qbb_y;
      REAL variance;
      REAL porosity;
      REAL rho;
      REAL c_s;
      REAL lambda_s;
      REAL kappa;
      REAL sigma0;
      Vector<REAL> correlation_lambda;
      REAL embedding_factor;

      void initialize(UINT dim_) {
        correlation_lambda.resize(dim_);
      }
    };



    /**********************************************************************************
     *
     *  CYFieldProperties is an object map of the XML-node <yfield_properties>
     *
     **********************************************************************************
     */
    struct CYFieldProperties{
    public:
      UINT random_seed;
      UINT nz;
      std::vector<CYFieldZone> zones;
      std::string well_type;

      void initialize(UINT nz_) {
        zones.resize(nz_);
        nz=nz_;
      }
    };


    class CStripe {
    public:
      std::vector<REAL> startpoint;
      std::vector<REAL> endpoint;
      REAL value;
      REAL regularization_factor;
      bool bfixedwidth;
    };

    class CBoundary {
    public:
      std::string orientation;
      int bctype;
      std::vector<CStripe> stripes; // Take care: This must be resized later!
    };

    class CBoundaryValueProblem {
    public:
      UINT nBoundaries;
      bool bRecycleMatrixHierarchy;
      REAL injection_time;
      //std::string active;
      std::vector<CBoundary> boundaries;
      CBoundaryValueProblem(){
        bRecycleMatrixHierarchy= false;
        injection_time=0.0;
      }

      void initialize(UINT dim_) {
        if( dim_ == 3 )
          nBoundaries = 6;
        else
          nBoundaries = 4;
        boundaries.resize( nBoundaries );
      }
    };



    /**********************************************************************************
     *
     *  CTransportParameters is an object map of the XML-node <transport_parameters>
     *
     **********************************************************************************
     */
    struct CTransportParameters {
      REAL a_longitudinal;
      REAL a_transversal;
      REAL D_molecular;
      REAL refined_diffusion; // extra diffusion for locally refined element
      REAL deltaSD_factor;
      REAL tolerance;  // tolerance for over- and undershoots
      REAL l2_diffusion; // extra diffusion for L2-projected solution
      REAL l2_diffusion_adjoint; // extra diffusion for L2-projected adjoint
      std::string dg_method;
      std::string dg_weights;
      REAL dg_gamma;
      int istl_max_iter;
      REAL istl_reduction;
      std::string istl_atOnceAccu;
      REAL point_smearing;    // regularization factor for point source

      /**********************************************************************************
       *
       *  This part is an object map of the XML-node <heat_transport_parameters>
       *
       **********************************************************************************
       */
      REAL heat_a_longitudinal;
      REAL heat_a_transversal;
      REAL heat_rho_w;
      REAL heat_c_w;
      REAL heat_lambda_w;
      REAL heat_deltaSD_factor;

    };


    enum MEASURE_TYPE {
      UNDEF_MEASURE_TYPE,
      PUMPING_DATA,
      WELL_DATA,
      LN_K_DATA,
      HEAD_DATA,
      M0_DATA,
      M1_DATA,
      HEAT_M0_DATA,
      HEAT_M1_DATA,
      GP_M0_DATA,
      GP_M1_DATA,
      GP_IN_OUT_DATA
    };


    struct SMeasurementTypes {
      std::string head;
    };



    // The PointType specifies input data format.
    enum PointType {
      UNDEFINED,
      SOURCE_POINT,   // XML expects the 2d(3d)-format: x, y, (z,) value
      WELL_POINT,     // XML expects the 2d(3d)-format: x, (y,) zTop, zBottom, k_well, flow_rate
      //  INJECTION_SOURCE_POINT,   // XML expects the 2d(3d)-format: x, y, (z,) value
      //  INJECTION_WELL_POINT,     // XML expects the 2d(3d)-format: x, (y,) zTop, zBottom, k_well, flow_rate
      MEASURE_POINT,   // XML expects the 2d(3d)-format: x, y, (z,) rel_error, abs_error, value
      DIFF_MEASURE_POINT,   // XML expects the 2d(3d)-format: x, y, (z,),x2,y2,(z2) rel_error, abs_error, value
      DIPOLE_POINT		// XML expects the 2d(3d)-format: x, y, (z,),x2,y2,(z2) value
    };


    /**********************************************************************************
     *
     *  CPointData can represent
     *  1.) a mesauring point (MEASURE_POINT)
     *  or
     *  2.) a point with pumping data (SOURCE_POINT)!
     *
     **********************************************************************************
     */
    class CPointData{
    public:
      REAL x,x2;
      REAL y,y2;
#ifdef DIMENSION3
      REAL z,z2;
#endif
      REAL value;

      REAL well_conductivity;
      REAL well_rate;
      REAL well_top;
      REAL well_bottom;
      REAL temperature;
      REAL temperature_injection_time;
      REAL concentration;
      REAL concentration_injection_time;

      REAL rel_error;
      REAL abs_error;
      MEASURE_TYPE measure_type;
      /*
      // Copy constructor needed ??
      CPointData( const CPointData & original_list )
      {
      x=111.0;
      }
      */
      CPointData()
      {
        x=0;
        x2=0;
        y=0;
        y=2;
#ifdef DIMENSION3
        z=0;
        z2=0;
#endif
        well_conductivity=0;
        well_rate=0;
        well_top=0;
        well_bottom=0;
        temperature=0.0;
        temperature_injection_time=0.0;
        concentration=0.0;
        concentration_injection_time=0.0;

        value=0;
        rel_error=0;
        abs_error=0;
        measure_type = UNDEF_MEASURE_TYPE;
      };

      void print(){
        logger << "  x=" << x
               << "  y=" << y
#ifdef DIMENSION3
               << "  z=" << z
#endif
               << "  value=" << value
               << "  rel_error=" << rel_error
               << "  abs_error=" << abs_error
               << std::endl;
      };

      friend std::ostream& operator<<(std::ostream& Stream, const CPointData& d){
        Stream << "  x=" << d.x
               << "  y=" << d.y
#ifdef DIMENSION3
               << "  z=" << d.z
#endif
               << "  value=" << d.value
               << "  rel_error=" << d.rel_error
               << "  abs_error=" << d.abs_error;
        return Stream;
      }
    };


    /**********************************************************************************
     *
     *  CListOfDataPoints is a collection of CPointData objects
     *
     **********************************************************************************
     */
    class CListOfDataPoints
    {
    public:
      PointType pointtype;
      UINT total;
      std::vector<CPointData> pointdata_vector; // This must be resized later!
      CListOfDataPoints()
      {
        pointtype=UNDEFINED;
        pointdata_vector.resize(0);
        total=0;  // IMPORTANT!!!
      };

      // Copy constructor needed to create a copy of this list!
      CListOfDataPoints( const CListOfDataPoints & original_list )
      {
        pointtype = original_list.pointtype;
        total = original_list.total;
        pointdata_vector = original_list.pointdata_vector;
      };
    };




    /**********************************************************************************
     *
     *  CInversionParameters is an object map of the XML-node <inversion_parameters>
     *
     **********************************************************************************
     */
    struct CInversionParameters
    {
      std::string L_prior;
      UINT max_iter;
      REAL lim_stepsize;
      REAL weighting_lim;
      REAL dY_lim;
      REAL dL_lim_fac;
      REAL L_accept_confidence_interval;
      REAL disturbance;    // add artifical disturbance to measurement data for test of stability of the inversion scheme
      REAL m1relerror;
    };


    /**********************************************************************************
     *
     *  CParallelFineTuning is an object map of the XML-node <parallel_fine_tuning>
     *
     **********************************************************************************
     */

    struct CParallelFineTuning{
      int maxNComm;
      UINT JQJ_max;
      UINT JQJ_slices;
      CParallelFineTuning():maxNComm(1),JQJ_max(4096),JQJ_slices(1){
      };


    };

    /**********************************************************************************
     *
     *  CCRParameters is an object map of the XML-node <CR_parameters>
     *
     **********************************************************************************
     */
    struct CCRParameters
    {
      UINT total;
      UINT start;
      UINT max_redo;


      //constructor
      CCRParameters():total(0),start(1),max_redo(1){

      };

    };



    /**********************************************************************************
     *
     *  lnk_Inversion_Data is an object map of the XML-node <lnk_inversion_data>
     *
     **********************************************************************************
     */
    struct Clnk_Inversion_Data
    {
    public:
      CListOfDataPoints mplist;  // list of measuring points for the conductivity
    };



    /**********************************************************************************
     *
     *  CHead_Inversion_Data is an object map of the XML-node <Head_Inversion_Data>
     *
     **********************************************************************************
     */
    struct CHead_Inversion_Data
    {
    public:
      //  CListOfDataPoints pdlist;  // list of pumping data for the head_inversion
      //  CListOfDataPoints wdlist;  // list of well data for the head_inversion
      CListOfDataPoints mplist;  // list of measuring points for the head_inversion
      REAL regularization_factor;
      bool bfixedwidth;
    };


    /**********************************************************************************
     *
     *  Cm0_Inversion_Data is an object map of the XML-node <m0_inversion_data>
     *
     **********************************************************************************
     */
    struct Cm0_Inversion_Data
    {
    public:
      CListOfDataPoints mplist;  // list of measuring points for the m0_inversion
      std::vector<REAL> sampling_volume;
      REAL SV_limit;
      std::vector<REAL> SV_max;
      std::vector<REAL> SV_step;
      REAL regularization_factor;
      bool bfixedwidth;
    };



    /**********************************************************************************
     *
     *  Cm1_Inversion_Data is an object map of the XML-node <m0_inversion_data>
     *
     **********************************************************************************
     */
    struct Cm1_Inversion_Data
    {
    public:
      CListOfDataPoints mplist;  // list of measuring points for the m0_inversion
      std::vector<REAL> sampling_volume;
    };

    /**********************************************************************************
     *
     *  CSoluteConcentration_Inversion_Data is an object map of the XML-node <SoluteConcentration_Inversion_Data>
     *
     **********************************************************************************
     */
    struct CSoluteConcentration_Inversion_Data
    {
      //CListOfDataPoints pdlist;  // list of pumping data for the head_inversion
      //CListOfDataPoints wdlist;  // list of well data for the head_inversion_data
      Cm0_Inversion_Data m0_inversion_data;
      Cm1_Inversion_Data m1_inversion_data;
      REAL regularization_factor;
      bool bfixedwidth;
    };

    /**********************************************************************************
     *
     *  CHeat_Mean_Arrival_Time_Inversion_Data is an object map of the XML-node <heat_mean_arrival_time_inversion_data>
     *
     **********************************************************************************
     */
    struct CHeat_Mean_Arrival_Time_Inversion_Data
    {
    public:
      CListOfDataPoints mpM0list;
      CListOfDataPoints mpM1list;
      CListOfDataPoints mplist;
      REAL regularization_factor;
      bool bfixedwidth;
    };

    /**********************************************************************************
     *
     *  CGeoelectrical_Potential_Inversion_Data is an object map of the XML-node <geoelectrical_potential_inversion_data>
     *
     **********************************************************************************
     */
    struct CGeoelectrical_Potential_Inversion_Data
    {
    public:
      //constructor setting nconfig to 0;
      CGeoelectrical_Potential_Inversion_Data():nconfig(0){};

      UINT nconfig;
      std::string kappa_hdf5_file;
      std::string sigma0_hdf5_file;
      std::vector<CListOfDataPoints> inoutlist;
      std::vector<CListOfDataPoints> mpM0list;
      std::vector<CListOfDataPoints> mpM1list;
      std::vector<CListOfDataPoints> mplist;
    };

    /**********************************************************************************
     *
     *  CSetupData is an object map of the whole inputfile XML!
     *
     **********************************************************************************
     */
    class CSetupData {
#ifdef DIMENSION3
      static const UINT dim = 3;
#else
      static const UINT dim = 2;
#endif
    public:
      CBoundaryValueProblem flow_equation;
      CBoundaryValueProblem transport_equation;
      CBoundaryValueProblem heat_transport_equation;
      CBoundaryValueProblem geoelectrical_potential_equation;


      Clnk_Inversion_Data lnk_inversion_data;
      CHead_Inversion_Data head_inversion_data;
      CListOfDataPoints pdlist;  // list of pumping data for the head_inversion
      CListOfDataPoints wdlist;
      CSoluteConcentration_Inversion_Data solute_concentration_inversion_data;
      CHeat_Mean_Arrival_Time_Inversion_Data heat_inversion_data;
      CGeoelectrical_Potential_Inversion_Data geoelectrical_potential_inversion_data;

      int index; // identifier of setup, needed for class MeasurementList

      // constructor
      CSetupData()
      {
        pdlist.pointtype = SOURCE_POINT;
        wdlist.pointtype = WELL_POINT;

        lnk_inversion_data.mplist.pointtype = MEASURE_POINT;
        head_inversion_data.mplist.pointtype = MEASURE_POINT;

        solute_concentration_inversion_data.m0_inversion_data.mplist.pointtype = MEASURE_POINT;
        solute_concentration_inversion_data.m0_inversion_data.sampling_volume.resize(dim,0.0);
        solute_concentration_inversion_data.m0_inversion_data.SV_limit = 0.0;
        solute_concentration_inversion_data.m0_inversion_data.SV_max.resize(dim,0.0);
        solute_concentration_inversion_data.m0_inversion_data.SV_step.resize(dim,0.0);

        solute_concentration_inversion_data.m1_inversion_data.mplist.pointtype = MEASURE_POINT;
        solute_concentration_inversion_data.m1_inversion_data.sampling_volume.resize(dim,0.0);

        heat_inversion_data.mpM0list.pointtype = MEASURE_POINT;
        heat_inversion_data.mpM1list.pointtype = MEASURE_POINT;
        heat_inversion_data.mplist.pointtype = MEASURE_POINT;

        index = 0;

      };

      //destructor
      ~CSetupData() {
      };
    };


    /**********************************************************************************
     *
     *  CInputData is an object map of the whole inputfile XML!
     *
     **********************************************************************************
     */
    class CInputData {

    private:
      const Dune::MPIHelper &helper;

    public:

      typedef CSetupData SDT;
      typedef CDomainData DDT;

      // global verbosity level for std::cout
      // 0=silent (ideally no logging at all),
      // 1=workflow logging (with ISTL logging off),
      // 2=more detailed logging (with ISTL logging off),
      // 3=more detailed logging (with ISTL verbosity level 1)
      // 4=more detailed logging (with ISTL verbosity level 2)
      // 9=debug logging
      int verbosity;

#ifdef DIMENSION3
      static const UINT dim = 3;
#else
      static const UINT dim = 2;
#endif
      CDomainData domain_data;
      CProblemTypes problem_types;
      CPlotOptions plot_options;
      CYFieldProperties yfield_properties;
      CTransportParameters transport_parameters;
      CInversionParameters inversion_parameters;
      CCRParameters CR_parameters;

      std::vector< CSetupData > setups;

      std::vector<      // how many setups?
        std::vector<    // how many wells per setup?
          std::vector<  // how many elements per well?
            std::pair<REAL,Dune::FieldVector<CTYPE,dim> >
            >
          >
        > listOfAllWellCenters;

      CParallelFineTuning parallel_fine_tuning;

      // The constructor
      CInputData( const Dune::MPIHelper &helper_ )
        : helper( helper_ ),verbosity(0)
      {
        listOfAllWellCenters.resize(0);
      };



      ~CInputData() {
      };




      void loglistOfAllWellCenters() const {
        logger << "=== loglistOfAllWellCenters: " << std::endl;
        for(UINT i=0;i<listOfAllWellCenters.size();i++){
          for(UINT j=0;j<listOfAllWellCenters[i].size();j++){
            for(UINT k=0;k<listOfAllWellCenters[i][j].size();k++){
              logger << "setup #" << i;
              logger << ": well #" << j;
              logger << ": element #" << k;
              logger << ": location = " << listOfAllWellCenters[i][j][k].second;
              logger << ", conductivity = " << listOfAllWellCenters[i][j][k].first << std::endl;
            }
          }
        }
      };












      int well2pumpingdata_for_simple_well_model(const CDomainData dd, CListOfDataPoints& well_list, CListOfDataPoints& pumping_list)
      {
        if( (well_list.pointtype!= WELL_POINT) || (pumping_list.pointtype != SOURCE_POINT ) ){
          logger<<"ERROR in: well2pumpingdata_for_simple_well_model; inputfile.hh: wrong list-type"<<std::endl;
          std::cout<<"ERROR in: well2pumpingdata_for_simple_well_model; inputfile.hh: wrong list-type"<<std::endl;
          return -1;
        }

        //loop over all wells and put them intopoint sources
        for(UINT iWell=0; iWell<well_list.pointdata_vector.size(); iWell++){
          int dim=2;
#ifdef DIMENSION3
          dim=3;
#endif

          double lower=well_list.pointdata_vector[iWell].well_bottom;
          double upper=well_list.pointdata_vector[iWell].well_top;

          UINT   l_int=(UINT)std::floor(lower/dd.virtual_gridsizes[dim-1]);
          UINT   u_int=(UINT)std::floor(upper/dd.virtual_gridsizes[dim-1]);


          for(UINT ii=l_int; ii<=u_int; ii++){

            pumping_list.pointdata_vector.push_back(well_list.pointdata_vector[iWell]);
            /*
                        #ifdef DIMENSION3
                        pumping_list.pointdata_vector[pumping_list.total].z=0.5*(pumping_list.pointdata_vector[pumping_list.total].well_top+pumping_list.pointdata_vector[pumping_list.total].well_bottom);
                        double tmp=pumping_list.pointdata_vector[pumping_list.total].z/dd.gridsizes[2];
                        pumping_list.pointdata_vector[pumping_list.total].z= (dd.gridsizes[2]*std::floor(tmp))+ (dd.gridsizes[2]*0.5);

                        #else
                        pumping_list.pointdata_vector[pumping_list.total].y=0.5*(pumping_list.pointdata_vector[pumping_list.total].well_top+pumping_list.pointdata_vector[pumping_list.total].well_bottom);
                        #endif
            */

            double tmp2;
#ifdef DIMENSION3
            pumping_list.pointdata_vector[pumping_list.total].z=(double)(dd.virtual_gridsizes[2]*ii)+(dd.virtual_gridsizes[2]*0.5);
            //double tmp=pumping_list.pointdata_vector[pumping_list.total].z/dd.virtual_gridsizes[2];
            //pumping_list.pointdata_vector[pumping_list.total].z= (dd.virtual_gridsizes[2]*std::floor(tmp))+ (dd.virtual_gridsizes[2]*0.5);
            tmp2=pumping_list.pointdata_vector[pumping_list.total].y/dd.virtual_gridsizes[1];
            pumping_list.pointdata_vector[pumping_list.total].y= (dd.virtual_gridsizes[1]*std::floor(tmp2))+ (dd.virtual_gridsizes[1]*0.5);
#else
            pumping_list.pointdata_vector[pumping_list.total].y=(double)(dd.virtual_gridsizes[1]*ii)+(dd.virtual_gridsizes[1]*0.5);
#endif

            tmp2=pumping_list.pointdata_vector[pumping_list.total].x/dd.virtual_gridsizes[0];
            pumping_list.pointdata_vector[pumping_list.total].x= (dd.virtual_gridsizes[0]*std::floor(tmp2))+ (dd.virtual_gridsizes[0]*0.5);


            pumping_list.pointdata_vector[pumping_list.total].value=pumping_list.pointdata_vector[pumping_list.total].well_rate/( (u_int-l_int)+1 );

            pumping_list.pointdata_vector[pumping_list.total].well_top=0.0;
            pumping_list.pointdata_vector[pumping_list.total].well_bottom=0.0;
            pumping_list.pointdata_vector[pumping_list.total].well_conductivity=0.0;
            pumping_list.pointdata_vector[pumping_list.total].well_rate=0.0;


            if(helper.rank()==0){

              std::cout<<"well #"<<iWell<<" pointsource added at : "<<pumping_list.pointdata_vector[pumping_list.total].x<<", "
                       << pumping_list.pointdata_vector[pumping_list.total].y
#ifdef DIMENSION3
                       <<", "<<pumping_list.pointdata_vector[pumping_list.total].z
#endif
                       <<"."<<std::endl;
            }
            pumping_list.total+=1;
          }
        }
        return well_list.pointdata_vector.size();

      }

      int convert_well_pumping_data(const REAL & dz, CListOfDataPoints& well_list, CListOfDataPoints& pumping_list)
      {
        std::vector<UINT> to_be_converted(0);

        if( (well_list.pointtype!= WELL_POINT) || (pumping_list.pointtype != SOURCE_POINT ) ){
          logger<<"ERROR in: convert_well_pumping_data; inputfile.hh: wrong list-type"<<std::endl;
          std::cout<<"ERROR in: convert_well_pumping_data; inputfile.hh: wrong list-type"<<std::endl;
          return -1;
        }

        for(UINT ii=0; ii< well_list.pointdata_vector.size(); ii++){
          if(std::fabs(well_list.pointdata_vector[ii].well_top-well_list.pointdata_vector[ii].well_bottom)<dz){
            logger<<"found a pumping location in the wells at index: "<<ii<<std::endl;
            to_be_converted.push_back(ii);
          }
        }

        for(UINT ii=0; ii< to_be_converted.size(); ii++){
          UINT last_pumping=pumping_list.total;

          pumping_list.pointdata_vector.push_back(well_list.pointdata_vector[to_be_converted[ii]]);

#ifdef DIMENSION3
          pumping_list.pointdata_vector[last_pumping].z=0.5*(pumping_list.pointdata_vector[last_pumping].well_top+pumping_list.pointdata_vector[last_pumping].well_bottom);
#else
          pumping_list.pointdata_vector[last_pumping].y=0.5*(pumping_list.pointdata_vector[last_pumping].well_top+pumping_list.pointdata_vector[last_pumping].well_bottom);
#endif
          pumping_list.pointdata_vector[last_pumping].value=pumping_list.pointdata_vector[last_pumping].well_rate;

          pumping_list.pointdata_vector[last_pumping].well_top=0.0;
          pumping_list.pointdata_vector[last_pumping].well_bottom=0.0;
          pumping_list.pointdata_vector[last_pumping].well_conductivity=0.0;
          pumping_list.pointdata_vector[last_pumping].well_rate=0.0;

          pumping_list.total+=1;

        }

        // delete from well list in reverse order!
        for(INT ii= to_be_converted.size()-1; ii>=0; ii-=1 ){
          well_list.total-=1;
          well_list.pointdata_vector.erase(well_list.pointdata_vector.begin()+to_be_converted[ii]);
        }
        return to_be_converted.size();
      }



      /** \brief reads measurement data from a string
       *
       *  \param[in]  nTotal  number of measuring points
       *  \param[in]  datastring row-wise listed set of measurement data
       *  \param[out] list Object containing pointdata vector
       *  \returns success(true) or failure(false)
       *  \note
       */
      template<typename LIST>
      bool readDataToList( int nTotal,
                           const std::string datastring,
                           const MEASURE_TYPE measure_type,
                           LIST& list
                           ){

        if(nTotal<1)
          return true;

        REAL x,x2, y,y2, value;
        REAL well_conductivity;
        REAL well_rate;
        REAL well_top;
        REAL well_bottom;
        REAL temperature;
        REAL temperature_injection_time;
        REAL concentration;
        REAL concentration_injection_time;
        REAL abs_error, rel_error;
#ifdef DIMENSION3
        REAL z,z2;
#endif

        std::stringstream instream; // Declare an input string stream
        instream.clear();                       // clear possible previous data
        //instream.str( datastring.c_str() ); // fill string stream with data from XML
        instream << datastring.c_str(); // fill string stream with data from XML

        list.pointdata_vector.resize(nTotal);
        list.total = nTotal;
        int iLineCounter = 0;
        while( !instream.eof() && iLineCounter < nTotal ) {
          x = 0;
          x2 = 0;
          y = 0;
          y2 = 0;
#ifdef DIMENSION3
          z = 0;
          z2 = 0;
#endif

          well_conductivity = 0.0;
          well_rate = 0.0;
          well_top = 0.0;
          well_bottom = 0.0;
          temperature=0.0;
          temperature_injection_time=0.0;
          concentration=0.0;
          concentration_injection_time=0.0;

          value = 0;
          rel_error = 0;
          abs_error = 0;

          // Reading data...
          if( list.pointtype == WELL_POINT ){
            // ...for well location only:
            instream >> x;
#ifdef DIMENSION3
            instream >> y;
#endif
            instream >> well_top;
            instream >> well_bottom;
            instream >> well_conductivity;
            instream >> well_rate;
            instream >> temperature;
            instream >> temperature_injection_time;
            instream >> concentration;
            instream >> concentration_injection_time;
          }
          else if( list.pointtype == MEASURE_POINT ){
            // ...for measure point only:
            instream >> x >> y;
#ifdef DIMENSION3
            instream >> z;
#endif
            instream >> abs_error >> rel_error;
            instream >> value;
          }
          else if( list.pointtype == SOURCE_POINT ){
            instream >> x >> y;
#ifdef DIMENSION3
            instream >> z;
#endif
            instream >> value;
            instream >> temperature;
            instream >> temperature_injection_time;
            instream >> concentration;
            instream >> concentration_injection_time;
          }else if( list.pointtype == DIFF_MEASURE_POINT ){
            // ...for measure point only:
            instream >> x >> y;
#ifdef DIMENSION3
            instream >> z;
#endif
            instream >> x2 >> y2;
#ifdef DIMENSION3
            instream >> z2;
#endif
            instream >> abs_error >> rel_error;
            instream >> value;
          }else if( list.pointtype == DIPOLE_POINT ){
            // ...for measure point only:
            instream >> x >> y;
#ifdef DIMENSION3
            instream >> z;
#endif
            instream >> x2 >> y2;
#ifdef DIMENSION3
            instream >> z2;
#endif
            instream >> value;
          }


          if( list.pointtype == WELL_POINT ){
            logger << "  x=" << x;
#ifdef DIMENSION3
            logger << "  y=" << y;
#endif
            logger << "  well_top= " << well_top;
            logger << "  well_bottom= " << well_bottom;
            logger << "  well_conductivity=" << well_conductivity;
            logger << "  well_rate= " << well_rate;
            logger << "  temperature= " << temperature;
            logger << "  temperature_injection time= " << temperature_injection_time;
            logger << "  concentration= " << concentration;
            logger << "  concentration_injection time= " << concentration_injection_time;
            logger<< std::endl;
          }
          else if( list.pointtype == MEASURE_POINT ){
            logger << "  x=" << x;
            logger << "  y=" << y;
#ifdef DIMENSION3
            logger << "  z=" << z;
#endif
            logger << "  value=" << value
                   << "  rel_error=" << rel_error
                   << "  abs_error=" << abs_error
                   << std::endl;
          }
          else if( list.pointtype == SOURCE_POINT ){
            logger << "  x=" << x;
            logger << "  y=" << y;
#ifdef DIMENSION3
            logger << "  z=" << z;
#endif
            logger << "  value=" << value;
            logger << "  temperature= " << temperature;
            logger << "  temperature_injection time= " << temperature_injection_time;
            logger << "  concentration= " << concentration;
            logger << "  concentration_injection time= " << concentration_injection_time
                   << std::endl;
          }else if( list.pointtype == DIFF_MEASURE_POINT ){
            logger << "  x=" << x;
            logger << "  y=" << y;
#ifdef DIMENSION3
            logger << "  z=" << z;
#endif
            logger << "  x2=" << x2;
            logger << "  y2=" << y2;
#ifdef DIMENSION3
            logger << "  z2=" << z2;
#endif
            logger << "  value=" << value
                   << "  rel_error=" << rel_error
                   << "  abs_error=" << abs_error
                   << std::endl;
          }else if( list.pointtype == DIPOLE_POINT ){
            logger << "  x=" << x;
            logger << "  y=" << y;
#ifdef DIMENSION3
            logger << "  z=" << z;
#endif
            logger << "  x2=" << x2;
            logger << "  y2=" << y2;
#ifdef DIMENSION3
            logger << "  z2=" << z2;
#endif
            logger << "  value=" << value
                   << std::endl;
          }
          /*
            if(std::fabs(std::floor(x+0.5)-x)<GEO_EPSILON)
            x=std::floor(x+0.5);
            if(std::fabs(std::floor(y+0.5)-y)<GEO_EPSILON)
            y=std::floor(y+0.5);
            if(std::fabs(std::floor(x2+0.5)-x2)<GEO_EPSILON)
            x2=std::floor(x2+0.5);
            if(std::fabs(std::floor(y2+0.5)-y2)<GEO_EPSILON)
            y2=std::floor(y2+0.5);
            #ifdef DIMENSION3
            if(std::fabs(std::floor(z+0.5)-z)<GEO_EPSILON)
            z=std::floor(z+0.5);
            if(std::fabs(std::floor(z2+0.5)-z2)<GEO_EPSILON)
            z2=std::floor(z2+0.5);
            #endif
            if(std::fabs(std::floor(well_top+0.5)-well_top)<GEO_EPSILON)
            well_top=std::floor(well_top+0.5);
            if(std::fabs(std::floor(well_bottom+0.5)-well_bottom)<GEO_EPSILON)
            well_bottom=std::floor(well_bottom+0.5);
          */

          list.pointdata_vector[iLineCounter].x = x;
          list.pointdata_vector[iLineCounter].y = y;
#ifdef DIMENSION3
          list.pointdata_vector[iLineCounter].z = z;
#endif
          list.pointdata_vector[iLineCounter].x2 = x2;
          list.pointdata_vector[iLineCounter].y2 = y2;
#ifdef DIMENSION3
          list.pointdata_vector[iLineCounter].z2 = z2;
#endif
          list.pointdata_vector[iLineCounter].value = value;

          list.pointdata_vector[iLineCounter].well_top = well_top;
          list.pointdata_vector[iLineCounter].well_bottom = well_bottom;
          list.pointdata_vector[iLineCounter].well_conductivity = well_conductivity;
          list.pointdata_vector[iLineCounter].well_rate = well_rate;
          list.pointdata_vector[iLineCounter].temperature = temperature;
          list.pointdata_vector[iLineCounter].temperature_injection_time = temperature_injection_time;
          list.pointdata_vector[iLineCounter].concentration = concentration;
          list.pointdata_vector[iLineCounter].concentration_injection_time = concentration_injection_time;


          list.pointdata_vector[iLineCounter].abs_error = abs_error;
          list.pointdata_vector[iLineCounter].rel_error = rel_error;

          list.pointdata_vector[iLineCounter].measure_type = measure_type;

          iLineCounter++;

          char c=instream.peek();
          if(c=='\n'){
            instream.readsome(&c,1);
            c=instream.peek();
            while (c==' '){
              instream.readsome(&c,1);
              c=instream.peek();
            }
          }
        } // end while

        return true;
      }





      /** \brief fill a standard vector with numbers read from a data string
       *
       *  \param[in]  dataline a string containing values separated by white space
       *  \param[in]  defaultvector set of default values in case dataline is empty
       *  \returns vec vector to be filled with numbers from the dataline
       */
      template<typename NUMBER_TYPE>
      void readVectorFromString( const std::string& dataline,
                                 std::vector<NUMBER_TYPE>& vec,
                                 const std::vector<NUMBER_TYPE>& defaultvector=std::vector<NUMBER_TYPE>(0) )
      {
        vec.clear();
        std::istringstream instream; // Declare an input string stream
        instream.clear(); // Reset from possible previous errors.
        instream.str( dataline );

        NUMBER_TYPE zahl;
        while(!instream.eof()) {
          if (instream >> zahl) {
            //std::cout << zahl << std::endl;
            vec.push_back( zahl );
          }
        }
        if(vec.size()==0)
          vec = defaultvector;

      }





      /** \brief reads a list of measurement data from a string
       *
       *  \param[in]  pt property tree object holding the data element
       *  \param[out] list of pointdata
       *  \returns success(true) or failure(false)
       *
       *  \note The data string is supposed to be consisting of at least one line.
       *        Each line has several values which are separated by a white space.
       *        Each line has identical format.
       */
      template<typename PT, typename VECTOR>
      bool readDataListToVectorList( const PT& pt, VECTOR& list, const MEASURE_TYPE measure_type ){
        logger << "Reading ... " << pt.first << std::endl;
        int nTotal = pt.second.get("<xmlattr>.total",0);
        list.total = nTotal;
        if(nTotal<1){
          list.pointdata_vector.resize(0);
          return true;
        }
        std::string datafile = pt.second.get("<xmlattr>.datafile","");
        std::string datastring = "";
        if(datafile=="")
          datastring = pt.second.get("","");
        else{
          char *textbuffer;
          std::filebuf *pFileBuffer;
          std::ifstream filestream ( datafile.c_str() );
          if( filestream.is_open() ){
            pFileBuffer = filestream.rdbuf();
            // get file size using buffer's members
            int filesize=pFileBuffer->pubseekoff (0,std::ios::end,std::ios::in);
            pFileBuffer->pubseekpos (0,std::ios::in);
            // allocate memory to contain file data
            textbuffer=new char[filesize];
            // get file data
            pFileBuffer->sgetn (textbuffer,filesize);
            datastring = std::string( textbuffer );
            delete[] textbuffer;
            filestream.close();
          }
          else
            std::cout << "WARNING: cannot open datafile " << datafile << std::endl;
        }
        if( datastring != "" ){
          readDataToList(nTotal,datastring,measure_type,list);
          return true;
        }
        else
          return false;
      }




      /** \brief reads a list of measurement data from a string
       *
       *  \param[in]  pt Property tree object the XML element <boundary...></boundary>
       *  \param[out] boundary Boundary object
       *  \returns success(true) or failure(false)
       *
       */
      template<typename PT, typename BoundaryObject>
      bool readBoundaryToBoundaryObject( const PT& pt, BoundaryObject& boundary ){

        try {
          boundary.bctype
            = pt.second.get("<xmlattr>.bctype",1); // Dirichlet B.C. is default

          boundary.orientation
            = pt.second.get("<xmlattr>.orientation","west");

          UINT nAvailableStripes = pt.second.get_child("").size();
          if( nAvailableStripes < 1 ){
            std::cerr << "Error: Stripe is missing!" << std::endl;
            return false;
          }
          UINT iStripes = 0;

          BOOST_FOREACH( PT v, pt.second.get_child("") ){
            if( v.first == "stripe" ) {
              CStripe newstripe;
              newstripe.value
                = v.second.get("<xmlattr>.value",0.0);
              newstripe.regularization_factor
                = v.second.get("<xmlattr>.regularization_factor",1.0);
              newstripe.bfixedwidth =
                ("fixed"==v.second.get("<xmlattr>.regularization_type","dynamic")) ? true : false;

              readVectorFromString<REAL>( v.second.get("startpoint","0 0"), newstripe.startpoint );
              readVectorFromString<REAL>( v.second.get("endpoint", "1 1"), newstripe.endpoint );
              boundary.stripes.push_back( newstripe );
              iStripes++;
            }
          }
          return true;
        }

        catch(std::exception const& ex) {
          std::cout << "WARNING in readBoundaryToBoundaryObject(): " << ex.what() << std::endl;
          return false;
        }


      }



      /** \brief reads a list of measurement data from a string
       *
       *  \param[in]  pt Property tree object holding the XML node for an equation
       *  \param[out] equation Equation object
       *  \returns success(true) or failure(false)
       *
       *  \note Only the information common to all types of equations are retrieved here.
       *        Data specific to one equation type (e.g. injection type) are retrieved separately.
       */
      template<typename PT,typename EQ>
      bool readCommonEquationData( const PT& pt, EQ& equation ){

        try{
          equation.initialize( dim );
          equation.bRecycleMatrixHierarchy =
            ("yes" == pt.second.get( "<xmlattr>.recycle_matrix_hierarchy","no")) ? true : false;

          // per default, a 2d rectangle has 4 boudaries sides, a 3d hexahedron has 6
          // UINT nBoundaries = pt.second.get_child( "" ).size();
          UINT iBoundary=0;
          BOOST_FOREACH( PT boundary_node, pt.second.get_child( "" ) ){
            if( boundary_node.first == "boundary" ) {
              if(readBoundaryToBoundaryObject( boundary_node, equation.boundaries[iBoundary] ))
                iBoundary++;
              else
                return false;
            }
          }
          equation.nBoundaries = iBoundary;
          return true;
        }

        catch(std::exception const& ex) {
          std::cout << "WARNING in readCommonEquationData(): " << ex.what() << std::endl;
          return false;
        }

      }




      /** \brief reads a list of measurement data from a string
       *
       *  \param[in]  filename Name of the inputfile (XML format)
       *  \returns success(true) or failure(false)
       *
       *  \note The XML file is read as one single file stream into a text buffer on the root process.
       *        The contents of the textbuffer is then being sent to all other processes.
       *        Afterwards a copy of the textbuffer is moved to a string stream used as input
       *        for the XML parser's input function read_xml().
       */
      bool readInputFileXml(const std::string filename)
      {

        Dune::Timer watch;

        std::istringstream instream; // Declare an input string stream
        // -------------------------------------------------------------------------------
        // Read XML document as a stringstream only on one processor.
        // Then, send this stream to all other processes using a dynamic character buffer.
        // -------------------------------------------------------------------------------
        std::string xmlstring;
        int filesize = 0;

        if( helper.rank() == 0 ) {

          if( verbosity==9 )
            std::cout << "Read XML document as string stream on P0..." << std::endl;

          char *textbuffer;
          std::filebuf *pFileBuffer;
          std::ifstream filestream ( filename.c_str() );
          if( filestream.is_open() ){
            pFileBuffer = filestream.rdbuf();
            // get file size using buffer's members
            filesize=pFileBuffer->pubseekoff (0,std::ios::end,std::ios::in);
            pFileBuffer->pubseekpos (0,std::ios::in);

            // allocate memory to contain file data, reserve one more char for ending null
            textbuffer=new char[filesize+1];

            // get file data
            pFileBuffer->sgetn (textbuffer,filesize);

            // Important: append null character to string makes sure that the string ends here!
            textbuffer[filesize] = '\0';

            xmlstring = std::string( textbuffer );

            filestream.close();

          }
          else {
            std::cout << "Error: Unable to open inputfile " << filename.c_str() << std::endl;
            return false;
          }

          if( helper.size() > 1 ) {
            // Sender:
            // prepare to broadcast data to other processes
            // Send filesize
            MPI_Bcast( &filesize, 1, MPI_INT, 0, helper.getCommunicator());
            // Send text and do not forget the trailing '\0'
            MPI_Bcast( textbuffer, filesize+1, MPI_CHAR, 0, helper.getCommunicator());

            //textbuffer = "";    // A must-do thing before deleting a char pointer!
            delete[] textbuffer;
          }

        }

        else {

          if( verbosity==9 ){
            std::cout << "Broadcast XML stream to process " << helper.rank() << std::endl;
          }

          // Receive length of character buffer:
          MPI_Bcast( &filesize, 1, MPI_INT, 0, helper.getCommunicator());

          char textbuffer[filesize]; // Do not use char* with new and delete here!

          // Receive character buffer with trailing '\0':
          MPI_Bcast( textbuffer, filesize+1, MPI_CHAR, 0, helper.getCommunicator());
          //std::cout.write (textbuffer,filesize);

          xmlstring = std::string( textbuffer );

          // write content to stdout
          //std::cout << "receiver start" << std::endl;
          //std::cout << textbuffer << std::endl;
          //std::cout << "receiver end" << std::endl;

        }



        // Create empty property tree object
        using boost::property_tree::ptree;
        ptree pt;

        //std::stringstream xmlstream;
        //xmlstream << xmlstring;

        std::istringstream xmlstream;
        xmlstream.str( xmlstring.c_str() );

        //std::cout << xmlstring << std::endl;
        //std::cout << xmlstream.str() << std::endl;


        // starting try-block for XML-Parsing
        try{

          //read_xml( xmlstream, pt );
          read_xml( xmlstream, pt, boost::property_tree::xml_parser::no_comments );

          // verbosity level: (This shows how to read an integer attribute of an XML-element.)
          verbosity = pt.get("inputfile.<xmlattr>.verbosity",0);
          logger << "verbosity = " << verbosity << std::endl;


          // "<plot_options>":
          // (This shows how to read a string attribute of an XML-element with a given default value "no".)
          plot_options.vtk_plot_yfield
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_yfield","no")) ? true : false;
          plot_options.vtk_plot_head
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_head","no")) ? true : false;
          plot_options.vtk_plot_q
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_q","no")) ? true : false;
          plot_options.vtk_plot_element_ordering
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_element_ordering","no")) ? true : false;
          plot_options.vtk_plot_m0
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_m0","no")) ? true : false;
          plot_options.vtk_plot_m1
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_m1","no")) ? true : false;
          plot_options.vtk_plot_heat
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_heat","no")) ? true : false;
          plot_options.vtk_plot_wells
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_wells","no")) ? true : false;
          plot_options.vtk_plot_el_potential_field
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_el_potential_field","no")) ? true : false;
          plot_options.vtk_plot_adjoint_head
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_adjoint_head","no")) ? true : false;
          plot_options.vtk_plot_adjoint_m0
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_adjoint_m0","no")) ? true : false;
          plot_options.vtk_plot_adjoint_m1
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_adjoint_m1","no")) ? true : false;
          plot_options.vtk_plot_heat_adjoint
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_heat_adjoint","no")) ? true : false;
          plot_options.vtk_plot_gp_adjoint
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_gp_adjoint","no")) ? true : false;
          plot_options.vtk_plot_sensitivities
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_sensitivities","no")) ? true : false;
          plot_options.vtk_plot_cross_covariances
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_cross_covariances","no")) ? true : false;
          plot_options.vtk_plot_y_smooth
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_y_smooth","no")) ? true : false;
          plot_options.vtk_plot_y_try
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_y_try","no")) ? true : false;
          plot_options.vtk_plot_y_old
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_y_old","no")) ? true : false;
          plot_options.vtk_plot_h_old
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_h_old","no")) ? true : false;
          plot_options.vtk_plot_m0_old
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_m0_old","no")) ? true : false;
          plot_options.vtk_plot_trialfields
            = ("yes" == pt.get("inputfile.plot_options.<xmlattr>.vtk_plot_trialfields","no")) ? true : false;





          // "<problem_types>":

          // (This shows how to read a string attribute of an XML-element with a given default value "no".)
          problem_types.new_YField
            = ("yes" == pt.get("inputfile.problem_types.<xmlattr>.new_YField","no")) ? true : false;
          problem_types.new_Eigenvalues
            = ("yes" == pt.get("inputfile.problem_types.<xmlattr>.new_Eigenvalues","no")) ? true : false;

          // (This shows how to read a string attribute of an XML-element with a given default value "yes".)
          problem_types.synthetic
            = ("yes" == pt.get("inputfile.problem_types.<xmlattr>.synthetic","yes")) ? true : false;

          problem_types.using_existing_Yold
            = ("yes" == pt.get("inputfile.problem_types.<xmlattr>.using_existing_Yold","no")) ? true : false;

          problem_types.refine_estimate
            = pt.get("inputfile.problem_types.<xmlattr>.refine_estimate",0);

          // (This shows how to read a string attribute of an XML-element with a given default value "off".)
          problem_types.head_forward
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.head_forward","off")) ? true : false;
          problem_types.transport_forward
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.transport_forward","off")) ? true : false;
          problem_types.heat_forward
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.heat_forward","off")) ? true : false;
          problem_types.moments_geoeletric_potential_forward
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.moments_geoeletric_potential_forward","off")) ? true : false;
          problem_types.lnK_inversion
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.lnK_inversion","off")) ? true : false;
          problem_types.head_inversion
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.head_inversion","off")) ? true : false;
          problem_types.heat_mean_arrival_time_inversion
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.heat_mean_arrival_time_inversion","off")) ? true : false;
          problem_types.transport_inversion_m0
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.transport_inversion_m0","off")) ? true : false;
          problem_types.transport_inversion_m1
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.transport_inversion_m1","off")) ? true : false;
          problem_types.moments_geoeletric_potential_inversion
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.moments_geoeletric_potential_inversion","off")) ? true : false;
          problem_types.CR
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.conditional_realization","off")) ? true : false;
          problem_types.estimation_variance_only
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.estimation_variance_only","off")) ? true : false;
          problem_types.generate_measurement_data
            = ("on" == pt.get("inputfile.problem_types.<xmlattr>.generate_measurement_data","off")) ? true : false;

          // ================
          // "<domain_data>":
          // ================
          domain_data.initialize(dim);

          readVectorFromString<REAL>( pt.get("inputfile.domain_data.extensions","1   1"), domain_data.extensions );
          readVectorFromString<UINT>( pt.get("inputfile.domain_data.virtual_nCells","1   1"), domain_data.nCells );
          readVectorFromString<UINT>( pt.get("inputfile.domain_data.chunks","1   1"), domain_data.nChunks );

          readVectorFromString<REAL>( pt.get("inputfile.domain_data.lineplot.startpoint","0 0"), domain_data.lineplot.startpoint );
          readVectorFromString<REAL>( pt.get("inputfile.domain_data.lineplot.endpoint","0 0"), domain_data.lineplot.endpoint );

          readVectorFromString<UINT>( pt.get("inputfile.domain_data.partitions","1   1"), domain_data.nPartitions );

#ifdef USE_YASP
          readVectorFromString<UINT>( pt.get("inputfile.domain_data.yasp_grid","1   1"), domain_data.yasp_nCells, domain_data.nCells );
          domain_data.yasp_baselevel =
            pt.get("inputfile.domain_data.yasp_grid.<xmlattr>.baselevel",0);
          domain_data.yasp_maxlevel =
            pt.get("inputfile.domain_data.yasp_grid.<xmlattr>.maxlevel",0);
          domain_data.yasp_overlap =
            pt.get("inputfile.domain_data.yasp_grid.<xmlattr>.overlap",1);
          for(UINT i = 0; i < dim; i++){
            UINT nYaspGridCells_i = pow( 2, this->domain_data.yasp_baselevel ) * this->domain_data.yasp_nCells[i];
            this->domain_data.yasp_gridsizes[i]
              = this->domain_data.extensions[i] / (double) nYaspGridCells_i;
          }
          domain_data.baselevel_gridsizes = domain_data.yasp_gridsizes;
#endif


#ifdef USE_ALUGRID
          domain_data.alugrid_globalrefine =
            pt.get("inputfile.domain_data.alu_grid.<xmlattr>.globalrefine",0);
          domain_data.alugrid_baselevel =
            pt.get("inputfile.domain_data.alu_grid.<xmlattr>.baselevel",0);
          domain_data.alugrid_maxsteps =
            pt.get("inputfile.domain_data.alu_grid.<xmlattr>.maxsteps",0);
          domain_data.alugrid_refinementfraction =
            pt.get("inputfile.domain_data.alu_grid.<xmlattr>.refinementfraction",0.1);
          domain_data.alugrid_coarseningfraction =
            pt.get("inputfile.domain_data.alu_grid.<xmlattr>.coarseningfraction",0.0);
          domain_data.alugrid_tolerance =
            pt.get("inputfile.domain_data.alu_grid.<xmlattr>.TOL",0.1);
          domain_data.alugrid_strategy =
            pt.get("inputfile.domain_data.alugrid_grid.<xmlattr>.strategy",1);
          // Compute gridsizes of alugrid along directions parallel to the coordinate axis
          for(UINT i = 0; i < dim; i++){
            UINT nAluGridCells_i = std::pow( 2, this->domain_data.alugrid_globalrefine ) * this->domain_data.nCells[i];
            this->domain_data.alugrid_gridsizes[i]
              = this->domain_data.extensions[i] / (double) nAluGridCells_i;
          }
          domain_data.baselevel_gridsizes = domain_data.alugrid_gridsizes;
#endif

#ifdef USE_UG
          domain_data.ug_heapsize =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.heapsize",1000);
          domain_data.ug_globalrefine =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.globalrefine",0);
          domain_data.ug_baselevel =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.baselevel",0);
          domain_data.ug_maxsteps =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.maxsteps",0);
          domain_data.ug_refinementtype =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.refinementtype","nonconforming");
          domain_data.ug_refinementfraction =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.refinementfraction",0.1);
          domain_data.ug_coarseningfraction =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.coarseningfraction",0.0);
          domain_data.ug_tolerance =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.TOL",0.1);
          domain_data.ug_strategy =
            pt.get("inputfile.domain_data.ug_grid.<xmlattr>.strategy",1);
          for(UINT i = 0; i < dim; i++) {
            // Take care: This here works only for a rectangular domain:
            UINT nUGGridCells_i = std::pow( 2, this->domain_data.ug_globalrefine ) * this->domain_data.nCells[i];
            this->domain_data.ug_gridsizes[i] = this->domain_data.extensions[i] / (double) nUGGridCells_i;
          }
          domain_data.baselevel_gridsizes = domain_data.ug_gridsizes;
#endif

          // Calculate the virtual gridsizes for each dimension:
          for (UINT index = 0; index < dim; index++) {
            domain_data.virtual_gridsizes[index]
              = domain_data.extensions[index] / (CTYPE) domain_data.nCells[index];
          }


          // yfield properties:
          UINT nz = pt.get("inputfile.yfield_properties.<xmlattr>.vertical_zones",1);
          UINT nz_available = pt.get_child("inputfile.yfield_properties").size();
          // Be careful though! The number of children could be more than the real number of available zones.

          logger << "active zones = " << nz << std::endl;
          logger << "available zones = " << nz_available << std::endl;

          std::string well_type =
            pt.get("inputfile.yfield_properties.<xmlattr>.well_type","line");

          if( nz_available < nz ){
            std::cerr << "Missing model for conductivity field!" << std::endl;
            return false;
          }
          yfield_properties.initialize( nz_available );

          // Reading all children of an element
          int count=0;
          BOOST_FOREACH( ptree::value_type &v, pt.get_child("inputfile.yfield_properties") ){
            if( v.first == "zone" ) {
              yfield_properties.zones[count].model = v.second.get("<xmlattr>.model","gaussian");
              yfield_properties.zones[count].bottom = v.second.get("<xmlattr>.bottom",0.0);
              yfield_properties.zones[count].variance = v.second.get("<xmlattr>.variance",1.0);
              yfield_properties.zones[count].beta = v.second.get("<xmlattr>.beta",-6.0);
              yfield_properties.zones[count].qbb_y = v.second.get("<xmlattr>.qbb_y",1.0);
              yfield_properties.zones[count].embedding_factor = v.second.get("<xmlattr>.embedding_factor",2.0);
              yfield_properties.zones[count].porosity = v.second.get("<xmlattr>.porosity",0.3);
              yfield_properties.zones[count].rho = v.second.get("<xmlattr>.rho",2560);
              yfield_properties.zones[count].kappa = v.second.get("<xmlattr>.kappa",0.045);
              yfield_properties.zones[count].sigma0 = v.second.get("<xmlattr>.sigma0",0.03);
              yfield_properties.zones[count].c_s = v.second.get("<xmlattr>.heat_capacity",800);
              yfield_properties.zones[count].lambda_s = v.second.get("<xmlattr>.thermal_conductivity",1.5);

              readVectorFromString<REAL>(v.second.get("correlationlengths","1   1"),
                                         yfield_properties.zones[count].correlation_lambda,
                                         domain_data.extensions);
              count++;
            }
          }
          yfield_properties.nz = count; // This is the real number of available zones!

          yfield_properties.random_seed = pt.get("inputfile.yfield_properties.<xmlattr>.random_seed",0);




          // =========================
          // "<transport_parameters>":
          // =========================
          transport_parameters.a_longitudinal
            = pt.get("inputfile.transport_parameters.<xmlattr>.a_l",1e-1);

          transport_parameters.a_transversal
            = pt.get("inputfile.transport_parameters.<xmlattr>.a_t",1e-2);

          transport_parameters.D_molecular
            = pt.get("inputfile.transport_parameters.<xmlattr>.D_m",1.0);

          // SDFEM parameter:
          transport_parameters.deltaSD_factor
            = pt.get("inputfile.transport_parameters.<xmlattr>.deltaSD_factor",1.0);

          // amplifying factor for D_m during adaptive grid refinement:
          transport_parameters.refined_diffusion
            = pt.get("inputfile.transport_parameters.<xmlattr>.refined_diffusion",1.0);

          transport_parameters.tolerance
            = pt.get("inputfile.transport_parameters.<xmlattr>.tolerance",5.0);

          // DG parameters:
          transport_parameters.dg_method
            = pt.get("inputfile.transport_parameters.<xmlattr>.dg_method","SIPG");

          transport_parameters.dg_weights
            = pt.get("inputfile.transport_parameters.<xmlattr>.dg_weights","on");

          transport_parameters.dg_gamma
            = pt.get("inputfile.transport_parameters.<xmlattr>.dg_gamma",2.0);

          // diffusive L2 Projection parameter:
          transport_parameters.l2_diffusion
            = pt.get("inputfile.transport_parameters.<xmlattr>.l2_diffusion",1.0);

          // diffusive L2 Projection parameter for adjoint problem:
          transport_parameters.l2_diffusion_adjoint
            = pt.get("inputfile.transport_parameters.<xmlattr>.l2_diffusion_adjoint",1.0);

          // linear solver parameters
          transport_parameters.istl_max_iter
            = pt.get("inputfile.transport_parameters.<xmlattr>.istl_max_iter",1000);

          transport_parameters.istl_reduction
            = pt.get("inputfile.transport_parameters.<xmlattr>.istl_reduction",1e-10);

          transport_parameters.istl_atOnceAccu
            = pt.get("inputfile.transport_parameters.<xmlattr>.istl_atOnceAccu","on");

          // point source regularization
          transport_parameters.point_smearing
            = pt.get("inputfile.transport_parameters.<xmlattr>.point_smearing",0.1);



          // ==============================
          // "<heat_transport_parameters>":
          // ==============================
          transport_parameters.heat_a_longitudinal
            = pt.get("inputfile.heat_transport_parameters.<xmlattr>.a_l",1e-1);

          transport_parameters.heat_a_transversal
            = pt.get("inputfile.heat_transport_parameters.<xmlattr>.a_t",1e-2);

          transport_parameters.heat_deltaSD_factor
            = pt.get("inputfile.heat_transport_parameters.<xmlattr>.deltaSD_factor",1.0);

          transport_parameters.heat_rho_w
            = pt.get("inputfile.heat_transport_parameters.<xmlattr>.rho_water",1.0);

          transport_parameters.heat_lambda_w
            = pt.get("inputfile.heat_transport_parameters.<xmlattr>.thermal_conductivity_water",1.0);

          transport_parameters.heat_c_w
            = pt.get("inputfile.heat_transport_parameters.<xmlattr>.heat_capacity_water",1.0);




          /*
           * Reading parameters for the stopping criteria of the line search algorithm and the optimisation algorithm
           * <inversion_parameters>
           *
           */

          inversion_parameters.L_prior
            = pt.get("inputfile.inversion_parameters.<xmlattr>.L_prior","on");

          inversion_parameters.max_iter
            = pt.get("inputfile.inversion_parameters.<xmlattr>.max_iter",10);

          inversion_parameters.lim_stepsize
            = pt.get("inputfile.inversion_parameters.<xmlattr>.lim_stepsize",1.0);

          inversion_parameters.weighting_lim
            = pt.get("inputfile.inversion_parameters.<xmlattr>.weighting_lim",1.0);

          inversion_parameters.dY_lim
            = pt.get("inputfile.inversion_parameters.<xmlattr>.dY_lim",1.0);

          inversion_parameters.dL_lim_fac
            = pt.get("inputfile.inversion_parameters.<xmlattr>.dL_lim_fac",1.0);

          inversion_parameters.L_accept_confidence_interval
            = pt.get("inputfile.inversion_parameters.<xmlattr>.L_accept_confidence_interval",1.0);

          inversion_parameters.disturbance
            = pt.get("inputfile.inversion_parameters.<xmlattr>.disturbance",0.0);

          inversion_parameters.m1relerror
            = pt.get("inputfile.inversion_parameters.<xmlattr>.m1relerror",0.0);


          // =========================
          // "<parallel_fine_tuning>":
          // =========================

          parallel_fine_tuning.maxNComm
            = pt.get("inputfile.parallel_fine_tuning.<xmlattr>.maxNComm",1);

          if(parallel_fine_tuning.maxNComm>helper.size()){
            if( helper.rank()==0)
              std::cout<<"Set maximum possible maxNComm!!! ( "<<helper.size()<<" )"<<std::endl;
            logger<<"Set maximum possible maxNComm!!! ( "<<helper.size()<<" )"<<std::endl;
            parallel_fine_tuning.maxNComm=helper.size();
          }

          parallel_fine_tuning.JQJ_max
            = pt.get("inputfile.parallel_fine_tuning.<xmlattr>.JQJ_max",4096);

          parallel_fine_tuning.JQJ_slices
            = pt.get("inputfile.parallel_fine_tuning.<xmlattr>.JQJ_slices",1);



          // =========================
          // "<CR_parameters>":
          // =========================

          CR_parameters.total
            = pt.get("inputfile.CR_parameters.<xmlattr>.total",1);

          CR_parameters.start
            = pt.get("inputfile.CR_parameters.<xmlattr>.start_counting",1);

          CR_parameters.max_redo
            = pt.get("inputfile.CR_parameters.<xmlattr>.max_iterations",1);


          // =========================
          // Setups
          // =========================

          const int total_setups
            = pt.get("inputfile.Setups.<xmlattr>.total_setups",0);

          if( total_setups < 1 ){
            std::cerr << "ERROR: Number of setups is " << total_setups << std::endl;
          }
          else {

            logger << "-----------------------------------------" << std::endl;
            logger << "Reading Setups! total: " << total_setups << std::endl;
            logger << "-----------------------------------------" << std::endl;

            setups.resize(total_setups);

            int nSetups=0; // setup counter

            // loop over all setups
            BOOST_FOREACH( ptree::value_type &setup_node, pt.get_child("inputfile.Setups") ){

              if( nSetups >= total_setups ) {
                if( helper.rank()==0 )
                  std::cout << "Note: There are more setups in the input file. Here, we are using " << total_setups << std::endl;
                break;
              }

              if( setup_node.first == "Setup" ) {

                // First of all, count the items to be read later...
                int iFlowEquations=0;
                int iTransportEquations=0;
                int iHeatTransportEquations=0;
                int iGeoelectricalPotentialEquation=0;

                setups[nSetups].index = nSetups;

                BOOST_FOREACH( ptree::value_type &setup_child, setup_node.second.get_child("") ){

                  // ====================
                  // "<flow_equation>":
                  // ====================
                  if( setup_child.first == "flow_equation" ) {
                    iFlowEquations++;
                    readCommonEquationData( setup_child, setups[nSetups].flow_equation );
                  }

                  // =======================
                  // "<transport_equation>":
                  // =======================
                  if( setup_child.first == "transport_equation" ) {
                    iTransportEquations++;
                    readCommonEquationData( setup_child, setups[nSetups].transport_equation );
                    setups[nSetups].transport_equation.injection_time
                      = setup_child.second.get("<xmlattr>.injection_time",0.0);
                    //logger << "TPE: injection time = " << setups[nSetups].transport_equation.injection_time << std::endl;
                  }

                  // =============================
                  // "<head_transport_equation>":
                  // =============================
                  if( setup_child.first == "heat_transport_equation" ) {
                    iHeatTransportEquations++;
                    readCommonEquationData( setup_child, setups[nSetups].heat_transport_equation );
                    setups[nSetups].heat_transport_equation.injection_time
                      = setup_child.second.get("<xmlattr>.injection_time",0.0);
                  }

                  // ======================================
                  // "<geoelectrical_potential_equation>":
                  // ======================================
                  if( setup_child.first == "geoelectrical_potential_equation" ) {
                    iGeoelectricalPotentialEquation++;
                    readCommonEquationData( setup_child, setups[nSetups].geoelectrical_potential_equation );
                  }


                  //===============
                  // <pumping_data>
                  //===============
                  if( setup_child.first == "pumping_data" ) {
                    readDataListToVectorList( setup_child, setups[nSetups].pdlist, PUMPING_DATA );
                  }

                  //===============
                  // <well_data>
                  //===============
                  if( setup_child.first == "well_data" ) {
                    readDataListToVectorList( setup_child, setups[nSetups].wdlist, WELL_DATA );
                    //std::cout << "DEBUG: well_rate = " << setups[nSetups].wdlist.pointdata_vector[1].well_rate <<std::endl;

                  }


                  //====================
                  //<lnk_inversion_data>
                  //====================
                  if( setup_child.first == "lnk_inversion_data" ) {
                    BOOST_FOREACH( ptree::value_type measurements_node, setup_child.second.get_child("") ){
                      if( measurements_node.first == "lnk_measurements" ) {
                        readDataListToVectorList( measurements_node, setups[nSetups].lnk_inversion_data.mplist, LN_K_DATA );
                      }
                    }
                    //std::cout << "DEBUG: " << setups[nSetups].lnk_inversion_data.mplist.pointdata_vector[2].y << std::endl;
                  }


                  //======================
                  //<head_inversion_data>
                  //======================
                  if( setup_child.first == "head_inversion_data" ) {

                    setups[nSetups].head_inversion_data.regularization_factor
                      = setup_child.second.get("<xmlattr>.regularization_factor",1.0);

                    setups[nSetups].head_inversion_data.bfixedwidth =
                      ("fixed"==setup_child.second.get("<xmlattr>.regularization_type","dynamic")) ? true : false;

                    BOOST_FOREACH( ptree::value_type measurements_node, setup_child.second.get_child("") ){
                      if( measurements_node.first == "head_measurements" ) {
                        readDataListToVectorList( measurements_node, setups[nSetups].head_inversion_data.mplist, HEAD_DATA );

                      }
                    }
                    //std::cout << "DEBUG: " << setups[nSetups].head_inversion_data.mplist.pointdata_vector[1].value << std::endl;
                  }


                  //======================================
                  // <solute_concentration_inversion_data>
                  //======================================
                  if( setup_child.first == "solute_concentration_inversion_data" ) {



                    //====================
                    // <m0_inversion_data>
                    //====================
                    BOOST_FOREACH( ptree::value_type subtype_node, setup_child.second.get_child("") ){
                      if( subtype_node.first == "m0_inversion_data" ) {


                        readVectorFromString<REAL>( subtype_node.second.get("<xmlattr>.sampling_volume",""),
                                                    setups[nSetups].solute_concentration_inversion_data.m0_inversion_data.sampling_volume );

                        logger << "sampling volume = " << setups[nSetups].solute_concentration_inversion_data.m0_inversion_data.sampling_volume[0] << std::endl;
                        logger << "sampling volume = " << setups[nSetups].solute_concentration_inversion_data.m0_inversion_data.sampling_volume[1] << std::endl;


                        readVectorFromString<REAL>( subtype_node.second.get("<xmlattr>.SV_max",""),
                                                    setups[nSetups].solute_concentration_inversion_data.m0_inversion_data.SV_max );

                        readVectorFromString<REAL>( subtype_node.second.get("<xmlattr>.SV_step",""),
                                                    setups[nSetups].solute_concentration_inversion_data.m0_inversion_data.SV_step );

                        setups[nSetups].solute_concentration_inversion_data.m0_inversion_data.SV_limit
                          = subtype_node.second.get("<xmlattr>.SV_limit",1.0);



                        setups[nSetups].solute_concentration_inversion_data.regularization_factor
                          = setup_child.second.get("<xmlattr>.regularization_factor",1.0);

                        setups[nSetups].solute_concentration_inversion_data.bfixedwidth =
                          ("fixed"==setup_child.second.get("<xmlattr>.regularization_type","dynamic")) ? true : false;

                        BOOST_FOREACH( ptree::value_type measurements_node, subtype_node.second.get_child("") ){
                          //====================
                          // <m0_measurements>
                          //====================
                          if( measurements_node.first == "m0_measurements" ) {
                            readDataListToVectorList( measurements_node, setups[nSetups].solute_concentration_inversion_data.m0_inversion_data.mplist, M0_DATA );
                          }
                        }
                        //std::cout << "DEBUG: " << setups[nSetups].solute_concentration_inversion_data.m0_inversion_data.mplist.pointdata_vector[1].x << std::endl;

                      }
                    }
                    //====================
                    // <m1_inversion_data>
                    //====================
                    BOOST_FOREACH( ptree::value_type subtype_node, setup_child.second.get_child("") ){
                      if( subtype_node.first == "m1_inversion_data" ) {

                        readVectorFromString<REAL>( subtype_node.second.get("<xmlattr>.sampling_volume",""),
                                                    setups[nSetups].solute_concentration_inversion_data.m1_inversion_data.sampling_volume );

                        BOOST_FOREACH( ptree::value_type measurements_node, subtype_node.second.get_child("") ){
                          //====================
                          // <m1_measurements>
                          //====================
                          if( measurements_node.first == "m1_measurements" ) {
                            readDataListToVectorList( measurements_node, setups[nSetups].solute_concentration_inversion_data.m1_inversion_data.mplist, M1_DATA );
                          }
                        }
                      }
                    }
                  }


                  //======================================
                  // <heat_inversion_data>
                  //======================================
                  if( setup_child.first == "heat_inversion_data" ) {

                    setups[nSetups].heat_inversion_data.regularization_factor
                      = setup_child.second.get("<xmlattr>.regularization_factor",1.0);

                    setups[nSetups].heat_inversion_data.bfixedwidth =
                      ("fixed"==setup_child.second.get("<xmlattr>.regularization_type","dynamic")) ? true : false;

                    //====================
                    // <m0_measurements>
                    //====================
                    BOOST_FOREACH( ptree::value_type measurements_node, setup_child.second.get_child("") ){
                      if( measurements_node.first == "m0_measurements" ) {
                        readDataListToVectorList( measurements_node, setups[nSetups].heat_inversion_data.mpM0list, HEAT_M0_DATA );
                      }
                    }
                    //std::cout << "DEBUG: " << setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[1].x << std::endl;

                    //====================
                    // <m1_measurements>
                    //====================
                    BOOST_FOREACH( ptree::value_type measurements_node, setup_child.second.get_child("") ){
                      if( measurements_node.first == "m1_measurements" ) {
                        readDataListToVectorList( measurements_node, setups[nSetups].heat_inversion_data.mpM1list, HEAT_M1_DATA );
                      }
                    }
                    //std::cout << "DEBUG: " << setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector[1].x << std::endl;



                    // Ronnie fragen: What is this good for?

                    logger<<"-----------------------------------------"<<std::endl;
                    logger<<" calculating mean arrival time for temperature measurements!"<<std::endl;
                    if(this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector.size()!=this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector.size()){
                      logger<<std::endl<<"WARNING:  missmatch in m0 and m1 measurements for temperature mean arrival time measurements!!"<<std::endl<<std::endl;
                    }
                    UINT nmeas;
                    if(this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector.size()<=this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector.size()){
                      nmeas=this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector.size();
                      this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector.resize(nmeas);
                    }else{
                      nmeas=this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector.size();
                      this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector.resize(nmeas);
                    }

                    std::vector<UINT> ToBeErased;
                    for(UINT ii=0; ii<nmeas; ii++){
                      if( std::fabs(this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[ii].x-this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector[ii].x)< GEO_EPSILON
                          && std::fabs(this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[ii].y-this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector[ii].y)< GEO_EPSILON
#ifdef DIMENSION3
                          && std::fabs(this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[ii].z-this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector[ii].z)< GEO_EPSILON
#endif
                          ){
                        CPointData tmp_data_point=this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector[ii];

                        if(std::fabs(this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[ii].value)>GEO_EPSILON*0.5)
                          tmp_data_point.value/=this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[ii].value;
                        else{
                          logger<<"measurement of heat m0 too small in measurement #"<<ii<<" (use: "<<(GEO_EPSILON*0.5)<<")!"<<std::endl;
                          tmp_data_point.value/=(GEO_EPSILON*0.5);
                        }

                        this->setups[nSetups].heat_inversion_data.mplist.pointdata_vector.push_back(tmp_data_point);
                      }else{
                        logger<<"drop measurement: missmatch in m0 and m1 location! "<<std::endl;
                        logger<<"m0( "<<this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[ii].x;
                        logger<<", "<<this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[ii].y;
#ifdef DIMENSION3
                        logger<<", "<<this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector[ii].z;
#endif
                        logger<<" )"<<std::endl;
                        logger<<"m1( "<<this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector[ii].x;
                        logger<<", "<<this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector[ii].y;
#ifdef DIMENSION3
                        logger<<", "<<this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector[ii].z;
#endif
                        logger<<" )"<<std::endl;
                        ToBeErased.push_back(ii);
                      }
                    }

                    for(UINT ii=0; ii<ToBeErased.size(); ii++){
                      this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector.erase(this->setups[nSetups].heat_inversion_data.mpM0list.pointdata_vector.begin()+ToBeErased[ii]);
                      this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector.erase(this->setups[nSetups].heat_inversion_data.mpM1list.pointdata_vector.begin()+ToBeErased[ii]);
                    }


                  }

                  // ==============================
                  // <geoelectrical_potential_data>
                  // ==============================

                  if( setup_child.first == "geoelectrical_potential_data" ) {

                    setups[nSetups].geoelectrical_potential_inversion_data.kappa_hdf5_file =
                      setup_child.second.get("<xmlattr>.kappa_hdf5_file","");

                    setups[nSetups].geoelectrical_potential_inversion_data.sigma0_hdf5_file =
                      setup_child.second.get("<xmlattr>.sigma0_hdf5_file","");

                    int nGP=0;
                    int total_GP = setup_child.second.get("<xmlattr>.total_configs",0);

                    if( total_GP > 0 ){

                      //========================================
                      // <geoelectrical_potential_configuration>
                      //========================================
                      // start configuration loop
                      BOOST_FOREACH( ptree::value_type config_node, setup_child.second.get_child("") ){
                        if( config_node.first == "geoelectrical_potential_configuration" ) {

                          if(nGP>=total_GP){
                            break;
                          }
                          logger<<std::endl<<"geoelectrical_potential_configuration #"<<nGP+1<<std::endl;

                          CListOfDataPoints inoutlist_tmp;
                          CListOfDataPoints mpM0list_tmp;
                          CListOfDataPoints mpM1list_tmp;
                          CListOfDataPoints mplist_tmp;

                          inoutlist_tmp.pointtype = DIPOLE_POINT;
                          mpM0list_tmp.pointtype = DIFF_MEASURE_POINT;
                          mpM1list_tmp.pointtype = DIFF_MEASURE_POINT;
                          mplist_tmp.pointtype = DIFF_MEASURE_POINT;

                          // start loop over measurement types
                          BOOST_FOREACH( ptree::value_type measurement_type, config_node.second.get_child("") ){
                            //==================
                            // <inout_location>
                            //==================
                            if( measurement_type.first == "inout_location" ) {
                              readDataListToVectorList( measurement_type, inoutlist_tmp, GP_IN_OUT_DATA );
                              //std::cout << "DEBUG: size = " << inoutlist_tmp.pointdata_vector.size() << std::endl;
                              //std::cout << "DEBUG: x = " << inoutlist_tmp.pointdata_vector[0].x << std::endl;
                            }

                            //==================
                            // <m0_measurements>
                            //==================
                            if( measurement_type.first == "m0_measurements" ) {
                              readDataListToVectorList( measurement_type, mpM0list_tmp, GP_M0_DATA );
                              //std::cout << "DEBUG: size = " << mpM0list_tmp.pointdata_vector.size() << std::endl;
                              //std::cout << "DEBUG: x = " << mpM0list_tmp.pointdata_vector[0].x << std::endl;
                            }

                            //==================
                            // <m1_measurements>
                            //==================
                            if( measurement_type.first == "m1_measurements" ) {
                              readDataListToVectorList( measurement_type, mpM1list_tmp, GP_M1_DATA );
                              //std::cout << "DEBUG: size = " << mpM1list_tmp.pointdata_vector.size() << std::endl;
                              //std::cout << "DEBUG: x = " << mpM1list_tmp.pointdata_vector[0].x << std::endl;
                            }

                          } // end of loop over measurement types

                          // Ask Ronnie: What is this good for?
                          logger<<"-----------------------------------------"<<std::endl;

                          //correction if something directly on the domain boundary (only for the maximal extension)
                          for(UINT ii=0; ii<inoutlist_tmp.pointdata_vector.size(); ii++){

                            if( std::fabs(domain_data.extensions[0]-inoutlist_tmp.pointdata_vector[ii].x)<GEO_EPSILON )
                              inoutlist_tmp.pointdata_vector[ii].x=domain_data.extensions[0]-GEO_EPSILON;

                            if(std::fabs(domain_data.extensions[0]-inoutlist_tmp.pointdata_vector[ii].x2)<GEO_EPSILON )
                              inoutlist_tmp.pointdata_vector[ii].x2=domain_data.extensions[0]-GEO_EPSILON;

                            if(std::fabs(domain_data.extensions[1]-inoutlist_tmp.pointdata_vector[ii].y)<GEO_EPSILON )
                              inoutlist_tmp.pointdata_vector[ii].y=domain_data.extensions[1]-GEO_EPSILON;

                            if(std::fabs(domain_data.extensions[1]-inoutlist_tmp.pointdata_vector[ii].y2)<GEO_EPSILON )
                              inoutlist_tmp.pointdata_vector[ii].y2=domain_data.extensions[1]-GEO_EPSILON;

#ifdef DIMENSION3
                            if(std::fabs(domain_data.extensions[2]-inoutlist_tmp.pointdata_vector[ii].z)<GEO_EPSILON )
                              inoutlist_tmp.pointdata_vector[ii].z=domain_data.extensions[2]-GEO_EPSILON;

                            if(std::fabs(domain_data.extensions[2]-inoutlist_tmp.pointdata_vector[ii].z2)<GEO_EPSILON )
                              inoutlist_tmp.pointdata_vector[ii].z2=domain_data.extensions[2]-GEO_EPSILON;
#endif
                          }


                          logger<<"-----------------------------------------"<<std::endl;
                          logger<<" calculating mean arrival time (of geoelectrical potential)"<<std::endl;
                          if(mpM1list_tmp.pointdata_vector.size()!=mpM0list_tmp.pointdata_vector.size()){
                            logger<<std::endl<<"WARNING:  missmatch in m0 and m1 measurements for mean arrival time measurements in configuration setup #"<<nGP<<" (smaller number will be used for number of measurements)"<<std::endl<<std::endl;
                          }
                          UINT nmeas;
                          if(mpM1list_tmp.pointdata_vector.size()<=mpM0list_tmp.pointdata_vector.size()){
                            nmeas=mpM1list_tmp.pointdata_vector.size();
                            mpM0list_tmp.pointdata_vector.resize(nmeas);
                          }else{
                            nmeas=mpM0list_tmp.pointdata_vector.size();
                            mpM1list_tmp.pointdata_vector.resize(nmeas);
                          }
                          std::vector<UINT> ToBeErased;
                          for(UINT ii=0; ii<nmeas; ii++){
                            if(std::fabs(mpM0list_tmp.pointdata_vector[ii].x-mpM1list_tmp.pointdata_vector[ii].x)< GEO_EPSILON
                               && std::fabs(mpM0list_tmp.pointdata_vector[ii].y-mpM1list_tmp.pointdata_vector[ii].y)< GEO_EPSILON
                               && std::fabs(mpM0list_tmp.pointdata_vector[ii].x2-mpM1list_tmp.pointdata_vector[ii].x2)< GEO_EPSILON
                               && std::fabs(mpM0list_tmp.pointdata_vector[ii].y2-mpM1list_tmp.pointdata_vector[ii].y2)< GEO_EPSILON
#ifdef DIMENSION3
                               && std::fabs(mpM0list_tmp.pointdata_vector[ii].z-mpM1list_tmp.pointdata_vector[ii].z)< GEO_EPSILON
                               && std::fabs(mpM0list_tmp.pointdata_vector[ii].z2-mpM1list_tmp.pointdata_vector[ii].z2)< GEO_EPSILON
#endif
                               ){
                              CPointData tmp_data_point=mpM1list_tmp.pointdata_vector[ii];

                              if(std::fabs(mpM0list_tmp.pointdata_vector[ii].value)>GEO_EPSILON*0.5)
                                tmp_data_point.value/=mpM0list_tmp.pointdata_vector[ii].value;
                              else{
                                logger<<"measurement of m0 too small in measurement #"<<ii<<" (use: "<<(GEO_EPSILON*0.5)<<")!"<<std::endl;
                                tmp_data_point.value/=(GEO_EPSILON*0.5);
                              }

                              mplist_tmp.pointdata_vector.push_back(tmp_data_point);
                            }else{
                              logger<<"drop measurement: missmatch in m0 and m1 location! "<<std::endl;
                              logger<<"m0_1( "<<mpM0list_tmp.pointdata_vector[ii].x;
                              logger<<", "<<mpM0list_tmp.pointdata_vector[ii].y;
#ifdef DIMENSION3
                              logger<<", "<<mpM0list_tmp.pointdata_vector[ii].z;
#endif
                              logger<<" )"<<std::endl;
                              logger<<"m0_2( "<<mpM0list_tmp.pointdata_vector[ii].x2;
                              logger<<", "<<mpM0list_tmp.pointdata_vector[ii].y2;
#ifdef DIMENSION3
                              logger<<", "<<mpM0list_tmp.pointdata_vector[ii].z2;
#endif
                              logger<<" )"<<std::endl;
                              logger<<"m1_1( "<<mpM1list_tmp.pointdata_vector[ii].x;
                              logger<<", "<<mpM1list_tmp.pointdata_vector[ii].y;
#ifdef DIMENSION3
                              logger<<", "<<mpM1list_tmp.pointdata_vector[ii].z;
#endif
                              logger<<" )"<<std::endl;
                              logger<<"m1_2( "<<mpM1list_tmp.pointdata_vector[ii].x2;
                              logger<<", "<<mpM1list_tmp.pointdata_vector[ii].y2;
#ifdef DIMENSION3
                              logger<<", "<<mpM1list_tmp.pointdata_vector[ii].z2;
#endif
                              logger<<" )"<<std::endl;
                              ToBeErased.push_back(ii);
                            }
                          }

                          for(UINT ii=0; ii<ToBeErased.size(); ii++){
                            mpM0list_tmp.pointdata_vector.erase(mpM0list_tmp.pointdata_vector.begin()+ToBeErased[ii]);
                            mpM1list_tmp.pointdata_vector.erase(mpM1list_tmp.pointdata_vector.begin()+ToBeErased[ii]);
                          }


                          this->setups[nSetups].geoelectrical_potential_inversion_data.inoutlist.push_back(inoutlist_tmp);
                          this->setups[nSetups].geoelectrical_potential_inversion_data.mpM0list.push_back(mpM0list_tmp);
                          this->setups[nSetups].geoelectrical_potential_inversion_data.mpM1list.push_back(mpM1list_tmp);
                          this->setups[nSetups].geoelectrical_potential_inversion_data.mplist.push_back(mplist_tmp);
                          nGP+=1;

                        } // end if "geoelectrical_potential_configuration"

                      } // end of configuration loop

                      if(total_GP>nGP){ // warning more HT setups requested then in the input file
                        logger << std::endl
                               << "WARNING: " << total_GP
                               << " geoelectrical potential configurations requested but only "
                               << nGP << " configurations found!!"
                               << std::endl << std::endl;
                      }
                      setups[nSetups].geoelectrical_potential_inversion_data.nconfig = nGP;

                    }
                    else {
                      logger << std::endl
                             << "No geoelectrical potential data in the InputFile-XML."
                             << std::endl << std::endl;
                      setups[nSetups].geoelectrical_potential_inversion_data.nconfig = 0;
                    } // end if ( total_GP > 0 )

                  }

                }

                if( iFlowEquations<1 ){
                  std::cerr << "Missing Flow Equation in Setup number " << nSetups << std::endl;
                }

                nSetups += 1; // increment setup index

              } // end if(setup_node.first == "Setup")

            } // end of loop over all setups

          }

        } catch(std::exception const&  ex) {

          std::cout << "WARNING: " << ex.what() << std::endl;

        }


        std::stringstream jobtitle;
        jobtitle << "readInputFileXml: reading " << filename;
        Dune::Gesis::General::log_elapsed_time( watch.elapsed(),
                                                helper.getCommunicator(),
                                                verbosity,
                                                "IO",
                                                jobtitle.str() );

        return true;

      }; // end of readInputFileXml(const std::string filename)






      void writeToScreen() {

        if( helper.rank()==0 ){

          std::cout << "----------------------" << std::endl;
          std::cout << "Summary of input file:" << std::endl;
          std::cout << "----------------------" << std::endl;


          std::cout << " Physical domain sizes in meters :   ";
          for( UINT index=0; index<dim; index++ )
            std::cout << domain_data.extensions[index] << "   ";
          std::cout << std::endl;
          std::cout << " Initial grid: number of elements:   ";
          for( UINT index=0; index<dim; index++ )
            std::cout << domain_data.nCells[index] << "   ";
          std::cout << std::endl;

#ifdef USE_YASP
          std::cout << " YASP element sizes in meters    :   ";
          for( UINT index=0; index<dim; index++ )
            std::cout << domain_data.yasp_gridsizes[index] << "   ";
          std::cout << std::endl;
#endif

#ifdef USE_ALUGRID
          std::cout << " ALUGRID element sizes in meters :   ";
          for( UINT index=0; index<dim; index++ )
            std::cout << domain_data.alugrid_gridsizes[index] << "   ";
          std::cout << std::endl;
#endif

          std::cout << "Scheidegger Dispersion Parameters: " << std::endl;

          std::cout
            << " a_l = " << transport_parameters.a_longitudinal
            << std::endl;
          std::cout
            << " a_t = " << transport_parameters.a_transversal
            << std::endl;
          std::cout
            << " D_m = " << transport_parameters.D_molecular
            << std::endl;

          std::cout
            << "Total number of setups specified in the inputfile = "
            << setups.size()
            << std::endl;

          std::cout << "For more info see logfiles in the OUTPUT...LOG sub-directory!" << std::endl;
          std::cout << std::endl;
        }



        logger << "----------------------" << std::endl;
        logger << "Summary of input file:" << std::endl;
        logger << "----------------------" << std::endl;

        logger << "extensions = ";
        for (UINT index = 0; index < dim; index++)
          logger << this->domain_data.extensions[index] << "   ";
        logger << std::endl;

        logger << "virtual_nCells = ";
        for (UINT index = 0; index < dim; index++)
          logger << this->domain_data.nCells[index] << "   ";
        logger << std::endl;

        logger << "chunks = ";
        for (UINT index = 0; index < dim; index++)
          logger << this->domain_data.nChunks[index] << "   ";
        logger << std::endl;

        logger << "partitions = ";
        for (UINT index = 0; index < dim; index++)
          logger << this->domain_data.nPartitions[index] << "   ";
        logger << std::endl;

#ifdef USE_YASP
        logger << "yasp_nCells = ";
        for (UINT index = 0; index < dim; index++)
          logger << this->domain_data.yasp_nCells[index] << "   ";
        logger << std::endl;
        logger << "yasp_baselevel = " << this->domain_data.yasp_baselevel << std::endl;
        logger << "yasp_maxlevel = " << this->domain_data.yasp_maxlevel << std::endl;
        logger << "yasp_overlap = " << this->domain_data.yasp_overlap << std::endl;
        logger << std::endl;
#endif //USE_YASP

#ifdef USE_UG
        logger << "ug heapsize = " << this->domain_data.ug_heapsize << std::endl;
        logger << "ug globalrefine = " << this->domain_data.ug_globalrefine << std::endl;
#endif //USE_UG

        logger << "Y field properties with random_seed = "
               << this->yfield_properties.random_seed
               << " "
               << std::endl;
        logger<< "Y field properties with nz = "<<this->yfield_properties.nz<<" Zones: "<<std::endl;
        for(UINT ii=0;ii<this->yfield_properties.nz; ii++){
          logger <<" ZONE #"<<ii+1<<std::endl;
          logger << "model = " << this->yfield_properties.zones[ii].model.c_str() << std::endl;
          logger << "beta = " << this->yfield_properties.zones[ii].beta << std::endl;
          logger << "qbb_y = " << this->yfield_properties.zones[ii].qbb_y << std::endl;
          logger << "variance = " << this->yfield_properties.zones[ii].variance << std::endl;
          logger << "embedding_factor = " << this->yfield_properties.zones[ii].embedding_factor << std::endl;
          logger << "porosity = " << this->yfield_properties.zones[ii].porosity << std::endl;

          logger << "correlation lengths  = " << std::endl;
          for (UINT index = 0; index < dim; index++)
            logger << this->yfield_properties.zones[ii].correlation_lambda[index] << "   ";
          logger << std::endl;
        }

      };


    }; // end of class CInputData

  }

}


#endif	/* DUNE_GESIS_INPUTFILE_HH */
