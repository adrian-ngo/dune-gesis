#ifndef GW_OP_CCFV
#define GW_OP_CCFV

//#include<dune/grid/common/genericreferenceelements.hh>
//#include<dune/grid/common/quadraturerules.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>


#include "../wells/wellposition.hh"

/** a local operator for solving the equation
 *
 *   - \Delta u + a*u = f   in \Omega
 *                  u = g   on \Gamma_D\subseteq\partial\Omega
 *  -\nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * with cell-centered finite volumes on axiparallel, structured grids
 *
 * \tparam B a function indicating the type of boundary condition
 * \tparam G a function for the values of the Dirichlet boundary condition
 */

namespace Dune {
  namespace Gesis {

    template< typename GWP
              , typename IDT
              , typename SDT
              , typename BCType
              , typename BCExt
              , typename SourceType
              >
    class GroundWaterOperatorCCFV :    // implement jacobian evaluation in base classes
#ifdef USE_NUMDIFF
      public Dune::PDELab::NumericalJacobianVolume<GroundWaterOperatorCCFV<GWP,IDT,SDT,BCType,BCExt,SourceType>>,
      public Dune::PDELab::NumericalJacobianSkeleton<GroundWaterOperatorCCFV<GWP,IDT,SDT,BCType,BCExt,SourceType>>,
      public Dune::PDELab::NumericalJacobianBoundary<GroundWaterOperatorCCFV<GWP,IDT,SDT,BCType,BCExt,SourceType>>,
#endif // USE_NUMDIFF
      public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags
    {

    private:

      const GWP& gwp;
      const IDT& inputdata;
      const SDT& setupdata;
      const BCType& bctype;
      const BCExt& g_Dirichlet;

      SourceType& sourceterm;
      const EQ::Mode equationMode;


    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };                             // assemble skeleton term
      enum { doAlphaBoundary  = true };

      enum { doLambdaVolume = true };
      enum { doLambdaSkeleton = true }; // -> void lambda_skeleton(...)
      enum { doLambdaBoundary = true }; // -> void lambda_boundary(...)


      template<typename VCType>
        void set_rhs( const VCType& vc1 )
      {
        sourceterm.set_rhs( vc1 );
      }

      /*
        template<typename UType,typename VCType>
        void set_rhs( const UType& vc1,  // GW: eqn 66: sp_ccfv( m_k )
        const VCType& vc2  // GW: eqn 66: phi0
        )
        {
        sourceterm.set_rhs( vc1, vc2 );
        }
      */
      template<typename VCType>
        void set_rhs( const VCType& vc1, const VCType& vc2 )
      {
        sourceterm.set_rhs( vc1, vc2 );
      }

      template<typename VCType>
        void set_rhs( const VCType& vc1, const VCType& vc2, const VCType& vc3, const VCType& vc4 )
      {
        sourceterm.set_rhs( vc1, vc2, vc3, vc4 );
      }

      template<typename COORDINATES>
        void set_PointSourceLocation( const COORDINATES& x )
      {
        sourceterm.set_PointSourceLocation( x );
      }


      template<typename COORDINATES>
        void set_PointSourceLocation( const COORDINATES& x1, const COORDINATES& x2)
      {
        sourceterm.set_PointSourceLocations( x1,x2 );
      }






      // Constructor used to initialize private members:
      GroundWaterOperatorCCFV( const GWP& gwp_,
                               const IDT& inputdata_,
                               const SDT& setupdata_,
                               const BCType& bctype_,
                               const BCExt& g_Dirichlet_,
                               SourceType& sourceterm_,
                               const EQ::Mode equationMode_=EQ::forward
                               ) :
#ifdef USE_NUMDIFF
        NumericalJacobianVolume<GroundWaterOperatorCCFV<GWP,IDT,SDT,BCType,BCExt,SourceType>>(1e-7),
        NumericalJacobianSkeleton<GroundWaterOperatorCCFV<GWP,IDT,SDT,BCType,BCExt,SourceType>>(1e-7),
        NumericalJacobianBoundary<GroundWaterOperatorCCFV<GWP,IDT,SDT,BCType,BCExt,SourceType>>(1e-7),
#endif // USE_NUMDIFF
        gwp( gwp_ ),
        inputdata(inputdata_),
        setupdata(setupdata_),
        bctype( bctype_ ),
        g_Dirichlet( g_Dirichlet_ ),
        sourceterm( sourceterm_ ),
        equationMode( equationMode_ )
        {}




      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
        void alpha_volume( const EG& eg,
                           const LFSU& lfsu,
                           const X& x,
                           const LFSV& lfsv,
                           R& r
                           ) const
      {
        // domain and range field type
        //typedef typename LFSV::Traits::FiniteElementType::
        //  Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        //const int dim = EG::Geometry::dimension;

        // evaluate reaction term
        //Dune::FieldVector<DF,dim> center = eg.geometry().center();
        RF a = 0.0;
        RF f = 0.0;

        r.accumulate(lfsv,0,(a*x(lfsu,0)-f)*eg.geometry().volume());
      }



#ifndef USE_NUMDIFF
      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
        void jacobian_volume( const EG& eg,
                              const LFSU& lfsu,
                              const X& x,
                              const LFSV& lfsv,
                              M& mat
                              ) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // dimensions
        //const int dim = EG::Geometry::dimension;

        RF a = 0.0;
        //RF f = 0.0;

        // residual contribution
        mat.accumulate(lfsu,0,lfsu,0,a*eg.geometry().volume());
      }
#endif // USE_NUMDIFF




      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
        void alpha_skeleton (const IG& ig,
                             const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                             const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                             R& r_s, R& r_n) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = IG::dimension;

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,dim> distance_vector = ig.inside()->geometry().center();
        Dune::FieldVector<DF,dim> outside_global = ig.outside()->geometry().center();
        distance_vector -= outside_global;
        RF distance = distance_vector.two_norm();
        distance_vector /= distance;
        //std::cout << "distance_vector = " << distance_vector << std::endl;

        // evaluate Y-field tensor at cell center, assume it is constant over elements
        const Dune::FieldVector<DF,dim>& localcenter_inside
          = Dune::ReferenceElements<DF,dim>::general( ig.inside()->type() ).position(0,0);

        const Dune::FieldVector<DF,dim>& localcenter_outside
          = Dune::ReferenceElements<DF,dim>::general( ig.outside()->type() ).position(0,0);

        typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(*(ig.inside()),localcenter_inside));
        typename GWP::Traits::DiffusionTensorType tensor_outside(gwp.DiffusionTensor(*(ig.outside()),localcenter_outside));


        // Evaluate Diffusiontensor dependent on the direction
        Dune::FieldVector<DF,dim> kvector_s;
        tensor_inside.mv(distance_vector,kvector_s);
        double A_s = kvector_s.infinity_norm();
        //std::cout << "A_s = " << A_s << std::endl;

        Dune::FieldVector<DF,dim> kvector_n;
        tensor_outside.mv(distance_vector,kvector_n);
        double A_n = kvector_n.infinity_norm();
        //std::cout << "A_n = " << A_n << std::endl;

        // effective conductivity:

        RF element_volume_s = ig.inside()->geometry().volume();
#ifdef DIMENSION3
        RF element_length_s = std::pow( element_volume_s, 1.0/3.0 );
#else
        RF element_length_s = std::sqrt( element_volume_s );
#endif

        RF element_volume_n = ig.outside()->geometry().volume();
#ifdef DIMENSION3
        RF element_length_n = std::pow( element_volume_n, 1.0/3.0 );
#else
        RF element_length_n = std::sqrt( element_volume_n );
#endif


        //Dune::FieldVector<DF,dim>
        //globalcenter_inside = ig.inside()->geometry().global(localcenter_inside);
        //std::cout << "inside = " << globalcenter_inside << std::endl;

        //Dune::FieldVector<DF,dim>
        //globalcenter_outside = ig.outside()->geometry().global(localcenter_outside);
        //std::cout << "outside = " << globalcenter_outside << std::endl;


        double A_eff = 0.0;
        if( element_length_s - element_length_n<1E-12 )
          A_eff = Dune::Gesis::General::harmonicAverage(A_s,A_n);
        else
          A_eff = Dune::Gesis::General::harmonicAverageWeightedByDistance( A_s,
                                                                           A_n,
                                                                           element_length_s,
                                                                           element_length_n );

        // face geometry
        RF face_volume = ig.geometry().volume();

        // diffusive flux for both sides, with tensor:

        r_s.accumulate(lfsu_s,0,-A_eff*(x_n(lfsu_n,0)-x_s(lfsu_s,0))*face_volume/distance);
      }



      template<typename IG, typename LFSU, typename M>
        void addWellTo_Jacobian_Skeleton( const UINT iWell,
                                          const IG& ig,
                                          const LFSU& lfsu_s,
                                          const LFSU& lfsu_n,
                                          M& mat_ss, M& mat_sn,
                                          M& mat_ns, M& mat_nn
                                          ) const {

        return; // testwise switch off

        // dimensions
        const int dim = IG::dimension;

        // extract well-data:
        REAL well_conductivity = setupdata.wdlist.pointdata_vector[iWell].well_conductivity;

        // Check for elements on both sides of ig if they are penetrated by the well!
        // Only if both sides are penetrated and if the well goes vertically
        // they are considered as connected via a 1d PDE.


        // distance between cell centers in global coordinates
        Dune::FieldVector<CTYPE,dim> inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<CTYPE,dim> outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        REAL distance = inside_global.two_norm();

        typedef Dune::PDELab::ElementGeometry<typename IG::Entity> EG;

        if( w_outside !=
            isElementPenetratedByWell(EG(*(ig.inside())),iWell,inputdata,setupdata)
            &&
            w_outside !=
            isElementPenetratedByWell(EG(*(ig.outside())),iWell,inputdata,setupdata) ){


          REAL c1 = -well_conductivity / distance;
          REAL c2 = well_conductivity / distance;

          mat_ss.accumulate( lfsu_s,0, lfsu_s,0 , -c1 );
          mat_sn.accumulate( lfsu_s,0, lfsu_n,0 , c1 );
          mat_ns.accumulate( lfsu_n,0, lfsu_s,0 , -c2 );
          mat_nn.accumulate( lfsu_n,0, lfsu_n,0 , c2 );

        }
      }



#ifndef USE_NUMDIFF

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
        void jacobian_skeleton (const IG& ig,
                                const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                M& mat_ss, M& mat_sn,
                                M& mat_ns, M& mat_nn) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // dimensions
        const int dim = IG::dimension;

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,dim> distance_vector = ig.inside()->geometry().center();
        Dune::FieldVector<DF,dim> outside_global = ig.outside()->geometry().center();
        distance_vector -= outside_global;
        RF distance = distance_vector.two_norm();
        distance_vector /= distance;

        // evaluate Y-field tensor at cell center, assume it is constant over elements
        const Dune::FieldVector<DF,dim>& localcenter_inside
          = Dune::ReferenceElements<DF,dim>::general( ig.inside()->type() ).position(0,0);

        const Dune::FieldVector<DF,dim>& localcenter_outside
          = Dune::ReferenceElements<DF,dim>::general( ig.outside()->type() ).position(0,0);

        typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(*(ig.inside()),localcenter_inside));
        typename GWP::Traits::DiffusionTensorType tensor_outside(gwp.DiffusionTensor(*(ig.outside()),localcenter_outside));


        // Evaluate Diffusiontensor dependent on the direction
        Dune::FieldVector<DF,dim> kvector_s;
        tensor_inside.mv(distance_vector,kvector_s);
        double A_s = kvector_s.infinity_norm();
        //std::cout << "A_s = " << A_s << std::endl;

        Dune::FieldVector<DF,dim> kvector_n;
        tensor_outside.mv(distance_vector,kvector_n);
        double A_n = kvector_n.infinity_norm();
        //std::cout << "A_n = " << A_n << std::endl;


        // effective conductivity:

        RF element_volume_s = ig.inside()->geometry().volume();
#ifdef DIMENSION3
        RF element_length_s = std::pow( element_volume_s, 1.0/3.0 );
#else
        RF element_length_s = std::sqrt( element_volume_s );
#endif

        RF element_volume_n = ig.outside()->geometry().volume();
#ifdef DIMENSION3
        RF element_length_n = std::pow( element_volume_n, 1.0/3.0 );
#else
        RF element_length_n = std::sqrt( element_volume_n );
#endif


        //Dune::FieldVector<DF,dim>
        //globalcenter_inside = ig.inside()->geometry().global(localcenter_inside);
        //std::cout << "inside = " << globalcenter_inside << std::endl;

        //Dune::FieldVector<DF,dim>
        //globalcenter_outside = ig.outside()->geometry().global(localcenter_outside);
        //std::cout << "outside = " << globalcenter_outside << std::endl;


        double A_eff = 0.0;
        if( element_length_s - element_length_n<1E-12 )
          A_eff = Dune::Gesis::General::harmonicAverage(A_s,A_n);
        else
          A_eff = Dune::Gesis::General::harmonicAverageWeightedByDistance( A_s,
                                                                           A_n,
                                                                           element_length_s,
                                                                           element_length_n );

        // face geometry
        RF face_volume = ig.geometry().volume();
        RF c1 = -A_eff*face_volume/distance;
        RF c2 = A_eff*face_volume/distance;

        mat_ss.accumulate( lfsu_s,0, lfsu_s,0 , -c1 );
        mat_sn.accumulate( lfsu_s,0, lfsu_n,0 , c1 );
        mat_ns.accumulate( lfsu_n,0, lfsu_s,0 , -c2 );
        mat_nn.accumulate( lfsu_n,0, lfsu_n,0 , c2 );


        /*
        // Add inflow or outlow rates for all the wells
        UINT nWells = setupdata.wdlist.total;
        for( UINT iWell=0; iWell<nWells; iWell++ ) {
        addWellTo_Jacobian_Skeleton( iWell,
        ig,
        lfsu_s,
        lfsu_n,
        mat_ss,
        mat_sn,
        mat_ns,
        mat_nn
        );
        }
        */

      }

#endif // USE_NUMDIFF


      // skeleton integral depending on test and ansatz functions
      // Here Dirichlet and Neumann boundary conditions are evaluated
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
        void alpha_boundary (const IG& ig,
                             const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                             R& r_s) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = IG::dimension;

        // face geometry
        Dune::FieldVector<DF,dim> face_center = ig.geometry().center();
        RF face_volume = ig.geometry().volume();

        const Dune::FieldVector<DF,dim-1>& faceLocal
          = Dune::ReferenceElements<DF,dim-1>::general( ig.geometry().type() ).position(0,0);


        // do things depending on boundary condition type
        // evaluate boundary condition type
        if( bctype.isDirichlet( ig, faceLocal ) ){
          // Dirichlet boundary:

          Dune::FieldVector<DF,dim> distance_vector = ig.inside()->geometry().center();
          distance_vector -= face_center;
          RF distance = distance_vector.two_norm();
          distance_vector /= distance;

          const Dune::FieldVector<DF,dim>& localcenter_inside
            = Dune::ReferenceElements<DF,dim>::general( ig.inside()->type() ).position(0,0);

          typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(*(ig.inside()),localcenter_inside));

          Dune::FieldVector<DF,dim> kvector_s;
          tensor_inside.mv(distance_vector,kvector_s);
          double A_s = kvector_s.infinity_norm();

          // evaluate boundary condition function
          typename BCExt::Traits::DomainType x
            = ig.geometryInInside().global( faceLocal );

          typename BCExt::Traits::RangeType g;
          g_Dirichlet.evaluate(*(ig.inside()),x,g);

          //r_s.accumulate(lfsu_s,0,-(g-x_s(lfsu_s,0))*face_volume/distance);
          r_s.accumulate(lfsu_s,0,-A_s*(g-x_s(lfsu_s,0))*face_volume/distance);
          return;

        }
        else{
          // Neumann boundary:
          RF j=0.0;
          r_s.accumulate(lfsu_s,0,j*face_volume);
          return;
        }
      }



#ifndef USE_NUMDIFF
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
        void jacobian_boundary (const IG& ig,
                                const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                M& mat_ss) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = IG::dimension;

        // face geometry
        Dune::FieldVector<DF,dim> face_center = ig.geometry().center();
        RF face_volume = ig.geometry().volume();

        const Dune::FieldVector<DF,dim-1>& faceLocal
          = Dune::ReferenceElements<DF,dim-1>::general( ig.geometry().type() ).position(0,0);


        // do things depending on boundary condition type
        // evaluate boundary condition type
        if( bctype.isDirichlet( ig, faceLocal ) ){
          // Dirichlet boundary:

          Dune::FieldVector<DF,dim> distance_vector = ig.inside()->geometry().center();
          distance_vector -= face_center;
          RF distance = distance_vector.two_norm();
          distance_vector /= distance;

          const Dune::FieldVector<DF,dim>& localcenter_inside
            = Dune::ReferenceElements<DF,dim>::general( ig.inside()->type() ).position(0,0);

          typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(*(ig.inside()),localcenter_inside));

          Dune::FieldVector<DF,dim> kvector_s;
          tensor_inside.mv(distance_vector,kvector_s);
          double A_s = kvector_s.infinity_norm();

          // evaluate boundary condition function
          typename BCExt::Traits::DomainType x
            = ig.geometryInInside().global( faceLocal );

          typename BCExt::Traits::RangeType g;
          g_Dirichlet.evaluate(*(ig.inside()),x,g);

          //r_s.accumulate(lfsu_s,0,-(g-x_s(lfsu_s,0))*face_volume/distance);
          mat_ss.accumulate(lfsu_s,0,lfsu_s,0,A_s*face_volume/distance);
          return;

        }
        else{
          // Neumann boundary:
          //RF j=0.0;
          mat_ss.accumulate(lfsu_s,0,lfsu_s,0,0.0);
          return;
        }

      }
#endif // USE_NUMDIFF













      template<typename EG, typename LFSV, typename R>
        void addWellTo_Lambda_Volume( const UINT iWell,
                                      const EG& eg,
                                      const LFSV& lfsv,
                                      R& residual
                                      ) const {

        // return; // testwise switch off

        // should not be used in the GPE!!!
        if( sourceterm.source_nature == GEOELECTRICAL_POTENTIAL_SOURCE ||
            sourceterm.source_nature == GEP_FUNCTIONAL_SOURCE ){
          return;
        }

        // extract well-data:
        const REAL well_rate = setupdata.wdlist.pointdata_vector[iWell].well_rate;

        if(std::abs(well_rate)<1e-12)
          return; // no contribution

        const int dim = EG::Geometry::dimension;

        REAL well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
        REAL well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;

        if( w_outside !=
            isElementPenetratedByWell(eg,iWell,inputdata,setupdata) ){

          const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1];
          const REAL screening_fraction = dz / ( well_top - well_bottom );

          // Note: the well_rate is given by (-K*gradu)*n,
          // i.e. it contains the well_conductivity implicitly!

          //const REAL j = well_rate * eg.geometry().volume();
          REAL fval = well_rate * screening_fraction;// good, verified with Dirchlet BC and zero well tracer

          // TODO: Maybe to be removed! "refinement_factor" not needed for global refinement:
          //REAL refinement_factor = std::pow(2.0,eg.entity().level()); // dim is no needed here for CCFV!
          //fval *= refinement_factor;

          residual.accumulate( lfsv,0,-fval );
        }
      }



      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
        void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;


        // Add to the residual the contributions of the pumping sources or point source in the given locations if needed

        std::vector<RangeType> shapefunctions( lfsv.size() );

        if( sourceterm.source_nature == PUMPING_SOURCES ){
          sourceterm.evaluate_residual_on_element( eg.entity(),
                                                   lfsv,
                                                   shapefunctions,
                                                   r );
        }

        if( sourceterm.source_nature == POINT_SOURCE || sourceterm.source_nature == GEOELECTRICAL_POINT_SOURCE){

          sourceterm.evaluate_residual_on_element( eg.entity(),
                                                   lfsv,
                                                   shapefunctions,
                                                   r );
        }

        if( equationMode == EQ::forward &&
            sourceterm.source_nature == GEOELECTRICAL_POTENTIAL_SOURCE ){
          sourceterm.evaluate_residual_on_element( eg.entity(),
                                                   lfsv,
                                                   shapefunctions,
                                                   r );
        }


        // Add inflow or outlow rates for all the wells
        UINT nWells = setupdata.wdlist.total;
        for( UINT iWell=0; iWell<nWells; iWell++ ) {
          addWellTo_Lambda_Volume( iWell,
                                   eg,
                                   lfsv,
                                   r
                                   );
        }



      }



      template<typename IG, typename LFSV, typename R>
        void lambda_skeleton (const IG& ig,
                              const LFSV& lfsv_s, const LFSV& lfsv_n,
                              R& r_s, R& r_n) const
      {
        if( ( equationMode == EQ::adjoint &&
              sourceterm.source_nature == FUNCTIONAL_SOURCE ) ||
            sourceterm.source_nature == GEP_FUNCTIONAL_SOURCE )
          sourceterm.evaluate_residual_on_skeleton( ig,
                                                    lfsv_s, lfsv_n,
                                                    r_s, r_n );
      }

      template<typename IG, typename LFSV, typename R>
        void lambda_boundary (const IG& ig,
                              const LFSV& lfsv_s,
                              R& r_s) const
      {
        if( ( equationMode == EQ::adjoint &&
              sourceterm.source_nature == FUNCTIONAL_SOURCE ) ||
            sourceterm.source_nature == GEP_FUNCTIONAL_SOURCE )
          sourceterm.evaluate_residual_on_boundary( ig,
                                                    lfsv_s,
                                                    r_s );
      }


    }; // end of class


  }//namespace PDELab

}// namespace Dune

#endif // GW_OP_CCFV
