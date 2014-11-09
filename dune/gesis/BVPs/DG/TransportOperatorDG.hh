// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef TRANSPORT_OPERATOR_DG_HH
#define TRANSPORT_OPERATOR_DG_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>


#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
//#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

//#include"convectiondiffusionparameter.hh"

#include "dune/gesis/common/eval.hh"
#include "dune/gesis/BVPs/wells/wellposition.hh"


#define HOUSTON_h_F

namespace Dune {
  namespace Gesis {

    struct ConvectionDiffusionDGMethod
    {
      enum Type { NIPG, SIPG };
    };

    struct ConvectionDiffusionDGWeights
    {
      enum Type { weightsOn, weightsOff };
    };

    /** a local operator for solving the convection-diffusion equation with discontinuous Galerkin
     *  
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                              u &=& g \mbox{ on } \partial\Omega_D \\
     *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template< 
      typename TP,
      typename FEM,
      typename SOURCE_TYPE,
      typename IDT,
      typename SDT
      >
    class ConvectionDiffusionDG
      : 
#ifdef USE_NUMDIFF
      public Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE,IDT,SDT> >,
      public Dune::PDELab::NumericalJacobianApplySkeleton<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE,IDT,SDT> >,
      public Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE,IDT,SDT> >,
      public Dune::PDELab::NumericalJacobianVolume<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE,IDT,SDT> >,
      public Dune::PDELab::NumericalJacobianSkeleton<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE,IDT,SDT> >,
      public Dune::PDELab::NumericalJacobianBoundary<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE,IDT,SDT> >,
#endif // USE_NUMDIFF
      public Dune::PDELab::FullSkeletonPattern, 
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {

    private:
      enum { dim = TP::Traits::GridViewType::dimension };
      typedef typename TP::Traits::RangeFieldType Real;
      typedef typename Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
      
      
      const IDT& inputdata;
      const SDT& setupdata;
      const TP& tp;
      SOURCE_TYPE& sourceterm;

      const Passenger::Type passenger;
      const EQ::Mode equationMode;

      ConvectionDiffusionDGMethod::Type method;
      ConvectionDiffusionDGWeights::Type weights;
      Real alpha;
      int intorderadd;
      int quadrature_factor;
      Real theta;
      typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasisType;

      typedef Dune::PDELab::LocalBasisCache<LocalBasisType> Cache;

      // In theory it is possible that one and the same local operator is
      // called first with a finite element of one type and later with a
      // finite element of another type.  Since finite elements of different
      // type will usually produce different results for the same local
      // coordinate they cannot share a cache.  Here we use a vector of caches
      // to allow for different orders of the shape functions, which should be
      // enough to support p-adaptivity.  (Another likely candidate would be
      // differing geometry types, i.e. hybrid meshes.)

      std::vector<Cache> cache;

      double epsilon;


    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };

      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume  = true };



      void reset_rhs()
      {
        sourceterm.reset_rhs();
      }

      template<typename VCType>
      void set_rhs( const VCType& vc_solution )
      {
        sourceterm.set_rhs( vc_solution );
      }
      
      template<typename VCType>
      void set_rhs( const VCType& vc_1,const VCType& vc_2 )
      {
        sourceterm.set_rhs( vc_1,vc_2 );
      }

      template<typename COORDINATES>
      void set_PointSourceLocation( const COORDINATES& x )
      {
        sourceterm.set_PointSourceLocation( x );
      }



      //! constructor: pass parameter object
      ConvectionDiffusionDG(
                            const IDT& inputdata_,
                            const SDT& setupdata_,
                            const TP& tp_,
                            SOURCE_TYPE& sourceterm_,
                            const Passenger::Type passenger_,
                            const EQ::Mode equationMode_,
                            ConvectionDiffusionDGMethod::Type method_=ConvectionDiffusionDGMethod::SIPG, 
                            ConvectionDiffusionDGWeights::Type weights_=ConvectionDiffusionDGWeights::weightsOn,
                            Real alpha_=2.0, 
                            int intorderadd_=0
                             )
        : 
#ifdef USE_NUMDIFF
        Dune::PDELab::NumericalJacobianApplyVolume<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE> >(1.0e-7),
        Dune::PDELab::NumericalJacobianApplySkeleton<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE> >(1.0e-7),
        Dune::PDELab::NumericalJacobianApplyBoundary<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE> >(1.0e-7),
        Dune::PDELab::NumericalJacobianVolume<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE> >(1.0e-7),
        Dune::PDELab::NumericalJacobianSkeleton<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE> >(1.0e-7),
        Dune::PDELab::NumericalJacobianBoundary<ConvectionDiffusionDG<TP,FEM,SOURCE_TYPE> >(1.0e-7),
#endif // USE_NUMDIFF
        inputdata( inputdata_ ),
        setupdata( setupdata_ ),
        tp(tp_), 
        sourceterm( sourceterm_ ),
        passenger( passenger_ ),
        equationMode( equationMode_ ),
        method(method_), 
        weights(weights_),
        alpha(alpha_), 
        intorderadd(intorderadd_), 
        quadrature_factor(2),
        cache(20),
        epsilon( inputdata.transport_parameters.D_molecular )
      {
        theta = 1.0;
        if( method==ConvectionDiffusionDGMethod::SIPG )
          theta = -1.0;
      }







      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        //const int dimw = EG::Geometry::dimensionworld;
        const int order = lfsu.finiteElement().localBasis().order();
        const int intorder = intorderadd + quadrature_factor * order;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            evalFunctionBasis( it->position(), lfsu, cache, phi );

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector< Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            evalGradientBasis( eg, 
                               it->position(),
                               lfsu,
                               cache,
                               gradphi );

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphi[i]);
            
            // evaluate diffusion tensor
            typename TP::Traits::DispersionTensorType A;
            A = tp.DispersionTensor( eg.entity(), it->position(), equationMode );

            // compute K * gradient of u
            Dune::FieldVector<RF,dim> Agradu(0.0);
            A.umv(gradu,Agradu);

            // evaluate velocity field on a qPoint in the interior of an element
            typename TP::Traits::RangeType beta 
              = tp.velocity( eg.entity(),it->position(),equationMode);

            // evaluate reaction term
            typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),it->position());

            // integrate (K grad u - bu)*grad phi_i + a*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsv.size(); i++)
              r.accumulate(lfsv,i,( Agradu*gradphi[i] - u*(beta*gradphi[i]) + c*u*phi[i] )*factor);
          }


      }


#ifndef USE_NUMDIFF



      template<
        typename EG,
        typename LFSU, 
        typename X,
        typename M
        >
      void addWellTo_Jacobian_Volume( const int& iWell, 
                                      const EG& eg, 
                                      const LFSU& lfsu, 
                                      const X& x, 
                                      M& mat,
                                      bool isExtractionWell
                                      ) const
      {

        // return; // testwise 

        isExtractionWell = false;

		// domain and range field type
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;

        REAL well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
        REAL well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;

        REAL well_rate 
          = setupdata.wdlist.pointdata_vector[iWell].well_rate;

        if( equationMode == EQ::adjoint )
          well_rate *= -1;  // adjoint case: injection becomes extraction!

        if( equationMode == EQ::forward &&
            well_rate > 1e-12 )
          return; // inflow is taken care of by addWellTo_Lambda_Volume()
        
        RF concentration = 1;
        RF injectiontime = 1;


        if( w_outside != 
            isElementWithinWellZone( eg, iWell, inputdata, setupdata )
            ) {

          isExtractionWell = true;

#ifdef USE_YASP
          UINT baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_ALUGRID
          UINT baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
#ifdef USE_UG
          UINT baselevel = inputdata.domain_data.ug_baselevel;
#endif


          // dim is important here: one level =  four cells in 2D, or 8 cells in 3D
          REAL baselevel_factor = std::pow(2.0,baselevel);
          REAL level_diff = eg.entity().level() - baselevel;

          //REAL refinement_factor_dim = std::pow(2.0,level_diff*(dim-1)); // dim is important here: one level =  four cells in 2D, or 8 cells in 3D
          REAL refinement_factor = std::pow(2.0,level_diff);
          REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1] / refinement_factor;
          //Dune::FieldVector<REAL,dim> cellcenter = eg.geometry().center();
          //REAL z=cellcenter[dim-1];
          REAL screening_length = well_top - well_bottom;
          REAL screening_fraction = dz / screening_length;

          if( equationMode == EQ::adjoint )
            screening_fraction = 1;

          const int order = lfsu.finiteElement().localBasis().order();
          const int qorder = 2 * order;


          // select quadrature rule
          Dune::GeometryType gt = eg.geometry().type();
          const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);


        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            evalFunctionBasis( it->position(), lfsu, cache, phi );

            //RF u0 = concentration * injectiontime / eg.geometry().volume();
            RF u0 = 1.0;
            // scaling with 1./eg.geometry().volume() is important due to refinement

            RF fval = well_rate * u0 * screening_fraction / eg.geometry().volume();


            //fval *= baselevel_factor; // works for outflow well with tracer coming from Dirichlet BC
            if( eg.entity().level() > (int)baselevel ){
              // fval *= (refinement_factor_dim / baselevel_factor);
              fval /= refinement_factor / baselevel_factor;  // for refinement of gv_tp
              //if( dim>2 )
              //fval /= refinement_factor / baselevel_factor;  // 3D: 1 block refined gives 4 blocks in the horizontal
            }
            //if( eg.entity().level() > baselevel )
            //fval *= (refinement_factor_dim / baselevel_factor);  // for refinement of gv_tp


            // integrate (K grad u - bu)*grad phi_i + a*u*phi_i

            RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++) {
                //mat.accumulate(lfsu,i,lfsu,j,( Agradphi[j]*gradphi[i] - phi[j]*(beta*gradphi[i]) + c*phi[j]*phi[i] )*factor);
                //RF extra_well_contribution = -u0*well_rate*phi[j]*phi[i]*factor *screening_fraction; // *screening_fraction
                //mat.accumulate( lfsu, i, lfsu, j, extra_well_contribution );
                mat.accumulate( lfsu, i, lfsu, j, -fval*phi[j]*phi[i]*factor );
                // addwellto_: jacobian_volume (element integral)
              }
          }

        }

      }


      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            M& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;
        
        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;
        const int order = lfsu.finiteElement().localBasis().order();
        const int intorder = intorderadd + quadrature_factor * order;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // transformation
        Dune::FieldMatrix<DF,dimw,dim> jac;


        //bool elementLiesInsideExtractionWell = false;
        UINT nWells = setupdata.wdlist.total;
        for( UINT iWell=0; iWell<nWells; iWell++ ) {
          bool isExtractionWell = false;
          addWellTo_Jacobian_Volume( iWell, 
                                     eg, 
                                     lfsu, 
                                     x, 
                                     mat,
                                     isExtractionWell );
          //if( isExtractionWell )
          //elementLiesInsideExtractionWell = true;
        }

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions
            std::vector<RangeType> phi(lfsu.size());
            evalFunctionBasis( it->position(), lfsu, cache, phi );

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            evalGradientBasis( eg,
                               it->position(), // expects element local coordinates!
                               lfsu,
                               cache,
                               gradphi );

            // evaluate diffusion tensor
            typename TP::Traits::DispersionTensorType A;
            A = tp.DispersionTensor( eg.entity(), it->position(), equationMode );

            std::vector<Dune::FieldVector<RF,dim> > Agradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++) {
              A.mv(gradphi[i],Agradphi[i]);
            }

            // evaluate velocity field on a qPoint in an interior of an element
            typename TP::Traits::RangeType beta 
              = tp.velocity(eg.entity(),it->position(),equationMode);

            // evaluate reaction term
            typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),it->position());

            // integrate (K grad u - bu)*grad phi_i + a*u*phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type j=0; j<lfsu.size(); j++)
              for (size_type i=0; i<lfsu.size(); i++)
                mat.accumulate(lfsu,i,lfsu,j,( Agradphi[j]*gradphi[i] - phi[j]*(beta*gradphi[i]) + c*phi[j]*phi[i] )*factor);

          }

      }

#endif // not USE_NUMDIFF




      // skeleton integral depending on test and ansatz functions
      // doSkeletonTwoSided = true
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s, 
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        const int intorder = intorderadd + quadrature_factor * std::max(lfsu_s.finiteElement().localBasis().order(),lfsu_n.finiteElement().localBasis().order());
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>& 
          outside_local = Dune::ReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);


        typename TP::Traits::DispersionTensorType A_s(0.0), A_n(0.0);
        A_s = tp.DispersionTensor(*(ig.inside()),inside_local,equationMode);
        A_n = tp.DispersionTensor(*(ig.outside()),outside_local,equationMode);

        // face diameter; this should be revised for anisotropic meshes?
#ifdef HOUSTON_h_F
        RF h_F = std::min( ig.inside()->geometry().volume(),
                        ig.outside()->geometry().volume() ) / ig.geometry().volume(); // Houston!
#else
        DF h_s, h_n;
        DF hmax_s = 0.;
        DF hmax_n = 0.;
        General::element_size(ig.inside()->geometry(),h_s,hmax_s);
        General::element_size(ig.outside()->geometry(),h_n,hmax_n);
        RF h_F = std::min(h_s,h_n);
#endif// HOUSTON_h_F

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule 
          = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // tensor times normal
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        Dune::FieldVector<RF,dim> An_F_n;
        A_n.mv(n_F,An_F_n);

        // compute weights
        RF omega_s;
        RF omega_n;
        RF harmonic_average(0.0);
        if (weights==ConvectionDiffusionDGWeights::weightsOn) {
          RF delta_s = (An_F_s*n_F);
          RF delta_n = (An_F_n*n_F);
          omega_s = delta_n/(delta_s+delta_n+1e-20);
          omega_n = delta_s/(delta_s+delta_n+1e-20);
          harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
        }
        else{
          omega_s = omega_n = 0.5;
          harmonic_average = 1.0;
        }

        // get polynomial degree
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        const int order_n = lfsu_n.finiteElement().localBasis().order();
        int degree = std::max( order_s, order_n );

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);


        // loop over quadrature points and integrate normal flux
        for( typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); 
             it!=rule.end(); ++it) {

          // exact normal
          const Dune::FieldVector<DF,dim> n_F_local 
            = ig.unitOuterNormal(it->position());

          // position of quadrature point in local coordinates of elements 
          Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
          Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());

          // evaluate basis functions on both sides of an intersection

          std::vector<RangeType> phi_s(lfsu_s.size());
          evalFunctionBasis( iplocal_s, lfsu_s, cache, phi_s );

          std::vector<RangeType> phi_n(lfsu_n.size());
          evalFunctionBasis( iplocal_n, lfsu_n, cache, phi_n );

          // evaluate u
          RF u_s=0.0;
          for (size_type i=0; i<lfsu_s.size(); i++)
            u_s += x_s(lfsu_s,i)*phi_s[i];

          RF u_n=0.0;
          for (size_type i=0; i<lfsu_n.size(); i++)
            u_n += x_n(lfsu_n,i)*phi_n[i];

          
          // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
          std::vector< Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
          evalGradientBasis( *(ig.inside()),
                             iplocal_s,
                             lfsu_s,
                             cache,
                             tgradphi_s );

          std::vector< Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
          evalGradientBasis( *(ig.outside()),
                             iplocal_n,
                             lfsu_n,
                             cache,
                             tgradphi_n );

          // compute gradient of u
          Dune::FieldVector<RF,dim> gradu_s(0.0);
          for (size_type i=0; i<lfsu_s.size(); i++)
            gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);

          Dune::FieldVector<RF,dim> gradu_n(0.0);
          for (size_type i=0; i<lfsu_n.size(); i++)
            gradu_n.axpy(x_n(lfsu_n,i),tgradphi_n[i]);

          // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side

          typename TP::Traits::RangeType beta 
            = tp.velocity( *(ig.inside()), iplocal_s, equationMode );

          RF normalflux = beta*n_F_local;
          RF omegaup_s, omegaup_n;
          if (normalflux>=0.0) {
            omegaup_s = 1.0;
            omegaup_n = 0.0;
          }
          else {
            omegaup_s = 0.0;
            omegaup_n = 1.0;
          }

          // integration factor
          RF factor = it->weight() * ig.geometry().integrationElement(it->position());

          // convection term (outflow part u_s AND inflow part u_n)
          RF term1 = (omegaup_s*u_s + omegaup_n*u_n) * normalflux *factor;
          
          for (size_type i=0; i<lfsu_s.size(); i++)
            r_s.accumulate(lfsu_s,i,term1 * phi_s[i]); // o.k. verified with Georgoulis paper

          for (size_type i=0; i<lfsu_n.size(); i++)
            r_n.accumulate(lfsu_n,i,-term1 * phi_n[i]);

          // diffusion term 
          RF term2 = -(omega_s*(An_F_s*gradu_s) + omega_n*(An_F_n*gradu_n)) * factor;
          for (size_type i=0; i<lfsu_s.size(); i++)
            r_s.accumulate(lfsu_s,i,term2 * phi_s[i]); // o.k. verified with Georgoulis paper
          // Adrian's note: 
          // ( - <n*Agradu> * [v] ), remember: <w> = 0.5 w_s + 0.5 w_n, [v] = v_s - v_n
          // phi_s contributes to v_s now, but once we are coming from the other side of the intersection
          // phi_s will automatically contribute to -v_n, 
          // because there, the outer normals will have the opposite sign

          for (size_type i=0; i<lfsu_n.size(); i++)
            r_n.accumulate(lfsu_n,i,-term2 * phi_n[i]);

          // (non-)symmetric IP term
          RF term3 = (u_s-u_n) * factor;
          for (size_type i=0; i<lfsu_s.size(); i++) 
            r_s.accumulate(lfsu_s,i,term3 * theta * omega_s * (An_F_s * tgradphi_s[i]));
          // o.k. verified with Georgoulis paper
          // Adrian's note:
          // (+- <n*Agradv> * [u] )

          for (size_type i=0; i<lfsu_n.size(); i++) 
            r_n.accumulate(lfsu_n,i,term3 * theta * omega_n * (An_F_n*tgradphi_n[i]));

          // standard IP term integral
          RF term4 = penalty_factor * (u_s - u_n) * factor;
          for (size_type i=0; i<lfsu_s.size(); i++)
            r_s.accumulate(lfsu_s,i,term4 * phi_s[i]); // o.k. verified with Georgoulis paper

          // Adrian's note: 
          // ( penalty_factor*[u]*[v] ), remember: [u] = u_s - u_n, [v] = v_s - v_n
          // phi_s contributes to v_s now, but once we are coming from the other side of the intersection
          // phi_s will automatically contribute to -v_n, 
          // because there, [u] will have the opposite sign

          for (size_type i=0; i<lfsu_n.size(); i++)
            r_n.accumulate(lfsu_n,i,-term4 * phi_n[i]);
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
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        const int intorder = intorderadd + quadrature_factor * std::max(lfsu_s.finiteElement().localBasis().order(),lfsu_n.finiteElement().localBasis().order());
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>& 
          outside_local = Dune::ReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);

        typename TP::Traits::DispersionTensorType A_s(0.0), A_n(0.0);
        A_s = tp.DispersionTensor(*(ig.inside()),inside_local,equationMode);
        A_n = tp.DispersionTensor(*(ig.outside()),outside_local,equationMode);

        // face diameter; this should be revised for anisotropic meshes?
#ifdef HOUSTON_h_F
        RF h_F = std::min( ig.inside()->geometry().volume(), 
                           ig.outside()->geometry().volume() ) / ig.geometry().volume(); // Houston!
#else
        DF h_s, h_n;
        DF hmax_s = 0., hmax_n = 0.;
        General::element_size(ig.inside()->geometry(),h_s,hmax_s);
        General::element_size(ig.outside()->geometry(),h_n,hmax_n);
        RF h_F = std::min(h_s,h_n);
#endif

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule 
          = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // tensor times normal
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        Dune::FieldVector<RF,dim> An_F_n;
        A_n.mv(n_F,An_F_n);

        // compute weights
        RF omega_s;
        RF omega_n;
        RF harmonic_average(0.0);
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          {
            RF delta_s = (An_F_s*n_F);
            RF delta_n = (An_F_n*n_F);
            omega_s = delta_n/(delta_s+delta_n+1e-20);
            omega_n = delta_s/(delta_s+delta_n+1e-20);
            harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
          }
        else
          {
            omega_s = omega_n = 0.5;
            harmonic_average = 1.0;
          }

        // get polynomial degree
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        const int order_n = lfsu_n.finiteElement().localBasis().order();
        int degree = std::max( order_s, order_n );

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points and integrate normal flux
        for( typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); 
             it!=rule.end(); ++it )
          {
            // exact normal
            const Dune::FieldVector<DF,dim> n_F_local 
              = ig.unitOuterNormal(it->position());

            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
            Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());

            // evaluate basis functions on both sides of an intersection

            std::vector<RangeType> phi_s(lfsu_s.size());
            evalFunctionBasis( iplocal_s, lfsu_s, cache, phi_s );

            std::vector<RangeType> phi_n(lfsu_n.size());
            evalFunctionBasis( iplocal_n, lfsu_n, cache, phi_n );


            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector< Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            evalGradientBasis( *(ig.inside()),
                               iplocal_s,
                               lfsu_s,
                               cache,
                               tgradphi_s );

            std::vector< Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
            evalGradientBasis( *(ig.outside()),
                               iplocal_n,
                               lfsu_n,
                               cache,
                               tgradphi_n );

            // evaluate velocity field and upwinding, assume H(div) velocity field => may choose any side
            typename TP::Traits::RangeType beta 
              = tp.velocity(*(ig.inside()),iplocal_s,equationMode);

            RF normalflux = beta * n_F_local;
            RF omegaup_s, omegaup_n;
            if (normalflux>=0.0) 
              {
                omegaup_s = 1.0;
                omegaup_n = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_n = 1.0;
              }

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());
            RF ipfactor = penalty_factor * factor;

            // do all terms in the order: I convection, II diffusion, III consistency, IV ip
            for (size_type j=0; j<lfsu_s.size(); j++) {
              RF temp1 = -(An_F_s*tgradphi_s[j])*omega_s*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,omegaup_s * phi_s[j] * normalflux *factor * phi_s[i]);
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,temp1 * phi_s[i]);
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * factor * theta * omega_s * (An_F_s*tgradphi_s[i]));
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * ipfactor * phi_s[i]);
              }
            }

            for (size_type j=0; j<lfsu_n.size(); j++) {
              RF temp1 = -(An_F_n*tgradphi_n[j])*omega_n*factor;
              for (size_type i=0; i<lfsu_s.size(); i++) {
                mat_sn.accumulate(lfsu_s,i,lfsu_n,j,omegaup_n * phi_n[j] * normalflux *factor * phi_s[i]);
                mat_sn.accumulate(lfsu_s,i,lfsu_n,j,temp1 * phi_s[i]);
                mat_sn.accumulate(lfsu_s,i,lfsu_n,j,-phi_n[j] * factor * theta * omega_s * (An_F_s*tgradphi_s[i]));
                mat_sn.accumulate(lfsu_s,i,lfsu_n,j,-phi_n[j] * ipfactor * phi_s[i]);
              }
            }

            for (size_type j=0; j<lfsu_s.size(); j++) {
              RF temp1 = -(An_F_s*tgradphi_s[j])*omega_s*factor;
              for (size_type i=0; i<lfsu_n.size(); i++) {
                mat_ns.accumulate(lfsu_n,i,lfsu_s,j,-omegaup_s * phi_s[j] * normalflux *factor * phi_n[i]);
                mat_ns.accumulate(lfsu_n,i,lfsu_s,j,-temp1 * phi_n[i]);
                mat_ns.accumulate(lfsu_n,i,lfsu_s,j,phi_s[j] * factor * theta * omega_n * (An_F_n*tgradphi_n[i]));
                mat_ns.accumulate(lfsu_n,i,lfsu_s,j,-phi_s[j] * ipfactor * phi_n[i]);
              }
            }

            for (size_type j=0; j<lfsu_n.size(); j++) {
              RF temp1 = -(An_F_n*tgradphi_n[j])*omega_n*factor;
              for (size_type i=0; i<lfsu_n.size(); i++) {
                mat_nn.accumulate(lfsu_n,i,lfsu_n,j,-omegaup_n * phi_n[j] * normalflux *factor * phi_n[i]);
                mat_nn.accumulate(lfsu_n,i,lfsu_n,j,-temp1 * phi_n[i]);
                mat_nn.accumulate(lfsu_n,i,lfsu_n,j,-phi_n[j] * factor * theta * omega_n * (An_F_n*tgradphi_n[i]));
                mat_nn.accumulate(lfsu_n,i,lfsu_n,j,phi_n[j] * ipfactor * phi_n[i]);
              }
            }

          }

      }

#endif // USE_NUMDIFF




      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        const int intorder = intorderadd+quadrature_factor*lfsu_s.finiteElement().localBasis().order();
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);

        typename TP::Traits::DispersionTensorType A_s(0.0);
        A_s = tp.DispersionTensor(*(ig.inside()),inside_local,equationMode);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // face diameter
#ifdef HOUSTON_h_F
        RF h_F = ig.inside()->geometry().volume()/ig.geometry().volume(); // Houston!
#else
        DF h_s;
        DF hmax_s = 0.;
        General::element_size(ig.inside()->geometry(),h_s,hmax_s);
        RF h_F = h_s;
#endif
        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // evaluate boundary condition
        const Dune::FieldVector<DF,dim-1> 
          face_local = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        BCType bctype = tp.bctype(ig.intersection(),face_local);

        // compute weights
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        RF harmonic_average;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_average = An_F_s*n_F;
        else
          harmonic_average = 1.0;

        // get polynomial degree
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        int degree = order_s;

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

            // local normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

            // evaluate basis functions

            std::vector<RangeType> phi_s(lfsu_s.size());
            evalFunctionBasis( iplocal_s, lfsu_s, cache, phi_s );

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
              {
                // evaluate flux boundary condition
                RF j = tp.j(ig.intersection(),it->position());
                
                // integrate
                for (size_type i=0; i<lfsv_s.size(); i++) 
                  r_s.accumulate(lfsu_s,i,j * phi_s[i] * factor);

                continue;
              }

            // evaluate u
            RF u_s=0.0; 
            for (size_type i=0; i<lfsu_s.size(); i++) 
              u_s += x_s(lfsu_s,i)*phi_s[i];

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            typename TP::Traits::RangeType beta 
              = tp.velocity(*(ig.inside()),iplocal_s,equationMode);

            RF normalflux = beta * n_F_local;

            if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow) {
              if( normalflux < -1e-10 ) {

                if( inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){                  
                  std::cout << "at " << ig.inside()->geometry().global( iplocal_s ) << std::endl;
                  std::cout << "with center " << ig.inside()->geometry().center() << std::endl;
                  std::cout << "beta " << beta << std::endl;
                  if( equationMode == EQ::adjoint )
                    std::cout << "adjoint " << std::endl;
                  std::cout << "normalflux " << normalflux << ", but bctype is outflow!!!" << std::endl;
                  //DUNE_THROW( Dune::Exception, "alpha_boundary: Outflow boundary condition on inflow!" );
                  std::cout << "WARNING: alpha_boundary: bctype=outflow but we have inflow! Set inflow=0." << std::endl;
                  //exit(99);
                }
                
                normalflux = 0.0; // behaves as if gDirichlet = 0 was given
              }
              // convection term
              RF term1 = u_s * normalflux *factor;
              for (size_type i=0; i<lfsu_s.size(); i++)
                r_s.accumulate(lfsu_s,i,term1 * phi_s[i]);

              // evaluate flux boundary condition
              RF o = tp.o(ig.intersection(),it->position());

              // integrate
              for (size_type i=0; i<lfsv_s.size(); i++)
                r_s.accumulate(lfsu_s,i,o * phi_s[i] * factor);

              continue;
            }


            // Additional check to prevent later problems during adjoint case:
            // This check should be done during the solution of the forward problem.
            // It is enough to have it here in jacobian_boundary
            if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet) {
              if( normalflux > 1e-10 ) {

                if( inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
                  std::cout << "at " << ig.inside()->geometry().global( iplocal_s ) << std::endl;
                  std::cout << "with center " << ig.inside()->geometry().center() << std::endl;
                  std::cout << "beta " << beta << std::endl;
                  if( equationMode == EQ::forward )
                    std::cout << "forward " << std::endl;
                  std::cout << "normalflux " << normalflux << ", but bctype is inflow!!!" << std::endl;
                  //DUNE_THROW( Dune::Exception, "jacobian_boundary: Outflow boundary condition on inflow!" );
                  std::cout << "WARNING: jacobian_boundary: bctype=inflow but we have outflow. Treat as outflow." << std::endl;
                  //exit(99);
                  //normalflux = 0.0;
                }

                // test this from above:
                // convection term
                RF term1 = u_s * normalflux *factor;
                for (size_type i=0; i<lfsu_s.size(); i++)
                  r_s.accumulate(lfsu_s,i,term1 * phi_s[i]);

                continue;
              }
            }



            

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector< Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            evalGradientBasis( *(ig.inside()),
                               iplocal_s,
                               lfsu_s,
                               cache,
                               tgradphi_s );

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu_s(0.0);
            for (size_type i=0; i<lfsu_s.size(); i++)
              gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);

            // evaluate Dirichlet boundary condition
            RF gDirichlet = tp.g_function(*(ig.inside()),iplocal_s);

            // upwind
            RF omegaup_s, omegaup_n;
            if (normalflux>=0.0) 
              {
                omegaup_s = 1.0;
                omegaup_n = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_n = 1.0;
              }

            // convection term
            RF term1 = (omegaup_s * u_s + omegaup_n * gDirichlet) * normalflux *factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term1 * phi_s[i]);

            // diffusion term
            RF term2 =  (An_F_s*gradu_s) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,-term2 * phi_s[i]);

            // (non-)symmetric IP term
            RF term3 = ( u_s - gDirichlet ) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term3 * theta * (An_F_s*tgradphi_s[i]));

            // standard IP term
            RF term4 = penalty_factor * ( u_s - gDirichlet ) * factor;
            for (size_type i=0; i<lfsu_s.size(); i++)
              r_s.accumulate(lfsu_s,i,term4 * phi_s[i]);
          }

        
      }


#ifndef USE_NUMDIFF

      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_ss) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;
        
        // dimensions
        const int dim = IG::dimension;
        const int intorder = intorderadd+quadrature_factor*lfsu_s.finiteElement().localBasis().order();
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);

        typename TP::Traits::DispersionTensorType A_s(0.0);
        A_s = tp.DispersionTensor(*(ig.inside()),inside_local,equationMode);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // face diameter
#ifdef HOUSTON_h_F
        RF h_F = ig.inside()->geometry().volume()/ig.geometry().volume(); // Houston!
#else
        DF h_s;
        DF hmax_s = 0.;
        General::element_size(ig.inside()->geometry(),h_s,hmax_s);
        RF h_F = h_s;
#endif
        // transformation
        Dune::FieldMatrix<DF,dim,dim> jac;

        // evaluate boundary condition
        const Dune::FieldVector<DF,dim-1> 
          face_local = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        BCType bctype = tp.bctype(ig.intersection(),face_local);

        // compute weights
        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
        Dune::FieldVector<RF,dim> An_F_s;
        A_s.mv(n_F,An_F_s);
        RF harmonic_average;
        if (weights==ConvectionDiffusionDGWeights::weightsOn)
          harmonic_average = An_F_s*n_F;
        else
          harmonic_average = 1.0;

        // get polynomial degree
        const int order_s = lfsu_s.finiteElement().localBasis().order();
        int degree = order_s;

        // penalty factor
        RF penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);

        // Neumann boundary makes no contribution to boundary
        if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann) return;

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // position of quadrature point in local coordinates of elements 
            Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());

            // local normal
            const Dune::FieldVector<DF,dim> n_F_local = ig.unitOuterNormal(it->position());

            // evaluate basis functions
            std::vector<RangeType> phi_s(lfsu_s.size());
            evalFunctionBasis( iplocal_s, lfsu_s, cache, phi_s );

            // integration factor
            RF factor = it->weight() * ig.geometry().integrationElement(it->position());

            // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
            typename TP::Traits::RangeType beta 
              = tp.velocity(*(ig.inside()),iplocal_s,equationMode);

            RF normalflux = beta * n_F_local;

            if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow) {
              if (normalflux<-1e-10){

                if( inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
                  std::cout << "at " << ig.inside()->geometry().global(iplocal_s) << std::endl;
                  std::cout << "with center " << ig.inside()->geometry().center() << std::endl;
                  std::cout << "beta " << beta << std::endl;
                  if( equationMode == EQ::adjoint )
                    std::cout << "adjoint " << std::endl;
                  std::cout << "normalflux " << normalflux << ", but bctype is outflow!!!" << std::endl;
                  //DUNE_THROW(Dune::Exception,"jacobian_boundary: Outflow boundary condition on inflow!");
                  std::cout << "WARNING: jacobian_boundary: bctype=outflow but we have inflow! Set inflow = 0." << std::endl;
                  //exit(99);
                }

                normalflux = 0.0;// behaves as if gDirichlet = 0 was given
              }
              // convection term
              for (size_type j=0; j<lfsu_s.size(); j++) 
                for (size_type i=0; i<lfsu_s.size(); i++) 
                  mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * normalflux * factor * phi_s[i]);

              continue;
            }


            // Additional check to prevent later problems during adjoint case:
            // This check should be done during the solution of the forward problem.
            // It is enough to have it here in jacobian_boundary
            if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet) {

              if( normalflux > 1e-10 ) {
                if( inputdata.verbosity >= VERBOSITY_DEBUG_LEVEL ){
                  std::cout << "at " << ig.inside()->geometry().global( iplocal_s ) << std::endl;
                  std::cout << "with center " << ig.inside()->geometry().center() << std::endl;
                  std::cout << "beta " << beta << std::endl;
                  if( equationMode == EQ::forward )
                    std::cout << "forward " << std::endl;
                  std::cout << "normalflux " << normalflux << ", but bctype is inflow!!!" << std::endl;
                  //DUNE_THROW( Dune::Exception, "jacobian_boundary: Outflow boundary condition on inflow!" );
                  std::cout << "WARNING: jacobian_boundary: bctype=inflow but we have outflow. Treat as outflow." << std::endl;
                  //exit(99);
                  //normalflux = 0.0;
                }
                // test this from above:
                // convection term
                for (size_type j=0; j<lfsu_s.size(); j++) 
                  for (size_type i=0; i<lfsu_s.size(); i++) 
                    mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * normalflux * factor * phi_s[i]);

                continue;
              }
            }



            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            std::vector< Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
            evalGradientBasis( *(ig.inside()),
                               iplocal_s,
                               lfsu_s,
                               cache,
                               tgradphi_s );

            // upwind
            RF omegaup_s, omegaup_n;
            if (normalflux>=0.0) 
              {
                omegaup_s = 1.0;
                omegaup_n = 0.0;
              }
            else
              {
                omegaup_s = 0.0;
                omegaup_n = 1.0;
              }

            // convection term
            for (size_type j=0; j<lfsu_s.size(); j++) 
              for (size_type i=0; i<lfsu_s.size(); i++) 
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,omegaup_s * phi_s[j] * normalflux * factor * phi_s[i]);

            // diffusion term
            for (size_type j=0; j<lfsu_s.size(); j++) 
              for (size_type i=0; i<lfsu_s.size(); i++) 
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,-(An_F_s*tgradphi_s[j]) * factor * phi_s[i]);

            // (non-)symmetric IP term
            for (size_type j=0; j<lfsu_s.size(); j++) 
              for (size_type i=0; i<lfsu_s.size(); i++) 
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * factor * theta * (An_F_s*tgradphi_s[i]));

            // standard IP term
            for (size_type j=0; j<lfsu_s.size(); j++) 
              for (size_type i=0; i<lfsu_s.size(); i++) 
                mat_ss.accumulate(lfsu_s,i,lfsu_s,j,penalty_factor * phi_s[j] * phi_s[i] * factor);
          }
      }

#endif // USE_NUMDIFF







      template<typename EG, typename LFSV, typename R>
      void addWellTo_Lambda_Volume( const UINT iWell, 
                                    const EG& eg, 
                                    const LFSV& lfsv,
                                    R& residual
                                    ) const {
        
        //return; // testwise

        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSV::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;

        const int dim = EG::Geometry::dimension;
        
        REAL well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
        REAL well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;

        const int order = lfsv.finiteElement().localBasis().order();
        const int qorder = 2 * order;

        REAL well_rate 
          = setupdata.wdlist.pointdata_vector[iWell].well_rate;


        if(well_rate < -1e-12)
          return; // extraction is taken care of by addWellTo_Jacobian_Volume()

        RF concentration = 0;
        if( passenger == Passenger::heat )
          concentration = setupdata.wdlist.pointdata_vector[iWell].temperature;
        else
          concentration = setupdata.wdlist.pointdata_vector[iWell].concentration;

        if( equationMode == EQ::adjoint )
          concentration = 0;  // adjoint case: extraction well becomes injection well with zero concentration/temperature!

        if(concentration<1e-12)
          return; // no inflow contribution


        RF injectiontime = 0;
        if( passenger == Passenger::heat )
          injectiontime = setupdata.wdlist.pointdata_vector[iWell].temperature_injection_time;
        else
          injectiontime = setupdata.wdlist.pointdata_vector[iWell].concentration_injection_time;
        if(injectiontime<1e-12)
          return; // no inflow contribution


        
        if( w_outside != 
            isElementWithinWellZone( eg, iWell, inputdata, setupdata )
            ){

#ifdef USE_YASP
          UINT baselevel = inputdata.domain_data.yasp_baselevel;
#endif
#ifdef USE_ALUGRID
          UINT baselevel = inputdata.domain_data.alugrid_baselevel;
#endif
#ifdef USE_UG
          UINT baselevel = inputdata.domain_data.ug_baselevel;
#endif

          REAL baselevel_factor = std::pow(2.0,baselevel ); // good without *dim
          REAL level_diff = eg.entity().level() - baselevel;

          //REAL refinement_factor_dim = std::pow(2.0,level_diff*(dim-1)); // dim is important here: one level =  four cells in 2D, or 8 cells in 3D
          REAL refinement_factor = std::pow(2.0,level_diff); 
          const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1] / refinement_factor;

          //Dune::FieldVector<REAL,dim> cellcenter = eg.geometry().center();
          //REAL z=cellcenter[dim-1];
          REAL screening_length = well_top - well_bottom;
          REAL screening_fraction = dz / screening_length;
          
          // Note: the well_rate is given by (-K*gradu)*n, 
          // i.e. it contains the well_conductivity implicitly!
          
          //residual.accumulate( lfsv,0,-j );

          // element integral along well elements:
          Dune::GeometryType gt = eg.geometry().type();
          const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

          // loop over quadrature points
          for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {

              //typename TP::Traits::RangeType beta 
              //    = tp.velocity(eg.entity(),it->position(),equationMode);

            // get gradient basis
            std::vector<Dune::FieldVector<RF,dim> > gradphi( lfsv.size() );
            evalGradientBasis( eg,
                               it->position(), // expects element local coordinates!
                               lfsv,
                               cache,
                               gradphi );

            // evaluate basis functions 
            std::vector<RangeType> phi(lfsv.size());
            evalFunctionBasis( it->position(), lfsv, cache, phi );

            RF u0 = concentration * injectiontime;
            //RF u0 = concentration * injectiontime; //*eg.geometry().volume();
            // scaling with 1./eg.geometry().volume() is important due to refinement

            RF factor = it->weight() * eg.geometry().integrationElement(it->position());

            RF fval = well_rate * u0 / eg.geometry().volume() * screening_fraction; // good!


            //fval *= baselevel_factor;
            if( eg.entity().level() > (int)baselevel ){
              //fval /= (refinement_factor_dim / baselevel_factor);

              fval /= refinement_factor / baselevel_factor;  // for refinement of gv_tp
              //if( dim>2 )
              //fval /= refinement_factor / baselevel_factor;  // 3D: 1 block refined gives 4 blocks in the horizontal
            }

            for (size_type i=0; i<lfsv.size(); i++){
              residual.accumulate(lfsv,i,-fval*phi[i]*factor); // active
              // addwellto_: lambda_volume (element integral)
            }
          }

        }

      }




      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSV::Traits::SizeType size_type;
        

        // TODO: Adding Point-sources only! Check implementation!
        if( sourceterm.source_nature == POINT_SOURCE
            || sourceterm.source_nature == SOLUTE_PUMPING_SOURCES
            )
          {
            std::vector<RangeType> shapefunctions( lfsv.size() );
            sourceterm.evaluate_residual_on_element( eg.entity(),
                                                     lfsv,
                                                     shapefunctions,
                                                     r ); // <--- sourceterm adds contribution to the residual
          }


        // TODO: Adding Point-sources only! Check implementation!
        //if(!bAdjoint && sourceterm.source_nature == FUNCTIONAL_SOURCE){
        //  std::vector<RangeType> shapefunctions( lfsv.size() );
        //  UINT flag_source=0;
        //  sourceterm.evaluate_residual_on_element( eg.entity(),
        //                                           lfsv,
        //                                           shapefunctions,
        //                                           r  // <--- sourceterm adds contribution to the residual
        //                                           );
        //}


        // dimensions
        const int dim = EG::Geometry::dimension;
        const int order = lfsv.finiteElement().localBasis().order();
        const int intorder = intorderadd + 2 * order * 1; // <-- wegen der Gaussglockenfunktion

        // forward eqn: (63)
        // adjoint eqn: (150)
        if( sourceterm.source_nature == FUNCTIONAL_SOURCE || 
            sourceterm.source_nature == GEP_FUNCTIONAL_SOURCE ){

          // select quadrature rule
          Dune::GeometryType gt = eg.geometry().type();
          const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

          // loop over quadrature points
          for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
            {
              // evaluate basis functions 
              std::vector<RangeType> phi(lfsv.size());
              evalFunctionBasis( it->position(), lfsv, cache, phi );

              // evaluate right hand side parameter function
              //Real f;
              //f = tp.f(eg.entity(),it->position());
              // evaluate source term
              typename TP::Traits::RangeFieldType fval = 0.0; 
              sourceterm.evaluate_function( eg.entity(),
                                            it->position(),
                                            fval );

              // integrate fval
              RF factor = it->weight() * eg.geometry().integrationElement(it->position());
              for (size_type i=0; i<lfsv.size(); i++)
                r.accumulate(lfsv,i,-fval*phi[i]*factor);
            }

        }


        // Add inflow rates for all injection wells
        UINT nWells = setupdata.wdlist.total;
        for( UINT iWell=0; iWell<nWells; iWell++ ) {
          addWellTo_Lambda_Volume( iWell, 
                                   eg, 
                                   lfsv,
                                   r );
        }

      } // void lambda_volume()

      //! set time in parameter class
      void setTime (double t)
      {
        tp.setTime(t);
      }
    };
  }
}


#endif // TRANSPORT_OPERATOR_DG_HH
