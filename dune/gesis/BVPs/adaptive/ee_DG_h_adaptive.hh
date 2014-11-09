// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_GESIS_DG_RESIDUAL_ESTIMATOR_HH
#define DUNE_GESIS_DG_RESIDUAL_ESTIMATOR_HH

/* 
   Internal edges are visited only once per default.
   This gives a better for performance. (doSkeletonTwoSided = false)
*/

#include <cstddef>
#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
//#include<dune/common/static_assert.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>


#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
//#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
//#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>


#include <dune/gesis/common/eval.hh>

namespace Dune {
  namespace Gesis {

    /** a local operator for residual-based error estimation
     *  
     * A call to residual() of a grid operator space will assemble
     * the quantity \eta_T^2 for each cell. Note that the squares
     * of the cell indicator \eta_T is stored. To compute the global
     * error estimate sum up all values and take the square root.
     *
     * Assumptions and limitations:
     * - Assumes that LFSU is P_1/Q_1 finite element space
     *   and LFSV is a P_0 finite element space (one value per cell).
     * - Convection term is ignored (but reaction term is included)
     *
     * \tparam TP model of ConvectionDiffusionParameterInterface
     */
    template<typename IDT,typename SDT,typename TP>
    class DGResidualEstimator 
      : public Dune::PDELab::LocalOperatorDefaultFlags
    {
      enum { dim = TP::Traits::GridViewType::dimension };
 
      typedef typename TP::Traits::RangeFieldType Real;

      typedef typename Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

    public:
      // pattern assembly flags
      enum { doPatternVolume = false };
      enum { doPatternSkeleton = false };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };


      //! constructor: pass parameter object
      DGResidualEstimator( const IDT& inputdata_,
                           const SDT& setupdata_,
                           const TP& param_, 
                           const Passenger::Type passenger_=Passenger::solute,
                           const EQ::Mode equationMode_=EQ::forward,
                           ConvectionDiffusionDGMethod::Type method_
                           = ConvectionDiffusionDGMethod::NIPG, 
                           ConvectionDiffusionDGWeights::Type weights_
                           = ConvectionDiffusionDGWeights::weightsOff,
                           Real gamma_=2.0
                           )
        : inputdata(inputdata_),
          setupdata(setupdata_),
          param(param_), 
          passenger(passenger_),
          equationMode(equationMode_),
          method(method_), 
          weights(weights_),
          gamma(gamma_) // interior penalty parameter, same as alpha in TransportOperatorDG
      {}










      template<typename EG, typename LFSV, typename FVAL>
      void addWellTo_Alpha_Volume( const UINT iWell, 
                                   const EG& eg, 
                                   const LFSV& lfsv,
                                   FVAL& fval
                                   ) const {

        const int dim = EG::Geometry::dimension;
        
        REAL well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
        REAL well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;
        REAL well_rate 
          = setupdata.wdlist.pointdata_vector[iWell].well_rate;

        if(well_rate < -1e-12)
          return; // extraction is taken care of by addWellTo_Jacobian_Volume()

        REAL concentration = 0;
        if( passenger == Passenger::heat )
          concentration = setupdata.wdlist.pointdata_vector[iWell].temperature;
        else
          concentration = setupdata.wdlist.pointdata_vector[iWell].concentration;

        if( equationMode == EQ::adjoint )
          concentration = 0;  // adjoint case: extraction well becomes injection well with zero concentration/temperature!

        if(concentration<1e-12)
          return; // no inflow contribution


        REAL injectiontime = 0;
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

          REAL refinement_factor_dim = std::pow(2.0,level_diff*(dim-1)); // dim is important here: one level =  four cells in 2D, or 8 cells in 3D
          REAL refinement_factor = std::pow(2.0,level_diff); 
          const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1] / refinement_factor;

          //Dune::FieldVector<REAL,dim> cellcenter = eg.geometry().center();
          REAL screening_length = well_top - well_bottom;
          REAL screening_fraction = dz / screening_length;


          REAL u0 = concentration * injectiontime;
          REAL fval = well_rate * u0 / eg.geometry().volume() * screening_fraction; // good!

          if( eg.entity().level() > baselevel ){
            fval /= (refinement_factor_dim / baselevel_factor);
          }
          // addwellto_: alpha_volume (element integral)

        } // end if inside well zone

      }








      // volume integral depending on test and ansatz functions
      template< typename EG, 
                typename LFSU, 
                typename X, 
                typename LFSV, 
                typename R >
      void alpha_volume( const EG& eg, 
                         const LFSU& lfsu, 
                         const X& x, 
                         const LFSV& lfsv, 
                         R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        
        // dimensions
        const int dim = EG::Geometry::dimension;

        // pOrder is constant on all grid elements (h-adaptive scheme).
        const int pOrder = lfsu.finiteElement().localBasis().order();
        const int intorder = 2 * pOrder;

        // Get geometry of element
        Dune::GeometryType gt = eg.geometry().type();

        // Diffusion tensor at cell center
        typename TP::Traits::DispersionTensorType A;

        //Dune::FieldVector<DF,dim> localcenter = Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);

        A = param.DispersionTensor(eg.entity(),localcenter,equationMode);
        
#ifdef TEST_MAX
        RF epsilon = std::max( A[0][0], A[1][1]);
        if( dim>2 ) epsilon = std::max( A[2][2], epsilon );
#else
        RF epsilon = std::min( A[0][0], A[1][1]);
        if( dim>2 ) epsilon = std::min( A[2][2], epsilon );
#endif

        // select quadrature rule
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
        
        // loop over quadrature points
        RF sum(0.0);
        for( typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin()
               ; it!=rule.end(); ++it ){

          // evaluate u_h
          RF u=0.0;
          evalFunction( it->position(), lfsu, x, u );

          // evaluate reaction term
          typename TP::Traits::RangeFieldType c = param.c(eg.entity(),it->position());


          
          // Add inflow rates for all injection wells
          typename TP::Traits::RangeFieldType fval(0);
          UINT nWells = setupdata.wdlist.total;
          for( UINT iWell=0; iWell<nWells; iWell++ ) {
            addWellTo_Alpha_Volume( iWell,
                                    eg,
                                    lfsv,
                                    fval );
          }



          // evaluate right hand side parameter function (TODO: sourceterm?)
          typename TP::Traits::RangeFieldType f0 = param.f(eg.entity(),it->position());


          // evaluate convection term
          typename TP::Traits::RangeType beta 
            = param.velocity(eg.entity(),it->position(),equationMode);

          // compute gradient of u_h
          Dune::FieldVector<RF,dim> gradu(0.0);
          evalGradient( it->position(), eg.entity(), lfsu, x, gradu );

          // integrate f^2
          RF factor = it->weight() * eg.geometry().integrationElement(it->position());

          RF square = f0 + fval - (beta*gradu) - c*u; // + eps * Laplacian_u (TODO for pMax=2)
          square *= square;
          sum += square * factor;

        }
        
        // accumulate cell indicator 
        DF h_T = diameter(eg.geometry());

        // Liang Zhu: First term, interior residual squared
        RF eta_RK = h_T * h_T / pOrder / pOrder / epsilon * sum;
        r.accumulate( lfsv, 0, eta_RK );

#ifdef USE_SHOW_PECLET
        typename TP::Traits::RangeType beta 
          = param.velocity(eg.entity(),localcenter,equationMode);
        RF mesh_peclet = 0.5 * beta.two_norm() * h_T / epsilon;
        r.accumulate( lfsv, 1, mesh_peclet ); 
        /* 
           If mesh_peclet goes down to the value 2.0, 
           the problem should become diffusive. 
           In that case, there should be no more oscillations and we are done!
           Otherwise, the number of iterations required for the iterative solver
           will become proportional to the problem size.
        */
#endif

      }


      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG,
               typename LFSU,
               typename X,
               typename LFSV,
               typename R>
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
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = IG::dimension;
        
        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>& 
          outside_local = Dune::ReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);

        typename TP::Traits::DispersionTensorType A_s, A_n;
        A_s = param.DispersionTensor(*(ig.inside()),inside_local,equationMode);
        A_n = param.DispersionTensor(*(ig.outside()),outside_local,equationMode);

#ifdef TEST_MAX
        RF epsilon_s = std::max( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::max( A_s[2][2], epsilon_s );

        RF epsilon_n = std::max( A_n[0][0], A_n[1][1]);
        if( dim>2 ) epsilon_n = std::max( A_n[2][2], epsilon_n );
#else
        RF epsilon_s = std::min( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::min( A_s[2][2], epsilon_s );

        RF epsilon_n = std::min( A_n[0][0], A_n[1][1]);
        if( dim>2 ) epsilon_n = std::min( A_n[2][2], epsilon_n );
#endif
        
        // select quadrature rule
        const int pOrder_s = lfsu_s.finiteElement().localBasis().order();
        const int pOrder_n = lfsu_n.finiteElement().localBasis().order();
        const int intorder = 2 * std::max( pOrder_s, pOrder_n ); // 4;

        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        // loop over quadrature points and integrate normal flux

        RF flux_jump_L2normSquare(0.0);
        RF uh_jump_L2normSquare(0.0);

        // loop over quadrature points and integrate normal flux

        for( typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin()
               ; it!=rule.end(); ++it ){

          // position of quadrature point in local coordinates of elements 
          Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(it->position());
          Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(it->position());


          // Diffusion tensor at quadrature point
          typename TP::Traits::DispersionTensorType A_s = param.DispersionTensor( *(ig.inside()), iplocal_s, equationMode );
          typename TP::Traits::DispersionTensorType A_n = param.DispersionTensor( *(ig.outside()), iplocal_n, equationMode );

          Dune::FieldVector<RF,dim> An_F_s;
          A_s.mv(n_F,An_F_s);
          Dune::FieldVector<RF,dim> An_F_n;
          A_n.mv(n_F,An_F_n);
          
          /**********************/
          /* Evaluate Functions */
          /**********************/

          // evaluate uDG_s, uDG_n
          RF uDG_s=0.0;
          evalFunction( iplocal_s, lfsu_s, x_s, uDG_s );
          
          RF uDG_n=0.0;
          evalFunction( iplocal_n, lfsu_n, x_n, uDG_n );            
          

          /**********************/
          /* Evaluate Gradients */
          /**********************/

          Dune::FieldVector<RF,dim> gradu_s(0.0);
          evalGradient( iplocal_s, *(ig.inside()), lfsu_s, x_s, gradu_s );

          Dune::FieldVector<RF,dim> gradu_n(0.0);
          evalGradient( iplocal_n, *(ig.outside()), lfsu_n, x_n, gradu_n );

          // integrate
          RF factor = it->weight() * ig.geometry().integrationElement(it->position());

          // evaluate flux jump term
          RF flux_jump = (An_F_s*gradu_s)-(An_F_n*gradu_n);
          flux_jump_L2normSquare += flux_jump * flux_jump * factor;

          // evaluate jump term
          RF jump_uDG = (uDG_s - uDG_n);
          uh_jump_L2normSquare += jump_uDG * jump_uDG * factor;

        } // loop over quadrature points on the face

        // accumulate indicator
        DF h_face = diameter(ig.geometry());

        // Liang Zhu: second term, edge residual
        RF eta_Ek_s = 0.5 * h_face * flux_jump_L2normSquare;
        RF eta_Ek_n = eta_Ek_s; // equal on both sides of the intersection!

        // Liang Zhu: third term, edge jumps
        RF eta_Jk_s = 0.5 * ( gamma / h_face + h_face / epsilon_s) * uh_jump_L2normSquare;
        RF eta_Jk_n = 0.5 * ( gamma / h_face + h_face / epsilon_n) * uh_jump_L2normSquare;
        
        // add contributions from both sides of the intersection
        r_s.accumulate( lfsv_s, 0, eta_Ek_s + eta_Jk_s );
        r_n.accumulate( lfsv_n, 0, eta_Ek_n + eta_Jk_n );
        
      } // void alpha_skeleton





      // boundary integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG,
               typename LFSU,
               typename X,
               typename LFSV,
               typename R>
      void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        
        // dimensions
        const int dim = IG::dimension;

        // evaluate permeability tensors
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);

        const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

        typename TP::Traits::RangeType beta 
          = param.velocity(*(ig.inside()),inside_local,equationMode);

        RF normalflux = beta * n_F;

        if( normalflux > 1e-10 )
          return; // No need to compare against gDirichlet in this case.

        typename TP::Traits::DispersionTensorType A_s;
        A_s = param.DispersionTensor(*(ig.inside()),inside_local,equationMode);
        
#ifdef TEST_MAX
        RF epsilon_s = std::max( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::max( A_s[2][2], epsilon_s );
#else
        RF epsilon_s = std::min( A_s[0][0], A_s[1][1]);
        if( dim>2 ) epsilon_s = std::min( A_s[2][2], epsilon_s );
#endif        

        // select quadrature rule
        const int pOrder_s = lfsu_s.finiteElement().localBasis().order();
        const int intorder = 2 * pOrder_s;
        
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule 
          = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

        // evaluate boundary condition
        const Dune::FieldVector<DF,dim-1> 
          face_local = Dune::ReferenceElements<DF,dim-1>::general(gtface).position(0,0);
        BCType bctype = param.bctype(ig.intersection(),face_local);
        if (bctype != Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet)
          return;

        // loop over quadrature points and integrate normal flux
        RF uh_jump_L2normSquare(0.0);

        for( typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin()
               ; it!=rule.end(); ++it ){

          // position of quadrature point in local coordinates of elements 
          Dune::FieldVector<DF,dim> iplocal_s 
            = ig.geometryInInside().global(it->position());
            
          // evaluate Dirichlet boundary condition
          RF gDirichlet = param.g( *(ig.inside()), iplocal_s );

            
          /**********************/
          /* Evaluate Functions */
          /**********************/

          // evaluate uDG_s, uDG_n
          RF uDG_s=0.0;
          evalFunction( iplocal_s, lfsu_s, x_s, uDG_s );
          RF jump_uDG = uDG_s - gDirichlet;
                
          // integrate
          RF factor = it->weight() * ig.geometry().integrationElement(it->position());

          uh_jump_L2normSquare += jump_uDG * jump_uDG * factor;

        }

        // accumulate indicator
        DF h_face = diameter(ig.geometry());

        // Liang Zhu: third term, edge jumps on the Dirichlet boundary
        RF eta_Jk_s = (gamma / h_face + h_face / epsilon_s) * uh_jump_L2normSquare;
        // TODO: This term is too big! Why??

        r_s.accumulate( lfsv_s, 0, eta_Jk_s );  // boundary edge
        
      }

    private:
      const IDT& inputdata;
      const SDT& setupdata;
      const TP& param;  // two phase parameter class
      const Passenger::Type passenger;
      const EQ::Mode equationMode;
      ConvectionDiffusionDGMethod::Type method;
      ConvectionDiffusionDGWeights::Type weights;
      Real gamma; // interior penalty parameter, same as alpha in class TransportOperatorDG

      template<class GEO>
      typename GEO::ctype diameter (const GEO& geo) const
      {
        typedef typename GEO::ctype DF;
        DF hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        for (int i=0; i<geo.corners(); i++)
          {
            Dune::FieldVector<DF,dim> xi = geo.corner(i);
            for (int j=i+1; j<geo.corners(); j++)
              {
                Dune::FieldVector<DF,dim> xj = geo.corner(j);
                xj -= xi;
                hmax = std::max(hmax,xj.two_norm());
              }
          }
        return hmax;
      }

    };

  }

}


#endif // DUNE_GESIS_DG_RESIDUAL_ESTIMATOR_HH
