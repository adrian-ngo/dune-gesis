// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GESIS_GRIDFUNCTIONSPACEUTILITIES_HH
#define DUNE_GESIS_GRIDFUNCTIONSPACEUTILITIES_HH

#include <math.h>
#include <cstdlib>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

//#include <dune/pdelab/common/jacobiantocurl.hh>
#include <dune/pdelab/constraints/common/constraintstransformation.hh> // backward compatibility

//#include"dune/pdelab/common/countingptr.hh"
//#include"dune/pdelab/common/multitypetree.hh"
//#include"dune/pdelab/common/cpstoragepolicy.hh"

#include"dune/pdelab/common/function.hh"
#include"dune/pdelab/gridfunctionspace/gridfunctionspace.hh"
#include"dune/pdelab/gridfunctionspace/localfunctionspace.hh"

#include<dune/localfunctions/raviartthomas/raviartthomascube.hh>

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

namespace Dune {
  namespace Gesis {

    //! \brief Convert the pressure h (scalar grid function) into a 
    //! vector-valued grid function (with dim components) 
    //! representing the Darcy flux q = -K grad h
    /**
     * The function values should be single-component vectors.  
     * The gradient will be a dim-component function.
     *
     * \tparam GWP Problem parameters: Diffusion-tensor and Dirichlet BC type and values
     * \tparam GFS Grid function space
     * \tparam IDT Input data type: location of wells
     */
    template<typename GWP,
             typename GFS>
    class DiscreteGridFunctionDarcy
      : public TypeTree::LeafNode,
        public Dune::PDELab::GridFunctionInterface<
      Dune::PDELab::GridFunctionTraits<
        typename GFS::Traits::GridViewType,
        typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
        GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain,
        FieldVector<
          typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
          GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::dimDomain> >,
      DiscreteGridFunctionDarcy<GWP,GFS> >
    {

    public:
      typedef typename GFS::Traits::GridViewType GV;
      typedef typename GV::ctype DF;      
      typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;
      typedef typename LBTraits::RangeFieldType RF;
      enum{ dim = GV::dimension };
      typedef Dune::FieldVector<DF,dim> Domain;
      typedef Dune::FieldVector<RF,dim> VectorRange;

      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type VCType;

    public:
      //typedef typename LBTraits::dimDomain dim;
      //typedef GridFunctionTraits< GV, RF, LBTraits::dimDomain, Dune::FieldVector<RF,LBTraits::dimDomain> > Traits;
      typedef Dune::PDELab::GridFunctionTraits< GV, RF, dim, Dune::FieldVector<RF,dim> > Traits;

    private:

      typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
      typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
      typedef typename VCType::template ConstLocalView<LFSCache> XView;

      const GWP& gwp;
	  shared_ptr<GFS const> pgfs;

      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<typename Traits::RangeFieldType> xl;

	  // const VCType& xg;

      const int baselevel;
      const bool withoutK;  // "false" evaluates darcy flux, "true" evaluates gradient
      bool well_included;

    public:
      /** \brief Construct a DiscreteGridFunctionDarcy
       *
       * \param gfs The GridFunctionsSpace
       * \param x_  The coefficients vector
       * \param k_  The conductivity field
       */
      DiscreteGridFunctionDarcy( const GWP& gwp_
                                 , const GFS& gfs
                                 , const VCType& x_
                                 , const int baselevel_=0
                                 , const bool withoutK_=false
                                 , bool well_included_=false // TODO: Ask Ronnie
                                 )
        : gwp( gwp_ )
        , pgfs( stackobject_to_shared_ptr( gfs ) )
          // xg( x_ ),
        , lfs(gfs)
        , lfs_cache(lfs)
        , x_view(x_)
        , xl(lfs.size())
        //
        , baselevel(baselevel_)
        , withoutK (withoutK_)
        , well_included( well_included_ )
      { }
      


      inline int getLocalBasisOrder() const {

        Dune::PDELab::LocalFunctionSpace<GFS> lfs(*pgfs);
        return lfs.finiteElement().localBasis().order();
        
      }


      inline void evaluate_on_root( const typename Traits::ElementType& e, 
                                    const typename Traits::DomainType& x, 
                                    typename Traits::RangeType& flux ) const {

        // convert to global coordinate wrt to element e
        typename Traits::DomainType global = e.geometry().global(x);
        if(e.level()>baselevel){
          // get father element
          typedef typename GFS::Traits::GridViewType::Grid::template Codim<0>::EntityPointer ElementPointer;
          ElementPointer pAncestor = e.father();
          while( pAncestor->level() > baselevel )
            pAncestor = (*pAncestor).father();
          // convert to local coordinate wrt to element *pAncestor
          typename Traits::DomainType xx = (*pAncestor).geometry().local( global );
          this->evaluate( *pAncestor, xx, flux );
        }
        else{
          this->evaluate(e, x, flux);
        }

        return;

      }
      
      
      
      // Evaluate
      inline void evaluate_on_leaf( const int maxGridLevel,
                                    const int baselevel,
                                    const typename Traits::ElementType& e,
                                    const typename Traits::DomainType& x,
                                    typename Traits::RangeType& y
                                   ) const
      {
          
        if( maxGridLevel == baselevel ){
          evaluate( e, x, y );
          //std::cout << "DEBUG: DiscreteGridFunctionDarcy::evaluate_on_leaf() with baselevel = " 
          //<< baselevel << std::endl;
          return;
        }

        const typename Traits::DomainType xglobal = e.geometry().global(x);

        //std::cout << "x = " << x << std::endl;
        //std::cout << "xglobal = " << xglobal << std::endl;

        const typename Traits::ElementType::HierarchicIterator& hbegin 
          = e.hbegin(maxGridLevel);
        const typename Traits::ElementType::HierarchicIterator& hend 
          = e.hend(maxGridLevel);
        int counter=0;
        for (typename Traits::ElementType::HierarchicIterator hit = hbegin;
             hit != hend; ++hit) {
          // only evaluate on entities with data
          if ((*hit).isLeaf()) {
            //std::cout << "a leaf... volume = " << (*hit).geometry().volume() << std::endl;
            const typename Traits::DomainType hit_local =
              (*hit).geometry().local( xglobal );
            //std::cout << "hit_local = " << hit_local << std::endl;
            bool bInside = true;
            for( int i=0; i<dim; i++ ){
              if( hit_local[i]<0.0 || hit_local[i] > 1.0 ){
                bInside=false;
                break;
              }
            }
            if(bInside){
              counter++;

              //std::cout << "point is inside cell with center " 
              //          << (*hit).geometry().center()
              //          << std::endl;

              typename Traits::RangeType y_contribution;
              evaluate( (*hit), hit_local, y_contribution );

              //std::cout << "y_contribution = " << y_contribution << std::endl;

              y += y_contribution;
            }
          }
        }
        if( counter>1 ){
          std::cout << "Taking mean for y: counter = " 
                    << counter 
                    << " y = " << y 
                    << std::endl;
          y /= (REAL) counter;
        }
        //std::cout << "counter = " << counter << std::endl;
        //std::cout << "y = " << y << std::endl;

        return;
      }


      // Evaluate
      inline void evaluate(
                           const typename Traits::ElementType& e,
                           const typename Traits::DomainType& x,
                           typename Traits::RangeType& y
                           ) const
      {
        typedef FiniteElementInterfaceSwitch<typename Dune::PDELab::LocalFunctionSpace<GFS>::Traits::FiniteElementType> FESwitch;
        
        // get and bind local functions space
        //LocalFunctionSpace<GFS> lfs(*pgfs);
        lfs.bind( e );
        lfs_cache.update();
        x_view.bind(lfs_cache);
        // get local coefficients
        xl.resize(lfs.size());
        x_view.read( xl );
        x_view.unbind();

        // get local coefficients
        // std::vector<typename Traits::RangeFieldType> xl( lfs.size() );
        // lfs.vread( xg, xl );

        // get local Jacobians/gradients of the shape functions
        std::vector<typename LBTraits::JacobianType> js( lfs.size() );
        
        FESwitch::basis( lfs.finiteElement() ).evaluateJacobian( x, js );

        // get Jacobian of geometry
        // const typename Traits::ElementType::Geometry::Jacobian&
        const Dune::FieldMatrix<REAL,dim,dim>
          jac = e.geometry().jacobianInverseTransposed( x );
        
        //Dune::GeometryType gt = e.geometry().type();

        typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(e,x));
        //typename Traits::RangeFieldType kk = tensor_inside[0][0];

        typename Traits::RangeType gradphi;
        y = 0;
        for(unsigned i = 0; i < lfs.size(); ++i) {
          // compute global gradient of shape function i
          gradphi = 0;
          //jac.umv( js[i][0], gradphi );  // A.umv(x,y) bedeutet: y += Ax
          jac.mv( js[i][0], gradphi );  // A.mv(x,y) bedeutet: y = Ax

          if( withoutK ){
            // sum up global gradients, weighting them with the appropriate coeff
            y.axpy( xl[i], gradphi );
          }
          else{
            // sum up global gradients, weighting them with the appropriate coeff
            //typename Traits::RangeFieldType ku_i = - kk * xl[i];

            Dune::FieldVector<RF,dim> Kgradphi(0.0);
            tensor_inside.umv(gradphi,Kgradphi);
            y.axpy( -xl[i], Kgradphi ); // y += -xl[i] * Kgradphi
          }
        }

        return;

        // SDFEM seems to be very sensitive on fluctuations in beta, even at 16 digits after the comma!

        
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return pgfs->gridview();
      }

    };










    /** @brief Provide velocity/flux field
     *
     *  Uses RT0 interpolation on a cell.
     *
     *  @tparam Parameter     provides FlowEquationParameterInterface
     *  @tparam ScalarDGF P0 function for scalar
     */
    template<typename GWP,typename GFS>
    class VectorDGF
      : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename GFS::Traits::GridViewType,
        typename GWP::Traits::RangeFieldType, 
        GWP::Traits::dimDomain, // == GV::dimension
        Dune::FieldVector<typename GWP::Traits::RangeFieldType,GWP::Traits::dimDomain>
        >,
      VectorDGF<GWP,GFS> >
    {
      
      typedef typename GFS::Traits::GridViewType GV;

      enum {dim = GV::dimension};
      typedef typename GV::ctype DF;      

      typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;
      typedef typename LBTraits::RangeFieldType RF;

      typedef Dune::FieldVector<RF,dim> VectorRange;
      typedef Dune::FieldVector<DF,dim> Domain;
      typedef Dune::FieldVector<DF,dim-1> IntersectionDomain;

      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename GV::Intersection Intersection;
      
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type VCType;
      typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> ScalarDGF;
      
      typedef typename Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

	  shared_ptr<GFS const> pgfs;
      const GWP& gwp;
      ScalarDGF scalar;
      const int baselevel;
      const bool withoutK;  // "false" evaluates darcy flux, "true" evaluates gradient

#ifdef DIMENSION3
      typedef Dune::RT0Cube3DLocalFiniteElement<DF,RF> RT0CubeLocalFiniteElement;
#else
      typedef Dune::RT0Cube2DLocalFiniteElement<DF,RF> RT0CubeLocalFiniteElement;
#endif
      RT0CubeLocalFiniteElement rt0fe;
      typedef typename RT0CubeLocalFiniteElement::Traits::LocalBasisType::Traits::RangeType RT0RangeType;


      public:
      
      // No copy constructor allowed, this class contains reference members. 
      // Use pointers instead.

      // overloaded constructor, version 1 of 2 (stationary equation):
      VectorDGF ( const GWP& gwp_,
                  const GFS& gfs_,
                  const VCType& vc_,
                  const int baselevel_=0,
                  const bool withoutK_=false
                  ) :
        pgfs( stackobject_to_shared_ptr( gfs_ ) ),
        gwp(gwp_),
        scalar(gfs_,vc_), 
        baselevel(baselevel_),
        withoutK (withoutK_)
      {}


      inline int getLocalBasisOrder() const {

        Dune::PDELab::LocalFunctionSpace<GFS> lfs(*pgfs);
        return lfs.finiteElement().localBasis().order();
        
      }


      inline void evaluate_on_root( const Element& e, 
                                    const Domain& x, 
                                    //const int maxGridLevel,
                                    VectorRange& flux ) const{

        // Map each cell to unique id
        typedef typename GFS::Traits::GridViewType GV;
        const GV &gv = pgfs->gridView();
        //Dune::PDELab::MultiGeomUniqueIDMapper<GV> cell_mapper(gv);

        // convert to global coordinate wrt to element e
        Domain global = e.geometry().global(x);

          
        if(e.level()>baselevel){
          // get father element
          typedef typename GFS::Traits::GridViewType::Grid::template Codim<0>::EntityPointer ElementPointer;
          ElementPointer pAncestor = e.father();
          while( pAncestor->level() > baselevel )
            pAncestor = (*pAncestor).father();
          // convert to local coordinate wrt to element *pAncestor
          Domain xx = (*pAncestor).geometry().local( global );

          //const typename GV::IndexSet::IndexType ids = cell_mapper.map(*pAncestor);
          //std::cout << "Coarse Level a: I am cell number " << ids 
          //          << " with center " 
          //          << global
          //          << std::endl;
          this->evaluate( *pAncestor, xx, flux );
        }
        else{
          //const typename GV::IndexSet::IndexType ids = cell_mapper.map(e);
          //std::cout << "Coarse Level b: I am cell number " << ids 
          //          << " with center " 
          //          << global
          //          << std::endl;
          this->evaluate(e, x, flux);
        }

        //if( bAdjoint )
        //  flux*=-1.0;   // Reverse velocity field for adjoint case!
        return;
      }






      inline void evaluate( const Element& e, 
                            const Domain& x, 
                            VectorRange& y,
                            bool bPiola=true
                            ) const
      {

#ifdef CONSTANT_VELOCITY
        y[0] = 1.0;
        y[1] = 0.3;
        return;
#endif

        // cell geometry
        const Domain& insideCellCenterLocal  = Dune::ReferenceElements<DF,dim>::general(e.type()).position(0,0);
        Domain        insideCellCenterGlobal = e.geometry().global(insideCellCenterLocal);

        typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(e,insideCellCenterLocal));

        // absolute permeability
        RF K_inside;
        if(withoutK)
          K_inside = -1.;
        else{
          // Evaluation of K_inside is depening on the direction of the 
          // face to the neighbor.
          // Therefore it must be done in the intersection loop below.
        }

        // scalar evaluation
        Dune::FieldVector<RF,1> pInside;
        scalar.evaluate(e,insideCellCenterLocal,pInside);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        RF coeff[2*dim]; // RT0 coefficient
        Dune::FieldMatrix<DF,dim,dim>
          B = e.geometry().jacobianInverseTransposed(x); // the transformation. Assume it is linear
        RF determinant = B.determinant();

        // loop over cell neighbors
        IntersectionIterator endit = scalar.getGridView().iend(e);
        for (IntersectionIterator iit = scalar.getGridView().ibegin(e); iit!=endit; ++iit)
        {

          // set to zero for processor boundary
          vn[iit->indexInInside()] = 0.0;

          // face geometry
          const IntersectionDomain& faceLocal = Dune::ReferenceElements<DF,dim-1>::general(iit->geometry().type()).position(0,0);

          //std::cout << "DEBUG: iit->indexInInside() = " 
          //          << iit->indexInInside()
          //          << std::endl;
          
          //Domain faceGlobal = iit->geometryInInside().global(faceLocal);
          //std::cout << "DEBUG: faceGlobal = " 
          //          << faceGlobal
          //          << std::endl;


          RF element_volume_s = e.geometry().volume();
#ifdef DIMENSION3
          RF element_length_s = std::pow( element_volume_s, 1.0/3.0 );
#else
          RF element_length_s = std::sqrt( element_volume_s );
#endif

          // interior face
          if (iit->neighbor())
            {

              RF element_volume_n = iit->outside()->geometry().volume();
#ifdef DIMENSION3
              RF element_length_n = std::pow( element_volume_n, 1.0/3.0 );
#else
              RF element_length_n = std::sqrt( element_volume_n );
#endif

              const Domain& outsideCellCenterLocal = 
                Dune::ReferenceElements<DF,dim>::general(iit->outside()->type()).position(0,0);
              Domain distance_vector = 
                iit->outside()->geometry().global(outsideCellCenterLocal); 
              
              // distance of cell centers
              distance_vector -= insideCellCenterGlobal;
              RF distance = distance_vector.two_norm();
              distance_vector /= distance;

              //std::cout << "DEBUG: distance = " 
              //          << distance
              //          << std::endl;

              // absolute permeability
              RF K_outside;
              if(withoutK)
                K_outside = -1.;
              else{

                Dune::FieldVector<DF,dim> kvector_s;
                tensor_inside.mv(distance_vector,kvector_s);
                K_inside = kvector_s.infinity_norm();

                typename GWP::Traits::DiffusionTensorType tensor_outside(gwp.DiffusionTensor(*(iit->outside()),outsideCellCenterLocal));
                Dune::FieldVector<DF,dim> kvector_n;
                tensor_outside.mv(distance_vector,kvector_n);
                K_outside = kvector_n.infinity_norm();


              }
              // scalar evaluation
              Dune::FieldVector<RF,1> pOutside;
              scalar.evaluate(*(iit->outside()),outsideCellCenterLocal,pOutside);

              // liquid phase calculation
              RF w = (pOutside-pInside)/distance; // determines direction

              // set coefficient
              if( element_length_s - element_length_n<1E-12 )
                vn[iit->indexInInside()] = 
                  - Dune::Gesis::General::harmonicAverage( K_inside, K_outside ) * w;
              else
                vn[iit->indexInInside()] = 
                  - Dune::Gesis::General::harmonicAverageWeightedByDistance( K_inside, 
                                                       K_outside, 
                                                       element_length_s, 
                                                       element_length_n ) * w;
              
              //std::cout << "DEBUG: vn[iit->indexInInside()] = " 
              //          << vn[iit->indexInInside()]
              //          << std::endl;


            }

          // boundary face
          if (iit->boundary())
          {
            // distance of cell center to boundary
            Domain distance_vector = iit->geometry().global(faceLocal);
            distance_vector -= insideCellCenterGlobal;
            RF distance = distance_vector.two_norm();
            distance_vector /= distance;

            // evaluate boundary condition type
            //int bc = parameter.bc(*iit,faceLocal,time);
            // liquid phase Dirichlet boundary
            //if (bc==1) 

            BCType bctype = gwp.bctype( *iit, faceLocal );
            if( !withoutK && 
                bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet ) {
              
              // evaluate Dirichlet boundary condition
              typename GWP::Traits::RangeFieldType g = gwp.g( *(iit->inside()), iit->geometry().global(faceLocal) );

              Dune::FieldVector<DF,dim> kvector_s;
              tensor_inside.mv(distance_vector,kvector_s);
              K_inside = kvector_s.infinity_norm();

              RF w = (g-pInside)/distance;
              vn[iit->indexInInside()] = - K_inside * w;
            }
            
            else {
              // Neumann:
              RF j = 0.0; // parameter.j(*iit,faceLocal,time);
              vn[iit->indexInInside()] = j;
            }

          }

          // compute coefficient
          Domain vstar=iit->unitOuterNormal(faceLocal); // normal on transformed element
          vstar *= vn[iit->indexInInside()];

          Domain vstarhat(0);
          B.umtv(vstar,vstarhat); // Piola backward transformation

          vstarhat *= determinant;

          VectorRange normalhat(0); // normal on reference element
          if (iit->indexInInside()%2==0)
            normalhat[iit->indexInInside()/2] = -1.0;
          else
            normalhat[iit->indexInInside()/2] =  1.0;

          coeff[iit->indexInInside()] = vstarhat*normalhat;
        }

        // compute velocity on reference element
        std::vector<RT0RangeType> rt0vectors(rt0fe.localBasis().size());
        rt0fe.localBasis().evaluateFunction(x,rt0vectors);
        VectorRange yhat(0);

        //Domain xglobal = e.geometry().global(x);
        //std::cout << "DEBUG:  xglobal = " << xglobal
        //          << "      x = " << x
        //          << std::endl;

        for (unsigned int i=0; i<rt0fe.localBasis().size(); i++){
          
          //std::cout << "DEBUG:  i = " << i
          //          << "  coeff[i] = " << coeff[i] 
          //          << "  rt0vectors[i] = " << rt0vectors[i]
          //          << std::endl;

          yhat.axpy(coeff[i],rt0vectors[i]);
        }

        if( bPiola ){
          // apply Piola transformation
          B.invert();
          y = 0;
          B.umtv(yhat,y);
          y /= determinant;
        }
        else{
          y = yhat;
        }
      }


      //inline const typename Traits::GridView& getGridView ()
      //{
      //  return scalar.getGridView();
      //}

      private:

      template<typename T>
        T aavg (T a, T b) const
        {
          return 0.5*(a+b);
        }

      template<typename T>
      T havg (T a, T b) const
      {
        T eps = 1E-20;
        //return 2.0/(1.0/(a+eps) + 1.0/(b+eps));
        return T(2.0)*a*b / ( a + b + eps );
      }
    }; // class VectorDGF 



    // 
    // This class is inherited from DiscreteGridFunction.
    // Extra feature:
    // New method evaluate_on_leaf( maxGridLevel,...)
    // used to evaluate solution of TPE for the source term of the 
    // combi-adjoint of the GWE. Background: The solution of the TPE lives on a 
    // refined grid whereas the solution of the GWE lives on the coarsest grid level 0.
    // 
	template<typename GFS, typename VCType>
	class DiscreteRefinedGridFunction: public Dune::PDELab::DiscreteGridFunction<GFS,VCType>{

    private:
      const int baselevel;

    public:
      typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> base_type;
	  typedef typename base_type::Traits Traits;
      enum {dim = GFS::Traits::GridViewType::dimension};

      DiscreteRefinedGridFunction( const GFS& gfs_, const VCType& x_, const int baselevel_=0 )
        : base_type( gfs_, x_ ), baselevel(baselevel_) {}


      inline void evaluate_on_root( const typename Traits::ElementType& e,
                                    const typename Traits::DomainType& x,
                                    //const int maxGridLevel,
                                    typename Traits::RangeType& y) const{

        // Map each cell to unique id
        typedef typename GFS::Traits::GridViewType GV;
        const GV &gv = this->getGridView();
        //Dune::PDELab::MultiGeomUniqueIDMapper<GV> cell_mapper(gv);

        // convert to global coordinate wrt to element e
        typename Traits::DomainType global = e.geometry().global(x);
          
        if(e.level()>baselevel){
          // get father element
          typedef typename GFS::Traits::GridViewType::Grid::template Codim<0>::EntityPointer ElementPointer;
          ElementPointer pAncestor = e.father();
          while( pAncestor->level() > baselevel )
            pAncestor = (*pAncestor).father();
          // convert to local coordinate wrt to element *pAncestor
          typename Traits::DomainType xx = (*pAncestor).geometry().local( global );
          this->evaluate( *pAncestor, xx, y );
        }
        else{
          this->evaluate(e, x, y);
        }

        return;

      }


	  inline void evaluate_on_leaf ( const int maxGridLevel,
                                     const int baselevel,
                                     const typename Traits::ElementType& e,
                                     const typename Traits::DomainType& x,
                                     typename Traits::RangeType& y) const
      {
        y = 0;
        
        if( maxGridLevel == baselevel ){
          this->evaluate( e, x, y );
          //std::cout << "DEBUG: DiscreteRefinedGridFunction::evaluate_on_leaf() with baselevel = " << baselevel << std::endl;
          return;
        }

        // The following is actually not needed if we are using
        // the DG solution projected onto the CCFV space!

        const typename Traits::DomainType xglobal = e.geometry().global(x);

        //std::cout << "DEBUG: xglobal = " << xglobal << std::endl;
        //std::cout << "DEBUG: xlocal = " << x << std::endl;

        const typename Traits::ElementType::HierarchicIterator& hbegin 
          = e.hbegin(maxGridLevel);
        const typename Traits::ElementType::HierarchicIterator& hend 
          = e.hend(maxGridLevel);
        int counter=0;
        for (typename Traits::ElementType::HierarchicIterator hit = hbegin;
             hit != hend; ++hit) {
          // only evaluate on entities with data
          if ((*hit).isLeaf()) {
            
            //std::cout << "DEBUG: leaf level = " 
            //          << (*hit).level() << std::endl;

            const typename Traits::DomainType hit_local =
              (*hit).geometry().local( xglobal );
            
            //std::cout << "DEBUG: hit_local = " << hit_local << std::endl;

            bool bInside = true;
            for( int i=0; i<dim; i++ ){
              if( hit_local[i]<0.0 || hit_local[i] > 1.0 ){
                bInside=false;
                break;
              }
            }
            if(bInside){
              counter++;

              //std::cout << "DEBUG: point is inside sub-cell with center " 
              //          << (*hit).geometry().center()
              //          << std::endl;

              typename Traits::RangeType y_contribution;
              this->evaluate( (*hit), hit_local, y_contribution );

              //std::cout << "DEBUG: y_contribution = " << y_contribution << std::endl;

              y += y_contribution;
            }
          }
        }
        if( counter>0 )
          y /= (typename Traits::RangeType) counter;

        //std::cout << "DEBUG: counter = " << counter << std::endl;
        //std::cout << "DEBUG: y = " << y << std::endl;

        return;
      }

	};  // class DiscreteRefinedGridFunction

  } // namespace Gesis

} // namespace Dune

#endif // DUNE_GESIS_GRIDFUNCTIONSPACEUTILITIES_HH
