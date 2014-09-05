#ifndef DUNE_GEO_INVERSION_RT0_FLUX_HH
#define DUNE_GEO_INVERSION_RT0_FLUX_HH

#include<dune/localfunctions/raviartthomas/raviartthomascube.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

namespace Dune {
  namespace GeoInversion {

    /** @brief Provide velocity/flux field
     *
     *  Uses RT0 interpolation on a cell.
     *
     *  @tparam Parameter     provides FlowEquationParameterInterface
     *  @tparam ScalarDGF P0 function for scalar
     */
    template<typename GWP,typename GFS>
    class RT0FluxDGF
      : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename GFS::Traits::GridViewType,
        typename GWP::Traits::RangeFieldType, 
        GWP::Traits::dimDomain, // == GV::dimension
        Dune::FieldVector<typename GWP::Traits::RangeFieldType,GWP::Traits::dimDomain>
        >,
          RT0FluxDGF<GWP,GFS> >
    {
    public:
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

    private:
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
      RT0FluxDGF ( const GWP& gwp_,
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

      /*
      inline int getLocalBasisOrder() const {

        LocalFunctionSpace<GFS> lfs(*pgfs);
        return lfs.finiteElement().localBasis().order();
        
      }
      */

      inline void evaluate_on_root( const Element& e, 
                                    const Domain& x, 
                                    //const int maxGridLevel,
                                    VectorRange& flux ) const{

        // Map each cell to unique id
        //typedef typename GFS::Traits::GridViewType GV;
        //const GV &gv = pgfs->gridView();
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
                            const Domain& xlocal,
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


        //std::cout << "DEBUG: insideCellCenterGlobal = " 
        //          << insideCellCenterGlobal
        //          << std::endl;


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
        //RF coeff[2*dim]; // RT0 coefficient
        Dune::FieldMatrix<DF,dim,dim>
          B = e.geometry().jacobianInverseTransposed(xlocal); // the transformation. Assume it is linear

#ifdef OLD_RT0
        RF determinant = B.determinant();
#endif

        // intersection loop:
        // loop over cell neighbors
        IntersectionIterator endit = scalar.getGridView().iend(e);
        for (IntersectionIterator iit = scalar.getGridView().ibegin(e); iit!=endit; ++iit)
          {

            // set to zero for processor boundary
            vn[iit->indexInInside()] = 0.0;

            // face geometry
            const IntersectionDomain& faceLocal = Dune::ReferenceElements<DF,dim-1>::general(iit->geometry().type()).position(0,0);

            //std::cout << "DEBUG: iit->indexInInside() = " 
            //<< iit->indexInInside()
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
                    - Dune::GeoInversion::General::harmonicAverage( K_inside, K_outside ) * w;
                else
                  vn[iit->indexInInside()] = 
                    - Dune::GeoInversion::General::harmonicAverageWeightedByDistance( K_inside, 
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

#ifdef OLD_RT0

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
#endif
          } // end of intersection loop

#ifdef OLD_RT0

        // compute velocity on reference element
        std::vector<RT0RangeType> rt0vectors(rt0fe.localBasis().size());
        rt0fe.localBasis().evaluateFunction(xlocal,rt0vectors);
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

        //std::cout << "flux = " << y << std::endl;
#else

        // Note:
        // You have to be sure that the interfaces are indexed like this:
        // In 2D:
        // 
        // y^
        //  |  +-----3-----+
        //  |  |           |
        //  |  0           1
        //  |  |           |
        //  |  +-----2-----+
        //  |
        //  +-----------------> x
        // This works on YASP.

        y[0] = -vn[0] + (vn[1]+vn[0]) * xlocal[0];
        y[1] = -vn[2] + (vn[3]+vn[2]) * xlocal[1];
        if(dim>2)
          y[dim-1] = -vn[2*dim-2] + (vn[2*dim-1]+vn[2*dim-2]) * xlocal[dim-1];
        //    y[2] = -vn[4]       + (vn[5]      +vn[4] ) * xlocal[2];
        //std::cout << "y = " << y << std::endl;
#endif
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
    }; // class RT0FluxDGF 

  } // GeoInversion

} // Dune
#endif // DUNE_GEO_INVERSION_RT0_FLUX_HH
