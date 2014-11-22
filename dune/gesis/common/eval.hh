#ifndef DUNE_GESIS_EVAL_HH
#define DUNE_GESIS_EVAL_HH


namespace Dune {
  namespace Gesis {


    //
    // Be careful: This makes sense only for the Lagrangian basis!
    //
    template<class UType>
    class EvalU{

      typedef typename UType::ElementType RF;
      typedef typename UType::const_iterator Iterator;

      const UType& u_h;
    public:
      EvalU( const UType& u_h_ ) : u_h(u_h_){};
      void extremeValues(RF& maximum, RF& minimum){
        maximum = RF(0.0);
        minimum = RF(1e99);

        std::vector<RF> x;
        u_h.std_copy_to(x);

        std::cout << "u_h.flatsize() = " << u_h.flatsize() << std::endl;
        std::cout << "x.size() = " << x.size() << std::endl;

        for( std::size_t i=0; i<x.size(); i++ ){
          maximum = std::max( maximum, x[i] );
          minimum = std::min( minimum, x[i] );
        }
      }
    };




    //#ifndef USE_CACHE

    template< class DomainType,
              class LFS,
              class BasisType>
    void evalFunctionBasis( const DomainType& location,  // expects element local coordinates!
                            const LFS& lfs,
                            BasisType& phi
                            ){
      lfs.finiteElement().localBasis().evaluateFunction(location,phi);
    }


    template< class DomainType,
              class LFS,
              class LV,
              class REAL>
    void evalFunction( const DomainType& location,  // expects element local coordinates!
                       const LFS& lfs,
                       const LV& xlocal,
                       REAL &valU
                       ){

      typedef typename LFS::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RangeType;  // normally = Dune::FieldVector<double,1>

      std::vector<RangeType> phi(lfs.size());
      evalFunctionBasis(location,lfs,phi);
      // evaluate u
      valU = REAL(0.0);
      for( std::size_t i=0; i<lfs.size(); i++ )
        valU += xlocal( lfs, i ) * phi[i];
    }





    template< class DomainType,
              class EG,
              class LFS,
              class RangeTypeVector
              >
    void evalGradientBasis( const EG& eg,
                            const DomainType& location, // expects element local coordinates!
                            const LFS& lfs,
                            RangeTypeVector &gradphi
                            ){

      typedef typename LFS::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DF;

      typedef typename LFS::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;

      typedef typename LFS::Traits::SizeType size_type;
      //const int order = lfs.finiteElement().localBasis().order();

      std::vector<JacobianType> js(lfs.size());
      lfs.finiteElement().localBasis().evaluateJacobian( location, js );

      enum { dim = DomainType::dimension };
      // transformation
      Dune::FieldMatrix<DF,dim,dim> jac;

      // transform gradients of shape functions to real element
      jac = eg.geometry().jacobianInverseTransposed( location );
      for (size_type i=0; i<lfs.size(); i++)
        jac.mv( js[i][0], gradphi[i] );

    }




    template< class DomainType,
              class EG,
              class LFS,
              class LV,
              class RType
              >
    void evalGradient( const DomainType& location, // expects element local coordinates!
                       const EG& eg,
                       const LFS& lfs,
                       const LV& xlocal,
                       RType &gradu
                       ){

      typedef typename LFS::Traits::SizeType size_type;
      enum { dim = DomainType::dimension };
      std::vector< Dune::FieldVector<REAL,dim> > gradphi(lfs.size());
      evalGradientBasis( eg,
                         location,
                         lfs,
                         gradphi
                         );

      // compute gradient of u
      gradu = RType(0);
      for (size_type i=0; i<lfs.size(); i++)
        gradu.axpy( xlocal(lfs,i), gradphi[i] );

    }



    /*
      Numerical computation of laplacian of u by central differences:
    */

    template< class DomainType,
              class EG,
              class LFS,
              class LV,
              class REAL>
    void evalNumericalLaplacian( const DomainType& location,
                                 const EG& eg,
                                 const LFS& lfs,
                                 const LV& xlocal,
                                 REAL &laplacian
                                 ){

      typedef typename LFS::Traits::SizeType size_type;

      REAL delta=1e-7; // sqrt(eps of double)

      // For DG, the central differences can be computed numerically
      // only outside a certain distance from the element border!
      enum { dim = DomainType::dimension };
      REAL minimum = 100.0;
      REAL maximum = 0.0;
      for (size_type i=0; i<dim ; i++){
        minimum = std::min( location[i], minimum );
        maximum = std::min( location[i], maximum );
      }

      if( (minimum > 1e-6) && ( maximum < 1.0-1e-6 ) ) {

        std::vector<DomainType> location_plus_delta(dim, location);
        std::vector<DomainType> location_minus_delta(dim, location);

        std::vector<REAL> u_plus(dim, 0.0);
        std::vector<REAL> u_minus(dim, 0.0);

        for (size_type j=0; j<dim ; j++){
          location_plus_delta[j][j] += delta;
          location_minus_delta[j][j] -= delta;
        }


        // evaluate u

        std::cout << "location = " << location << std::endl;
        REAL u=0.0;
        evalFunction( location,
                      lfs,
                      xlocal,
                      u
                      );

        REAL sum(0.0);

        for(size_type j=0; j<dim ; j++){
          /*
            std::cout << "location_plus_delta = "
            << std::setprecision(9)
            << location_plus_delta[j]
            << std::endl;
          */
          evalFunction( location_plus_delta[j],
                        lfs,
                        xlocal,
                        u_plus[j]
                        );
          /*
            std::cout << "location_minus_delta = "
            << std::setprecision(9)
            << location_minus_delta[j]
            << std::endl;
          */
          evalFunction( location_minus_delta[j],
                        lfs,
                        xlocal,
                        u_minus[j]
                        );

          std::cout << j << " : " << u_plus[j] << " ... " << u_minus[j]
                    << std::endl;

          std::cout << j << " => " << u_plus[j] - 2.0 * u + u_minus[j]
                    << std::endl;

          sum += u_plus[j] - 2.0 * u + u_minus[j];

        } // end of for j

        laplacian = sum/delta/delta; // central difference approximation

      } // end if within element?

    }



    //#else  // if USE_CACHE:




    template< class DomainType,
              class LFS,
              class CACHE,
              class BasisType>
    void evalFunctionBasis( const DomainType& location,  // expects element local coordinates!
                            const LFS& lfs,
                            const CACHE& cache,
                            BasisType& phi
                            ){
      const int order = lfs.finiteElement().localBasis().order();
      phi = cache[order].evaluateFunction(location,lfs.finiteElement().localBasis());
    }





    template< class DomainType,
              class LFS,
              class CACHE,
              class LV,
              class REAL>
    void evalFunction( const DomainType& location,  // expects element local coordinates!
                       const LFS& lfs,
                       const CACHE& cache,
                       const LV& xlocal,
                       REAL &valU
                       ){

      typedef typename LFS::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeType RangeType;  // normally = Dune::FieldVector<double,1>

      std::vector<RangeType> phi(lfs.size());
      evalFunctionBasis(location,lfs,cache,phi);
      // evaluate u
      valU = REAL(0.0);
      for( std::size_t i=0; i<lfs.size(); i++ )
        valU += xlocal( lfs, i ) * phi[i];
    }


    template< class DomainType,
              class LFS,
              class CACHE,
              class LV,
              class BasisType,
              class REAL
              >
    void evalFunctionAndBasis( const DomainType& location,  // expects element local coordinates!
                               const LFS& lfs,
                               const CACHE& cache,
                               const LV& xlocal,
                               BasisType& phi,
                               REAL &u
                               ){
      evalFunctionBasis(location,lfs,cache,phi);
      // evaluate u
      u=0.0;
      for( std::size_t i=0; i<lfs.size(); i++ )
        u += xlocal( lfs, i ) * phi[i];
    }




    template< class DomainType,
              class EG,
              class LFS,
              class Cache,
              class RangeTypeVector
              >
    void evalGradientBasis( const EG& eg,
                            const DomainType& location, // expects element local coordinates!
                            const LFS& lfs,
                            const Cache& cache,
                            RangeTypeVector &gradphi
                            ){

      typedef typename LFS::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::DomainFieldType DF;

      typedef typename LFS::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;

      typedef typename LFS::Traits::SizeType size_type;
      const int order = lfs.finiteElement().localBasis().order();

      const std::vector<JacobianType>& js = cache[order].evaluateJacobian( location, lfs.finiteElement().localBasis() );

      enum { dim = DomainType::dimension };
      // transformation
      Dune::FieldMatrix<DF,dim,dim> jac;

      // transform gradients of shape functions to real element
      jac = eg.geometry().jacobianInverseTransposed( location );
      for (size_type i=0; i<lfs.size(); i++)
        jac.mv( js[i][0], gradphi[i] );

    }



    template< class DomainType,
              class EG,
              class LFS,
              class CACHE,
              class LV,
              class RType
              >
    void evalGradient( const DomainType& location, // expects element local coordinates!
                       const EG& eg,
                       const LFS& lfs,
                       const CACHE& cache,
                       const LV& xlocal,
                       RType &gradu
                       ){

      typedef typename LFS::Traits::SizeType size_type;
      enum { dim = DomainType::dimension };
      std::vector< Dune::FieldVector<REAL,dim> > gradphi(lfs.size());
      evalGradientBasis( eg,
                         location, // expects element local coordinates!
                         lfs,
                         cache,
                         gradphi );

      // compute gradient of u
      gradu = RType(0);
      for (size_type i=0; i<lfs.size(); i++)
        gradu.axpy( xlocal(lfs,i), gradphi[i] );

    }

    //#endif // USE_CACHE
  } // Gesis
} // Dune

#endif // DUNE_GESIS_EVAL_HH
