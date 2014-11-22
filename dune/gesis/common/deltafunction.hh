#ifndef DUNE_GESIS_DELTAFUNCTION_HH
#define DUNE_GESIS_DELTAFUNCTION_HH

// If the smearing is below this level it does not make sense to use the Gaussian.
// The quadrature point might miss the peak of the Gaussian resulting in a zero solution!
#define PEAK_THRESHOLD 0.125

namespace Dune{
  namespace Gesis{

    /***********************************************************************************
     *
     * GridFunction used to define a Smoothed Pointsource via Gaussian function
     *
     ***********************************************************************************
     */
    template<typename GV, typename RF, typename COORDINATES, typename IDT>
    class DeltaFunction
      : public Dune::PDELab::GridFunctionBase
    <
      Dune::PDELab::GridFunctionTraits< GV,RF,1,Dune::FieldVector<RF,1> >
      , DeltaFunction<GV,RF,COORDINATES,IDT> >
    {
    public:
      typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
      typedef typename GV::ctype ctype;
      enum{ dim = GV::dimension }; //Traits::GridViewType::Grid::dimension;

    private:
      const GV& gv;
      const IDT& inputdata;
      const COORDINATES& global_location;
      const RF pointValue;
      RF pointSmearing;
      const bool bfixedwidth;

      Dune::FieldMatrix<CTYPE,dim,dim> Sigma;
      RF detV;
      RF scale_factor;

    public:

      //! construct from grid view
      DeltaFunction (const GV& gv_,
                     const IDT& inputdata_,
                     const COORDINATES& global_location_,
                     const RF pointValue_=1.0,
                     const RF pointSmearing_=1.0,
                     const bool bfixedwidth_=false)
        : gv(gv_)
        , inputdata(inputdata_)
        , global_location(global_location_)
        , pointValue(pointValue_)
        , pointSmearing(pointSmearing_)
        , bfixedwidth(bfixedwidth_)
        , Sigma(0) // better initialize before the computer does it for you!
      {

        if(pointSmearing < PEAK_THRESHOLD)
          pointSmearing = PEAK_THRESHOLD;  // pointSmearing must not be zero

#ifdef USE_YASP
        Vector<CTYPE> sigma ( inputdata_.domain_data.yasp_gridsizes );
#endif
#ifdef USE_ALUGRID
        Vector<CTYPE> sigma ( inputdata_.domain_data.alugrid_gridsizes );
#endif
#ifdef USE_UG
        Vector<CTYPE> sigma ( inputdata_.domain_data.ug_gridsizes );
#endif

        if( gv.grid().maxLevel() > 0 )
          sigma /= std::pow(2.0,gv.grid().maxLevel());

        if( bfixedwidth )
          sigma = pointSmearing;
        else
          sigma *= pointSmearing;

        for( int i=0;i<dim;i++ )
          Sigma[i][i] = sigma[i]*sigma[i];
        detV = Sigma.determinant();
        scale_factor = std::pow(std::sqrt(2.0*M_PI),dim)*std::sqrt(detV);
        Sigma.invert();
      }

      // http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Definition

      //! evaluate extended function on element
      inline void evaluate( const typename Traits::ElementType& e,
                            const typename Traits::DomainType& xlocal,
                            typename Traits::RangeType& y ) const
      {
        COORDINATES distance = e.geometry().global( xlocal );
        distance -= global_location;

        Dune::FieldVector<CTYPE,dim> Vdiff (0);
        Sigma.mv(distance,Vdiff);

        y = pointValue * std::exp( -0.5 * (distance * Vdiff) ) / scale_factor;

        return;
      }

      //! get a reference to the grid view
      inline const GV& getGridView ()
      {
        return gv;
      }
    };

  } // Gesis
} // Dune

#endif // DUNE_GESIS_DELTAFUNCTION_HH
