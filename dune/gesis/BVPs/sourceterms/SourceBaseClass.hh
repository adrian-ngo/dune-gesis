#ifndef DUNE_GESIS_SOURCE_BASE_CLASS_HH
#define DUNE_GESIS_SOURCE_BASE_CLASS_HH

namespace Dune {
  namespace Gesis {

    /**********************************************************
     *
     *  SourceNature is an enum-type specifying what kind
     *  of source-terms are possible
     *
     *
     **********************************************************
     */
    enum SourceNature {
      NULL_SOURCE,
      POINT_SOURCE,
      PUMPING_SOURCES,  // several point sources with values
      SOLUTE_PUMPING_SOURCES,
      HEAT_PUMPING_SOURCES,
      GEOELECTRICAL_POTENTIAL_SOURCE,
      GEOELECTRICAL_POINT_SOURCE,
      FUNCTIONAL_SOURCE,
      GEP_FUNCTIONAL_SOURCE,
      TP_FUNCTIONAL_SOURCE
    };


    //============================================================
    //! Define a traits class for both flow and transport problems
    //============================================================
    template<typename GV, typename RF>
    struct SourceTermTraits
    {
      // brief Export type for range field
      typedef RF RangeFieldType;

      //! \brief the grid view
      typedef GV GridViewType;

      //! \brief Enum for domain dimension

      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      typedef typename GV::Grid::ctype DomainFieldType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;


      //! \brief range type
      typedef Dune::FieldVector<RF, GV::dimensionworld> RangeType;
      //! grid types
      typedef typename GV::Traits::template Codim< 0 >::Entity ElementType;

    };

    //=================================
    //! base class for the source classes
    //=================================
    template<class T, class Imp>
    class SourceTermInterface {
    public:
      typedef T Traits;
      typedef typename T::RangeFieldType RF;

      enum{dim=Traits::dimDomain};

      template<
        typename EG
        , typename LFSV
        , typename SF
        , typename R
        >
      bool evaluate_residual_on_element( const EG& eg
                                         , const LFSV& lfsv
                                         , SF& shapefunctions
                                         , R& residual
                                         , const bool bAdjoint=false
                                         ) const {
        return true;
        //return asImp().evaluate( eg, lfsv, shapefunctions, residual );
      }


      template<typename IG, typename LFSV, typename R>
      void evaluate_residual_on_skeleton (const IG& ig,
                                          const LFSV& lfsv_s, const LFSV& lfsv_n,
                                          R& r_s, R& r_n
                                          ) const {
        return;
      }


      template<typename IG, typename LFSV, typename R>
      void evaluate_residual_on_boundary (const IG& ig,
                                          const LFSV& lfsv_s,
                                          R& r_s
                                          ) const {
        return;
      }



      template<typename EG, typename COORDINATES>
      bool evaluate_function( const EG& e
                              , const COORDINATES& xlocal
                              , REAL& fval
                              ) const
      {
        return true;
        //return asImp().evaluate( e, xlocal, fval);
      }

      // base-class
      template<typename EG, typename COORDINATES>
      bool evaluate_vector( const EG& e
                            , const COORDINATES& xlocal
                            , Dune::FieldVector<RF,dim>& fvec
                            ) const
      {
        return true;
      }


      template<typename COORDINATES>
      void set_PointSourceLocation( const COORDINATES& pointsource_global )
      {
        return;
      }

      template<typename COORDINATES>
      void set_PointSourceLocations( const COORDINATES& pointsource1_global,
                                     const COORDINATES& pointsource2_global )
      {
        return;
      }


      void reset_rhs() {
        return;
      }

      template<typename VCType>
      void set_rhs( const VCType& vc_solution )
      {
        return;
      }

      template<typename VCType>
      void set_rhs( const VCType& vc_1, const VCType& vc_2 ) {
        return;
      }



    private:

      Imp& asImp() {
        return static_cast<Imp &> (*this);
      }

      const Imp& asImp() const {
        return static_cast<const Imp &> (*this);
      }

    };

  }
}
#endif // DUNE_GESIS_SOURCE_BASE_CLASS_HH
