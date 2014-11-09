// -*- tab-width: 4; indent-tabs-mode: nil -*-

#ifndef DIFFERENCEADAPTERS_HH
#define DIFFERENCEADAPTERS_HH




namespace Dune {
  namespace PDELab {

    template<typename GF>
    void egCenterMaximumOfGridFunction(const GF& gf,
                                       typename GF::Traits::RangeType& maximum,
                                       unsigned qorder = 1) {
      typedef typename GF::Traits::GridViewType GV;
      typedef typename GV::template Codim<0>::
        template Partition<Interior_Partition>::Iterator EIterator;
      typedef typename GV::template Codim<0>::Geometry Geometry;
      typedef typename GF::Traits::RangeType Range;
      typedef typename GF::Traits::DomainFieldType DF;
      static const int dimD = GF::Traits::dimDomain;
      //typedef Dune::QuadratureRule<DF,dimD> QR;
      //typedef Dune::QuadratureRules<DF,dimD> QRs;
      //typedef typename QR::const_iterator QIterator;

      maximum = 0;
      Range val;
      const EIterator eend = gf.getGridView().template end<0,
        Interior_Partition>();
      for(EIterator eit = gf.getGridView().template begin<0,
            Interior_Partition>(); eit != eend; ++eit) {
        const Geometry& geo = eit->geometry();
        Dune::GeometryType gt = geo.type();
        //const QR& rule = QRs::rule(gt,qorder);
        //const QIterator qend = rule.end();

        Dune::FieldVector<DF,dimD> localcenter = Dune::ReferenceElements<DF,dimD>::general(gt).position(0,0);

        gf.evaluate(*eit,localcenter,val);

        maximum = std::max( val.infinity_norm(), maximum.infinity_norm() );

      }
    }

  }
}

/*! \brief Adapter returning f1(x)-f2(x) for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceAdapter
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                   typename T1::Traits::RangeFieldType,
                                   1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
  ,DifferenceAdapter<T1,T2> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor
  DifferenceAdapter (const T1& t1_, 
                     const T2& t2_,
                     const double neglectRadius_=1E-12
                     ) : t1(t1_), t2(t2_), neglectRadius(neglectRadius_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    const typename Traits::DomainType xglobal = e.geometry().global(x);
    if( xglobal.two_norm() < neglectRadius ) {
      y = 0;
      return;
    }

    typename Traits::RangeType y1;
    t1.evaluate(e,x,y1);
    typename Traits::RangeType y2;
    t2.evaluate(e,x,y2);
    y1 -= y2;
    y = -y1;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }

private:
  const T1& t1;
  const T2& t2;
  const double neglectRadius;
};









/*! \brief Adapter returning ||f1(x)-f2(x)||^2 for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceSquaredAdapter
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                   typename T1::Traits::RangeFieldType,
                                   1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
  ,DifferenceSquaredAdapter<T1,T2> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor
  DifferenceSquaredAdapter( const T1& t1_, 
                            const T2& t2_,
                            const double neglectRadius_=1E-12
                            ) : t1(t1_),
                                t2(t2_),
                                neglectRadius(neglectRadius_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    const typename Traits::DomainType xglobal = e.geometry().global(x);
    if( xglobal.two_norm() < neglectRadius ) {
      y = 0;
      return;
    }
    
    typename T1::Traits::RangeType y1;
    t1.evaluate(e,x,y1);
    typename T2::Traits::RangeType y2;
    t2.evaluate(e,x,y2);
    y1 -= y2;
    y = y1.two_norm2();
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }

private:
  const T1& t1;
  const T2& t2;
  const double neglectRadius;

};



#endif // DIFFERENCEADAPTERS_HH

