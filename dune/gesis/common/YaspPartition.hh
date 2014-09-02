#ifndef DUNE_GEO_INVERSION_YASP_PARTITION_HH
#define DUNE_GEO_INVERSION_YASP_PARTITION_HH

namespace Dune {

  namespace GeoInversion {

    /*
      With this class you can specify how to distribute the total number of
      processes to the YASP grid by passing an integer vector
      to the constructor.
    */
    template<typename IntegerVector,int dim>
    class YaspPartition : public Dune::YLoadBalance<dim>
    {
    private:
      typedef Dune::FieldVector<int,dim> iTupel;
      const IntegerVector& yasppartitions;
      
    public:
      //constructor:
      YaspPartition( const IntegerVector& yasppartitions_ )
        : yasppartitions( yasppartitions_ )
      {
      }
      
      void loadbalance (const iTupel& size, int P, iTupel& dims) const
      {
        for(UINT i=0;i<dim;i++)
          dims[i] = yasppartitions[i];
      }

    };

  }

}


#endif // DUNE_GEO_INVERSION_YASP_PARTITION_HH
