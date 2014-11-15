#ifndef DUNE_GESIS_ORDERING_SCHEMES_HH
#define DUNE_GESIS_ORDERING_SCHEMES_HH


namespace Dune{

  namespace Gesis{


    /** *********************************************************** *
     * Structure for comparing two head values
     * ************************************************************ 
     */
    struct PressureLikeOrdering {
      PressureLikeOrdering()
      {};
      template <typename PressureType>
      bool operator()(const PressureType& value1,
                      const PressureType& value2) {
        //
        // This is very tricky:
        // ====================
        // Sort the grid elements in the same ascending order as the
        // hydraulic head h! This is essential for writing and reading 
        // the m0 or m1 solution in parallel mode to the HDF5 file
        // using different numbers of processes for writing and retrieving
        // forward solutions m0 or m1 as is required during the inversion 
        // with more than one MPI_Pools where a repartitioning of the
        // grid is necessary.
        //
        if( value1 < value2 )
          return true;
        else
          return false;

      }

    }; // used by PositionIndexCompare


    /** *********************************************************** *
     * Structure for comparing two head values
     * To the linear solver is does not matter whether we are sorting
     * in ascending (against the flow) or descending (in flow direction)
     * order, because BiCGStab is designed to cope with 
     * non-sysmmetric matrices!
     * ************************************************************ 
     */
    struct ReversePressureLikeOrdering {
      ReversePressureLikeOrdering()
      {};
      template <typename PressureType>
      bool operator()(const PressureType& value1,
                      const PressureType& value2) {

        if( value1 > value2 )
          return true;
        else
          return false;

      }
    }; // used by PositionIndexCompare



    struct HorizontalOrderingCartesian {
      template <typename V>
      bool operator()(const V& v1,
                      const V& v2) {
        assert(v1.dim() == v2.dim());
        for ( typename V::size_type i=v1.dim()-1; i>=0; i-- ) {
          if (v1[i] < v2[i])
            return true;
          else if (v1[i] > v2[i])
            return false;
        }
        return false;
      }
    };



    /* ************************************************************ *
     * Structure for comparing two vectors, 
     * such that highest preference is given to the first index
     * ************************************************************ */
    struct VerticalOrderingCartesian {
      template <typename V>    // <--- V = LocationAndVelocity
      bool operator()(const V& v1,
                      const V& v2) {

        assert(v1.first.dim() == v2.first.dim());
        // v1.first = x
        // v2.first = y
        // v1.second = beta(x)
        // v2.second = beta(y)

        for( UINT i=0;i<v1.first.dim();i++ ) {
          if (v1.first[i] < v2.first[i])
            return true;
          else if (v1.first[i] > v2.first[i])
            return false;
        }
        return false;
      }

    };

  } // Gesis

} // Dune

#endif // DUNE_GESIS_ORDERING_SCHEMES_HH
