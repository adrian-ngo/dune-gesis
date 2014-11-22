#ifndef DUNE_GESIS_RGV_HH
#define DUNE_GESIS_RGV_HH

/*
 * Classes for renumbering grid cells by the pressure field
 *
 * Authors: J. Fahlke, E. Müller, A. Ngo
 * Date   : 12/2012
 *
 */

#include <iostream>
#include<list>
#include<vector>

#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/grid.hh>

#include <dune/pdelab/common/elementmapper.hh>
//#include "dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh"

#include "dune/gesis/common/general.hh"

#include "orderingschemes.hh"
#include "reorderedindexset.hh"

namespace Dune{

  namespace Gesis{




    /* ************************************************************ *
     * Grid view with reordered index set.
     *
     * The ordering of the index set is defined by the Compare class
     * ************************************************************ */
    template< typename GridView,
              typename RT0_PF,
              typename Compare >
    class ReorderedGridView : public GridView {

    private:
      Dune::shared_ptr<RT0_PF> rt0_pf;
      //const Compare& comparefunctor;
      Compare comparefunctor;

    public:

      /* ************************************************************ *
       *  subclass for reordered index set for elements and vertices.
       *  Reordering is realised with the (std::vector) mapping tables
       *
       *    reorderedindicesElements[] and
       *    reorderedindicesVertices[]
       *
       *  which return the reordered index for an index returned
       *  by the base class index function.
       * ************************************************************ */
      // Internal index set with vertical index running fastest
      typedef ReorderedIndexSet<typename GridView::Grid,
                                typename GridView::IndexSet,
                                RT0_PF
                                > IndexSet;
      typedef typename IndexSet::IndexType IndexType;

      Dune::shared_ptr<IndexSet> reorderedindexset;

      // Copy (from base class) constructor
      ReorderedGridView(const ReorderedGridView& other) :
        GridView(other.getGridView()),
        rt0_pf(other.get_rt0pf()),
        comparefunctor( other.getCompare() ),
        reorderedindexset( other.reorderedindexset ) {}


      // constructor of class ReorderedGridView
      ReorderedGridView( const GridView& other,
                         Dune::shared_ptr<RT0_PF> rt0_pf_,
                         const Compare& comparefunctor_
                         ) :
        GridView(other),
        rt0_pf(rt0_pf_),
        comparefunctor( comparefunctor_ ),
        reorderedindexset(new IndexSet( this->getGridView(),
                                        rt0_pf,
                                        comparefunctor ))
      {
      }

      ReorderedGridView & operator=(const ReorderedGridView& other){

        if(this != &other)
          GridView::operator=(other.getGridView());
        rt0_pf = other.get_rt0pf();
        comparefunctor = other.getCompare();
        reorderedindexset = other.reorderedindexset;

      }

      template<typename GFS_GW,typename VCType_GW>
      void update(
                  const GFS_GW& gfs_gw,
                  const VCType_GW& vc_h,
                  const GridView& other
                  ) {
#ifndef USE_DGF_PressureField
        rt0_pf->assemble(gfs_gw,vc_h);
#endif
        reorderedindexset->update( static_cast<const GridView&>(*this),
                                   rt0_pf,
                                   comparefunctor );
      }

      const IndexSet &indexSet () const {
        return *reorderedindexset;
      }

      const GridView& getGridView() const {
        return *static_cast<const GridView*>(this);
      }

      Dune::shared_ptr<RT0_PF> get_rt0pf() const {
        return rt0_pf;
        // return *static_cast<Dune::shared_ptr<const RT0_PF>> (rt0_pf);
      }

      Compare getCompare() const {
        return comparefunctor;
      }


    };

  } // Gesis
} // Dune
#endif
