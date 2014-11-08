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


namespace Dune{

  namespace Gesis{

template<typename GV, typename IDT>
void outputGridviewIndexToDGF( const GV& gv,
                               const IDT& inputdata,
                               const std::string filename,
                               const std::string meshtype="cube"
                               ){
  Dune::Timer watch;

  // (1) Construct finite element space
  const int dim = GV::dimension;
  Dune::GeometryType gt;

  if( meshtype == "cube" )
    gt.makeCube(dim);
  else
    gt.makeSimplex(dim);
  
  typedef Dune::PDELab::P0LocalFiniteElementMap<double,double,dim> P0FEM;
  P0FEM p0fem( gt );

  // (2) Set up grid function space
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;

  typedef Dune::PDELab::NoConstraints NOCON;
  NOCON nocon;

  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,NOCON,VBE> P0GFS;
  P0GFS p0gfs(gv, p0fem, nocon );

  typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,double>::Type UType;    
  UType u0(p0gfs, 0.0);


  for (unsigned int i=0;i<u0.flatsize();i++) {
    u0.base()[i] = REAL(i)/(u0.flatsize()-1);
    // plot the processor rank for debugging purposes:
    // u0[i] = gv.comm().rank();
  }

  /*
  // Alternative loop implementation:
  typedef typename GV::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
  int iElement = 0;
  for (ElementIterator it=gv.template begin<0,Dune::All_Partition>();
       it!=gv.template end<0,Dune::All_Partition>();++it) {

    // plotting element level
    //u0[gv.indexSet().index(*it)] = (*it).level(); 
    // plotting element index
    u0[gv.indexSet().index(*it)] = gv.indexSet().index(*it);
    iElement++;
  }
  */
    
  typedef Dune::PDELab::DiscreteGridFunction<P0GFS,UType> P0DGF;
  P0DGF u0dgf(p0gfs,u0);

  std::stringstream jobtitle;
  jobtitle << "outputGridviewIndexToDGF: copying ";
  General::log_elapsed_time( watch.elapsed(),
                             gv.comm(),
                             General::verbosity,
                             "IO",
                             jobtitle.str() );

  Dune::Gesis::VTKPlot::output2vtu( //gv, 
                                           p0gfs, 
                                           u0, 
                                           filename, 
                                           "elementorder", 
                                           General::verbosity, 
                                           true, 
                                           0
                                           );

}


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


/* ************************************************************ *
 * Structure for comparing two pairs
 * ************************************************************ */
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
    // To the linear solver is does not matter whether we are sorting
    // in ascending (against the flow) or descending (in flow direction)
    // order, because BiCGStab is designed to cope with 
    // non-sysmmetric matrices!
    //
    if( value1 < value2 )
      return true;
    else
      return false;

  }
}; // used by PositionIndexCompare



/* ************************************************************ *
 * Structure for comparing two pairs
 * ************************************************************ */
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



template<typename Grid, 
         typename BackendIndexSet, 
         typename RT0_PF
         >
class ReorderedIndexSet :
  public Dune::IndexSet<Grid,
                        ReorderedIndexSet<Grid, 
                                          BackendIndexSet, 
                                          RT0_PF>,
                        typename BackendIndexSet::IndexType>
{
  typedef BackendIndexSet Backend;
public:
  // Type of base class index type
  typedef typename Backend::IndexType IndexType;
private:
  typedef Dune::IndexSet<Grid, ReorderedIndexSet, IndexType> Base;
  // Extract dimension
  enum {dim = Grid::dimension };   

  template<typename PressureType, typename Compare>
  struct PositionIndexCompare {
    PositionIndexCompare(const Compare& compare_) : compare(compare_) {}
    bool operator()(const typename std::pair<PressureType,IndexType>& v1i,
                    const typename std::pair<PressureType,IndexType>& v2i) {
      return compare(v1i.first,v2i.first);
    }
    Compare compare;
  };



  // Type of mapping tables
  typedef std::vector<IndexType> ReorderedIndices;

  // Generate the mapping table for a given codim.
  template<int codim, typename GridView, typename Compare>
  void generateReorderedIndices (const GridView& gridview,
                                 ReorderedIndices& reorderedindices,
                                 Dune::shared_ptr<RT0_PF> rt0_pf,
                                 const Compare &compare)
  {
    //std::cout << "DEBUG: Generating reordered indices for codim = " << codim << std::endl;
    //typedef typename GridView::template Codim<0>::Geometry::GlobalCoordinate CoordType;

    typedef REAL PressureType;
    typedef typename std::pair<PressureType,IndexType> PositionIndex;
    //typedef typename std::list<PositionIndex> IndexList;
    typedef typename std::vector<PositionIndex> IndexList;

    // Loop over all entities it and generate a list of pairs 
    // {v(it),index(it)} where v(it) is the center of the entity
    // and index(it) its index 
    IndexList indexlist;

    // Map each cell to unique id
    Dune::PDELab::ElementMapper<GridView> cell_mapper(gridview);


    if( codim==dim ) {
      // std::cout << "DEBUG: Generating reordered indices for codim = " << codim << std::endl;
      typedef typename GridView::template Codim<dim>::template Partition<Dune::All_Partition>::Iterator Iterator;
      typedef typename GridView::template Codim<dim>::Geometry::GlobalCoordinate CoordType;
      for (Iterator it=gridview.template begin<dim,Dune::All_Partition>();
           it!=gridview.template end<dim,Dune::All_Partition>();++it) {

        CoordType xglobal(0.0);
        xglobal = it->geometry().corner(0);

        IndexType idx = gridview.indexSet().index(*it);
        PressureType pressure = (PressureType) idx; // vertices order remain unchanged!
        indexlist.push_back( std::make_pair(pressure,idx) );
      }
    }
    else if( codim==0 ) {
      // std::cout << "DEBUG: Generating reordered indices for codim = " << codim << std::endl;

      typedef typename GridView::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator Iterator;
      typedef typename GridView::Traits::template Codim<0>::Geometry::GlobalCoordinate CoordType;

      for (Iterator it=gridview.template begin<0,Dune::All_Partition>();
           it!=gridview.template end<0,Dune::All_Partition>();++it) {

        CoordType xlocal(0.5); // evaluate pressure on the center of the element
        Dune::FieldVector<REAL,1> rt0head(0.0);
        rt0_pf->evaluate_on_root( *it, xlocal, rt0head );

        PressureType pressure = rt0head[0]; // <-- we use pressure reconstructed from RT0 velocity
        
        IndexType idx = gridview.indexSet().index(*it);
        indexlist.push_back(std::make_pair(pressure,idx));
      } // for
    } // if(codim==0)

    PositionIndexCompare<PressureType,Compare> comparePairs(compare);

    // Sort list according to first entry in each pair
    // indexlist.sort(comparePairs);
    std::sort( indexlist.begin(), indexlist.end(), comparePairs );

    UINT indexListSize = indexlist.size();
    reorderedindices.resize( indexListSize );

    int i=0;
    // Generate the mapping table using the index information in the
    // reordered list.
    for(typename IndexList::const_iterator it=indexlist.begin();
        it!=indexlist.end();it++) {
      reorderedindices[it->second] = i;  // new index [old index]
      i++;
    }
  }

public:
  // Constructor. Initialise with a given gridview
  ReorderedIndexSet(const ReorderedIndexSet& other) :
    backend(other.backend)
    , reorderedindicesElements(other.reorderedindicesElements)
    , reorderedindicesVertices(other.reorderedindicesVertices)
  {} 

  /*
  ReorderedIndexSet & operator=(const ReorderedIndexSet& other){
    backend = other.backend;
    reorderedindicesElements = other.reorderedindicesElements;
    reorderedindicesVertices = other.reorderedindicesVertices;
  }
  */

  // called from constructor of class ReorderedGridView:
  template<typename GridView, typename Compare>
  ReorderedIndexSet( const GridView& gridview, 
                     Dune::shared_ptr<RT0_PF> rt0_pf,
                     const Compare &compare )
    : backend(&gridview.indexSet())
  {
    Dune::Timer watch;

    generateReorderedIndices<0>(gridview,
                                reorderedindicesElements,
                                rt0_pf,
                                compare);

    generateReorderedIndices<dim>(gridview,
                                  reorderedindicesVertices,
                                  rt0_pf,
                                  compare);
    
    std::stringstream jobtitle;
    jobtitle << "Reordering elements on initial grid.";
    General::log_elapsed_time( watch.elapsed(),
                               gridview.comm(),
                               General::verbosity,
                               "RGV",
                               jobtitle.str() );
  }


  template<typename GridView, typename Compare>
  void update( const GridView& gridview, 
               Dune::shared_ptr<RT0_PF> rt0_pf,
               const Compare &compare
               ) {
    Dune::Timer watch;

    backend = &gridview.indexSet();
    generateReorderedIndices<0>(gridview,
                                reorderedindicesElements,
                                rt0_pf,
                                compare);

    generateReorderedIndices<dim>(gridview,
                                  reorderedindicesVertices,
                                  rt0_pf,
                                  compare);

    std::stringstream jobtitle;
    jobtitle << "Reordering elements on updated grid.";
    General::log_elapsed_time( watch.elapsed(),
                               gridview.comm(),
                               General::verbosity,
                               "RGV",
                               jobtitle.str() );
  }



  template<int cc>
  IndexType index (const typename Dune::remove_const<Grid>::type::
                   Traits::template Codim<cc>::Entity& e) const
  {
    //std::cout << " this->index(e) = " 
    //          << this->index(e) 
    //          << std::endl;
    return this->index(e);
  }

  // Return index of entity
  template<typename Entity>
  IndexType index (const Entity& e) const {
    int cc = Entity::codimension;
    assert( cc == 0 || cc == Backend::dimension );
    if( cc==0 ) {

      //std::cout << " new index = " 
      //          << reorderedindicesElements[backend->index(e)]
      //          << " old index = "
      //          << backend->index(e)
      //          << std::endl;

      return reorderedindicesElements[backend->index(e)]; // new index [old index]

    }

    else if( cc==dim ) {
      return reorderedindicesVertices[backend->index(e)]; 
    }

    else{
        DUNE_THROW(Dune::RangeError,
                   "rgv.index() for codim " << Entity::codimension
                   << " should never be "
                   "called.");
        return 0;
    }
        
  }

  template< int cc >
  IndexType subIndex (
                      const typename Dune::remove_const<Grid>::type::
                      Traits::template Codim<cc>::Entity &e,
                      //const typename Backend::Traits::template Codim< cc >::Entity &e,
                      int i, unsigned int codim ) const
  {
    std::cout << "RGV: cc this->subIndex(e," << i << "," << codim << ") = " 
              << this->subIndex(e, i, codim) 
              << std::endl;

    return this->subIndex(e, i, codim);
  }

  // Return index of subentity
  template<typename Entity>
  IndexType subIndex (const Entity &e,
                      int i, unsigned int codim ) const {
    int cc = Entity::codimension;
    assert( cc == 0 || cc == Backend::dimension );
    int cc_codim = cc+codim;
    if ( cc_codim == 0 )
      return reorderedindicesElements[backend->subIndex(e,i,codim)]; 

    else if (cc_codim == dim) 
      return reorderedindicesVertices[backend->subIndex(e,i,codim)];

    else {
      DUNE_THROW(Dune::RangeError,
                 "subIndex for entity with codim " << Entity::codimension
                 << " should never be "
                 "called.");
      std::cout << "RGV: ReorderedIndexSet ERROR" << std::endl;
    }

    std::cout << "RGV: Entity this->subIndex(e," << i << "," << codim << ") = " 
              << this->subIndex(e, i, codim) 
              << std::endl;

    exit(-1);
  }

  const std::vector<Dune::GeometryType>& geomTypes (int codim) const
  {
    return backend->geomTypes(codim);
  }

  IndexType size (Dune::GeometryType type) const
  {
    return backend->size(type);
  }

  IndexType size (int codim) const
  {
    return backend->size(codim);
  }
  
  template<typename EntityType>
  bool contains (const EntityType& e) const
  {
    return backend->contains(e);
  }

private:
  const Backend* backend;

  enum { ncodim = dim+1 };
  ReorderedIndices reorderedindicesVertices;
  ReorderedIndices reorderedindicesElements;
};


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
