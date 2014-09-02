#ifndef COARSE_GRID_P0_DATA_HANDLE_HH
#define COARSE_GRID_P0_DATA_HANDLE_HH

template<typename GV,typename ContainerType>
class CoarseGridP0Datahandle
  : public Dune::CommDataHandleIF<CoarseGridP0Datahandle<GV,ContainerType>,typename ContainerType::mapped_type>
{
private:
  const GV& gv;
  ContainerType& data_container;

public:
  // Note: 
  // 1.) 
  // For adaptive refinement with load balancing, 
  // we do not use IndexType idx = gv.indexSet().index(e) here,
  // because the index of a cell might get changed after 
  // adaptive mesh refinement + load balancing,
  // whereas the LocalIdSet keeps the numbering of unchanged elements.
  // 
  // LocalIdSet is good enough for the purpose of MPI communications.
  // Its structure is much simpler than that of the GlobalIdSet and it is just enough
  // for an element index to be locally unique.
  //
  // 2.)
  // For uniform refinement on YASP grid, the LocalIdSet index type
  // relies on bigunsignedint<k>::touint() which has produced reoccuring 
  // and therefore non-unique indices!!!
  // Since we are do not have to deal with load-balancing here, all coarse elements
  // will remain on the partition and won't get changed indices.
  // Here we can use gv.indexSet().index(e).
  // 
#ifdef USE_YASP
  typedef typename GV::IndexSet IdSet;
  typedef typename IdSet::IndexType Idx;
#else
  typedef typename GV::Grid::LocalIdSet IdSet; // reuired for ALUGRID
  typedef typename IdSet::IdType Idx;
#endif

  /* constructor */
  CoarseGridP0Datahandle( const GV& gv_,
                          ContainerType& data_container_ )
    : 
    gv(gv_),
    data_container(data_container_)
  {}
  
  bool contains(int dim, int codim) const{
    return (codim==0);
  }

  bool fixedsize(int dim, int codim) const{
    return true;
  }

  template<typename EntityType>
  size_t size( EntityType& e ) const{
    return 1;
  }

  /* Sender */
  template<typename MessageBuffer, typename EntityType>
  void gather( MessageBuffer& buff, const EntityType& e) const{

#ifdef USE_YASP
    Idx idx = gv.indexSet().index(e);
#else
    Idx idx = gv.grid().localIdSet().id(e);
#endif


#ifdef USE_YASP
    typename ContainerType::const_iterator it =
      data_container.find( idx );
#else
    typename ContainerType::const_iterator it =
      data_container.find( idx );
#endif
    if(it != data_container.end())
      buff.write(it->second);
  }

  /* Receiver */
  template<typename MessageBuffer, typename EntityType>
  void scatter( MessageBuffer& buff, const EntityType& e, size_t n){

#ifdef USE_YASP
    Idx idx = gv.indexSet().index(e);
#else
    Idx idx = gv.grid().localIdSet().id(e);
#endif
    REAL x;
    buff.read(x);
#ifdef USE_YASP
    data_container[ idx ] = x;
#else
    data_container[ idx ] = x;
#endif
  }

};




#endif // COARSE_GRID_P0_DATA_HANDLE_HH
