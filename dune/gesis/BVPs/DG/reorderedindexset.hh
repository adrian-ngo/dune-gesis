#ifndef DUNE_GESIS_REORDERED_INDEX_SET_HH
#define DUNE_GESIS_REORDERED_INDEX_SET_HH

namespace Dune{

  namespace Gesis{

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

  } // Gesis

} // Dune

#endif // DUNE_GESIS_REORDERED_INDEX_SET_HH

