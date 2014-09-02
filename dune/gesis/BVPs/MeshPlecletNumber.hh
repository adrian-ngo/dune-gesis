#ifndef MESH_PECLET_NUMBER_HH
#define MESH_PECLET_NUMBER_HH

#include "../projectors/CoarseGridP0Datahandle.hh"

namespace Dune{
  namespace GeoInversion{

    template<typename GV,typename TP>
    class MeshPlecletNumber{

    private:
      const GV& gv;
      const TP& tp;
      REAL maxPeclet_l;
      REAL maxPeclet_t;
      REAL max_deltaSDFEM;
      std::vector<REAL> peclet_field;

    public:

      MeshPlecletNumber( const GV& gv_,
                         const TP& tp_,
                         const int step = 0 ) 
        :
        gv(gv_),
        tp(tp_),
        maxPeclet_l(0),
        maxPeclet_t(0),
        max_deltaSDFEM(0),
        peclet_field(gv.size(0)) {

        Dune::Timer watch;

#ifdef USE_YASP
        typedef typename GV::IndexSet IdSet;
        typedef typename IdSet::IndexType Idx;
#else
        typedef typename GV::Grid::LocalIdSet IdSet;  // required for ALUGRID
        typedef typename IdSet::IdType Idx;
#endif

        typedef typename GV::IndexSet IndexSet;
        typedef typename IndexSet::IndexType Index;

        typedef std::map<UINT,REAL> ContainerType;
        ContainerType data_container;

        // types and constants
        const int dim = GV::dimension;

        typedef typename GV::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
        //const typename GV::IndexSet& is = gv.indexSet();
  
        // loop over grid elements of codim 0
        for( ElementIterator eit=gv.template begin<0,Dune::All_Partition>()
               ; eit != gv.template end<0,Dune::All_Partition>(); ++eit) {

#ifdef USE_YASP
          Idx idx = gv.indexSet().index(*eit);
#else
          Idx idx = gv.grid().localIdSet().id(*eit);
#endif

          if(eit->partitionType()==Dune::InteriorEntity){

            Dune::FieldVector<CTYPE,dim> xglobal = eit->geometry().center();
            Dune::FieldVector<CTYPE,dim> xlocal = eit->geometry().local( xglobal );

            REAL Pe_l = 1.0;
            REAL Pe_t = 1.0;
            REAL L = 1.0; // characteristic length of grid element
            Dune::FieldVector<REAL,dim> beta;
            tp.getMeshPecletNumbers( *eit, xlocal, beta, L, Pe_l, Pe_t );

            if( eit->level() < gv.grid().maxLevel() )
              Pe_t = 0.0;
            
            data_container[idx] = Pe_t;

            // get maximal Mesh Peclet Number of all elements
            maxPeclet_l = std::max( maxPeclet_l, Pe_l );
            maxPeclet_t = std::max( maxPeclet_t, Pe_t );

            // get maximal SDFEM artifical diffusion
            max_deltaSDFEM = std::max( max_deltaSDFEM, tp.deltaSDFEM( *eit, xlocal ) );

          }

          else {

            data_container[idx] = 0;

          }

        }  // end of element loop



        if(gv.comm().size()>1){
          
          // MPI communicate
          logger << "MeshPecletNumber: Start data exchange ..." << std::endl;
          typedef CoarseGridP0Datahandle<GV,ContainerType> DataHandleType;
          DataHandleType datahandle(gv,data_container);
          //logger << "datahandle created." << std::endl;
    
          gv.communicate( datahandle, 
                          Dune::InteriorBorder_All_Interface, 
                          Dune::ForwardCommunication );
          //logger << "gv_gw.communicate() done." << std::endl;

        }


        // Re-assignment to Vector container
        for( ElementIterator eit=gv.template begin<0,Dune::All_Partition>()
               ; eit!=gv.template end<0,Dune::All_Partition>()
               ; ++eit) {

#ifdef USE_YASP
          Idx idx = gv.indexSet().index(*eit);
#else
          Idx idx = gv.grid().localIdSet().id(*eit);
#endif
          Index index = gv.indexSet().index(*eit);
          peclet_field[index] = data_container[idx];
        }


        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   General::verbosity,
                                   "REDIST",
                                   "Get maximal Mesh Peclet Number" );

      };



      void plot2vtu( const std::string filename ){
        VTKPlot::output_vector_to_vtu( gv, 
                                       peclet_field,
                                       filename,
                                       "Peclet"
                                       );

        
      };


      REAL maximum_l(){ return maxPeclet_l; };
      REAL maximum_t(){ return maxPeclet_t; };
      REAL maximum_deltaSDFEM(){ return max_deltaSDFEM; };
      
    }; // class MeshPlecletNumber
  }
}

#endif // MESH_PECLET_NUMBER_HH
