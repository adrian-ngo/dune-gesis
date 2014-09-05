#ifndef COARSE_GRID_P0_PROJECTOR_HH
#define COARSE_GRID_P0_PROJECTOR_HH

/*
  This class should work for adaptive refinement of gv_tp as well
  because we are using the check isLeaf() to filter out those
  elements in the level hierarchy which do not contribute!
*/


#include "CoarseGridP0Datahandle.hh"

namespace Dune{
  namespace GeoInversion{

template<typename GFS_TP,
         typename GFS_GW,
         typename IDT,
         int dim>
class CoarseGridP0Projector{

private:
  
  typedef Dune::FieldVector<REAL,dim> DomainType;

  typedef typename GFS_GW::Traits::GridViewType GV_GW;
  typedef typename GFS_TP::Traits::GridViewType GV_TP;

#ifdef USE_YASP
  typedef typename GV_GW::IndexSet IdSet;
  typedef typename IdSet::IndexType Idx;
#else
  typedef typename GV_GW::Grid::LocalIdSet IdSet;  // required for ALUGRID
  typedef typename IdSet::IdType Idx;
#endif

  typedef typename Dune::PDELab::BackendVectorSelector<GFS_TP,REAL>::Type VCType_TP;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;

  typedef Dune::PDELab::DiscreteGridFunction<GFS_TP,VCType_TP> DGF_TP;

  const GFS_TP& gfs_tp;
  const GFS_GW& gfs_gw;
  const IDT& inputdata;

  typedef std::map<UINT,REAL> ContainerType;

  ContainerType data_container;
  
public:  
  CoarseGridP0Projector( const GFS_TP& gfs_tp_,
                         const GFS_GW& gfs_gw_,
                         const IDT& inputdata_ ):
    gfs_tp(gfs_tp_),
    gfs_gw(gfs_gw_),
    inputdata(inputdata_)
  {
    /*
    if( gfs_gw.gridView().comm().rank()==0 && inputdata.verbosity > 3 ){
      std::cout << "CoarseGridP0Projector: source gfs size() = " 
                << gfs_tp.size() << std::endl;
      std::cout << "CoarseGridP0Projector: target gfs_gw.size() = " 
                << gfs_gw.size() << std::endl;
    }
    */
  }

  
  void apply( const VCType_TP& vc_tp, 
              const int maxGridLevel, 
              const int baselevel, 
              VCType_GW& vc_gw ){

    Dune::Timer watch;

    //std::cout << "DEBUG: maxGridLevel = " << maxGridLevel << std::endl;
    //std::cout << "DEBUG: baselevel = " << baselevel << std::endl;

    data_container.clear(); // clear contents leaving it with size 0

    DGF_TP dgf_tp(gfs_tp,vc_tp);
    
    typedef typename GV_GW::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
    typedef typename GV_GW::Traits::template Codim<0>::Entity ElementType;
    
    GV_GW gv_gw = gfs_gw.gridView();
    
    int iCounter=0;
    // Loop over elements of gv_gw on internal partions only and fill a data container of size gv_gw.size(0)
    for( ElementIterator eit=gv_gw.template begin<0,Dune::All_Partition>()
           ; eit!=gv_gw.template end<0,Dune::All_Partition>()
           ; ++eit) {

#ifdef USE_YASP
      Idx idx = gv_gw.indexSet().index(*eit);
#else
      Idx idx = gv_gw.grid().localIdSet().id(*eit);
#endif

      //std::cout << "DEBUG: cell center = " << (*eit).geometry().center() << std::endl;
      //std::cout << "DEBUG: idx = " << idx << std::endl;
      //std::cout << "DEBUG: idx.touint() = " << idx.touint() << std::endl;

      if(eit->partitionType()==Dune::InteriorEntity){

        REAL value = 0;

        if( maxGridLevel-baselevel == 0 ) {
          // case: no refinement at all
          const DomainType localcenter(0.5);
          Dune::FieldVector<REAL,1> contribution;
          dgf_tp.evaluate( *eit, localcenter, contribution );
          value = contribution[0];
          iCounter++;
        }
        else {
          // case: go through refined children
          // Loop over subgrid hierarchy
          for( int iLevel=1; iLevel<=maxGridLevel; iLevel ++ ){

            const typename ElementType::HierarchicIterator& hbegin 
              = eit->hbegin(iLevel);
            const typename ElementType::HierarchicIterator& hend 
              = eit->hend(iLevel);

            for (typename ElementType::HierarchicIterator hit = hbegin;
                 hit != hend; ++hit) {
              if ((*hit).isLeaf()) {
                //std::cout << "DEBUG: sub cell center = " 
                //          << (*hit).geometry().center() << std::endl;

                // Evaluate DG solution on the cellcenters of the refined grid!
                const DomainType localcenter(0.5);
                Dune::FieldVector<REAL,1> contribution;
                dgf_tp.evaluate( *hit, localcenter, contribution );
                value += contribution[0] * (*hit).geometry().volume();

                //std::cout << "DEBUG: sub cell volume = " 
                //          << (*hit).geometry().volume() << std::endl;
                //std::cout << "DEBUG: iLevel = " 
                //          << iLevel << std::endl;
                iCounter++;
              }
            }
          }
          value /= (*eit).geometry().volume();  // Take the volume averaged mean value of the DG solution.
        }

        data_container[idx] = value;
      }
      else {
        data_container[idx] = 0;
      }

    } // loop over elements of gv_gw

    General::log_elapsed_time( watch.elapsed(),
                               gv_gw.comm(),
                               General::verbosity,
                               "EVAL",
                               "CoarseGridP0Projector - part 1" );
    watch.reset();

    //std::cout << "DEBUG: iElementCounter = " << iCounter << std::endl;

    // MPI communicate
    logger << "Start data exchange." << std::endl;
    typedef CoarseGridP0Datahandle<GV_GW,ContainerType> DataHandleType;
    DataHandleType datahandle(gv_gw,data_container);
    //logger << "datahandle created." << std::endl;
    
    gv_gw.communicate( datahandle, 
                       Dune::InteriorBorder_All_Interface, 
                       Dune::ForwardCommunication );
    //logger << "gv_gw.communicate() done." << std::endl;

    General::log_elapsed_time( watch.elapsed(),
                               gv_gw.comm(),
                               General::verbosity,
                               "REDIST",
                               "CoarseGridP0Projector" );
    watch.reset();


    // Re-assignment to Vector container
    for( ElementIterator eit=gv_gw.template begin<0,Dune::All_Partition>()
           ; eit!=gv_gw.template end<0,Dune::All_Partition>()
           ; ++eit) {

#ifdef USE_YASP
      Idx idx = gv_gw.indexSet().index(*eit);
#else
      Idx idx = gv_gw.grid().localIdSet().id(*eit);
#endif
      Idx index = gv_gw.indexSet().index(*eit);
      
      //std::cout << "DEBUG: cell center = " << (*eit).geometry().center() << std::endl;
      //std::cout << "DEBUG: idx = " << idx << std::endl;
      //std::cout << "DEBUG: idx.touint() = " << idx.touint() << std::endl;
      //std::cout << "DEBUG: index = " << index << std::endl;

      vc_gw.base()[index] = data_container[idx];
    }

    data_container.clear(); // clear contents leaving it with size 0

    General::log_elapsed_time( watch.elapsed(),
                               gv_gw.comm(),
                               General::verbosity,
                               "EVAL",
                               "CoarseGridP0Projector - part 2" );
    
  }

};

  }
}
#endif // COARSE_GRID_P0_PROJECTOR_HH
