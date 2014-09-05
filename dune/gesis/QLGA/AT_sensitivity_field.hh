#ifndef AT_SENSITIVITY_FIELD_HH
#define AT_SENSITIVITY_FIELD_HH

#include "dune/gesis/BVPs/projectors/CoarseGridP0Datahandle.hh"

namespace Dune{
  namespace GeoInversion{

    template<typename GV_GW, 
             typename VECTOR>
    void AT_sensitivity_field( // input:
                              const GV_GW& gv_gw        // the grid view
                              , const REAL & heatM0     // heat M0 at measuring point x=x_i
                              , const VECTOR & sensM0   // sensitivity of heat M0
                              , const REAL & heatM1     // heat M1 at measuring point x=x_i
                              , const VECTOR & sensM1   // sensitivity of heat M1
                              // result of this calculation:
                              , VECTOR& sensitivity				// the claculated sensitivity of the mean arrival time d(heatM0/heatM1)/dY	
                               ) {


#ifdef USE_YASP
      typedef typename GV_GW::IndexSet IdSet;
      typedef typename IdSet::IndexType Idx;
#else
      typedef typename GV_GW::Grid::LocalIdSet IdSet;  // required for ALUGRID
      typedef typename IdSet::IdType Idx;
#endif

      typedef typename GV_GW::IndexSet IndexSet;
      typedef typename IndexSet::IndexType Index;

      typedef std::map<UINT,REAL> ContainerType;
      ContainerType data_container;

      typedef typename GV_GW::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;

      //const typename GV_GW::IndexSet& is = gv_gw.indexSet();

      logger << "AT_sensitivity_field(): heatM0 = " << heatM0 << std::endl;
      logger << "AT_sensitivity_field(): heatM1 = " << heatM1 << std::endl;

      // loop over the grid
      for( ElementIterator eit=gv_gw.template begin<0,Dune::All_Partition>()
             ; eit!=gv_gw.template end<0,Dune::All_Partition>()
             ; ++eit) {
#ifdef USE_YASP
        Idx idx = gv_gw.indexSet().index(*eit);
#else
        Idx idx = gv_gw.grid().localIdSet().id(*eit);
#endif
        if(eit->partitionType()==Dune::InteriorEntity){
          data_container[idx] = sensM1[idx] / heatM0 -  heatM1 * sensM0[idx] / heatM0 / heatM0;
        } // end if(eit->partitionType()==Dune::InteriorEntity)
        else {
          data_container[idx] = 0;
        }
      }

      // MPI communicate
      logger << "mo_sensitivity_field(): Start data exchange ..." << std::endl;
      typedef CoarseGridP0Datahandle<GV_GW,ContainerType> DataHandleType;
      DataHandleType datahandle(gv_gw,data_container);
      //logger << "datahandle created." << std::endl;
    
      gv_gw.communicate( datahandle, 
                         Dune::InteriorBorder_All_Interface, 
                         Dune::ForwardCommunication );
      //logger << "gv_gw.communicate() done." << std::endl;
      
      
      // Re-assignment to Vector container
      for( ElementIterator eit=gv_gw.template begin<0,Dune::All_Partition>()
             ; eit!=gv_gw.template end<0,Dune::All_Partition>()
             ; ++eit) {
        
#ifdef USE_YASP
        Idx idx = gv_gw.indexSet().index(*eit);
#else
        Idx idx = gv_gw.grid().localIdSet().id(*eit);
#endif
        Index index = gv_gw.indexSet().index(*eit);
        sensitivity[index] = data_container[idx];
      }

    } //END: void AT_sensitivity_field()
      
  }

}

#endif // AT_SENSITIVITY_FIELD_HH
