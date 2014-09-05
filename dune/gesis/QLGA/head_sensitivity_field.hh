#ifndef HEAD_SENSITIVITY_FIELD_HH
#define HEAD_SENSITIVITY_FIELD_HH


/*

  The sensitivity field is calculated cell wise by integration over the
  same grid on which the forward and adjoint groundwater equation was solved.
  Note 1:
  The calculation need to be done on the Dune::InteriorEntity only.
  The ovelap entities get their values synchronized via communication
  using a datahandle. The information is stored in a map because this 
  guarantees a precise location of the data.
  Note 2:
  The groundwater equation might be solve on a globally refined grid with
  baselevel>0.
  But the sensitivity field is ALWAYS stored on the level 0 grid.
  It must have the same resolution as the parameter field Y.

*/

#include "dune/gesis/BVPs/projectors/CoarseGridP0Datahandle.hh"

namespace Dune {
  namespace GeoInversion {

    template<typename GV_GW,
             typename DARCY_BASE,
             typename DARCY_DGF,
             typename VECTOR,
             typename IDVECTOR
             >
    void head_sensitivity_field(
                                const GV_GW& gv_0
                                , const GV_GW& gv_gw                 // the grid view
                                , const DARCY_DGF& darcyflux_dgf   // the darcy flux dgf
                                , const DARCY_BASE& gradient_dgf // the gradient of the adjoint
                                , const int qorderP_Q1             // integration order
                                // result of this calculation:
                                , VECTOR& sensitivity // the calculated sensitivity
                                , VECTOR& iSensitivity // only interior partition!
                                , IDVECTOR& iGlobalGridIndexVector // only interior partition!
                                ){
      Dune::Timer watch;


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

      // types and constants
      typedef typename GV_GW::Grid::ctype DF;
      const int dim = GV_GW::dimension;

      typedef typename GV_GW::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;


      Idx idx;
      // loop over the Level0-grid
      for( ElementIterator e0it=gv_0.template begin<0,Dune::All_Partition>()
             ; e0it!=gv_0.template end<0,Dune::All_Partition>()
             ; ++e0it) {
#ifdef USE_YASP
        Idx idx = gv_0.indexSet().index(*e0it);
#else
        Idx idx = gv_0.grid().localIdSet().id(*e0it);
#endif
        data_container[ idx ] = 0.0; // initialize with zero
      }


      // loop over the grid
      for( ElementIterator eit=gv_gw.template begin<0,Dune::All_Partition>()
             ; eit!=gv_gw.template end<0,Dune::All_Partition>()
             ; ++eit) {
    
#ifdef USE_YASP
        idx = gv_gw.indexSet().index(*eit);
#else
        idx = gv_gw.grid().localIdSet().id(*eit);
#endif

        if( eit->level() > 0 ){
          typedef typename GV_GW::Grid::template Codim<0>::EntityPointer ElementPointer;
          ElementPointer pAncestor = eit->father();
          while( pAncestor->level() > 0 )
            pAncestor = (*pAncestor).father();
#ifdef USE_YASP
          idx = gv_gw.indexSet().index(*pAncestor);
#else
          idx = gv_gw.grid().localIdSet().id(*pAncestor);
#endif
        }


        if(eit->partitionType()==Dune::InteriorEntity){

          REAL value = 0;

          // select quadrature rule
          Dune::GeometryType geometrytype = eit->geometry().type();
          const int integration_order = qorderP_Q1;
          const Dune::QuadratureRule<DF, dim>& rule = Dune::QuadratureRules<DF, dim>::rule(geometrytype, integration_order);
      
          // quadrature points
          Dune::FieldVector<REAL, dim> minus_k_gradh_QP(0.0);
          Dune::FieldVector<REAL, dim> gradpsi_QP(0.0);

          // loop over quadrature points 
          for( typename Dune::QuadratureRule<DF, dim>::const_iterator iterQP = rule.begin()
                 ; iterQP != rule.end()
                 ; ++iterQP) {

            // Floating Point Exception in 3D parallel is caused here, when exception handler is activated!
        
            // evaluate the darcy flux
            minus_k_gradh_QP = 0.0;
            darcyflux_dgf.evaluate(*eit, iterQP->position(), minus_k_gradh_QP);
      
            //evaluate the gradient of the adjoint
            gradpsi_QP = 0.0;
            gradient_dgf.evaluate(*eit, iterQP->position(), gradpsi_QP);

            // integration factor
            REAL factor = iterQP->weight() * eit->geometry().integrationElement(iterQP->position());
            
            // add to the sensitivity
            value += minus_k_gradh_QP * gradpsi_QP * factor;

          } // quadrature loop

          if( eit->level() > 0 ){
            data_container[idx] += value;
          }
          else{
            data_container[idx] += value;
          }

          iSensitivity.push_back( data_container[ idx ] );
          auto gidx = gv_gw.grid().globalIdSet().id(*eit);
          /*
          std::cout << "QQ" << gv_gw.comm().rank()
                    << " >> gidx: " << gidx
                    << " >> size: " << sizeof(gidx)
                    << std::endl;
          */
          iGlobalGridIndexVector.push_back( gidx );

        } // end if(eit->partitionType()==Dune::InteriorEntity)
        else {
          data_container[idx] = 0;
        }


      } // element loop

      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 General::verbosity,
                                 "EVAL",
                                 "head_sensitivity_field - part 1" );
      watch.reset();

      // MPI communicate
      logger << "head_sensitivity_field(): Start data exchange ..." << std::endl;
      typedef CoarseGridP0Datahandle<GV_GW,ContainerType> DataHandleType;
      DataHandleType datahandle(gv_0,data_container);
      //logger << "datahandle created." << std::endl;
    
      gv_0.communicate( datahandle,
                        Dune::InteriorBorder_All_Interface,
                        Dune::ForwardCommunication );
      //logger << "gv_gw.communicate() done." << std::endl;

      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 General::verbosity,
                                 "REDIST",
                                 "head_sensitivity_field" );
      watch.reset();

      // Re-assignment to Vector container
      for( ElementIterator e0it=gv_0.template begin<0,Dune::All_Partition>()
             ; e0it!=gv_0.template end<0,Dune::All_Partition>()
             ; ++e0it) {

#ifdef USE_YASP
        Idx idx = gv_0.indexSet().index(*e0it);
#else
        Idx idx = gv_0.grid().localIdSet().id(*e0it);
#endif
        Index index = gv_0.indexSet().index(*e0it);
        sensitivity[index] = data_container[ idx ];
      }

      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 General::verbosity,
                                 "EVAL",
                                 "head_sensitivity_field  - part 2" );

    } // void head_sensitivity_field()

  } // GeoInversion

} // Dune


#endif // HEAD_SENSITIVITY_FIELD_HH
