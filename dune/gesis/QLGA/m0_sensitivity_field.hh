#ifndef DUNE_GESIS_M0_SENSITIVITY_FIELD_HH
#define DUNE_GESIS_M0_SENSITIVITY_FIELD_HH


#include "dune/gesis/BVPs/projectors/CoarseGridP0Datahandle.hh"

namespace Dune {
  namespace Gesis {

    template<typename GV_GW,
             typename DGF,
             typename DARCY_BASE,
             typename DARCY_DGF,
             typename DGF_GRADIENT,
             typename VECTOR>
    void m0_sensitivity_field(
                              const GV_GW& gv_0  // level 0 grid view
                              , const GV_GW& gv_gw // baselevel grid view
                              , const DGF& m0_adj_dgf // dgf of m0 adjoint
                              , const DARCY_DGF& darcyflux_dgf // dgf of the darcy flux
                              , const DARCY_BASE& gradient_h_adj_m0_dgf // dgf of: head adjoint
                              , const DGF_GRADIENT& gradientM0_dgf // dgf of: m0 adjoint
                              , const int maxGridLevel
                              , const int baselevel
                              , const int qorderP_Q1 // integration order
                              // result of this calculation:
                              , VECTOR& sensitivity
                              , VECTOR& iSensitivity
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
      typedef typename GV_GW::Traits::template Codim<0>::Entity ElementType;


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
        Idx idx = gv_gw.indexSet().index(*eit);
#else
        Idx idx = gv_gw.grid().localIdSet().id(*eit);
#endif

        if( eit->level() > 0 ){
          typedef typename GV_GW::Grid::template Codim<0>::EntityPointer ElementPointer;
          ElementPointer pAncestor = eit->father();
          while( pAncestor->level() > 0 )
            pAncestor = (*pAncestor).father();
#ifdef USE_YASP
          idx = gv_0.indexSet().index(*pAncestor);
#else
          idx = gv_0.grid().localIdSet().id(*pAncestor);
#endif
        }

        if(eit->partitionType()==Dune::InteriorEntity){

          // variables for the evaluation at the quadrature points:
          REAL value = 0;
          Dune::FieldVector<REAL,dim> minus_k_gradh_QP(0.0);
          Dune::FieldVector<REAL,dim> gradpsi_h_QP(0.0);
          Dune::FieldVector<REAL,dim> gradm0_QP(0.0);
          Dune::FieldVector<REAL,1> psi_m0_QP(0.0);

          const int integration_order = qorderP_Q1;

          // Integration on the coarse level:
          // select quadrature rule
          Dune::GeometryType gt = eit->geometry().type();
          const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,integration_order);

          // loop over quadrature points
          for( typename Dune::QuadratureRule<DF, dim>::const_iterator iterQP = rule.begin()
                 ; iterQP != rule.end()
                 ; ++iterQP) {

            // evaluate the darcy flux
            minus_k_gradh_QP = 0.0;
            darcyflux_dgf.evaluate(*eit, iterQP->position(), minus_k_gradh_QP);
      
            // evaluate gradient of h adjoint
            gradpsi_h_QP = 0.0;
            gradient_h_adj_m0_dgf.evaluate(*eit, iterQP->position(), gradpsi_h_QP);

            // integration factor
            REAL factor = iterQP->weight() * eit->geometry().integrationElement(iterQP->position());

            value -= (minus_k_gradh_QP * gradpsi_h_QP ) * factor ;

            if(maxGridLevel==baselevel) {

              // evaluate gradient of m0
              gradm0_QP= 0.0;
              gradientM0_dgf.evaluate(*eit, iterQP->position(), gradm0_QP);

              // evaluate m0
              psi_m0_QP = 0.0;
              m0_adj_dgf.evaluate(*eit, iterQP->position(), psi_m0_QP);

              //value -= ((gradm0_QP*minus_k_gradh_QP)*psi_m0_QP) * factor;
              value += ((gradm0_QP*minus_k_gradh_QP)*psi_m0_QP) * factor;
            }

          } // end of coarse quadratur loop


          if( maxGridLevel > baselevel ) {

            // Integration on the fine levels:
            // Loop over subgrid hierarchy
            for( int iLevel=1; iLevel<=maxGridLevel; iLevel ++ ){

              const typename ElementType::HierarchicIterator& hbegin 
                = eit->hbegin(iLevel);
              const typename ElementType::HierarchicIterator& hend 
                = eit->hend(iLevel);
            
              for (typename ElementType::HierarchicIterator hit = hbegin;
                   hit != hend; ++hit) {
                if ((*hit).isLeaf()) {

                  Dune::GeometryType hgt = hit->geometry().type();
                  const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(hgt,integration_order);

                  // loop over quadrature points of the subcell *hit
                  for( typename Dune::QuadratureRule<DF,dim>::const_iterator itQP = rule.begin()
                         ; itQP != rule.end()
                         ; ++itQP) {
                  
                    // evaluate the darcy flux on root
                    minus_k_gradh_QP = 0.0;
                    darcyflux_dgf.evaluate_on_root(*hit, itQP->position(), minus_k_gradh_QP);
      
                    // evaluate gradient of h adjoint on root
                    gradpsi_h_QP = 0.0;
                    gradient_h_adj_m0_dgf.evaluate_on_root(*hit, itQP->position(), gradpsi_h_QP);

                    // evaluate gradient of m0
                    gradm0_QP= 0.0;
                    gradientM0_dgf.evaluate(*hit, itQP->position(), gradm0_QP);

                    // evaluate m0
                    psi_m0_QP = 0.0;
                    m0_adj_dgf.evaluate(*hit, itQP->position(), psi_m0_QP);
            
                    // integration factor
                    REAL factor = itQP->weight() * hit->geometry().integrationElement(itQP->position());
                    // add the data to the sensitivity
                    // Remark: without division by porosity because we use q in the advection-dispersion and not the seepage velocity (v)

                    //value -= ((gradm0_QP*minus_k_gradh_QP)*psi_m0_QP) * factor ;
                    value += ((gradm0_QP*minus_k_gradh_QP)*psi_m0_QP) * factor ;
                  
                  } // quadrature loop

                } // if isLeaf

              } // loop over refined grid

            } // loop over refinement levels

          } // end if( maxGridLevel > baselevel )

          if( eit->level() > 0 ){
            data_container[idx] += value;
          }
          else{
            data_container[idx] += value;
          }

          iSensitivity.push_back( data_container[ idx ] );

        } // end if(eit->partitionType()==Dune::InteriorEntity)
        else {
          data_container[idx] = 0;
        }
        
      } // end of loop over the grid


      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 General::verbosity,
                                 "EVAL",
                                 "m0_sensitivity_field - part 1" );
      watch.reset();

      // MPI communicate
      logger << "mo_sensitivity_field(): Start data exchange ..." << std::endl;
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
                                 "m0_sensitivity_field" );
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
        sensitivity[index] = data_container[idx];

      }

      General::log_elapsed_time( watch.elapsed(),
                                 gv_0.comm(),
                                 General::verbosity,
                                 "EVAL",
                                 "m0_sensitivity_field - part2" );

    } // void m0_sensitivity_field()





    template<typename GV_GW,
             typename DARCY_BASE,
             typename DARCY_DGF,
             typename GFS_TP,
             typename VCType_TP,
             typename V>
    void get_m0_Sensitivities(
                              const GV_GW& gv_0,  // sensitivity grid level = level 0
                              const GV_GW& gv_gw, // head grid level = baselevel
                              const DARCY_DGF& darcyflux_dgf,
                              const DARCY_BASE& grad_hadj_m0_dgf,
                              const GFS_TP& gfs_tp,
                              const VCType_TP& vc_m0adj,
                              const VCType_TP& vc_m0,
                              const int maxGridLevel,
                              const int baselevel,
                              const int pOrder,
                              V& solution,
                              V& iSolution
                              ){
    
      typedef Dune::PDELab::DiscreteGridFunction<GFS_TP,VCType_TP> DGF_TP;

      DGF_TP m0adj_dgf( gfs_tp, vc_m0adj );

      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_TP,VCType_TP> GRADIENT_DGF_TP;
      GRADIENT_DGF_TP grad_m0_dgf( gfs_tp, vc_m0 );
    
      m0_sensitivity_field( gv_0
                            , gv_gw
                            , m0adj_dgf
                            , darcyflux_dgf
                            , grad_hadj_m0_dgf
                            , grad_m0_dgf
                            , maxGridLevel
                            , baselevel
                            , pOrder 
                            , solution
                            , iSolution
                            );

    }


  } // namespace Gesis

} // namespace Dune
#endif // DUNE_GESIS_M0_SENSITIVITY_FIELD_HH
