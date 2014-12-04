#ifndef DUNE_GESIS_RT0_PRESSUREFIELD_HH
#define DUNE_GESIS_RT0_PRESSUREFIELD_HH


#define USE_LOCAL_COORD


#include<dune/grid/common/datahandleif.hh>

extern CLogfile logger;


namespace Dune{
  namespace Gesis{

    template<
      typename GV,
      typename ContainerType
      >
    class RT0_PressureField_Exchange
      : public Dune::CommDataHandleIF< RT0_PressureField_Exchange<GV,ContainerType>, typename ContainerType::mapped_type >
    {
    private:
      const GV& gv;
      ContainerType& data_container;
      const int blocksize;

    public:
      // Note:
      // 1.) We do not use IndexType idx = gv.indexSet().index(e) here,
      // because the index of a cell might get changed after local mesh refinement + load balancing,
      // whereas the LocalIdSet keeps the numbering of unchanged elements.
      //
      // 2.) LocalIdSet is good enough for the purpose of MPI communications.
      // Its structure is much simpler than that of the GlobalIdSet and it is just enough
      // for an element index to be locally unique.
      //

#ifndef USE_ALUGRID
      typedef typename GV::IndexSet IdSet;
      typedef typename IdSet::IndexType Idx;
#else
      typedef typename GV::Grid::LocalIdSet IdSet; // reuired for ALUGRID
      typedef typename IdSet::IdType Idx;
#endif

      /* constructor */
      RT0_PressureField_Exchange(
                                 const GV& gv_,
                                 ContainerType& data_container_,
                                 const int blocksize_
                                  )
        :
        gv(gv_),
        data_container(data_container_),
        blocksize(blocksize_)
      {
      }

      bool contains(int dim, int codim) const{
        return (codim==0);
      }

      bool fixedsize(int dim, int codim) const{
        return true;
      }

      template<typename EntityType>
      size_t size( EntityType& e ) const{
        return blocksize;
      }

      /* Sender */
      template<typename MessageBuffer, typename EntityType>
      void gather( MessageBuffer& buff, const EntityType& e) const{

#ifndef USE_ALUGRID
        Idx idx = gv.indexSet().index(e);
#else
        Idx idx = gv.grid().localIdSet().id(e);
#endif

        for(int i=0;i<blocksize;i++){
          typename ContainerType::const_iterator it =
            data_container.find( blocksize*idx+i );
          if(it != data_container.end())
            buff.write(it->second);
        }
      }

      /* Receiver */
      template<typename MessageBuffer, typename EntityType>
      void scatter( MessageBuffer& buff, const EntityType& e, size_t n){

#ifndef USE_ALUGRID
        Idx idx = gv.indexSet().index(e);
#else
        Idx idx = gv.grid().localIdSet().id(e);
#endif

        for(int i=0;i<blocksize;i++){
          REAL x;
          buff.read(x);
          data_container[ blocksize*idx+i ] = x;
        }
      }

    };


    template<
      typename GV,
      typename GFS_GW,
      typename GWP
      >
    class RT0_PressureField {

      typedef typename GV::Grid GRID;
      typedef typename GV::ctype DF;

#ifndef USE_ALUGRID
      typedef typename GV::IndexSet IdSet;
      typedef typename IdSet::IndexType Idx;
#else
      typedef typename GV::Grid::LocalIdSet IdSet;  // required for ALUGRID
      typedef typename IdSet::IdType Idx;
#endif

      typedef typename Dune::PDELab::BackendVectorSelector<GFS_GW,REAL>::Type VCType_GW;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_GW,VCType_GW> DGF;

      typedef Dune::Gesis::VectorDGF<GWP,GFS_GW> DARCY_FLUX_DGF;


      typedef typename GV::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection Intersection;

      //typedef std::map<Dune::bigunsignedint<58>,REAL> ContainerType; // using a map is better than using a vector because idx is not necessarily successive

      //typedef std::map<UINT,REAL> ContainerType;
      typedef std::map<Idx,REAL> ContainerType;


      typedef typename GV::Traits::template Codim<0>::Entity ElementType;

    private:
      enum{ dim = GV::dimension };

      GV& gv; // gv may be updated after load balancing
      const GWP& gwp;
      const UINT baselevel;
      const int blocksize;
      ContainerType data_container;


      void appendNewElement(
                            const ElementType& e,
                            const Vector<REAL>& coefficients ){


#ifndef USE_ALUGRID
        Idx idx = gv.indexSet().index(e);
#else
        Idx idx = gv.grid().localIdSet().id(e);
#endif

        logger << "decltype(idx)=" << General::type_name<decltype(idx)>() << std::endl;

        for(UINT i=0;i<coefficients.size();i++)
          data_container[blocksize*idx+i] = coefficients[i];
      }


    public:
      RT0_PressureField( GV& gv_, // gv may be updated after load balancing
                         const GWP& gwp_,
                         const UINT baselevel_=0 )
        : gv(gv_),
          gwp(gwp_),
          baselevel(baselevel_),
          blocksize(2*dim+1)
      {
      }



      void assemble( const GFS_GW& gfs_gw_,
                     const VCType_GW& vc_h_ ){

        gv = gfs_gw_.gridView(); // gv is being updated here possibly after load balancing
        data_container.clear();

        Dune::Timer watch;

        watch.reset();

        if(gv.comm().rank()==0)
          std::cout << "Reconstruct rt0_pressurefield..." << std::endl;

        const DGF dgf(gfs_gw_,vc_h_);
        const DARCY_FLUX_DGF darcyflux_dgf( gwp, gfs_gw_, vc_h_, baselevel );

        //logger << "RT0_pf: Start element loop." << std::endl;
        // element loop
        for( ElementIterator eit=gv.template begin<0,Dune::All_Partition>()
               ; eit!=gv.template end<0,Dune::All_Partition>()
               ; ++eit) {

          // Evaluate pressure on cellcenter
          Dune::FieldVector<REAL,dim> insideCellCenterLocal =
            Dune::ReferenceElements<DF,dim>::general(eit->type()).position(0,0);

          Dune::FieldVector<REAL,dim> insideCellCenterGlobal =
            eit->geometry().global(insideCellCenterLocal);

          //logger << "insideCellCenterGlobal = " << insideCellCenterGlobal << std::endl;
          //logger << "insideCellCenterLocal = " << insideCellCenterLocal << std::endl;

          Dune::FieldVector<REAL,1> pCellcenter(0.0);
          dgf.evaluate( *eit, insideCellCenterLocal, pCellcenter );

          typename GWP::Traits::DiffusionTensorType tensor_inside(gwp.DiffusionTensor(*eit,insideCellCenterLocal));
          Vector<REAL> coeff( blocksize );

          if(eit->partitionType()==Dune::InteriorEntity){

            Dune::Gesis::DenseMatrix<REAL> XY( blocksize, blocksize, 0.0 );
            Vector<REAL> rhs( blocksize );

            // It suffices to traverse only two internal of the four intersections of element *eit
            // This is enough to get 2*dim equations for 2*dim unknows
            //
            int faceCounter=0;

            //intersection loop:
            for( IntersectionIterator iit = gv.ibegin(*eit)
                   ; iit!=gv.iend(*eit)
                   ; ++iit )
              {

                //std::cout << "face " << faceCounter << std::endl;

                // face geometry
                Dune::FieldVector<REAL,dim-1> faceLocal =
                  Dune::ReferenceElements<REAL,dim-1>::general(iit->geometry().type()).position(0,0);

                Dune::FieldVector<REAL,dim> normal = iit->unitOuterNormal( faceLocal );
                //std::cout << "DEBUG: normal = " << normal << std::endl;


                REAL element_volume_s = (*eit).geometry().volume();
#ifdef DIMENSION3
                REAL element_length_s = std::pow( element_volume_s, 1.0/3.0 );
#else
                REAL element_length_s = std::sqrt( element_volume_s );
#endif


                REAL K_inside = tensor_inside[0][0];
                REAL K_outside = tensor_inside[0][0];
                REAL K_effective = tensor_inside[0][0];


                if (iit->neighbor()){

                  //std::cout << "DEBUG: internal face... " << std::endl;
                  // interior face:

                  REAL element_volume_n = iit->outside()->geometry().volume();
#ifdef DIMENSION3
                  REAL element_length_n = std::pow( element_volume_n, 1.0/3.0 );
#else
                  REAL element_length_n = std::sqrt( element_volume_n );
#endif

                  const Dune::FieldVector<REAL,dim>& outsideCellCenterLocal =
                    Dune::ReferenceElements<DF,dim>::general(iit->outside()->type()).position(0,0);

                  Dune::FieldVector<REAL,dim> distance_vector =
                    iit->outside()->geometry().global(outsideCellCenterLocal);

                  // distance of cell centers
                  distance_vector -= insideCellCenterGlobal;
                  REAL distance = distance_vector.two_norm();
                  distance_vector /= distance;

                  Dune::FieldVector<DF,dim> kvector_s;
                  tensor_inside.mv(distance_vector,kvector_s);
                  K_inside = kvector_s.infinity_norm();

                  typename GWP::Traits::DiffusionTensorType tensor_outside(gwp.DiffusionTensor(*(iit->outside()),outsideCellCenterLocal));
                  Dune::FieldVector<DF,dim> kvector_n;
                  tensor_outside.mv(distance_vector,kvector_n);
                  K_outside = kvector_n.infinity_norm();

                  if( element_length_s - element_length_n<1E-12 )
                    K_effective =
                      Dune::Gesis::General::harmonicAverage( K_inside, K_outside );
                  else
                    K_effective =
                      Dune::Gesis::General::harmonicAverageWeightedByDistance( K_inside,
                                                                               K_outside,
                                                                               element_length_s,
                                                                               element_length_n );

                }
                else {
                  //std::cout << "DEBUG: boundary face... " << std::endl;
                  // distance of cell center to boundary
                  Dune::FieldVector<DF,dim> distance_vector = iit->geometry().global(faceLocal);
                  distance_vector -= insideCellCenterGlobal;
                  REAL distance = distance_vector.two_norm();
                  distance_vector /= distance;

                  Dune::FieldVector<DF,dim> kvector_s;
                  tensor_inside.mv(distance_vector,kvector_s);
                  K_inside = kvector_s.infinity_norm();
                  K_effective = K_inside;
                }

                //std::cout << "DEBUG: setup matrix... " << std::endl;

                Dune::FieldVector<REAL,dim> faceElementLocal 
                  = iit->geometryInInside().global(faceLocal);
                //std::cout << "DEBUG: faceElementLocal = " << faceElementLocal << std::endl;

                Dune::FieldVector<REAL,dim> flux;
                darcyflux_dgf.evaluate( *eit, faceElementLocal, flux );
                //std::cout << "DEBUG: flux = " << flux << std::endl;

                // This maybe not such a bad idea:
                //for(int i=0;i<dim;i++)
                //  normal[i] *= normal[i];
                //std::cout << "DEBUG: normalized normal = " << normal << std::endl;


                rhs[faceCounter] = - (flux * normal) / K_effective;

                Dune::FieldVector<REAL,dim> indexvector;
                for(int i=0;i<dim;i++)
                  indexvector[i] = (REAL) i + 1E-3; // This +1E-3 is very important due to the (int) casting in the next line!

                int component_index = (int) std::abs( normal * indexvector );

#ifdef USE_LOCAL_COORD
                //logger << "faceCounter = " << faceCounter << std::endl;
                //logger << "normal = " << normal << std::endl;
                //logger << "indexvector = " << indexvector << std::endl;
                //logger << "component_index = " << component_index << std::endl;

                XY(faceCounter,2*component_index) 
                  = faceElementLocal[component_index] 
                  - insideCellCenterLocal[component_index];
#else
                Dune::FieldVector<REAL,dim> faceGlobal = eit->geometry().global(faceElementLocal);
                //std::cout << "DEBUG: faceGlobal = " << faceGlobal << std::endl;
                XY(faceCounter,2*component_index) 
                  = faceGlobal[component_index];
#endif

                XY(faceCounter,2*component_index+1) = 1.0;

                faceCounter++;

                //logger << std::endl;


              } // intersection loop

            // The very last row of the RHS:
            rhs[2*dim] = pCellcenter;

            // The very last row of the matrix:


#ifdef USE_LOCAL_COORD
            REAL x0 = insideCellCenterLocal[0];
            REAL x1 = insideCellCenterLocal[1];
            REAL x2 = insideCellCenterLocal[2];
#else
            REAL x0 = insideCellCenterGlobal[0];
            REAL x1 = insideCellCenterGlobal[1];
            REAL x2 = insideCellCenterGlobal[2];
#endif

            XY(2*dim,0) = 0.5 * x0 * x0;
            XY(2*dim,1) = x0;
            XY(2*dim,2) = 0.5 * x1 * x1;
            XY(2*dim,3) = x1;
            if(dim==3){
              XY(2*dim,2*dim-2) = 0.5 * x2 * x2;
              XY(2*dim,2*dim-1) = x2;
            }

            XY(2*dim,2*dim) = 1.0;

            //logger << "XY = " << XY << std::endl;
            //logger << "rhs = " << rhs << std::endl;

            XY.gauss_solver( coeff, rhs );

          } // if(eit->partitionType()==Dune::InteriorEntity)

          //std::cout << "  rank = " << gv.comm().rank() 
          //          << " coeff = " << coeff << std::endl;

          this->appendNewElement( *eit, coeff );  // data_container gets filled here!

        } // element loop

        //logger << "RT0_pf: Element loop done." << std::endl;

        // std::cout << "rt0_pressurefield calculated on interior partitions." << std::endl;
        // std::cout << "Start MPI communications..." << std::endl;
        // Do the MPI communication here...

        logger << "RT0_pf: Start data exchange." << std::endl;
        typedef RT0_PressureField_Exchange<GV,ContainerType> DataHandleType;
        DataHandleType datahandle(gv,data_container,blocksize);
        gv.communicate( datahandle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
        logger << "RT0_pf: gv.communicate() done." << std::endl;


        std::stringstream jobtitle;
        jobtitle << "Reconstructing Q2 potential field from RT0 flux on (locally refined) grid.";
        General::log_elapsed_time( watch.elapsed(),
                                   gv.comm(),
                                   General::verbosity,
                                   "RGV",
                                   jobtitle.str() );
      }



      void evaluate( const ElementType& e,
                     const Dune::FieldVector<REAL,dim>& x, // element-local coordinate
                     Dune::FieldVector<REAL,1>& y
                     ) const {

#ifndef USE_ALUGRID
        Idx idx = gv.indexSet().index(e);
#else
        Idx idx = gv.grid().localIdSet().id(e);
#endif

        Vector<REAL> coefficients(blocksize,0.0);
        for(int i=0;i<blocksize;i++){
          //coefficients[i] = data_container[ blocksize * idx + i ];
          typename ContainerType::const_iterator it =
            data_container.find( blocksize*idx + i );
          if(it != data_container.end())
            coefficients[i] = it->second;
          else
            logger << "WARNING: rt0 pressure coefficients not found for element idx = "
                   << idx << std::endl;
        }

#ifdef USE_LOCAL_COORD
        REAL x0 = x[0]-0.5;   // Very important! Must be relative to the center.
        REAL x1 = x[1]-0.5;
        y[0] = coefficients[2*dim];
        y[0] += 0.5*coefficients[0]*x0*x0 + coefficients[1]*x0;
        y[0] += 0.5*coefficients[2]*x1*x1 + coefficients[3]*x1;
        if( dim>2 ){
          REAL x2 = x[dim-1]-0.5;
          y[0] += 0.5*coefficients[2*dim-2]*x2*x2 + coefficients[2*dim-1]*x2;
        }
#else
        Dune::FieldVector<REAL,dim> xglobal = e.geometry().global(x);
        REAL x0 = xglobal[0];
        REAL x1 = xglobal[1];
        y[0] = coefficients[2*dim];
        y[0] += 0.5*coefficients[0]*x0*x0 + coefficients[1]*x0;
        y[0] += 0.5*coefficients[2]*x1*x1 + coefficients[3]*x1;
        if( dim>2 ){
          REAL x2 = xglobal[dim-1];
          y[0] += 0.5*coefficients[2*dim-2]*x2*x2 + coefficients[2*dim-1]*x2;
        }
#endif

      }

      void evaluate_on_root( const ElementType& e,
                             const Dune::FieldVector<REAL,dim>& x, // element-local coordinate
                             Dune::FieldVector<REAL,1>& y
                             ) const {

        // convert to global coordinate wrt to element e
        Dune::FieldVector<REAL,dim> global = e.geometry().global(x);

        if(e.level()>baselevel){
          // get father element
          typedef typename GV::Grid::template Codim<0>::EntityPointer ElementPointer;
          ElementPointer pAncestor = e.father();
          while( pAncestor->level() > baselevel )
            pAncestor = (*pAncestor).father();
          // convert to local coordinate wrt to element *pAncestor
          Dune::FieldVector<REAL,dim> xx = (*pAncestor).geometry().local( global );
          this->evaluate( *pAncestor, xx, y );
        }
        else{
          this->evaluate(e, x, y);
        }
        return;

      }


      template<typename GV_TP>
      void plot2vtu( const GV_TP& gv_tp, const std::string filename ) const {

        // Get the conductivity field on the grid for the vtkwriter
        typedef typename GV_TP::Grid GRIDTYPE; // get the grid-type out of the gridview-type
        // typedef typename GV_TP::Grid::template Codim < 0 > ::Entity Entity;
        // typedef typename GV_TP::Grid::template Codim < 0 > ::EntityPointer EntityPointer;

        //typedef typename GV::Traits::template Codim < 0 > ::Iterator ElementLeafIterator;
        typedef typename GV_TP::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementLeafIterator;
        // typedef typename ElementLeafIterator::Entity::Geometry LeafGeometry;

        // make a mapper for codim 0 entities in the leaf grid
        Dune::LeafMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout>
          mapper(gv_tp.grid()); // get the underlying hierarchical grid out ot the gridview

        // make a mapper for codim 0 entities in the level(0) grid
        //Dune::LevelMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout>
        //  mapper(gv.grid(),baselevel); // get the underlying hierarchical grid out ot the gridview
        const int nGridCells = mapper.size();

        std::vector<REAL> hydraulic_head(nGridCells);

        if(gv_tp.comm().rank()==0)
          std::cout << "RT0_pf: Plot reconstructed hydraulic head: Loop over leaf elements..." << std::endl;

        const typename GV_TP::IndexSet& indexset = gv_tp.indexSet();

        for (ElementLeafIterator it = gv_tp.template begin<0,Dune::All_Partition> ()
               ; it != gv_tp.template end<0,Dune::All_Partition> (); ++it) {
          int index = indexset.index(*it);

          Dune::FieldVector<REAL,dim> xglobal = it->geometry().center();

          Dune::FieldVector<REAL,1> head(0);

          Dune::FieldVector<REAL,dim> xlocal(0.5);// = it->geometry().local(xglobal);
          this->evaluate_on_root( *it, xlocal, head );


          hydraulic_head[index] = head[0];
        }

        Dune::Gesis::VTKPlot::output_vector_to_vtu( gv_tp,
                                                    hydraulic_head,
                                                    filename,
                                                    "head on leafView",
                                                    0,
                                                    true,
                                                    0
                                                    );


      } // end of plot2vtu()



    }; // class RT0_PressureField

  }
}

#endif // DUNE_GESIS_RT0_PRESSUREFIELD_HH






/*

  Herleitung: Rekonstruktion des Potentialfeldes phi aus dem Darcyfluss q=(qx,qy,qz)

  Ansatz für das Potential:
  phi(x,y,z) = 0.5*a(x-x0)² + bx + 0.5*c(y-y0)² + dy + 0.5*e(z-z0)² + fz + g

  1 Gleichung, 1 Unbekante g an der Stelle center = (x0,y0,z0):
  h(center) = phi(center)

  3 Gleichungen mit 6 Unbekannten a,...,f:
  - K * grad phi (x,y,z) = -K * ( ax+b, cy+d, ez+f ) = (qx, qy, qz)

  An zwei Stellen face1 und face2 evaluiert, ergibt 6 Gleichungen:
  face1   (east): x1,y1,z1 with normal (1,0,0)
  face2  (north): x2,y2,z2 with normal (0,1,0)
  face3    (top): x3,y3,z3 with normal (0,0,1)
  face4   (west): x4,y4,z4 with normal (-1,0,0)
  face5  (south): x5,y5,z5 with normal (0,-1,0)
  face6 (bottom): x6,y6,z6 with normal (0,0,-1)

  1.) -K1 (a*x1 + b) = qx (face1)
  =>
  (a*x1 + b) = - qx(face1) / K1


  2.) -K2 (c*y1 + d) = qy (face2)
  =>
  (c*y1 + d) = - qy(face2) / K2


  3.) -K3 (e*z1 + f) = qz (face3)
  =>
  (e*z1 + f) = - qz(face3) / K3


  4.) -K4 (a*x2 + b) = qx (face4)
  =>
  (a*x2 + b) = - qx(face4) / K4


  5.) -K5 (c*y2 + d) = qy (face5)
  =>
  (c*y2 + d) = - qy(face5) / K5


  6.) -K6 (e*z2 + f) = qz (face6)
  =>
  (e*z2 + f) = - qz(face6) / K6


  7.) phi(x0,y0,z0) 
  = 0.5*a*x0² + b*x0 
  + 0.5*c*y0² + d*y0 
  + 0.5*e*z0² + f*z0 + g 
  = h(x0,y0,z0)



  ==> Matrix-Schreibweise: XX * coeff = rhs

  a        b       c   d       e   f   g

  x1       1       0   0       0   0   0       a         - qx(face1) / K1
  0        0      y1   1       0   0   0       b         - qy(face2) / K2
  0        0       0   0      z1   1   0       c         - qz(face3) / K3
  x2       1       0   0       0   0   0   *   d    =    - qx(face4) / K4
  0        0      y2   1       0   0   0       e         - qy(face5) / K5
  0        0       0   0      z2   1   0       f         - qz(face6) / K6
  0.5*x0² x0  0.5*y0² y0  0.5*z0² z0   1       g             h(center)


  Evaluierung:
  phi(x,y,z) = 0.5*a*(x-x0)² + b*(x-x0) + 0.5*c*(y-y0)² + d*(y-y0) + 0.5*e*(z-z0)² + f*(z-z0) + g






*/
