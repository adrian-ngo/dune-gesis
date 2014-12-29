#ifndef DUNE_GESIS_DGF_PressureField_HH
#define DUNE_GESIS_DGF_PressureField_HH

/*
  Note: This class is an extension of DiscreteGridFunction used to evaluate the CCFV pressure.
  It has the method evaluate_on_root() used to evaluate the pressure on the baselevel element.
  Since the head solution of CCFV is P0, of course it will return only one value h=h(center)
  no matter how many sub-cells the baselevel reference element has.
  But for global refinenment on YASP this seems to offer
  the best element order for the solution of the transport equation.
*/


template<typename GFS,typename VCType>
class DGF_PressureField
  : public Dune::PDELab::DiscreteGridFunction<GFS,VCType> {

private:
  typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> base_type;

  typedef typename GFS::Traits::GridViewType GV;
  enum{ dim = GV::dimension };
  typedef typename GV::Traits::template Codim<0>::Entity ElementType;

  const int baselevel;

public:
  DGF_PressureField( const GFS& gfs, const VCType& vc, const int baselevel_=0 )
    : base_type( gfs, vc ),
      baselevel(baselevel_)
  {}

  void evaluate_on_root( const ElementType& e,
                         const Dune::FieldVector<REAL,dim>& x, // element-local coordinate
                         Dune::FieldVector<REAL,1>& y // result
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
    // typedef typename GV_TP::Grid GRIDTYPE; // get the grid-type out of the gridview-type
    // typedef typename GV_TP::Grid::template Codim < 0 > ::Entity Entity;
    // typedef typename GV_TP::Grid::template Codim < 0 > ::EntityPointer EntityPointer;

    //typedef typename GV::Traits::template Codim < 0 > ::Iterator ElementLeafIterator;
    typedef typename GV_TP::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementLeafIterator;
    //typedef typename ElementLeafIterator::Entity::Geometry LeafGeometry;

    // make a mapper for codim 0 entities in the leaf grid
    //Dune::LeafMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout>
    //mapper(gv_tp.grid()); // get the underlying hierarchical grid out ot the gridview

    // make a mapper for codim 0 entities in the level(0) grid
    //Dune::LevelMultipleCodimMultipleGeomTypeMapper<GRIDTYPE, P0Layout>
    //  mapper(gv.grid(),baselevel); // get the underlying hierarchical grid out ot the gridview
    //const int nGridCells = mapper.size();

    const int nGridCells = gv_tp.size(0);

#ifdef DEBUG_LOG
    logger << "DEBUG: leafView mapper.size() = " << nGridCells << std::endl;
#endif

    std::vector<REAL> hydraulic_head(nGridCells);

    if(gv_tp.comm().rank()==0)
      std::cout << "Plot hydraulic head: Loop over leaf elements..." << std::endl;

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
                                                1
                                                );

  } // end of plot2vtu()

};

#endif // DUNE_GESIS_DGF_PressureField_HH
