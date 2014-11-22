// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GESIS_ADAPTIVITY_HH
#define DUNE_GESIS_ADAPTIVITY_HH

#include<dune/common/exceptions.hh>

#include<limits>
#include<vector>
#include<map>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>

// for CTLFA
#include<dune/pdelab/common/function.hh>
// for InterpolateBackendStandard
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
// for intersectionoperator
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "RedistributeDataHandle.hh"

namespace Dune {
  namespace Gesis {

    /*! @class CoeffsToLocalFunctionAdapter
     *
     *  @brief A wrapper class interpreting coefficients as a LocalFunction
     *
     *         A wrapper class interpreting a vector of coefficients (and a local finite element)
     *         as a LocalFunction living on an Elem. The function may be evaluated on another Elem,
     *         e.g. a child or grandchild. Useful if the GridFunctionSpace the coefficients originate from is
     *         no longer available (think adaptation).
     *
     *  @tparam CoeffType Type of the interpreted coefficients
     *  @tparam DGF       The DiscreteGridFunction this should mimic
     *  @tparam FEM       The FiniteElementMap used
     *  @tparam E         Type for Elems
     */
    template<class CoeffType, class DGF, class FEM, class E>
    class CoeffsToLocalFunctionAdapter : public Dune::PDELab::FunctionInterface<typename DGF::Traits,
                                                                  CoeffsToLocalFunctionAdapter<CoeffType,DGF,FEM,E> >
    {

    public:

      typedef typename DGF::Traits Traits;
      typedef typename FEM::Traits::FiniteElementType FiniteElement;

      /** @brief The constructor.
       *
       *  @param[in] coeff_ A vector of coefficients
       *  @param[in] fem_ A FiniteElementMap used for interpretation of the coefficients
       *  @param[in] from_ The Elem on which the LocalFunction lives
       *  @param[in] to_ The Elem on which the LocalFunction is evaluated (should be e.g. a child
       *             of from_ to make the values of the shape functions meaningful)
       */
      CoeffsToLocalFunctionAdapter (const std::vector<CoeffType>& coeff_, const FEM& fem_, const E& from_, const E& to_)
        : coeff(coeff_), fem(fem_), from(from_), to(to_), fe(fem.find(from)) {}

      /** @brief Special version of the constructor with to = from.
       *
       *  @param[in] coeff_ A vector of coefficients
       *  @param[in] fem_ A FiniteElementMap used for interpretation of the coefficients
       *  @param[in] elem_ The Elem on which the LocalFunction lives
       */
      CoeffsToLocalFunctionAdapter (const std::vector<CoeffType>& coeff_, const FEM& fem_, const E& elem_)
        : coeff(coeff_), fem(fem_), from(elem_), to(elem_), fe(fem.find(from)) {}

      /** @brief Evaluate the LocalFunction at the given position
       *
       *  @param[in]  x The position in local coordinates
       *  @param[out] y The result of the evaluation
       */
      inline void evaluate (const typename Traits::DomainType& x, typename Traits::RangeType& y) const
      {
        std::vector<typename Traits::RangeType> yVector;
        fe.localBasis().evaluateFunction(from.geometry().local(to.geometry().global(x)),yVector);
        if (yVector.size() != coeff.size()) DUNE_THROW(Dune::Exception,"Coefficient vector has wrong length in CoeffsToLocalFunctionAdapter");
        typename Traits::RangeType sum=0.0;
        for (unsigned int i = 0; i < yVector.size(); ++i)
          sum += yVector[i]*coeff[i];
        y = sum;
      }

    private:
      const std::vector<CoeffType>& coeff;
      const FEM& fem;
      const E& from;
      const E& to;
      const FiniteElement& fe;
    };


    /*! @class L2Projection
     *
     * @brief @todo
     *
     * @tparam GFSU Type of ansatz space
     * @tparam U    Container class for the solution
     */
    template<class GFSU, class U>
    class L2Projection
    {
      typedef typename GFSU::Traits::GridViewType::Grid Grid;
      typedef typename Grid::template Codim<0>::Entity Element;
      typedef Dune::PDELab::LocalFunctionSpace<GFSU> LFSU;
      typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;

    public:

      /*! @brief The constructor.
       *
       * @todo Doc params!
       */
      explicit L2Projection(int intorder_=2) : intorder(intorder_), haveMatrix(), matrix() {}

      /*! @brief Calculate the L2 scalar product of functions X and Y on e, but in the geometry of its ancestor
       *
       * @todo Doc template params
       * @todo params
       */
      template<typename F1, typename F2>
      void apply (const Element& father, const Element& e, F1& X, F2& Y, typename U::ElementType& u) const
      {
        Dune::GeometryType gt = e.geometry().type();

        const int dim  = Element::Geometry::dimension;
        const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);
        typename F1::Traits::RangeType x;
        typename F2::Traits::RangeType y;

        // iterate over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            X.evaluate(it->position(),x);
            Y.evaluate(it->position(),y);
            const typename U::ElementType factor = it->weight()
              * e.geometry().integrationElement(it->position())
              / father.geometry().integrationElement(father.geometry().local(e.geometry().global(it->position())));

            u += x * y * factor;
          }
      }

      /*! @brief Calculate the L2 norm of X - Y on e, but in the geometry of its ancestor
       *
       * @todo Doc template params
       * @todo Doc params
       */
      template<typename F1, typename F2>
      void error (const Element& father, const Element& e, F1& X, F2& Y, typename U::ElementType& u) const
      {
        Dune::GeometryType gt = e.geometry().type();

        const int dim  = Element::Geometry::dimension;
        const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);
        typename F1::Traits::RangeType x;
        typename F2::Traits::RangeType y;

        typename U::ElementType temp = 0.;

        // iterate over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            X.evaluate(it->position(),x);
            Y.evaluate(it->position(),y);
            const typename U::ElementType factor = it->weight();
            //  * e.geometry().integrationElement(it->position());
            //  / father.geometry().integrationElement(father.geometry().local(e.geometry().global(it->position())));

            temp += (x-y) * (x-y) * factor;
          }

        u += sqrt(temp);
      }

      /*! @brief Calculate the inverse local mass matrix, used in the local L2 projection
       *
       * @todo Doc template params
       * @todo Doc params
       */
      template <class CTLFA, class FEM>
      const std::vector<double>& inverseMassMatrix(const Element& e, const FEM& fem, int k)
      {
        const Dune::GeometryType gt = e.geometry().type();
        if (haveMatrix.find(gt) != haveMatrix.end())
          {
            return matrix[gt][k];
          }

        const int localSize = (fem.find(e)).localBasis().size();
        typename CTLFA::Traits::RangeType x;
        typename CTLFA::Traits::RangeType y;

        std::vector<std::vector<typename U::ElementType> > massMatrix(localSize,std::vector<typename U::ElementType>(localSize,0.));
        std::vector<std::vector<typename U::ElementType> > inverseMatrix(localSize,std::vector<typename U::ElementType>(localSize,0.));

        std::vector<double> phi_i(localSize,0.), phi_j(localSize,0.);
        CTLFA ctlfa_i(phi_i,fem,e,e);
        CTLFA ctlfa_j(phi_j,fem,e,e);

        for (int i = 0; i < localSize; ++i)
          {
            ++phi_i[i];
            for (int j = 0; j < localSize; ++j)
              {
                ++phi_j[j];
                (*this).template apply<CTLFA,CTLFA>(e,e,ctlfa_i,ctlfa_j,massMatrix[i][j]);
                --phi_j[j];
              }
            --phi_i[i];
          }

        for (int i = 0; i < localSize; ++i)
          {
            inverseMatrix[i][i] = 1.;
          }

        for (int i = 0; i < localSize; ++i)
          {
            for (int j = 0; j < localSize; ++j)
              {
                if (i != j)
                  {
                    const typename U::ElementType factor = massMatrix[j][i]/massMatrix[i][i];
                    for (int l = 0; l < localSize; ++l)
                      {
                        massMatrix[j][l] -= factor * massMatrix[i][l];
                        inverseMatrix[j][l] -= factor * inverseMatrix[i][l];
                      }
                  }
              }
          }
        for (int i = localSize - 1; i >= 0; --i)
          {
            for (int j = localSize - 1; j >= 0; --j)
              {
                if (i != j)
                  {
                    const typename U::ElementType factor = massMatrix[j][i]/massMatrix[i][i];
                    for (int l = localSize - 1; l >= 0; --l)
                      {
                        massMatrix[j][l] -= factor * massMatrix[i][l];
                        inverseMatrix[j][l] -= factor * inverseMatrix[i][l];
                      }
                  }
              }
          }
        for (int i = 0; i < localSize; ++i)
          {
            const typename U::ElementType factor = massMatrix[i][i];
            for (int j = 0; j < localSize; ++j)
              {
                massMatrix[i][j] /= factor;
                inverseMatrix[i][j] /= factor;
              }
          }

        matrix.insert(std::pair<Dune::GeometryType,std::vector<std::vector<typename U::ElementType> > >(gt,inverseMatrix));
        haveMatrix.insert(gt);
        return matrix[gt][k];
      }

    private:

      const int intorder;

      // only for inverseMassMatrix
      std::set<Dune::GeometryType> haveMatrix;
      std::map<Dune::GeometryType, std::vector<std::vector<typename U::ElementType> > > matrix;
    };




    /*! @class GridAdaptor
     *
     * @brief Class for automatic adaptation of the grid.
     *
     *        The GridAdaptor capsules the act of deciding which Elems to refine and coarsen,
     *        adapting the grid, and transfering the solution from the old grid to the new one.
     *        Currrently this only works for scalar solutions.
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFSU       Type of ansatz space, we need to update it after adaptation
     * @tparam U          Container class of the solution
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid, class GFSU, class U, class Projection>
    class GridAdaptor
    {
      typedef typename Grid::LeafGridView LeafGridView;
      typedef typename LeafGridView::template Codim<0>
      ::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
      typedef typename Grid::template Codim<0>::Entity Element;
      typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
      typedef typename Grid::GlobalIdSet IdSet;
      typedef typename IdSet::IdType IdType;
      typedef typename Grid::LeafIndexSet IndexSet;
      typedef typename IndexSet::IndexType IndexType;
      typedef Dune::PDELab::LocalFunctionSpace<GFSU> LFSU;
      typedef Dune::PDELab::DiscreteGridFunction<GFSU, U> DGF;
      typedef typename GFSU::Traits::FiniteElementMapType FEM;
      typedef Dune::PDELab::InterpolateBackendStandard IB;
      typedef CoeffsToLocalFunctionAdapter<typename U::ElementType,DGF,FEM,Element> CTLFA;
      typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;

    public:
      typedef std::map<IdType,std::vector<typename U::ElementType> > MapType;


      /*! @brief The constructor.
       *
       * @param grid_       The grid we want to adapt
       * @param gfsu_       The ansatz space, we need to update it
       * @param projection_ The Projection used when Elems vanish
       */
      GridAdaptor() {}

      /* @brief @todo
       *
       * @param[in]  u           The solution that will be saved
       * @param[out] transferMap The map containing the solution during adaptation
       */
      void backupData(Grid& grid, GFSU& gfsu, Projection& projection, U& u, MapType& transferMap)
      {
        const IdSet& idset = grid.globalIdSet();
        LFSU lfsu(gfsu);
        DGF dgf(gfsu,u);
        const FEM& fem = gfsu.finiteElementMap();
        // IB ib = IB();
        std::vector<typename U::ElementType> ul;

        // iterate over all elems
        LeafGridView leafView = grid.leafView();
        for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
             it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
          {
            const Element& e = *it;

            // save local coeffs in map
            lfsu.bind(e);
            //lfsu.vread(u,transferMap[e.level()][idset.id(e)]);
            lfsu.vread(u,transferMap[idset.id(e)]);

            // save local coeffs of father in map
            if (e.mightVanish())
              {
                ElementPointer father = e.father();

                projection.template inverseMassMatrix<CTLFA,FEM>(e,fem,0);
                {
                  const int localSize = (fem.find(*father)).localBasis().size();
                  std::vector<typename U::ElementType> uCoarse(localSize,0.);
                  std::vector<typename U::ElementType> coarseBasis(localSize,0.);
                  const typename Element::HierarchicIterator& hbegin = (*father).hbegin(grid.maxLevel());
                  const typename Element::HierarchicIterator& hend   = (*father).hend(grid.maxLevel());

                  for (typename Element::HierarchicIterator hit = hbegin; hit != hend; ++hit)
                    {
                      // only evaluate on entities with data
                      if ((*hit).isLeaf())
                        {
                          typedef Dune::PDELab::GridFunctionToLocalFunctionAdapter<DGF> GFTLFA;
                          GFTLFA gftlfa(dgf,*hit);
                          CTLFA  ctlfa(coarseBasis,fem,*father,*hit);

                          // iterate over canonical basis vectors
                          for (unsigned int i = 0; i < coarseBasis.size(); ++i)
                            {
                              coarseBasis = projection.template inverseMassMatrix<CTLFA,FEM>(e,fem,i);
                              projection.template apply<GFTLFA,CTLFA>(*father,*hit,gftlfa,ctlfa,uCoarse[i]);
                            }
                        }
                    }
                  //transferMap[(*father).level()][idset.id(*father)] = uCoarse;
                  transferMap[idset.id(*father)] = uCoarse;
                }

                // conforming grids may need this
                while ((*father).mightVanish())
                  {
                    father = (*father).father();
                    {
                      const int localSize = (fem.find(*father)).localBasis().size();
                      std::vector<typename U::ElementType> uCoarse(localSize,0.);
                      std::vector<typename U::ElementType> coarseBasis(localSize,0.);
                      const typename Element::HierarchicIterator& hbegin = (*father).hbegin(grid.maxLevel());
                      const typename Element::HierarchicIterator& hend   = (*father).hend(grid.maxLevel());

                      for (typename Element::HierarchicIterator hit = hbegin; hit != hend; ++hit)
                        {
                          // only evaluate on entities with data
                          if ((*hit).isLeaf())
                            {
                              typedef Dune::PDELab::GridFunctionToLocalFunctionAdapter<DGF> GFTLFA;
                              GFTLFA gftlfa(dgf,*hit);
                              CTLFA  ctlfa(coarseBasis,fem,*father,*hit);

                              // iterate over canonical basis vectors
                              for (unsigned int i = 0; i < coarseBasis.size(); ++i)
                                {
                                  coarseBasis = projection.template inverseMassMatrix<CTLFA,FEM>(e,fem,i);
                                  projection.template apply<GFTLFA,CTLFA>(*father,*hit,gftlfa,ctlfa,uCoarse[i]);
                                }
                            }
                        }
                      //transferMap[(*father).level()][idset.id(*father)] = uCoarse;
                      transferMap[idset.id(*father)] = uCoarse;
                    }
                  }
              }
          }
      }

      /* @brief @todo
       *
       * @param[out] u           The solution after adaptation
       * @param[in]  transferMap The map that contains the information for the rebuild of u
       */
      void replayData(Grid& grid, GFSU& gfsu, Projection& projection, U& u, MapType& transferMap)
      {
        const IdSet& idset = grid.globalIdSet();
        LFSU lfsu(gfsu);
        const FEM& fem = gfsu.finiteElementMap();
        IB ib = IB();
        std::vector<typename U::ElementType> ul;
        std::vector<typename U::ElementType> ulc;
        U uc(gfsu,0.0);
        const IndexSet& indexset = grid.leafIndexSet();
        std::vector<typename U::ElementType> ug(indexset.size(IndexSet::dimension),0.);
        std::vector<typename U::ElementType> ugc(indexset.size(IndexSet::dimension),0.);
        Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGElementLayout> mapper(grid);

        // iterate over all elems
        LeafGridView leafView = grid.leafView();
        for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
             it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
          {
            const Element& e = *it;
            lfsu.bind(e);
            const IdType& id = idset.id(e);
            const int level = e.level();

            if (e.isNew()) // id is not in map, we have to interpolate
              {
                // find ancestor with data
                int levelAncestor = level - 1;
                ElementPointer pAncestor = e.father();
                //while (transferMap[levelAncestor].find(idset.id(*pAncestor)) == transferMap[levelAncestor].end())
                while (transferMap.find(idset.id(*pAncestor)) == transferMap.end())
                  {
                    //this ancestor does not have data, check next one
                    if (levelAncestor == 0)
                      DUNE_THROW(Dune::Exception,
                                 "transferMap of GridAdaptor didn't contain ancestor of element with id " << id);
                    pAncestor = (*pAncestor).father();
                    --levelAncestor;
                  }

                // this one has the data
                const Element& Ancestor = *pAncestor;
                //CTLFA ctlfa(transferMap[levelAncestor][idset.id(Ancestor)],fem,Ancestor,e);
                CTLFA ctlfa(transferMap[idset.id(Ancestor)],fem,Ancestor,e);

                ul.clear();
                ib.interpolate(fem.find(e),ctlfa,ul);

                lfsu.vadd(ul,u);
              }
            else // this entity is not new and should have data
              {
                //lfsu.vadd(transferMap[level][id],u);
                lfsu.vadd(transferMap[id],u);
              }

            ulc = std::vector<typename U::ElementType>(lfsu.size(),1.0);
            lfsu.vadd(ulc,uc);
          }

        typedef Dune::PDELab::AddDataHandle<GFSU,U> Handle;
        Handle addHandle1(gfsu,u);
        leafView.communicate (addHandle1,
                              Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
        Handle addHandle2(gfsu,uc);
        leafView.communicate (addHandle2,
                              Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

        for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
             it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
          {
            const Element& e = *it;
            lfsu.bind(e);

            ul = std::vector<typename U::ElementType>(lfsu.size(),0.0);
            ulc = std::vector<typename U::ElementType>(lfsu.size(),0.0);
            lfsu.vread(u,ul);
            lfsu.vread(uc,ulc);

            for (unsigned int i = 0; i<ul.size();++i)
              {
                if(ulc[i]>1)
                  {
                    ul[i] /= ulc[i];
                    ulc[i] = 1.;
                  }
              }
            lfsu.vwrite(ul,u);
            lfsu.vwrite(ulc,uc);
          }
      }

    };


    /*! grid adaptation as a function
     *
     * @brief adapt a grid, corresponding function space and solution vectors
     *
     * Assumes that the grid's elements have been marked for refinement and coarsening appropriately before
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFS        Type of ansatz space, we need to update it after adaptation
     * @tparam X          Container class for DOF vectors
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid, class GFS, class X, class Projection=L2Projection<GFS,X> >
    void adapt_grid (Grid& grid, GFS& gfs, X& x1, Projection projection=L2Projection<GFS,X>() )
    {
      GridAdaptor<Grid,GFS,X,Projection> grid_adaptor;

      // prepare the grid for refinement
      grid.preAdapt();
      
      // save u
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap1;
      grid_adaptor.backupData(grid,gfs,projection,x1,transferMap1);
      
      // adapt the grid
      grid.adapt();
      
      // update the function spaces
      gfs.update();
      
      // reset u
      x1 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x1,transferMap1);
      
      // clean up
      grid.postAdapt();
    }

    /*! grid adaptation as a function
     *
     * @brief adapt a grid, corresponding function space and solution vectors
     *
     * Assumes that the grid's elements have been marked for refinement and coarsening appropriately before
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFS        Type of ansatz space, we need to update it after adaptation
     * @tparam X          Container class for DOF vectors
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid, class GFS, class X, class Projection=L2Projection<GFS,X> >
    void adapt_grid (Grid& grid, GFS& gfs, X& x1, X& x2, Projection projection=L2Projection<GFS,X>() )
    {
      GridAdaptor<Grid,GFS,X,Projection> grid_adaptor;

      // prepare the grid for refinement
      grid.preAdapt();
      
      // save solution
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap1;
      grid_adaptor.backupData(grid,gfs,projection,x1,transferMap1);
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap2;
      grid_adaptor.backupData(grid,gfs,projection,x2,transferMap2);
      
      // adapt the grid
      grid.adapt();
      
      // update the function spaces
      gfs.update();
      
      // interpolate solution
      x1 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x1,transferMap1);
      x2 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x2,transferMap2);
      
      // clean up
      grid.postAdapt();
    }





    /*! grid adaptation as a function
     *
     * @brief adapt a grid, corresponding function space and solution vectors
     *
     * Assumes that the grid's elements have been marked for refinement and coarsening appropriately before
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GV         Gridview with some ordering system, must offer update() method
     * @tparam GFS        Type of ansatz space, we need to update it after adaptation
     * @tparam X          Container class for DOF vectors
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid,
             class GFS_GW,
             class VCType_GW,
             class GV,
             class GFS,
             class YFG,
             class IDT,
             class DIR
#ifdef BACKUP_REPLAY
             , class X
             , class Projection=L2Projection<GFS,X> 
#endif
             >
    void adapt_grid_rgv( Grid& grid,
                         GFS_GW& gfs_gw,
                         VCType_GW& vc_h,
                         GV& gv,
                         GFS& gfs,
                         YFG& yfg,
                         const IDT& inputdata,
                         const DIR& dir,
                         const UINT step
#ifdef BACKUP_REPLAY
                         , X& x1
                         , Projection projection=L2Projection<GFS,X>() 
#endif
                         )
    {
      /*
      std::cout << "before adapt: grid.maxLevel() = "
                << grid.maxLevel() << std::endl;

      std::cout << "before adapt: grid.levelView( grid.maxLevel() ).size(0) = "
                << grid.levelView( grid.maxLevel() ).size(0) << std::endl;

      std::cout << "before adapt: grid.leafView().size(0) = "
                << grid.leafView().size(0) << std::endl;
      */

      // prepare the grid for refinement
      //grid.preAdapt();

#ifdef BACKUP_REPLAY
      //GridAdaptor<Grid,GFS,X,Projection> grid_adaptor;
      // save u
      //typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap1;
      //grid_adaptor.backupData(grid,gfs,projection,x1,transferMap1);
#endif

      // adapt the grid
      std::cout << "grid.adapt()..." << std::endl;
      grid.adapt();
      std::cout << "grid.adapt() done." << std::endl;

      
      //#ifdef USE_ALUGRID

      if( gv.comm().size() > 1 ) {

        typedef typename Grid::LocalIdSet::IdType Id;
        typedef std::map<Id,REAL> TagMap;
        TagMap elementTagMap;
        const typename Grid::LocalIdSet &lis = grid.localIdSet();

        typedef typename GFS_GW::Traits::GridViewType GV_GW;
        const GV_GW &gv_gw = gfs_gw.gridView();
        //typedef typename Grid::LevelGridView GV_GW;
        //const int baselevel = 0;
        //const GV_GW &gv_gw = grid.levelView( baselevel );

        for( typename GV_GW::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator eit = gv_gw.template begin<0,Dune::Interior_Partition>()
               ; eit != gv_gw.template end<0,Dune::Interior_Partition>()
               ; ++eit ) {
          elementTagMap[ lis.id( *eit ) ] =
            vc_h.base()[ gv_gw.indexSet().index( *eit ) ];
        }

        // loadbalance the grid
        typedef RedistributeDataHandle<Grid,VCType_GW> DataHandle;
        //typedef CoarseGridP0Datahandle<Grid,VCType_GW> DataHandle;

        DataHandle dh(grid,elementTagMap /* vc_h */);
        grid.loadBalance
          (static_cast<Dune::CommDataHandleIF<DataHandle, int>&>(dh));

        gv_gw.communicate( dh, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);

        // yfg must be re-loaded for the new processor distribution!

        Vector<REAL> local_Yfield_vector;
        Vector<UINT> local_count;
        Vector<UINT> local_offset;

        Dune::Gesis
          ::HDF5Tools::read_parallel_from_HDF5( gv_gw
                                                , inputdata
                                                , local_Yfield_vector
                                                , "/YField"
                                                , local_count
                                                , local_offset
                                                , yfg.getFilename()
                                                , 1
                                                , FEMType::DG // P0
                                                , 0 //baselevel
                                                );
        yfg.parallel_import_from_local_vector( local_Yfield_vector, local_count, local_offset );

        logger << "adaptivity: loglistOfAllWellCenters: " << std::endl;
        inputdata.loglistOfAllWellCenters();
        
        logger << "adaptivity: inputdata all well data: " << std::endl;
        for( int i=0; i<inputdata.setups[0].wdlist.pointdata_vector.size(); ++i ){
          logger << " x=" << inputdata.setups[0].wdlist.pointdata_vector[i].x
                 << " y=" << inputdata.setups[0].wdlist.pointdata_vector[i].y
                 << " z=" << inputdata.setups[0].wdlist.pointdata_vector[i].y
                 << " w=" << inputdata.setups[0].wdlist.pointdata_vector[i].well_conductivity
                 << " c=" << inputdata.setups[0].wdlist.pointdata_vector[i].concentration
                 << std::endl;
        }

        yfg.setWellConductivities( gv_gw );


#ifdef DEBUG_PLOT
        std::stringstream vtu_yfg_step;
        vtu_yfg_step << dir.vtudir << "/test_y_step" << step;
        yfg.plot2vtu( gv_gw, vtu_yfg_step.str(), "Y_orig", 0 );
#endif // DEBUG_PLOT

        gfs_gw.update();  // gfs_gw gets new size
        vc_h = VCType_GW(gfs_gw);  // resizing done implicitly

        for( typename GV_GW::template Codim<0>::template Partition<Dune::All_Partition>::Iterator eit = gv_gw.template begin<0,Dune::All_Partition>()
               ; eit != gv_gw.template end<0,Dune::All_Partition>()
               ; ++eit )
          vc_h.base()[gv_gw.indexSet().index(*eit)] = elementTagMap[lis.id(*eit)];
      
        elementTagMap.clear();

#ifdef DEBUG_PLOT
        std::stringstream vtu_h_step;
        vtu_h_step << dir.vtudir << "/test_h_step" << step;
        Dune::Gesis::VTKPlot::output2vtu( gfs_gw,
                                          vc_h,
                                          vtu_h_step.str(),
                                          "h_orig",
                                          inputdata.verbosity,
                                          true,
                                          0 );
#endif // DEBUG_PLOT
      }
      //#endif // USE_ALUGRID



      /*
      std::cout << "after adapt: grid.maxLevel() = "
                << grid.maxLevel() << std::endl;

      std::cout << "after adapt: grid.levelView( grid.maxLevel() ).size(0) = "
                << grid.levelView( grid.maxLevel() ).size(0) << std::endl;

      std::cout << "after adapt: grid.leafView().size(0) = "
                << grid.leafView().size(0) << std::endl;
      */

      // reorder updated gridview
      gv.update( gfs_gw,
                 vc_h,
                 grid.leafGridView()
                 //grid.levelView( grid.maxLevel() )
                 );

      std::cout << "gv.size(0) = "
                << gv.size(0) << std::endl;

      // update the function spaces
      gfs.update();
      
#ifdef BACKUP_REPLAY
      // reset u
      //x1 = X(gfs,0.0);
      //grid_adaptor.replayData(grid,gfs,projection,x1,transferMap1);
#endif
      // clean up
      //grid.postAdapt();
    }


    template<typename T,typename GV>
    void error_fraction( const T& x,
                         const GV& gv,
                         typename T::ElementType alpha,
                         typename T::ElementType beta,
                         typename T::ElementType& eta_alpha,
                         typename T::ElementType& eta_beta,
                         int verbose=0 )
    {
      if (verbose>=VERBOSITY_DEBUG_LEVEL) 
        std::cout << "+++ error fraction: alpha=" << alpha
                  << " beta=" << beta << std::endl;
      const int steps=20; // max number of bisection steps
      typedef typename T::ElementType NumberType;

      NumberType total_error = x.base().two_norm2();
      //NumberType total_error = x.one_norm();

      NumberType max_error = x.infinity_norm();

      if (verbose >= VERBOSITY_DEBUG_LEVEL){
        logger << "error_fraction: total_error squared = " 
               << total_error
               << std::endl
               << "error_fraction: max_error = " 
               << max_error
               << std::endl;
      }

      total_error = gv.comm().sum( total_error );
      max_error = gv.comm().max( max_error );

      if (verbose >= VERBOSITY_DEBUG_LEVEL){
        logger << "error_fraction: comm: total_error squared = " 
               << total_error
               << std::endl
               << "error_fraction: comm: max_error = " 
               << max_error
               << std::endl;
        std::cout << "error_fraction: comm: total_error squared = " 
                  << total_error
                  << std::endl
                  << "error_fraction: comm: max_error = " 
                  << max_error
                  << std::endl;
      }

      NumberType eta_alpha_left = 0.0;
      NumberType eta_alpha_right = max_error;
      NumberType eta_beta_left = 0.0;
      NumberType eta_beta_right = max_error;
      for (int j=1; j<=steps; j++)
        {
          logger << "bisection step j = " << j << std::endl;

          eta_alpha = 0.5*(eta_alpha_left+eta_alpha_right);
          eta_beta = 0.5*(eta_beta_left+eta_beta_right);
          NumberType sum_alpha=0.0;
          NumberType sum_beta=0.0;
          unsigned int alpha_count = 0;
          unsigned int beta_count = 0;
          for (unsigned int i=0; i<x.N(); i++) {
            REAL val = x.base()[i];

            if (val>=eta_alpha) { sum_alpha += val*val; alpha_count++;}
            if (val< eta_beta) { sum_beta += val*val; beta_count++;}
            //if (val>=eta_alpha) { sum_alpha += val; alpha_count++;}
            //if (val< eta_beta) { sum_beta += val; beta_count++;}
          }
          if (verbose >= VERBOSITY_DEBUG_LEVEL)
            {
              logger << j << " eta_alpha=" << eta_alpha << " alpha_fraction=" << sum_alpha/total_error 
                        << " elements: " << alpha_count << " of " << x.N() << std::endl;
              logger << j << " eta_beta=" << eta_beta << " beta_fraction=" << sum_beta/total_error 
                        << " elements: " << beta_count << " of " << x.N() << std::endl;
            }


          if (verbose >= VERBOSITY_DEBUG_LEVEL) {
            logger << "sum_alpha = " << sum_alpha << std::endl;
          }
          sum_alpha = gv.comm().sum( sum_alpha );
          sum_beta = gv.comm().sum( sum_beta );
          if (verbose >= VERBOSITY_DEBUG_LEVEL) {
            logger << "total_error(all) = " << total_error << std::endl;
            logger << "sum_alpha(all) = " << sum_alpha << std::endl;
            logger << "alpha_fraction(all) = " << sum_alpha/total_error << std::endl;
            logger << "sum_beta(all) = " << sum_beta << std::endl;
            logger << "beta_fraction(all) = " << sum_beta/total_error << std::endl;
          }

          if (std::abs(alpha-sum_alpha/total_error) <= 0.01 
              && std::abs(beta-sum_beta/total_error) <= 0.01) break;
          
          if (sum_alpha > alpha*total_error)
            eta_alpha_left = eta_alpha;
          else
            eta_alpha_right = eta_alpha;
          
          if (sum_beta > beta*total_error) 
            eta_beta_right = eta_beta;
          else
            eta_beta_left = eta_beta;

        }

      if (verbose >= VERBOSITY_LS_DETAILS) {
        logger << " eta_alpha = " << eta_alpha 
               << " eta_beta  = " << eta_beta << std::endl;
      }

      // eta_alpha = gv.comm().max( eta_alpha );
      // eta_beta = gv.comm().min( eta_beta );

      if( beta < 1E-12 )
        eta_beta = 0.0;

      if( alpha < 1E-12 )
        eta_alpha = 0.0;

      if (verbose >= VERBOSITY_LS_DETAILS) {
        logger << "comm: eta_alpha = " << eta_alpha 
               << " eta_beta  = " << eta_beta << std::endl;
      }
    }


    template<typename T>
    void element_fraction( const T& x,
                           typename T::ElementType alpha,
                           typename T::ElementType beta,
                           typename T::ElementType& eta_alpha,
                           typename T::ElementType& eta_beta,
                           int verbose=0 )
    {
      const int steps=20; // max number of bisection steps
      typedef typename T::ElementType NumberType;
      NumberType total_error =x.N();
      NumberType max_error = x.infinity_norm();
      NumberType eta_alpha_left = 0.0;
      NumberType eta_alpha_right = max_error;
      NumberType eta_beta_left = 0.0;
      NumberType eta_beta_right = max_error;
      for (int j=1; j<=steps; j++)
        {
          eta_alpha = 0.5*(eta_alpha_left+eta_alpha_right);
          eta_beta = 0.5*(eta_beta_left+eta_beta_right);
          NumberType sum_alpha=0.0;
          NumberType sum_beta=0.0;
          unsigned int alpha_count = 0;
          unsigned int beta_count = 0;
          for (unsigned int i=0; i<x.N(); i++)
            {
              if (x.base()[i]>=eta_alpha) { sum_alpha += 1.0; alpha_count++;}
              if (x.base()[i]< eta_beta) { sum_beta +=1.0; beta_count++;}
            }
          if (verbose>=VERBOSITY_DEBUG_LEVEL)
            {
              logger << j 
                     << " eta_alpha=" << eta_alpha 
                     << " alpha_fraction=" << sum_alpha/total_error 
                     << " elements: " << alpha_count << " of " << x.N() << std::endl;
              logger << j 
                     << " eta_beta=" << eta_beta 
                     << " beta_fraction=" << sum_beta/total_error 
                     << " elements: " << beta_count << " of " << x.N() << std::endl;
            }
          if (std::abs(alpha-sum_alpha/total_error) <= 0.01 && std::abs(beta-sum_beta/total_error) <= 0.01) break;
          if (sum_alpha>alpha*total_error)
            eta_alpha_left = eta_alpha;
          else
            eta_alpha_right = eta_alpha;
          if (sum_beta>beta*total_error) 
            eta_beta_right = eta_beta;
          else
            eta_beta_left = eta_beta;
        }

      if( beta < 1E-12 )
        eta_beta = 0.0;

      if( alpha < 1E-12 )
        eta_alpha = 0.0;

      if (verbose >= VERBOSITY_LS_DETAILS)
        {
          logger << "+++ refine_threshold=" << eta_alpha 
                 << " coarsen_threshold=" << eta_beta << std::endl;
          std::cout << "+++ refine_threshold=" << eta_alpha 
                    << " coarsen_threshold=" << eta_beta << std::endl;
        }

    }

    /** Compute error distribution
     */
    template<typename T>
    void error_distribution(const T& x, int bins)
    {
      const int steps=30; // max number of bisection steps
      typedef typename T::ElementType NumberType;
      NumberType total_error = x.one_norm();
      NumberType total_elements = x.N();
      NumberType max_error = x.infinity_norm();
      std::vector<NumberType> left(bins,0.0);
      std::vector<NumberType> right(bins,max_error*(1.0+1e-8));
      std::vector<NumberType> eta(bins);
      std::vector<NumberType> target(bins);
      for (unsigned int k=0; k<bins; k++)
        target[k]= (k+1)/((NumberType)bins);
      for (int j=1; j<=steps; j++)
        {
          for (unsigned int k=0; k<bins; k++)
            eta[k]= 0.5*(left[k]+right[k]);
          std::vector<NumberType> sum(bins,0.0);
          std::vector<int> count(bins,0);
          for (unsigned int i=0; i<x.N(); i++)
            for (unsigned int k=0; k<bins; k++)
              if (x[i]<=eta[k])
                {
                  sum[k] += x[i];
                  count[k] += 1;
                }
          // std::cout << std::endl;
          // std::cout << "// step " << j << std::endl;
          // for (unsigned int k=0; k<bins; k++)
          //    std::cout << k+1 << " " << count[k] << " " << eta[k] << " " << right[k]-left[k] 
          //          << " " << sum[k]/total_error << " " << target[k] << std::endl;
          for (unsigned int k=0; k<bins; k++)
            if (sum[k]<=target[k]*total_error) 
              left[k] = eta[k];
            else
              right[k] = eta[k];
        }
      std::vector<NumberType> sum(bins,0.0);
      std::vector<int> count(bins,0);
      for (unsigned int i=0; i<x.N(); i++)
        for (unsigned int k=0; k<bins; k++)
          if (x[i]<=eta[k])
            {
              sum[k] += x[i];
              count[k] += 1;
            }
      std::cout << "// error distribution" << std::endl;
      std::cout << " number of elements: " << x.N() << std::endl;
      std::cout << " max element error:  " << max_error << std::endl;
      std::cout << " total error:        " << total_error << std::endl;
      std::cout << " bin #elements eta sum/total " << std::endl;
      for (unsigned int k=0; k<bins; k++)
        std::cout << k+1 << " " << count[k] << " " << eta[k] << " " << sum[k]/total_error << std::endl;
    }

    template<typename Grid, typename X>
    void mark_grid (Grid &grid, const X& x, typename X::ElementType refine_threshold, 
                    typename X::ElementType coarsen_threshold, int verbose=0)
    {
      typedef typename Grid::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator 
        Iterator;
      typedef typename Grid::LeafGridView GV;
      typedef typename GV::IndexSet IndexSet;

      const GV& gv=grid.leafView();
      const IndexSet& is(gv.indexSet());
      Iterator it = grid.template leafbegin<0,Dune::All_Partition>();
      Iterator eit = grid.template leafend<0,Dune::All_Partition>();

      unsigned int refine_cnt=0;
      unsigned int coarsen_cnt=0;

      for(;it!=eit;++it)
        {
          typename IndexSet::IndexType myid = is.template index<0>(*it);
          if (x[myid]>=refine_threshold)
            {
              grid.mark(1,*(it));
              refine_cnt++;
            }
          if (x[myid]<=coarsen_threshold)
            {
              grid.mark(-1,*(it));
              coarsen_cnt++;
            }
        }
      if (verbose>0)
        std::cout << "+++ mark_grid: " << refine_cnt << " marked for refinement, " 
                  << coarsen_cnt << " marked for coarsening" << std::endl;
    }

    template<typename Grid, typename GV, typename X>
    void mark_grid_rgv( Grid &grid,
                        const GV &gv,
                        const X& x,
                        typename X::ElementType refine_threshold,
                        typename X::ElementType coarsen_threshold,
                        int verbose=0
                        )
    {
      typedef typename Grid::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator 
        Iterator;
      typedef typename GV::IndexSet IndexSet;

      const IndexSet& is(gv.indexSet());
      Iterator it = grid.template leafbegin<0,Dune::All_Partition>();
      Iterator eit = grid.template leafend<0,Dune::All_Partition>();

      unsigned int refine_cnt=0;
      unsigned int coarsen_cnt=0;

      for(;it!=eit;++it) {
        typename IndexSet::IndexType myid = is.template index<0>(*it);
        if (x.base()[myid]>=refine_threshold){
          grid.mark(1,*(it));
          refine_cnt++;
        }
        if (x.base()[myid]<=coarsen_threshold && it->level()>0 ){
          grid.mark(-1,*(it));
          coarsen_cnt++;
        }
      }

      if (verbose >= VERBOSITY_LS_DETAILS){
        std::cout << "mark_grid_rgv: refine_cnt = " << refine_cnt << std::endl;
        std::cout << "mark_grid_rgv: coarsen_cnt = " << coarsen_cnt << std::endl;

        logger << "mark_grid_rgv: refine_cnt = " << refine_cnt << std::endl;
        logger << "mark_grid_rgv: coarsen_cnt = " << coarsen_cnt << std::endl;

      }
    }


    template<typename Grid>
    void mark_all_coarsen( Grid &grid ) {
      typedef typename Grid::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator 
        Iterator;
      Iterator it = grid.template leafbegin<0,Dune::All_Partition>();
      Iterator eit = grid.template leafend<0,Dune::All_Partition>();
      int maxLevel = grid.maxLevel();
      if( maxLevel > 0 ){
        for(;it!=eit;++it) {
          if( it->level() == maxLevel )
            grid.mark(-1,*(it));
        }
      }
    }


  } // namespace Gesis
} // namespace Dune

#endif
