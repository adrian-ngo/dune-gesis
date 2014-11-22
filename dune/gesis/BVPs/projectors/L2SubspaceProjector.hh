#ifndef DUNE_GESIS_L2_SUBSPACE_PROJECTOR_HH
#define DUNE_GESIS_L2_SUBSPACE_PROJECTOR_HH



#include "L2ProjectionOperator.hh"
#include "L2Projection.hh"

namespace Dune{
  namespace Gesis{

    template<typename GFS_DG, // source solution space
             typename GFS_CG, // target solution space
             typename IDT,
             int dim>
    class L2SubspaceProjector{

    private:
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_DG,REAL>::Type VCType_DG;
      typedef Dune::PDELab::DiscreteGridFunction<GFS_DG,VCType_DG> DGF_DG;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS_CG,REAL>::Type VCType_CG;

      const GFS_DG& gfs_dg;
      const GFS_CG& gfs_cg;
      const IDT& inputdata;

    public:
      L2SubspaceProjector(const GFS_DG& gfs_dg_,
                          const GFS_CG& gfs_cg_,
                          const IDT& inputdata_):
        gfs_dg(gfs_dg_),
        gfs_cg(gfs_cg_),
        inputdata(inputdata_)
      {
        //std::cout << "DEBUG LOG: L2SubspaceProjector: gfs_cg.maxLocalSize() = " << gfs_cg.maxLocalSize() << std::endl;
        //std::cout << "DEBUG LOG: L2SubspaceProjector: gfs_cg.size() = " << gfs_cg.size() << std::endl;
      }

      void apply( const VCType_DG& v,
                  VCType_CG& u,
                  const REAL l2_diffusion = 1.0
                  ) {
        Dune::Timer watch;
#ifdef L2ProjectionOfM0
        u = 0;
        DGF_DG dgf_v(gfs_dg,v);
        L2Projection<GFS_CG,DGF_DG,IDT> l2proj(gfs_dg.gridView(),gfs_cg,dgf_v,l2_diffusion,inputdata);

        General::log_elapsed_time( watch.elapsed(),
                                   gfs_dg.gridView().comm(),
                                   inputdata.verbosity,
                                   "L2Proj",
                                   "Matrix Pattern + Assembly" );

        l2proj.solve_forward( u );
        if( inputdata.verbosity>=VERBOSITY_EQ_SUMMARY && gfs_dg.gridView().comm().rank()==0 )
          std::cout << l2proj.show_ls_result() << std::endl;
        return;
#else
        u = v;
        return;
#endif

      }

    };

  }
}

#endif // DUNE_GESIS_L2_SUBSPACE_PROJECTOR_HH
