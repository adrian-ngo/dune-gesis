#ifndef TOTAL_MASS_HH
#define	TOTAL_MASS_HH

#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <assert.h>
#include <sstream>

#include <dune/common/timer.hh>


extern CLogfile logger;

namespace Dune {
  namespace Gesis {

    class GridFunctionTools {

    private:

      GridFunctionTools(){};

    public:

      // Loop over grd cells and compute integral of the solution.
      
      template<typename GFS,typename VCType>
      static void totalMass( const GFS& gfs
                             , const VCType& xSolution
                             , REAL& negMass
                             , REAL& posMass
                             , REAL& totMass
                             ) {
        Dune::Timer watch;

        typedef typename GFS::Traits::GridViewType GV;
        typedef typename GV::Grid::ctype DF;
        const int dim = GV::dimension;
        //typedef typename GV::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementIterator;
        //typedef typename GV::Traits::template Codim<0>::Entity Element;
        //typedef typename GV::IntersectionIterator IntersectionIterator;
        //typedef typename IntersectionIterator::Intersection Intersection;

        GV gv = gfs.gridView();
        typedef Dune::PDELab::DiscreteGridFunction<GFS,VCType> DGF;
        DGF dgf(gfs,xSolution);

        REAL value = 0;
        REAL posvalue = 0;
        REAL negvalue = 0;

        // element loop
        for( auto eit=gv.template begin<0,Dune::All_Partition>()
               ; eit!=gv.template end<0,Dune::All_Partition>()
               ; ++eit) {

          if(eit->partitionType()==Dune::InteriorEntity){

            // select quadrature rule
            Dune::GeometryType geometrytype = eit->geometry().type();
            const int integration_order = 2*pMAX;
            const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF, dim>::rule(geometrytype, integration_order);
      
            // quadrature points
            Dune::FieldVector<REAL,1> fval(0.0);

            // loop over quadrature points 
            for( typename Dune::QuadratureRule<DF,dim>::const_iterator iterQP = rule.begin()
                   ; iterQP != rule.end()
                   ; ++iterQP) {

              // Floating Point Exception in 3D parallel is caused here, when exception handler is activated!
        
              // evaluate the darcy flux
              fval = 0.0;
              dgf.evaluate(*eit, iterQP->position(), fval);

              // integration factor
              REAL factor = iterQP->weight() * eit->geometry().integrationElement(iterQP->position());
            
              // add to the sensitivity
              if(fval>1E-12)
                posvalue += fval * factor;
              else
                negvalue += fval * factor;

              value += fval * factor;

            } // quadrature loop

          } // end if interior

        } // element loop

        posMass = gv.comm().sum( posvalue );
        negMass = gv.comm().sum( negvalue );
        totMass = gv.comm().sum( value );

        if( General::verbosity >= VERBOSITY_DEBUG_LEVEL ){
          std::cout << "P" << gv.comm().rank() 
                    << ": pos./neg. mass: " 
                    << posvalue << "/" << negvalue 
                    << std::endl;
          std::cout << "P" << gv.comm().rank() << ": tot. mass: " << value << std::endl;
        }
        
        return;

      } // totalMass

    }; // class VTKPlot

  }

}

#endif // TOTAL_MASS_HH
