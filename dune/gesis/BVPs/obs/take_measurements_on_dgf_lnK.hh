/*
 * This function is to take simulated measurement of lnK
 */

// Inputs:
// gv: grid view
// K: current estimate of ln K field (YFieldGenerator)
// measuring_points: container for the simulated measurements
// helper: MPI utility

#ifndef DUNE_GESIS_TAKE_MEASUREMENTS_ON_DGF_LNK_HH
#define	DUNE_GESIS_TAKE_MEASUREMENTS_ON_DGF_LNK_HH

template<typename GV, typename YFieldGenerator, typename VEC>
void take_measurements_on_dgf_lnK(const GV& gv,
				  const YFieldGenerator& yfg,
				  VEC& measuring_points,
				  const Dune::MPIHelper& helper) {

  const UINT dim = GV::dimension;
  const UINT nPoints = measuring_points.size();
  typedef typename GV::template Codim < 0 > ::Iterator LeafIterator;

  //  logger << "Taking measurements (lnK)... " << std::endl;


  //const typename GV::IndexSet& is = gv.indexSet();

  // loop over all elements of the grid!
  for (LeafIterator it = gv.template begin < 0 > ()
	 ; it != gv.template end < 0 > ()
	 ; ++it) {

    //int id = is.index(*it);

    // Loop over list of measuring points
    for (UINT index = 0; index < nPoints; index++) {

      Dune::FieldVector<double, dim> currentpoint(0.0);
      currentpoint[0] = measuring_points[index].x;
      currentpoint[1] = measuring_points[index].y;
#ifdef DIMENSION3
      currentpoint[2] = measuring_points[index].z;
#endif

      // Get the local coordinate:
      Dune::FieldVector<double, dim> currentpoint_local = it->geometry().local(currentpoint);

      for(UINT ii=0; ii<dim; ii++){
        if(std::fabs(std::floor(currentpoint_local[ii]+0.5)- currentpoint_local[ii])<1e-10){
          currentpoint_local[ii]=std::floor(currentpoint_local[ii]+0.5);
        }
      }

      // Check if the current point belongs to that cell
      if ( it->partitionType()!=Dune::OverlapEntity && (currentpoint_local[0] >= 0.0 && currentpoint_local[0] < 1.0
                                                        && currentpoint_local[1] >= 0.0 && currentpoint_local[1] < 1.0
#ifdef DIMENSION3
                                                        && currentpoint_local[2] >= 0.0 && currentpoint_local[2] < 1.0
#endif
                                                        )) {


	// get the lnK-value at the measurement location

	yfg.evaluateY( currentpoint, measuring_points[index].value );
        //measuring_points[index].value = (double) K[id];


      }

    }  // end of loop over measurement points

  } // end of loop over leaf elements

}

#endif
