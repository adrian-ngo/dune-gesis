// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef WELLPOSITION_HH
#define WELLPOSITION_HH


enum WellPosition {w_outside=0,w_inside=1,w_top=2,w_bottom=3};

// Check whether this element is (partially) penetrated by one of the wells
template<typename EG,typename IDT,typename SDT>
WellPosition isElementPenetratedByWell( const EG& eg,
                                        const int iWell,
                                        const IDT& inputdata,
                                        const SDT& setupdata
                                        ) {

  const int dim = EG::Geometry::dimension;
  Dune::FieldVector<REAL,dim> cellcenter = eg.geometry().center();
  //const int refinementlevel = eg.entity().level();

  Dune::FieldVector<REAL,dim-1> well_position;
  well_position[0] = setupdata.wdlist.pointdata_vector[iWell].x; // + GEO_EPSILON; // Avoid being directly on an intersection! This could cause ambiguities!
#ifdef DIMENSION3
  well_position[1] = setupdata.wdlist.pointdata_vector[iWell].y;// + GEO_EPSILON; // Avoid being directly on an intersection! This could cause ambiguities!
#endif

  REAL well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
  REAL well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;


  // define element zone and check if well is within element zone

  const REAL dx = inputdata.domain_data.baselevel_gridsizes[0]; //   / std::pow(2.0,refinementlevel);
  REAL x1 = cellcenter[0] - 0.5 * dx;
  REAL x2 = cellcenter[0] + 0.5 * dx;
#ifdef DIMENSION3
  const REAL dy = inputdata.domain_data.baselevel_gridsizes[1]; //   / std::pow(2.0,refinementlevel);
  REAL y1 = cellcenter[1] - 0.5 * dy;
  REAL y2 = cellcenter[1] + 0.5 * dy;
#endif
  const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1]; //   / std::pow(2.0,refinementlevel);
  REAL z1 = cellcenter[dim-1] - 0.5 * dz;
  REAL z2 = cellcenter[dim-1] + 0.5 * dz;

  if( well_position[0] > x1
      &&
      well_position[0] <= x2
#ifdef DIMENSION3
      &&
      well_position[1] > y1
      &&
      well_position[1] <= y2
#endif
      &&
      well_top > z1 + GEO_EPSILON  // This GEO_EPSILON is important here!
      &&
      well_bottom < z2 - GEO_EPSILON // This GEO_EPSILON is important here!
      ) {

    if( well_top <= z2 )
      return w_top;
    else if( well_bottom >= z1 )
      return w_bottom;
    else
      return w_inside;
  }
  else
    return w_outside;
}


template<typename EG,typename IDT,typename SDT>
WellPosition isElementWithinWellZone( const EG& eg,
                                      const int iWell,
                                      const IDT& inputdata,
                                      const SDT& setupdata
                                      ) {

  const int dim = EG::Geometry::dimension;
  Dune::FieldVector<REAL,dim> cellcenter = eg.geometry().center();
  //const int refinementlevel = eg.entity().level();

  Dune::FieldVector<REAL,dim-1> well_center;
  well_center[0] = setupdata.wdlist.pointdata_vector[iWell].x;

  // define well zone and check if element center is within well zone

  const REAL dx = inputdata.domain_data.baselevel_gridsizes[0];
  REAL x1 = well_center[0] - 0.5 * dx;
  REAL x2 = well_center[0] + 0.5 * dx;
#ifdef DIMENSION3
  well_center[1] = setupdata.wdlist.pointdata_vector[iWell].y;
  const REAL dy = inputdata.domain_data.baselevel_gridsizes[1];
  REAL y1 = well_center[1] - 0.5 * dy;
  REAL y2 = well_center[1] + 0.5 * dy;
#endif
  const REAL dz = inputdata.domain_data.baselevel_gridsizes[dim-1];

  REAL well_top =  setupdata.wdlist.pointdata_vector[iWell].well_top;
  REAL well_bottom =  setupdata.wdlist.pointdata_vector[iWell].well_bottom;

  if( cellcenter[0] > x1
      &&
      cellcenter[0] < x2
#ifdef DIMENSION3
      &&
      cellcenter[1] > y1
      &&
      cellcenter[1] < y2
#endif
      &&
      cellcenter[dim-1] > well_bottom
      &&
      cellcenter[dim-1] < well_top
      ) {
    if( well_top - cellcenter[dim-1] < 0.6*dz )
      return w_top;
    else if( cellcenter[dim-1] - well_bottom < 0.6*dz )
      return w_bottom;
    else
      return w_inside;
  }
  else
    return w_outside;

}



#endif // WELLPOSITION_HH
