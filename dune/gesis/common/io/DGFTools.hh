/* 
 * File:   DGFTools.hh
 * Author: ngo
 *
 * Created on Jan 26, 2013
 */

#ifndef DUNE_GESIS_DGF_TOOLS_HH
#define	DUNE_GESIS_DGF_TOOLS_HH

#include <dune/pdelab/common/geometrywrapper.hh>
#include <assert.h>
#include <sstream>


extern CLogfile logger;


namespace Dune {
  namespace Gesis {


    /*
      DGF - Dune Grid Format 
    */
    class DGFTools{

    private:
      DGFTools(){};
    public:
      
      /* 
         Creates a dgf (Dune Grid Format) file for initial generation of a unstructured grid 
         over a cuboid or rectangular domain. 
         On the one hand, a minimal number of cells given by 
         npfactors[0] * npfactors[1] * npfactors[2] 
         is needed to serve all nProcesses.
         On the other hand, the number of cells per dimension, npfactors[i],
         must fit the number of virtual cells nCells[i], 
         i.e. npfactors[i] must be an integer multiple of nCells[i].
         Otherwise, the Y-field data would not fit into the unstructured grid smoothly!
      */
      static void output_Domain_DGF( 
                             const std::string& filename
                             , const Vector<REAL>& extensions
                             , const Vector<UINT>& nCells
                             , const std::string& elementshape
                             , const UINT heapSize
                             , const UINT nProcesses
                              )
      {

#ifdef DIMENSION3
        UINT dim = 3;
        Vector<UINT> npfactors( 3, 1 );
        General::factorize_three_optimally( nProcesses, npfactors[0], npfactors[1], npfactors[2] );
#else
        UINT dim = 2;
        Vector<UINT> npfactors( 2, 1 );
        General::factorize_two_optimally( nProcesses, npfactors[0], npfactors[1] );
#endif


        UINT nAllVirtualCells = nCells[ 0 ] * nCells[ 1 ];
#ifdef DIMENSION3
        nAllVirtualCells *= nCells[ 2 ];
#endif


        /* nAllVirtualCells >> nProcesses should be the normal case */
        if( nAllVirtualCells > nProcesses )
          {
            std::cout << "nAllVirtualCells = " << nAllVirtualCells  << std::endl;
            for( UINT i=0; i<dim; i++ )
              {
                /* If npfactors is not an integer multiple of nCells, make it at least the same number! 
                   This is necessary to make the unstructured grid fit into the virtual grid.
                   This is to prevent cubes (containing two triangular grid cells) from getting 
                   two different values of the Y-field!
                */
                if( !General::isFactorOf(  nCells[i], npfactors[i] ) )
                  {
                    npfactors[ i ] = nCells[ i ];
                  }
              }
          }
  
        logger << "create DGF file for the domain" << std::endl;
        std::ofstream outfile( filename.c_str(), std::ios::out );
        if( outfile.is_open() )
          {

            outfile << "DGF" << std::endl;

            outfile << "Interval" << std::endl;

            outfile << "0.0 0.0" 
#ifdef DIMENSION3
                    << " 0.0"
#endif
                    << std::endl;

            outfile << extensions[ 0 ] << " " << extensions[ 1 ] 
#ifdef DIMENSION3
                    << " " << extensions[ 2 ]
#endif
                    << std::endl;


            outfile << npfactors[0] << " " << npfactors[1]
#ifdef DIMENSION3
                    << " " << npfactors[2]
#endif
                    << std::endl;

            outfile << "#" << std::endl;

            if( elementshape=="cube" )
              outfile << "Cube" << std::endl;
            else
              outfile << "Simplex" << std::endl;

            outfile << "#" << std::endl;

            outfile << std::endl;

            outfile << "GridParameter" << std::endl;
            outfile << "% set overlap to 1" << std::endl;
            outfile << "overlap 1" << std::endl;
            outfile << "% set closure to NONE" << std::endl;
            outfile << "closure NONE" << std::endl;

            outfile << "% set copies to YES" << std::endl;
            outfile << "copies yes" << std::endl;
            //outfile << "% set copies to NO" << std::endl;
            //outfile << "copies no" << std::endl;

            outfile << "% set heapsize " << std::endl;
            outfile << "heapsize " << heapSize << std::endl;
            outfile << "#" << std::endl;
	  
            outfile << std::endl;
            outfile << "Boundarydomain" << std::endl;
            outfile << "default 1" << std::endl;
            outfile << "#" << std::endl;
            outfile << "#" << std::endl;
            outfile << "# " << filename.c_str() << std::endl;

            outfile.close();
            logger << "DGF file written to " << filename.c_str() << std::endl;
          }
  
      } // static void output_Domain_DGF

    }; // class DGFTools






  } // Gesis
} // Dune








#endif	/* DUNE_GESIS_DGF_TOOLS_HH */

