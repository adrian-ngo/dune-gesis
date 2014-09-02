/* 
 * File:   ALUTools.hh
 * Author: ngo
 *
 * Created on Jan 23, 2013
 */

#ifndef DUNE_GEOINVERSION_ALU_TOOLS_HH
#define	DUNE_GEOINVERSION_ALU_TOOLS_HH

#include <dune/pdelab/common/geometrywrapper.hh>
#include <assert.h>
#include <sstream>


extern CLogfile logger;


namespace Dune {
  namespace GeoInversion {

    class ALUTools{
    private:
      ALUTools(){};
    public:
      // Write a simple 2D Alugrid input file
      static void output_Domain_ALU( const std::string filename,
                                     double domain_length_x,
                                     double domain_length_y ){

        std::ofstream outfile( filename.c_str(), std::ios::out );
        if( outfile.is_open() ){
          outfile << "!Triangles" << std::endl;
          outfile << "4" << std::endl;
          outfile << "0 0 " << std::endl;
          outfile << "0 " << domain_length_y << " " << std::endl;
          outfile << domain_length_x << " 0 " << std::endl;
          outfile << domain_length_x << " " << domain_length_y << " " << std::endl;
          outfile << "2" << std::endl;
          outfile << " 2  3  0 " << std::endl;
          outfile << " 1  0  3 " << std::endl;
          outfile << "4" << std::endl;
          outfile << "-1 1 0 " << std::endl;
          outfile << "-1 0 2 " << std::endl;
          outfile << "-1 3 1 " << std::endl;
          outfile << "-1 2 3 " << std::endl;
          outfile.close();
        }
      };

      
      static void create_ALU_Hexahedra_Cuboid( const std::string filename,
                                               double domain_length_x,
                                               double domain_length_y, 
                                               double domain_length_z ){


        std::ofstream outfile( filename.c_str(), std::ios::out );
        if( outfile.is_open() ){
          outfile << "!Hexahedra" << std::endl;
          outfile << "8" << std::endl;

          outfile << "0.0  0.0  0.0" << std::endl;
          outfile << domain_length_x << " 0.0 0.0 " << std::endl;
          outfile << domain_length_x << " " << domain_length_y << " 0.0 " << std::endl;
          outfile << "0.0 " << domain_length_y << " 0.0 " << std::endl;
          outfile << "0.0 0.0 " << domain_length_z << std::endl;
          outfile << domain_length_x << " 0.0 " << domain_length_z << std::endl;
          outfile << domain_length_x << " " << domain_length_z << " " << domain_length_z << std::endl;
          outfile << "0.0 " << domain_length_y << " " << domain_length_z << std::endl;

          outfile << "1" << std::endl;
          outfile << " 0 1 2 3 4 5 6 7 " << std::endl;

          outfile << "6" << std::endl;
		  outfile << "-2 4 0 3 7 4" << std::endl; 
		  outfile << "-3 4 1 5 6 2" << std::endl; 
		  outfile << "-1 4 0 4 5 1" << std::endl; 
		  outfile << "-1 4 3 2 6 7" << std::endl; 
		  outfile << "-1 4 0 1 2 3" << std::endl; 
		  outfile << "-1 4 5 4 7 6" << std::endl; 

          outfile.close();
        }


        /*

1
0  1  2  3  4  5  6  7 

6
-2 4 0 3 7 4
-3 4 1 5 6 2
-1 4 0 4 5 1 
-1 4 3 2 6 7
-1 4 0 1 2 3 
-1 4 5 4 7 6 

0 -1
1 -1
2 -1
3 -1
4 -1
5 -1
6 -1
7 -1
        */

      }


    }; // class ALUTools




  } // GeoInversion
} // Dune








#endif	/* DUNE_GEOINVERSION_ALU_TOOLS_HH */

