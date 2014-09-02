#ifndef REGULARIZE_HH
#define REGULARIZE_HH



extern CLogfile logger; // declared + initalized in the main function!

namespace Dune {
  namespace GeoInversion {

    class Regularizer {

    private:
      Regularizer(){};

    public:
      /* 1d function using Gaussian-like exponential function to smoothen the 
         transition from zero up to a certain non-zero value.
         
         x : variable of the function
         y1: start position of non-zero section
         y2: end position of non-zero section
         gDelta : transition zone
         concentration : the non-zero value
         
         Note: 
         This function is supposed to be used for a 2D problem 
         with a 1D inflow boundary!

      */
      template<typename RF>
      static RF gRegular1( RF x, 
                           RF y1, 
                           RF y2, 
                           RF delta_y, 
                           RF gDelta, 
                           RF concentration ) {
  
        RF gSigma=0.5*gDelta;
        RF gP=2.8;

        if( x >= y1 && x <= y2 )
          return concentration;
    
#ifdef REGULARIZED_DIRICHLET
        else if( x > y2 && x < y2+gDelta ){
          x -= y2;
          return concentration * exp(-std::pow(std::abs(x)/gSigma,gP));
        }
        else if( x < y1 && x > y1-gDelta ){
          x -= y1;
          return concentration * exp(-std::pow(std::abs(x)/gSigma,gP));
        }
#endif
        else
          return RF(0.0);

      } // gRegular1()






      template<typename RF>
      RF gRegular2( RF x, RF yStart, RF yEnd, RF delta_y, RF concentration ) {
        if( 
           x > yStart - 1e-6
           && 
           x < yEnd + 1e-6
            )
          {
            return concentration;
          }
#ifdef REGULARIZED_DIRICHLET
        else if( 
                x >= yStart - 0.9 * delta_y
                &&
                x <= yEnd + 0.9 * delta_y
                 )
          {
            return concentration * 5.0 / 6.0;
          }
        else if( 
                x >= yStart - 1.9 * delta_y
                &&
                x <= yEnd + 1.9 * delta_y
                 )
          {
            return concentration * 0.50;
          }
        else if( 
                x >= yStart - 2.9 * delta_y
                &&
                x <= yEnd + 2.9 * delta_y
                 )
          {
            return concentration / 6.0;
          }
#endif // REGULARIZED_DIRICHLET
        else 
          {
            return 0.0;
          }
      } // gRegular2()




      template<typename RF>
      static RF gRegularYZ( RF y, RF z, 
                            RF y1, RF z1, 
                            RF y2, RF z2, 
                            RF gDeltaY, RF gDeltaZ, 
                            RF concentration,
                            bool bfixedwidth,
                            RF regularization_factor
                            ) {
        
        if( regularization_factor < 1E-12 ){
          if( y >= y1 && y <= y2 
              &&
              z >= z1 && z <= z2 ){
            return concentration;
          }
          return 0;            
        }

        
        if( bfixedwidth ){
          gDeltaY=regularization_factor;
          gDeltaZ=regularization_factor;
        }

        RF gSigmaY=0.5*gDeltaY;
        RF gSigmaZ=0.5*gDeltaZ;

        // share the smearing to both sides of the jump:
        y1 += 0.5*gDeltaY;
        y2 -= 0.5*gDeltaY;

        z1 += 0.5*gDeltaZ;
        z2 -= 0.5*gDeltaZ;


        RF gP=2.8;
        if( y >= y1 && y <= y2 ){
          if( z >= z1 && z <= z2 ){
            return concentration;
          }
          else if( z<z1 && z>z1-gDeltaZ  ){
            z -= z1;
            return concentration * exp(-std::pow(std::abs(z)/gSigmaZ,gP));
          }
          else if( z>z2 && z<z2+gDeltaZ  ){
            z -= z2;
            return concentration * exp(-std::pow(std::abs(z)/gSigmaZ,gP));
          }
          else
            return 0;
        }
        else if( y<y1 && y>y1-gDeltaY ){
          if( z >= z1 && z <= z2 ){
            y -= y1;
            return concentration * exp(-std::pow(std::abs(y)/gSigmaY,gP));
          }
          else if( z<z1 && z>z1-gDeltaY  ){
            y -= y1;
            z -= z1;
            return concentration * exp(-std::pow(std::abs(y)/gSigmaY,gP)) * exp(-std::pow(std::abs(z)/gSigmaZ,gP));
          }
          else if( z>z2 && z<z2+gDeltaZ  ){
            y -= y1;
            z -= z2;
            return concentration * exp(-std::pow(std::abs(y)/gSigmaY,gP)) * exp(-std::pow(std::abs(z)/gSigmaZ,gP));
          }
          else
            return 0;
        }
        else if( y>y2 && y<y2+gDeltaY ){
          if( z >= z1 && z <= z2 ){
            y -= y2;
            return concentration * exp(-std::pow(std::abs(y)/gSigmaY,gP));
          }
          else if( z<z1 && z>z1-gDeltaZ  ){
            y -= y2;
            z -= z1;
            return concentration * exp(-std::pow(std::abs(y)/gSigmaY,gP)) * exp(-std::pow(std::abs(z)/gSigmaZ,gP));
          }
          else if( z>z2 && z<z2+gDeltaZ  ){
            y -= y2;
            z -= z2;
            return concentration * exp(-std::pow(std::abs(y)/gSigmaY,gP)) * exp(-std::pow(std::abs(z)/gSigmaZ,gP));
          }
          else
            return 0;
        }  
        else
          return RF(0.0);
      }
      



      /*
        Note: In the inputfile, you can specify
        
        regularization_type="dynamic" or
        regularization_type="fixed".
        
        With regularization_factor you can set the regularization zone
        to a fixed value (if regularization_type="fixed") or to a 
        multiple of the gridsize in that dimension (if regularization_type="dynamic").

       */


      template<typename RF>
      static RF gRegularY( RF y,
                           RF y1,
                           RF y2, 
                           RF gDeltaY, 
                           RF concentration,
                           bool bfixedwidth,
                           RF regularization_factor
                           ) {

        if( regularization_factor < 1E-12 ){
          if( y >= y1 && y <= y2 ){
            return concentration;
          }
          return 0;            
        }
        
        if( bfixedwidth ){
          gDeltaY=regularization_factor;
        }

        RF gSigmaY=0.5*gDeltaY;

        // share the smearing to both sides of the jump:
        y1 += 0.5*gDeltaY;
        y2 -= 0.5*gDeltaY;

        RF gP=2.8;
        if( y >= y1 && y <= y2 ){
          return concentration;
        }
        else if( y<y1 && y>y1-gDeltaY ){
          y -= y1;
          return concentration * exp(-std::pow(std::abs(y)/gSigmaY,gP));
        }
        else if( y>y2 && y<y2+gDeltaY ){
          y -= y2;
          return concentration * exp(-std::pow(std::abs(y)/gSigmaY,gP));
        }  
        else
          return RF(0.0);
      }






    }; // class Regularizer

  } // GeoInversion
  
} // Dune






#endif // REGULARIZE_HH

