#ifndef DUNE_GESIS_DELTA_SD_HH
#define DUNE_GESIS_DELTA_SD_HH

// Elman, Silvester, Wathen (2005)
REAL zeta_1( const REAL Pe)
{
  if(Pe > 1.0){
    return 1.0 - 1.0/Pe ;
  } else {
    return 0;
  }
}

// Cirpka(2000) bzw. Brooks&Hughes(1982)
REAL zeta_11( const REAL Pe)
{
  return 1.0/std::tanh(Pe) - 1.0/Pe;
}

REAL zeta_12( const REAL Pe)
{
  return ( 1.0 / tanh( Pe ) - 1.0 / Pe );
}



// The following variations of the zeta functions are taken from the dissertation of "Areti Papastavrou" (1998)

// doppelt asymtotisch:
REAL zeta_4( const REAL Peclet )
{
  return std::min( 1.0, Peclet / 3.0 );
}

// Mizukami:
REAL zeta_5( const REAL Peclet )
{
  return Peclet / (1.0 + Peclet);
}

// kritisch:
REAL zeta_6( const REAL Peclet )
{
  return std::max( 0.0, 1.0 - 1.0 / Peclet );
}

// Johnson:
REAL zeta_7( const REAL Peclet )
{
  return std::max( 0.0, 1.0 - 0.5 / Peclet );
}

// Hughes, Tezduyar:
REAL zeta_8( const REAL Peclet )
{
  REAL a = Peclet / 3.0;
  return a * 1.0 / sqrt( 1.0 + a*a );
}


// Tezduyar:
REAL zeta_9( const REAL Peclet )
{
  REAL a = Peclet;
  return a * 1.0 / sqrt( 1.0 + a*a );
}



#endif // DUNE_GESIS_DELTA_SD_HH

