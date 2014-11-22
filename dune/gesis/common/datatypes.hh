#ifndef DUNE_GESIS_DATATYPES_HH
#define DUNE_GESIS_DATATYPES_HH

// define some general data types

#ifdef _USE_FLOAT_
typedef float CTYPE;  // Be careful: This was never tested!
typedef float REAL;   // Be careful: This was never tested!
#else
typedef double CTYPE;
typedef double REAL;
#endif

typedef unsigned int UINT;
typedef  int INT;

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
  bool contains (Dune::GeometryType gt)
  {
    if (gt.dim()==dim) return true;
    return false;
  }
};


// enum type to distinguish between forward and adjoint equation:
struct EQ {
  enum Mode{ adjoint=199, forward=911 };
};

// enum type to distinguish between solute and heat transport:
struct Passenger {
  enum Type{ solute=123, heat=321 };
};

struct FEMType {
  enum Type{ DG, CG };
};

#endif // DUNE_GESIS_DATATYPES_HH
