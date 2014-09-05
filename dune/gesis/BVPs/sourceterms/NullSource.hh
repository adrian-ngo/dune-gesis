#ifndef DUNE_GESIS_NULL_SOURCE_HH
#define DUNE_GESIS_NULL_SOURCE_HH

#include "SourceBaseClass.hh"

namespace Dune {

  namespace Gesis {

    /**
     * The NullSource is supposed to be a dummy class adding 
     * no contribution to source term.
     */

    template<typename GV,typename RF>
    class NullSource 
      : public SourceTermInterface<
      SourceTermTraits<GV,RF>,
      NullSource<GV,RF>
      >
    {
    private:
      enum{dim=GV::dimension};

    public:
      const SourceNature source_nature;
      typedef SourceTermTraits<GV,RF> Traits;
      
      NullSource()
        :
        source_nature( NULL_SOURCE )
      {}
      
    }; // class NullSource
  }
}

#endif
