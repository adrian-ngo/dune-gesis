// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_GEOINVERSION_DARCYVELOCITYCACHE_HH
#define DUNE_GEOINVERSION_DARCYVELOCITYCACHE_HH

#include<vector>
#include<map>

#include<dune/common/exceptions.hh>
#include<dune/common/static_assert.hh>

#include "../DG/RT0flux.hh"

namespace Dune {
  namespace GeoInversion {

    //! \brief store values of basis functions and gradients in a cache
    template<typename GWP,typename GFS,typename DarcyFluxType>
    class DarcyVelocityCache
    {
      typedef typename DarcyFluxType::Element Element;
      typedef typename DarcyFluxType::VCType VCType;
      typedef typename DarcyFluxType::DF DomainFieldType;
      enum {dim = DarcyFluxType::dim};
      typedef typename DarcyFluxType::Domain DomainType;
      typedef typename DarcyFluxType::VectorRange RangeType;

      struct less_than
      {
        bool operator() (const DomainType& v1, const DomainType& v2) const
        {
          for (typename DomainType::size_type i=0; i<DomainType::dimension; i++)
            {
              if ( v1[i] < v2[i]-1e-5 ) return true;   // is less than
              if ( v1[i] > v2[i]+1e-5 ) return false;  // is greater than
            }
          return false; // is equal
        }
      };

      typedef std::map<DomainType,RangeType,less_than> FluxCache;

    public:

      //! \brief constructor
      DarcyVelocityCache( const GWP& gwp_,
                          const GFS& gfs_,
                          const VCType& vchead_,
                          const int baselevel_=0,
                          const bool withoutK_=false ) :
        dgf( gwp_,
             gfs_,
             vchead_,
             baselevel_,
             withoutK_ )
      {
        if( !fluxcache.empty() )
          fluxcache.clear();
      }
      
      const DarcyFluxType& exportDGF(){
        return dgf;
      }

      //! evaluate flux at a point on an element
      void evaluate( const Element& e,
                     const DomainType& localposition,
                     RangeType& fluxVector ) const
      {
        evaluate_on_root( e, localposition, fluxVector );
        return;
      }

      void evaluate_on_root( const Element& e,
                             const DomainType& localposition,
                             RangeType& fluxVector ) const
      {
        DomainType globalposition = 
          e.geometry().global( localposition );

        // Retrieve fluxVector from cache, if available:
        typename FluxCache::iterator it = fluxcache.find( globalposition );
        if( it != fluxcache.end() ){
          fluxVector = it->second;
          //std::cout << "huhu" << std::endl;
          return;
        }
        
        // Otherwise, compute the fluxVector and add it to the fluxcache:
        dgf.evaluate_on_root( e, localposition, fluxVector );
        it = fluxcache.insert( fluxcache.begin(), std::pair<DomainType,RangeType>( globalposition, fluxVector ) );
        //std::cout << "haha haha size = " << fluxcache.size() << std::endl;
        return;
      }

      void emptyCache(){
        fluxcache.clear();
      }

    private:
      DarcyFluxType dgf;
      mutable FluxCache fluxcache;
    };

  }
}

#endif // DUNE_GEOINVERSION_DARCYVELOCITYCACHE_HH

