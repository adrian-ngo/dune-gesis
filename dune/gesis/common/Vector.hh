/* 
 * File:   Vector.hh
 * Author: A. Ngo
 *
 * 2010-2014
 */

#ifndef DUNE_GESIS_VECTOR_HH
#define	DUNE_GESIS_VECTOR_HH
#include <vector>
#include <cmath>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <iomanip>


namespace Dune {
  namespace Gesis {

    //! Extension of the STL vector class implementing a finite-dimensional vector.
    //! \tparam ComponentType Number type of the vector components.
    template<typename ComponentType>
    class Vector : public std::vector<ComponentType> {

      typedef std::vector<ComponentType> basetype;

    public:

      //! \brief Default Constructor
      Vector() : std::vector<ComponentType>() {
      }

      //! \brief Constructor initializing Vector from an std::vector.
      Vector( const std::vector<ComponentType>&x ) 
        : std::vector<ComponentType>( x ) {
      }

      //! \brief Constructor initializing Vector from another Vector.
      Vector( const Vector& orig ) : basetype(orig) {
      }

      //! \brief Constructor
      //! \param[in] size The dimension of the vector.
      //! \param[in] defaultvalue All components of the vector are set to defaultvalue.
      Vector( const int size, const ComponentType defaultvalue = 0 )
        : std::vector<ComponentType>( size, defaultvalue ) {
      }

      //! \brief sets all components of the vector to zero.
      //! \param[in] dim Dimension of the vector
      void reset( const int dim ){
        Vector &self = *this;
        for( int i=0;i<dim;++i )
          self[i] = 0;
      }

      //! \brief inner product with another vector:
      //! \param[in] x The other vector.
      ComponentType operator*(Vector & x) {
        assert( x.size() == this->size() );
        ComponentType sum( 0 );
        Vector & self = *this;
        for( size_t i = 0; i < this->size(); ++i )
          sum += self[i] * x[i];
        return sum;
      }

      // get the product of all the components
      ComponentType componentsproduct() const {
        ComponentType product( 1 );
        const Vector & self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          product *= self[i];
        return product;
      }

      // square of the euclidean norm:
      ComponentType two_norm_sqr() const {
        ComponentType sum( 0 );
        const Vector & self = *this;
        for (int i = 0; i < (int) this->size(); ++i)
          sum += self[i] * self[i];
        return sum;
      }

      // euclidean norm:
      ComponentType two_norm() const {
        return sqrt(two_norm_sqr());
      }

      // Manhattan norm:
      ComponentType one_norm() const {
        ComponentType sum( 0 );
        const Vector & self = *this;
        for (int i = 0; i < (int) this->size(); ++i)
          sum += std::abs( self[i] );
        return sum;
      }
  
      // Maximum norm:
      ComponentType maximum_norm() const {
        ComponentType maximum( 0 );
        const Vector & self = *this;
        for (int i = 0; i < (int) this->size(); ++i)
          maximum = std::max( std::abs( self[i] ), maximum );
        return maximum;
      }

      // Arithmetic mean:
      ComponentType mean() const {
        ComponentType sum( 0 );
        const Vector & self = *this;
        int N = (int) this->size();
        for (int i = 0; i<N; ++i)
          sum += self[i];

        return REAL(sum) / REAL(N);
      }

      // Variance:
      void mean_and_variance( REAL& mean_value, REAL& variance) const {
        ComponentType sum( 0 );
        const Vector & self = *this;
        mean_value = mean();
        int N = (int) this->size();
        for (int i = 0; i<N; ++i){
          ComponentType summand = self[i] - mean_value;
          sum += summand * summand;
        }
        variance = REAL(sum) / REAL(N-1);
        return;
      }


      // x[i] = value for all components
      Vector &operator=(const ComponentType value) {
        Vector & self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          self[i] = value;
        return *this;
      }

      // x[i] *= value
      Vector &operator*=(const ComponentType value) {
        Vector& self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          self[i] *= value;
        return *this;
      }

      // x[i] /= value
      Vector &operator/=(const ComponentType value) {
        Vector& self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          self[i] /= value;
        return *this;
      }

      // x[i] += value
      Vector &operator+=(const ComponentType value) {
        Vector& self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          self[i] += value;
        return *this;
      }

      // x += y
      Vector &operator+=(const Vector& b) {
        Vector& self = *this;
        assert( self.size() == b.size() );
        for (size_t i = 0; i < this->size(); ++i)
          self[i] += b[i];
        return *this;
      }

      // x -= y
      Vector &operator-=(const Vector& b) {
        Vector& self = *this;
        assert( self.size() == b.size() );
        for (size_t i = 0; i < this->size(); ++i)
          self[i] -= b[i];
        return *this;
      }

      // x[i] -= value
      Vector &operator-=(const ComponentType value) {
        Vector& self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          self[i] -= value;
        return *this;
      }

      Vector operator-(const Vector& b){    
        const Vector &a = *this;
        assert( a.size() == b.size() );
        Vector c(a);
        assert( c.size() == b.size() );
        c.axpy(-1,b);
        return c;
      };

      // x = x + alpha * y :
      Vector & axpy(const ComponentType alpha, const Vector & y) {
        assert( this->size() == y.size());

        Vector & self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          self[i] += alpha * y[i];
        return *this;
      }

      // x = x + alpha * y
      Vector & scaled_add(const ComponentType alpha, const Vector & y) {
        assert( this->size() == y.size());

        Vector & self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          self[i] += alpha * y[i];
        return *this;
      }

      // x = alpha * x + y

      Vector & a_times_x_plus_y(const ComponentType alpha, const Vector & y) {
        assert( this->size() == y.size());

        Vector & self = *this;
        for (size_t i = 0; i < this->size(); ++i)
          self[i] = alpha * self[i] + y[i];
        return *this;
      }

      void print(void) const {
        const Vector &self = *this;
        std::cout << "( ";
        for (size_t i = 0; i < self.size(); i++) {
          std::cout << self[i] << " ";
        }
        std::cout << " )" << std::endl;
      }

    };


    template<typename ComponentType>
    std::ostream & operator <<(std::ostream & os, const Vector<ComponentType> & x) {
      os << "( ";
      for (int r = 0; r < (int) x.size(); ++r) {
        os << x[r] << " ";
      }
      os << ")" << std::endl;

      return os;
    }


  } // Gesis
} // Dune

#endif	/* DUNE_GESIS_VECTOR_HH */

