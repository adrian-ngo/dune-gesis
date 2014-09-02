/* 
 * File:   Vector.hh
 * Author: ngo
 *
 * Created on July 15, 2010, 2:32 PM
 */

#ifndef _VECTOR_HH
#define	_VECTOR_HH
#include <vector>
#include <cmath>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <iomanip>

template<typename ComponentType>
class Vector : public std::vector<ComponentType> {

  typedef std::vector<ComponentType> basetype;

public:

  Vector() : std::vector<ComponentType>() {
  }

  Vector( const std::vector<ComponentType>&x ) 
    : std::vector<ComponentType>( x ) {
  }

  Vector( const Vector& orig ) : basetype(orig) {
  }

  Vector(const int size, const ComponentType defaultvalue_ = 0)
    : std::vector<ComponentType>( size, defaultvalue_ ) {
  }

  /*
  Vector(const int dim, const ComponentType* array_){
    for(size_t i=0; i<dim; i++){
      ComponentType a = array_[i];
      this->push_back(a);
    }
  }
  */

  void reset( const int dim ){
    Vector &self = *this;
    for( int i=0;i<dim;++i )
      self[i] = 0;
  }

  // inner product with another vector:
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




#endif	/* _VECTOR_HH */

