#ifndef tws_cg_hpp
#define tws_cg_hpp

#include <cassert>
#include <type_traits>
#include "vector.hpp"

namespace tws {
  //(C) Copyright Karl Meerbergen, Gowri Suryanarayana, Yuya Suzuki & Joris Tavernier, 2016.
  //
  // This is code for the conjugate gradient Krylov method for solving a symmetric linear system (that is guaranteed to work for positive matrices)
  // 
  // Stopcriterion:
  //
  //   number of iterations <= maximum number of iterations max_it
  //   ||r|| <= tolerance
  //
  // Input
  //   op:
  //     binary functor:
  //         op( x, y ) computes y = A * x, i.e., op is a linear operator
  //     where x and y are Vectors.
  //
  //   r takes the right-hand side.
  //   is_vector<R> is equivalent to std::true_type
  //
  //   x takes the initial solution.
  //   is_vector<X> is equivalent to std::true_type
  //
  //   tolerance: is the residual tolerance of the linear solver, see further.
  //   maxit: is the maximum number of iterations.
  //
  // Output:
  //   r is the residual of the final solution: r = b - A * x
  //   x is the solution.
  //

  template <typename X, typename R, typename Op>
  void cg( Op const& op, X& x, R& r, double const& tolerance, int max_it ) {
    typedef typename X::value_type value_type ;
    typedef vector<typename X::value_type> container_type ;

    assert( x.size() == r.size() ) ;
    assert(tolerance>0.);
    assert(max_it>0);

    container_type  p( x.size() ) ;
    container_type  q( x.size() ) ;
    container_type  x_initial( x.size(),0. ) ;
    // Compute residual

    auto res_norm_0 = norm_2( r ) ;
    op( x, q ) ;
    r -= q ;
    value_type rho2 ;
    value_type rho1 ;

    for ( int it=0 ; it<max_it; ++it ) {
      if ( norm_2( r ) < tolerance*res_norm_0 ) break ;
      rho1 = inner_product( r, r ) ;

      if (it==0) {
        p = r ;
      } else {
        p = r + (rho1/rho2) * p ;
      }
      op( p, q ) ;
      value_type alpha = rho1 / inner_product( p, q ) ;
      x_initial += alpha * p ;
      r -= alpha * q ;
      rho2 = rho1 ;
    }
    x+=x_initial;
  } // cg
}

#endif
