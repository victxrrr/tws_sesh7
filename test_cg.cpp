#include "vector.hpp"
#include "cg.hpp"
#include <iostream>
#include <typeinfo>
#include <type_traits>
#include <algorithm>
void matvec ( tws::vector<double> const& x, tws::vector<double>& y ) {
    assert( x.size()==y.size() ) ;

    for (decltype(x.size()) i=0; i<x.size(); ++i) {
      y(i) = x(i) / (i+1) ;
    }
}

template <typename scalar>
void matvec1( tws::vector<scalar> const& x, tws::vector<scalar>& y ) {
    assert( x.size()==y.size() ) ;

    for (decltype(x.size()) i=0; i<x.size(); ++i) {
      y(i) = x(i) / (i+1) ;
    }
}
	
struct matvec2{
  void operator()(tws::vector<float> const& x, tws::vector<float>& y) const{
    assert( x.size()==y.size() ) ;

    for (decltype(x.size()) i=0; i<x.size(); ++i) {
      y(i) = x(i) / (i+1) ;
    }
  }
};

struct matvec3{
  matvec3(float m)
  :m_(m)
  {}
  void operator()(tws::vector<float> const& x, tws::vector<float>& y) const{
    assert( x.size()==y.size() ) ;

    for (decltype(x.size()) i=0; i<x.size(); ++i) {
      y(i) = x(i) / ((i+1) + m_) ;
    }
  }
  float m_;
};

struct matvec4{
  void operator()(tws::vector<float> const& x, tws::vector<float>& y, float m) const{
    assert( x.size()==y.size() ) ;

    for (decltype(x.size()) i=0; i<x.size(); ++i) {
      y(i) = x(i) / ((i+1) + m) ;
    }
  }
};

int main() {
  int n=100;
  tws::vector<double> b(n) ;
  tws::vector<double> sol(n) ;
  tws::vector<double> x(n) ;
  tws::vector<double> b_ex(n) ;

  //x random between 0 and 1
  x.randomize(0,1);

  matvec( x, b ) ;

  b_ex=b ;

  //x zero vector
  std::fill(x.begin(),x.end(),0.);
  tws::cg( matvec, x, b, 1.e-10, n ) ;
  matvec ( x, sol ) ;

  std::cout<<"relative error: "<<tws::norm_2(sol-b_ex)/tws::norm_2(b_ex)<<std::endl;

  tws::vector<float> b_f(n), sol_f(n), x_f(n), b_ex_f(n);
  x_f.randomize(0,1);
  matvec1( x_f, b_f );

  b_ex_f = b_f;
  std::fill(x_f.begin(),x_f.end(),0.f);
  tws::cg(matvec1<float>, x_f, b_f, 1.e-10, n);
  matvec1(x_f, sol_f);
  std::cout<<"relative error: "<<tws::norm_2(sol_f-b_ex_f)/tws::norm_2(b_ex_f)<<std::endl;
  
  x_f.randomize(0,1);
  matvec2 mv;
  mv(x_f, b_f);
  b_ex_f = b_f;
  std::fill(x_f.begin(),x_f.end(),0.f);
  tws::cg(mv, x_f, b_f, 1.e-10, n);
  mv(x_f, sol_f);
  std::cout<<"relative error: "<<tws::norm_2(sol_f-b_ex_f)/tws::norm_2(b_ex_f)<<std::endl;
  
  x_f.randomize(0,1);
  float m = 3;
  matvec3 mv3(m);
  mv3(x_f, b_f);
  b_ex_f = b_f;
  std::fill(x_f.begin(),x_f.end(),0.f);
  tws::cg(mv3, x_f, b_f, 1.e-10, n);
  mv3(x_f, sol_f);
  // for (decltype(x.size()) i=0; i<x.size(); ++i) {
  //     b_ex_f(i) = b_ex_f(i) / (i+1) * ((i+1) + m) ;
  //     sol_f(i) = sol_f(i) / (i+1) * ((i+1) + m) ;
  // }
  std::cout<<"relative error: "<<tws::norm_2(sol_f-b_ex_f)/tws::norm_2(b_ex_f)<<std::endl;
  
  x_f.randomize(0,1);
  matvec4 mv4;
  auto apply = [&mv4, &m] (auto const& x, auto & y) -> void {mv4(x,y,m);};
  apply(x_f, b_f);
  b_ex_f = b_f;
  std::fill(x_f.begin(),x_f.end(),0.f);
  tws::cg(apply, x_f, b_f, 1.e-10, n);
  apply(x_f, sol_f);
  std::cout<<"relative error: "<<tws::norm_2(sol_f-b_ex_f)/tws::norm_2(b_ex_f)<<std::endl;

  return 0 ;
} 
