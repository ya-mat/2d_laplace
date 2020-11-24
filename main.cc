// -*- coding: utf-8 -*-

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <chrono>
#include <assert.h>


void test(const int &a){
  std::cout << "test a " << a << std::endl;
}

void test2(Eigen::MatrixXi &a){
  std::cout << "test2 a " << a << std::endl;
}


double slp_laplace(const Eigen::Vector2d &xm, const Eigen::Vector2d &y1,
		   const Eigen::Vector2d &y2, const double &hy, const Eigen::Vector2d &yn){
//  std::cout << "xm " << xm << std::endl; //dbg
//  std::cout << "y1 " << y1 << std::endl; //dbg
//  std::cout << "y2 " << y2 << std::endl; //dbg
//  std::cout << "hy " << hy << std::endl; //dbg
//  std::cout << "yn " << yn << std::endl; //dbg

  Eigen::Vector2d ym;
  Eigen::Vector2d rr1;
  Eigen::Vector2d rr2;
  double r1;
  double r2;
  Eigen::Vector2d yt;
  double uk1, uk2, uk3, uk4, uk5;
  constexpr double arctwopi = 1.0/(2.0*std::acos(-1.0));

  ym = 0.5*(y1 + y2);

  if((xm - ym).squaredNorm() < 1.0e-8){
    return hy*(1.0 - std::log(hy*0.5))*arctwopi; //log
  }
  else{
    rr1 = xm - y1;
    rr2 = xm - y2;
    r1 = rr1.norm();
    r2 = rr2.norm();
    yt(1) = -yn(2);
    yt(2) = yn(1);

    uk1 = rr1.dot(yn);
    uk2 = rr1.dot(yt);
    uk3 = rr2.dot(yn);
    uk4 = rr2.dot(yt);
    uk5 = std::atan(uk3/uk4) - std::atan(uk1/uk2);
    return (uk4*std::log(r2) - uk2*std::log(r1) + hy - uk1*uk5)*arctwopi;
  }
}

int main(){
  //std::cout << std::fixed << setprecision(15);

  int n;
  double rad;
  std::ifstream fr("input");
  fr >> n >> rad;
  std::cout << "n " << n << std::endl;
  std::cout << "rad " << rad << std::endl;

  Eigen::MatrixXd x(2, n);
  Eigen::MatrixXd xn(2, n);
  Eigen::MatrixXi edge(2, n);
  Eigen::VectorXd hs(n);
  Eigen::MatrixXd slp(n, n);
  Eigen::MatrixXd dlp(n, n);
  Eigen::VectorXd u(n);
  Eigen::VectorXd kai(n);
  Eigen::Vector2d xm;

  constexpr double pi = std::acos(-1.0);

  int i;
  int j;
  double th;
  int i0;
  int i1;
  int j0;
  int j1;

  for(i = 0; i < n; i++){
    th = 2.0*pi*(i + 0.5)/n;
    x(0, i) = rad*cos(th);
    x(1, i) = rad*sin(th);
    edge(0, i) = i;
    edge(1, i) = (i + 1) % n;

//    std::cout << "i " << i << std::endl;
//    std::cout << "th " << th << std::endl;
//    std::cout << "x(0, i) " << x(0, i) << std::endl;
//    std::cout << "x(1, i) " << x(1, i) << std::endl;
//    std::cout << "edge(0, i) " << edge(0, i) << std::endl;
//    std::cout << "edge(1, i) " << edge(1, i) << std::endl;
  }

  for(i = 0; i < n; i++){
    i0 = edge(0, i);
    i1 = edge(1, i);

    if((i0 >= 0 && i0 < n) && (i1 >= 0 && i1 < n)){
      xn(0, i) = x(1, i1) - x(1, i0);
      xn(1, i) = x(0, i0) - x(0, i1);

      hs(i) = std::sqrt(std::pow(xn(0, i), 2.0) + std::pow(xn(1, i), 2.0));

      xn(0, i) = xn(0, i)/hs(i);
      xn(1, i) = xn(1, i)/hs(i);
//      std::cout << "i, xn(0, i), xn(1, i) " << i << ' ' << xn(0, i) << ' '<< xn(1, i) << std::endl; //dbg
//      std::cout << "hs(i) " << hs(i) << std::endl; //dbg
    }
    else{
      assert(false);
    }
    // x**3*y - x*y**3
    u(i) = std::pow(x(0, i), 3.0)*x(1, i) - x(0, i)*std::pow(x(1, i), 3.0);

    // (3x**2*y - y**3)nx + (x**3 - 3xy**2)ny
    kai(i) = (3.0*std::pow(x(0, i), 2.0)*x(1, i) - std::pow(x(1, i), 3.0))*xn(0, i)
      + (std::pow(x(0, i), 3.0) - 3.0*x(0, i)*std::pow(x(1, i), 2.0))*xn(1, i);

    //std::cout << "i, u(i), kai(i) " << i << ' ' << u(i) << ' '<< kai(i) << std::endl; //dbg
  }

  for(j = 0; j < n; j++){
     j0 = edge(0, j);
     j1 = edge(1, j);
     for(i = 0; i < n; i++){
       i0 = edge(0, i);
       i1 = edge(1, i);
       xm = 0.5*(x.col(i0) + x.col(i1));

//       std::cout << "xm " << xm << std::endl; //dbg
//       std::cout << "x.col(j0) " << x.col(j0) << std::endl; //dbg
//       std::cout << "x.col(j1) " << x.col(j1) << std::endl; //dbg
//       std::cout << "hs(j) " << hs(j) << std::endl; //dbg
//       std::cout << "xn.col(j) " << xn.col(j) << std::endl; //dbg

       slp(i, j) = slp_laplace(xm, x.col(j0), x.col(j1), hs(j), xn.col(j));
       std::exit(1);
//        lp2(i, j) = dble(dlp_laplace(xm, x(:,j0), x(:,j1), hs(j), xn(:,j), exterior))
     }
  }
}
