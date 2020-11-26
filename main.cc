// -*- coding: utf-8 -*-

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <chrono>
#include <assert.h>

double slp_laplace(const Eigen::Vector2d &xm, const Eigen::Vector2d &ym, const Eigen::Vector2d &y1,
		   const Eigen::Vector2d &y2, const double &hy, const Eigen::Vector2d &yn){
//  Eigen::Vector2d ym;
  Eigen::Vector2d rr1;
  Eigen::Vector2d rr2;
  double r1;
  double r2;
  Eigen::Vector2d yt;
  double uk1, uk2, uk3, uk4, uk5;
  constexpr double arctwopi = 1.0/(2.0*std::acos(-1.0));

//  ym = 0.5*(y1 + y2);

  if((xm - ym).squaredNorm() < 1.0e-8){
    return hy*(1.0 - std::log(hy*0.5))*arctwopi;
  }
  else{
    rr1 = xm - y1;
    rr2 = xm - y2;
    r1 = rr1.norm();
    r2 = rr2.norm();
    yt(0) = -yn(1);
    yt(1) = yn(0);

    uk1 = rr1.dot(yn);
    uk2 = rr1.dot(yt);
    uk3 = rr2.dot(yn);
    uk4 = rr2.dot(yt);
    uk5 = std::atan2(uk3, uk4) - std::atan2(uk1, uk2);
    return (uk4*std::log(r2) - uk2*std::log(r1) + hy - uk1*uk5)*arctwopi;
//    std::cout << "tes3 " << std::endl; //dbg
  }
}

double dlp_laplace(const Eigen::Vector2d &xm, const Eigen::Vector2d &ym, const Eigen::Vector2d &y1,
		   const Eigen::Vector2d &y2, const double &hy, const Eigen::Vector2d &yn,
		   const int &exterior){

//  Eigen::Vector2d ym;
  Eigen::Vector2d rr1;
  Eigen::Vector2d rr2;
  double r1;
  double r2;
  Eigen::Vector2d yt;
  double uk1, uk2, uk3, uk4;
  constexpr double arctwopi = 1.0/(2.0*std::acos(-1.0));

  //  ym = 0.5*(y1 + y2);

  if((xm - ym).squaredNorm() < 1.0e-8){
    if(exterior == 0){
      return 0.5;
    }else if(exterior == 1){
      return -0.5;
    }else{
      assert(false);
    }
  }
  else{
    rr1 = xm - y1;
    rr2 = xm - y2;
    r1 = rr1.norm();
    r2 = rr2.norm();
    yt(0) = -yn(1);
    yt(1) = yn(0);

    uk1 = rr1.dot(yn);
    uk2 = rr1.dot(yt);
    uk3 = rr2.dot(yn);
    uk4 = rr2.dot(yt);
    return arctwopi*(std::atan2(uk3, uk4) - std::atan2(uk1, uk2));
  }
}

int main(){
  int n;
  double rad;
  std::ifstream fr("input");
  fr >> n >> rad;
  std::cout << "n " << n << std::endl;
  std::cout << "rad " << rad << std::endl;

  Eigen::MatrixXd x = Eigen::MatrixXd::Zero(2, n);
  Eigen::MatrixXd xn = Eigen::MatrixXd::Zero(2, n);
  Eigen::MatrixXi edge = Eigen::MatrixXi::Zero(2, n);
  Eigen::VectorXd hs = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd slp = Eigen::MatrixXd::Zero(n, n);
  Eigen::MatrixXd dlp = Eigen::MatrixXd::Zero(n, n);
  Eigen::VectorXd u = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd kai = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd xm = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd ym = Eigen::VectorXd::Zero(2);

  constexpr double pi = std::acos(-1.0);

  int i;
  int j;
  double th;
  int i0;
  int i1;
  int j0;
  int j1;
  int exterior = 0;

  for(i = 0; i < n; i++){
    th = 2.0*pi*(i + 1.5)/n;
    x(0, i) = rad*cos(th);
    x(1, i) = rad*sin(th);
    edge(0, i) = i;
    edge(1, i) = (i + 1) % n;
  }

  Eigen::Vector2d rr1; //dbg
  Eigen::Vector2d rr2; //dbg
  double r1; //dbg
  double r2; //dbg
  Eigen::Vector2d yt; //dbg
  double uk1, uk2, uk3, uk4, uk5; //dbg
  constexpr double arctwopi = 1.0/(2.0*std::acos(-1.0)); //dbg

  for(i = 0; i < n; i++){
    i0 = edge(0, i);
    i1 = edge(1, i);
    xm = 0.5*(x.col(i0) + x.col(i1));

    if((i0 >= 0 && i0 < n) && (i1 >= 0 && i1 < n)){
      xn(0, i) = x(1, i1) - x(1, i0);
      xn(1, i) = x(0, i0) - x(0, i1);

      hs(i) = std::sqrt(std::pow(xn(0, i), 2.0) + std::pow(xn(1, i), 2.0));

      xn(0, i) = xn(0, i)/hs(i);
      xn(1, i) = xn(1, i)/hs(i);
    }
    else{
      assert(false);
    }
    // x**3*y - x*y**3
    u(i) = std::pow(xm(0), 3.0)*xm(1) - xm(0)*std::pow(xm(1), 3.0);

    // (3x**2*y - y**3)nx + (x**3 - 3xy**2)ny
    kai(i) = (3.0*std::pow(xm(0), 2.0)*xm(1) - std::pow(xm(1), 3.0))*xn(0, i)
      + (std::pow(xm(0), 3.0) - 3.0*xm(0)*std::pow(xm(1), 2.0))*xn(1, i);
  }

  for(j = 0; j < n; j++){
     j0 = edge(0, j);
     j1 = edge(1, j);
     ym = 0.5*(x.col(j0) + x.col(j1));
     for(i = 0; i < n; i++){
       i0 = edge(0, i);
       i1 = edge(1, i);
       xm = 0.5*(x.col(i0) + x.col(i1));

       slp(i, j) = slp_laplace(xm, ym, x.col(j0), x.col(j1), hs(j), xn.col(j));
       dlp(i, j) = dlp_laplace(xm, ym, x.col(j0), x.col(j1), hs(j), xn.col(j), exterior);
     }
  }

  u = slp.partialPivLu().solve(dlp*u);

  std::cout << "relative error " << std::sqrt(((u - kai).dot(u - kai))/kai.dot(kai)) << std::endl;
}
