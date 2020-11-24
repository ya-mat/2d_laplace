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

  int i;
  double th;
  int i0;
  int i1;

  for(i = 0; i < n; i++){
    th = 2.0*M_PI*(i + 0.5)/n;
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

    std::cout << "i, u(i), kai(i) " << i << ' ' << u(i) << ' '<< kai(i) << std::endl; //dbg
  }

  for(j = 0; j < n; j++){
     j0 = edge(0, j);
     j1 = edge(1, j);
     for(i = 0; i < n; i++){
       i0 = edge(0, i);
       i1 = edge(1, i);
//        xm(:) = 0.5d0*(x(:,i0) + x(:,i1))
//        lp1(i, j) = dble(slp_laplace(xm, x(:,j0), x(:,j1), hs(j), xn(:,j)))
//        lp2(i, j) = dble(dlp_laplace(xm, x(:,j0), x(:,j1), hs(j), xn(:,j), exterior))
     }
  }
}

//  exterior = 0
// 
//  do j = 1, n
//     j0 = edge(1, j)
//     j1 = edge(2, j)
//     do i = 1, n
//        i0 = edge(1, i)
//        i1 = edge(2, i)
//        xm(:) = 0.5d0*(x(:,i0) + x(:,i1))
//        lp1(i, j) = dble(slp_laplace(xm, x(:,j0), x(:,j1), hs(j), xn(:,j)))
//        lp2(i, j) = dble(dlp_laplace(xm, x(:,j0), x(:,j1), hs(j), xn(:,j), exterior))
//     end do
//  end do
