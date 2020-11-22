// -*- coding: utf-8 -*-

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <chrono>


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
    test(edge(1, i));
    test2(edge);
  }
}


//  do i = 1, n
//     i0 = edge(1, i)
//     i1 = edge(2, i)
//     if((i0 .ge. 1 .and. i0 .lt. n+1) .and. (i1 .ge. 1 .and. i1 .lt. n+1)) then
//        xn(1, i) = x(2, i1) - x(2, i0)
//        xn(2, i) = x(1, i0) - x(1, i1)
// 
//        hs(i) = sqrt(xn(1, i)**2 + xn(2, i)**2)
// 
//        xn(1, i) = xn(1, i)/hs(i)
//        xn(2, i) = xn(2, i)/hs(i)
//     else
//        call force_raise()
//     end if
// 
//     !x**3*y - x*y**3
//     u(i) = x(1, i)**3*x(2, i) - x(1, i)*x(2, i)**3
// 
//     !(3x**2*y - y**3)nx + (x**3 - 3xy**2)ny
//     kai(i) = (3*x(1, i)**2*x(2, i) - x(2, i)**3)*xn(1, i) + (x(1, i)**3 - 3*x(1, i)*x(2, i)**2)*xn(2, i)
//  end do
