! -*- coding: utf-8 -*-
program main
  use lp_lap
  use mod_force_raise
  implicit none

  real(8) :: rad
  real(8),allocatable :: x(:, :)
  real(8),allocatable :: xn(:, :)
  integer,allocatable :: edge(:, :)
  real(8),allocatable :: hs(:)
  real(8),allocatable :: lp1(:, :)
  real(8),allocatable :: lp2(:, :)
  real(8),allocatable :: u(:)
  real(8),allocatable :: kai(:)
  integer,allocatable :: ipiv(:)
  integer :: n
  integer :: info
  real(8),parameter :: pi = 3.1415926535897932384626433832795028841971d0
  integer :: exterior

  real(8) :: xm(2)
  real(8) :: th
  integer :: i
  integer :: j
  integer :: i0
  integer :: i1
  integer :: j0
  integer :: j1

  write(6,*) '# n, rad'
  open(1000, file='input', status='old')
  read(1000,*) n, rad
  close(1000)

  allocate(x(2, n))
  allocate(xn(2, n))
  allocate(edge(2, n))
  allocate(hs(n))
  allocate(lp1(n, n))
  allocate(lp2(n, n))
  allocate(u(n))
  allocate(kai(n))
  allocate(ipiv(n))

  do i = 1, n
     th = 2d0*pi*(dble(i)+0.5d0)/dble(n)

     x(1, i) = rad*cos(th)
     x(2, i) = rad*sin(th)
     edge(1, i) = i
     edge(2, i) = mod(i, n) + 1
  end do

  do i = 1, n
     i0 = edge(1, i)
     i1 = edge(2, i)
     if((i0 .ge. 1 .and. i0 .lt. n+1) .and. (i1 .ge. 1 .and. i1 .lt. n+1)) then
        xn(1, i) = x(2, i1) - x(2, i0)
        xn(2, i) = x(1, i0) - x(1, i1)

        hs(i) = sqrt(xn(1, i)**2 + xn(2, i)**2)

        xn(1, i) = xn(1, i)/hs(i)
        xn(2, i) = xn(2, i)/hs(i)
     else
        call force_raise()
     end if

     !x**3*y - x*y**3
     u(i) = x(1, i)**3*x(2, i) - x(1, i)*x(2, i)**3

     !(3x**2*y - y**3)nx + (x**3 - 3xy**2)ny
     kai(i) = (3*x(1, i)**2*x(2, i) - x(2, i)**3)*xn(1, i) + (x(1, i)**3 - 3*x(1, i)*x(2, i)**2)*xn(2, i)
  end do

  exterior = 0

  do j = 1, n
     j0 = edge(1, j)
     j1 = edge(2, j)
     do i = 1, n
        i0 = edge(1, i)
        i1 = edge(2, i)
        xm(:) = 0.5d0*(x(:,i0) + x(:,i1))
        lp1(i, j) = slp_laplace(xm, x(:,j0), x(:,j1), hs(j), xn(:,j))
        lp2(i, j) = dlp_laplace(xm, x(:,j0), x(:,j1), hs(j), xn(:,j), exterior)
     end do
  end do

  u = matmul(lp2, u)
  deallocate(lp2)

  call DGESV(N, 1, lp1, n, ipiv, u, n, info)
  if(info.ne.0) then
     write(6,*) 'dgesv, info=', info
     call force_raise()
  endif
  deallocate(lp1)
  deallocate(ipiv)

  write(*,*) 'number of dof', n
  write(*,*) 'rad', rad
  write(*,*) 'relative error', sqrt(dot_product(u - kai, u - kai)/dot_product(kai, kai))

  deallocate(u)
  deallocate(kai)
  deallocate(x)
  deallocate(xn)
  deallocate(edge)
  deallocate(hs)

end program main

