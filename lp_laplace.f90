!-------------------------------------------
module lp_lap
  implicit none
contains
!----------------------------------------------
  complex(c_double_complex) function slp_laplace(xm,y1,y2,hy,n) bind(c)
    use iso_c_binding
    implicit none

    real(c_double),intent(in) :: xm(2),y1(2),y2(2),hy,n(2)

    real(c_double) :: t(2), r1, r2
    real(c_double) :: rr1(2), rr2(2)
    real(c_double) :: uk1,uk2,uk3,uk4,uk5
    real(c_double) :: dnrm2 !blas
    real(c_double),parameter :: arctwopi = 0.159154943091895d0

    real(c_double) :: ym(2)

    ym = (y1 + y2)*0.5d0

    if(abs(xm(1)-ym(1))+abs(xm(2)-ym(2)) .lt. 1d-8) then
       slp_laplace = hy*(1d0-log(hy*0.5d0))*arctwopi
    else
       rr1 = xm - y1
       rr2 = xm - y2
       r1 = dnrm2(2, rr1(1), 1)
       r2 = dnrm2(2, rr2(1), 1)
       t(1) = -n(2)
       t(2) = n(1)

       uk1 = dot_product(rr1,n)
       uk2 = dot_product(rr1,t)
       uk3 = dot_product(rr2,n)
       uk4 = dot_product(rr2,t)
       uk5 = atan2(uk3,uk4)-atan2(uk1,uk2)
       slp_laplace = dcmplx((uk4*log(r2)-uk2*log(r1)+hy-uk1*uk5)*arctwopi, 0d0)
    end if

  end function slp_laplace
!----------------------------------------
  complex(c_double_complex) function dlp_laplace(xm,y1,y2,hy,n,exterior) bind(c)
    use iso_c_binding
    implicit none

    real(c_double),intent(in) :: xm(2),y1(2),y2(2),hy,n(2)
    integer(c_int),optional,intent(in) :: exterior

    real(c_double) :: t(2),r1,r2
    real(c_double) :: rr1(2), rr2(2)
    real(c_double) :: uk1,uk2,uk3,uk4,uk5
    real(c_double) :: dnrm2 !blas
    real(c_double),parameter :: arctwopi = 0.159154943091895d0

    real(c_double) :: ym(2)

    ym = (y1 + y2)*0.5d0

    if(abs(xm(1)-ym(1))+abs(xm(2)-ym(2)) .lt. 1d-8) then
       if(present(exterior)) then
          if(exterior == 0) dlp_laplace = dcmplx(0.5d0, 0d0)
          if(exterior == 1) dlp_laplace = dcmplx(-0.5d0, 0d0)
       else
          stop 'if(present(exterior)) is false in dlp_laplace'
       end if
    else
       rr1 = xm - y1
       rr2 = xm - y2
       r1 = dnrm2(2, rr1(1), 1)
       r2 = dnrm2(2, rr2(1), 1)
       t(1) = -n(2)
       t(2) = n(1)

       uk1 = dot_product(rr1, n)
       uk2 = dot_product(rr1, t)
       uk3 = dot_product(rr2, n)
       uk4 = dot_product(rr2, t)
       uk5 = atan2(uk3,uk4)-atan2(uk1,uk2)
       dlp_laplace = dcmplx(uk5*arctwopi, 0d0)
    end if

    return
  end function dlp_laplace
!----------------------------------------
end module lp_lap
