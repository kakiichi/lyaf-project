! Usage of linear interpolation subroutine
!
!  Input:  double precision array of size n: xx(n), yy(n)
!          double precision: x, point of evaluation
!  Output: double precision: y, interpolated value at x
! 
!  xx(1:n) is an array monotically increasing or decreasing
!  such as grid xx=1,2,3,4,....,10, with n grid point
!  yy(1:n) is an array. Discretely sampled function y=y(x)
!  at the point of xx(1:n) array
!
!  Example:
!
!  integer :: j  <- internal value to indicate the array index
!
!  call locate(xx,n,x,j)
!  call linint(xx(j),xx(j+1),yy(j),yy(j+1),x,y)
!
!  Now we have the value of y at x !!



!==== Numerical Recipes code ==================================
!==== Linear interpolation   ==================================
!  
!  Suppose y=y(x), and we want the value y at x given 
!  array y_i=y(x_i). x1 is lower value x2 is upper value
!  near to x i.e. x1<x<x2. y1 and y2 are the known values of
!  array at x1 and x2. linint returns the value y at x. 

!  linint : input, double precision
!             x1  lower value of x
!             x2  upper value of x
!             y1  array value at x1
!             y2  array value at x2
!             x   value we want at which y(x) is evaluated
!          output, double precision
!             y   interpolated value at x
!
               

subroutine linint(x1,x2,y1,y2,x,y)
  implicit none
  double precision :: x1,x2,y1,y2,x,y
  double precision :: A,B
  
  A=(x2-x)/(x2-x1)
  B=1.d0-A
  y=A*y1+B*y2

  return
end subroutine linint

!==== Bisection search =======================================
!  'locate' subroutine from Numerical recipe
!  Given an array xx(1:n) and given a value x, returns a value
!  j such that x is between xx(j) and xx(j+1). x(1:n) must be
!  monotonic, either increasing or decreasing. j=0 or j=n
!  is returned to indicate that x is out of range.
 
subroutine locate(xx,n,x,j)

  integer :: j,n
  double precision :: x, xx(n)
  integer :: jl,jm,ju

  jl=0
  ju=n+1
10 if (ju-jl .gt. 1) then
     jm=(ju+jl)/2
     if((xx(n) .ge. xx(1)) .eqv. (x .ge. xx(jm))) then
        jl=jm
     else
        ju=jm
     endif
     goto 10
  endif
  ! set the output
  if (x .eq. xx(1)) then
     j=1
  else if (x .eq. xx(n)) then
     j=n-1
  else
     j=jl
  endif

  ! error call
  if ((j .eq. n) .or. (j .eq. 0)) then
     stop "Error: Interpolation beyond the known range."
  endif
  
  return
end subroutine locate
