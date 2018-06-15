module mfunciones
public :: a, indices, inverse, error_inverse


contains

function a (teta, phi)
! Vector unitario para la representacion espacial

real, intent(in)    :: teta, phi
real, dimension (3) :: a

a (1) = sin(teta) * cos(phi)
a (2) = sin(teta) * sin(phi)
a (3) = cos(teta)

endfunction a


! ================================================================



subroutine indices (i, j, m, factor)
! Coeficientes para el cambio de notacion reducida a tensorial en el tensor S_{ij}

integer, intent(in) :: i, j
integer, intent(out) :: m
real, intent(out)   :: factor

if (i == j) then
factor = 1.0
else
factor = 0.5
endif

if (i == 1) then
  if (j == 1) then
    m = 1
  endif
  if (j == 2) then
    m = 6
  endif
  if (j == 3) then
    m = 5
  endif
endif

if (i == 2) then
  if (j == 1) then
    m = 6
  endif
  if (j == 2) then
    m = 2
  endif
  if (j == 3) then
    m = 4
  endif
endif

if (i == 3) then
  if (j == 1) then
    m = 5
  endif
  if (j == 2) then
    m = 4
  endif
  if (j == 3) then
    m = 3
  endif
endif

endsubroutine indices



! ================================================================



subroutine inverse(a, n, c)
! Inversion de una matriz mediante el descomposicion LU.

implicit none

integer, intent (in) :: n
integer :: i, j, k
real, dimension (:,:), intent(in)   :: a
real, dimension (n, n), intent(out) :: c
real, dimension (n, n)              :: L, U, aa
real, dimension (n)                 :: b, d, x
real                                :: coeff

do i = 1, n
  do j = 1, n
    aa(i, j) = a(i, j)
  enddo
enddo

L=0.0
U=0.0
b=0.0

do k=1, n-1
   do i=k+1,n
      coeff=aa(i,k)/aa(k,k)
      L(i,k) = coeff
      do j=k+1,n
         aa(i,j) = aa(i,j)-coeff*aa(k,j)
      end do
   end do
end do

do i=1,n
  L(i,i) = 1.0
end do
do j=1,n
  do i=1,j
    U(i,j) = aa(i,j)
  end do
end do

do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse



! ================================================================



subroutine error_inverse (eps, M_inv, eps_inv)
! Propagacion de errores en la inversion matricial.

real, dimension (:,:), intent (in)  :: eps, M_inv
real, dimension (6,6), intent (out) :: eps_inv

integer                             :: i, j, k, l
real                                :: aux

do k = 1, 6
  do l = 1, 6
    aux = 0.0
    do i = 1, 6
      do j = 1, 6
        aux = aux + ( (M_inv(k,i)) * (eps(i,j)) * (M_inv(j,l)) )**2
      enddo
    enddo
    eps_inv(k,l) = sqrt(aux)
  enddo
enddo

endsubroutine error_inverse




endmodule mfunciones
