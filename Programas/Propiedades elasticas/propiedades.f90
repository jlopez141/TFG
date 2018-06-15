module m_young_compresibilidad
use mfunciones

real, parameter, public :: pi = acos(-1.0)
public :: young, comp_lineal, average_voight, average_reuss



contains

function young (S, teta, phi, a, indices)
! Valor del modulo de Young en una direccion (teta,phi) dada.

real, intent(in)                 :: teta, phi
real, intent(in), dimension(:,:) :: S

interface
  function a (teta, phi)
    real, intent (in)   :: teta, phi
    real, dimension (3) :: a
  endfunction a
end interface

interface
  subroutine indices (i, j, m, factor)
    integer, intent(in)  :: i, j
    integer, intent(out) :: m
    real, intent (out)   :: factor
  endsubroutine indices
end interface

real, dimension(3) :: aux
real :: young, fact1, fact2
integer :: i, j, k, l, m, n

young = 0.0
aux = a(teta, phi)

do i = 1, 3
  do j = 1, 3
    do k = 1, 3
      do l = 1, 3
        call indices (i, j, m, fact1)
        call indices (k, l, n, fact2)
        young = young + S(m, n)*fact1*fact2*aux(i)*aux(j)*aux(k)*aux(l)
      enddo
    enddo
  enddo
enddo

young = 1.0/young

endfunction young



! ================================================================



function comp_lineal (S, teta, phi, a, indices)
! Valor del coeficiente de compresibilidad lineal en una direccion dada.

real, intent(in)                 :: teta, phi
real, intent(in), dimension(:,:) :: S

interface
  function a (teta, phi)
    real, intent (in)   :: teta, phi
    real, dimension (3) :: a
  endfunction a
end interface

interface
  subroutine indices (i, j, m, factor)
    integer, intent(in)  :: i, j
    integer, intent(out) :: m
    real, intent (out)   :: factor
  endsubroutine indices
end interface

real, dimension(3) :: aux
real :: comp_lineal, fact1
integer :: i, j, k, l, m, n

comp_lineal = 0.0
aux = a(teta, phi)

do i = 1, 3
  do j = 1, 3
    do k = 1, 3  
      call indices (i, j, m, fact1)
      comp_lineal = comp_lineal + S(m, k)*fact1 * aux(i) * aux(j)
    enddo
  enddo
enddo

endfunction comp_lineal



! ================================================================



subroutine average_voight (C_tensor, d_tensor, K, G, E, nu, d_K, d_g, d_E, d_nu)
! Propiedades elasticas y sus errores. ESQUEMA DE VOIGT

real, dimension (:, :), intent(in) :: C_tensor, d_tensor
real, intent(out)                  :: K, G, E, nu
real, intent(out)                  :: d_K, d_G, d_E, d_nu


real                               :: A, B, C
real                               :: k1, k2, k3


A = (C_tensor(1,1) + C_tensor(2,2) + C_tensor(3,3)) / 3.0
B = (C_tensor(2,3) + C_tensor(1,3) + C_tensor(1,2)) / 3.0
C = (C_tensor(4,4) + C_tensor(5,5) + C_tensor(6,6)) / 3.0

k1 = (d_tensor(1,1))**2 + (d_tensor(2,2))**2 + (d_tensor(3,3))**2
k2 = (d_tensor(1,2))**2 + (d_tensor(2,3))**2 + (d_tensor(1,3))**2
k3 = (d_tensor(4,4))**2 + (d_tensor(5,5))**2 + (d_tensor(6,6))**2

K = (A + 2*B) / 3.0
G = (A - B + 3*C) / 5.0

E = 1.0 / ( 1.0/(3*G) + 1.0/(9*K) )
nu = (1.0 - 3*G/(3*K + G)) / 2.0

d_K = sqrt(k1 + 4*k2) / 9
d_G = sqrt(k1 + k2 + 9*k3) /15

d_E = sqrt((d_G/G**2)**2 + (d_K/(3*K**2))**2) * E**2/3.0
d_nu = sqrt((G/(3*K+G)-1.0)**2*d_G**2 + (3*G/(3*K+G))*2*d_K**2) * 3/(6*K+2*G)

endsubroutine average_voight



! ================================================================



subroutine average_reuss (S, d_tensor, K, G, E, nu, d_K, d_g, d_E, d_nu)
! Propiedades elasticas y sus errores. ESQUEMA DE REUSS

real, dimension (:, :), intent(in) :: S, d_tensor
real, intent(out)                  :: K, G, E, nu
real, intent(out)                  :: d_K, d_G, d_E, d_nu


real                               :: a, b, c
real                               :: k1, k2, k3


a = (S(1,1) + S(2,2) + S(3,3)) / 3.0
b = (S(2,3) + S(1,3) + S(1,2)) / 3.0
c = (S(4,4) + S(5,5) + S(6,6)) / 3.0

k1 = (d_tensor(1,1))**2 + (d_tensor(2,2))**2 + (d_tensor(3,3))**2
k2 = (d_tensor(1,2))**2 + (d_tensor(2,3))**2 + (d_tensor(1,3))**2
k3 = (d_tensor(4,4))**2 + (d_tensor(5,5))**2 + (d_tensor(6,6))**2

K = 1.0 / (3*a + 6*b)
G = 5.0 / (4*a - 4*b + 3*c)

E = 1.0 / ( 1.0/(3*G) + 1.0/(9*K) )
nu = (1.0 - 3*G/(3*K + G)) / 2.0


d_K = K**2 * sqrt(k1 + 4*k2)
d_G = (4.0/15) * G**2 * sqrt(k1 + k2 + 9*k3)

d_E = sqrt((d_G/G**2)**2 + (d_K/(3*K**2))**2) * E**2/3.0
d_nu = sqrt((G/(3*K+G)-1.0)**2*d_G**2 + (3*G/(3*K+G))*2*d_K**2) * 3/(6*K+2*G)


endsubroutine average_reuss

endmodule m_young_compresibilidad
