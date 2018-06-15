program GENERAL

real, parameter                    :: escala = 5.0E-3
character (len = 4), dimension (6) :: aux_def
real, dimension (6,6)              :: defor_lista
real, dimension (6)                :: sigma

character (len = 8), dimension (6) :: out_lista
real, dimension (36, 13)           :: MAT
real, dimension (36)               :: B_1
real, dimension (6)                :: deformacion
real, dimension (6, 13)            :: mat_aux

real, dimension (36, 1)         :: B
real, dimension (13, 36)        :: MAT_T, pseudoinverse
real, dimension (13, 13)        :: aux1, aux2
real, dimension (13, 1)         :: z, resultado
character(len=3), dimension(13) :: constantes

integer :: i, kont, kont2



! 1. SE CREAN LOS 6 VECTORES DE DEFORMACION GENERAL QUE SE UTILIZAN EN LAS SIMULACIONES


aux_def = (/"def1", "def2", "def3", "def4", "def5", "def6"/)

defor_lista (1, 1:6) = (/-4.0, 3.0, 2.0, -3.0, -2.0, 1.0/)
defor_lista (2, 1:6) = (/-2.0, 1.0, 4.0, -3.0, 6.0, -5.0/)
defor_lista (3, 1:6) = (/3.0, -5.0, -1.0, 6.0, 2.0, -4.0/)
defor_lista (4, 1:6) = (/-4.0, -6.0, 5.0, 1.0, -3.0, 2.0/)
defor_lista (5, 1:6) = (/5.0, 4.0, 6.0, -2.0, -1.0, -3.0/)
defor_lista (6, 1:6) = (/5.0, -6.0, 1.0, -2.0, -3.0, 4.0/)

defor_lista = defor_lista * escala



! 2. SE LEEN LOS VECTORES DE ESFUERZO DE LAS 6 SIMULACIONES DE LOS ARCHIVOS 'def1.out',...,'def6.out'
!    SE CREAN LA MATRIZ 'MAT' Y EL VECTOR  'B'
!        EN LA MATRIZ 'MAT' SE ALMACENAN LAS DEFORMACIONES (CON LA FORMA PARA EL SISTEMA LINEAL)
!        EN LA MATRIZ 'B' SE ALMACENAN LOS ESFUERZOS
!    PARA QUE ASI EL SISTEMA SEA MAT*X = B


kont = 1
do i = 1, 6
  kont2 = kont + 5
  
  deformacion = defor_lista(i, :)
  sigma       = get_sigma(aux_def(i)//'.out')
  
  mat_aux = matrix_eps (deformacion)
  MAT (kont:kont2, :) = mat_aux

  B_1 (kont:kont2) = sigma
  kont = kont2 + 1
enddo



! 3. SE CALCULA LA PSEUDOINVERSA DE MOORE-PENROSE:  A^{+} = (A^{t} * A)^{-1} * A^{t}
!    Y SE RESUELVE EL SISTEMA. SOLUCION: z = A^{+} * B

B (:,1) = B_1(:)
MAT_T = transpose (MAT)
aux1 = matmul (MAT_T, MAT)
call inverse (aux1, 13, aux2)

pseudoinverse = matmul (aux2, MAT_T)

z = matmul(pseudoinverse, B)



! 4. SE GUARDAN LOS RESULTADOS EN EL ARCHIVO 'CONSTANTES.txt'


constantes = (/'C11', 'C12', 'C13', 'C16', 'C22', 'C23', 'C26', 'C33', 'C36', 'C44', 'C45', 'C55', 'C66'/)

open(unit = 3, action = 'write', file = 'CONSTANTES.txt', status = 'replace')
do i = 1, 13
  write (unit = 3, fmt = '(a10,10f12.6)'), constantes(i)//' = ', z(i,1)*160.21766208
enddo
close(unit = 3)





contains

function get_sigma(archivo)
! Dado el nombre del archivo *.out devuelve un vector con los componentes del vector de esfuerzo

character(len = *), intent(in) :: archivo
real, dimension (3,3) :: A
real, dimension (6) :: get_sigma
integer :: i

call system ("tail -n 14 "//archivo//" | head -n 3 | cut -c 10- > b.txt")

open (unit = 1, action = "read", status = "old", file = "b.txt")
do i = 1, 3
  read(unit = 1, fmt = '(3f12.6)'), A(i,1:3)
enddo
close (unit = 1)

call system ("rm b.txt")

do i = 1, 3
  get_sigma(i) = A(i,i)
enddo

get_sigma(4) = (A(2,3) + A(3,2))/2.0
get_sigma(5) = (A(1,3) + A(3,1))/2.0
get_sigma(6) = (A(2,1) + A(1,2))/2.0

endfunction get_sigma



! ================================================================================



function matrix_eps (epsilon_A)
! Devuelve una matriz de 6x13 con los coeficientes de la matriz total correspondientes a la deformacion epsilon_A

real, dimension (:), intent (in) :: epsilon_A
real, dimension (6, 13)          :: matrix_eps

matrix_eps = 0.0

matrix_eps (1, 1:3)   = epsilon_A (1:3)
matrix_eps (1, 4)     = epsilon_A (6)

matrix_eps (2, 2)     = epsilon_A (1)
matrix_eps (2, 5:6)   = epsilon_A (2:3)
matrix_eps (2, 7)     = epsilon_A (6)

matrix_eps (3, 3)     = epsilon_A (1)
matrix_eps (3, 6)     = epsilon_A (2)
matrix_eps (3, 8)     = epsilon_A (3)
matrix_eps (3, 9)     = epsilon_A (6)

matrix_eps (4, 10:11) = epsilon_A (4:5)

matrix_eps (5, 11:12) = epsilon_A (4:5)

matrix_eps (6, 4)     = epsilon_A (1)
matrix_eps (6, 7)     = epsilon_A (2)
matrix_eps (6, 9)     = epsilon_A (3)
matrix_eps (6, 13)    = epsilon_A (6)

endfunction matrix_eps



! ================================================================================



subroutine inverse(a, n, c)
! Se calcula la inversa de una matriz mediante descomposicion LU

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


endprogram GENERAL
