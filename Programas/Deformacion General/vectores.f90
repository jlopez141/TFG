program vectores
real, parameter :: pi = acos(-1.0), escala = 5.0E-3
real, dimension (6) :: vec, resultado
real, dimension (6, 6) :: defor_lista

integer :: i, j

vec = (/6.812110, 7.541201, 22.543466, 89.9949, 89.9954, 123.2664/)

defor_lista (1, 1:6) = (/-4.0, 3.0, 2.0, -3.0, -2.0, 1.0/)
defor_lista (2, 1:6) = (/-2.0, 1.0, 4.0, -3.0, 6.0, -5.0/)
defor_lista (3, 1:6) = (/3.0, -5.0, -1.0, 6.0, 2.0, -4.0/)
defor_lista (4, 1:6) = (/-4.0, -6.0, 5.0, 1.0, -3.0, 2.0/)
defor_lista (5, 1:6) = (/5.0, 4.0, 6.0, -2.0, -1.0, -3.0/)
defor_lista (6, 1:6) = (/5.0, -6.0, 1.0, -2.0, -3.0, 4.0/)


defor_lista = defor_lista * escala


do i = 1, 6
  resultado = deformacion (vec, defor_lista (i, 1:6))
  write (unit = 1, fmt = "(i5, a, 3f12.6, 3f10.4)"), i, "  --->", resultado(1:6)
enddo

close (unit = 1)

contains

function deformacion (vector, eps)
! eps = porcentaje de la deformacion a aplicar en cada una de las direcciones
! La funcion devuelve un vector con los vectores de base deformados.

real, dimension (:), intent (in) :: vector, eps

real, dimension (size(vector)) :: deformacion
integer :: i

do i = 1, 3
  deformacion (i)   = vector (i) * (1.0 + eps (i))
  deformacion (i+3) = vector (i+3) - eps (i+3) * 180.0/pi
enddo

endfunction deformacion

endprogram vectores
