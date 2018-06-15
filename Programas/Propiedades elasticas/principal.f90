! PARA COMPILAR:
! f95 funciones.f90 propiedades.f90 principal.f90

! INPUT:
! 'matriz.txt' :: archivo con el tensor de constante elasticas
! 'error.txt'  :: archivo con la matriz de errores

program PRINCIPAL
use mfunciones
use m_young_compresibilidad


! PRIMERO ELEGIR EL NOMBRE DEL ARCHIVO CON EL TENSOR C_ij

! Magnitudes a elegir para representar :: (nombre de la funcion)
! Modulo de Young :: young
! Compresibilidad lineal :: comp_lineal


! subroutine average_voight (C, error_c, K_v, G_v, E_v, nu_v)
! subroutine average_reuss (S, error_s, K_r, G_r, E_r, nu_r)


real, dimension (6,6) :: C, S, error_C, error_S
real                  :: K_v, G_v, E_v, nu_v, K_r, G_r, E_r, nu_r
real                  :: d_K_v, d_G_v, d_E_v, d_nu_v, d_K_r, d_G_r, d_E_r,d_nu_r

integer :: i, j
real, dimension (6) ::kk




C = read_matrix ("matriz.txt", 101)
error_C = read_matrix ("error_C.txt", 102)


call inverse (C, 6, S)
call error_inverse (error_C, S, error_S)


call average_voight (C, error_C, K_v, G_v, E_v, nu_v, d_K_v, d_G_v, d_E_v, d_nu_v)
call average_reuss (S, error_S, K_r, G_r, E_r, nu_r, d_K_r, d_G_r, d_E_r, d_nu_r)

write(unit = *, fmt = *), "Bulk (K)  |  Young (E)  |  Shear (G)  |  Poisson (nu)"
write(unit = *, fmt = *), "VOIGHT"
write(unit = *, fmt = "(10f12.6)"), k_v, e_v, G_v, nu_v

write(unit = *, fmt = *), "REUSS"
write(unit = *, fmt = "(10f12.6)"), k_r, e_r, G_r, nu_r

write(unit = *, fmt = *), "HILL"
write(unit = *, fmt = "(10f12.6)"), (k_v+k_r)/2, (e_v+e_r)/2, (g_v+g_r)/2, (nu_v+nu_r)/2

write(unit = *, fmt = *), ""
write(unit = *, fmt = *), "-----------------------------"
write(unit = *, fmt = *), ""
write(unit = *, fmt = *), "ERRORES:"

write(unit = *, fmt = *), "Bulk (K)  |  Young (E)  |  Shear (G)  |  Poisson (nu)"
write(unit = *, fmt = *), "VOIGHT"
write(unit = *, fmt = "(10f12.6)"), d_k_v, d_e_v, d_G_v, d_nu_v

write(unit = *, fmt = *), "REUSS"
write(unit = *, fmt = "(10f12.6)"), d_k_r, d_e_r, d_G_r, d_nu_r

write(unit = *, fmt = *), "HILL"
write(unit = *, fmt = "(10f12.6)"),(d_k_v+d_k_r)/2,(d_e_v+d_e_r)/2,(d_g_v+d_g_r)/2,(d_nu_v+d_nu_r)/2





call plot_2D (young, S, 2, "young")
call plot_3D (young, S, 3, "young")



call plot_2D (comp_lineal, S, 4, "compr")
call plot_3D (comp_lineal, S, 5, "compr")





contains


function read_matrix (archivo, N)
! Leer una matriz de 6x6 de un archivo.

character (len = *), intent(in) :: archivo
integer, intent(in)              :: N              ! unit del que leer
real, dimension(6,6)             :: read_matrix
integer                          :: i

open (unit = N, action = "read", status = "old", file = trim(archivo))

do i = 1, 6
  read (unit = N, fmt = *) read_matrix(i,1:6)
enddo

close (unit = N)

end function read_matrix



! ================================================================



subroutine plot_2D (funcion, MAT, n_unit, nombre)
! Guarda en el archivo 2D_plot.txt los datos para representar las proyecciones xy, xz y yz
! PROYECCION XY: Columnas 1 y 2  -->  x   y
! PROYECCION XZ: Columnas 3 y 4  -->  x   z
! PROYECCION YZ: columnas 5 y 6  -->  y   z

integer, intent (in)             :: n_unit
real, dimension(:,:), intent(in) :: MAT
character (len=*), intent(in)    :: nombre

interface
  function funcion (MAT, teta, phi, a, indices)
    real, intent(in)                 :: teta, phi
    real, intent(in), dimension(:,:) :: MAT
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
  endfunction
endinterface

real               :: teta_xy, phi_xy, teta_xz_yz, phi_xz, phi_yz, d_ang
real               :: mag_xy, mag_xz, mag_yz
integer            :: N, i
real, dimension(3) :: a_xy, a_xz, a_yz

N = 300
d_ang = 2*pi/(N-1)

open (unit = n_unit, action = "write", status = "replace", file = "2D_plot_"//trim(nombre)//".txt")

teta_xy = pi/2
phi_xy = 0.0

teta_xz_yz = 0.0
phi_xz = 0.0
phi_yz = pi/2

do i = 1, N
  a_xy = a(teta_xy, phi_xy)
  mag_xy = funcion (MAT, teta_xy, phi_xy, a, indices)

  a_xz = a(teta_xz_yz, phi_xz)
  mag_xz = funcion (MAT, teta_xz_yz, phi_xz, a, indices)

  a_yz = a(teta_xz_yz, phi_yz)
  mag_yz = funcion (MAT, teta_xz_yz, phi_yz, a, indices)

  write(unit=n_unit,fmt="(10f12.6)"),mag_xy*a_xy(1),mag_xy*a_xy(2),mag_xz*a_xz(1),mag_xz*a_xz(3),mag_yz*a_yz(2),mag_yz*a_yz(3)

  teta_xz_yz = teta_xz_yz + d_ang
  phi_xy = phi_xy + d_ang
enddo

close (unit = n_unit)

end subroutine plot_2D



! ================================================================



subroutine plot_3D (funcion, MAT, n_unit, nombre)
! Guarda en el archivo 3D_plot.txt los datos para representar en tres dimensiones la magnitud.

integer, intent (in)             :: n_unit
real, dimension(:,:), intent(in) :: MAT
character (len = *), intent(in)  :: nombre

interface
  function funcion (MAT, teta, phi, a, indices)
    real, intent(in)                 :: teta, phi
    real, intent(in), dimension(:,:) :: MAT
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
  endfunction
endinterface

real, dimension(3) :: aux
real               :: teta, phi, d_teta, d_phi, magnitud
integer            :: N, i, j

N = 250

d_teta = pi / (N-1)
d_phi = 2*pi / (N-1)

teta = 0.0
phi = 0.0

open(unit = n_unit, action = "write", status = "replace", file = "3D_plot_"//trim(nombre)//".txt")

do i = 1, N
  do j = 1, N
    aux = a(teta, phi)
    magnitud = funcion (MAT, teta, phi, a, indices)
!    write (unit = n_unit, fmt = "(10f12.6)"), magnitud * aux(1), magnitud * aux(2), magnitud * aux(3)
    write (unit = n_unit, fmt = "(2i5, 10f12.6)"), i, j, magnitud * aux(1), magnitud * aux(2), magnitud * aux(3)
    phi = phi + d_phi
  enddo
  teta = teta + d_teta
enddo

close (unit = n_unit)

end subroutine plot_3D




end program PRINCIPAL
