!=============================================================================!
!  SOLUCIÓN NUMÉRICA: ECUACIONES DE POISSON Y LAPLACE                        !
!  Método de Diferencias Finitas — Iteración de Gauss-Seidel/Jacobi          !
!  Equivalente al notebook: Poisson_Laplace_Mallado.ipynb                    !
!                                                                             !
!  Compilar: gfortran -O2 -o poisson_laplace poisson_laplace.f90 -lm        !
!  Ejecutar: ./poisson_laplace                                                !
!=============================================================================!

program poisson_laplace
  implicit none

  write(*,*) "========================================="
  write(*,*) "  EJERCICIO 1: Poisson  f=(x²+y²)e^(xy) "
  write(*,*) "  Dominio: [0,2]x[0,1]                  "
  write(*,*) "  Solución analítica: V = e^(xy)         "
  write(*,*) "========================================="
  call ejercicio1()

  write(*,*)
  write(*,*) "========================================="
  write(*,*) "  EJERCICIO 2: Laplace (f=0)             "
  write(*,*) "  Dominio: [1,2]x[0,1]                  "
  write(*,*) "  Solución analítica: V = ln(x²+y²)      "
  write(*,*) "========================================="
  call ejercicio2()

  write(*,*)
  write(*,*) "========================================="
  write(*,*) "  EJERCICIO 3: Poisson  f=4              "
  write(*,*) "  Dominio: [1,2]x[0,2]                  "
  write(*,*) "  Solución analítica: V = (x-y)²         "
  write(*,*) "========================================="
  call ejercicio3()

  write(*,*)
  write(*,*) "========================================="
  write(*,*) "  EJERCICIO 4: Poisson  f=x/y + y/x      "
  write(*,*) "  Dominio: [1,2]x[1,2]                  "
  write(*,*) "  Solución analítica: V = xy·ln(xy)      "
  write(*,*) "========================================="
  call ejercicio4()

end program poisson_laplace


!=============================================================================!
! SUBRUTINA AUXILIAR: Imprime estadísticas de error y muestra la grilla      !
!=============================================================================!
subroutine print_stats(T_num, T_ana, M, N, label)
  implicit none
  integer, intent(in) :: M, N
  real(8), intent(in) :: T_num(0:N, 0:M), T_ana(0:N, 0:M)
  character(len=*), intent(in) :: label
  real(8) :: err_max, err_mean, err
  integer :: i, j

  err_max  = 0.0d0
  err_mean = 0.0d0
  do i = 0, N
    do j = 0, M
      err = abs(T_num(i,j) - T_ana(i,j))
      if (err > err_max) err_max = err
      err_mean = err_mean + err
    end do
  end do
  err_mean = err_mean / real((M+1)*(N+1), 8)

  write(*,'(A,A)')       "  Resultado: ", label
  write(*,'(A,ES12.4)')  "  Error máximo    : ", err_max
  write(*,'(A,ES12.4)')  "  Error promedio  : ", err_mean
end subroutine print_stats


!=============================================================================!
! SUBRUTINA AUXILIAR: Escribe la grilla completa en un archivo .dat           !
!=============================================================================!
subroutine write_grid(filename, X, Y, T_num, T_ana, M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), intent(in) :: X(0:N,0:M), Y(0:N,0:M)
  real(8), intent(in) :: T_num(0:N,0:M), T_ana(0:N,0:M)
  character(len=*), intent(in) :: filename
  integer :: i, j, funit
  funit = 20
  open(funit, file=filename, status='replace')
  write(funit,'(A)') "# x   y   V_numerica   V_analitica   |error|"
  do i = 0, N
    do j = 0, M
      write(funit,'(5ES15.6)') X(i,j), Y(i,j), T_num(i,j), T_ana(i,j), &
                                abs(T_num(i,j)-T_ana(i,j))
    end do
    write(funit,*)   ! línea en blanco por filas (para gnuplot splot)
  end do
  close(funit)
  write(*,'(A,A)') "  Datos escritos en: ", filename
end subroutine write_grid


!=============================================================================!
! EJERCICIO 1                                                                 !
! ∇²V = (x²+y²)·e^(xy)   en [0,2]x[0,1]                                    !
! CC:  V(0,y)=1, V(2,y)=e^(2y), V(x,0)=1, V(x,1)=e^x                       !
! Analítica: V(x,y) = e^(xy)                                                 !
!=============================================================================!
subroutine ejercicio1()
  implicit none
  integer, parameter :: M = 50, N = 50
  real(8), parameter :: x0=0.0d0, xf=2.0d0, y0=0.0d0, yf=1.0d0
  real(8), parameter :: tol = 1.0d-6
  real(8) :: h, k, h2, k2, denom
  real(8) :: X(0:N,0:M), Y(0:N,0:M)
  real(8) :: T(0:N,0:M), T_new(0:N,0:M), T_ana(0:N,0:M), source(0:N,0:M)
  integer :: i, j, it
  real(8) :: delta

  h = (xf - x0) / M
  k = (yf - y0) / N
  h2 = h*h;  k2 = k*k
  denom = 2.0d0*(h2 + k2)

  ! Construir grilla
  do i = 0, N
    do j = 0, M
      X(i,j) = x0 + j*h
      Y(i,j) = y0 + i*k
    end do
  end do

  ! Inicializar T con ceros
  T = 0.0d0

  ! Condiciones de contorno
  do i = 0, N
    T(i, 0)  = 1.0d0                       ! V(0,y) = 1
    T(i, M)  = exp(2.0d0 * Y(i,M))         ! V(2,y) = e^(2y)
  end do
  do j = 0, M
    T(0, j)  = 1.0d0                        ! V(x,0) = 1
    T(N, j)  = exp(X(N,j))                  ! V(x,1) = e^x
  end do

  ! Término fuente: f(x,y) = (x²+y²)·e^(xy)
  do i = 0, N
    do j = 0, M
      source(i,j) = (X(i,j)**2 + Y(i,j)**2) * exp(X(i,j)*Y(i,j))
    end do
  end do

  ! ---- Iteración de Jacobi ----
  T_new = T
  do it = 1, 100000
    delta = 0.0d0
    do i = 1, N-1
      do j = 1, M-1
        T_new(i,j) = ( (T(i+1,j) + T(i-1,j))*k2 &
                     + (T(i,j+1) + T(i,j-1))*h2  &
                     - source(i,j)*h2*k2          ) / denom
        delta = max(delta, abs(T_new(i,j) - T(i,j)))
      end do
    end do
    T = T_new
    if (delta < tol) then
      write(*,'(A,I6,A)') "  Convergió en ", it, " iteraciones."
      exit
    end if
  end do

  ! Solución analítica
  do i = 0, N
    do j = 0, M
      T_ana(i,j) = exp(X(i,j)*Y(i,j))
    end do
  end do

  call print_stats(T, T_ana, M, N, "Ejercicio 1 — Poisson")
  call write_grid("ejercicio1.dat", X, Y, T, T_ana, M, N)
end subroutine ejercicio1


!=============================================================================!
! EJERCICIO 2                                                                 !
! ∇²V = 0   (Laplace)  en [1,2]x[0,1]                                       !
! CC:  V(1,y)=ln(y²+1), V(2,y)=ln(y²+4), V(x,0)=2ln(x), V(x,1)=ln(x²+1)  !
! Analítica: V(x,y) = ln(x²+y²)                                              !
!=============================================================================!
subroutine ejercicio2()
  implicit none
  integer, parameter :: M = 50, N = 50
  real(8), parameter :: x0=1.0d0, xf=2.0d0, y0=0.0d0, yf=1.0d0
  real(8), parameter :: tol = 1.0d-6
  real(8) :: h, k, h2, k2, denom
  real(8) :: X(0:N,0:M), Y(0:N,0:M)
  real(8) :: T(0:N,0:M), T_new(0:N,0:M), T_ana(0:N,0:M)
  integer :: i, j, it
  real(8) :: delta

  h = (xf - x0) / M
  k = (yf - y0) / N
  h2 = h*h;  k2 = k*k
  denom = 2.0d0*(h2 + k2)

  do i = 0, N
    do j = 0, M
      X(i,j) = x0 + j*h
      Y(i,j) = y0 + i*k
    end do
  end do

  T = 0.0d0

  ! Condiciones de contorno
  do i = 0, N
    T(i, 0)  = log(Y(i,0)**2 + 1.0d0)   ! V(1,y) = ln(y²+1)
    T(i, M)  = log(Y(i,M)**2 + 4.0d0)   ! V(2,y) = ln(y²+4)
  end do
  do j = 0, M
    T(0, j)  = 2.0d0*log(X(0,j))         ! V(x,0) = 2ln(x)
    T(N, j)  = log(X(N,j)**2 + 1.0d0)    ! V(x,1) = ln(x²+1)
  end do

  ! Laplace: fuente = 0  → stencil estándar con h=k simplificado
  ! Para h≠k se mantiene la fórmula general (f=0)
  T_new = T
  do it = 1, 100000
    delta = 0.0d0
    do i = 1, N-1
      do j = 1, M-1
        T_new(i,j) = ( (T(i+1,j) + T(i-1,j))*k2 &
                     + (T(i,j+1) + T(i,j-1))*h2  ) / denom
        delta = max(delta, abs(T_new(i,j) - T(i,j)))
      end do
    end do
    T = T_new
    if (delta < tol) then
      write(*,'(A,I6,A)') "  Convergió en ", it, " iteraciones."
      exit
    end if
  end do

  do i = 0, N
    do j = 0, M
      T_ana(i,j) = log(X(i,j)**2 + Y(i,j)**2)
    end do
  end do

  call print_stats(T, T_ana, M, N, "Ejercicio 2 — Laplace")
  call write_grid("ejercicio2.dat", X, Y, T, T_ana, M, N)
end subroutine ejercicio2


!=============================================================================!
! EJERCICIO 3                                                                 !
! ∇²V = 4   en [1,2]x[0,2]                                                  !
! CC consistentes con la solución analítica V(x,y) = (x-y)²                 !
!=============================================================================!
subroutine ejercicio3()
  implicit none
  integer, parameter :: M = 50, N = 50
  real(8), parameter :: x0=1.0d0, xf=2.0d0, y0=0.0d0, yf=2.0d0
  real(8), parameter :: tol = 1.0d-6
  real(8) :: h, k, h2, k2, denom
  real(8) :: X(0:N,0:M), Y(0:N,0:M)
  real(8) :: T(0:N,0:M), T_new(0:N,0:M), T_ana(0:N,0:M), source(0:N,0:M)
  integer :: i, j, it
  real(8) :: delta

  h = (xf - x0) / M
  k = (yf - y0) / N
  h2 = h*h;  k2 = k*k
  denom = 2.0d0*(h2 + k2)

  do i = 0, N
    do j = 0, M
      X(i,j) = x0 + j*h
      Y(i,j) = y0 + i*k
    end do
  end do

  T = 0.0d0

  ! Condiciones de contorno: V = (x-y)²
  do i = 0, N
    T(i, 0)  = (x0 - Y(i,0))**2    ! pared izquierda x=1
    T(i, M)  = (xf - Y(i,M))**2    ! pared derecha   x=2
  end do
  do j = 0, M
    T(0, j)  = (X(0,j) - y0)**2    ! pared inferior  y=0
    T(N, j)  = (X(N,j) - yf)**2    ! pared superior  y=2
  end do

  ! Fuente: f = 4
  source = 4.0d0

  T_new = T
  do it = 1, 100000
    delta = 0.0d0
    do i = 1, N-1
      do j = 1, M-1
        T_new(i,j) = ( (T(i+1,j) + T(i-1,j))*k2 &
                     + (T(i,j+1) + T(i,j-1))*h2  &
                     - source(i,j)*h2*k2          ) / denom
        delta = max(delta, abs(T_new(i,j) - T(i,j)))
      end do
    end do
    T = T_new
    if (delta < tol) then
      write(*,'(A,I6,A)') "  Convergió en ", it, " iteraciones."
      exit
    end if
  end do

  do i = 0, N
    do j = 0, M
      T_ana(i,j) = (X(i,j) - Y(i,j))**2
    end do
  end do

  call print_stats(T, T_ana, M, N, "Ejercicio 3 — Poisson f=4")
  call write_grid("ejercicio3.dat", X, Y, T, T_ana, M, N)
end subroutine ejercicio3


!=============================================================================!
! EJERCICIO 4                                                                 !
! ∇²V = x/y + y/x   en [1,2]x[1,2]                                          !
! CC:  V(1,y)=y·ln(y), V(2,y)=2y·ln(2y), V(x,1)=x·ln(x), V(x,2)=x·ln(4x²)!
! Analítica: V(x,y) = xy·ln(xy)                                              !
!=============================================================================!
subroutine ejercicio4()
  implicit none
  integer, parameter :: M = 50, N = 50
  real(8), parameter :: x0=1.0d0, xf=2.0d0, y0=1.0d0, yf=2.0d0
  real(8), parameter :: tol = 1.0d-6
  real(8) :: h, k, h2, k2, denom
  real(8) :: X(0:N,0:M), Y(0:N,0:M)
  real(8) :: T(0:N,0:M), T_new(0:N,0:M), T_ana(0:N,0:M), source(0:N,0:M)
  integer :: i, j, it
  real(8) :: delta

  h = (xf - x0) / M
  k = (yf - y0) / N
  h2 = h*h;  k2 = k*k
  denom = 2.0d0*(h2 + k2)

  do i = 0, N
    do j = 0, M
      X(i,j) = x0 + j*h
      Y(i,j) = y0 + i*k
    end do
  end do

  T = 0.0d0

  ! Condiciones de contorno
  do i = 0, N
    T(i, 0)  = Y(i,0) * log(Y(i,0))             ! V(1,y)  = y·ln(y)
    T(i, M)  = 2.0d0*Y(i,M)*log(2.0d0*Y(i,M))  ! V(2,y)  = 2y·ln(2y)
  end do
  do j = 0, M
    T(0, j)  = X(0,j) * log(X(0,j))             ! V(x,1)  = x·ln(x)
    T(N, j)  = X(N,j) * log(4.0d0*X(N,j)**2)    ! V(x,2)  = x·ln(4x²)
  end do

  ! Término fuente: f = x/y + y/x
  do i = 0, N
    do j = 0, M
      source(i,j) = X(i,j)/Y(i,j) + Y(i,j)/X(i,j)
    end do
  end do

  T_new = T
  do it = 1, 100000
    delta = 0.0d0
    do i = 1, N-1
      do j = 1, M-1
        T_new(i,j) = ( (T(i+1,j) + T(i-1,j))*k2 &
                     + (T(i,j+1) + T(i,j-1))*h2  &
                     - source(i,j)*h2*k2          ) / denom
        delta = max(delta, abs(T_new(i,j) - T(i,j)))
      end do
    end do
    T = T_new
    if (delta < tol) then
      write(*,'(A,I6,A)') "  Convergió en ", it, " iteraciones."
      exit
    end if
  end do

  do i = 0, N
    do j = 0, M
      T_ana(i,j) = X(i,j)*Y(i,j)*log(X(i,j)*Y(i,j))
    end do
  end do

  call print_stats(T, T_ana, M, N, "Ejercicio 4 — Poisson f=x/y+y/x")
  call write_grid("ejercicio4.dat", X, Y, T, T_ana, M, N)
end subroutine ejercicio4
