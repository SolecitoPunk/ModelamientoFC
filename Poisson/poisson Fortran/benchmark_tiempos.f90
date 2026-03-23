!=============================================================================!
!  BENCHMARK: Rendimiento Computacional — Poisson & Laplace                  !
!  Mide el tiempo de ejecución para mallas 400x400, 700x700, 1000x1000       !
!  en los 4 ejercicios y escribe los resultados en tiempos_fortran.dat        !
!                                                                             !
!  Compilar:  gfortran -O2 -o benchmark benchmark_tiempos.f90                !
!  Ejecutar:  ./benchmark                                                     !
!=============================================================================!

program benchmark_tiempos
  implicit none

  ! Tamaños de malla a probar (igual que el notebook Python)
  integer, parameter :: N_MALLAS = 3
  integer :: mallas(N_MALLAS)
  data mallas / 400, 700, 1000 /

  ! Matriz de tiempos: fila=ejercicio, col=malla
  real(8) :: tiempos(4, N_MALLAS)

  integer :: im, M, N
  real(8) :: t_inicio, t_fin

  write(*,'(A)') "=================================================="
  write(*,'(A)') "  BENCHMARK — Diferencias Finitas (Fortran)"
  write(*,'(A)') "  Ejercicios 1-4 | Mallas: 400, 700, 1000"
  write(*,'(A)') "=================================================="
  write(*,*)

  do im = 1, N_MALLAS
    M = mallas(im)
    N = mallas(im)
    write(*,'(A,I4,A,I4)') "  >> Malla ", M, " x ", N

    ! ── Ejercicio 1: Poisson  f=(x²+y²)e^(xy)  ─────────────────────────────
    call cpu_time(t_inicio)
    call ej1(M, N)
    call cpu_time(t_fin)
    tiempos(1, im) = t_fin - t_inicio
    write(*,'(A,F10.4,A)') "     Ej1: ", tiempos(1,im), " s"

    ! ── Ejercicio 2: Laplace  f=0  ──────────────────────────────────────────
    call cpu_time(t_inicio)
    call ej2(M, N)
    call cpu_time(t_fin)
    tiempos(2, im) = t_fin - t_inicio
    write(*,'(A,F10.4,A)') "     Ej2: ", tiempos(2,im), " s"

    ! ── Ejercicio 3: Poisson  f=4  ──────────────────────────────────────────
    call cpu_time(t_inicio)
    call ej3(M, N)
    call cpu_time(t_fin)
    tiempos(3, im) = t_fin - t_inicio
    write(*,'(A,F10.4,A)') "     Ej3: ", tiempos(3,im), " s"

    ! ── Ejercicio 4: Poisson  f=x/y+y/x  ───────────────────────────────────
    call cpu_time(t_inicio)
    call ej4(M, N)
    call cpu_time(t_fin)
    tiempos(4, im) = t_fin - t_inicio
    write(*,'(A,F10.4,A)') "     Ej4: ", tiempos(4,im), " s"

    write(*,*)
  end do

  ! ── Escribir .dat para el script Python de gráficos ──────────────────────
  call escribir_dat(tiempos, mallas, N_MALLAS)

  ! ── Tabla resumen en consola ──────────────────────────────────────────────
  call tabla_resumen(tiempos, mallas, N_MALLAS)

end program benchmark_tiempos


!=============================================================================!
!  EJERCICIO 1 — Poisson: f = (x²+y²)·e^(xy)                                !
!  Dominio [0,2]×[0,1]    Analítica: e^(xy)                                  !
!=============================================================================!
subroutine ej1(M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), parameter  :: x0=0.0d0, xf=2.0d0, y0=0.0d0, yf=1.0d0, tol=1.0d-6
  real(8), allocatable :: X(:,:), Y(:,:), T(:,:), Tn(:,:), src(:,:)
  real(8) :: h, k, h2, k2, denom, delta
  integer :: i, j, it

  allocate(X(0:N,0:M), Y(0:N,0:M), T(0:N,0:M), Tn(0:N,0:M), src(0:N,0:M))

  h = (xf-x0)/M;  k = (yf-y0)/N
  h2 = h*h;  k2 = k*k;  denom = 2.0d0*(h2+k2)

  do i=0,N; do j=0,M
    X(i,j) = x0+j*h;  Y(i,j) = y0+i*k
  end do; end do

  T = 0.0d0
  do i=0,N
    T(i,0) = 1.0d0
    T(i,M) = exp(2.0d0*Y(i,M))
  end do
  do j=0,M
    T(0,j) = 1.0d0
    T(N,j) = exp(X(N,j))
  end do

  do i=0,N; do j=0,M
    src(i,j) = (X(i,j)**2 + Y(i,j)**2)*exp(X(i,j)*Y(i,j))
  end do; end do

  Tn = T
  do it=1,200000
    delta = 0.0d0
    do i=1,N-1; do j=1,M-1
      Tn(i,j) = ((T(i+1,j)+T(i-1,j))*k2 + (T(i,j+1)+T(i,j-1))*h2 &
                 - src(i,j)*h2*k2) / denom
      delta = max(delta, abs(Tn(i,j)-T(i,j)))
    end do; end do
    T = Tn
    if (delta < tol) exit
  end do

  deallocate(X, Y, T, Tn, src)
end subroutine ej1


!=============================================================================!
!  EJERCICIO 2 — Laplace: f = 0                                              !
!  Dominio [1,2]×[0,1]    Analítica: ln(x²+y²)                              !
!=============================================================================!
subroutine ej2(M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), parameter  :: x0=1.0d0, xf=2.0d0, y0=0.0d0, yf=1.0d0, tol=1.0d-6
  real(8), allocatable :: X(:,:), Y(:,:), T(:,:), Tn(:,:)
  real(8) :: h, k, h2, k2, denom, delta
  integer :: i, j, it

  allocate(X(0:N,0:M), Y(0:N,0:M), T(0:N,0:M), Tn(0:N,0:M))

  h = (xf-x0)/M;  k = (yf-y0)/N
  h2 = h*h;  k2 = k*k;  denom = 2.0d0*(h2+k2)

  do i=0,N; do j=0,M
    X(i,j) = x0+j*h;  Y(i,j) = y0+i*k
  end do; end do

  T = 0.0d0
  do i=0,N
    T(i,0) = log(Y(i,0)**2 + 1.0d0)
    T(i,M) = log(Y(i,M)**2 + 4.0d0)
  end do
  do j=0,M
    T(0,j) = 2.0d0*log(X(0,j))
    T(N,j) = log(X(N,j)**2 + 1.0d0)
  end do

  Tn = T
  do it=1,200000
    delta = 0.0d0
    do i=1,N-1; do j=1,M-1
      Tn(i,j) = ((T(i+1,j)+T(i-1,j))*k2 + (T(i,j+1)+T(i,j-1))*h2) / denom
      delta = max(delta, abs(Tn(i,j)-T(i,j)))
    end do; end do
    T = Tn
    if (delta < tol) exit
  end do

  deallocate(X, Y, T, Tn)
end subroutine ej2


!=============================================================================!
!  EJERCICIO 3 — Poisson: f = 4                                              !
!  Dominio [1,2]×[0,2]    Analítica: (x-y)²                                 !
!=============================================================================!
subroutine ej3(M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), parameter  :: x0=1.0d0, xf=2.0d0, y0=0.0d0, yf=2.0d0, tol=1.0d-6
  real(8), allocatable :: X(:,:), Y(:,:), T(:,:), Tn(:,:)
  real(8) :: h, k, h2, k2, denom, delta
  integer :: i, j, it

  allocate(X(0:N,0:M), Y(0:N,0:M), T(0:N,0:M), Tn(0:N,0:M))

  h = (xf-x0)/M;  k = (yf-y0)/N
  h2 = h*h;  k2 = k*k;  denom = 2.0d0*(h2+k2)

  do i=0,N; do j=0,M
    X(i,j) = x0+j*h;  Y(i,j) = y0+i*k
  end do; end do

  T = 0.0d0
  do i=0,N
    T(i,0) = (x0 - Y(i,0))**2
    T(i,M) = (xf - Y(i,M))**2
  end do
  do j=0,M
    T(0,j) = (X(0,j) - y0)**2
    T(N,j) = (X(N,j) - yf)**2
  end do

  Tn = T
  do it=1,200000
    delta = 0.0d0
    do i=1,N-1; do j=1,M-1
      Tn(i,j) = ((T(i+1,j)+T(i-1,j))*k2 + (T(i,j+1)+T(i,j-1))*h2 &
                 - 4.0d0*h2*k2) / denom
      delta = max(delta, abs(Tn(i,j)-T(i,j)))
    end do; end do
    T = Tn
    if (delta < tol) exit
  end do

  deallocate(X, Y, T, Tn)
end subroutine ej3


!=============================================================================!
!  EJERCICIO 4 — Poisson: f = x/y + y/x                                     !
!  Dominio [1,2]×[1,2]    Analítica: xy·ln(xy)                              !
!=============================================================================!
subroutine ej4(M, N)
  implicit none
  integer, intent(in) :: M, N
  real(8), parameter  :: x0=1.0d0, xf=2.0d0, y0=1.0d0, yf=2.0d0, tol=1.0d-6
  real(8), allocatable :: X(:,:), Y(:,:), T(:,:), Tn(:,:), src(:,:)
  real(8) :: h, k, h2, k2, denom, delta
  integer :: i, j, it

  allocate(X(0:N,0:M), Y(0:N,0:M), T(0:N,0:M), Tn(0:N,0:M), src(0:N,0:M))

  h = (xf-x0)/M;  k = (yf-y0)/N
  h2 = h*h;  k2 = k*k;  denom = 2.0d0*(h2+k2)

  do i=0,N; do j=0,M
    X(i,j) = x0+j*h;  Y(i,j) = y0+i*k
  end do; end do

  T = 0.0d0
  do i=0,N
    T(i,0) = Y(i,0)*log(Y(i,0))
    T(i,M) = 2.0d0*Y(i,M)*log(2.0d0*Y(i,M))
  end do
  do j=0,M
    T(0,j) = X(0,j)*log(X(0,j))
    T(N,j) = X(N,j)*log(4.0d0*X(N,j)**2)
  end do

  do i=0,N; do j=0,M
    src(i,j) = X(i,j)/Y(i,j) + Y(i,j)/X(i,j)
  end do; end do

  Tn = T
  do it=1,200000
    delta = 0.0d0
    do i=1,N-1; do j=1,M-1
      Tn(i,j) = ((T(i+1,j)+T(i-1,j))*k2 + (T(i,j+1)+T(i,j-1))*h2 &
                 - src(i,j)*h2*k2) / denom
      delta = max(delta, abs(Tn(i,j)-T(i,j)))
    end do; end do
    T = Tn
    if (delta < tol) exit
  end do

  deallocate(X, Y, T, Tn, src)
end subroutine ej4


!=============================================================================!
!  ESCRIBIR .dat  para el script Python de visualización                      !
!=============================================================================!
subroutine escribir_dat(tiempos, mallas, nm)
  implicit none
  integer, intent(in) :: nm, mallas(nm)
  real(8), intent(in) :: tiempos(4, nm)
  integer :: im, ie

  open(10, file='tiempos_fortran.dat', status='replace')
  write(10,'(A)') "# malla   ej1       ej2       ej3       ej4"
  do im = 1, nm
    write(10,'(I6,4F12.4)') mallas(im), (tiempos(ie,im), ie=1,4)
  end do
  close(10)
  write(*,'(A)') "  >> Datos guardados en: tiempos_fortran.dat"
end subroutine escribir_dat


!=============================================================================!
!  TABLA RESUMEN en consola                                                   !
!=============================================================================!
subroutine tabla_resumen(tiempos, mallas, nm)
  implicit none
  integer, intent(in) :: nm, mallas(nm)
  real(8), intent(in) :: tiempos(4, nm)
  integer :: im, ie
  real(8) :: total(4)

  write(*,*)
  write(*,'(A)') "=================================================="
  write(*,'(A)') "  TABLA RESUMEN DE TIEMPOS (segundos)"
  write(*,'(A)') "=================================================="
  write(*,'(A6,4A12)') "Malla", "Ej1", "Ej2", "Ej3", "Ej4"
  write(*,'(A)') "--------------------------------------------------"
  do im = 1, nm
    write(*,'(I6,4F12.4)') mallas(im), (tiempos(ie,im), ie=1,4)
  end do
  write(*,'(A)') "--------------------------------------------------"

  total = 0.0d0
  do im = 1, nm
    do ie = 1, 4
      total(ie) = total(ie) + tiempos(ie,im)
    end do
  end do
  write(*,'(A6,4F12.4)') "TOTAL", (total(ie), ie=1,4)
  write(*,'(A)') "=================================================="
  write(*,*)

  ! Factor de aceleración respecto a los tiempos Python del notebook
  write(*,'(A)') "  SPEEDUP vs Python (tiempos del notebook):"
  write(*,'(A)') "--------------------------------------------------"
  block
    real(8) :: py(4,3)
    real(8) :: speedup
    integer :: col
    ! Tiempos Python del notebook original
    py(1,:) = [16.3305d0, 42.2586d0, 86.4030d0]
    py(2,:) = [ 5.6702d0, 14.6251d0, 35.9491d0]
    py(3,:) = [ 5.2973d0, 13.7404d0, 32.2162d0]
    py(4,:) = [ 5.8790d0, 18.3472d0, 36.4926d0]
    write(*,'(A14,3A12)') "Ejercicio", "400x400", "700x700", "1000x1000"
    do ie = 1, 4
      write(*,'(A,I1,A,3F12.1,A)') "  Ejercicio ", ie, ":", &
        (py(ie,col)/tiempos(ie,col), col=1,nm), "  x"
    end do
  end block
  write(*,'(A)') "=================================================="
end subroutine tabla_resumen
