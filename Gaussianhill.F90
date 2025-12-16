!======================================================================
!  D2Q5 MRT-LBM for Gaussian hill advection-diffusion 
!
!  PDE (physical):
!    ∂_t φ + u_phys · ∇φ = D_phys ( ∂_xx φ + ∂_yy φ ),
!
!  Analytic solution (Gaussian hill with constant convection):
!    σ_t^2 = σ0^2 + 2 D_phys t
!    x_c(t) = ux_phys t,  y_c(t) = uy_phys t
!    φ(x,y,t) = (σ0^2 / σ_t^2) *
!               exp( -[(x-x_c)^2 + (y-y_c)^2] / (2 σ_t^2) )
!
!  LBM model:
!    - D2Q5, MRT (C = 0 case of PRE-2020, with B = φ u).
!    - Equilibrium:
!        f_eq_j = w_j [ φ + (c_j · B) / cs2 ],
!      with B = φ u_LB (LB-unit velocity).
!    - First-moment source (C=0):
!        M1G = (I - S1/2) ∂_t B,
!        G_j = w_j (c_j · M1G) / cs2,
!      where ∂_t B ≈ B^n - B^{n-1} (Δt_LB = 1).
!    - Periodic BC.
!
!  Mapping:
!    D_LB = cs2 (tau_phi - 0.5)
!    D_phys = D_LB * dx^2 / dt_phys
!      ⇒ dt_phys = D_LB * dx^2 / D_phys
!    u_LB = u_phys * dt_phys / dx
!
!======================================================================
module commondata_gauss_adv_d2q5
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none

  !----------------- mesh & domain -----------------
  integer, parameter :: nx = 64, ny = 64
  real(wp), parameter :: Lx = 2.0_wp   ! x ∈ [-1,1]
  real(wp), parameter :: Ly = 2.0_wp   ! y ∈ [-1,1]
  real(wp), parameter :: dx = Lx / real(nx,wp)
  real(wp), parameter :: dy = Ly / real(ny,wp)

  !----------------- physical & model parameters -----------------
  real(wp), parameter :: pi      = acos(-1.0_wp)
  real(wp), parameter :: D_phys  = 2.0e-4_wp   ! physical diffusion
  real(wp), parameter :: sigma0  = 0.05_wp     ! initial Gaussian width

  ! physical convection velocity
  real(wp), parameter :: ux_phys = 0.05_wp
  real(wp), parameter :: uy_phys = 0.05_wp

  real(wp), parameter :: cs2     = 1.0_wp/3.0_wp
  real(wp), parameter :: tau_phi = 0.62_wp      ! relaxation time for scalar

  ! LB diffusion coefficient
  real(wp), parameter :: D_LB = cs2 * (tau_phi - 0.5_wp)

  ! physical time step mapping
  real(wp), parameter :: dt_phys   = D_LB * dx*dx / D_phys
  real(wp), parameter :: tEnd_phys = 1.0_wp

  ! LB velocity (in lattice units)
  real(wp), parameter :: ux_LB = ux_phys * dt_phys / dx
  real(wp), parameter :: uy_LB = uy_phys * dt_phys / dy

  !----------------- D2Q5 discrete velocities -----------------
  integer, parameter :: q = 5
  integer, parameter :: ex(0:4) = (/ 0,  1,  0, -1,  0 /)
  integer, parameter :: ey(0:4) = (/ 0,  0,  1,  0, -1 /)
  real(wp), parameter :: w (0:4) = (/ 1.0_wp/3.0_wp, &
                                     1.0_wp/6.0_wp, 1.0_wp/6.0_wp, &
                                     1.0_wp/6.0_wp, 1.0_wp/6.0_wp /)

  !----------------- MRT relaxation rates -----------------
  real(wp), parameter :: s0 = 0.0_wp
  real(wp), parameter :: s1 = 1.0_wp/tau_phi   ! x-flux mode
  real(wp), parameter :: s2 = 1.0_wp/tau_phi   ! y-flux mode
  real(wp), parameter :: s3 = 1.0_wp
  real(wp), parameter :: s4 = 1.0_wp

  !----------------- fields & coordinates -----------------
  real(wp), allocatable :: phi(:,:)                 ! scalar φ
  real(wp), allocatable :: f(:,:,:), f_post(:,:,:)  ! distributions

  ! for explicit time derivative of B = φ u_LB
  real(wp), allocatable :: Bx_prev(:,:), By_prev(:,:)

  real(wp), allocatable :: x(:), y(:)               ! physical cell centers

  ! time & step control
  integer :: nSteps
  integer :: it
  real(wp) :: t_phys

contains

  !------------------------------------------------------------------
  pure function phi_exact_phys(xp, yp, tp) result(val)
    ! Analytic solution: Gaussian hill with convection
    real(wp), intent(in) :: xp, yp, tp
    real(wp) :: val
    real(wp) :: sig2, rx, ry

    sig2 = sigma0*sigma0 + 2.0_wp*D_phys*tp

    ! center moves with physical velocity (ux_phys, uy_phys)
    rx = xp - ux_phys*tp
    ry = yp - uy_phys*tp

    val = (sigma0*sigma0 / sig2) * &
          exp( -0.5_wp * (rx*rx + ry*ry) / sig2 )
  end function phi_exact_phys

end module commondata_gauss_adv_d2q5

!======================================================================
program lbm_gauss_adv_d2q5
  use omp_lib
  use commondata_gauss_adv_d2q5
  implicit none
  real(wp) :: t0, t1
  real(wp) :: errL2, errL1

  call omp_set_num_threads(omp_get_max_threads())

  call initial()

  call cpu_time(t0)

  do it = 1, nSteps

    call collision_mrt_adv()   ! MRT+convection (C=0 scheme, M1G)
    call streaming()           ! periodic BC
    call macro_phi()           ! update phi

    t_phys = t_phys + dt_phys

  end do

  call cpu_time(t1)

  call compute_L2_error(errL2)
  call compute_L1_error(errL1)

  write(*,'(A,ES16.8)') 'Final physical time t_phys = ', t_phys
  write(*,'(A,ES16.8)') 'Final relative L2 error   = ', errL2
  write(*,'(A,ES16.8)') 'Final relative L1 error   = ', errL1
  write(*,'(A,F12.4)')  'CPU time (s)              = ', real(t1 - t0,wp)

contains

  !---------------------------------
  subroutine initial()
    implicit none
    integer :: i,j
    real(wp) :: xc, yc, phi_loc
    real(wp) :: Bx_loc, By_loc

    allocate(phi(nx,ny))
    allocate(f(nx,ny,0:4))
    allocate(f_post(nx,ny,0:4))
    allocate(Bx_prev(nx,ny), By_prev(nx,ny))
    allocate(x(nx), y(ny))

    ! physical coordinates: x ∈ [-1,1], y ∈ [-1,1]
    do i = 1, nx
      x(i) = -0.5_wp*Lx + (real(i,wp)-0.5_wp)*dx
    end do
    do j = 1, ny
      y(j) = -0.5_wp*Ly + (real(j,wp)-0.5_wp)*dy
    end do

    t_phys = 0.0_wp

    nSteps = int( tEnd_phys / dt_phys + 0.5_wp )

    write(*,'(A,ES16.8)') 'User-specified t_end_phys = ', tEnd_phys
    write(*,'(A,ES16.8)') 'dt_phys (per LB step)     = ', dt_phys
    write(*,'(A,I10)')    'Number of LB steps        = ', nSteps
    write(*,'(A,ES16.8)') 'Actual t_end* = nSteps*dt = ', &
                          real(nSteps,wp)*dt_phys
    write(*,*) 'ux_phys, uy_phys = ', ux_phys, uy_phys
    write(*,*) 'ux_LB,   uy_LB   = ', ux_LB,   uy_LB
    write(*,*) '---------------------------------------------'

    ! initial field: analytic solution at t=0
!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc,phi_loc,Bx_loc,By_loc)
    do j = 1, ny
      do i = 1, nx
        xc = x(i)
        yc = y(j)

        phi_loc   = phi_exact_phys(xc,yc,t_phys)
        phi(i,j)  = phi_loc

        ! B = φ u_LB (LB units)
        Bx_loc = ux_LB * phi_loc
        By_loc = uy_LB * phi_loc

        Bx_prev(i,j) = Bx_loc
        By_prev(i,j) = By_loc

        ! equilibrium distribution: f_eq_j = w_j [ φ + (c_j·B)/cs2 ]
        f(i,j,0) = w(0) * ( phi_loc )
        f(i,j,1) = w(1) * ( phi_loc + Bx_loc/cs2 )
        f(i,j,2) = w(2) * ( phi_loc + By_loc/cs2 )
        f(i,j,3) = w(3) * ( phi_loc - Bx_loc/cs2 )
        f(i,j,4) = w(4) * ( phi_loc - By_loc/cs2 )

        f_post(i,j,0) = f(i,j,0)
        f_post(i,j,1) = f(i,j,1)
        f_post(i,j,2) = f(i,j,2)
        f_post(i,j,3) = f(i,j,3)
        f_post(i,j,4) = f(i,j,4)
      end do
    end do
!$omp end parallel do

  end subroutine initial

  !---------------------------------
  subroutine collision_mrt_adv()
    ! MRT collision in moment space + C=0 convection correction (M1G, G_j).
    !
    ! Moments:
    !   m0 = f0 + f1 + f2 + f3 + f4
    !   m1 = f1 - f3
    !   m2 = f2 - f4
    !   m3 = f1 + f3 - 2 f0
    !   m4 = f2 + f4 - 2 f0
    !
    ! Equilibrium moments (with B = φ u_LB):
    !   m0^eq = φ
    !   m1^eq = Bx / (3 cs2) = Bx   (cs2=1/3)
    !   m2^eq = By / (3 cs2) = By
    !   m3^eq = -φ/3
    !   m4^eq = -φ/3
    !
    ! C=0 scheme:
    !   M1G = (I - S1/2) ∂_t B,   ∂_t B ≈ B^n - B^{n-1}
    !   G_j = w_j (c_j · M1G) / cs2
    implicit none
    integer :: i,j,a
    real(wp) :: f0,f1,f2,f3,f4
    real(wp) :: m0,m1,m2,m3,m4
    real(wp) :: me0,me1,me2,me3,me4
    real(wp) :: ms0,ms1,ms2,ms3,ms4
    real(wp) :: phi_loc
    real(wp) :: Bx_loc, By_loc
    real(wp) :: dBx, dBy
    real(wp) :: M1Gx, M1Gy
    real(wp) :: Gj

!$omp parallel do collapse(2) default(shared) private(i,j,a, &
!$omp& f0,f1,f2,f3,f4, &
!$omp& m0,m1,m2,m3,m4, &
!$omp& me0,me1,me2,me3,me4, &
!$omp& ms0,ms1,ms2,ms3,ms4, &
!$omp& phi_loc,Bx_loc,By_loc,dBx,dBy,M1Gx,M1Gy,Gj)
    do j = 1, ny
      do i = 1, nx

        phi_loc = phi(i,j)

        ! local B in LB units
        Bx_loc = ux_LB * phi_loc
        By_loc = uy_LB * phi_loc

        ! ∂_t B ≈ B^n - B^{n-1} (LB time step = 1)
        dBx   = Bx_loc - Bx_prev(i,j)
        dBy   = By_loc - By_prev(i,j)
        Bx_prev(i,j) = Bx_loc
        By_prev(i,j) = By_loc

        ! M1G = (I - S1/2) ∂_t B
        M1Gx  = (1.0_wp - 0.5_wp*s1) * dBx
        M1Gy  = (1.0_wp - 0.5_wp*s2) * dBy

        !------ moments m = M f -------
        f0 = f(i,j,0)
        f1 = f(i,j,1)
        f2 = f(i,j,2)
        f3 = f(i,j,3)
        f4 = f(i,j,4)

        m0 = f0 + f1 + f2 + f3 + f4
        m1 = f1 - f3
        m2 = f2 - f4
        m3 = f1 + f3 - 2.0_wp*f0
        m4 = f2 + f4 - 2.0_wp*f0

        !------ equilibrium moments ----
        me0 = phi_loc
        me1 = Bx_loc            ! = Bx / (3 cs2)
        me2 = By_loc
        me3 = -phi_loc / 3.0_wp
        me4 = -phi_loc / 3.0_wp

        !------ MRT relaxation ----------
        ms0 = m0
        ms1 = m1 - s1*(m1 - me1)
        ms2 = m2 - s2*(m2 - me2)
        ms3 = m3 - s3*(m3 - me3)
        ms4 = m4 - s4*(m4 - me4)

        !------ back to velocity space: f_post = M^{-1} m* ----
        f_post(i,j,0) = (1.0_wp/5.0_wp)*ms0 &
                      - (1.0_wp/5.0_wp)*ms3 &
                      - (1.0_wp/5.0_wp)*ms4

        f_post(i,j,1) = (1.0_wp/5.0_wp)*ms0 &
                      + 0.5_wp*ms1 &
                      + 0.3_wp*ms3 &
                      - (1.0_wp/5.0_wp)*ms4

        f_post(i,j,2) = (1.0_wp/5.0_wp)*ms0 &
                      + 0.5_wp*ms2 &
                      - (1.0_wp/5.0_wp)*ms3 &
                      + 0.3_wp*ms4

        f_post(i,j,3) = (1.0_wp/5.0_wp)*ms0 &
                      - 0.5_wp*ms1 &
                      + 0.3_wp*ms3 &
                      - (1.0_wp/5.0_wp)*ms4

        f_post(i,j,4) = (1.0_wp/5.0_wp)*ms0 &
                      - 0.5_wp*ms2 &
                      - (1.0_wp/5.0_wp)*ms3 &
                      + 0.3_wp*ms4

        !------ add G_j (C=0 scheme) ----
        do a = 0, 4
          Gj = w(a) * ( ex(a)*M1Gx + ey(a)*M1Gy ) / cs2
          f_post(i,j,a) = f_post(i,j,a) + Gj
        end do

      end do
    end do
!$omp end parallel do

  end subroutine collision_mrt_adv

  !---------------------------------
  subroutine streaming()
    ! periodic streaming (pull scheme)
    implicit none
    integer :: i,j,a,ip,jp

!$omp parallel do collapse(2) default(shared) private(i,j,a,ip,jp)
    do j = 1, ny
      do i = 1, nx
        do a = 0, 4
          ip = i - ex(a)
          jp = j - ey(a)

          if (ip < 1)  ip = ip + nx
          if (ip > nx) ip = ip - nx
          if (jp < 1)  jp = jp + ny
          if (jp > ny) jp = jp - ny

          f(i,j,a) = f_post(ip,jp,a)
        end do
      end do
    end do
!$omp end parallel do
  end subroutine streaming

  !---------------------------------
  subroutine macro_phi()
    implicit none
    integer :: i,j
!$omp parallel do collapse(2) default(shared) private(i,j)
    do j = 1, ny
      do i = 1, nx
        phi(i,j) = f(i,j,0) + f(i,j,1) + f(i,j,2) + f(i,j,3) + f(i,j,4)
      end do
    end do
!$omp end parallel do
  end subroutine macro_phi

  !---------------------------------
  subroutine compute_L2_error(err)
    implicit none
    real(wp), intent(out) :: err
    integer :: i,j
    real(wp) :: num, den, xc, yc, exact

    num = 0.0_wp
    den = 0.0_wp

!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc,exact) reduction(+:num,den)
    do j = 1, ny
      do i = 1, nx
        xc = x(i)
        yc = y(j)
        exact = phi_exact_phys(xc,yc,t_phys)
        num = num + (phi(i,j) - exact)**2
        den = den + exact**2
      end do
    end do
!$omp end parallel do

    if (den > 0.0_wp) then
      err = sqrt(num/den)
    else
      err = 0.0_wp
    end if
  end subroutine compute_L2_error

  !---------------------------------
  subroutine compute_L1_error(err)
    implicit none
    real(wp), intent(out) :: err
    integer :: i,j
    real(wp) :: num, den, xc, yc, exact

    num = 0.0_wp
    den = 0.0_wp

!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc,exact) reduction(+:num,den)
    do j = 1, ny
      do i = 1, nx
        xc = x(i)
        yc = y(j)
        exact = phi_exact_phys(xc,yc,t_phys)
        num = num + abs(phi(i,j) - exact)
        den = den + abs(exact)
      end do
    end do
!$omp end parallel do

    if (den > 0.0_wp) then
      err = num/den
    else
      err = 0.0_wp
    end if
  end subroutine compute_L1_error

end program lbm_gauss_adv_d2q5
