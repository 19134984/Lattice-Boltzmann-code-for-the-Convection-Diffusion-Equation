!======================================================================
!  PRE2020 Case-1 MRT-LBM (D2Q5) for Example 3.1 (periodic BC)
!  - 2D convection–diffusion with source term
!  - Case 1: beta=1, D = phi I, C = 0, M1F = 0, gamma = 0
!  - D2Q5 + MRT for scalar, OpenMP parallel
!  - Domain: [0,2] x [0,2], nx=ny=256
!  - Pe = 100, u_phys = (15,15), tau_phi = 0.56
!======================================================================
module commondata
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none

  !----------------- mesh & domain -----------------
  integer, parameter :: nx = 256, ny = 256
  integer, parameter :: maxThreads = 24

  real(wp), parameter :: Lx = 2.0_wp
  real(wp), parameter :: Ly = 2.0_wp
  real(wp), parameter :: dx = Lx / real(nx,wp)
  real(wp), parameter :: dy = Ly / real(ny,wp)

  !----------------- physical parameters -----------------
  real(wp), parameter :: pi      = acos(-1.0_wp)
  real(wp), parameter :: Pe      = 100.0_wp
  real(wp), parameter :: u_phys  = 0.1_wp        
  real(wp), parameter :: cs2     = 1.0_wp/3.0_wp  ! c_s^2
  real(wp), parameter :: tau_phi = 0.56_wp

  ! LB scalar diffusivity: K_LB = cs2*(tau_phi - 0.5)
  real(wp), parameter :: alphaLB = cs2*(tau_phi - 0.5_wp)

  ! physical diffusivity: alpha_phys = 1/Pe
  real(wp), parameter :: alpha_phys = 1.0_wp/Pe

  ! time step mapping:
  ! alpha_phys = alphaLB * dx^2 / dt_phys  ⇒  dt_phys = alphaLB * Pe * dx^2
  real(wp), parameter :: dt_phys = alphaLB * Pe * dx*dx

  ! LB velocity corresponding to u_phys (used in B = uLB * rho)
  real(wp), parameter :: uLB  = u_phys * dt_phys / dx
  real(wp), parameter :: uLBx = uLB
  real(wp), parameter :: uLBy = uLB

  ! growth rate in analytic solution
  real(wp), parameter :: lambda = 1.0_wp - 2.0_wp*pi*pi/Pe

  ! user-specified physical end time
  real(wp), parameter :: tEnd_phys = 1.0_wp

  !----------------- D2Q5 discrete velocities -----------------
  integer, parameter :: q = 5
  integer, parameter :: ex(0:4) = (/ 0,  1, -1,  0,  0 /)
  integer, parameter :: ey(0:4) = (/ 0,  0,  0,  1, -1 /)
  real(wp), parameter :: w (0:4) = (/ 1.0_wp/3.0_wp, &
                                     1.0_wp/6.0_wp, 1.0_wp/6.0_wp, &
                                     1.0_wp/6.0_wp, 1.0_wp/6.0_wp /)

  !----------------- MRT relaxation rates -----------------
  ! m0: conserved → s0 = 0
  ! m1,m2: first moments → s_q = 1/tau_phi
  ! m3,m4: higher moments
  real(wp), parameter :: s0 = 0.0_wp
  real(wp), parameter :: s_q = 1.0_wp/tau_phi
  real(wp), parameter :: s_e = 1.0_wp
  real(wp), parameter :: s_p = 1.0_wp

  !----------------- fields & coordinates -----------------
  real(wp), allocatable :: rho(:,:), rho_prev(:,:)
  real(wp), allocatable :: f(:,:,:), f_post(:,:,:)
  real(wp), allocatable :: S_prev(:,:)   ! S_LB^{n-1}(i,j)
  real(wp), allocatable :: x(:), y(:)    ! physical cell centers

  ! time & step control
  integer :: nSteps
  integer :: it
  real(wp) :: t_phys

contains

  pure function rho_exact_phys(xp, yp, tp) result(val)
    real(wp), intent(in) :: xp, yp, tp
    real(wp) :: val
    ! rho_exact(x,y,t) = exp[(1 - 2π^2/Pe) t] * sin[π (x+y)]
    val = exp( lambda * tp ) * sin( pi*(xp+yp) )
  end function rho_exact_phys

  pure function F_phys(xp, yp, tp) result(val)
    real(wp), intent(in) :: xp, yp, tp
    real(wp) :: val
    real(wp) :: phase
    ! F(x,y,t) = exp(λ t) * { π(u1+u2) cos[π(x+y)] + sin[π(x+y)] }
    phase = pi*(xp+yp)
    val = exp( lambda * tp ) * ( pi*(2.0_wp*u_phys)*cos(phase) + sin(phase) )
  end function F_phys

end module commondata

!======================================================================
program lbm_case1_example31
  use omp_lib
  use commondata
  implicit none
  real(wp) :: t0, t1
  real(wp) :: errL1, errL2

  call omp_set_num_threads(maxThreads)

  call initial()

  call cpu_time(t0)

  do it = 1, nSteps

    call collision_case1_mrt()   ! MRT + G_j + Scheme 1 source
    call apply_periodic_ghosts() ! x,y periodic ghosts
    call streaming()             ! D2Q5 streaming
    call macro_rho()             ! update rho & rho_prev

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
    real(wp) :: xc, yc

    allocate(rho(nx,ny), rho_prev(nx,ny))
    allocate(f(0:nx+1,0:ny+1,0:4))
    allocate(f_post(0:nx+1,0:ny+1,0:4))
    allocate(S_prev(nx,ny))
    allocate(x(nx), y(ny))

    ! cell-center physical coordinates
    do i = 1, nx
      x(i) = (real(i,wp)-0.5_wp) * dx
    end do
    do j = 1, ny
      y(j) = (real(j,wp)-0.5_wp) * dy
    end do

    ! initial physical time
    t_phys = 0.0_wp

    ! compute LB steps from user-specified physical time
    nSteps = int( tEnd_phys / dt_phys + 0.5_wp )

    write(*,'(A,ES16.8)') 'User-specified t_end_phys = ', tEnd_phys
    write(*,'(A,ES16.8)') 'dt_phys (per LB step)     = ', dt_phys
    write(*,'(A,I10)')    'Number of LB steps        = ', nSteps
    write(*,'(A,ES16.8)') 'Actual t_end* = nSteps*dt = ', &
                          real(nSteps,wp)*dt_phys
    write(*,*) '---------------------------------------------'

    ! initial scalar field: exact solution at t=0
!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc)
    do j = 1, ny
      do i = 1, nx
        yc = y(j)
        xc = x(i)
        rho(i,j) = rho_exact_phys(xc,yc,t_phys)
      end do
    end do
!$omp end parallel do

    ! rho_prev initialized as rho^0
    rho_prev = rho

    ! initial source term S_prev = F_LB^0 = dt_phys * F_phys(x,y,t=0)
!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc)
    do j = 1, ny
      do i = 1, nx
      yc = y(j)
        xc = x(i)
        S_prev(i,j) = dt_phys * F_phys(xc,yc,t_phys)
      end do
    end do
!$omp end parallel do

    ! initialize distribution function to equilibrium f^eq
!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc)
    do j = 1, ny
      do i = 1, nx
      yc = y(j)
        xc = x(i)
        f(i,j,0) = w(0)*rho(i,j) * ( 1.0_wp + (ex(0)*uLBx + ey(0)*uLBy)/cs2 )
        f(i,j,1) = w(1)*rho(i,j) * ( 1.0_wp + (ex(1)*uLBx + ey(1)*uLBy)/cs2 )
        f(i,j,2) = w(2)*rho(i,j) * ( 1.0_wp + (ex(2)*uLBx + ey(2)*uLBy)/cs2 )
        f(i,j,3) = w(3)*rho(i,j) * ( 1.0_wp + (ex(3)*uLBx + ey(3)*uLBy)/cs2 )
        f(i,j,4) = w(4)*rho(i,j) * ( 1.0_wp + (ex(4)*uLBx + ey(4)*uLBy)/cs2 )
      end do
    end do
!$omp end parallel do

    ! initialize ghosts as periodic
    call apply_periodic_ghosts()

  end subroutine initial

  !---------------------------------
  subroutine collision_case1_mrt()
    implicit none
    integer :: i,j,a
    real(wp) :: m0,m1,m2,m3,m4
    real(wp) :: me0,me1,me2,me3,me4
    real(wp) :: ms0,ms1,ms2,ms3,ms4
    real(wp) :: Bx_curr,By_curr,Bx_prev,By_prev
    real(wp) :: dBx_dtLB, dBy_dtLB
    real(wp) :: M1Gx, M1Gy, Gj
    real(wp) :: xc,yc
    real(wp) :: S_curr, S_eff

!$omp parallel do collapse(2) default(shared) &
!$omp& private(i,j,a,m0,m1,m2,m3,m4, &
!$omp& me0,me1,me2,me3,me4, &
!$omp& ms0,ms1,ms2,ms3,ms4, &
!$omp& Bx_curr,By_curr,Bx_prev,By_prev, &
!$omp& dBx_dtLB,dBy_dtLB,M1Gx,M1Gy,Gj, &
!$omp& xc,yc,S_curr,S_eff)
    do j = 1, ny
      do i = 1, nx
      yc = y(j)
        xc = x(i)

        ! 1) moments m = M f
        m0 =  f(i,j,0) + f(i,j,1) + f(i,j,2) + f(i,j,3) + f(i,j,4)
        m1 =               f(i,j,1) - f(i,j,2)
        m2 =                                  f(i,j,3) - f(i,j,4)
        m3 = 4.0_wp*f(i,j,0) - ( f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4) )
        m4 =               f(i,j,1) + f(i,j,2) - f(i,j,3) - f(i,j,4)

        ! 2) equilibrium moments meq (Case 1: B = uLB * rho)
        Bx_curr = uLBx * rho(i,j)
        By_curr = uLBy * rho(i,j)

        me0 = 0.0_wp
        me1 = Bx_curr
        me2 = By_curr
        me3 = (2.0_wp/3.0_wp) * rho(i,j)
        me4 = 0.0_wp

        ! 3) MRT relaxation: m* = m - S (m - meq)
        ms0 = m0 - s0*(m0 - me0)
        ms1 = m1 - s_q*(m1 - me1)
        ms2 = m2 - s_q*(m2 - me2)
        ms3 = m3 - s_e*(m3 - me3)
        ms4 = m4 - s_p*(m4 - me4)

        ! 4) back to velocity space: f_coll = M^{-1} m*
        f_post(i,j,0) = 0.2_wp*ms0 + 0.2_wp*ms3
        f_post(i,j,1) = 0.2_wp*ms0 + 0.5_wp*ms1 - 0.05_wp*ms3 + 0.25_wp*ms4
        f_post(i,j,2) = 0.2_wp*ms0 - 0.5_wp*ms1 - 0.05_wp*ms3 + 0.25_wp*ms4
        f_post(i,j,3) = 0.2_wp*ms0 + 0.5_wp*ms2 - 0.05_wp*ms3 - 0.25_wp*ms4
        f_post(i,j,4) = 0.2_wp*ms0 - 0.5_wp*ms2 - 0.05_wp*ms3 - 0.25_wp*ms4

        ! 5) Case 1: M1G = (I - S1/2) * ∂t B,  ∂tB ≈ B^n - B^{n-1}
        Bx_prev  = uLBx * rho_prev(i,j)
        By_prev  = uLBy * rho_prev(i,j)
        dBx_dtLB = Bx_curr - Bx_prev
        dBy_dtLB = By_curr - By_prev

        M1Gx = (1.0_wp - 0.5_wp*s_q) * dBx_dtLB
        M1Gy = (1.0_wp - 0.5_wp*s_q) * dBy_dtLB

        ! 6) Scheme 1 source: Ω_F = 1.5 F^n - 0.5 F^{n-1}
        S_curr = dt_phys * F_phys(xc,yc,t_phys)
        S_eff  = 1.5_wp*S_curr - 0.5_wp*S_prev(i,j)
        S_prev(i,j) = S_curr

        ! 7) add source in velocity space: f_post += G_j + Ω_Fj
        do a = 0, 4
          Gj = w(a) * ( ex(a)*M1Gx + ey(a)*M1Gy ) / cs2
          f_post(i,j,a) = f_post(i,j,a) + Gj + w(a)*S_eff
        end do

      end do
    end do
!$omp end parallel do

  end subroutine collision_case1_mrt

  !---------------------------------
  subroutine apply_periodic_ghosts()
    implicit none
    integer :: i,j,a

    ! x periodic: i=0 ← nx,  i=nx+1 ← 1
!$omp parallel do collapse(2) default(shared) private(j,a)
    do j = 1, ny
      do a = 0, 4
        f_post(0    ,j,a) = f_post(nx,j,a)
        f_post(nx+1 ,j,a) = f_post(1 ,j,a)
      end do
    end do
!$omp end parallel do

    ! y periodic: j=0 ← ny,  j=ny+1 ← 1
!$omp parallel do collapse(2) default(shared) private(i,a)
    do i = 1, nx
      do a = 0, 4
        f_post(i,0    ,a) = f_post(i,ny,a)
        f_post(i,ny+1 ,a) = f_post(i,1 ,a)
      end do
    end do
!$omp end parallel do

  end subroutine apply_periodic_ghosts

  !---------------------------------
  subroutine streaming()
    implicit none
    integer :: i,j

!$omp parallel do collapse(2) default(shared) private(i,j)
    do j = 1, ny
      do i = 1, nx
        f(i,j,0) = f_post(i  ,j  ,0)
        f(i,j,1) = f_post(i-1,j  ,1)   ! E ← from left
        f(i,j,2) = f_post(i+1,j  ,2)   ! W ← from right
        f(i,j,3) = f_post(i  ,j-1,3)   ! N ← from bottom
        f(i,j,4) = f_post(i  ,j+1,4)   ! S ← from top
      end do
    end do
!$omp end parallel do

  end subroutine streaming

  !---------------------------------
  subroutine macro_rho()
    implicit none
    integer :: i,j

!$omp parallel do collapse(2) default(shared) private(i,j)
    do j = 1, ny
      do i = 1, nx
        rho_prev(i,j) = rho(i,j)
        rho(i,j) = f(i,j,0) + f(i,j,1) + f(i,j,2) + f(i,j,3) + f(i,j,4)
      end do
    end do
!$omp end parallel do

  end subroutine macro_rho

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
      yc = y(j)
        xc = x(i)
        exact = rho_exact_phys(xc,yc,t_phys)
        num = num + (rho(i,j) - exact)**2
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
  use commondata
  implicit none
  real(wp), intent(out) :: err
  integer :: i,j
  real(wp) :: num, den, xc, yc, exact

  num = 0.0_wp
  den = 0.0_wp

!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc,exact) reduction(+:num,den)
  do j = 1, ny
    do i = 1, nx
      yc = y(j)
      xc    = x(i)
      exact = rho_exact_phys(xc,yc,t_phys)
      num   = num + abs( rho(i,j) - exact )
      den   = den + abs( exact )
    end do
  end do
!$omp end parallel do

  if (den > 0.0_wp) then
    err = num / den
  else
    err = 0.0_wp
  end if

end subroutine compute_L1_error


end program lbm_case1_example31
