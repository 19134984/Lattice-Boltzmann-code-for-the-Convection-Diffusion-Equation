!======================================================================
!  D2Q9 MRT-LBM for Example 3.2 (nonlinear heat conduction, Dirichlet BC)
!
!  PDE (Phys. Rev. E 97, 043310, Eq. (58)):
!    u_t = alpha_phys * ( ∂_xx (u^δ) + ∂_yy (u^δ) ) + u - u^δ
!        = alpha_phys * ∇^2(u^δ) + S(u),
!    where S(u) = u - u^δ,  δ > 1.
!
!  Analytic solution (Eq. (59)):
!    u(x,y,t) = [ 1/2 + 1/2 * tanh( (δ-1)/(2δ√(2α)) * (x + y + √(2α)t) ) ]^{1/(1-δ)}
!             = inner^{ -1/(δ-1) }.
!
!  PRE-2020 MRT-LBM model (Case 1) with D2Q9:
!    - Case 1: beta = 1,  D = u^δ I, B = 0, C = 0, M1F = 0, γ = 0.
!    - Equilibrium distribution f_eq:
!        Eq. (77a), B = 0, isotropic D:
!        f_eq_j = w_j [ u + ((D - u) (|c_j|^2 - d c_s^2)) / (2 c_s^2) ],
!        with d = 2, c_s^2 = 1/3,  |c_j|^2 = ex_j^2 + ey_j^2.
!      This is algebraically equivalent to J. Sci. Comput. Eq. (2.4)/(3.23).
!
!  MRT part (moment space):
!    - Moment basis M given by J. Sci. Comput. Eq. (2.3).
!    - Equilibrium moments (Case 1, B = 0) [Eq. (2.11)]:
!        m0^eq = u
!        m1^eq = 2 D - 4 u
!        m2^eq = 3 u - 2 D
!        m3^eq = m4^eq = m5^eq = m6^eq = m7^eq = m8^eq = 0
!
!    - Relaxation matrix S = diag(s0,...,s8), with
!        s0 = 0  (conserved scalar)
!        s3 = s5 = 1/τ_phi  (diffusive modes, control alpha_LB)
!        other non-conserved modes set to O(1) (here = 1).
!
!    - Diffusion coefficient (Eq. (2.18) / (42)-(43)):
!        alpha_LB = c_s^2 (1/s3 - 1/2) Δt_LB
!                 = c_s^2 (τ_phi - 1/2)   (Δt_LB = 1).
!
!  Physical–lattice mapping (保持你原来的思路)：
!    alpha_phys = alpha_LB * dx^2 / dt_phys
!      ⇒ dt_phys = alpha_LB * dx^2 / alpha_phys
!
!  Boundary conditions:
!    - Dirichlet BC from analytic solution at all four sides.
!    - Implemented by half-way anti-bounce-back (ABB) for D2Q9:
!        e.g. left wall (x=0), for incoming directions 1,5,8:
!          f_1(ghost) = - f_3(fluid) + 2 f_1^eq(u_w),
!          f_5(ghost) = - f_7(fluid) + 2 f_5^eq(u_w),
!          f_8(ghost) = - f_6(fluid) + 2 f_8^eq(u_w),
!        similarly for right/bottom/top.
!
!  Time integration of source term S(u) = u - u^δ:
!    - Scheme 1 (Guo PRE 2018, Example 3.1):
!        S_LB^n   = dt_phys * S_phys(u^n),
!        S_eff^n  = 1.5 S_LB^n - 0.5 S_LB^{n-1},
!        f_j ← f_j + w_j * S_eff^n.
!
!======================================================================
module commondata_nl_d2q9
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none

  !----------------- mesh & domain -----------------
  integer, parameter :: nx = 100, ny = 100
  integer, parameter :: maxThreads = 24

  real(wp), parameter :: Lx = 1.0_wp
  real(wp), parameter :: Ly = 1.0_wp
  real(wp), parameter :: dx = Lx / real(nx,wp)
  real(wp), parameter :: dy = Ly / real(ny,wp)

  !----------------- physical & model parameters -----------------
  real(wp), parameter :: pi         = acos(-1.0_wp)
  real(wp), parameter :: alpha_phys = 0.01_wp       ! physical diffusivity α
  real(wp), parameter :: delta_n    = 1.2_wp        ! nonlinearity exponent δ

  real(wp), parameter :: cs2     = 1.0_wp/3.0_wp    ! c_s^2 for D2Q9
  real(wp), parameter :: tau_phi = 0.62_wp          ! diffusion "relaxation time"

  ! LB diffusion coefficient (Eq. α_LB = c_s^2 (τ_phi - 1/2))
  real(wp), parameter :: alphaLB = cs2 * (tau_phi - 0.5_wp)

  ! physical time step mapping (保持原来的写法):
  !   alpha_phys = alphaLB * dx^2 / dt_phys
  !   ⇒ dt_phys = alphaLB * dx^2 / alpha_phys
  real(wp), parameter :: dt_phys   = alphaLB * dx*dx / alpha_phys
  real(wp), parameter :: tEnd_phys = 1.0_wp

  !----------------- D2Q9 discrete velocities -----------------
  integer, parameter :: q = 9
  integer, parameter :: ex(0:8) = (/ 0,  1,  0, -1,  0,  1, -1, -1,  1 /)
  integer, parameter :: ey(0:8) = (/ 0,  0,  1,  0, -1,  1,  1, -1, -1 /)
  real(wp), parameter :: w (0:8) = (/ 4.0_wp/9.0_wp, &
                                     1.0_wp/9.0_wp, 1.0_wp/9.0_wp, 1.0_wp/9.0_wp, 1.0_wp/9.0_wp, &
                                     1.0_wp/36.0_wp,1.0_wp/36.0_wp,1.0_wp/36.0_wp,1.0_wp/36.0_wp /)

  !----------------- MRT relaxation rates (S = diag(s0,...,s8)) -----------------
  !  Moment ordering (J. Sci. Comput. Eq. (2.3)):
  !    m0: scalar u     (conserved)
  !    m1,m2: energy-like modes
  !    m3,m5: "flux" modes (control diffusion) → s3,s5 = 1/τ_phi
  !    others: higher-order modes (set to O(1), here 1.0)
  real(wp), parameter :: s0 = 0.0_wp
  real(wp), parameter :: s1 = 1.0_wp
  real(wp), parameter :: s2 = 1.0_wp
  real(wp), parameter :: s3 = 1.0_wp/tau_phi   ! controls α_LB (together with s5)
  real(wp), parameter :: s4 = 1.0_wp
  real(wp), parameter :: s5 = 1.0_wp/tau_phi   ! controls α_LB
  real(wp), parameter :: s6 = 1.0_wp
  real(wp), parameter :: s7 = 1.0_wp
  real(wp), parameter :: s8 = 1.0_wp

  !----------------- fields & coordinates -----------------
  real(wp), allocatable :: phi(:,:)                 ! scalar u
  real(wp), allocatable :: f(:,:,:), f_post(:,:,:)  ! distributions
  real(wp), allocatable :: S_prev(:,:)              ! S_LB^{n-1}(i,j)
  real(wp), allocatable :: x(:), y(:)               ! physical cell centers

  ! time & step control
  integer :: nSteps
  integer :: it
  real(wp) :: t_phys

contains

  !------------------------------------------------------------------
  pure function u_exact_phys(xp, yp, tp) result(val)
    ! Analytic solution (Eq. (59) in PRE 97, 043310):
    !   u(x,y,t) = [ 1/2 + 1/2 tanh( (δ-1)/(2δ√(2α)) * (x + y + √(2α)t) ) ]^{1/(1-δ)}
    real(wp), intent(in) :: xp, yp, tp
    real(wp) :: val
    real(wp) :: theta, inner

    theta = (delta_n - 1.0_wp) / ( 2.0_wp * delta_n * sqrt( 2.0_wp * alpha_phys ) ) * &
            ( xp + yp + sqrt( 2.0_wp * alpha_phys ) * tp )

    inner = 0.5_wp + 0.5_wp * tanh(theta)

    ! exponent 1/(1-δ) = -1/(δ-1)
    val = inner ** ( 1.0_wp / ( 1.0_wp - delta_n ) )
  end function u_exact_phys

  !------------------------------------------------------------------
  pure function S_phys_from_u(u) result(val)
    ! Reaction term S(u) = u - u^δ (PDE RHS, Eq. (58))
    real(wp), intent(in) :: u
    real(wp) :: val
    val = u - u**delta_n
  end function S_phys_from_u

  !------------------------------------------------------------------
  pure function feq_dir(u, D_u, a) result(val)
    ! Equilibrium distribution f_eq_j, Case 1, B=0 (PRE-2020 Eq. (77a) special case):
    !
    !   f_eq_j = w_j [ u + ((D - u)(|c_j|^2 - d c_s^2)) / (2 c_s^2) ],
    !
    ! with d = 2, c_s^2 = 1/3, |c_j|^2 = ex_j^2 + ey_j^2.
    ! This is algebraically equivalent to Eq. (2.4)/(3.23) in J. Sci. Comput.
    real(wp), intent(in) :: u, D_u
    integer , intent(in) :: a
    real(wp) :: val
    real(wp) :: c2, d

    d  = 2.0_wp
    c2 = real(ex(a)*ex(a) + ey(a)*ey(a), wp)

    val = w(a) * ( u + ( (D_u - u) * ( c2 - d*cs2 ) ) / ( 2.0_wp * cs2 ) )
  end function feq_dir

end module commondata_nl_d2q9

!======================================================================
program lbm_case1_example32_d2q9
  use omp_lib
  use commondata_nl_d2q9
  implicit none
  real(wp) :: t0, t1
  real(wp) :: errL2, errL1

  call omp_set_num_threads(maxThreads)

  call initial()

  call cpu_time(t0)

  do it = 1, nSteps

    call collision_case1_mrt()   ! MRT collision + Scheme 1 source
    call apply_dirichlet_bc()    ! Dirichlet BC via half-way ABB
    call streaming()             ! D2Q9 streaming
    call macro_phi()             ! update phi

    t_phys = t_phys + dt_phys

    ! 可选：每步监控误差或极值
    !call compute_L2_error(errL2)
    !write(*,'("it=",I6,"  L2_err=",ES12.4)') it, errL2

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
    real(wp) :: phi_loc, D_loc

    allocate(phi(nx,ny))
    allocate(f(0:nx+1,0:ny+1,0:8))
    allocate(f_post(0:nx+1,0:ny+1,0:8))
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

    ! compute number of LB steps from user-specified physical time
    nSteps = int( tEnd_phys / dt_phys + 0.5_wp )

    write(*,'(A,ES16.8)') 'User-specified t_end_phys = ', tEnd_phys
    write(*,'(A,ES16.8)') 'dt_phys (per LB step)     = ', dt_phys
    write(*,'(A,I10)')    'Number of LB steps        = ', nSteps
    write(*,'(A,ES16.8)') 'Actual t_end* = nSteps*dt = ', &
                          real(nSteps,wp)*dt_phys
    write(*,*) '---------------------------------------------'

    ! initial scalar field: exact solution at t=0
!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc,phi_loc,D_loc)
    do j = 1, ny
      do i = 1, nx
      yc = y(j)
        xc = x(i)
        phi_loc    = u_exact_phys(xc,yc,t_phys)
        phi(i,j)   = phi_loc
        D_loc      = phi_loc**delta_n

        ! initialize distributions to equilibrium f^eq (Case 1, Eq. (77a))
        f(i,j,0) = feq_dir( phi_loc, D_loc, 0 )
        f(i,j,1) = feq_dir( phi_loc, D_loc, 1 )
        f(i,j,2) = feq_dir( phi_loc, D_loc, 2 )
        f(i,j,3) = feq_dir( phi_loc, D_loc, 3 )
        f(i,j,4) = feq_dir( phi_loc, D_loc, 4 )
        f(i,j,5) = feq_dir( phi_loc, D_loc, 5 )
        f(i,j,6) = feq_dir( phi_loc, D_loc, 6 )
        f(i,j,7) = feq_dir( phi_loc, D_loc, 7 )
        f(i,j,8) = feq_dir( phi_loc, D_loc, 8 )
      end do
    end do
!$omp end parallel do

    ! initial source term S_prev = S_LB^0 = dt_phys * S_phys(u^0)
!$omp parallel do collapse(2) default(shared) private(i,j)
    do j = 1, ny
      do i = 1, nx
        S_prev(i,j) = dt_phys * S_phys_from_u( phi(i,j) )
      end do
    end do
!$omp end parallel do

  end subroutine initial

  !---------------------------------
  subroutine collision_case1_mrt()
    ! MRT collision in moment space + Scheme-1 source term.
    ! Moment basis M and inverse M^{-1} correspond to Eq. (2.3) in J. Sci. Comput.
    implicit none
    integer :: i,j
    integer :: a
    real(wp) :: f0,f1,f2,f3,f4,f5,f6,f7,f8
    real(wp) :: m0,m1,m2,m3,m4,m5,m6,m7,m8
    real(wp) :: me0,me1,me2,me3,me4,me5,me6,me7,me8
    real(wp) :: ms0,ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8
    real(wp) :: phi_loc, D_loc
    real(wp) :: S_phys_loc, S_curr, S_eff
    real(wp), parameter :: phi_floor = 1.0e-12_wp
    real(wp) :: phi_min_loc, phi_max_loc
    real(wp) :: phi_min_glb, phi_max_glb

    phi_min_glb =  1.0e30_wp
    phi_max_glb = -1.0e30_wp

!$omp parallel do collapse(2) default(shared) private(i,j,a, &
!$omp& f0,f1,f2,f3,f4,f5,f6,f7,f8, &
!$omp& m0,m1,m2,m3,m4,m5,m6,m7,m8, &
!$omp& me0,me1,me2,me3,me4,me5,me6,me7,me8, &
!$omp& ms0,ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8, &
!$omp& phi_loc,D_loc,S_phys_loc,S_curr,S_eff, &
!$omp& phi_min_loc,phi_max_loc) reduction(min:phi_min_glb) reduction(max:phi_max_glb)
    do j = 1, ny
      do i = 1, nx

        ! local scalar and nonlinearity
        phi_loc = phi(i,j)
        !if (phi_loc < phi_floor) phi_loc = phi_floor
        D_loc   = phi_loc**delta_n

        ! Scheme 1 source term (Guo PRE 2018, Example 3.1)
        S_phys_loc = S_phys_from_u( phi_loc )
        S_curr     = dt_phys * S_phys_loc
        S_eff      = 1.5_wp*S_curr - 0.5_wp*S_prev(i,j)
        S_prev(i,j)= S_curr

        ! store local min/max BEFORE clamping in macroscopic field (for监控)
        !phi_min_loc = phi_loc
        !phi_max_loc = phi_loc
        !if (phi_min_loc < phi_min_glb) phi_min_glb = phi_min_loc
        !if (phi_max_loc > phi_max_glb) phi_max_glb = phi_max_loc

        !---------------- moments m = M f  (J. Sci. Comput. Eq. (2.3)) -----------
        f0 = f(i,j,0);  f1 = f(i,j,1);  f2 = f(i,j,2)
        f3 = f(i,j,3);  f4 = f(i,j,4);  f5 = f(i,j,5)
        f6 = f(i,j,6);  f7 = f(i,j,7);  f8 = f(i,j,8)

        m0 = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8

        m1 = -4.0_wp*f0 - f1 - f2 - f3 - f4 + 2.0_wp*(f5+f6+f7+f8)
        m2 =  4.0_wp*f0 - 2.0_wp*(f1+f2+f3+f4) + (f5+f6+f7+f8)

        m3 =  f1 - f3 + f5 - f6 - f7 + f8
        m4 = -2.0_wp*f1 + 2.0_wp*f3 + f5 - f6 - f7 + f8

        m5 =  f2 - f4 + f5 + f6 - f7 - f8
        m6 = -2.0_wp*f2 + 2.0_wp*f4 + f5 + f6 - f7 - f8

        m7 =  f1 - f2 + f3 - f4
        m8 =  f5 - f6 + f7 - f8

        !---------------- equilibrium moments m^eq (Case 1, Eq. (2.11)) ----------
        me0 = phi_loc
        me1 = 2.0_wp*D_loc - 4.0_wp*phi_loc
        me2 = 3.0_wp*phi_loc - 2.0_wp*D_loc
        me3 = 0.0_wp
        me4 = 0.0_wp
        me5 = 0.0_wp
        me6 = 0.0_wp
        me7 = 0.0_wp
        me8 = 0.0_wp

        !---------------- MRT relaxation: m* = m - S (m - m^eq) ------------------
        ms0 = m0                         ! s0 = 0 → conserved
        ms1 = m1 - s1*(m1 - me1)
        ms2 = m2 - s2*(m2 - me2)
        ms3 = m3 - s3*(m3 - me3)        ! s3 = 1/τ_phi (controls α_LB)
        ms4 = m4 - s4*(m4 - me4)
        ms5 = m5 - s5*(m5 - me5)        ! s5 = 1/τ_phi (controls α_LB)
        ms6 = m6 - s6*(m6 - me6)
        ms7 = m7 - s7*(m7 - me7)
        ms8 = m8 - s8*(m8 - me8)

        !---------------- back to velocity space: f_post = M^{-1} m* -------------
        !  M^{-1} rows (from symbolic inversion), equivalent to J. Sci. Comput.:
        !
        !  f0 =  1/9 m0 - 1/9 m1 + 1/9 m2
        !  f1 =  1/9 m0 - 1/36 m1 - 1/18 m2 + 1/6 m3 - 1/6 m4 + 0 m5 + 0 m6 + 1/4 m7 + 0 m8
        !  f2 =  1/9 m0 - 1/36 m1 - 1/18 m2 + 0 m3   + 0 m4   + 1/6 m5 - 1/6 m6 - 1/4 m7 + 0 m8
        !  f3 =  1/9 m0 - 1/36 m1 - 1/18 m2 - 1/6 m3 + 1/6 m4 + 0 m5 + 0 m6 + 1/4 m7 + 0 m8
        !  f4 =  1/9 m0 - 1/36 m1 - 1/18 m2 + 0 m3   + 0 m4   - 1/6 m5 + 1/6 m6 - 1/4 m7 + 0 m8
        !  f5 =  1/9 m0 + 1/18 m1 + 1/36 m2 + 1/6 m3 + 1/12 m4 + 1/6 m5 + 1/12 m6 + 0 m7 + 1/4 m8
        !  f6 =  1/9 m0 + 1/18 m1 + 1/36 m2 - 1/6 m3 - 1/12 m4 + 1/6 m5 + 1/12 m6 + 0 m7 - 1/4 m8
        !  f7 =  1/9 m0 + 1/18 m1 + 1/36 m2 - 1/6 m3 - 1/12 m4 - 1/6 m5 - 1/12 m6 + 0 m7 + 1/4 m8
        !  f8 =  1/9 m0 + 1/18 m1 + 1/36 m2 + 1/6 m3 + 1/12 m4 - 1/6 m5 - 1/12 m6 + 0 m7 - 1/4 m8

        f_post(i,j,0) = (1.0_wp/9.0_wp)*ms0 + (-1.0_wp/9.0_wp)*ms1 + (1.0_wp/9.0_wp)*ms2

        f_post(i,j,1) = (1.0_wp/9.0_wp)*ms0 + (-1.0_wp/36.0_wp)*ms1 + (-1.0_wp/18.0_wp)*ms2 + &
                        (1.0_wp/6.0_wp)*ms3 + (-1.0_wp/6.0_wp)*ms4 + &
                        0.0_wp*ms5 + 0.0_wp*ms6 + (1.0_wp/4.0_wp)*ms7 + 0.0_wp*ms8

        f_post(i,j,2) = (1.0_wp/9.0_wp)*ms0 + (-1.0_wp/36.0_wp)*ms1 + (-1.0_wp/18.0_wp)*ms2 + &
                        0.0_wp*ms3 + 0.0_wp*ms4 + &
                        (1.0_wp/6.0_wp)*ms5 + (-1.0_wp/6.0_wp)*ms6 + (-1.0_wp/4.0_wp)*ms7 + 0.0_wp*ms8

        f_post(i,j,3) = (1.0_wp/9.0_wp)*ms0 + (-1.0_wp/36.0_wp)*ms1 + (-1.0_wp/18.0_wp)*ms2 + &
                        (-1.0_wp/6.0_wp)*ms3 + (1.0_wp/6.0_wp)*ms4 + &
                        0.0_wp*ms5 + 0.0_wp*ms6 + (1.0_wp/4.0_wp)*ms7 + 0.0_wp*ms8

        f_post(i,j,4) = (1.0_wp/9.0_wp)*ms0 + (-1.0_wp/36.0_wp)*ms1 + (-1.0_wp/18.0_wp)*ms2 + &
                        0.0_wp*ms3 + 0.0_wp*ms4 + &
                        (-1.0_wp/6.0_wp)*ms5 + (1.0_wp/6.0_wp)*ms6 + (-1.0_wp/4.0_wp)*ms7 + 0.0_wp*ms8

        f_post(i,j,5) = (1.0_wp/9.0_wp)*ms0 + (1.0_wp/18.0_wp)*ms1 + (1.0_wp/36.0_wp)*ms2 + &
                        (1.0_wp/6.0_wp)*ms3 + (1.0_wp/12.0_wp)*ms4 + &
                        (1.0_wp/6.0_wp)*ms5 + (1.0_wp/12.0_wp)*ms6 + 0.0_wp*ms7 + (1.0_wp/4.0_wp)*ms8

        f_post(i,j,6) = (1.0_wp/9.0_wp)*ms0 + (1.0_wp/18.0_wp)*ms1 + (1.0_wp/36.0_wp)*ms2 + &
                        (-1.0_wp/6.0_wp)*ms3 + (-1.0_wp/12.0_wp)*ms4 + &
                        (1.0_wp/6.0_wp)*ms5 + (1.0_wp/12.0_wp)*ms6 + 0.0_wp*ms7 + (-1.0_wp/4.0_wp)*ms8

        f_post(i,j,7) = (1.0_wp/9.0_wp)*ms0 + (1.0_wp/18.0_wp)*ms1 + (1.0_wp/36.0_wp)*ms2 + &
                        (-1.0_wp/6.0_wp)*ms3 + (-1.0_wp/12.0_wp)*ms4 + &
                        (-1.0_wp/6.0_wp)*ms5 + (-1.0_wp/12.0_wp)*ms6 + 0.0_wp*ms7 + (1.0_wp/4.0_wp)*ms8

        f_post(i,j,8) = (1.0_wp/9.0_wp)*ms0 + (1.0_wp/18.0_wp)*ms1 + (1.0_wp/36.0_wp)*ms2 + &
                        (1.0_wp/6.0_wp)*ms3 + (1.0_wp/12.0_wp)*ms4 + &
                        (-1.0_wp/6.0_wp)*ms5 + (-1.0_wp/12.0_wp)*ms6 + 0.0_wp*ms7 + (-1.0_wp/4.0_wp)*ms8

        !---------------- add Scheme-1 source term in velocity space --------------
        do a = 0, 8
          f_post(i,j,a) = f_post(i,j,a) + w(a)*S_eff
        end do

      end do
    end do
!$omp end parallel do

    ! 每一步输出一次全场的 phi_min, phi_max 方便监控
    !write(*,'("it=",I6,"  phi_min_loc≈",ES12.4,"  phi_max_loc≈",ES12.4)') &
         !it, phi_min_glb, phi_max_glb

  end subroutine collision_case1_mrt

  !---------------------------------
  subroutine apply_dirichlet_bc()
    ! Half-way anti-bounce-back (ABB) for Dirichlet BC using analytic solution.
    ! For each boundary:
    !   f_q(ghost) = - f_qbar(fluid) + 2 f_q^eq(u_w),
    ! where q and qbar are opposite directions.
    implicit none
    integer :: i,j
    real(wp) :: xw, yw, u_w, D_w
    real(wp) :: t_bc

    t_bc = t_phys + dt_phys   ! boundary value at time level n+1 (保持你现在的用法)

    !---------------- Left boundary (x = 0): unknown incoming q = 1,5,8 -----------
!$omp parallel do default(shared) private(j,xw,yw,u_w,D_w)
    do j = 1, ny
      xw = 0.0_wp
      yw = y(j)
      u_w = u_exact_phys(xw,yw,t_bc)
      D_w = u_w**delta_n

      f_post(0 ,j,1) = - f_post(1 ,j,3) + 2.0_wp * feq_dir(u_w, D_w, 1)
      f_post(0 ,j,5) = - f_post(1 ,j,7) + 2.0_wp * feq_dir(u_w, D_w, 5)
      f_post(0 ,j,8) = - f_post(1 ,j,6) + 2.0_wp * feq_dir(u_w, D_w, 8)
    end do
!$omp end parallel do

    !---------------- Right boundary (x = Lx): incoming q = 3,6,7 -----------------
!$omp parallel do default(shared) private(j,xw,yw,u_w,D_w)
    do j = 1, ny
      xw = Lx
      yw = y(j)
      u_w = u_exact_phys(xw,yw,t_bc)
      D_w = u_w**delta_n

      f_post(nx+1 ,j,3) = - f_post(nx ,j,1) + 2.0_wp * feq_dir(u_w, D_w, 3)
      f_post(nx+1 ,j,6) = - f_post(nx ,j,8) + 2.0_wp * feq_dir(u_w, D_w, 6)
      f_post(nx+1 ,j,7) = - f_post(nx ,j,5) + 2.0_wp * feq_dir(u_w, D_w, 7)
    end do
!$omp end parallel do

    !---------------- Bottom boundary (y = 0): incoming q = 2,5,6 -----------------
!$omp parallel do default(shared) private(i,xw,yw,u_w,D_w)
    do i = 1, nx
      xw = x(i)
      yw = 0.0_wp
      u_w = u_exact_phys(xw,yw,t_bc)
      D_w = u_w**delta_n

      f_post(i ,0,2) = - f_post(i ,1,4) + 2.0_wp * feq_dir(u_w, D_w, 2)
      f_post(i ,0,5) = - f_post(i ,1,7) + 2.0_wp * feq_dir(u_w, D_w, 5)
      f_post(i ,0,6) = - f_post(i ,1,8) + 2.0_wp * feq_dir(u_w, D_w, 6)
    end do
!$omp end parallel do

    !---------------- Top boundary (y = Ly): incoming q = 4,7,8 -------------------
!$omp parallel do default(shared) private(i,xw,yw,u_w,D_w)
    do i = 1, nx
      xw = x(i)
      yw = Ly
      u_w = u_exact_phys(xw,yw,t_bc)
      D_w = u_w**delta_n

      f_post(i ,ny+1,4) = - f_post(i ,ny,2) + 2.0_wp * feq_dir(u_w, D_w, 4)
      f_post(i ,ny+1,7) = - f_post(i ,ny,5) + 2.0_wp * feq_dir(u_w, D_w, 7)
      f_post(i ,ny+1,8) = - f_post(i ,ny,6) + 2.0_wp * feq_dir(u_w, D_w, 8)
    end do
!$omp end parallel do

 !====== 角点：只处理对角方向 ======
    ! 左下角 (0,0)，fluid (1,1)，方向 5 ↔ 7
    xw = 0.0_wp
    yw = 0.0_wp
    u_w = u_exact_phys(xw,yw,t_bc)
    !if (u_w < phi_floor) u_w = phi_floor
    D_w = u_w**delta_n
    f_post(0,0,5) = - f_post(1,1,7) + 2.0_wp * feq_dir(u_w, D_w, 5)

    ! 左上角 (0,ny+1)，fluid (1,ny)，方向 8 ↔ 6
    xw = 0.0_wp
    yw = Ly
    u_w = u_exact_phys(xw,yw,t_bc)
    !if (u_w < phi_floor) u_w = phi_floor
    D_w = u_w**delta_n
    f_post(0,ny+1,8) = - f_post(1,ny,6) + 2.0_wp * feq_dir(u_w, D_w, 8)

    ! 右下角 (nx+1,0)，fluid (nx,1)，方向 6 ↔ 8
    xw = Lx
    yw = 0.0_wp
    u_w = u_exact_phys(xw,yw,t_bc)
    !if (u_w < phi_floor) u_w = phi_floor
    D_w = u_w**delta_n
    f_post(nx+1,0,6) = - f_post(nx,1,8) + 2.0_wp * feq_dir(u_w, D_w, 6)

    ! 右上角 (nx+1,ny+1)，fluid (nx,ny)，方向 7 ↔ 5
    xw = Lx
    yw = Ly
    u_w = u_exact_phys(xw,yw,t_bc)
    !if (u_w < phi_floor) u_w = phi_floor
    D_w = u_w**delta_n
    f_post(nx+1,ny+1,7) = - f_post(nx,ny,5) + 2.0_wp * feq_dir(u_w, D_w, 7)

  end subroutine apply_dirichlet_bc


  !---------------------------------
  subroutine streaming()
    implicit none
    integer :: i,j

!$omp parallel do collapse(2) default(shared) private(i,j)
    do j = 1, ny
      do i = 1, nx
        ! streaming: f(i,j,a) at new time comes from neighbor along -c_a
        f(i,j,0) = f_post(i    , j    , 0)
        f(i,j,1) = f_post(i-1  , j    , 1)   ! from left
        f(i,j,2) = f_post(i    , j-1  , 2)   ! from bottom
        f(i,j,3) = f_post(i+1  , j    , 3)   ! from right
        f(i,j,4) = f_post(i    , j+1  , 4)   ! from top
        f(i,j,5) = f_post(i-1  , j-1  , 5)   ! from bottom-left
        f(i,j,6) = f_post(i+1  , j-1  , 6)   ! from bottom-right
        f(i,j,7) = f_post(i+1  , j+1  , 7)   ! from top-right
        f(i,j,8) = f_post(i-1  , j+1  , 8)   ! from top-left
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
        phi(i,j) = f(i,j,0) + f(i,j,1) + f(i,j,2) + f(i,j,3) + f(i,j,4) + &
                   f(i,j,5) + f(i,j,6) + f(i,j,7) + f(i,j,8)
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
      yc = y(j)
        xc = x(i)
        exact = u_exact_phys(xc,yc,t_phys)
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
        yc = y(j)
        xc = x(i)
        exact = u_exact_phys(xc,yc,t_phys)
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

end program lbm_case1_example32_d2q9
