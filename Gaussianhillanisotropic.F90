!======================================================================
!  D2Q5 MRT-LBM for Gaussian hill advection–diffusion (anisotropic K)
!
!  PDE (physical):
!    ∂_t φ + u_phys · ∇φ = ∇·( K_phys ∇φ ),
!
!    K_phys =
!      [ Kxx_phys   Kxy_phys ;
!        Kxy_phys   Kyy_phys ].
!
!  Time mapping:
!    tEnd_phys   = 1.0  (物理时间 1 s)
!    nSteps      = 200
!    dt_phys     = tEnd_phys / nSteps
!    每一步 LBM 对应 dt_phys 的物理时间。
!
!  LBM mapping (C = 0, Δt_LB = 1, β = 1):
!    K_phys = K_LB * dx^2 / dt_phys,
!    K_LB   = cs2 ( S1^{-1} - 0.5 I ),
!    S1     = [s11 s12; s12 s22]  (一阶矩的弛豫矩阵).
!
!  Convection:
!    B = φ u_LB,  u_LB = u_phys * dt_phys / dx,
!    C=0 scheme:
!      M1G = (I - S1/2) ∂_t B,  ∂_t B ≈ B^n - B^{n-1},
!      G_j = w_j ( c_j · M1G ) / cs2.
!
!  Analytic solution: anisotropic Gaussian hill with convection:
!    Σ_t = Σ_0 + 2 K_phys t,
!    Σ_0 = σ0^2 I,
!    r = [ x - ux_phys t ; y - uy_phys t ],
!    detΣ = Σ_xx Σ_yy - Σ_xy^2,
!    φ(x,y,t) = σ0^2 / sqrt(detΣ)
!               * exp( -0.5 * r^T Σ_t^{-1} r ).
!
!======================================================================
module commondata_gauss_adv_aniso_d2q5
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none

  !----------------- mesh & domain -----------------
  integer, parameter :: nx = 128, ny = 128
  real(wp), parameter :: Lx = 2.0_wp      ! x ∈ [-1, 1]
  real(wp), parameter :: Ly = 2.0_wp      ! y ∈ [-1, 1]
  real(wp), parameter :: dx = Lx / real(nx,wp)
  real(wp), parameter :: dy = Ly / real(ny,wp)

  !----------------- physical time control -----------------
  real(wp), parameter :: tEnd_phys = 1.0_wp   ! 1 s
  integer              :: nSteps

  !----------------- physical diffusion tensor K_phys -----------------
  ! 你可以在这里设置各向异性的扩散系数：
  !
  !
  ! 1) 轴向各向异性: K * 1e4 = [1 0; 0 4]
  ! real(wp), parameter :: Kxx_phys = 1.0e-4_wp
  ! real(wp), parameter :: Kxy_phys = 0.0e-4_wp
  ! real(wp), parameter :: Kyy_phys = 4.0e-4_wp
  !
  ! 2) 完全各向异性: K * 1e4 = [1 1; 1 4]
  real(wp), parameter :: Kxx_phys = 1.0e-4_wp
  real(wp), parameter :: Kxy_phys = 1.0e-4_wp
  real(wp), parameter :: Kyy_phys = 4.0e-4_wp

  !----------------- physical convection velocity -----------------
  real(wp), parameter :: ux_phys = 0.05_wp
  real(wp), parameter :: uy_phys = 0.05_wp

  !----------------- Gaussian initial width -----------------
  real(wp), parameter :: sigma0 = 0.05_wp

  !----------------- D2Q5 discrete velocities -----------------
  integer, parameter :: q = 5
  integer, parameter :: ex(0:4) = (/ 0,  1,  0, -1,  0 /)
  integer, parameter :: ey(0:4) = (/ 0,  0,  1,  0, -1 /)
  real(wp), parameter :: w (0:4) = (/ 1.0_wp/3.0_wp, &
                                     1.0_wp/6.0_wp, 1.0_wp/6.0_wp, &
                                     1.0_wp/6.0_wp, 1.0_wp/6.0_wp /)

  real(wp), parameter :: cs2 = 1.0_wp/3.0_wp

    !----------------- 以 λ_max 和 KLB_target 决定 dt_phys -----------------
  ! 1) 物理扩散张量的最大特征值 λ_max_phys
  real(wp), parameter :: traceK          = Kxx_phys + Kyy_phys
  real(wp), parameter :: deltaK          = Kxx_phys - Kyy_phys
  real(wp), parameter :: discK           = sqrt( deltaK*deltaK + 4.0_wp*Kxy_phys*Kxy_phys )
  real(wp), parameter :: lambda_max_phys = 0.5_wp*( traceK + discK )

  ! 2) 目标格子扩散数 (类似各向同性时的 D_LB)，其实就是cs2*(tau_phi - 0.5_wp)
  real(wp), parameter :: KLB_target = 0.10_wp

  ! 3) 时间步长: dt_phys = KLB_target * dx^2 / lambda_max_phys
  real(wp), parameter :: dt_phys = KLB_target * dx*dx / lambda_max_phys


  !----------------- LB diffusion tensor from K_phys -----------------
  ! K_LB = K_phys * dt_phys / dx^2  (假设 dx = dy)
  real(wp), parameter :: Kxx_LB = Kxx_phys * dt_phys / (dx*dx)
  real(wp), parameter :: Kxy_LB = Kxy_phys * dt_phys / (dx*dx)
  real(wp), parameter :: Kyy_LB = Kyy_phys * dt_phys / (dx*dx)

  !----------------- MRT S1 matrix from K_LB ------------------
  ! S1^{-1} = K_LB / cs2 + 0.5 I
  real(wp), parameter :: a11 = Kxx_LB/cs2 + 0.5_wp
  real(wp), parameter :: a22 = Kyy_LB/cs2 + 0.5_wp
  real(wp), parameter :: a12 = Kxy_LB/cs2
  real(wp), parameter :: detA = a11*a22 - a12*a12

  real(wp), parameter :: s11 =  a22/detA
  real(wp), parameter :: s22 =  a11/detA
  real(wp), parameter :: s12 = -a12/detA

  ! 其它模式的弛豫率：取 O(1) 即可
  real(wp), parameter :: s0 = 0.0_wp
  real(wp), parameter :: s3 = 1.0_wp
  real(wp), parameter :: s4 = 1.0_wp

  !----------------- LB convection velocity (lattice units) ----------
  ! u_LB = u_phys * dt_phys / dx
  real(wp), parameter :: ux_LB = ux_phys * dt_phys / dx
  real(wp), parameter :: uy_LB = uy_phys * dt_phys / dy

  !----------------- fields & coordinates -----------------
  real(wp), allocatable :: phi(:,:)                 ! scalar φ
  real(wp), allocatable :: f(:,:,:), f_post(:,:,:)  ! distributions
  real(wp), allocatable :: Bx_prev(:,:), By_prev(:,:) ! 存储上一时间步的 B

  real(wp), allocatable :: x(:), y(:)               ! cell centers

  real(wp) :: t_phys
  integer :: it

contains

  !------------------------------------------------------------------
  pure function phi_exact_phys(xp, yp, tp) result(val)
    ! Anisotropic Gaussian hill with convection
    real(wp), intent(in) :: xp, yp, tp
    real(wp) :: val
    real(wp) :: rx, ry
    real(wp) :: s_xx, s_xy, s_yy
    real(wp) :: detS, quad

    ! 移动中心
    rx = xp - ux_phys*tp
    ry = yp - uy_phys*tp

    ! Σ_t = Σ_0 + 2 K_phys t,  Σ_0 = σ0^2 I
    s_xx = sigma0*sigma0 + 2.0_wp*Kxx_phys*tp
    s_yy = sigma0*sigma0 + 2.0_wp*Kyy_phys*tp
    s_xy = 2.0_wp*Kxy_phys*tp

    detS = s_xx*s_yy - s_xy*s_xy

    quad = ( s_yy*rx*rx - 2.0_wp*s_xy*rx*ry + s_xx*ry*ry ) / detS

    ! 归一化使 t=0 时 φ(0,0,0) = 1
    val = ( sigma0*sigma0 / sqrt(detS) ) * exp( -0.5_wp*quad )
  end function phi_exact_phys

end module commondata_gauss_adv_aniso_d2q5

!======================================================================
program lbm_gauss_adv_aniso_d2q5
  use omp_lib
  use commondata_gauss_adv_aniso_d2q5
  implicit none
  real(wp) :: t0, t1
  real(wp) :: errL2, errL1

  call omp_set_num_threads(omp_get_max_threads())

  call initial()

  call cpu_time(t0)

  do it = 1, nSteps

    call collision_mrt_adv_aniso()  ! MRT + anisotropic diffusion + convection
    call streaming()                ! periodic
    call macro_phi()                ! update phi

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

    ! physical coordinates: x ∈ [-Lx/2, Lx/2], y ∈ [-Ly/2, Ly/2]
    do i = 1, nx
      x(i) = -0.5_wp*Lx + (real(i,wp)-0.5_wp)*dx
    end do
    do j = 1, ny
      y(j) = -0.5_wp*Ly + (real(j,wp)-0.5_wp)*dy
    end do

    t_phys = 0.0_wp
    ! 根据 tEnd_phys 和 dt_phys 计算步数
    nSteps = int( tEnd_phys / dt_phys + 0.5_wp )

    write(*,'(A,ES16.8)') 'tEnd_phys  = ', tEnd_phys
    write(*,'(A,ES16.8)') 'dt_phys    = ', dt_phys
    write(*,'(A,I10)')    'nSteps     = ', nSteps
    write(*,'(A,ES16.8)') 'dx         = ', dx
    write(*,'(A,ES16.8)') 'Kxx_phys   = ', Kxx_phys
    write(*,'(A,ES16.8)') 'Kxy_phys   = ', Kxy_phys
    write(*,'(A,ES16.8)') 'Kyy_phys   = ', Kyy_phys
    write(*,'(A,ES16.8)') 'Kxx_LB     = ', Kxx_LB
    write(*,'(A,ES16.8)') 'Kxy_LB     = ', Kxy_LB
    write(*,'(A,ES16.8)') 'Kyy_LB     = ', Kyy_LB
    write(*,'(A,ES16.8)') 'ux_LB      = ', ux_LB
    write(*,'(A,ES16.8)') 'uy_LB      = ', uy_LB
    write(*,*) '---------------------------------------------'

!$omp parallel do collapse(2) default(shared) private(i,j,xc,yc,phi_loc,Bx_loc,By_loc)
    do j = 1, ny
      do i = 1, nx
        xc = x(i)
        yc = y(j)

        phi_loc  = phi_exact_phys(xc,yc,t_phys)
        phi(i,j) = phi_loc

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
  subroutine collision_mrt_adv_aniso()
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
    real(wp) :: dm1, dm2
    real(wp) :: Gj

!$omp parallel do collapse(2) default(shared) private(i,j,a, &
!$omp& f0,f1,f2,f3,f4, &
!$omp& m0,m1,m2,m3,m4, &
!$omp& me0,me1,me2,me3,me4, &
!$omp& ms0,ms1,ms2,ms3,ms4, &
!$omp& phi_loc,Bx_loc,By_loc,dBx,dBy,M1Gx,M1Gy,dm1,dm2,Gj)
    do j = 1, ny
      do i = 1, nx

        phi_loc = phi(i,j)

        Bx_loc = ux_LB * phi_loc
        By_loc = uy_LB * phi_loc

        dBx   = Bx_loc - Bx_prev(i,j)
        dBy   = By_loc - By_prev(i,j)
        Bx_prev(i,j) = Bx_loc
        By_prev(i,j) = By_loc

        ! M1G = (I - S1/2) dB
        M1Gx = dBx - 0.5_wp*( s11*dBx + s12*dBy )
        M1Gy = dBy - 0.5_wp*( s12*dBx + s22*dBy )

        ! moments
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

        ! equilibrium moments
        me0 = phi_loc
        me1 = Bx_loc           ! Bx / (3 cs2) = Bx
        me2 = By_loc           ! By / (3 cs2) = By
        me3 = -phi_loc / 3.0_wp
        me4 = -phi_loc / 3.0_wp

        ! MRT relaxation (anisotropic for m1,m2)
        ms0 = m0

        dm1 = m1 - me1
        dm2 = m2 - me2
        ms1 = m1 - ( s11*dm1 + s12*dm2 )
        ms2 = m2 - ( s12*dm1 + s22*dm2 )

        ms3 = m3 - s3*(m3 - me3)
        ms4 = m4 - s4*(m4 - me4)

        ! inverse transform (D2Q5 M^{-1})
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

        ! add G_j
        do a = 0, 4
          Gj = w(a) * ( ex(a)*M1Gx + ey(a)*M1Gy ) / cs2
          f_post(i,j,a) = f_post(i,j,a) + Gj
        end do

      end do
    end do
!$omp end parallel do

  end subroutine collision_mrt_adv_aniso

  !---------------------------------
  subroutine streaming()
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
        den = den +  exact**2
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

end program lbm_gauss_adv_aniso_d2q5
