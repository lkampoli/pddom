module halfspace
  !-----------------------------------------------------------------------------
  ! Module: halfspace
  !  This modules contains subroutines specific to problems with a halfspace
  !  geometry.
  !
  ! Subroutines:
  !  driver_reflectance (public)
  !    Compute the Fourier modes of the polarized bidirectional reflectance
  !    from a uniform halfspace.
  !
  !  vrte0_iq (public)
  !    Compute the upper-left block of the discretize system corresponding to
  !    the reduced form of the zeroth Fourier mode of the vector radiative
  !    transport equation (vRTE).
  !
  !  vrte0_uv (public)
  !    Compute the discretize system corresponding to the bottom-right block of
  !    the zeroth Fourier mode of the reduced form of the vector radiative
  !    transport equation (vRTE).
  !
  !  source (public)
  !    Compute the discretized source terms for the reduced problem
  !    resulting from the direct portion of the Stokes vector.
  !
  !  source0_iq (public)
  !    Compute the I and Q parts of the discretized source terms for the
  !    zeroth Fourier mode of the reduced problem resulting from the direct
  !    portion of the Stokes vector.
  !
  !  source0_uv (public)
  !    Compute the U and V parts of the discretized source terms for the
  !    zeroth Fourier mode of the reduced problem resulting from the direct
  !    portion of the Stokes vector.
  !
  ! Useful Citations:
  !  Siewert, C. E. "A discrete-ordinates solution for radiative-transfer
  !  models that include polarization effects." JQRST, 2000.
  !---
  ! Developed by J.P. Dark
  ! Contact: email@jpdark.com
  ! Last modified June 2017
  !-----------------------------------------------------------------------------

  implicit none

  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.d0)

  private
  public :: driver_reflectance, dirichlet, evaluate

contains


  subroutine driver_reflectance(N, L, kmax, tol, albedo, mu0, mu, wt, &
       FL_coeff, Rk)
    !!--------------------------------------------------------------------------
    !! Name: driver_reflectance
    !!  Compute the Fourier modes of the polarized bidirectional reflectance
    !!  from a uniform halfspace.
    !!
    !! Synopsis
    !!  driver_reflectance(N, L, kmax, tol, albedo, mu0, mu, wt, FL_coeff, Rk)
    !!
    !!    INPUT
    !!      Integer ------------- N, L, kmax
    !!      Double precision ---- tol, albedo, mu0, mu(N), wt(N),
    !!                            FL_coeff(6, L+1)
    !!
    !!    OUTPUT
    !!      Double precision ---- Rk(4N, 4, kmax+1)
    !!
    !! Purpose
    !!   Compute the Fourier modes, Rk(mu, mu0), of the reflection matrix,
    !!   R(mu, mu0, phi-phi0), using the discrete ordinate method where the
    !!   reflection matrix gives the diffuse intensity reflected from the
    !!   scattering halfspace,
    !!
    !!     Id(0, -mu, phi) = mu0/pi * R(mu, mu0, phi-phi0) I0.
    !!
    !! Arguments
    !!  N (input) integer
    !!   The number of ordinate quadrature points.
    !!
    !!  L (input) integer
    !!   The order of the expansion of the scattering matrix.
    !!
    !!  kmax (input) integer
    !!   Maximum Fourier mode. Must have kmax <= L.
    !!
    !!  tol (input) double precision
    !!   Numerical precision threshold for comparing eigenvalues.
    !!
    !!  albedo (input) double precision
    !!   Single scattering albedo.
    !!
    !!  mu0 (input) double precision
    !!   Ordinate direction of incident light.
    !!
    !!  mu (input) double precision, array(N)
    !!   Quadrature of ordinate directions.
    !!
    !!  wt (input) double precision, array(M)
    !!   Associated weights of the quadrature for ordinates.
    !!
    !!  FL_coeff (input) double precision, array(6, L+1)
    !!   Expansion coefficients of the scattering matrix.
    !!
    !!  Rk (output) double precision, array(4N, 4, kmax+1)
    !!   Reflectance matrix for the kth Fourier mode.
    !!--------------------------------------------------------------------------
    use special, only: auxiliary
    use discretize, only: vrte, source, vrte0_iq, vrte0_uv, source0_iq, &
         source0_uv
    use solver, only: vdisord, vdisord0_iq

    implicit none
    ! INPUT
    integer(8) :: N, L, kmax
    real(dp), intent(in) :: tol
    real(dp), intent(in) :: albedo
    real(dp), intent(in) :: mu0
    real(dp), intent(in), dimension(:) :: mu, wt ! dim(N)
    real(dp), intent(in), dimension(:, :) :: FL_coeff ! dim(6, L+1)
    ! OUTPUT
    real(dp), intent(out), dimension(4*N, 4, kmax+1) :: Rk
    ! Workspace: general
    integer(8) :: mm
    real(dp), dimension(N+1) :: x
    integer(8) :: k, nsource, nrows, nreal, nzero
    ! Workspace: zeroth Fourier mode
    real(dp), dimension(2*N, 2*N) :: A0, B0, AB0
    real(dp), dimension(2*N, 2) :: Psip0, Gammap0, Ip_pos0, Ip_neg0, reflmat0
    real(dp), dimension(2*N, 2*N) :: Psih0, Gammah0, Ih_pos0, Ih_neg0
    real(dp), dimension(2*N) :: lambda0
    real(dp), dimension(2*N, 2) :: coeff0
    ! Workspace: higher Fourier modes
    real(dp), dimension(L+1, N+1) :: plk, rlk, tlk
    real(dp), dimension(4*N, 4*N) :: A, B, AB
    real(dp), dimension(4*N, 4) :: Psi_p, Gamma_p
    real(dp), dimension(4*N, 4*N) :: Psi_h, Gamma_h
    real(dp), dimension(4*N, 4) :: coeff_pos
    real(dp), dimension(4*N, 4*N) :: Ih_pos, Ih_neg
    real(dp), dimension(4*N, 4) :: Ip_pos, Ip_neg
    real(dp), dimension(4*N) :: lambda
    real(dp), dimension(4*N, 4) :: reflmat

    ! Create extended ordinate array.
    x(1:N) = mu
    x(N+1) = mu0

    !---------------------------------------------------------------------------
    ! Compute the discretized system for k=0.
    !
    k = 0
    nsource = 2
    nrows = 2*N
    call auxiliary(k, N+1, L, x, plk, rlk, tlk)
    Rk(:, :, 1) = 0.0_dp
    !--
    ! 1. Compute the discretized system for k=0, m=1.
    !
    call vrte0_iq(L, N, nrows, albedo, mu, wt, plk, rlk, FL_coeff, A0, B0, AB0)
    call source0_iq(L, N, albedo, mu, plk, rlk, FL_coeff, Psip0, Gammap0)
    !
    ! 2. Compute Psi and Gamma from the discrete system.
    !
    call vdisord0_iq(nrows, nsource, mu0, A0, B0, AB0, tol, &
         Psip0, Gammap0, Psih0, Gammah0, lambda0, nreal, nzero)
    if (nzero == 1) then
       do mm = 0, N-1
          Psih0(2*mm+1, 1) = 1.0_dp
          Psih0(2*mm+2, 1) = 0.0_dp
       enddo
    endif
    !
    ! 3. Compute the Stokes vector from Psi and Gamma.
    ! The Stokes vector for decaying planewave solutions is given by
    !
    !  I = [[Psi + Gamma],
    !       [D34 * (Psi - Gamma)]].
    !
    if (nzero == 1) then
       Ih_pos0(:, 1) = Psih0(:, 1)
       Ih_neg0(:, 1) = Psih0(:, 1)
    endif
    Ih_pos0(:, nzero+1 : nrows) = &
         Psih0(:, nzero+1 : nrows) + Gammah0(:, nzero+1 : nrows)
    Ih_neg0(:, nzero+1 : nrows) = &
         Psih0(:, nzero+1 : nrows) - Gammah0(:, nzero+1 : nrows)
    Ip_pos0 = Psip0 + Gammap0
    Ip_neg0 = Psip0 - Gammap0
    !
    ! 4. Use the boundary conditions to compute the expansion coefficients.
    !
    coeff0 = -Ip_pos0
    call dirichlet(nrows, nsource, Ih_pos0, coeff0)
    !
    ! 5.  Compute the Fourier modes of the reflectance matrix.
    ! The reflectance matrix is given by
    !
    !   R0 = Ih_neg0 * coeff0 + Ip_neg0.
    !
    reflmat0 = Ip_neg0
    call dgemm('N', 'N', nrows, nsource, nrows, 1.0_dp, Ih_neg0, nrows, &
         coeff0, nrows, 1.0_dp, reflmat0, nrows)
    do mm = 0 , N-1
       Rk(4*mm+1:4*mm+2, 1:2, k+1) = reflmat0(2*mm+1:2*mm+2,:)
    enddo
    !---
    !
    ! 1. Compute the discretized system for k=0, m=1.
    !
    call vrte0_uv(L, N, nrows, albedo, mu, wt, plk, rlk, FL_coeff, A0, B0, AB0)
    call source0_uv(L, N, albedo, mu, plk, rlk, FL_coeff, Psip0, Gammap0)
    !
    ! 2. Compute Psi and Gamma from the discrete system.
    !
    call vdisord(nrows, nsource, mu0, A0, B0, AB0, tol, &
         Psip0, Gammap0, Psih0, Gammah0, lambda0, nreal, nzero)
    !
    ! 3. Compute the Stokes vector from Psi and Gamma.
    ! The Stokes vector for decaying planewave solutions is given by
    !
    !  I = [[Psi + Gamma],
    !       [D34 * (Psi - Gamma)]].
    !
    Ih_pos0 = Psih0 + Gammah0
    Ih_neg0 = -Psih0 + Gammah0
    Ip_pos0 = Psip0 + Gammap0
    Ip_neg0 = -Psip0 + Gammap0
    !
    ! 4. Use the boundary conditions to compute the expansion coefficients.
    !
    coeff0 = -Ip_pos0
    call dirichlet(nrows, nsource, Ih_pos0, coeff0)
    !
    ! 5.  Compute the Fourier modes of the reflectance matrix.
    ! The reflectance matrix is given by
    !
    !   R0 = Ih_neg0 * coeff0 + Ip_neg0.
    !
    reflmat0 = Ip_neg0
    call dgemm('N', 'N', nrows, nsource, nrows, 1.0_dp, Ih_neg0, nrows, &
         coeff0, nrows, 1.0_dp, reflmat0, nrows)
    do mm = 0 , N-1
       Rk(4*mm+3:4*mm+4, 3:4, k+1) = reflmat0(2*mm+1:2*mm+2,:)
    enddo
    !---------------------------------------------------------------------------


    nrows = 4*N
    nsource = 4
    do k = 1, kmax
       !
       ! 1. Compute discretized system.
       !
       call auxiliary(k, N+1, L, x, plk, rlk, tlk)
       call vrte(k, L, N, nrows, albedo, mu, wt, plk, rlk, tlk, FL_coeff, &
            A, B, AB)
       call source(k, L, N, albedo, mu, plk, rlk, tlk, FL_coeff, Psi_p, Gamma_p)
       !
       ! 2. Compute components of the general solution for the Fourier mode.
       !
       !      I_pos(tau, mu) = Psi(tau, mu) + Gamma(tau, mu),
       !      I_neg(tau, mu) = D34 * [Psi(tau, mu) - Gamma(tau, mu)],
       !
       ! where D34 = (1, 1, -1, -1).
       !
       call vdisord(nrows, nsource, mu0, A, B, AB, tol, &
            Psi_p, Gamma_p, Psi_h, Gamma_h, lambda, nreal, nzero)
       Ip_pos = Psi_p + Gamma_p
       Ip_neg = Psi_p - Gamma_p
       do mm = 0, N-1
          Ip_neg(4*mm+1 : 4*mm+2, :) = Psi_p(4*mm+1 : 4*mm+2, :) &
               - Gamma_p(4*mm+1 : 4*mm+2, :)
          Ip_neg(4*mm+3 : 4*mm+4, :) = -Psi_p(4*mm+3 : 4*mm+4, :) &
               + Gamma_p(4*mm+3 : 4*mm+4, :)
       enddo
       Ih_pos = Psi_h + Gamma_h
       do mm = 0, N-1
          Ih_neg(4*mm+1 : 4*mm+2, :) = Psi_h(4*mm+1 : 4*mm+2, :) &
               - Gamma_h(4*mm+1 : 4*mm+2, :)
          Ih_neg(4*mm+3 : 4*mm+4, :) = -Psi_h(4*mm+3 : 4*mm+4, :) &
               + Gamma_h(4*mm+3 : 4*mm+4, :)
       enddo
       !
       ! 3. Compute boundary conditions.
       !
       coeff_pos = -Ip_pos
       call dirichlet(nrows, nsource, Ih_pos, coeff_pos)
       !
       ! 4. Compute the Fourier modes of the reflectance matrix.
       ! The reflectance matrix is given by
       !
       !    Rk = Ih_neg * coeff_pos + Ip_neg
       !
       reflmat = Ip_neg
       call dgemm('N', 'N', nrows, nsource, nrows, 1.0_dp, Ih_neg, nrows, &
            coeff_pos, nrows, 1.0_dp, reflmat, nrows)
       Rk(:, :, k+1) = reflmat
    enddo

  end subroutine driver_reflectance


  subroutine dirichlet(nrows, nsource, Ih_pos, coeff_pos)
    !!--------------------------------------------------------------------------
    !! Name: dirichlet
    !!  Compute the expansion coefficients for scattering on a halfspace with
    !!  Dirichlet boundary conditions.
    !!
    !! Synopsis
    !!  subroutine dirichlet(nrows, nsource, Ih_pos, coeff_pos)
    !!
    !!     INPUT
    !!       Integer ------------- nrows, nsource
    !!       Double precision ---- Ipos(nrows, nrows)
    !!
    !!     IN/OUT
    !!       Double precsion ---- coeff_pos(nrows, nsource)
    !!
    !! Purpose
    !!  Solve for the expansion coefficients, coeff_pos, for Dirichlet
    !!  boundary conditions,
    !!
    !!    sum[j=1,4N] coeff_pos[j] * Ih_pos[j](0)  = I1 - Ip_pos(0).
    !!
    !! Arguments:
    !!  nrows (input) integer
    !!    The number of rows in the discretized system.
    !!
    !!  nsource (input) integer
    !!    The number of source terms.
    !!
    !!  Ih_pos (input) double precision, array(nrows, nrows)
    !!    The homogeneous solution at positive quadrature points
    !!
    !!  coeff_pos(in/out) double precision, array(nrows, source)
    !!    On input coeff_pos is the right-hand side of the boundary conditions,
    !!
    !!       coeff_pos = I1 - Ip_pos(0),
    !!
    !!    and on output coeff_pos is the expansion coefficients.
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8), intent(in) :: nrows, nsource
    real(dp), intent(in), dimension(:, :) :: Ih_pos ! dim(nrows, nrows)
    ! OUTPUT
    real(dp), intent(inout), dimension(:, :) :: coeff_pos ! dim(nrows, nsource)
    ! Workspace
    real(dp), dimension(nrows, nrows) :: workspace
    integer(8), dimension(nrows) :: IPIV
    integer(8) :: info
    !---------------------------------------------------------------------------
    ! Solve for coeff_pos using
    !
    !  sum[j=1,4N] coeff_pos[j] * Ipos[j](0)  = I1 - Ip(0).
    !
    workspace = Ih_pos
    info = 0
    call dgesv(nrows, nsource, workspace, nrows, IPIV, coeff_pos, nrows, info)
    if (info > 0) then
       write (*, '(A)') "> PROGRAM STOPPED BEFORE COMPLETION."
       write (*, '(A)') "> ERROR in __solver_MOD_dirichlet."
       write (*, '(A)') "> PLU decomposition of the boundary conditions failed."
       write (*, '(A, I4, A, I4, A)') "> U(", INFO, "," ,INFO,") = 0."
       write (*, '(A, I4, A)') "> Dimension of system is ", nrows, "."
       stop
    endif
    !---------------------------------------------------------------------------
  end subroutine dirichlet


  subroutine evaluate(z, nrows, ncols, nsource, nreal, tol, &
       Ih, Ip, mu0, lambda, coeff, Stokes)
    !!--------------------------------------------------------------------------
    !! Name: evaluate
    !!  Compute a Fourier mode of the Stokes vector evaluated at z.
    !!
    !! Synopsis
    !!  subroutine evaluate(z, nrows, ncols, nsource, nreal, tol, Ih, Ip, mu0,
    !!                      lambda, coeff, Stokes)
    !!
    !!     INPUT
    !!       Integer ------------- nrows, ncols, nsource, nreal
    !!       Double precision ---- z, tol, Ih, Ip, mu0, Ipos(nrows, ncols),
    !!                             lambda(ncols), coeff(ncols)
    !!
    !!     OUTPUT
    !!       Double precsion ---- Stokes(nrows, nsource)
    !!
    !! Purpose
    !!  Compute a Fourier mode of the Stokes vector at z from the general
    !!  solution using the general solution,
    !!
    !!   I = Ip * exp(-z/mu0)
    !!    + sum[j=1,nreal] coeff[j] * Ih[j] * exp(-lambda[j] * z)
    !!    + sum[j=nreal+1,nrows] (coeff(j) * Re[Ih[j] * exp(-lambda[j] * z)]
    !!                        + coeff(j+1) * Im[Ih[j] * exp(-lambda[j] * z)])
    !!
    !!
    !! Arguments:
    !!  z (input) double precisions
    !!   The depth at which to evaluate the Stokes vector.
    !!
    !!  nrows (input) integer
    !!    The number of rows in the discretized system.
    !!
    !!  ncols (input) integer
    !!    The number of planewave solutions in the homogeneous solution.
    !!
    !!  nsource (input) integer
    !!    The number of source terms.
    !!
    !!  nreal (input) integer
    !!    The number of eigenvalues, lambda, which are strictly real.
    !!
    !!  Ih (input) double precision, array(nrows, ncols)
    !!    The eigenvectors that make up the homogeneous solution.
    !!
    !!  Ip (input) double precision, array(nrows, nsource)
    !!    The Stokes vector for particular solution for each source.
    !!
    !!  mu0 (input) double precision
    !!    The incident polar direction.
    !!
    !!  coeff(in/out) double precision, array(nrows, source)
    !!    The expansion coefficients for the homogeneous solution.
    !!
    !!  Stokes (output) double precision
    !!    A Fourier mode of the Stokes vector evaluated at z.
    !!--------------------------------------------------------------------------
    implicit none
    ! Input
    real(dp) :: z
    integer(8), intent(in) :: nrows, ncols, nsource, nreal
    real(dp) :: tol
    real(dp), dimension(:, :) :: Ih ! dim(nrows, ncols)
    real(dp), dimension(:, :) :: Ip ! dim(nrows, nsource)
    real(dp) :: mu0
    real(dp), dimension(:) :: lambda ! dim(ncols)
    real(dp), dimension(:, :) :: coeff ! dim(ncols, nsource)
    ! Output
    real(dp), dimension(:, :) :: Stokes ! dim(nrows, nsource)
    ! Workspace
    real(dp), dimension(nrows, ncols) :: work
    integer(8) :: jj

    !---------------------------------------------------------------------------
    ! Compute a Fourier mode of the Stokes vector at z from the general
    ! solution,
    !
    !   I = Ip * exp(-z/mu0)
    !    + sum[j=1,nreal] coeff[j] * Ih[j] * exp(-lambda[j] * z)
    !    + sum[j=nreal+1,nrows] (coeff(j) * Re[Ih[j] * exp(-lambda[j] * z)]
    !                           + coeff(j+1) * Im[Ih[j] * exp(-lambda[j] * z)])
    if (abs(z) < tol) then
       !
       work = Ih
    else
       do jj = 1, nreal
          work(:, jj) = Ih(:, jj) * exp(-lambda(jj) * z)
       enddo
       do jj = nreal+1, ncols, 2
          !
          ! Re[Ih * exp(-lambda *z)] =
          !   Re[Ih] * cos(Im[lambda]*z) * exp(-lambda*z)
          !   - Im[Ih] * sin(-Im[lambda]*z) * exp(-lambda*z)
          !
          ! Im[Ih * exp(-lambda *z)] =
          !   Im[Ih] * cos(Im[lambda]*z) * exp(-lambda*z)
          !   + Re[Ih] * sin(-Im[lambda]*z) * exp(-lambda*z)
          work(:, jj) = exp(-lambda(jj) * z) &
               * (cos(lambda(jj+1) * z) * Ih(:, jj) &
               + sin(lambda(jj+1) * z) * Ih(:, jj+1))
          work(:, jj+1) = exp(-lambda(jj) * z) &
               * (cos(lambda(jj+1) * z) * Ih(:, jj+1) &
               - sin(lambda(jj+1) * z) * Ih(:, jj))
       enddo
    endif
    Stokes = Ip * exp(-z / mu0)
    call dgemm('N', 'N', nrows, nsource, ncols, 1.0_dp, work, nrows, &
         coeff, ncols, 1.0_dp, Stokes, nrows)
    !---------------------------------------------------------------------------
  end subroutine evaluate


end module halfspace
