module discretize
  !-----------------------------------------------------------------------------
  ! Module: discretize
  !
  !  This module contains subroutines that discretize the radiative transport
  !  equation (RTE) over the polar cosine (ordinate) directions.
  !
  !  The discretization uses the inherent symmetries in the RTE to reduce the
  !  dimensions of the problem. This method is often referred to as the
  !  double-Gauss method.
  !
  ! Subroutines:
  !  vrte (public)
  !    Compute the discretize system corresponding to the reduced form of the
  !    kth Fourier mode of the vector radiative transport equation (vRTE).
  !
  !   vrte0_iq (public)
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
  !  models that include polarization effects." JQSRT 64.3 (2000): 227-254.
  !---
  ! Developed by J.P. Dark
  ! Contact: email@jpdark.com
  ! Last modified June 2017
  !-----------------------------------------------------------------------------

  implicit none

  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: PI = 3.141592653589793_dp

  private
  public :: vrte, vrte0_iq, vrte0_uv, source, source0_iq, source0_uv


contains


  subroutine vrte(k, L, N, nrows, albedo, mu, wt, plk, rlk, tlk, FL_coeff, &
       A, B, AB)
    !! --------------------------------------------------------------------------
    !! Name: vrte
    !!  Compute the discretize system corresponding to the reduced form of the
    !!  kth Fourier mode of the vector radiative transport equation (vRTE).
    !!
    !! Synopsis
    !!  subroutine vrte(k, L, N, nrows, albedo, mu, wt, plk, rlk, tlk,
    !!                   FL_coeff, A, B, AB)
    !!
    !!    INPUT
    !!      integer ----------- k, L, N, nrows
    !!      double precision -- albedo, mu(N), wt(N), plk(L+1, N),
    !!                           rlk(L+1, N), tlk(L+1, N), FL_coeff(6, L+1),
    !!    OUTPUT
    !!      double precision -- A(nrows, nrows), B(nrows, nrows)
    !!
    !! Purpose
    !!  Compute the discretized radiative transport matrix for the kth
    !!  Fourier mode. That is, compute
    !!
    !!    A = M^{-1} * (-I + 0.5 * albedo (Zpos - D34 * Zneg) * W)
    !!      = M^{-1} * (-I + 0.5 * alebdo * sum[l=k,L] Aux[l,k]
    !!                  * (1 - (-1)^{j-k} * D_{34}) * F[j] * Aux[l,k]^T * W)
    !!
    !!  and
    !!
    !!    B = M^{-1} * (-I + 0.5 * albedo (Zpos + D34 * Zneg) * W)
    !!      = M^{-1} * (-I + 0.5 * alebdo * sum[l=k,L] Aux[l,k]
    !!           * (1 + (-1)^{j-k} * D_{34}) * F[j] * Aux[l,k]^T * W)
    !!
    !!  so that Psi and Gamma solve the system
    !!
    !!      Psi'(tau) = A * Gamma(tau) + src1 * exp(-tau/mu0),
    !!    Gamma'(tau) = B *   Psi(tau) + src2 * exp(-tau/mu0).
    !!
    !!  Note: D34 is a 4Nx4N matrix with (1,1,-1,-1) repeated along the
    !!  diagonal N times,
    !!
    !!    D34 = kronecker_product(Identity(N, N), (1,1,-1,-1))
    !!
    !!  and Aux[l,k] is a 4N x 4 matrix of the auxiliary functions.
    !!
    !! Arguments
    !!  k (input) integer
    !!    The Fourier mode of the azimuth.
    !!
    !!  L (input) integer
    !!    The order of the expansion in auxiliary rotation functions.
    !!
    !!  N (input) integer
    !!    The number of ordinate quadrature points.
    !!
    !!  nrows (input) integer
    !!    The number of rows in the discretized system, nrows = 4*N.
    !!
    !!  albedo (input) double precision
    !!    The single scattering albedo.
    !!
    !!  mu (input) double precision, array(N)
    !!    Quadrature of ordinate directions.
    !!
    !!  wt (input) double precision, array(M)
    !!    Associated weights of the quadrature for ordinates.
    !!
    !!  FL_coeff (input) double precision, array(L+1, 4*N, 4*N)
    !!    Expansion coefficients of the scattering matrix.
    !!
    !!  A (output) double precision, array(4*N, 4*N)
    !!    Discretization of reduced system,
    !!      A = (-I + 0.5 * albedo (Zpos - D34 * Zneg) * W) M^{-1}.
    !!
    !!  B (output) double precision array, dimension (4*N, 4*N)
    !!    Discretization of reduced system,
    !!      B = (-I + 0.5 * albedo (Zpos + D34 * Zneg) * W) M^{-1}.
    !!
    !!  AB (output) double precision array, dimension (4*N, 4*N)
    !!    Discretization of reduced system,
    !!      AB = A * B.
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8) :: k, L, N, nrows
    real(dp) :: albedo
    real(dp), intent(in), dimension(:) :: mu, wt ! dim(N)
    real(dp), intent(in), dimension(:, :) :: FL_coeff ! dim(6, L+1)
    real(dp), intent(in), dimension(:, :) :: plk, rlk, tlk ! dim(L+1, N)
    ! OUTPUT
    real(dp), intent(out), dimension(:, :) :: A, B, AB ! dim(4N, 4N)
    ! Workspace
    integer(8) :: nn, ii, jj
    integer(8), dimension(4) :: indx
    real(dp), dimension(4*N, 2) :: AUX1, AUX2
    real(dp), dimension(2, 4*N) :: F1, F2
    real(dp) :: tmp1

    !---------------------------------------------------------------------------
    ! Compute
    !  A = Zpos - D34 * Zneg
    ! and
    !  B = Zpos + D34 * Zneg.
    !
    A = 0.0_dp
    B = 0.0_dp
    F1 = 0.0_dp
    F2 = 0.0_dp
    AUX1 = 0.0_dp
    AUX2 = 0.0_dp
    do jj = k, L, 2
       ! Construct the slices of the matrix of auxiliary functions
       do nn = 0, N-1
          AUX1(4*nn+1, 1) = plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          AUX1(4*nn+2, 2) = rlk(jj+1, nn+1)
          AUX1(4*nn+3, 2) = -tlk(jj+1, nn+1)
       enddo
       ! Construct the slices of the matrix of auxiliary functions
       do nn = 0, N-1
          AUX2(4*nn+2, 1) = -tlk(jj+1, nn+1)
          AUX2(4*nn+3, 1) = rlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          AUX2(4*nn+4, 2) = plk(jj+1, nn+1)
       enddo
       ! Construct F1 = (F[l] * Aux[l,k])(1:2, :)
       do nn = 0, N-1
          ! First column.
          F1(1, 4*nn+1) = FL_coeff(1, jj+1) * plk(jj+1, nn+1)
          F1(2, 4*nn+1) = FL_coeff(5, jj+1) * plk(jj+1, nn+1)
          ! Second column.
          F1(1, 4*nn+2) = FL_coeff(5, jj+1) * rlk(jj+1, nn+1)
          F1(2, 4*nn+2) = FL_coeff(2, jj+1) * rlk(jj+1, nn+1)
          ! Third column.
          F1(1, 4*nn+3) = -FL_coeff(5, jj+1) * tlk(jj+1, nn+1)
          F1(2, 4*nn+3) = -FL_coeff(2, jj+1) * tlk(jj+1, nn+1)
       enddo
       ! Construct F2 = (F[l] * Aux[l,k])(3:4, :)
       do nn = 0, N-1
          ! Second column.
          F2(1, 4*nn+2) = -FL_coeff(3, jj+1) * tlk(jj+1, nn+1)
          F2(2, 4*nn+2) =  FL_coeff(6, jj+1) * tlk(jj+1, nn+1)
          ! Third column.
          F2(1, 4*nn+3) =  FL_coeff(3, jj+1) * rlk(jj+1, nn+1)
          F2(2, 4*nn+3) = -FL_coeff(6, jj+1) * rlk(jj+1, nn+1)
          ! Fourth column.
          F2(1, 4*nn+4) = FL_coeff(6, jj+1) * plk(jj+1, nn+1)
          F2(2, 4*nn+4) = FL_coeff(4, jj+1) * plk(jj+1, nn+1)
       enddo
       ! j-k is even.
       ! A += Aux[j,k] * (1 - D_{34}) * F[j] * Aux[j,k]^T
       ! B += Aux[j,k] * (1 + D_{34}) * F[j] * Aux[j,k]^T
       call dgemm('N', 'N', nrows, nrows, 2, 2.0_dp, AUX2, nrows, &
            F2, 2, 1.0_dp, A, nrows)
       call dgemm('N', 'N', nrows, nrows, 2, 2.0_dp, AUX1, nrows, &
            F1, 2, 1.0_dp, B, nrows)
    enddo
    do jj = k+1, L, 2
       ! Construct the slices of the matrix of auxiliary functions
       do nn = 0, N-1
          AUX1(4*nn+1, 1) = plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          AUX1(4*nn+2, 2) = rlk(jj+1, nn+1)
          AUX1(4*nn+3, 2) = -tlk(jj+1, nn+1)
       enddo
       ! Construct the slices of the matrix of auxiliary functions
       do nn = 0, N-1
          AUX2(4*nn+2, 1) = -tlk(jj+1, nn+1)
          AUX2(4*nn+3, 1) = rlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          AUX2(4*nn+4, 2) = plk(jj+1, nn+1)
       enddo
       ! Construct F1 = (F[l] * Aux[l,k])(1:2, :)
       do nn = 0, N-1
          ! First column.
          F1(1, 4*nn+1) = FL_coeff(1, jj+1) * plk(jj+1, nn+1)
          F1(2, 4*nn+1) = FL_coeff(5, jj+1) * plk(jj+1, nn+1)
          ! Second column.
          F1(1, 4*nn+2) = FL_coeff(5, jj+1) * rlk(jj+1, nn+1)
          F1(2, 4*nn+2) = FL_coeff(2, jj+1) * rlk(jj+1, nn+1)
          ! Third column.
          F1(1, 4*nn+3) = -FL_coeff(5, jj+1) * tlk(jj+1, nn+1)
          F1(2, 4*nn+3) = -FL_coeff(2, jj+1) * tlk(jj+1, nn+1)
       enddo
       ! Construct F2 = (F[l] * Aux[l,k])(3:4, :)
       do nn = 0, N-1
          ! Second column.
          F2(1, 4*nn+2) = -FL_coeff(3, jj+1) * tlk(jj+1, nn+1)
          F2(2, 4*nn+2) =  FL_coeff(6, jj+1) * tlk(jj+1, nn+1)
          ! Third column.
          F2(1, 4*nn+3) =  FL_coeff(3, jj+1) * rlk(jj+1, nn+1)
          F2(2, 4*nn+3) = -FL_coeff(6, jj+1) * rlk(jj+1, nn+1)
          ! Fourth column.
          F2(1, 4*nn+4) = FL_coeff(6, jj+1) * plk(jj+1, nn+1)
          F2(2, 4*nn+4) = FL_coeff(4, jj+1) * plk(jj+1, nn+1)
       enddo
       ! j-k is odd.
       ! A += Aux[j,k] * (1 + D_{34}) * F[j] * Aux[j,k]^T
       ! B += Aux[j,k] * (1 - D_{34}) * F[j] * Aux[j,k]^T
       call dgemm('N', 'N', nrows, nrows, 2, 2.0_dp, AUX1, nrows, &
            F1, 2, 1.0_dp, A, nrows)
       call dgemm('N', 'N', nrows, nrows, 2, 2.0_dp, AUX2, nrows, &
            F2, 2, 1.0_dp, B, nrows)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = 0.5 * albedo (Zpos - D34 * Zneg) * W
    ! and
    !  B = 0.5 * albedo (Zpos + D34 * Zneg) * W.
    !
    do nn = 0, N-1
       indx = (/ (4*nn + ii, ii = 1, 4) /)
       tmp1 = 0.5_dp * albedo * wt(nn+1)
       A(:, indx) = tmp1 * A(:, indx)
       B(:, indx) = tmp1 * B(:, indx)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = -I + 0.5 * albedo (Zpos - D34 * Zneg) * W
    ! and
    !  B = -I + 0.5 * albedo (Zpos + D34 * Zneg) * W.
    !
    do nn = 1, nrows
       A(nn, nn) = A(nn, nn) - 1.0_dp
       B(nn, nn) = B(nn, nn) - 1.0_dp
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = M^{-1} * (-I + 0.5 * albedo (Zpos - D34 * Zneg) * W)
    ! and
    !  B = M^{-1} * (-I + 0.5 * albedo (Zpos + D34 * Zneg) * W).
    !
    do nn = 0, N-1
       indx = (/ (4*nn + ii, ii = 1, 4) /)
       A(indx, :) = A(indx, :) / mu(nn+1)
       B(indx, :) = B(indx, :) / mu(nn+1)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  AB = A * B.
    call dgemm('N', 'N', nrows, nrows, nrows, 1.0_dp, A, nrows, B, nrows, &
         0.0_dp, AB, nrows)
    !---------------------------------------------------------------------------

  end subroutine vrte


  subroutine vrte0_iq(L, N, nrows, albedo, mu, wt, plk, rlk, FL_coeff, &
       A, B, AB)
    !!--------------------------------------------------------------------------
    !! Name: vrte0_iq
    !!  Compute the upper-left block of the discretize system corresponding to
    !!  the reduced form of the zeroth Fourier mode of the vector radiative
    !!  transport equation (vRTE).
    !!
    !! Synopsis
    !!  subroutine vrte(L, N, nrows, albedo, mu, wt, plk, rlk, tlk, FL_coeff,
    !!                  A, B, AB)
    !!
    !!     INPUT
    !!       integer ----------- L, N, nrows
    !!       double precision -- albedo, mu(N), wt(N), plk(L+1, N),
    !!                           rlk(L+1, N), tlk(L+1, N), FL_coeff(6, L+1),
    !!     OUTPUT
    !!       double precision -- A(nrows, nrows), B(nrows, nrows)
    !!
    !! Purpose
    !!   Compute the upper-left block of the discretized vRTE of the zeroth
    !!   Fourier mode. That is, compute the upper block:
    !!
    !!     A = (-I + 0.5 * albedo (Zpos_UL - Zneg_UL) * W) M^{-1}
    !!       = (-I + 0.5 * alebdo * sum[l=k,L] Aux[l,k]
    !!            * (1 - (-1)^{j-k}) * F[j] * Aux[l,k]^T * W) M^{-1}
    !!
    !!   and
    !!
    !!     B = (-I + 0.5 * albedo (Zpos_UL + Zneg_UL) * W) M^{-1}
    !!       = (-I + 0.5 * alebdo * sum[l=k,L] Aux[l,k]
    !!            * (1 + (-1)^{j-k}) * F[j] * Aux[l,k]^T * W) M^{-1}
    !!
    !!   so that Psi and Gamma solve the system
    !!
    !!       Psi'(tau) = A * Gamma(tau) + src1 * exp(-tau/mu0),
    !!     Gamma'(tau) = B *   Psi(tau) + src2 * exp(-tau/mu0),.
    !!
    !! Arguments
    !!   L (input) integer
    !!     The order of the expansion in auxiliary rotation functions.
    !!
    !!   N (input) integer
    !!     The number of ordinate quadrature points.
    !!
    !!   nrows (input) integer
    !!     The number of rows in the discretized system, nrows = 4*N.
    !!
    !!   albedo (input) double precision
    !!     The single scattering albedo.
    !!
    !!   mu (input) double precision, array(N)
    !!     Quadrature of ordinate directions.
    !!
    !!   wt (input) double precision, array(M)
    !!     Associated weights of the quadrature for ordinates.
    !!
    !!   FL_coeff (input) double precision, array(L+1, 4*N, 4*N)
    !!     Expansion coefficients of the scattering matrix.
    !!
    !!   A (output) double precision, array(2*N, 2*N)
    !!     Discretization of the upper left-block of the reduced system,
    !!       A = (-I + 0.5 * albedo (Zpos_UL - Zneg_UL) * W) M^{-1}.
    !!
    !!   B (output) double precision array, dimension (2*N, 2*N)
    !!     Discretization of the upper left-block of the reduced system,
    !!       B = (-I + 0.5 * albedo (Zpos_UL + Zneg_UL) * W) M^{-1}.
    !!
    !!   AB (output) double precision array, dimension (2*N, 2*N)
    !!     Discretization of the upper left-block reduced system,
    !!       AB = A * B.
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8) :: L, N, nrows
    integer(8), parameter :: k = 0
    real(dp) :: albedo
    real(dp), intent(in), dimension(:) :: mu, wt ! dim(N)
    real(dp), intent(in), dimension(:, :) :: FL_coeff ! dim(6, L+1)
    real(dp), intent(in), dimension(:, :) :: plk, rlk ! dim(L+1, N)
    ! OUTPUT
    real(dp), intent(out), dimension(:, :) :: A, B, AB ! dim(2N, 2N)
    ! Workspace
    integer(8) :: nn, jj
    real(dp), dimension(2*N, 2) :: AUX1
    real(dp), dimension(2, 2*N) :: F1
    real(dp) :: tmp1

    !---------------------------------------------------------------------------
    ! Compute the upper-left blocks
    !   A = Zpos_UL - D34 * Zneg_UL
    ! and
    !   B = Zpos_UL + D34 * Zneg_UL.
    !
    A = 0.0_dp
    B = 0.0_dp
    F1 = 0.0_dp
    AUX1 = 0.0_dp
    do jj = k, L, 2
       ! Construct the slices of the matrix of auxiliary functions
       do nn = 0, N-1
          AUX1(2*nn+1, 1) = plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          AUX1(2*nn+2, 2) = rlk(jj+1, nn+1)
       enddo
       ! Construct F1 = (F[l] * Aux[l,k])(1:2, :)
       do nn = 0, N-1
          ! First column.
          F1(1, 2*nn+1) = FL_coeff(1, jj+1) * plk(jj+1, nn+1)
          F1(2, 2*nn+1) = FL_coeff(5, jj+1) * plk(jj+1, nn+1)
          ! Second column.
          F1(1, 2*nn+2) = FL_coeff(5, jj+1) * rlk(jj+1, nn+1)
          F1(2, 2*nn+2) = FL_coeff(2, jj+1) * rlk(jj+1, nn+1)
       enddo
       ! j-k is even.
       ! B += 2 * Aux[j,k] * F[j] * Aux[j,k]^T
       call dgemm('N', 'N', nrows, nrows, 2, 2.0_dp, AUX1, nrows, &
            F1, 2, 1.0_dp, B, nrows)
    enddo
    do jj = k+1, L, 2
       ! Construct the slices of the matrix of auxiliary functions
       do nn = 0, N-1
          AUX1(2*nn+1, 1) = plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          AUX1(2*nn+2, 2) = rlk(jj+1, nn+1)
       enddo
       ! Construct F1 = (F[l] * Aux[l,k])(1:2, :)
       do nn = 0, N-1
          ! First column.
          F1(1, 2*nn+1) = FL_coeff(1, jj+1) * plk(jj+1, nn+1)
          F1(2, 2*nn+1) = FL_coeff(5, jj+1) * plk(jj+1, nn+1)
          ! Second column.
          F1(1, 2*nn+2) = FL_coeff(5, jj+1) * rlk(jj+1, nn+1)
          F1(2, 2*nn+2) = FL_coeff(2, jj+1) * rlk(jj+1, nn+1)
       enddo
       ! j-k is odd.
       ! A += 2 * Aux[j,k] * F[j] * Aux[j,k]^T
       call dgemm('N', 'N', nrows, nrows, 2, 2.0_dp, AUX1, nrows, &
            F1, 2, 1.0_dp, A, nrows)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = 0.5 * albedo (Zpos_UL - Zneg_UL) * W
    ! and
    !  B = 0.5 * albedo (Zpos_UL + Zneg_UL) * W.
    !
    do nn = 0, N-1
       tmp1 = 0.5_dp * albedo * wt(nn+1)
       A(:, 2*nn+1 : 2*nn+2) = tmp1 * A(:, 2*nn+1 : 2*nn+2)
       B(:, 2*nn+1 : 2*nn+2) = tmp1 * B(:, 2*nn+1 : 2*nn+2)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = -I + 0.5 * albedo (Zpos_UL - Zneg_UL) * W
    ! and
    !  B = -I + 0.5 * albedo (Zpos_UL + Zneg_UL) * W.
    !
    do nn = 1, nrows
       A(nn, nn) = A(nn, nn) - 1.0_dp
       B(nn, nn) = B(nn, nn) - 1.0_dp
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = (-I + 0.5 * albedo (Zpos_UL - Zneg_UL) * W) * M^{-1}
    ! and
    !  B = (-I + 0.5 * albedo (Zpos_UL + Zneg_UL) * W) * M^{-1}.
    !
    do nn = 0, N-1
       A(2*nn+1 : 2*nn+2, :) = A(2*nn+1 : 2*nn+2, :) / mu(nn+1)
       B(2*nn+1 : 2*nn+2, :) = B(2*nn+1 : 2*nn+2, :) / mu(nn+1)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  AB = A * B.
    call dgemm('N', 'N', nrows, nrows, nrows, 1.0_dp, A, nrows, B, nrows, &
         0.0_dp, AB, nrows)
    !---------------------------------------------------------------------------

  end subroutine vrte0_iq


  subroutine vrte0_uv(L, N, nrows, albedo, mu, wt, plk, rlk, FL_coeff, &
       A, B, AB)
    !!--------------------------------------------------------------------------
    !! Name: vrte0_uv
    !!  Compute the discretize system corresponding to the bottom-right block of
    !!  the zeroth Fourier mode of the reduced form of the vector radiative
    !!  transport equation (vRTE).
    !!
    !! Synopsis
    !!   subroutine vrte(L, N, nrows, albedo, mu, wt, plk, rlk, tlk,
    !!                   FL_coeff, A, B, AB)
    !!
    !!     INPUT
    !!       integer ----------- L, N, nrows
    !!       double precision -- albedo, mu(N), wt(N), plk(L+1, N),
    !!                           rlk(L+1, N), tlk(L+1, N), FL_coeff(6, L+1)
    !!     OUTPUT
    !!       double precision -- A(nrows, nrows), B(nrows, nrows)
    !!
    !! Purpose
    !!   Compute the bottom-right block of the zeroth Fourier mode of the
    !!   discretized vRTE. That is, compute
    !!
    !!     A = (-I + 0.5 * albedo (Zpos_BR + Zneg_BR) * W) M^{-1}
    !!       = (-I + 0.5 * alebdo * sum[l=k,L] Aux[l,k]
    !!            * (1 + (-1)^{j-k}) * F[j] * Aux[l,k]^T * W) M^{-1}
    !!
    !!   and
    !!
    !!     B = (-I + 0.5 * albedo (Zpos_BR - Zneg_BR) * W) M^{-1}
    !!       = (-I + 0.5 * alebdo * sum[l=k,L] Aux[l,k]
    !!            * (1 - (-1)^{j-k}) * F[j] * Aux[l,k]^T * W) M^{-1}
    !!
    !!   so that Psi and Gamma solve the system
    !!
    !!       Psi'(tau) = A * Gamma(tau) + src1 * exp(-tau/mu0),
    !!     Gamma'(tau) = B *   Psi(tau) + src2 * exp(-tau/mu0),.
    !!
    !! Arguments
    !!   k (input) integer
    !!     The Fourier mode of the azimuth.
    !!
    !!   L (input) integer
    !!     The order of the expansion in auxiliary rotation functions.
    !!
    !!   N (input) integer
    !!     The number of ordinate quadrature points.
    !!
    !!   nrows (input) integer
    !!     The number of rows in the discretized system, nrows = 4*N.
    !!
    !!   albedo (input) double precision
    !!     The single scattering albedo.
    !!
    !!   mu (input) double precision, array(N)
    !!     Quadrature of ordinate directions.
    !!
    !!   wt (input) double precision, array(M)
    !!     Associated weights of the quadrature for ordinates.
    !!
    !!   FL_coeff (input) double precision, array(L+1, 4*N, 4*N)
    !!     Expansion coefficients of the scattering matrix.
    !!
    !!   A (output) double precision, array(2*N, 2*N)
    !!     Bottom-right corner of discretization of reduced system,
    !!       A = (-I + 0.5 * albedo (Zpos_BR + Zneg_BR) * W) M^{-1}.
    !!
    !!   B (output) double precision array, dimension (2*N, 2*N)
    !!     Bottom-right corner of discretization of reduced system,
    !!       B = (-I + 0.5 * albedo (Zpos_BR - Zneg_BR) * W) M^{-1}.
    !!
    !!   AB (output) double precision array, dimension (2*N, 2*N)
    !!     Bottom-right corner of discretization of reduced system,
    !!       AB = A * B.
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8) :: L, N, nrows
    integer(8), parameter :: k=0
    real(dp) :: albedo
    real(dp), intent(in), dimension(:) :: mu, wt ! dim(N)
    real(dp), intent(in), dimension(:, :) :: FL_coeff ! dim(6, L+1)
    real(dp), intent(in), dimension(:, :) :: plk, rlk ! dim(L+1, N)
    ! OUTPUT
    real(dp), intent(out), dimension(:, :) :: A, B, AB ! dim(2N, 2N)
    ! Workspace
    integer(8) :: nn, jj
    real(dp), dimension(2*N, 2) :: AUX2
    real(dp), dimension(2, 2*N) :: F2
    real(dp) :: tmp1

    !---------------------------------------------------------------------------
    ! Compute
    !  A = Zpos_BR - D34 * Zneg_BR
    ! and
    !  B = Zpos_BR + D34 * Zneg_BR.
    !
    A = 0.0_dp
    B = 0.0_dp
    F2 = 0.0_dp
    AUX2 = 0.0_dp
    do jj = k, L, 2
       ! Construct the slices of the matrix of auxiliary functions
       do nn = 0, N-1
          AUX2(2*nn+1, 1) = rlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          AUX2(2*nn+2, 2) = plk(jj+1, nn+1)
       enddo
       ! Construct F2 = (F[l] * Aux[l,k])(3:4, :)
       do nn = 0, N-1
          ! Third column.
          F2(1, 2*nn+1) =  FL_coeff(3, jj+1) * rlk(jj+1, nn+1)
          F2(2, 2*nn+1) = -FL_coeff(6, jj+1) * rlk(jj+1, nn+1)
          ! Fourth column.
          F2(1, 2*nn+2) = FL_coeff(6, jj+1) * plk(jj+1, nn+1)
          F2(2, 2*nn+2) = FL_coeff(4, jj+1) * plk(jj+1, nn+1)
       enddo
       ! j-k is even.
       ! A += 2 * Aux[j,k] * F[j] * Aux[j,k]^T
       call dgemm('N', 'N', nrows, nrows, 2, 2.0_dp, AUX2, nrows, &
            F2, 2, 1.0_dp, A, nrows)
    enddo
    do jj = k+1, L, 2
       ! Construct the slices of the matrix of auxiliary functions
       do nn = 0, N-1
          AUX2(2*nn+1, 1) = rlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          AUX2(2*nn+2, 2) = plk(jj+1, nn+1)
       enddo
       ! Construct F2 = (F[l] * Aux[l,k])(3:4, :)
       do nn = 0, N-1
          ! Third column.
          F2(1, 2*nn+1) =  FL_coeff(3, jj+1) * rlk(jj+1, nn+1)
          F2(2, 2*nn+1) = -FL_coeff(6, jj+1) * rlk(jj+1, nn+1)
          ! Fourth column.
          F2(1, 2*nn+2) = FL_coeff(6, jj+1) * plk(jj+1, nn+1)
          F2(2, 2*nn+2) = FL_coeff(4, jj+1) * plk(jj+1, nn+1)
       enddo
       ! j-k is odd.
       ! B += 2 * Aux[j,k] * F[j] * Aux[j,k]^T
       call dgemm('N', 'N', nrows, nrows, 2, 2.0_dp, AUX2, nrows, &
            F2, 2, 1.0_dp, B, nrows)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = 0.5 * albedo (Zpos_BR + Zneg_BR) * W
    ! and
    !  B = 0.5 * albedo (Zpos_BR - Zneg_BR) * W.
    !
    do nn = 0, N-1
       tmp1 = 0.5_dp * albedo * wt(nn+1)
       A(:, 2*nn+1 : 2*nn+2) = tmp1 * A(:, 2*nn+1 : 2*nn+2)
       B(:, 2*nn+1 : 2*nn+2) = tmp1 * B(:, 2*nn+1 : 2*nn+2)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = -I + 0.5 * albedo (Zpos_BR + Zneg_BR) * W
    ! and
    !  B = -I + 0.5 * albedo (Zpos_BR - Zneg_BR) * W.
    !
    do nn = 1, nrows
       A(nn, nn) = A(nn, nn) - 1.0_dp
       B(nn, nn) = B(nn, nn) - 1.0_dp
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  A = (-I + 0.5 * albedo (Zpos_BR - Zneg_BR) * W) * M^{-1}
    ! and
    !  B = (-I + 0.5 * albedo (Zpos_BR + Zneg_BR) * W) * M^{-1}.
    !
    do nn = 0, N-1
       A(2*nn+1 : 2*nn+2, :) = A(2*nn+1 : 2*nn+2, :) / mu(nn+1)
       B(2*nn+1 : 2*nn+2, :) = B(2*nn+1 : 2*nn+2, :) / mu(nn+1)
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute
    !  AB = A * B.
    call dgemm('N', 'N', nrows, nrows, nrows, 1.0_dp, A, nrows, B, nrows, &
         0.0_dp, AB, nrows)
    !---------------------------------------------------------------------------
  end subroutine vrte0_uv


  subroutine source(k, L, N, albedo, mu, plk, rlk, tlk, FL_coeff, src1, src2)
    !!--------------------------------------------------------------------------
    !! Name: source
    !!  Compute the discretized source terms for the reduced problem
    !!  resulting from the direct portion of the Stokes vector.
    !!
    !! Synopsis
    !!   subroutine source(k, L, N, albedo, mu, I0, plk, rlk, tlk, FL_coeff,
    !!                     src1, src2)
    !!
    !!     INPUT
    !!       integer ------------- k, L, N
    !!       double precision ---- albedo, mu, plk(L+1, N+1), rlk(L+1, N+1),
    !!                             tlk(L+1, N+1), FL_coeff(6, L+1)
    !!
    !!     OUTPUT
    !!       double precision ---- src1(nrows, 2), src2(nrows, 2)
    !!
    !!
    !! Purpose
    !!  Compute the reduced source terms for the kth Fourier mode resulting
    !!  from the direct portion of the Stokes vector for a general source
    !!  term (e.g. do not multiply by I_inc at this point). The sources are
    !!  given by
    !!
    !!     src1 = 0.5 * albedo * M^{-1} * (Zpos_{0} - D34 * Zneg_{0})
    !!          = 0.5 * albedo * M^{-1} * sum[l=k,L] Aux[l,k]
    !!            * (1 - (-1)^{l-k} * D34) * F[l] * Aux[l,k](mu0)
    !!
    !!  and
    !!
    !!     src2 = 0.5 * albedo * M^{-1} * (Zpos_{0} + D34 * Zneg_{0})
    !!          = 0.5 * albedo * M^{-1} * sum[l=k,L] Aux[l,k]
    !!            * (1 + (-1)^{l-k} * D34) * F[l] * Aux[l,k](mu0).
    !!
    !! Arguments
    !!  k (input) integer
    !!    The Fourier mode of the azimuth.
    !!
    !!  L (input) integer
    !!    The order of the expansion in auxiliary rotation functions.
    !!
    !!  N (input) integer
    !!    The number of ordinate quadrature points.
    !!
    !!  albedo (input) double precision
    !!    The single scattering albedo.
    !!
    !!  mu (input) double precision, array(N)
    !!   Quadrature of ordinate directions.
    !!
    !!  plk (output) double precision, array(L+1, N)
    !!   Auxiliary functions as defined by Siewert et.al. 2000.
    !!
    !!  rlk (output) double precision, array(L+1, N)
    !!   Auxiliary functions as defined by Siewert et.al. 2000.
    !!
    !!  tlk (output) double precision, array(L+1, N)
    !!   Auxiliary functions as defined by Siewert et.al. 2000.
    !!
    !!  FL_coeff (input) double precision, array(L+1, 4*N, 4*N)
    !!    Expansion coefficients of the scattering matrix.
    !!
    !!  src1 (output) double precision, array(4*N, 4*N)
    !!    Source term given by
    !!     src1 = 0.5 * albedo * (Zpos_{0} - D34 * Zneg_{0})
    !!
    !!  src2 (output) double precision array, dimension (4*N, 4*N)
    !!    Source term govem by
    !!     src2 = 0.5 * albedo * (Zpos_{0} + D34 * Zneg_{0})
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8), intent(in) :: k, L, N
    real(dp), intent(in) :: albedo
    real(dp), intent(in), dimension(:) :: mu
    real(dp), intent(in), dimension(:, :) :: FL_coeff ! dim(6, L+1)
    real(dp), intent(in), dimension(:, :) :: plk, rlk, tlk ! dim(L+1, N+1)
    ! OUTPUT
    real(dp), intent(inout), dimension(:, :) ::  src1, src2 ! dim(4N, 4)
    ! Workspace
    integer(8) :: jj, nn
    real(dp) :: tmp1

    !---------------------------------------------------------------------------
    ! Compute discretized vector associated with the source term.
    ! ==========================================================
    !   S[a,l,k] = (albedo/2) * F_l * P[l,k](mu0) * I_init
    !
    src1 = 0.0_dp
    src2 = 0.0_dp
    do jj = k, L, 2
       !------------------------------------------------------------------------
       ! j-k is even,
       !   src1 += Aux[l,k] * diag(0,0,1,1) * F[l] * Aux0
       !
       do nn = 0, N-1
          ! I0 = [0, 1, 0, 0]
          src1(4*nn+2, 2) = src1(4*nn+2, 2) &
               + FL_coeff(3, jj+1) * tlk(jj+1, N+1) * tlk(jj+1, nn+1)
          src1(4*nn+3, 2) = src1(4*nn+3, 2) &
               - FL_coeff(3, jj+1) * tlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src1(4*nn+4, 2) = src1(4*nn+4, 2) &
               + FL_coeff(6, jj+1) * tlk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 0, 1, 0]
          src1(4*nn+2, 3) = src1(4*nn+2, 3) &
               - FL_coeff(3, jj+1) * rlk(jj+1, N+1) * tlk(jj+1, nn+1)
          src1(4*nn+3, 3) = src1(4*nn+3, 3) &
               + FL_coeff(3, jj+1) * rlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src1(4*nn+4, 3) = src1(4*nn+4, 3) &
               - FL_coeff(6, jj+1) * rlk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 0, 0, 1]
          src1(4*nn+2, 4) = src1(4*nn+2, 4) &
               - FL_coeff(6, jj+1) * plk(jj+1, N+1) * tlk(jj+1, nn+1)
          src1(4*nn+3, 4) = src1(4*nn+3, 4) &
               + FL_coeff(6, jj+1) * plk(jj+1, N+1) * rlk(jj+1, nn+1)
          src1(4*nn+4, 4) =  src1(4*nn+4, 4) &
               + FL_coeff(4, jj+1) * plk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       !------------------------------------------------------------------------

       !------------------------------------------------------------------------
       ! j-k is even,
       !   src2 += Aux[l,k] * diag(1,1,0,0) * F[l] * Aux0.
       !
       do nn = 0, N-1
          ! I0 = [1, 0, 0, 0]
          src2(4*nn+1, 1) = src2(4*nn+1, 1) &
               + FL_coeff(1, jj+1) * plk(jj+1, N+1) * plk(jj+1, nn+1)
          src2(4*nn+2, 1) = src2(4*nn+2, 1) &
               + FL_coeff(5, jj+1) * plk(jj+1, N+1) * rlk(jj+1, nn+1)
          src2(4*nn+3, 1) = src2(4*nn+3, 1) &
               - FL_coeff(5, jj+1) * plk(jj+1, N+1) * tlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 1, 0, 0]
          src2(4*nn+1, 2) = src2(4*nn+1, 2) &
               + FL_coeff(5, jj+1) * rlk(jj+1, N+1) * plk(jj+1, nn+1)
          src2(4*nn+2, 2) = src2(4*nn+2, 2) &
               + FL_coeff(2, jj+1) * rlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src2(4*nn+3, 2) = src2(4*nn+3, 2) &
               - FL_coeff(2, jj+1) * rlk(jj+1, N+1) * tlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 0,1, 0]
          src2(4*nn+1, 3) = src2(4*nn+1, 3) &
               - FL_coeff(5, jj+1) * tlk(jj+1, N+1) * plk(jj+1, nn+1)
          src2(4*nn+2, 3) = src2(4*nn+2, 3) &
               - FL_coeff(2, jj+1) * tlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src2(4*nn+3, 3) = src2(4*nn+3, 3) &
               + FL_coeff(2, jj+1) * tlk(jj+1, N+1) * tlk(jj+1, nn+1)
       enddo
       !-----------------------------------------------------------------------
    enddo
    do jj = k+1, L, 2
       !------------------------------------------------------------------------
       ! j-k is odd,
       !   src1 += Aux[l,k] * diag(1,1,0,0) * F[l] * Aux0.
       !
       do nn = 0, N-1
          ! I0 = [1, 0, 0, 0]
          src1(4*nn+1, 1) = src1(4*nn+1, 1) &
               + FL_coeff(1, jj+1) * plk(jj+1, N+1) * plk(jj+1, nn+1)
          src1(4*nn+2, 1) = src1(4*nn+2, 1) &
               + FL_coeff(5, jj+1) * plk(jj+1, N+1) * rlk(jj+1, nn+1)
          src1(4*nn+3, 1) = src1(4*nn+3, 1) &
               - FL_coeff(5, jj+1) * plk(jj+1, N+1) * tlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 1, 0, 0]
          src1(4*nn+1, 2) = src1(4*nn+1, 2) &
               + FL_coeff(5, jj+1) * rlk(jj+1, N+1) * plk(jj+1, nn+1)
          src1(4*nn+2, 2) = src1(4*nn+2, 2) &
               + FL_coeff(2, jj+1) * rlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src1(4*nn+3, 2) = src1(4*nn+3, 2) &
               - FL_coeff(2, jj+1) * rlk(jj+1, N+1) * tlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 0, 1, 0]
          src1(4*nn+1, 3) = src1(4*nn+1, 3) &
               - FL_coeff(5, jj+1) * tlk(jj+1, N+1) * plk(jj+1, nn+1)
          src1(4*nn+2, 3) = src1(4*nn+2, 3) &
               - FL_coeff(2, jj+1) * tlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src1(4*nn+3, 3) = src1(4*nn+3, 3) &
               + FL_coeff(2, jj+1) * tlk(jj+1, N+1) * tlk(jj+1, nn+1)
       enddo
       !------------------------------------------------------------------------

       !------------------------------------------------------------------------
       ! j-k is odd,
       !   src2 += Aux[l,k] * diag(0,0,1,1) * F[l] * Aux0.
       !
       do nn = 0, N-1
          ! I0 = [0, 1, 0, 0]
          src2(4*nn+2, 2) = src2(4*nn+2, 2) &
               + FL_coeff(3, jj+1) * tlk(jj+1, N+1) * tlk(jj+1, nn+1)
          src2(4*nn+3, 2) = src2(4*nn+3, 2) &
               - FL_coeff(3, jj+1) * tlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src2(4*nn+4, 2) = src2(4*nn+4, 2) &
               + FL_coeff(6, jj+1) * tlk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 0, 1, 0]
          src2(4*nn+2, 3) = src2(4*nn+2, 3) &
               - FL_coeff(3, jj+1) * rlk(jj+1, N+1) * tlk(jj+1, nn+1)
          src2(4*nn+3, 3) = src2(4*nn+3, 3) &
               + FL_coeff(3, jj+1) * rlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src2(4*nn+4, 3) = src2(4*nn+4, 3) &
               - FL_coeff(6, jj+1) * rlk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 0, 0, 1]
          src2(4*nn+2, 4) = src2(4*nn+2, 4) &
               - FL_coeff(6, jj+1) * plk(jj+1, N+1) * tlk(jj+1, nn+1)
          src2(4*nn+3, 4) = src2(4*nn+3, 4) &
               + FL_coeff(6, jj+1) * plk(jj+1, N+1) * rlk(jj+1, nn+1)
          src2(4*nn+4, 4) =  src2(4*nn+4, 4) &
               + FL_coeff(4, jj+1) * plk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       !------------------------------------------------------------------------
    enddo

    !---------------------------------------------------------------------------
    ! Multiply by 0.5 * albedo and divide by mu.
    !
    do nn = 0, N-1
       tmp1 = 0.5_dp * albedo / mu(nn+1)
       src1(4*nn+1 : 4*nn+4, :) = tmp1 * src1(4*nn+1 : 4*nn+4, :)
       src2(4*nn+1 : 4*nn+4, :) = tmp1 * src2(4*nn+1 : 4*nn+4, :)
    enddo
    !---------------------------------------------------------------------------

  end subroutine source


  subroutine source0_iq(L, N, albedo, mu, plk, rlk, FL_coeff, src1, src2)
    !!--------------------------------------------------------------------------
    !! Name: source0_iq
    !!  Compute the I and Q parts of the discretized source terms for the
    !!  zeroth Fourier mode of the reduced problem resulting from the direct
    !!  portion of the Stokes vector.
    !!
    !! Synopsis
    !!  subroutine source0_iq(L, N, albedo, mu, plk, rlk, tlk, FL_coeff,
    !!                        src1, src2 )
    !!
    !!     INPUT
    !!       integer ------------- L, N
    !!       double precision ---- albedo, mu,  I0(4), plk(L+1, N+1),
    !!                             rlk(L+1, N+1), FL_coeff(6, L+1)
    !!
    !!     OUTPUT
    !!       double precision ---- src1(nrows, 2), src2(nrows, 2)
    !!
    !!
    !! Purpose
    !!  Compute the I and Q components of the reduced source term for the
    !!  zeroth Fourier mode resulting from the direct portion of the Stokes
    !!  vector for a general source term.
    !!
    !!  The source terms are given by
    !!
    !!       src1 = 0.5 * albedo * M^{-1} * (Zpos_{0}_UL - D34 * Zneg_{0}_UL)
    !!            = 0.5 * albedo * M^{-1} * sum[l=k,L] Aux[l,k]
    !!                  * (1 - (-1)^{l-k}) * F[l] * Aux[l,k](mu0)
    !!
    !!     and
    !!
    !!       src2 = 0.5 * albedo * M^{-1} * (Zpos_{0}_UL + D34 * Zneg_{0}_UL)
    !!            = 0.5 * albedo * M^{-1} * sum[l=k,L] Aux[l,k]
    !!                  * (1 + (-1)^{l-k}) * F[l] * Aux[l,k](mu0).
    !!
    !! Arguments
    !!  L (input) integer
    !!    The order of the expansion in auxiliary rotation functions.
    !!
    !!  N (input) integer
    !!    The number of ordinate quadrature points.
    !!
    !!  albedo (input) double precision
    !!    The single scattering albedo.
    !!
    !!  mu (input) double precision, array(N)
    !!   Quadrature of ordinate directions.
    !!
    !!  plk (output) double precision, array(L+1, N)
    !!   Auxiliary functions as defined by Siewert et.al. 2000.
    !!
    !!  rlk (output) double precision, array(L+1, N)
    !!   Auxiliary functions as defined by Siewert et.al. 2000.
    !!
    !!  tlk (output) double precision, array(L+1, N)
    !!   Auxiliary functions as defined by Siewert et.al. 2000.
    !!
    !!  FL_coeff (input) double precision, array(L+1, 4*N, 4*N)
    !!    Expansion coefficients of the scattering matrix.
    !!
    !!  src1 (output) double precision, array(2*N, 2)
    !!    I and Q components of the source term given by
    !!     src1 = 0.5 * albedo * M^{-1} * (Zpos_{0}_UL - Zneg_{0}_UL)
    !!
    !!  src2 (output) double precision array, dimension (2*N, 2)
    !!    I and Q components of the source term given by
    !!     src2 = 0.5 * albedo * M^{-1} * (Zpos_{0}_UL + Zneg_{0}_UL)
    !!--------------------------------------------------------------------------
    implicit none
    integer(8), parameter :: k = 0
    ! INPUT
    integer(8), intent(in) :: L, N
    real(dp), intent(in) :: albedo
    real(dp), intent(in), dimension(:) :: mu ! dim(N)
    real(dp), intent(in), dimension(:, :) :: FL_coeff ! dim(6, L+1)
    real(dp), intent(in), dimension(:, :) :: plk, rlk ! dim(L+1, N+1)
    ! OUTPUT
    real(dp), intent(inout), dimension(:, :) ::  src1, src2 ! dim(2N, 2)
    ! Workspace
    integer(8) :: jj, nn
    real(dp) :: tmp1

    !---------------------------------------------------------------------------
    ! Compute discretized vector associated with the source term.
    ! ==========================================================
    !   S[a,l,k] = (albedo/2) * F_l * P[l,k](mu0) * I_init
    !
    src1 = 0.0_dp
    src2 = 0.0_dp
    do jj = k, L, 2
       !------------------------------------------------------------------------
       ! j-k is even,
       !   src2 += Aux[l,k] * diag(1,1,0,0) * F[l] * Aux0.
       !
       do nn = 0, N-1
          ! I0 = [1, 0, 0, 0]
          src2(2*nn+1, 1) = src2(2*nn+1, 1) &
               + FL_coeff(1, jj+1) * plk(jj+1, N+1) * plk(jj+1, nn+1)
          src2(2*nn+2, 1) = src2(2*nn+2, 1) &
               + FL_coeff(5, jj+1) * plk(jj+1, N+1) * rlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 1, 0, 0]
          src2(2*nn+1, 2) = src2(2*nn+1, 2) &
               + FL_coeff(5, jj+1) * rlk(jj+1, N+1) * plk(jj+1, nn+1)
          src2(2*nn+2, 2) = src2(2*nn+2, 2) &
               + FL_coeff(2, jj+1) * rlk(jj+1, N+1) * rlk(jj+1, nn+1)
       enddo
       !-----------------------------------------------------------------------
    enddo
    do jj = k+1, L, 2
       !------------------------------------------------------------------------
       ! j-k is odd,
       !   src1 += Aux[l,k] * diag(1,1,0,0) * F[l] * Aux0.
       !
       do nn = 0, N-1
          ! I0 = [1, 0, 0, 0]
          src1(2*nn+1, 1) = src1(2*nn+1, 1) &
               + FL_coeff(1, jj+1) * plk(jj+1, N+1) * plk(jj+1, nn+1)
          src1(2*nn+2, 1) = src1(2*nn+2, 1) &
               + FL_coeff(5, jj+1) * plk(jj+1, N+1) * rlk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 1, 0, 0]
          src1(2*nn+1, 2) = src1(2*nn+1, 2) &
               + FL_coeff(5, jj+1) * rlk(jj+1, N+1) * plk(jj+1, nn+1)
          src1(2*nn+2, 2) = src1(2*nn+2, 2) &
               + FL_coeff(2, jj+1) * rlk(jj+1, N+1) * rlk(jj+1, nn+1)
       enddo
       !------------------------------------------------------------------------
    enddo

    !---------------------------------------------------------------------------
    ! Multiply by 0.5 * albedo and divide by mu.
    !
    do nn = 0, N-1
       tmp1 = 0.5_dp * albedo / mu(nn+1)
       src1(2*nn+1 : 2*nn+2, :) = tmp1 * src1(2*nn+1 : 2*nn+2, :)
       src2(2*nn+1 : 2*nn+2, :) = tmp1 * src2(2*nn+1 : 2*nn+2, :)
    enddo
    !---------------------------------------------------------------------------

  end subroutine source0_iq


  subroutine source0_uv(L, N, albedo, mu, plk, rlk, FL_coeff, src1, src2)
    !!--------------------------------------------------------------------------
    !! Name: source0_uv
    !!  Compute the U and V parts of the discretized source terms for the
    !!  zeroth Fourier mode of the reduced problem resulting from the direct
    !!  portion of the Stokes vector.
    !!
    !! Synopsis
    !!  subroutine source0_iq(L, N, albedo, mu, I0, plk, rlk, tlk, FL_coeff,
    !!                        src1, src2 )
    !!
    !!     INPUT
    !!       integer ------------- L, N
    !!       double precision ---- albedo, mu0,  I0(4), plk(L+1, N+1),
    !!                             rlk(L+1, N+1), FL_coeff(6, L+1)
    !!
    !!     OUTPUT
    !!       double precision ---- src1(nrows, 2), src2(nrows, 2)
    !!
    !!
    !! Purpose
    !!  Compute the I and Q components of the reduced source term for the
    !!  zeroth Fourier mode resulting from the direct portion of the Stokes
    !!  vector for a general source term.
    !!
    !!  The source terms are given by
    !!
    !!       src1 = 0.5 * albedo * M^{-1} * (Zpos_{0}_UL - D34 * Zneg_{0}_UL)
    !!                 = 0.5 * albedo * M^{-1} * sum[l=k,L] Aux[l,k]
    !!                  * (1 - (-1)^{l-k}) * F[l] * Aux[l,k](mu0)
    !!
    !!     and
    !!
    !!       src2 = 0.5 * albedo * (Zpos_{0}_UL + D34 * Zneg_{0}_UL)
    !!            = 0.5 * albedo * sum[l=k,L] Aux[l,k]
    !!                  * (1 + (-1)^{l-k}) * F[l] * Aux[l,k](mu0).
    !!
    !! Arguments
    !!   L (input) integer
    !!     The order of the expansion in auxiliary rotation functions.
    !!
    !!   N (input) integer
    !!     The number of ordinate quadrature points.
    !!
    !!   albedo (input) double precision
    !!     The number of rows in the discretized system, nrows = 4*N.
    !!
    !!   mu (input) double precision, array(N)
    !!     Incident ordinate direction.
    !!
    !!   FL_coeff (input) double precision, array(L+1, 4*N, 4*N)
    !!     Expansion coefficients of the scattering matrix.
    !!
    !!   src1 (output) double precision, array(2*N, 2)
    !!     I and Q components of the source term given by
    !!      src1 = 0.5 * albedo * (Zpos_{0}_UL - Zneg_{0}_UL)
    !!
    !!   src2 (output) double precision array, dimension (2*N, 2)
    !!     I and Q components of the source term given by
    !!      src2 = 0.5 * albedo * (Zpos_{0}_UL + Zneg_{0}_UL)
    !!--------------------------------------------------------------------------
    implicit none
    integer(8), parameter :: k = 0
    ! INPUT
    integer(8), intent(in) :: L, N
    real(dp), intent(in) :: albedo
    real(dp), intent(in), dimension(:, :) :: FL_coeff ! dim(6, L+1)
    real(dp), intent(in), dimension(:) :: mu ! dim(N)
    real(dp), intent(in), dimension(:, :) :: plk, rlk ! dim(N+1, L+1)
    ! OUTPUT
    real(dp), intent(inout), dimension(:, :) ::  src1, src2 ! dim(4N, 4)
    ! Workspace
    integer(8) :: jj, nn
    real(dp) :: tmp1

    src1 = 0.0_dp
    src2 = 0.0_dp
    do jj = k, L, 2
       !------------------------------------------------------------------------
       ! j-k is even,
       !   src1 += Aux[l,k] * diag(0,0,1,1) * F[l] * Aux0
       !
       do nn = 0, N-1
          ! I0 = [0, 0, 1, 0]
          src1(2*nn+1, 1) = src1(2*nn+1, 1) &
               + FL_coeff(3, jj+1) * rlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src1(2*nn+2, 1) = src1(2*nn+2, 1) &
               - FL_coeff(6, jj+1) * rlk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 0, 0, 1]
          src1(2*nn+1, 2) = src1(2*nn+1, 2) &
               + FL_coeff(6, jj+1) * plk(jj+1, N+1) * rlk(jj+1, nn+1)
          src1(2*nn+2, 2) =  src1(2*nn+2, 2) &
               + FL_coeff(4, jj+1) * plk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       !------------------------------------------------------------------------
    enddo

    do jj = k+1, L, 2
       !------------------------------------------------------------------------
       ! j-k is odd,
       !   src2 += Aux[l,k] * diag(0,0,1,1) * F[l] * Aux0.
       !
       do nn = 0, N-1
          ! I0 = [0, 0, 1, 0]
          src2(2*nn+1, 1) = src2(2*nn+1, 1) &
               + FL_coeff(3, jj+1) * rlk(jj+1, N+1) * rlk(jj+1, nn+1)
          src2(2*nn+2, 1) = src2(2*nn+2, 1) &
               - FL_coeff(6, jj+1) * rlk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       do nn = 0, N-1
          ! I0 = [0, 0, 0, 1]
          src2(2*nn+1, 2) = src2(2*nn+1, 2) &
               + FL_coeff(6, jj+1) * plk(jj+1, N+1) * rlk(jj+1, nn+1)
          src2(2*nn+2, 2) =  src2(2*nn+2, 2) &
               + FL_coeff(4, jj+1) * plk(jj+1, N+1) * plk(jj+1, nn+1)
       enddo
       !------------------------------------------------------------------------
    enddo

    !---------------------------------------------------------------------------
    ! Multiply by 0.5 * albedo and divide by mu.
    !
    do nn = 0, N-1
       tmp1 = 0.5_dp * albedo / mu(nn+1)
       src1(2*nn+1 : 2*nn+2, :) = tmp1 * src1(2*nn+1 : 2*nn+2, :)
       src2(2*nn+1 : 2*nn+2, :) = tmp1 * src2(2*nn+1 : 2*nn+2, :)
    enddo
    !---------------------------------------------------------------------------
  end subroutine source0_uv


end module discretize
