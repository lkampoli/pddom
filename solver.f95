module solver
  !-----------------------------------------------------------------------------
  ! Module: solver
  !  This module contains subroutines that compute the general solution of a
  !  given Fourier mode of the discrete radiative transport equation (DRTE).
  !
  ! Subroutines:
  !  vdisord (public)
  !    Compute the general solution of a Fourier mode of the DRTE.
  !
  !  particular (public)
  !    Compute the particular solution of a Fourier mode of the DRTE.
  !
  !  eigensystem (public)
  !    Compute the eigenvalues and eigenvectors of a Fourier mode of the DRTE.
  !
  !  source (private)
  !    Compute the discretized source terms for the reduced problem
  !    resulting from the direct portion of the Stokes vector.
  !
  !  vdisord0_iq (public)
  !vdisord0_iq
  !!   Compute the general solution for the upper-right block of the zeroth
  !!   Fourier mode of the DRTE in the form of planewave solutions.
  !  eigensystem0_iq (public)
  !    Compute the eigenvalues and eigenvectors of the upper-right block of the
  !    zeroth Fourier mode of the reduced DRTE.
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
  public :: vdisord, vdisord0_iq, eigensystem

contains


  subroutine vdisord( nrows, nsource, mu0, A, B, AB, tol, &
       Psi_p, Gamma_p, Psi_h, Gamma_h, lambda, nreal, nzero )
    !!--------------------------------------------------------------------------
    !! NAME
    !!   vdisord: Compute the general solution for the reduced discretized
    !!   radiative transport equation in the form of planewave solutions.
    !!
    !! SYNOPSIS
    !!   subroutine vdisord( N, nrows, nsource, mu, mu0, kappa, A, B, AB, tol,
    !!       Psi_p, Gamma_p, Psi_h, Gamma_h, lambda, nreal, nzero )
    !!
    !!   Integer ------------- N, nrows, nsource, nreal, nzero
    !!   Double precision ---- mu( N ), mu0, kappa, A( nrows, nrows ),
    !!                      B( nrows, nrows ), AB( nrows, nrows ), tol,
    !!                      Psi_p( nrows, nsource ), Gamma_h( nrows, nsource )
    !!                      Psi_h( nrows, nrows ), Gamma_h( nrows, nrows ),
    !!                      lambda( nrows )
    !!
    !! PURPOSE
    !!   Compute the general solution of the first order system of differential
    !!   equations
    !!
    !!         Psi'(tau) = A * Gamma(tau) + S1 * exp(-tau/mu0),
    !!       Gamma'(tau) = B *   Psi(tau) + S2 * exp(-tau/mu0).
    !!
    !!    where A and B are linear combinations of the discretized radiative
    !!    transport operator given by
    !!
    !!       A = M^{-1} * (-I + 0.5 * albedo * (Zpos - D34 * Zneg) * W),
    !!       B = M^{-1} * (-I + 0.5 * albedo * (Zpos + D34 * Zneg) * W),
    !!
    !!    and S1 and S2 are the source terms
    !!
    !!       S1 = 0.5 * M^{-1} * (S_pos - D34 * S_neg) / mu0,
    !!       S2 = 0.5 * M^{-1} * (S_pos + D34 * S_neg) / mu0.
    !!
    !! ARGUMENTS
    !!   N (input) integer
    !!     The number of ordinate quadrature points.
    !!
    !!   nsource (input) integer
    !!     The number of source vectors.
    !!
    !!   mu (input) double precision, array(N)
    !!     Quadrature of ordinate directions.
    !!
    !!   mu0 (input) double precision
    !!     Ordinate direction of incident light.
    !!
    !!   A (input) double precision, array(nrows, nrows)
    !!     Discretization of the modified radiative transport operator
    !!      A = (-I + 0.5 * albedo (Zpos - D34 * Zneg) * W) M^{-1}.
    !!
    !!   B (input) double precision, array(nrows, nrows)
    !!     Discretization of the modified radiative transport operator
    !!       B = (-I + 0.5 * albedo (Zpos + D34 * Zneg) * W) M^{-1}.
    !!
    !!   AB (input) double precision, array(nrows, nrows)
    !!     Matrix system AB = A*B.
    !!
    !!   tol (input)
    !!     Tolerance for numerical precision.
    !!
    !!   Psi_p (input/output) double precision, array( nrows, nrows )
    !!     On entry, the source term S1.
    !!     On exit, the particular solution Psi_p to the discrete system.
    !!
    !!   Gamma_p (input/output) double precision, array( nrows, nrows )
    !!     On entry, the source term S2.
    !!     On exit, the particular solution Gamma_p to the discrete system.
    !!
    !!   Psi_h (output) double precision, array ( nrows, nrows )
    !!     A matrix with row vectors corresponding to the planewave solutions.
    !!     Psi_h are the eigenvectors of the matrix (AB).
    !!
    !!   Gamma_h (output) double precision, array ( nrows, nrows )
    !!     A matrix with row vectors corresponding to the planewave solutions.
    !!     Gamma_h are the eigenvectors of the matrix (BA).
    !!
    !!   lambda (output) double precision, array ( nrows )
    !!     An array of eigenvalues/decay of the planewave solutions.
    !!
    !!   nreal (output) integer
    !!     The number of strictly real eigenvalues.
    !!
    !!   nzero (output) integer
    !!     The number of identically zero eigenvalues.
    !!---------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8), intent(in) :: nrows, nsource
    real(dp), intent(in) :: mu0
    real(dp), intent(in), dimension(:, :) :: A, B, AB
    real(dp), intent(in) :: tol
    ! IN/OUT
    real(dp), intent(inout), dimension(:, :) :: Psi_p, Gamma_p ! dim(n, nsource)
    ! OUTPUT
    real(dp), intent(out), dimension(:, :) :: Psi_h, Gamma_h ! dim(4N, 4N)
    real(dp), intent(out), dimension(:) :: lambda ! dim(4N)
    integer(8), intent(out) :: nreal, nzero
    ! Workspace
    real(dp), dimension(nrows)  :: wr, wi
    real(dp), dimension(nrows, nrows) :: workspace ! workspace
    real(dp) :: cond = 10D-7
    nzero = 0
    !
    ! 1. Store Psi_p = -1/mu0 * S1 + A * S2.
    !
    call dgemm('N', 'N', nrows, nsource, nrows, 1.0_dp, A, nrows, &
         Gamma_p, nrows, -1.0_dp / mu0, Psi_p, nrows)
    !
    ! 2. Compute spectrum of A and B.
    !
    call eigensystem(nrows, B, AB, workspace, Psi_h, Gamma_h, wr, wi, lambda, &
         nreal, tol)
    call particular(nrows, nreal, nsource, mu0, B, AB, wr, wi, Psi_h, &
         Psi_p, Gamma_p, workspace, cond)
  end subroutine vdisord


  subroutine eigensystem(nrows, B, AB, workspace, Psi_h, Gamma_h, wr, wi, &
       lambda, nreal, tol)
    !!--------------------------------------------------------------------------
    !! NAME: eigensystem
    !!  Compute the eigenvalues and eigenvectors of a Fourier mode of the DRTE.
    !!
    !! SYNOPSIS
    !!     subroutine eigensystem( N, nrows, mu, A, B, AB, Psi_h, Gamma_h, wr,
    !!         wi, lambda, nreal, nzero, tol )
    !!
    !!     Integer ------------- N, nrows, nreal, nzero
    !!     Double precision ----
    !!
    !! PURPOSE
    !!     Compute the eigenvalues and eigenvectors of the system
    !!
    !!        -lambda(j) * M *   Psi_h(j) = A * M * Gamma(j),
    !!        -lambda(j) * M * Gamma_h(j) = B * M *   Psi(j).
    !!
    !!      Return the eigenvalues lambda so that:
    !!        1. If there is a zero eigenvalue it is stored first.
    !!        2. The first nreal eigenvalues are strictly real.
    !!
    !!
    !! PARAMETERS
    !!     N (input) integer
    !!       Number of quadrature points.
    !!
    !!     nrows (input) integer
    !!       Dimension of the matrices A, B, and AB.
    !!
    !!     mu (input) real, double precision, array( N )
    !!       Quadrature points of the ordinate directions.
    !!
    !!     A (input) real, double precision, array(nrows, nrows)
    !!       Discretization of the modified radiative transport operator
    !!        A = (-I + 0.5 * albedo (Zpos - D34 * Zneg) * W) M^{-1}.
    !!
    !!     B (input) real, double precision, array(nrows, nrows)
    !!       Discretization of the modified radiative transport operator
    !!         B = (-I + 0.5 * albedo (Zpos + D34 * Zneg) * W) M^{-1}.
    !!
    !!     AB (input) real, double precision, array(nrows, nrows)
    !!       Matrix system AB = A*B.
    !!
    !!     Psi_h (output) real, double precision, array ( nrows, nrows )
    !!       The rows of Psi_h are the eigenvectors of the matrix AB. The
    !!       eigenvectors are stored so their sorting matches the eigenvalues
    !!       wr and wi, and if the jth and (j+1)th eigenvalue correspond to a
    !!       conjugate pair, then the corresponding eigenvector is stored so
    !!              Eigenvector = Psi_h(j) +/- i * Psi_h(j+1).
    !!
    !!     Gamma_h (output) real, double precision, array ( nrows, nrows )
    !!       The rows of Gamma_h are the eigenvectors of the matrix BA stored
    !!       in the same order as Psi_h.
    !!
    !!     wr (output) real, double precision, array ( nrows )
    !!     wi (output) real, double precision, array ( nrows )
    !!       wr and wi contain the real and imaginary parts, respectively, of
    !!       the eigenvalues of AB. The eigenvalues are sorted so that
    !!
    !!     lambda (output) real, double precision, array ( nrows )
    !!       Eigenvalues of [[0, A], [B, 0]] corresponding to the positive real
    !!       branch. Eigenvalues are stored so that if there is a zero
    !!       eigenvalue it is stored first, followed by the remaining strictly
    !!       real eigenvalues, and finally the complex eigenvalues which are
    !!       stored so complex conjugates occur in consecutive pairs.
    !!
    !!     nreal (output) integer
    !!       Number of strictly real eigenvalues.
    !!
    !!     tol (input) real, double precision
    !!       Tolerance for determining if an eigenvalue is zero or non-zero.
    !!
    !! ---
    !! Written by Julia CLark.
    !! University of California, Merced.
    !! Contact: jclark@ucmerced.edu
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8), intent(in) :: nrows
    real(dp), intent(in), dimension(:, :) :: B, AB
    real(dp), intent(in) :: tol
    ! Workspace
    real(dp), intent(inout), dimension(:, :) :: workspace
    real(dp), intent(out), dimension(:, :) :: Psi_h, Gamma_h
    real(dp), intent(out), dimension(:)  :: wr, wi, lambda
    integer(8), intent(out) :: nreal
    ! Workspace: DGEEVX
    integer(8) :: info, lwork
    real(dp), dimension(1) :: qwork
    real(dp), allocatable, dimension(:) :: work
    integer(8) :: ihi, ilo
    integer(8), dimension(2*nrows-2) :: iwork
    real(dp), dimension(nrows) :: rconde, rcondv
    real(dp), dimension(nrows) :: scale
    real(dp), dimension(nrows, nrows) :: Psi_adj
    real(dp) :: abnrm
    ! Preserve AB
    ! Workspace
    integer(8) :: ncomplex
    integer(8), dimension(nrows) :: eigindx
    integer(8) :: mm, ii
    real(dp) :: theta, rho
    real(dp), dimension(nrows) :: v1, v2
    real(dp) :: tmp1

    !---------------------------------------------------------------------------
    ! Compute the spectrum of the matrix AB.
    !
    workspace = AB
    ! Query work size.
    info = 0
    lwork = -1
    call dgeevx('B','V', 'V', 'B',&
         nrows, workspace, nrows, wr, wi, Psi_adj, nrows, Psi_h, nrows, &
         ilo, ihi, scale, abnrm, rconde, rcondv, &
         qwork, lwork, iwork, info)
    ! Solve for eigenvalues and eigenvectors.
    lwork = int(qwork(1), 8)
    allocate(work(lwork))
    call dgeevx('B','V', 'V', 'B',&
         nrows, workspace, nrows, wr, wi, Psi_adj, nrows, Psi_h, nrows, &
         ilo, ihi, scale, abnrm, rconde, rcondv, &
         work, lwork, iwork, info)
    if (info > 0) then
       write (*, '(A)') &
            "> Error in __solver_MOD_spectrum. "
       write (*, '(A)') &
            "> Lapack subroutine DGEEVX failed to compute all eigenvalues."
       write (*, '(A, I6)') "DGEEVX, info = ", info
       stop
    endif
    deallocate(work)
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Index Sort
    !   The planewave solutions are of the form
    !
    !      Psi(j) * exp(- lambda * z)
    !
    !   where the eigenvalues, lambda, are given by
    !
    !      lambda = sqrt( wr + i * wi )
    !             = sqrt((sqrt(wr^2 + wi^2 ) + wr) / 2)
    !               +/- i * sqrt((sqrt(wr^2 + wi^2) - wr) / 2)
    !
    !   Define:
    !      nreal = # of eigenvalues lambda that are strictly real,
    !
    !   and sort wr, wi, and Psi_h so that
    !
    !    1. The first nreal eigenvalues are strictly real.
    !    2. For complex conjugate pairs, store the eigenvectors in the form
    !          Psi = Psi(j) + i * Psi(j+1).
    !
    !
    ! *Stop Procedure:
    !    The program produces an error and halts if there is an eigenvalue with
    !    non-positive real part. This means this subroutine cannot handle the
    !    computation of the zeroth Fourier mode when albedo = 1.
    !
    nreal = nrows
    ncomplex = 0
    if (minval(wr) < tol) then
       write (*, '(A)') &
            "> Error in __solver_MOD_eigensort."
       write (*, '(A, I2, A)') &
            "> Computed eigenvalue with non-positive real part."
       stop
    endif
    !
    ! Resort any complex eigenvectors to the end of the array.
    !
    !  * Since there is no strictly imaginary eigenvalues, lambda has a
    !    nonzero imaginary part if and only if wi is nonzero.
    !
    do mm = 1, nrows
       eigindx(mm) = mm
    enddo
    if (maxval(wi) > tol) then
       mm = nrows
       do while (mm > 0)
          if (abs(wi(mm)) > tol) then
             ii = eigindx(mm)
             eigindx(mm) = eigindx(nrows-ncomplex)
             eigindx(nrows-ncomplex) = ii
             ncomplex = ncomplex + 1
             nreal = nreal - 1
          endif
          mm = mm -1
       enddo
    endif
    !
    ! Sort the eigenvalues and eigenvectors according to the computed indices.
    !
    wr = wr(eigindx)
    wi = wi(eigindx)
    workspace = Psi_h
    do mm = 1, nrows
       Psi_h(:, mm) = workspace(:, eigindx(mm))
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute lambda from wr and wi.
    !
    ! Complex eigenvalue with nonzero real part.
    !   lambda = rho * (cos(theta) + i * sin(theta))
    ! where
    !   (rho * exp(i(theta))^2 = wr + i * wi.
    !
    lambda(1 : nreal) = sqrt(wr(1 : nreal))
    do mm = nreal + 1, nrows, 2
       rho = (wr(mm) ** 2 + wi(mm) ** 2) ** 0.25_dp
       theta = 0.5_dp * atan2(wi(mm), wr(mm))
       lambda(mm) = rho * cos(theta)
       lambda(mm+1) = rho * sin(theta)
    enddo
    !--------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute (M * Gamma_h).
    !
    !   M * Gamma_h(j) = - B * M * Psi_h(j) / lambda(j),
    !
    !   *Note currently storing (M * Psi_h) in Psi_h.
    !
    call dgemm('N', 'N', nrows, nrows, nrows, -1.0_dp, B, nrows, Psi_h, nrows, &
         0.0_dp, Gamma_h, nrows) ! Set V = - B * U.
    do mm = 1, nreal
       Gamma_h(:, mm) = Gamma_h(:, mm) / lambda(mm)
    enddo
    do mm = nreal+1, nrows, 2
       !
       !   Re(Gamma_h) = -(Re[lambda] * (B * Re[Psi_h])
       !                    + Im[lambda] * (B * Im[Psi_h])) / |lambda|^2,
       !
       !   Im(Gamma_h) = -(Re[lambda] * (B * Im[Psi_h])
       !                   - Im[lambda] * (B * Re[Psi_h])) / |lambda|^2.
       !
       tmp1 = lambda(mm) ** 2 + lambda(mm+1) ** 2
       v1 = Gamma_h(:, mm)
       v2 = Gamma_h(:, mm+1)
       Gamma_h(:, mm) = (lambda(mm) * v1 + lambda(mm+1) * v2) / tmp1
       Gamma_h(:, mm+1) = (lambda(mm) * v2 - lambda(mm+1) * v1) / tmp1
    enddo
    !---------------------------------------------------------------------------
  end subroutine eigensystem


  subroutine particular(nrows, nreal, nsource, mu0, B, AB, wr, wi, Psi_h, &
       Psi_p, Gamma_p, workspace, cond)
    !!--------------------------------------------------------------------------
    !! NAME: particular
    !!  Compute the particular solution of a Fourier mode of the DRTE.
    !!
    !! SYNOPSIS
    !!     subroutine particular( nrows, nreal, mu0, AB, Psi_h, Psi_p,
    !!         workspace, cond )
    !!
    !!     Integer ------------- nrows, nreal
    !!     Double precision ---- mu0, AB(nrows, nrows), wr(nrows),
    !!                        Psi_h(nrows, nrows),
    !!                        Psi_p(nrows 2), workspace(nrows, nrows), cond
    !!
    !! PURPOSE
    !!   Compute the particular solution of the modified discritzed vRTE
    !!
    !!        Psi'(tau) = B * Gamma'(tau) + S1 * exp(-tau/mu0),
    !!        Gamma'(tau) = A * Psi'(tau) + S2 * exp(-tau/mu0).
    !!
    !!    The subroutine first solves the reduced second-order system
    !!
    !!        Psi_p''(tau) - AB * Psi_p(tau) = S * exp(-tau / mu0)
    !!
    !!    using the ansatz
    !!
    !!        Psi_p(tau) = Psi_p * exp(-tau/mu0)
    !!
    !!    and computes Gamma_p(tau) by using
    !!
    !!        Gamma_p(tau) = -mu0 * (A * Psi_p + S2) * exp(-tau/mu0).
    !!
    !! ARGUMENTS
    !!     nrows (input) integer
    !!       The number of equations in the system of equations.
    !!
    !!     nreal (input) integer
    !!       The number of real eigenvalues (1 < nreal <= nrows).
    !!
    !!     nsource (input) integer
    !!      The number of independent source terms.
    !!
    !!     mu0 (input) real, double precision
    !!      The incident polar direction of the light on the slab.
    !!
    !!     B ------- real, double precision, array(nrows, nrows)
    !!      Matrix of reduced discretized vRTE such that
    !!        Gamma'(tau) = B * Psi(tau) + S2 * exp(-tau/mu0).
    !!
    !!     AB ------- real, double precision, array(nrows, nrows)
    !!      Matrix of reduced discretized vRTE, AB = A *B.
    !!
    !!     wr ------- real, double precision, array(nrows)
    !!      Real eigenvalues of the reduced discretized vRTE, AB.
    !!
    !!     Psi_h ---- real, double precision, array(nrows, nrows)
    !!      Eigenvectors of the reduced discretzied vRTE, AB.
    !!
    !!     Psi_p (input/output) double precision, array( nrows, nrows )
    !!       On entry, the source term S to the reduced system,
    !!          S = - S1 / mu0 + A * S2.
    !!       On exit, the particular solution Psi_p to the discrete system.
    !!
    !!     Gamma_p (input/output) double precision, array( nrows, nrows )
    !!       On entry, the source term S2.
    !!       On exit, the particular solution Gamma_p to the discrete system.
    !!
    !!     workspace --- real, double precision, array(nrows, nrows)
    !!       Workspace for lapack routine dgesv.
    !!
    !!     cond ----- real, double precision
    !!      Conditioning number threshold for inverting (1/mu0^2 - lambda^2).
    !!
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8) :: nrows, nreal, nsource
    real(dp), intent(in) :: mu0
    real(dp), intent(in), dimension(:, :) :: Psi_h ! dim(4N, 4N)
    real(dp), intent(in), dimension(:, :) :: B, AB ! dim(4N, 4N)
    real(dp), intent(in), dimension(:) :: wr, wi ! dim(4N)
    real(dp), intent(in) :: cond
    ! OUTPUT
    real(dp), intent(inout), dimension(:, :) :: Psi_p, Gamma_p ! dim(8*N, nsource)
    ! In/Out workspace
    real(dp), intent(inout), dimension(:, :) :: workspace ! dim(4N, 4N)
    ! Workspace
    integer, parameter :: lwmax = 100000
    integer :: info
    integer, dimension(nrows) :: IPIV
    real(dp) :: tmp0, tmp1, tmp2
    real(dp), dimension(nsource) :: x1, x2
    real(dp), dimension(nrows, nsource) :: u0, u1
    real(dp), dimension(nrows) :: diag
    integer(8) :: mm
    !---------------------------------------------------------------------------
    ! Solve the system of equations
    !
    !  Psi_p''(tau) - AB * Psi_p(tau) = S * exp(-tau / mu0).
    !
    ! using AB = Psi_h * Diag( wr + i * wi ) * Psi_h^{-1}.
    !
    !
    !  Psi_p = sum[lambda(j) =/= 1/mu0] c[j] * Psi_h(j) * exp(-tau/mu0)
    !         + sum[lambda(j) = 1/mu0] c[j] * Psi_h(j) * tau * exp(-tau/mu0)
    !
    diag = 1.0_dp / mu0 ** 2 - wr
    if (minval(abs(diag(1:nreal))) > cond) then
       !
       ! System is not degenerate: 1/mu0^2 is not an eigenvalue of AB.
       !
       ! Solve the system
       !
       !   (1/mu0^2 * Identity - AB) * Psi_p = Source
       !
       ! directly using LAPACK subroutine DGESV.
       !
       u0 = Psi_p
       workspace = -AB
       tmp0 = 1.0 / mu0 ** 2
       do mm = 1, nrows
          workspace(mm, mm) = workspace(mm, mm) + tmp0
       enddo
       info = 0
       call dgesv(nrows, nsource, workspace, nrows, IPIV, Psi_p, nrows, info)
       if (info > 0) then
          write (*, '(A)') "> Error in __solver_MOD_particular."
          write (*, '(A)') "> The system (1/mu0^2 * Idenity -  AB) is singular."
          write(*, '(A)') "> Program halted."
          stop
       endif
    else
       !
       ! System is degenerate: 1/mu0^2 is close to an eigenvalue of AB.
       !
       ! Solve for the particular solution
       !
       !  Psi_p = (tau * u0 + u1) * exp(-tau/mu0)
       !
       ! where
       !
       !   1/mu0^2 * u0 - AB * u0 = 0,
       !   1/mu0^2 * u1 - AB * u1 = S + u0 /mu0,
       !
       ! using decomposition of AB using
       !
       !    AB = Psi_h * diag( wr + i * wi ) * Psi_h^{-1}.
       !
       u0 = Psi_p
       workspace = Psi_h
       info = 0
       ! Compute u0 = Psi_p^{-1} * Source.
       call dgesv(nrows, nsource, workspace, nrows, IPIV, Psi_p, nrows, info)
       if (info > 0) then
          write(*, '(A)') "> Error in __solver_MOD_particular."
          write(*, '(A)') "> Eigenvectors of system AB are defective."
          write(*, '(A)') "> Program halted."
          stop
       endif
       ! Compute Psi_p = (1/mu0^2 - (wr + i * wi))^{-1} * (Psi_h^{-1} * Source).
       do mm = 1, nreal
          if (abs(diag(mm)) > cond) then
             Psi_p(mm, :) = Psi_p(mm, :) / diag(mm)
          endif
       enddo
       do mm = nreal+1, nrows, 2
          ! Psi_h stored in block-diagonal form for complex eigenvalues.
          x1 = Psi_p(mm, :)
          x2 = Psi_p(mm+1, :)
          tmp1 = wr(mm) - 1.0_dp / mu0 ** 2
          tmp2 = wi(mm)
          tmp0 = tmp1 ** 2 + tmp2 ** 2
          tmp1 = tmp1 / tmp0
          tmp2 = tmp2 / tmp0
          Psi_p(mm, :) = -tmp1 * x1 + tmp2 * x2
          Psi_p(mm+1, :) = -tmp2 * x1 - tmp1 * x2
       enddo
       ! Compute the error.
       u1 = Psi_p
       call dgemm('N', 'N', nrows, nsource, nrows, 1.0_dp, Psi_h, nrows, &
            u1, nrows, 0.0_dp, Psi_p, nrows)
       write (*, '(A)') "Warning in _solver__MOD_particular."
       write (*, '(A)') "Particular solution is a resonant frequency."
       write (*, '(A,E12.6)') &
            "Error in particular solution is  ||Psi_p''(0) - AB * Psi_p(0)|| = ", &
            norm2(Psi_p / mu0 ** 2 - matmul(AB,  Psi_p) - u0)
    endif
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute (M * Gamma_p) where
    !   M * Gamma_p = -mu0 * (B * M * Psi_p + S2).
    !
    call dgemm('N', 'N', nrows, nsource, nrows, -mu0, B, nrows, Psi_p, nrows, &
         -mu0, Gamma_p, nrows)
    !---------------------------------------------------------------------------
  end subroutine particular


  subroutine vdisord0_iq( nrows, nsource, mu0, A, B, AB, tol, &
       Psi_p, Gamma_p, Psi_h, Gamma_h, lambda, nreal, nzero )
    !!--------------------------------------------------------------------------
    !! Name: vdisord0_iq
    !!   Compute the general solution for the upper-right block of the zeroth
    !!   Fourier mode of the DRTE in the form of planewave solutions.
    !!
    !! SYNOPSIS
    !!   subroutine vdisord( N, nrows, nsource, mu, mu0, kappa, A, B, AB, tol,
    !!       Psi_p, Gamma_p, Psi_h, Gamma_h, lambda, nreal, nzero)
    !!
    !!   Integer ------------- N, nrows, nsource, nreal, nzero
    !!   Double precision ---- mu( N ), mu0, kappa, A( nrows, nrows ),
    !!                      B( nrows, nrows ), AB( nrows, nrows ), tol,
    !!                      Psi_p( nrows, nsource ), Gamma_h( nrows, nsource )
    !!                      Psi_h( nrows, nrows ), Gamma_h( nrows, nrows ),
    !!                      lambda( nrows )
    !!
    !! PURPOSE
    !!   Compute the general solution of the first order system of differential
    !!   equations
    !!
    !!       M *   Psi'(tau) = A * M * Gamma(tau) + S1 * exp(-tau/mu0),
    !!       M * Gamma'(tau) = B * M *   Psi(tau) + S2 * exp(-tau/mu0).
    !!
    !!    where A and B are linear combinations of the discretized radiative
    !!    transport operator given by
    !!
    !!       A = (-I + 0.5 * albedo (Zpos - D34 * Zneg) * W) M^{-1},
    !!       B = (-I + 0.5 * albedo (Zpos + D34 * Zneg) * W) M^{-1},
    !!
    !!    and S1 and S2 are the source terms
    !!
    !!       S1 = 0.5 * (S_pos - D34 * S_neg),
    !!       S2 = 0.5 * (S_pos + D34 * S_neg).
    !!
    !! ARGUMENTS
    !!   N (input) integer
    !!     The number of ordinate quadrature points.
    !!
    !!   nsource (input) integer
    !!     The number of source vectors.
    !!
    !!   mu (input) double precision, array(N)
    !!     Quadrature of ordinate directions.
    !!
    !!   mu0 (input) double precision
    !!     Ordinate direction of incident light.
    !!
    !!   A (input) double precision, array(nrows, nrows)
    !!     Discretization of the modified radiative transport operator
    !!      A = (-I + 0.5 * albedo (Zpos - D34 * Zneg) * W) M^{-1}.
    !!
    !!   B (input) double precision, array(nrows, nrows)
    !!     Discretization of the modified radiative transport operator
    !!       B = (-I + 0.5 * albedo (Zpos + D34 * Zneg) * W) M^{-1}.
    !!
    !!   AB (input) double precision, array(nrows, nrows)
    !!     Matrix system AB = A*B.
    !!
    !!   tol (input)
    !!     Tolerance for numerical precision.
    !!
    !!   Psi_p (input/output) double precision, array( nrows, nrows )
    !!     On entry, the source term S to the reduced system,
    !!        S = - S1 / mu0 + A * S2.
    !!     On exit, the particular solution Psi_p to the discrete system.
    !!
    !!   Gamma_p (input/output) double precision, array( nrows, nrows )
    !!     On entry, the source term S2.
    !!     On exit, the particular solution Gamma_p to the discrete system.
    !!
    !!   Psi_h (output) double precision, array ( nrows, nrows )
    !!     A matrix with row vectors corresponding to the planewave solutions.
    !!     M * Psi_h are the eigenvectors of the matrix (AB).
    !!
    !!   Gamma_h (output) double precision, array ( nrows, nrows )
    !!     A matrix with row vectors corresponding to the planewave solutions.
    !!     M * Gamma_h are the eigenvectors of the matrix (BA).
    !!
    !!   lambda (output) double precision, array ( nrows )
    !!     An array of eigenvalues/decay of the planewave solutions.
    !!
    !!   nreal (output) integer
    !!     The number of strictly real eigenvalues.
    !!
    !!   nzero (output) integer
    !!     The number of identically zero eigenvalues.
    !!---------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8), intent(in) :: nrows, nsource
    real(dp), intent(in) :: mu0
    real(dp), intent(in), dimension(:, :) :: A, B, AB ! dim(nrows, nrows)
    real(dp), intent(in) :: tol
    ! IN/OUT
    real(dp), intent(inout), dimension(:, :) :: Psi_p ! dim(nrows, nsource)
    real(dp), intent(inout), dimension(:, :) :: Gamma_p ! dim(nrows, nsource)
    ! OUTPUT
    real(dp), intent(out), dimension(:, :) :: Psi_h, Gamma_h ! dim(nrows, nrows)
    real(dp), intent(out), dimension(:) :: lambda ! dim(4N)
    integer(8), intent(out) :: nreal, nzero
    ! Workspace
    real(dp), dimension(nrows)  :: wr, wi
    real(dp), dimension(nrows, nrows) :: workspace ! workspace
    real(dp) :: cond = 10D-7
    !
    ! 1. Compute Psi_h and Gamma_h.
    !
    call eigensystem0_iq(nrows, B, AB, Psi_h, Gamma_h, wr, wi, lambda, &
         nreal, nzero, tol)
    !
    ! 2. Store Psi_p = -1/mu0 * S1 + A * S2.
    !
    call dgemm('N', 'N', nrows, nsource, nrows, 1.0_dp, A, nrows, &
         Gamma_p, nrows, -1.0_dp / mu0, Psi_p, nrows)
    !
    ! 3. Compute Psi_p and Gamma_p.
    !
    call particular(nrows, nreal, nsource, mu0, B, AB, wr, wi, Psi_h, &
         Psi_p, Gamma_p, workspace, cond)
  end subroutine vdisord0_iq


  subroutine eigensystem0_iq(nrows, B, AB, Psi_h, Gamma_h, wr, wi, &
       lambda, nreal, nzero, tol)
    !!--------------------------------------------------------------------------
    !! Name: eigensystem0_iq
    !!  Compute the eigenvalues and eigenvectors of the upper-right block of the
    !!  zeroth Fourier mode of the reduced DRTE.
    !!
    !! SYNOPSIS
    !!   subroutine eigensystem( N, nrows, mu, A, B, AB, Psi_h, Gamma_h, wr,
    !!       wi, lambda, nreal, nzero, tol )
    !!
    !!   Integer ------------- N, nrows, nreal, nzero
    !!   Double precision ----
    !!
    !!
    !! PURPOSE
    !!   Compute the eigenvalues and eigenvectors of the system
    !!
    !!      -lambda(j) *   Psi_h(j) = A * Gamma(j),
    !!      -lambda(j) * Gamma_h(j) = B *   Psi(j).
    !!
    !!   Note: If nzero = 1, an analytical solution is known for the two
    !!    eigenvectors in associated with the zero eigenvalue.
    !!
    !!       Psi[0,+](tau, mu) = i0,
    !!     Gamma[0,+](tau, mu) = 0,
    !!
    !!   and
    !!
    !!       Psi[0,-](tau, mu) = tau * i0,
    !!     Gamma[0,-](tau, mu) = mu / (1-g)
    !!
    !!   where i0 is an isotropic unpolarized beam, i0 = [1,0,0,0]^T, and g is
    !!   the anisotropy factor.
    !!
    !!
    !! PARAMETERS
    !!   N (input) integer
    !!     Number of quadrature points.
    !!
    !!   nrows (input) integer
    !!     Dimension of the matrices A, B, and AB.
    !!
    !!   mu (input) real, double precision, array( N )
    !!     Quadrature points of the ordinate directions.
    !!
    !!   B (input) real, double precision, array(nrows, nrows)
    !!     Discretization of the modified radiative transport operator
    !!       B = (-I + 0.5 * albedo (Zpos + D34 * Zneg) * W) M^{-1}.
    !!
    !!   AB (input) real, double precision, array(nrows, nrows)
    !!     Matrix system AB = A*B.
    !!
    !!   Psi_h (output) real, double precision, array ( nrows, nrows )
    !!     The rows of Psi_h are the eigenvectors of the matrix AB. The
    !!     eigenvectors are stored so their sorting matches the eigenvalues
    !!     wr and wi, and if the jth and (j+1)th eigenvalue correspond to a
    !!     conjugate pair, then the corresponding eigenvector is stored so
    !!            Eigenvector = Psi_h(j) +/- i * Psi_h(j+1).
    !!
    !!   Gamma_h (output) real, double precision, array ( nrows, nrows )
    !!     The rows of Gamma_h are the eigenvectors of the matrix BA stored
    !!     in the same order as Psi_h.
    !!
    !!   wr (output) real, double precision, array ( nrows )
    !!   wi (output) real, double precision, array ( nrows )
    !!     wr and wi contain the real and imaginary parts, respectively, of
    !!     the eigenvalues of AB. The eigenvalues are sorted so that
    !!
    !!   lambda (output) real, double precision, array ( nrows )
    !!     Eigenvalues of [[0, A], [B, 0]] corresponding to the positive real
    !!     branch. Eigenvalues are stored so that if there is a zero
    !!     eigenvalue it is stored first.
    !!
    !!   nreal (output) integer
    !!     Number of strictly real eigenvalues, nreal = nrows.
    !!
    !!   nzero (output) integer
    !!     Number of zero eigenvalues.
    !!
    !!   tol (input) real, double precision
    !!     Tolerance for determining if an eigenvalue is zero or non-zero.
    !!
    !! ---
    !! Written by Julia CLark.
    !! University of California, Merced.
    !! Contact: jclark@ucmerced.edu
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8), intent(in) :: nrows
    real(dp), intent(in), dimension(:, :) :: B, AB
    real(dp), intent(in) :: tol
    ! Workspace
    real(dp), intent(out), dimension(:, :) :: Psi_h, Gamma_h
    real(dp), intent(out), dimension(:)  :: wr, wi, lambda
    integer(8), intent(out) :: nzero, nreal
    ! Workspace: DGEEVX
    integer(8) :: info, lwork
    real(dp), dimension(1) :: qwork
    real(dp), allocatable, dimension(:) :: work
    integer(8) :: ihi, ilo
    integer(8), dimension(2*nrows-2) :: iwork
    real(dp), dimension(nrows) :: rconde, rcondv
    real(dp), dimension(nrows) :: scale
    real(dp) :: abnrm
    ! Preserve AB
    real(dp), dimension(nrows, nrows) :: workspace, Psi_adj
    ! Workspace
    integer(8) :: mm

    !---------------------------------------------------------------------------
    ! Compute the spectrum of the matrix AB.
    !
    workspace = AB
    ! Query work size.
    info = 0
    lwork = -1
    call dgeevx('B','V', 'V', 'B',&
         nrows, workspace, nrows, wr, wi, Psi_adj, nrows, Psi_h, nrows, &
         ilo, ihi, scale, abnrm, rconde, rcondv, &
         qwork, lwork, iwork, info)
    ! Solve for eigenvalues and eigenvectors.
    lwork = int(qwork(1), 8)
    allocate(work(lwork))
    call dgeevx('B','V', 'V', 'B',&
         nrows, workspace, nrows, wr, wi, Psi_adj, nrows, Psi_h, nrows, &
         ilo, ihi, scale, abnrm, rconde, rcondv, &
         work, lwork, iwork, info)
    if (info > 0) then
       write (*, '(A)') &
            "> Error in __solver_MOD_spectrum. "
       write (*, '(A)') &
            "> Lapack subroutine DGEEVX failed to compute all eigenvalues."
       write (*, '(A, I6)') "DGEEVX, info = ", info
       write (*, '(A)') "Program halted."
       stop
    endif
    deallocate(work)
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! This subroutine is not intended for use with complex eigenvalues. Check
    ! for complex eigenvalues and stop program if any found.
    !
    nreal = nrows
    if (maxval(wi) > tol) then
       write (*, '(A)') &
            "> Error in __solver_MOD_eigensystem0."
       write (*, '(A)') &
            "> Complex eigenvalues computed. Do not use  __sovler_MOD_vdisord0."
       write (*, '(A)') 'Program halted.'
       stop
    endif
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Find any zero eigenvalues and stop program in nzero > 1.  Sort so that
    ! the zero eigenvalue is stored first.
    !
    nzero = 0
    do mm = 1, nrows
       if (wr(mm) < tol) then
          ! If Re[lambda^2] <= 0, then Re[lambda] < 0.
          if (wr(mm) < -tol .and. abs(wi(mm)) < tol) then
             write (*, '(A)') &
                  "> Error in __solver_MOD_eigensystem."
             write (*, '(A)') &
                  "> Spectrum computed strictly imaginary eigenvalues, "
             write (*, '(A, E12.6, A, E12.6, A)') &
                  "lambda^2 = ", wr(mm), " + i * ", wi(mm), "."
             stop
          endif
          ! Insertion sort.
          ! Don't need to save new wr(1) and Psi_h(:, 1) since we know
          !    wr(1) = 0
          !    Psi_h(4*mm+1 : 4*mm+4, 1) = (1, 0, 0, 0)^T.
          nzero = nzero+1
          if (nzero > 1) then
             write (*, '(A)') &
                  "> Error in __solver_MOD_eigensystem."
             write (*, '(A, I2, A)') &
                  "> The zero eigenvalue has multiplicity ", nzero, &
                  ". The system is degenerate."
             write (*, '(A)') 'Program halted.'
             stop
          endif
          wr(mm) = wr(1)
          Psi_h(:, mm) = Psi_h(:, 1)
       endif
    enddo
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute lambda = sqrt(wr).
    !
    if (nzero == 1) then
       lambda(1) = 0.0_dp
    endif
    lambda(nzero+1 : nreal) = sqrt(wr(nzero+1 : nreal))
    !--------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute Gamma_h(j) where
    !
    !   Gamma_h(j) = - B * Psi_h(j) / lambda(j),
    !
    call dgemm('N', 'N', nrows, nrows, nrows, -1.0_dp, B, nrows, Psi_h, nrows, &
         0.0_dp, Gamma_h, nrows)
    do mm = nzero+1, nreal
       Gamma_h(:, mm) = Gamma_h(:, mm) / lambda(mm)
    enddo
    !---------------------------------------------------------------------------
  end subroutine eigensystem0_iq


end module solver
