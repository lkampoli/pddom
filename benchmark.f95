program benchmark
  !-----------------------------------------------------------------------------
  ! Vector Discrete Ordinate Method Benchmark
  !
  !  This program runs the vdisord method and tested it using published results
  !  from Mischenko, et.al.
  !
  !--
  ! Developed by J.P. Dark
  ! Contact: email@jpdark.com
  ! Last modified June 2017
  !-----------------------------------------------------------------------------


  implicit none

  ! Basic parameters
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: PI = 3.141592653589793_dp
  real(dp), parameter :: tol = 10E-12


  call vbirefl_oblique()


contains

  subroutine vbirefl_oblique()
    !!---------------------------------------------------------------------------
    !! Name: vbirefl_oblique
    !!  Compare vector discrete ordinate method to a polarized bidirectional
    !!  reflectance benchmark results for an oblique beam.
    !!
    !! Model Description:
    !!  * Venutian atmosphere (Hansen).
    !!  * Non-absorbing, albedo = 1.
    !!  * Oblique beam, mu0 = 0.5.
    !!  * Gridpoints: Shifted Gauss-Legendre quadrature with N=49.
    !!
    !! References
    !! -----------
    !!  Benchmark results:
    !!   Mishchenko, M. I., Dlugach, J. M., Chowdhary, J., & Zakharova, N. T.
    !!   Polarized bidirectional reflectance of optically thick sparse
    !!   particulate layers: An efficient numerically exact radiative-transfer
    !!   solution. JQRST. 2015.
    !!
    !!  Bidirection reflectance code available at:
    !!   https://www.giss.nasa.gov/staff/mmishchenko/brf/
    !!--------------------------------------------------------------------------
    use halfspace, only: driver_reflectance
    use quadratures, only: gausslegendre_shifted
    implicit none
    ! Array sizes
    integer(8), parameter :: N = 49   ! Number of gridpoints
    integer(8), parameter :: kmax = 3 ! Number of Fourier modes computed
    integer(8) :: L                   ! Order of phase matrix expansion
    ! Expansion coefficients
    real(dp), allocatable, dimension(:, :) :: FL_coeff ! dim(6, L+1)
    ! Model parameters
    real(dp), parameter :: albedo = 1.0_dp
    real(dp), parameter :: mu0 = 0.5_dp
    ! Quadrature
    real(dp), dimension(N) :: mu, wt
    ! Results
    real(dp), dimension(4*N, 4, kmax+1) :: reflect ! Reflectance
    real(dp), dimension(N, 16, kmax+1) :: vdisord
    ! File output name
    character(len=50) :: reflout
    ! Workspace
    integer(8) :: k, ii, jj ! Loop index.
    ! Vector Bidirection Reflectance Benchmark Results for comparison.
    integer, parameter :: N_bm  = 49         ! Number of gridpoints.
    real(dp), dimension(N_bm) :: mu_bm          ! Array of gridpoints.
    real(dp), dimension(N_bm, 16, 4) :: vbirefl    ! Benchmark results for comparison.

    write (*, '(A)') "> Run benchmark test: halfspace reflectance from oblique beam."

    !---------------------------------------------------------------------------
    ! Get coefficients for the Venutian atmsophere model.
    !
    include 'coeff_venutian.f95'
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute the quadrature.
    !  *Benchmark is evaluated on a N=49 point *shifted* Gaussian quadrature.
    !
    call gausslegendre_shifted(N, mu, wt)
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Run the model.
    !
    call driver_reflectance(N, L, kmax, tol, albedo, mu0, mu, wt, &
       FL_coeff, reflect)
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Print results and errors to file.
    !
    do k=0, kmax
       if (kmax < 10) then
          write (reflout, "(A25,I1,A4)") &
               'results/venus/venus_refl_',k, '.dat'
       else if (kmax < 100) then
          write (reflout, "(A25,I2,A4)") &
               'results/venus/venus_refl_',k, '.dat'
       else
          write (reflout, "(A25,I3,A4)") &
               'results/venus/venus_refl_',k, '.dat'
       endif
       call print_refl(N, mu, reflect(:, :, k+1), reflout)
    enddo
    !---------------------------------------------------------------------------
    write (*, '(A)') "> Run ./venus.sh to print results."

    !---------------------------------------------------------------------------
    ! Get Vector Bidirectional Reflectance benchmark results.
    include 'vbirefl_results.f95'
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! Change reflectance array output to match format of vbirefl results.
    !
    do k = 1, kmax+1
       do ii = 0, N - 1
          do jj = 1, 4
             vdisord(ii+1, 4*(jj-1)+1:4*(jj-1)+4, k) = reflect(4*ii+jj, :, k)
          enddo
       enddo
       write(*, '(A,I1,A,E12.4)') 'k = ', k-1,', Max Error = ', &
            maxval(abs(vdisord(:, : , k) - vbirefl(:, :, k)))
       write(*, '(A,I1,A,E12.4)') 'k = ', k-1,',  L2 Error = ', &
            norm2(vdisord(:, : , k) - vbirefl(:, :, k))
    enddo
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! Clean allocated arrays.
    !
    deallocate(FL_coeff)
    !---------------------------------------------------------------------------

  end subroutine vbirefl_oblique


  subroutine print_refl(N, mu, R0, filename)
    !!--------------------------------------------------------------------------
    !! Name: print_refl
    !!  Print reflectance matrix to file.
    !!
    !! Synopsis
    !!  print_refl(N, mu, R0, filename)
    !!
    !!   INPUT
    !!     Integer ------------- N
    !!     Double precision ---- mu(N), R0(4*N, 4)
    !!     Character------------ filename(*)
    !!
    !! Purpose
    !!  Print the reflectance matrix in the same format as results from the
    !!  vector bidirection reflectance benchmark results from Mishchenko, et.al.
    !!
    !! Arguments
    !!  N (input) integer
    !!   The number of ordinate quadrature points.
    !!
    !!  mu (input) double precision, array(N)
    !!   Quadrature of ordinate directions.
    !!
    !!  R0 (input) double precision, array(4N, 4)
    !!   Fourier mode of the reflectance matrix.
    !!
    !!  filename (input) character, array(*)
    !!   Name of the file the reflectance matrix is printed to.
    !!
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8), intent(in) :: N
    real(dp), intent(in), dimension(:) :: mu ! dim(N)
    real(dp), intent(in), dimension(:, :) :: R0 ! dim(4*N, 4)
    character(*), intent(in) :: filename
    ! Workspace
    integer(8) :: stat, ii, jj
    open (unit=10, file=filename, action="write", status="replace", iostat=stat)
    if (stat > 0) then
       write (*, '(A)') "> PROGRAM STOPPED BEFORE COMPLETION."
       write (*,'(A)') "> ERROR in subroutine print_refl."
       write (*,*) &
            "> Could not open file ", filename
       stop
    endif
    do ii = 0, N - 1
       write(10, '(17F14.6)', advance='no') mu(ii+1)
       do jj = 1, 4
          write(10, '(17F14.6)', advance='no') R0(4*ii+jj, :)
       enddo
       write(10, '(A)') ""
    enddo
    close(10)
  end subroutine print_refl


end program benchmark
