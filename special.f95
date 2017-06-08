module special
  !-----------------------------------------------------------------------------
  ! Module: special
  !  This modules contains subroutines to compute special functions.
  !
  ! Subroutines:
  !  auxiliary (public)
  !    Compute the auxiliary rotation functions (also called generalized
  !    spherical functions) for a given k.
  !
  ! Useful Citations:
  !  Siewert, CE. "A discrete-ordinates solution for radiative-transfer
  !  models that include polarization effects." JQSRT, 2000.
  !
  !  Hovenier, JW, CVM van der Mee, and H Domke. Transfer of
  !  polarized light in planetary atmospheres: basic concepts and practical
  !  methods. Vol. 318. Springer Science & Business Media, 2014.
  !
  !  Varshalovich DA, Moskalev AN, Khersonskii VK. Quantum theory of angular
  !  momentum. World scientific, 1988.
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
  public :: auxiliary

contains

  subroutine auxiliary(k, N, L, mu, plk, rlk, tlk)
    !!--------------------------------------------------------------------------
    !! Name: auxiliary
    !!  Compute the auxiliary rotation functions (also called generalized
    !!  spherical functions) for a given k.
    !!
    !! Synopsis
    !!  auxiliary(k, N, L, mu, plk, rlk, tlk)
    !!
    !!     INPUT
    !!      Integer ------------- k, N, L
    !!      Double precision ---- mu(N)
    !!
    !!    OUTPUT
    !!      Double precision ---- plk(L+1, N), rlk(L+1, N), tlk(L+1, N)
    !!
    !!
    !! Purpose
    !!  Compute the auxiliary functions Plk, Rlk, and Tlk for a given k at
    !!  gridpoints mu (in the interval [-1,1]) using the initial values:
    !!
    !!    P[k,k](x) = sqrt((2k)!) * (1-x^2)^{k-2},
    !!    R[k,k](x) = 1/2^k * sqrt((2k)! / ((k+2)!(k-2)!)) * (1+x^2)
    !!                 * (1-x^2)^(1+l/2),
    !!    T[k,k](x) = 1/2^k * sqrt((2k)! / ((k+2)!(k-2)!)) * (1+x^2)
    !!                 * (1-x^2)^(1+l/2),
    !!
    !!
    !!  and the recurrence relations
    !!
    !!     sqrt((j+1)^2 - k^2) * P[j+1,k](x) = (2*j+1) * x * P[j,k](x)
    !!                                        - sqrt(j^2-k^2) P[j-1,k](x),
    !!
    !!     j * sqrt((j+1)^2 - k^2) * sqrt((j+1)^2 - 4) * R[j+1,k](x)
    !!        = j * (j+1) * (2*j+1) * x * R[j,k](x)
    !!          - (j+1) * sqrt(j^2 - 4) * sqrt(j^2 - k^2) * R[j-1,k](x)
    !!          - 2 * k * (2j+1) T[j,k](x),
    !!
    !!     j * sqrt((j+1)^2 - k^2) * sqrt((j+1)^2 - 4) * T[j+1,k](x)
    !!        = j * (j+1) * (2*j+1) * x * T[j,k](x)
    !!          - (j+1) * sqrt(j^2 - 4) * sqrt(j^2 - k^2) * T[j-1,k](x)
    !!          - 2 * k * (2j+1) R[j,k](x)
    !!
    !!
    !! Arguments
    !!  k (input) integer
    !!   Index for the auxiliary functions.
    !!
    !!  N (input) integer
    !!   The number of ordinate quadrature points.
    !!
    !!  L (input) integer
    !!   The order of the expansion.
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
    !!--------------------------------------------------------------------------
    implicit none
    ! INPUT
    integer(8), intent(in) :: k, N, L
    real(dp), dimension(:), intent(in) :: mu ! dim(N)
    ! OUTPUT
    real(dp), dimension(L+1, N), intent(out) :: plk, rlk, tlk
    ! Workspace
    real(dp), dimension(N) :: smu
    real(dp) :: a1, a2, c1, c2, d1
    real(dp) :: tmp
    integer(8) :: jj
    plk = 0.0_dp
    rlk = 0.0_dp
    tlk = 0.0_dp
    if (k > L .or. k < 0) then
       write (*,*) "> PROGRAM STOPPED BEFORE COMPLETING."
       write (*,*) "> ERROR in __special_MOD_auxiliary."
       write (*,*) "> Illegal value of k: must have 0 <= k <= L."
       stop
    endif
    smu = sqrt((1.0_dp - mu) * (1.0_dp + mu))
    if (k > 2) then
       !
       ! Compute initial value for Plk(x), Rlk(x), and Tlk(x).
       !
       tmp = 1.0_dp
       do jj = 1, k
          tmp = tmp * sqrt(0.25_dp + 0.25 * real(k, dp) / real(jj, dp))
       enddo
       plk(k+1, :) = tmp * smu **  k
       !
       ! Compute initial value for Plk(x), Rlk(x), and Tlk(x).
       !
       tmp = 0.5_dp ** k
       do jj = 1, k-2
          tmp = tmp * sqrt(1.0_dp + real(k + 2, dp) / real(jj, dp))
       enddo
       rlk(k+1, :) = tmp * (1.0_dp + mu ** 2) * smu ** (k - 2_dp)
       tlk(k+1, :) = 2.0_dp * mu * tmp * smu ** (k - 2_dp)
       !
       ! Compute higher order polynomials using recurrence relation.
       !
       if (L == k) then
          return
       endif
       plk(k+2, :) = sqrt(real(2 * k + 1, dp)) * mu * plk(k+1, :)
       tmp = sqrt(real((k+1) ** 2 - 4, dp))
       c1 = sqrt(real(2 *k + 1, dp)) * real((k+1), dp) / tmp
       d1 = 2.0_dp * sqrt(real(2 *k + 1, dp)) / tmp
       rlk(k+2, :) = c1 * mu * rlk(k+1, :) - d1 * tlk(k+1, :)
       tlk(k+2, :) = c1 * mu * tlk(k+1, :) - d1 * rlk(k+1, :)
       do jj = k+2, L
          tmp = sqrt(real(jj ** 2 - k ** 2, dp))
          a1 = real(2 * jj - 1, dp) / tmp
          a2 = sqrt(real((jj-1) ** 2 - k ** 2, dp)) / tmp
          plk(jj+1, :) = a1 * mu * plk(jj, :) - a2 * plk(jj-1, :)
          tmp = sqrt(real((jj ** 2 - 4) * (jj ** 2 - k ** 2), dp)) &
               * real(jj - 1, dp)
          c1 = real((2 * jj - 1) * jj * (jj - 1), dp) / tmp
          c2 = sqrt(real(((jj - 1) ** 2 - 4) * ((jj - 1) ** 2 - k ** 2), dp)) &
               * real(jj, dp) / tmp
          d1 = real(2 * k * (2 * jj - 1), dp) / tmp
          rlk(jj+1, :) = c1 * mu * rlk(jj, :) - c2 * rlk(jj-1, :) &
               - d1 * tlk(jj, :)
          tlk(jj+1, :) = c1 * mu * tlk(jj, :) - c2 * tlk(jj-1, :) &
               - d1 * rlk(jj, :)
       enddo

    else if (k == 2) then
       !
       ! Compute initial value for Plk(x).
       !
       plk(3, :) = sqrt(0.375_dp) * smu ** 2
       rlk(3, :) = 0.25_dp * (1.0_dp + mu ** 2)
       tlk(3, :) = 0.5_dp * mu
       if (L == 2) then
          return
       endif
       plk(4, :) = sqrt(5.0_dp) * mu * plk(3, :)
       rlk(4, :) = 3.0_dp * mu * rlk(3, :) - 2.0_dp * tlk(3, :)
       tlk(4, :) = 3.0_dp * mu * tlk(3, :) - 2.0_dp * rlk(3, :)
       !
       ! Compute higher order polynomials using recurrence relation.
       !
       do jj = 4, L
          ! Plk(x)
          tmp = sqrt(real(jj ** 2 - 4_dp, dp))
          a1 = real(2_dp * jj - 1_dp, dp) / tmp
          a2 = sqrt(real((jj-1) ** 2 - 4_dp, dp)) / tmp
          plk(jj+1, :) = a1 * mu * plk(jj, :) - a2 * plk(jj-1, :)
          ! Rlk(x), Tlk(x)
          tmp = real((jj ** 2 - 4) * (jj - 1), dp)
          c1 = real((2 * jj - 1) * jj * (jj - 1), dp) / tmp
          c2 = real(((jj - 1) ** 2 - 4) * jj, dp) / tmp
          d1 = real(4 * (2 * jj - 1), dp) / tmp
          rlk(jj+1, :) = c1 * mu * rlk(jj, :) - c2 * rlk(jj-1, :) &
               - d1 * tlk(jj, :)
          tlk(jj+1, :) = c1 * mu * tlk(jj, :) - c2 * tlk(jj-1, :) &
               - d1 * rlk(jj, :)
       enddo

    else if (k == 1) then
       !
       ! Compute initial value for Plk(x).
       !
       plk(2, :) = sqrt(0.5_dp) * smu ! Check sign.
       if (L == 1) then
          return
       endif
       plk(3, :) = sqrt(1.5_dp) * mu * smu
       tlk(3, :) = - 0.5_dp * smu
       rlk(3, :) = mu * tlk(3, :)
       if (L == 2) then
          return
       endif
       plk(4, :) = sqrt(3.125_dp) * mu * plk(3, :) &
            - sqrt(0.375_dp) * plk(2, :)
       rlk(4, :) = sqrt(5.625_dp) * mu * rlk(3, :) &
            - sqrt(0.625_dp) * tlk(3, :)
       tlk(4, :) = sqrt(5.625_dp) * mu * tlk(3, :) &
            - sqrt(0.625_dp) * rlk(3, :)
       !
       ! Compute higher order polynomials using recurrence relation.
       !
       do jj = 4, L
          !
          ! P[j+1,1](x) = a1(l) * x * P[j,1](x) - a2(l) * P[j-1,1](x)
          !
          tmp = sqrt(real(jj ** 2 - 1, dp))
          a1 = real(2 * jj - 1, dp) / tmp
          a2 = sqrt(real((jj - 1) ** 2 - 1, dp)) / tmp
          plk(jj+1, :) = a1 * mu * plk(jj, :) - a2 * plk(jj-1, :)
          !
          ! R[j+1,1](x) = c1 * mu R[j,1](x) - c2 * R[j-1,1](x)
          !               - d1 * T[j,1](x)
          !
          ! and
          !
          ! T[j+1,1](x) = c1 * mu T[j,1](x) - c2 * T[j-1,1](x)
          !               - d1 * R[j,1](x)
          !
          tmp = sqrt( real( (jj ** 2 - 4) * (jj ** 2 - 1), dp)) &
               * real(jj - 1, dp)
          c1 = real((2 * jj - 1) * jj * (jj - 1), dp) / tmp
          c2 = sqrt(real(((jj - 1) ** 2 - 4) * ((jj - 1) ** 2 - k), dp)) &
               * real(jj, dp) / tmp
          d1 = real(2 * k * (2 * jj - 1), dp) / tmp
          rlk(jj+1, :) = c1 * mu * rlk(jj, :) - c2 * rlk(jj-1, :) &
               - d1 * tlk(jj, :)
          tlk(jj+1, :) = c1 * mu * tlk(jj, :) - c2 * tlk(jj-1, :) &
               - d1 * rlk(jj, :)
       enddo

    else if (k == 0) then
       !
       ! Compute initial values for Plk(x) and Rlk(x). Tlk(x) = 0.
       !
       plk(1, :) = 1.0_dp
       if (L == 0) then
          return
       endif
       plk(2, :) = mu
       if (L == 1) then
          return
       endif
       plk(3, :) = 1.5_dp * mu ** 2 - 0.5_dp
       rlk(3, :) = sqrt(0.375_dp) * smu ** 2
       if (L == 2) then
          return
       endif
       plk(4, :) = 0.5 * (5.0_dp * mu ** 3 - 3.0_dp * mu)
       rlk(4, :) = sqrt(5.0_dp) * mu * rlk(3, :)
       !
       ! Compute higher order values.
       !
       do jj = 4, L
          !
          ! (j+1) * P[j+1,k](x) = (2j+1) * x * P[j,k](x) - j * P[j-1,k](x).
          !
          tmp = real(jj, dp)
          a1 = real(2_dp * jj - 1_dp, dp) / tmp
          a2 = real(jj - 1, dp) / tmp
          plk(jj+1, :) = a1 * mu * plk(jj, :) - a2 * plk(jj-1, :)
          !
          ! sqrt((j+1)^2 - 4) * R[j+1,k](x) = (2j+1) * R[j,k}(x)
          !     - sqrt(j^2 - 4) * R[j-1,k](x).
          !
          tmp = sqrt(real(jj ** 2 - 4_dp, dp))
          c1 = real(2_dp * jj - 1_dp, dp) / tmp
          c2 = sqrt(real((jj-1) ** 2 - 4_dp, dp)) / tmp
          rlk(jj+1, :) = c1 * mu * rlk(jj, :) - c2 * rlk(jj-1, :)
       enddo
    endif
    !---------------------------------------------------------------------------
    ! Correct for sign error in convention chosen?
    !
    rlk = -rlk
    tlk = -tlk
    !---------------------------------------------------------------------------
  end subroutine auxiliary

end module special
