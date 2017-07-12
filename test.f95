program test
  use quadratures, only: clencurt, polytest
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: sp = kind(1.0)
  real(dp), parameter :: PI = 3.141592653589793_dp
  integer(8) :: N
  real(dp), dimension(:), allocatable :: mu, wt
  character(len=12) :: arg

  if (command_argument_count() .ne. 1) then
     write (*, '(A)') &
          '> Runtime error: must supply the number of gridpoints as input.'
     write (*, '(A)') 'Example: ./test 16'
     stop
  else
     call get_command_argument(1, arg)
     read(arg, *) N
  endif
  allocate(mu(N), wt(N))


  call clencurt(N, mu, wt)
  call polytest(N, 0, -1.0_dp, 1.0_dp, mu, wt)
  call polytest(N, 1, -1.0_dp, 1.0_dp, mu, wt)
  call polytest(N, 2, -1.0_dp, 1.0_dp, mu, wt)
  call polytest(N, 3, -1.0_dp, 1.0_dp, mu, wt)
  call polytest(N, 4, -1.0_dp, 1.0_dp, mu, wt)
  call polytest(N, 16, -1.0_dp, 1.0_dp, mu, wt)
  call polytest(N, 32, -1.0_dp, 1.0_dp, mu, wt)
  call polytest(N, 50, -1.0_dp, 1.0_dp, mu, wt)
  call polytest(N, 0, -1.0_dp, 1.0_dp, mu, wt)

  deallocate(mu, wt)

end program test
