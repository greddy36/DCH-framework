      subroutine vamp_integrate(fct, xl, xu, d,                         &
     &numCallsGridOpt, numCallsIntEval, integral, integralErr)
    
      use kinds
      use exceptions
      use tao_random_numbers
      use vamp
      implicit none

      interface
          function fct (x, dim, params) result (f)
              real(kind=8), dimension(:), intent(in) :: x
              integer(kind=4), intent(in) :: dim
              real(kind=8), dimension(:), intent(in) :: params
              real(kind=8) :: f
          end function fct
      end interface
      real(kind=8), dimension(d), intent(in) :: xl
      real(kind=8), dimension(d), intent(in) :: xu
      integer(kind=4), intent(in) :: d
      integer(kind=4), intent(in) :: numCallsGridOpt
      integer(kind=4), intent(in) :: numCallsIntEval
      real(kind=8), intent(out) :: integral
      real(kind=8), intent(out) :: integralErr

      integer :: i

      type(exception) :: exc
      type(tao_random_state) :: rng
      type(vamp_grid) :: grid
      real(kind=default), dimension(2,d) :: domain
      real(kind=default) :: vamp_integral, vamp_error, vamp_chi2

!      print *, "<vamp_integrate>:"
!      print *, " numCallsGridOpt = ", numCallsGridOpt
!      print *, " numCallsIntEval = ", numCallsIntEval

      do i = 1, d
          domain(1,i) = xl(i)
          domain(2,i) = xu(i)
      end do

      call tao_random_create(rng, seed=0)
      call clear_exception(exc)
      call vamp_create_grid(grid, domain, num_calls=numCallsGridOpt,    &
     &exc=exc)
      call handle_exception(exc)

      call clear_exception(exc)
      call vamp_sample_grid(rng, grid, vamp_fct, NO_DATA, 6, exc=exc)
      call handle_exception(exc)

      call clear_exception(exc)
      call vamp_discard_integral(grid, num_calls=numCallsIntEval,       &
     &exc=exc)
      call handle_exception(exc)
      call clear_exception(exc)
      call vamp_sample_grid(rng, grid, vamp_fct, NO_DATA, 4,            &
     &vamp_integral, vamp_error, vamp_chi2, exc=exc)
      call handle_exception(exc)
      call clear_exception(exc)
      call vamp_delete_grid(grid)
      call handle_exception(exc)
!      print *, "integral = ", vamp_integral, " +/- ", vamp_error
!      print *, " (chi^2 = ", vamp_chi2, ")"

      integral = vamp_integral
      integralErr = vamp_error

      contains
          function vamp_fct (x, data, weights, channel, grids) result   &
     &(f)
              real(kind=default), dimension(:), intent(in) :: x
              class(vamp_data_t), intent(in) :: data
              real(kind=default), dimension(:), intent(in), optional :: &
     &weights
              integer, intent(in), optional :: channel
              type(vamp_grid), dimension(:), intent(in), optional ::    &
     &grids
              real(kind=default) :: f
              real(kind=8), dimension(0) :: dummy_params
!             print *, "<vamp_fct>:"
!             print *, " x = ", x
              f = fct(x, d, dummy_params) 
          end function vamp_fct   
     
      end subroutine vamp_integrate
