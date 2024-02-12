! vampi.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vamp_serial_mpi
  use vamp, &
       vamp0_create_grid => vamp_create_grid, &
       vamp0_discard_integral => vamp_discard_integral, &
       vamp0_reshape_grid => vamp_reshape_grid, &
       vamp0_sample_grid => vamp_sample_grid, &
       vamp0_delete_grid => vamp_delete_grid, &
       vamp0_print_history => vamp_print_history, &
       vamp0_grids => vamp_grids, &
       vamp0_create_grids => vamp_create_grids, &
       vamp0_discard_integrals => vamp_discard_integrals, &
       vamp0_update_weights => vamp_update_weights, &
       vamp0_refine_weights => vamp_refine_weights, &
       vamp0_delete_grids => vamp_delete_grids, &
       vamp0_sample_grids => vamp_sample_grids, &
       vamp0_warmup_grid => vamp_warmup_grid, &
       vamp0_warmup_grids => vamp_warmup_grids, &
       vamp0_next_event => vamp_next_event, &
       vamp0_write_grid => vamp_write_grid, &
       vamp0_read_grid => vamp_read_grid, &
       vamp0_write_grids => vamp_write_grids, &
       vamp0_read_grids => vamp_read_grids, &
       VAMP0_RCS_ID => VAMP_RCS_ID
  public
end module vamp_serial_mpi
module vamp_parallel_mpi
  use kinds
  use utils
  use tao_random_numbers
  use exceptions
  use mpi90
  use divisions
  use vamp_serial_mpi !NODEP!
  use iso_fortran_env
  implicit none
  private
  public :: vamp_create_grid
  public :: vamp_discard_integral
  public :: vamp_reshape_grid
  public :: vamp_sample_grid
  public :: vamp_delete_grid
  public :: vamp_print_history
  private :: vamp_print_one_history, vamp_print_histories
  public :: vamp_create_grids
  public :: vamp_discard_integrals
  public :: vamp_update_weights
  public :: vamp_refine_weights
  public :: vamp_delete_grids
  public :: vamp_sample_grids
  private :: object
  private :: schedule
  public :: vamp_warmup_grid
  public :: vamp_warmup_grids
  public :: vamp_next_event
  private :: vamp_next_event_single, vamp_next_event_multi
  public :: vamp_write_grid, vamp_read_grid
  private :: write_grid_unit, write_grid_name
  private :: read_grid_unit, read_grid_name
  public :: vamp_write_grids, vamp_read_grids
  private :: write_grids_unit, write_grids_name
  private :: read_grids_unit, read_grids_name
    public :: vamp_send_grid
    public :: vamp_receive_grid
    public :: vamp_broadcast_grid
    public :: vamp_broadcast_grids
    public :: vamp_send_history
    public :: vamp_receive_history
  interface vamp_print_history
     module procedure vamp_print_one_history, vamp_print_histories
  end interface
  interface vamp_next_event
     module procedure vamp_next_event_single, vamp_next_event_multi
  end interface
  interface vamp_write_grid
     module procedure write_grid_unit, write_grid_name
  end interface
  interface vamp_read_grid
     module procedure read_grid_unit, read_grid_name
  end interface
    interface vamp_write_grids
       module procedure write_grids_unit, write_grids_name
    end interface
    interface vamp_read_grids
       module procedure read_grids_unit, read_grids_name
    end interface

    interface vamp_broadcast_grid
       module procedure &
            vamp_broadcast_one_grid, vamp_broadcast_many_grids
    end interface
  integer, public, parameter :: VAMP_ROOT = 0
  real(kind=default), private, parameter :: VAMP_MAX_WASTE = 1.0
  ! real(kind=default), private, parameter :: VAMP_MAX_WASTE = 0.3
  integer, public, parameter :: &
       TAG_INTEGRAL = 1, &
       TAG_STD_DEV = 2, &
       TAG_GRID = 3, &
       TAG_HISTORY = 6, &
       TAG_NEXT_FREE = 9
  type, public :: vamp_grids
     !!! private
     type(vamp0_grids) :: g0
     logical, dimension(:), pointer :: active
     integer, dimension(:), pointer :: proc
     real(kind=default), dimension(:), pointer :: integrals, std_devs
  end type vamp_grids
  character(len=*), public, parameter :: VAMPI_RCS_ID = &
       "$Id: vampi.nw 314 2010-04-17 20:32:33Z ohl $"
contains
  subroutine vamp_create_grid &
       (g, domain, num_calls, num_div, &
        stratified, quadrupole, covariance, map, exc)
    type(vamp_grid), intent(inout) :: g
    real(kind=default), dimension(:,:), intent(in) :: domain
    integer, intent(in) :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole, covariance
    real(kind=default), dimension(:,:), intent(in), optional :: map
    type(exception), intent(inout), optional :: exc
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_create_grid &
            (g, domain, num_calls, num_div, &
             stratified, quadrupole, covariance, map, exc)
    else
       call vamp_create_empty_grid (g)
    end if
  end subroutine vamp_create_grid
  subroutine vamp_discard_integral &
       (g, num_calls, num_div, stratified, quadrupole, covariance, exc)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole, covariance
    type(exception), intent(inout), optional :: exc
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_discard_integral &
            (g, num_calls, num_div, stratified, quadrupole, covariance, exc)
    end if
  end subroutine vamp_discard_integral
  subroutine vamp_reshape_grid &
       (g, num_calls, num_div, stratified, quadrupole, covariance, exc)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole, covariance
    type(exception), intent(inout), optional :: exc
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_reshape_grid &
            (g, num_calls, num_div, stratified, quadrupole, covariance, exc)
    end if
  end subroutine vamp_reshape_grid
  subroutine vamp_sample_grid &
       (rng, g, func, iterations, integral, std_dev, avg_chi2, accuracy, &
        channel, weights, grids, exc, history)
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grid), intent(inout) :: g
    integer, intent(in) :: iterations
    real(kind=default), intent(out), optional :: integral, std_dev, avg_chi2
    real(kind=default), intent(in), optional :: accuracy
    integer, intent(in), optional :: channel
    real(kind=default), dimension(:), intent(in), optional :: weights
    type(vamp_grid), dimension(:), intent(inout), optional :: grids
    type(exception), intent(inout), optional :: exc
    type(vamp_history), dimension(:), intent(inout), optional :: history
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    character(len=*), parameter :: FN = "vamp_sample_grid"
    real(kind=default) :: local_integral, local_std_dev, local_avg_chi2
    type(vamp_grid), dimension(:), allocatable :: gs, gx
    integer, dimension(:,:), pointer :: d
    integer :: iteration, i
    integer :: num_proc, proc_id, num_workers
    nullify (d)
    call mpi90_size (num_proc)
    call mpi90_rank (proc_id)
    iterate: do iteration = 1, iterations
       if (proc_id == VAMP_ROOT) then
          call vamp_distribute_work (num_proc, vamp_rigid_divisions (g), d)
          num_workers = max (1, product (d(2,:)))
       end if
       call mpi90_broadcast (num_workers, VAMP_ROOT)
       if ((present (grids)) .and. (num_workers > 1)) then
          call vamp_broadcast_grid (grids, VAMP_ROOT)
       end if
       if (proc_id == VAMP_ROOT) then
          allocate (gs(num_workers), gx(vamp_fork_grid_joints (d)))
          call vamp_create_empty_grid (gs)
          call vamp_fork_grid (g, gs, gx, d, exc)
          do i = 2, num_workers
             call vamp_send_grid (gs(i), i-1, 0)
          end do
       else if (proc_id < num_workers) then
          call vamp_receive_grid (g, VAMP_ROOT, 0)
       end if
       if (proc_id == VAMP_ROOT) then
          if (num_workers > 1) then
             call vamp_sample_grid0 &
                  (rng, gs(1), func, channel, weights, grids, exc)
          else
             call vamp_sample_grid0 &
                  (rng, g, func, channel, weights, grids, exc)
          end if
       else if (proc_id < num_workers) then
          call vamp_sample_grid0 &
               (rng, g, func, channel, weights, grids, exc)
       end if
       if (proc_id == VAMP_ROOT) then
          do i = 2, num_workers
             call vamp_receive_grid (gs(i), i-1, 0)
          end do
          call vamp_join_grid (g, gs, gx, d, exc)
          call vamp0_delete_grid (gs)
          deallocate (gs, gx)
          call vamp_refine_grid (g)
          call vamp_average_iterations &
               (g, iteration, local_integral, local_std_dev, local_avg_chi2)
          if (present (history)) then
             if (iteration <= size (history)) then
                call vamp_get_history &
                     (history(iteration), g, &
                      local_integral, local_std_dev, local_avg_chi2)
             else
                call raise_exception (exc, EXC_WARN, FN, "history too short")
             end if
             call vamp_terminate_history (history(iteration+1:))
          end if
          if (present (accuracy)) then
             if (local_std_dev <= accuracy * local_integral) then
                call raise_exception (exc, EXC_INFO, FN, &
                     "requested accuracy reached")
                exit iterate
             end if
          end if
       else if (proc_id < num_workers) then
          call vamp_send_grid (g, VAMP_ROOT, 0)
       end if
    end do iterate
    if (proc_id == VAMP_ROOT) then
       deallocate (d)
       if (present (integral)) then
          integral = local_integral
       end if
       if (present (std_dev)) then
          std_dev = local_std_dev
       end if
       if (present (avg_chi2)) then
          avg_chi2 = local_avg_chi2
       end if
    end if
  end subroutine vamp_sample_grid
  subroutine vamp_delete_grid (g)
    type(vamp_grid), intent(inout) :: g
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_reshape_grid (g)
    end if
  end subroutine vamp_delete_grid
  subroutine vamp_print_one_history (h, tag)
    type(vamp_history), dimension(:), intent(in) :: h
    character(len=*), intent(in), optional :: tag
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_print_history (h, tag)
    end if
  end subroutine vamp_print_one_history
  subroutine vamp_print_histories (h, tag)
    type(vamp_history), dimension(:,:), intent(in) :: h
    character(len=*), intent(in), optional :: tag
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_print_history (h, tag)
    end if
  end subroutine vamp_print_histories
  subroutine vamp_create_grids (g, domain, num_calls, weights, maps, &
                                num_div, stratified, quadrupole, exc)
    type(vamp_grids), intent(inout) :: g
    real(kind=default), dimension(:,:), intent(in) :: domain
    integer, intent(in) :: num_calls
    real(kind=default), dimension(:), intent(in) :: weights
    real(kind=default), dimension(:,:,:), intent(in), optional :: maps
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    type(exception), intent(inout), optional :: exc
    integer :: proc_id, nch
    call mpi90_rank (proc_id)
    nch = size (weights)
    allocate (g%active(nch), g%proc(nch), g%integrals(nch), g%std_devs(nch))
    if (proc_id == VAMP_ROOT) then
       call vamp0_create_grids (g%g0, domain, num_calls, weights, maps, &
                                num_div, stratified, quadrupole, exc)
    else
       allocate (g%g0%grids(nch), g%g0%weights(nch), g%g0%num_calls(nch))
       call vamp_create_empty_grid (g%g0%grids)
    end if
  end subroutine vamp_create_grids
  subroutine vamp_discard_integrals &
       (g, num_calls, num_div, stratified, quadrupole, exc)
    type(vamp_grids), intent(inout) :: g
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    type(exception), intent(inout), optional :: exc
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_discard_integrals &
            (g%g0, num_calls, num_div, stratified, quadrupole, exc)
    end if
  end subroutine vamp_discard_integrals
  subroutine vamp_update_weights &
       (g, weights, num_calls, num_div, stratified, quadrupole, exc)
    type(vamp_grids), intent(inout) :: g
    real(kind=default), dimension(:), intent(in) :: weights
    integer, intent(in), optional :: num_calls
    integer, dimension(:), intent(in), optional :: num_div
    logical, intent(in), optional :: stratified, quadrupole
    type(exception), intent(inout), optional :: exc
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_update_weights &
            (g%g0, weights, num_calls, num_div, stratified, quadrupole, exc)
    end if
  end subroutine vamp_update_weights
  subroutine vamp_refine_weights (g, power)
    type(vamp_grids), intent(inout) :: g
    real(kind=default), intent(in), optional :: power
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_refine_weights (g%g0, power)
    end if
  end subroutine vamp_refine_weights
  subroutine vamp_delete_grids (g)
    type(vamp_grids), intent(inout) :: g
    character(len=*), parameter :: FN = "vamp_delete_grids"
    deallocate (g%active, g%proc, g%integrals, g%std_devs)
    call vamp0_delete_grids (g%g0)
  end subroutine vamp_delete_grids
  subroutine vamp_sample_grids &
       (rng, g, func, iterations, integral, std_dev, avg_chi2, &
        accuracy, history, histories, exc)
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grids), intent(inout) :: g
    integer, intent(in) :: iterations
    real(kind=default), intent(out), optional :: integral, std_dev, avg_chi2
    real(kind=default), intent(in), optional :: accuracy
    type(vamp_history), dimension(:), intent(inout), optional :: history
    type(vamp_history), dimension(:,:), intent(inout), optional :: histories
    type(exception), intent(inout), optional :: exc
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    character(len=*), parameter :: FN = "vamp_sample_grids"
    integer :: num_proc, proc_id, nch, ch, iteration
    real(kind=default), dimension(size(g%g0%weights)) :: weights
    real(kind=default) :: local_integral, local_std_dev, local_avg_chi2
    real(kind=default) :: current_accuracy, waste
    logical :: distribute_complete_grids
    call mpi90_size (num_proc)
    call mpi90_rank (proc_id)
    nch = size (g%g0%weights)
    if (proc_id == VAMP_ROOT) then
       g%active = (g%g0%num_calls >= 2)
       where (g%active)
          weights = g%g0%num_calls
       elsewhere
          weights = 0.0
       endwhere
       weights = weights / sum (weights)
       call schedule (weights, num_proc, g%proc, waste)
       distribute_complete_grids = (waste <= VAMP_MAX_WASTE)
    end if
    call mpi90_broadcast (weights, VAMP_ROOT)
    call mpi90_broadcast (g%active, VAMP_ROOT)
    call mpi90_broadcast (distribute_complete_grids, VAMP_ROOT)
    if (distribute_complete_grids) then
       call mpi90_broadcast (g%proc, VAMP_ROOT)
    end if
    iterate: do iteration = 1, iterations
       if (distribute_complete_grids) then
          call vamp_broadcast_grid (g%g0%grids, VAMP_ROOT)
          do ch = 1, nch
             if (g%active(ch)) then
                if (proc_id == g%proc(ch)) then
                   call vamp0_discard_integral (g%g0%grids(ch))
                   call vamp_sample_grid0 &
                        (rng, g%g0%grids(ch), func, ch, weights, g%g0%grids, exc)
                   call vamp_average_iterations &
                        (g%g0%grids(ch), iteration, g%integrals(ch), g%std_devs(ch), local_avg_chi2)
                   if (present (histories)) then
                      if (iteration <= ubound (histories, dim=1)) then
                         call vamp_get_history &
                              (histories(iteration,ch), g%g0%grids(ch), &
                               g%integrals(ch), g%std_devs(ch), local_avg_chi2)
                      else
                         call raise_exception (exc, EXC_WARN, FN, "history too short")
                      end if
                      call vamp`'_terminate_history (histories(iteration+1:,ch))
                   end if
                end if
             else
                call vamp_nullify_variance (g%g0%grids(ch))
                call vamp_nullify_covariance (g%g0%grids(ch))
             end if
          end do
          do ch = 1, nch
             if (g%active(ch) .and. (proc_id == g%proc(ch))) then
                call vamp_refine_grid (g%g0%grids(ch))
                if (proc_id /= VAMP_ROOT) then
                   call mpi90_send (g%integrals(ch), VAMP_ROOT, object (ch, TAG_INTEGRAL))
                   call mpi90_send (g%std_devs(ch), VAMP_ROOT, object (ch, TAG_STD_DEV))
                   call vamp_send_grid (g%g0%grids(ch), VAMP_ROOT, object (ch, TAG_GRID))
                   if (present (histories)) then
                      call vamp_send_history &
                           (histories(iteration,ch), VAMP_ROOT, object (ch, TAG_HISTORY))
                   end if
                end if
             end if
          end do
          if (proc_id == VAMP_ROOT) then
             do ch = 1, nch
                if (g%active(ch) .and. (g%proc(ch) /= proc_id)) then
                   call mpi90_receive (g%integrals(ch), g%proc(ch), object (ch, TAG_INTEGRAL))
                   call mpi90_receive (g%std_devs(ch), g%proc(ch), object (ch, TAG_STD_DEV))
                   call vamp_receive_grid (g%g0%grids(ch), g%proc(ch), object (ch, TAG_GRID))
                   if (present (histories)) then
                      call vamp_receive_history &
                           (histories(iteration,ch), g%proc(ch), object (ch, TAG_HISTORY))
                   end if
                end if
             end do
             call vamp_reduce_channels (g%g0, g%integrals, g%std_devs, g%active)
             call vamp_average_iterations &
                  (g%g0, iteration, local_integral, local_std_dev, local_avg_chi2)
             if (present (history)) then
                if (iteration <= size (history)) then
                   call vamp_get_history &
                        (history(iteration), g%g0, local_integral, local_std_dev, &
                         local_avg_chi2)
                else
                   call raise_exception (exc, EXC_WARN, FN, "history too short")
                end if
                call vamp_terminate_history (history(iteration+1:))
             end if
          end if
       else
          do ch = 1, size (g%g0%grids)
             if (g%active(ch)) then
                call vamp_discard_integral (g%g0%grids(ch))
                if (present (histories)) then
                   call vamp_sample_grid &
                        (rng, g%g0%grids(ch), func, 1, g%integrals(ch), g%std_devs(ch), &
                         channel = ch, weights = weights, grids = g%g0%grids, &
                         history = histories(iteration:iteration,ch))
                else      
                   call vamp_sample_grid &
                        (rng, g%g0%grids(ch), func, 1, g%integrals(ch), g%std_devs(ch), &
                         channel = ch, weights = weights, grids = g%g0%grids)
                end if
             else
                if (proc_id == VAMP_ROOT) then
                   call vamp_nullify_variance (g%g0%grids(ch))
                   call vamp_nullify_covariance (g%g0%grids(ch))
                end if
             end if
          end do
          if (proc_id == VAMP_ROOT) then
             call vamp_reduce_channels (g%g0, g%integrals, g%std_devs, g%active)
             call vamp_average_iterations &
                  (g%g0, iteration, local_integral, local_std_dev, local_avg_chi2)
             if (present (history)) then
                if (iteration <= size (history)) then
                   call vamp_get_history &
                        (history(iteration), g%g0, local_integral, local_std_dev, &
                         local_avg_chi2)
                else
                   call raise_exception (exc, EXC_WARN, FN, "history too short")
                end if
                call vamp`'_terminate_history (history(iteration+1:))
             end if
          end if
       end if
       if (present (accuracy)) then
          if (proc_id == VAMP_ROOT) then
             current_accuracy = local_std_dev / local_integral
          end if
          call mpi90_broadcast (current_accuracy, VAMP_ROOT)
          if (current_accuracy <= accuracy) then
             call raise_exception (exc, EXC_INFO, FN, &
                  "requested accuracy reached")
             exit iterate
          end if
       end if
    end do iterate
    if (present (integral)) then
       call mpi90_broadcast (local_integral, VAMP_ROOT)
       integral = local_integral
    end if
    if (present (std_dev)) then
       call mpi90_broadcast (local_std_dev, VAMP_ROOT)
       std_dev = local_std_dev
    end if
    if (present (avg_chi2)) then
       call mpi90_broadcast (local_avg_chi2, VAMP_ROOT)
       avg_chi2 = local_avg_chi2
    end if
  end subroutine vamp_sample_grids
  pure function object (ch, obj) result (tag)
    integer, intent(in) :: ch, obj
    integer :: tag
    tag = 100 * ch + obj
  end function object
  pure subroutine schedule (jobs, num_procs, assign, waste)
    real(kind=default), dimension(:), intent(in) :: jobs
    integer, intent(in) :: num_procs
    integer, dimension(:), intent(out) :: assign
    real(kind=default), intent(out), optional :: waste
    integer, dimension(size(jobs)) :: idx
    real(kind=default), dimension(size(jobs)) :: sjobs
    real(kind=default), dimension(num_procs) :: fill
    integer :: job, proc
    sjobs = jobs / sum (jobs) * num_procs
    idx = (/ (job, job = 1, size(jobs)) /)
    call sort (sjobs, idx, reverse = .true.)
    fill = 0.0
    fill(VAMP_ROOT+1) = 0.1
    do job = 1, size (sjobs)
       proc = sum (minloc (fill))
       fill(proc) = fill(proc) + sjobs(job)
       assign(idx(job)) = proc - 1
    end do
    if (present (waste)) then
       waste = 1.0 - sum (fill) / (num_procs * maxval (fill))
    end if
  end subroutine schedule
  subroutine vamp_next_event_single &
       (x, rng, g, func, weight, channel, weights, grids, exc)
    real(kind=default), dimension(:), intent(out) :: x
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grid), intent(inout) :: g
    real(kind=default), intent(out), optional :: weight
    integer, intent(in), optional :: channel
    real(kind=default), dimension(:), intent(in), optional :: weights
    type(vamp_grid), dimension(:), intent(in), optional :: grids
    type(exception), intent(inout), optional :: exc
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_next_event &
            (x, rng, g, func, weight, channel, weights, grids, exc)
    end if
  end subroutine vamp_next_event_single
  subroutine vamp_next_event_multi (x, rng, g, func, phi, weight, exc)
    real(kind=default), dimension(:), intent(out) :: x
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grids), intent(inout) :: g
    real(kind=default), intent(out), optional :: weight
    type(exception), intent(inout), optional :: exc
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    interface
       pure function phi (xi, channel) result (x)
         use kinds
         real(kind=default), dimension(:), intent(in) :: xi
         integer, intent(in) :: channel
         real(kind=default), dimension(size(xi)) :: x
       end function phi
    end interface
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_next_event (x, rng, g%g0, func, phi, weight, exc)
    end if
  end subroutine vamp_next_event_multi
  subroutine vamp_warmup_grid (rng, g, func, iterations, exc, history)
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grid), intent(inout) :: g
    integer, intent(in) :: iterations
    type(exception), intent(inout), optional :: exc
    type(vamp_history), dimension(:), intent(inout), optional :: history
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    call vamp_sample_grid &
       (rng, g, func, iterations - 1, exc = exc, history = history)
    call vamp_sample_grid0 (rng, g, func, exc = exc)
  end subroutine vamp_warmup_grid
  subroutine vamp_warmup_grids &
       (rng, g, func, iterations, history, histories, exc)
    type(tao_random_state), intent(inout) :: rng
    type(vamp_grids), intent(inout) :: g
    integer, intent(in) :: iterations
    type(vamp_history), dimension(:), intent(inout), optional :: history
    type(vamp_history), dimension(:,:), intent(inout), optional :: histories
    type(exception), intent(inout), optional :: exc
    interface
       function func (xi, data, weights, channel, grids) result (f)
         use kinds
         use vamp_grid_type !NODEP!
         import vamp_data_t
         real(kind=default), dimension(:), intent(in) :: xi
         class(vamp_data_t), intent(in) :: data
         real(kind=default), dimension(:), intent(in), optional :: weights
         integer, intent(in), optional :: channel
         type(vamp_grid), dimension(:), intent(in), optional :: grids
         real(kind=default) :: f
       end function func
    end interface
    integer :: ch
    call vamp0_sample_grids (rng, g%g0, func, iterations - 1, exc = exc, &
                             history = history, histories = histories)
    do ch = 1, size (g%g0%grids)
       ! if (g%g0%grids(ch)%num_calls >= 2) then
          call vamp_sample_grid0 (rng, g%g0%grids(ch), func, exc = exc)
       ! end if
    end do
  end subroutine vamp_warmup_grids
  subroutine write_grid_unit (g, unit)
    type(vamp_grid), intent(in) :: g
    integer, intent(in) :: unit
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_write_grid (g, unit)
    end if
  end subroutine write_grid_unit
  subroutine read_grid_unit (g, unit)
    type(vamp_grid), intent(inout) :: g
    integer, intent(in) :: unit
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_read_grid (g, unit)
    end if
  end subroutine read_grid_unit
  subroutine write_grid_name (g, name)
    type(vamp_grid), intent(inout) :: g
    character(len=*), intent(in) :: name
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_write_grid (g, name)
    end if
  end subroutine write_grid_name
  subroutine read_grid_name (g, name)
    type(vamp_grid), intent(inout) :: g
    character(len=*), intent(in) :: name
    integer :: proc_id
    call mpi90_rank (proc_id)
    if (proc_id == VAMP_ROOT) then
       call vamp0_read_grid (g, name)
    end if
  end subroutine read_grid_name
    subroutine write_grids_unit (g, unit)
      type(vamp_grids), intent(in) :: g
      integer, intent(in) :: unit
      integer :: proc_id
      call mpi90_rank (proc_id)
      if (proc_id == VAMP_ROOT) then
         call vamp0_write_grids (g%g0, unit)
      end if
    end subroutine write_grids_unit

    subroutine read_grids_unit (g, unit)
      type(vamp_grids), intent(inout) :: g
      integer, intent(in) :: unit
      integer :: proc_id
      call mpi90_rank (proc_id)
      if (proc_id == VAMP_ROOT) then
         call vamp0_read_grids (g%g0, unit)
      end if
    end subroutine read_grids_unit

    subroutine write_grids_name (g, name)
      type(vamp_grids), intent(inout) :: g
      character(len=*), intent(in) :: name
      integer :: proc_id
      call mpi90_rank (proc_id)
      if (proc_id == VAMP_ROOT) then
         call vamp0_write_grids (g%g0, name)
      end if
    end subroutine write_grids_name

    subroutine read_grids_name (g, name)
      type(vamp_grids), intent(inout) :: g
      character(len=*), intent(in) :: name
      integer :: proc_id
      call mpi90_rank (proc_id)
      if (proc_id == VAMP_ROOT) then
         call vamp0_read_grids (g%g0, name)
      end if
    end subroutine read_grids_name
    subroutine vamp_send_grid (g, target, tag, domain, error)
      type(vamp_grid), intent(in) :: g
      integer, intent(in) :: target, tag
      integer, intent(in), optional :: domain
      integer, intent(out), optional :: error
      integer, dimension(2) :: words
      integer, dimension(:), allocatable :: ibuf
      real(kind=default), dimension(:), allocatable :: dbuf
      call vamp_marshal_grid_size (g, words(1), words(2))
      allocate (ibuf(words(1)), dbuf(words(2)))
      call vamp_marshal_grid (g, ibuf, dbuf)
      call mpi90_send (words, target, tag, domain, error)
      call mpi90_send (ibuf, target, tag+1, domain, error)
      call mpi90_send (dbuf, target, tag+2, domain, error)
      deallocate (ibuf, dbuf)
    end subroutine vamp_send_grid

    subroutine vamp_receive_grid (g, source, tag, domain, status, error)
      type(vamp_grid), intent(inout) :: g
      integer, intent(in) :: source, tag
      integer, intent(in), optional :: domain
      type(mpi90_status), intent(out), optional :: status
      integer, intent(out), optional :: error
      integer, dimension(2) :: words
      integer, dimension(:), allocatable :: ibuf
      real(kind=default), dimension(:), allocatable :: dbuf
      call mpi90_receive (words, source, tag, domain, status, error)
      allocate (ibuf(words(1)), dbuf(words(2)))
      call mpi90_receive (ibuf, source, tag+1, domain, status, error)
      call mpi90_receive (dbuf, source, tag+2, domain, status, error)
      call vamp_unmarshal_grid (g, ibuf, dbuf)
      deallocate (ibuf, dbuf)
    end subroutine vamp_receive_grid

    subroutine vamp_broadcast_one_grid (g, root, domain, error)
      type(vamp_grid), intent(inout) :: g
      integer, intent(in) :: root
      integer, intent(in), optional :: domain
      integer, intent(out), optional :: error
      integer, dimension(:), allocatable :: ibuf
      real(kind=default), dimension(:), allocatable :: dbuf
      integer :: iwords, dwords, me
      call mpi90_rank (me)
      if (me == root) then
         call vamp_marshal_grid_size (g, iwords, dwords)
      end if
      call mpi90_broadcast (iwords, root, domain, error)
      call mpi90_broadcast (dwords, root, domain, error)
      allocate (ibuf(iwords), dbuf(dwords))
      if (me == root) then
         call vamp_marshal_grid (g, ibuf, dbuf)
      end if
      call mpi90_broadcast (ibuf, root, domain, error)
      call mpi90_broadcast (dbuf, root, domain, error)
      if (me /= root) then
         call vamp_unmarshal_grid (g, ibuf, dbuf)
      end if
      deallocate (ibuf, dbuf)
    end subroutine vamp_broadcast_one_grid

    subroutine vamp_broadcast_many_grids (g, root, domain, error)
      type(vamp_grid), dimension(:), intent(inout) :: g
      integer, intent(in) :: root
      integer, intent(in), optional :: domain
      integer, intent(out), optional :: error
      integer :: i
      do i = 1, size(g)
         call vamp_broadcast_one_grid (g(i), root, domain, error)
      end do
    end subroutine vamp_broadcast_many_grids

    subroutine vamp_broadcast_grids (g, root, domain, error)
      type(vamp0_grids), intent(inout) :: g
      integer, intent(in) :: root
      integer, intent(in), optional :: domain
      integer, intent(out), optional :: error
      integer :: nch, me
      call mpi90_broadcast (g%sum_chi2, root, domain, error)
      call mpi90_broadcast (g%sum_integral, root, domain, error)
      call mpi90_broadcast (g%sum_weights, root, domain, error)
      call mpi90_rank (me)
      if (me == root) then
         nch = size (g%grids)
      end if
      call mpi90_broadcast (nch, root, domain, error)
      if (me /= root) then
         if (associated (g%grids)) then
            if (size (g%grids) /= nch) then
               call vamp0_delete_grid (g%grids)
               deallocate (g%grids, g%weights, g%num_calls)
               allocate (g%grids(nch), g%weights(nch), g%num_calls(nch))
               call vamp_create_empty_grid (g%grids)
            end if
         else
            allocate (g%grids(nch), g%weights(nch), g%num_calls(nch))
            call vamp_create_empty_grid (g%grids)
         end if
      end if
      call vamp_broadcast_grid (g%grids, root, domain, error)
      call mpi90_broadcast (g%weights, root, domain, error)
      call mpi90_broadcast (g%num_calls, root, domain, error)
    end subroutine vamp_broadcast_grids

    subroutine vamp_send_history (g, target, tag, domain, error)
      type(vamp_history), intent(in) :: g
      integer, intent(in) :: target, tag
      integer, intent(in), optional :: domain
      integer, intent(out), optional :: error
      integer, dimension(2) :: words
      integer, dimension(:), allocatable :: ibuf
      real(kind=default), dimension(:), allocatable :: dbuf
      call vamp_marshal_history_size (g, words(1), words(2))
      allocate (ibuf(words(1)), dbuf(words(2)))
      call vamp_marshal_history (g, ibuf, dbuf)
      call mpi90_send (words, target, tag, domain, error)
      call mpi90_send (ibuf, target, tag+1, domain, error)
      call mpi90_send (dbuf, target, tag+2, domain, error)
      deallocate (ibuf, dbuf)
    end subroutine vamp_send_history

    subroutine vamp_receive_history (g, source, tag, domain, status, error)
      type(vamp_history), intent(inout) :: g
      integer, intent(in) :: source, tag
      integer, intent(in), optional :: domain
      type(mpi90_status), intent(out), optional :: status
      integer, intent(out), optional :: error
      integer, dimension(2) :: words
      integer, dimension(:), allocatable :: ibuf
      real(kind=default), dimension(:), allocatable :: dbuf
      call mpi90_receive (words, source, tag, domain, status, error)
      allocate (ibuf(words(1)), dbuf(words(2)))
      call mpi90_receive (ibuf, source, tag+1, domain, status, error)
      call mpi90_receive (dbuf, source, tag+2, domain, status, error)
      call vamp_unmarshal_history (g, ibuf, dbuf)
      deallocate (ibuf, dbuf)
    end subroutine vamp_receive_history
end module vamp_parallel_mpi
module vampi
  use vamp_serial_mpi !NODEP!
  use vamp_parallel_mpi !NODEP!
  public
end module vampi
