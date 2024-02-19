! mpi90.f90 --
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
module mpi90
  use kinds
  use mpi
  implicit none
  private
  public :: mpi90_init
  public :: mpi90_finalize
  public :: mpi90_abort
  public :: mpi90_print_error
  public :: mpi90_size
  public :: mpi90_rank
  public :: mpi90_send
  public :: mpi90_receive
  public :: mpi90_receive_pointer
  private :: decode_status
  public :: mpi90_broadcast
  interface mpi90_send
     module procedure &
          mpi90_send_integer, mpi90_send_double, &
          mpi90_send_integer_array, mpi90_send_double_array, &
          mpi90_send_integer_array2, mpi90_send_double_array2
  end interface
  interface mpi90_receive
     module procedure &
          mpi90_receive_integer, mpi90_receive_double, &
          mpi90_receive_integer_array, mpi90_receive_double_array, &
          mpi90_receive_integer_array2, mpi90_receive_double_array2
  end interface
  interface mpi90_receive_pointer
     module procedure &
          mpi90_receive_integer_pointer, mpi90_receive_double_pointer
  end interface
  interface mpi90_broadcast
     module procedure &
          mpi90_broadcast_integer, mpi90_broadcast_integer_array, &
          mpi90_broadcast_integer_array2, mpi90_broadcast_integer_array3, &
          mpi90_broadcast_double, mpi90_broadcast_double_array, &
          mpi90_broadcast_double_array2, mpi90_broadcast_double_array3, &
          mpi90_broadcast_logical, mpi90_broadcast_logical_array, &
          mpi90_broadcast_logical_array2, mpi90_broadcast_logical_array3
  end interface
  
  
  type, public :: mpi90_status
     integer :: count, source, tag, error
  end type mpi90_status
  character(len=*), public, parameter :: MPI90_RCS_ID = &
       "$Id: mpi90.nw 314 2010-04-17 20:32:33Z ohl $"
contains
  subroutine mpi90_init (error)
    integer, intent(out), optional :: error
    integer :: local_error
    character(len=*), parameter :: FN = "mpi90_init"
    external mpi_init
    call mpi_init (local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          stop
       end if
    end if
  end subroutine mpi90_init
  subroutine mpi90_finalize (error)
    integer, intent(out), optional :: error
    integer :: local_error
    character(len=*), parameter :: FN = "mpi90_finalize"
    external mpi_finalize
    call mpi_finalize (local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
  end subroutine mpi90_finalize
  subroutine mpi90_abort (code, domain, error)
    integer, intent(in), optional :: code, domain
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_abort"
    integer :: local_domain, local_code, local_error
    external mpi_abort
    if (present (code)) then
       local_code = code
    else
       local_code = MPI_ERR_UNKNOWN
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_abort (local_domain, local_code, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          stop
       end if
    end if
  end subroutine mpi90_abort
  subroutine mpi90_print_error (error, msg)
    integer, intent(in) :: error
    character(len=*), optional :: msg
    character(len=*), parameter :: FN = "mpi90_print_error"
    integer :: msg_len, local_error
    external mpi_error_string
    call mpi_error_string (error, msg, msg_len, local_error)
    if (local_error /= MPI_SUCCESS) then
       print *, "PANIC: even MPI_ERROR_STRING() failed!!!"
       call mpi90_abort (local_error)
    else if (present (msg)) then
       print *, trim (msg), ": ", trim (msg(msg_len+1:))
    else
       print *, "mpi90: ", trim (msg(msg_len+1:))
    end if
  end subroutine mpi90_print_error
  subroutine mpi90_size (sz, domain, error)
    integer, intent(out) :: sz
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_size"
    integer :: local_domain, local_error
    external mpi_comm_size
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_comm_size (local_domain, sz, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
  end subroutine mpi90_size
  subroutine mpi90_rank (rank, domain, error)
    integer, intent(out) :: rank
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_rank"
    integer :: local_domain, local_error
    external mpi_comm_rank
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_comm_rank (local_domain, rank, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
  end subroutine mpi90_rank
  subroutine mpi90_send_integer (value, target, tag, domain, error)
    integer, intent(in) :: value
    integer, intent(in) :: target, tag
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    call mpi90_send_integer_array ((/ value /), target, tag, domain, error)
  end subroutine mpi90_send_integer
  subroutine mpi90_send_double (value, target, tag, domain, error)
    real(kind=default), intent(in) :: value
    integer, intent(in) :: target, tag
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    call mpi90_send_double_array ((/ value /), target, tag, domain, error)
  end subroutine mpi90_send_double
  subroutine mpi90_send_integer_array (buffer, target, tag, domain, error)
    integer, dimension(:), intent(in) :: buffer
    integer, intent(in) ::  target, tag
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_send_integer_array"
    integer, parameter :: datatype = MPI_INTEGER
    integer :: local_domain, local_error
    external mpi_send
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_send (buffer, size (buffer), datatype, target, tag, &
                   local_domain, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
  end subroutine mpi90_send_integer_array
  subroutine mpi90_send_double_array (buffer, target, tag, domain, error)
    real(kind=default), dimension(:), intent(in) :: buffer
    integer, intent(in) :: target, tag
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_send_double_array"
    integer, parameter :: datatype = MPI_DOUBLE_PRECISION
    integer :: local_domain, local_error
    external mpi_send
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_send (buffer, size (buffer), datatype, target, tag, &
                   local_domain, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
  end subroutine mpi90_send_double_array
  subroutine mpi90_send_integer_array2 (value, target, tag, domain, error)
    integer, dimension(:,:), intent(in) :: value
    integer, intent(in) :: target, tag
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    integer, dimension(size(value)) :: buffer
    buffer = reshape (value, shape(buffer))
    call mpi90_send_integer_array (buffer, target, tag, domain, error)
  end subroutine mpi90_send_integer_array2
  subroutine mpi90_send_double_array2 (value, target, tag, domain, error)
    real(kind=default), dimension(:,:), intent(in) :: value
    integer, intent(in) :: target, tag
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    real(kind=default), dimension(size(value)) :: buffer
    buffer = reshape (value, shape(buffer))
    call mpi90_send_double_array (buffer, target, tag, domain, error)
  end subroutine mpi90_send_double_array2
  subroutine mpi90_receive_integer (value, source, tag, domain, status, error)
    integer, intent(out) :: value
    integer, intent(in), optional :: source, tag, domain
    type(mpi90_status), intent(out), optional :: status
    integer, intent(out), optional :: error
    integer, dimension(1)  :: buffer
    call mpi90_receive_integer_array (buffer, source, tag, domain, status, error)
    value = buffer(1)
  end subroutine mpi90_receive_integer
  subroutine mpi90_receive_double (value, source, tag, domain, status, error)
    real(kind=default), intent(out) :: value
    integer, intent(in), optional :: source, tag, domain
    type(mpi90_status), intent(out), optional :: status
    integer, intent(out), optional :: error
    real(kind=default), dimension(1) :: buffer
    call mpi90_receive_double_array (buffer, source, tag, domain, status, error)
    value = buffer(1)
  end subroutine mpi90_receive_double
  subroutine mpi90_receive_integer_array &
       (buffer, source, tag, domain, status, error)
    integer, dimension(:), intent(out) :: buffer
    integer, intent(in), optional :: source, tag, domain
    type(mpi90_status), intent(out), optional :: status
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_receive_integer_array"
    integer, parameter :: datatype = MPI_INTEGER
    integer :: local_source, local_tag, local_domain, local_error
    integer, dimension(MPI_STATUS_SIZE) :: local_status
    external mpi_receive, mpi_get_count
    if (present (source)) then
       local_source = source
    else
       local_source = MPI_ANY_SOURCE
    end if
    if (present (tag)) then
       local_tag = tag
    else
       local_tag = MPI_ANY_TAG
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_recv (buffer, size (buffer), datatype, local_source, local_tag, &
                   local_domain, local_status, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
    if (present (status)) then
       call decode_status (status, local_status, datatype)
    end if
  end subroutine mpi90_receive_integer_array
  subroutine decode_status (status, mpi_status, datatype)
    type(mpi90_status), intent(out) :: status
    integer, dimension(:), intent(in) :: mpi_status
    integer, intent(in), optional :: datatype
    integer :: ierror
    if (present (datatype)) then
       call mpi_get_count (mpi_status, datatype, status%count, ierror)
    else
       status%count = 0
    end if
    status%source = mpi_status(MPI_SOURCE)
    status%tag = mpi_status(MPI_TAG)
    status%error = mpi_status(MPI_ERROR)
  end subroutine decode_status
  subroutine mpi90_receive_double_array &
       (buffer, source, tag, domain, status, error)
    real(kind=default), dimension(:), intent(out) :: buffer
    integer, intent(in), optional :: source, tag, domain
    type(mpi90_status), intent(out), optional :: status
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_receive_double_array"
    integer, parameter :: datatype = MPI_DOUBLE_PRECISION
    integer :: local_source, local_tag, local_domain, local_error
    integer, dimension(MPI_STATUS_SIZE) :: local_status
    external mpi_receive, mpi_get_count
    if (present (source)) then
       local_source = source
    else
       local_source = MPI_ANY_SOURCE
    end if
    if (present (tag)) then
       local_tag = tag
    else
       local_tag = MPI_ANY_TAG
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_recv (buffer, size (buffer), datatype, local_source, local_tag, &
                   local_domain, local_status, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
    if (present (status)) then
       call decode_status (status, local_status, datatype)
    end if
  end subroutine mpi90_receive_double_array
  subroutine mpi90_receive_integer_array2 &
       (value, source, tag, domain, status, error)
    integer, dimension(:,:), intent(out) :: value
    integer, intent(in), optional :: source, tag, domain
    type(mpi90_status), intent(out), optional :: status
    integer, intent(out), optional :: error
    integer, dimension(size(value)) :: buffer
    call mpi90_receive_integer_array &
         (buffer, source, tag, domain, status, error)
    value = reshape (buffer, shape(value))
  end subroutine mpi90_receive_integer_array2
  subroutine mpi90_receive_double_array2 &
       (value, source, tag, domain, status, error)
    real(kind=default), dimension(:,:), intent(out) :: value
    integer, intent(in), optional :: source, tag, domain
    type(mpi90_status), intent(out), optional :: status
    integer, intent(out), optional :: error
    real(kind=default), dimension(size(value)) :: buffer
    call mpi90_receive_double_array &
          (buffer, source, tag, domain, status, error)
    value = reshape (buffer, shape(value))
  end subroutine mpi90_receive_double_array2
  subroutine mpi90_receive_integer_pointer &
       (buffer, source, tag, domain, status, error)
    integer, dimension(:), pointer :: buffer
    integer, intent(in), optional :: source, tag, domain
    type(mpi90_status), intent(out), optional :: status
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_receive_integer_pointer"
    integer, parameter :: datatype = MPI_INTEGER
    integer :: local_source, local_tag, local_domain, local_error, buffer_size
    integer, dimension(MPI_STATUS_SIZE) :: local_status
    integer :: ierror
    external mpi_receive, mpi_get_count
    if (present (source)) then
       local_source = source
    else
       local_source = MPI_ANY_SOURCE
    end if
    if (present (tag)) then
       local_tag = tag
    else
       local_tag = MPI_ANY_TAG
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_probe (local_source, local_tag, local_domain, &
                    local_status, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
    call mpi_get_count (local_status, datatype, buffer_size, ierror)
    if (associated (buffer)) then
       if (size (buffer) /= buffer_size) then
          deallocate (buffer)
          allocate (buffer(buffer_size))
       end if
    else
       allocate (buffer(buffer_size))
    end if
    call mpi_recv (buffer, size (buffer), datatype, local_source, local_tag, &
                   local_domain, local_status, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
    if (present (status)) then
       call decode_status (status, local_status, datatype)
    end if
  end subroutine mpi90_receive_integer_pointer
  subroutine mpi90_receive_double_pointer &
       (buffer, source, tag, domain, status, error)
    real(kind=default), dimension(:), pointer :: buffer
    integer, intent(in), optional :: source, tag, domain
    type(mpi90_status), intent(out), optional :: status
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_receive_double_pointer"
    integer, parameter :: datatype = MPI_DOUBLE_PRECISION
    integer :: local_source, local_tag, local_domain, local_error, buffer_size
    integer, dimension(MPI_STATUS_SIZE) :: local_status
    integer :: ierror
    external mpi_receive, mpi_get_count
    if (present (source)) then
       local_source = source
    else
       local_source = MPI_ANY_SOURCE
    end if
    if (present (tag)) then
       local_tag = tag
    else
       local_tag = MPI_ANY_TAG
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_probe (local_source, local_tag, local_domain, &
                    local_status, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
    call mpi_get_count (local_status, datatype, buffer_size, ierror)
    if (associated (buffer)) then
       if (size (buffer) /= buffer_size) then
          deallocate (buffer)
          allocate (buffer(buffer_size))
       end if
    else
       allocate (buffer(buffer_size))
    end if
    call mpi_recv (buffer, size (buffer), datatype, local_source, local_tag, &
                   local_domain, local_status, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
    if (present (status)) then
       call decode_status (status, local_status, datatype)
    end if
  end subroutine mpi90_receive_double_pointer
  subroutine mpi90_broadcast_integer (value, root, domain, error)
    integer, intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    integer, dimension(1) :: buffer
    buffer(1) = value
    call mpi90_broadcast_integer_array (buffer, root, domain, error)
    value = buffer(1)
  end subroutine mpi90_broadcast_integer
  subroutine mpi90_broadcast_double (value, root, domain, error)
    real(kind=default), intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    real(kind=default), dimension(1) :: buffer
    buffer(1) = value
    call mpi90_broadcast_double_array (buffer, root, domain, error)
    value = buffer(1)
  end subroutine mpi90_broadcast_double
  subroutine mpi90_broadcast_logical (value, root, domain, error)
    logical, intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    logical, dimension(1) :: buffer
    buffer(1) = value
    call mpi90_broadcast_logical_array (buffer, root, domain, error)
    value = buffer(1)
  end subroutine mpi90_broadcast_logical
  subroutine mpi90_broadcast_integer_array (buffer, root, domain, error)
    integer, dimension(:), intent(inout) :: buffer
    integer, intent(in) ::  root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    character(len=*), parameter :: FN = "mpi90_broadcast_integer_array"
    integer, parameter :: datatype = MPI_INTEGER
    integer :: local_domain, local_error
    external mpi_bcast
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_bcast (buffer, size (buffer), datatype, root, &
                    local_domain, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
  end subroutine mpi90_broadcast_integer_array
  subroutine mpi90_broadcast_double_array (buffer, root, domain, error)
    real(kind=default), dimension(:), intent(inout) :: buffer
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    integer, parameter :: datatype = MPI_DOUBLE_PRECISION
    character(len=*), parameter :: FN = "mpi90_broadcast_double_array"
    integer :: local_domain, local_error
    external mpi_bcast
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_bcast (buffer, size (buffer), datatype, root, &
                    local_domain, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
  end subroutine mpi90_broadcast_double_array
  subroutine mpi90_broadcast_logical_array (buffer, root, domain, error)
    logical, dimension(:), intent(inout) :: buffer
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    integer, parameter :: datatype = MPI_LOGICAL
    character(len=*), parameter :: FN = "mpi90_broadcast_logical_array"
    integer :: local_domain, local_error
    external mpi_bcast
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    if (present (domain)) then
       local_domain = domain
    else
       local_domain = MPI_COMM_WORLD
    end if
    call mpi_bcast (buffer, size (buffer), datatype, root, &
                    local_domain, local_error)
    if (present (error)) then
       error = local_error
    else
       if (local_error /= MPI_SUCCESS) then
          call mpi90_print_error (local_error, FN)
          call mpi90_abort (local_error)
          stop
       end if
    end if
  end subroutine mpi90_broadcast_logical_array
  subroutine mpi90_broadcast_integer_array2 (value, root, domain, error)
    integer, dimension(:,:), intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    integer, dimension(size(value)) :: buffer
    buffer = reshape (value, shape(buffer))
    call mpi90_broadcast_integer_array (buffer, root, domain, error)
    value = reshape (buffer, shape(value))
  end subroutine mpi90_broadcast_integer_array2
  subroutine mpi90_broadcast_double_array2 (value, root, domain, error)
    real(kind=default), dimension(:,:), intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    real(kind=default), dimension(size(value)) :: buffer
    buffer = reshape (value, shape(buffer))
    call mpi90_broadcast_double_array (buffer, root, domain, error)
    value = reshape (buffer, shape(value))
  end subroutine mpi90_broadcast_double_array2
  subroutine mpi90_broadcast_logical_array2 (value, root, domain, error)
    logical, dimension(:,:), intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    logical, dimension(size(value)) :: buffer
    buffer = reshape (value, shape(buffer))
    call mpi90_broadcast_logical_array (buffer, root, domain, error)
    value = reshape (buffer, shape(value))
  end subroutine mpi90_broadcast_logical_array2
  subroutine mpi90_broadcast_integer_array3 (value, root, domain, error)
    integer, dimension(:,:,:), intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    integer, dimension(size(value)) :: buffer
    buffer = reshape (value, shape(buffer))
    call mpi90_broadcast_integer_array (buffer, root, domain, error)
    value = reshape (buffer, shape(value))
  end subroutine mpi90_broadcast_integer_array3
  subroutine mpi90_broadcast_double_array3 (value, root, domain, error)
    real(kind=default), dimension(:,:,:), intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    real(kind=default), dimension(size(value)) :: buffer
    buffer = reshape (value, shape(buffer))
    call mpi90_broadcast_double_array (buffer, root, domain, error)
    value = reshape (buffer, shape(value))
  end subroutine mpi90_broadcast_double_array3
  subroutine mpi90_broadcast_logical_array3 (value, root, domain, error)
    logical, dimension(:,:,:), intent(inout) :: value
    integer, intent(in) :: root
    integer, intent(in), optional :: domain
    integer, intent(out), optional :: error
    logical, dimension(size(value)) :: buffer
    buffer = reshape (value, shape(buffer))
    call mpi90_broadcast_logical_array (buffer, root, domain, error)
    value = reshape (buffer, shape(value))
  end subroutine mpi90_broadcast_logical_array3
end module mpi90
