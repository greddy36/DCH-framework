 module kinds
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  implicit none
  private
 
  !!! available REAL kinds               ! prec.  ! ISO     ! C
  integer, parameter :: single    =  4   !  1.. 6 ! real32  ! c_float      
  integer, parameter :: double    =  8   !  7..15 ! real64  ! c_double     
  integer, parameter :: extended  = 10   ! 16..18 ! -1      ! c_long_double
  integer, parameter :: quadruple = 16   ! 19..33 ! real128 ! c_float128   
 
  !!! available INTEGER kinds            ! range  ! ISO     ! C
  integer, parameter :: dflt_int  =  4   !  5.. 9 ! int32   ! 
  integer, parameter :: range02   =  1   !  1.. 2 ! int8    ! 
  integer, parameter :: range04   =  2   !  3.. 4 ! int16   ! 
  integer, parameter :: range18   =  8   ! 10..18 ! int64   ! 
  integer, parameter :: range38   = 16   ! 19..38 ! -1      ! 
 
  !!! additional INTEGER kinds
  public :: i8, i16, i32, i64
  integer, parameter :: i8  = selected_int_kind (2)
  integer, parameter :: i16 = selected_int_kind (4)
  integer, parameter :: i32 = selected_int_kind (9)
  integer, parameter :: i64 = selected_int_kind (18)
  public :: TC
  integer, parameter :: TC = i32
 
  !!! default REAL kinds
  public :: single, double
  public :: default, c_default_float, c_default_complex
  integer, parameter :: default           = double   
  integer, parameter :: c_default_float   = c_double     
  integer, parameter :: c_default_complex = c_double_complex     
 
 end module kinds
