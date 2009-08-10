module timer_class
implicit none

integer, parameter :: dbl = selected_real_kind ( p = 14 )

type, public :: timer

    private
    real(kind=dbl) :: saved_time

    contains

    procedure, public :: start_timer => start_timer_sub
    procedure, public :: elapsed_time => elapsed_time_fn

end type timer

private :: start_timer_sub, elapsed_time_fn

contains

subroutine start_timer_sub (this)
implicit none
    type(timer) :: this
    integer, dimension(8) :: value
    
    call date_and_time ( values = value )
    this%saved_time = 86400.D0 * value(3) + 3600.D0 * value(5) &
        + 60.D0 * value(6) + value(7) + 0.001D0 * value(8)

end subroutine start_timer_sub

real function elapsed_time_fn (this)
implicit none
    type(timer) :: this
    integer, dimension(8) :: value
    real(kind=dbl) :: current_time

    call date_and_time ( values = value )
    current_time = 86400.D0 * value(3) + 3600.D0 * value(5) &
        + 60.D0 * value(6) + value(7) + 0.001D0 * value(8)

    elapsed_time_fn = current_time - this%saved_time

end function elapsed_time_fn

end module timer_class
