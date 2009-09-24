module possion_solvers

contains

function jacobi_iter ( a, b, x)

    implicit none
    
    real, intent(in) :: a(:,:), b(:)
    real, intent(inout), optional, allocatable :: x(:)
    real :: tol
    real, allocatable :: x0(:)
    integer :: jacobi_iter
    real :: norm0, sum_

    integer :: nPts, nIter, nIterMax
    integer :: i, j

    nPts    = size ( b )
    nIterMax   = 1000
    tol = 0.01

    allocate ( x0 ( nPts ) )
    x0  = 0

    if ( .not. allocated ( x ) ) allocate ( x ( nPts ) )
    x   = 0


    nIter = 0
    iteration: &
    do
        
        do i = 1, nPts

            sum_    = 0.0
            do j = 1, nPts

                if (i /= j ) sum_ = sum_ - ( a(i,j) * x0(j) )

            enddo
            x(i)    = ( sum_ + b(i) ) / a(i,i)

        enddo

        norm0   = abs ( sum ( x - x0 ) ) / abs ( sum ( x ) )
        
        !write(*,*) nIter, norm0 
    
        if ( norm0 < tol .and. nIter > 0 ) then
            jacobi_iter = 0
            !write(*,*) 'Convergence ... yeah baby!', nIter, abs ( sum ( x - x0 ) )
            deallocate ( x0 )
            exit
        endif

        if ( nIter > nIterMax ) then
            jacobi_iter = 1
            write(*,*) 'Soln did not converge, you are stupid'
            deallocate ( x0 )
            exit
        endif

        nIter   = nIter + 1

        x0  = x
        !write(*,*) x0

    enddo iteration

end function jacobi_iter

end module possion_solvers
