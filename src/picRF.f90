program picRF

use iso_fortran_env
use netcdf
use luxury
use random
use timer_class
use possion_solvers

implicit none

!   define types

type :: particle
    real :: x, y, z, vx, vy, vz, q, w
    real :: x_, y_, z_, vx_, vy_, vz_, q_
    integer :: ii 
end type

type :: species
    integer(kind=selected_int_kind(9)) :: nP
    type(particle), dimension(:), pointer :: p
end type species

type :: grid
    integer :: nBins
    real :: mn, mn_
    real :: mx, mx_
    real :: step, step_
    real :: rng, rng_
    real, allocatable :: binEdges(:), binEdges_(:)
    real, allocatable :: binCenters(:), binCenters_(:)
end type grid

!   variable list

type(species) :: H
type(grid) :: x

logical :: two_electron
logical :: two_stream

integer :: iStat, i, t
character(len=5) :: stepChar
real, allocatable :: rhoNGP(:) 
real, allocatable :: rhoNGP_cuda(:)
real, allocatable :: laplaceOp(:,:)
real, allocatable :: phi(:), Ex(:), Ex_cuda(:), phi_jacobi(:)

!   physical constants

real :: e, me, mi, e0, pi
real :: dt, wp, total_charge
real :: xConv, tConv, vConv, aConv, EConv, pConv, rConv
real :: ne, weight

!   lapack variables

integer :: lapack_info, lapack_lda, lapack_ldb, &
    lapack_n, lapack_nrhs
integer, allocatable :: lapack_ipiv(:)
real, allocatable :: lapack_a(:,:), lapack_b(:,:)
integer :: jStat

!   netcdf variables

character(len=30) :: outFName 
integer :: ncid, ncStat
integer :: x_id, nP_id, vx_id, t_id
integer :: nc_start(2), nc_count(2), nc_count_grid(2)
integer :: x_nBins_id, rho_id, phi_id, Ex_id
integer :: dx_id, scalar_id, weight_id
integer :: output_freq, outCnt

!   luxury random number generator vars

integer :: values(8)
character(len=8) :: date
character(len=10) :: time
character(len=5) :: zone
real, allocatable :: randPos(:)

!   random.f90 vars

real, allocatable :: randNormal(:)

!   timing

type(timer) :: clock1, clockHist, clockHistCuda
real(kind=dbl) :: tmpFileTime, tmpFileTime2, fileTime

!   cuda copies

real, allocatable :: cuda_x(:), cuda_vx(:)

e   = 1.602e-19
me  = 9.10938188e-31
mi  = 1.67262158e-27
e0  = 8.85418782e-12
pi  = 3.14

i=0
t=0

!   initialise variables

!   periodic grid

x%mn  = 0.0
x%mx  = 1.0e-06
x%nBins   = 256 
x%rng    = x%mx-x%mn
x%step  = x%rng/x%nBins
allocate(x%binCenters(x%nBins),x%binCenters_(x%nBins))
x%binCenters  = (/ (i*x%step,i=0,x%nBins) /) 

allocate(rhoNGP(x%nBins), stat = iStat)
allocate(rhoNGP_cuda(x%nBins), stat = iStat)

lapack_n    = x%nBins
lapack_nRhs = 1
lapack_lda  = x%nBins
lapack_ldb  = x%nBins

allocate ( laplaceOp(x%nBins,x%nBins) )
allocate ( lapack_a(lapack_lda,x%nBins))
allocate ( lapack_ipiv(lapack_n))
allocate ( lapack_b(lapack_ldb,lapack_nRhs))


allocate( phi(x%nBins), &
          Ex(x%nBins), &
          Ex_cuda(x%nBins), &
          phi_jacobi(x%nBins) )

rhoNGP  = 0.0
rhoNGP_cuda = 0
phi = 0.0
phi_jacobi  = 0.0
Ex  = 0.0
Ex_cuda = 0.0

!   benchmarking scenarios
!   ----------------------

two_electron = .false.
two_stream = .true.

two_electron_oscillation: &
if (two_electron) then

    H%nP    = 2 
    allocate ( H%p(H%nP), stat = iStat )
    if ( iStat > 0 ) stop
    
    H%p%x   = (/ (i * x%rng/4,i=0,H%nP-1) /) + x%mn + 0.22
    H%p%y   = 0.0
    H%p%z   = 0.0
    H%p%vx  = 0.0
    H%p%vy  = 0.0
    H%p%vz  = 0.0
    
    H%p%q    = -1.0 * e
    
    ne  = 1e12
    
    !H%p%w   = ne * x%rng / H%nP 
    weight   = ne * x%rng / H%nP 

    
    wp  = sqrt ( ne * e**2 / ( me * e0 ) )

    dt  = 0.15 / wp *1e-4
    output_freq = 1
    write(*,*) 'plasma freq: ', wp
    write(*,*) 'dt: ', dt
    
endif two_electron_oscillation


two_stream_instability: &
if (two_stream) then

    H%nP    = 25600
    output_freq = 1

    allocate ( H%p(H%nP), stat = iStat )
    if ( iStat > 0 ) stop

    allocate ( randPos(H%nP), stat = iStat )
    allocate ( randNormal(H%nP), stat = iStat )
    if ( iStat > 0 ) stop

    !   create seed
    call date_and_time ( date, time, zone, values )
    !   initialise generator with seed
    call rLuxGo ( 223, abs ( values(8)*values(7) ) + 1, 0, 0)
    !   use generator
    call ranLux ( randPos, H%nP )

    do i=1,H%nP
        randNormal(i)   = random_normal()
    enddo

    H%p%x   = randPos * x%rng 
    H%p%y   = 0.0
    H%p%z   = 0.0
    H%p(1:H%nP/2)%vx  = randNormal(1:H%nP/2)*5e2+5e3 
    H%p(H%nP/2+1:H%nP)%vx = randNormal(H%nP/2+1:H%nP)*5e2-5e3
    H%p%vy  = 0.0
    H%p%vz  = 0.0
    
    H%p%q    = -1.0 * e
    
    ne  = 1e19
    
    !H%p%w   = ne * x%rng / H%nP 
    weight   = ne * x%rng / H%nP 
    
    wp  = sqrt ( ne * e**2 / ( me * e0 ) )

    dt  = 1.0 / wp * 0.1
    
    write(*,*) 'plasma freq: ', wp
    write(*,*) 'dt: ', dt

endif two_stream_instability

!   allocate cuda copies

    allocate ( cuda_x ( H%nP ), cuda_vx ( H%nP ) )
    cuda_x  = H%p%x
    cuda_vx = H%p%vx

!   setup dimensionless conversions
!   -----------------------------

xConv   = 1.0 / x%step
tConv   = 1.0 / dt
vConv   = dt / x%step
aConv   = dt**2 / x%step
EConv   = e * dt**2 / ( me * x%step )
pConv   = e * dt**2 / ( 2.0 * me * x%step**2 )
rConv   = e * dt**2 / ( 2.0 * me * e0 )

!   setup laplace operator
!   ----------------------

laplaceOp(1,1)  =-2.0
laplaceOp(2,1)  = 1.0

laplaceOp(1,x%nBins)    = 1.0
laplaceOp(x%nBins-1,x%nBins)    = 1.0
laplaceOp(x%nBins,x%nBins)    =-2.0

do i=2,x%nBins-1
    laplaceOp(i-1,i)    = 1.0
    laplaceOp(i,i)      = -2.0
    laplaceOp(i+1,i)    = 1.0
enddo
lapack_a    = laplaceOp

!   convert inital conditions to dimensionless vars
!   -----------------------------------------------

H%p%x_   = H%p%x * xConv
H%p%y_   = H%p%y * xConv
H%p%z_   = H%p%z * xConv
H%p%vx_  = H%p%vx * vConv
H%p%vy_  = H%p%vy * vConv
H%p%vz_  = H%p%vz * vConv

x%binCenters_   = x%binCenters * xConv
x%rng_  = x%rng * xConv

!   initialise ouput file
!   ----------------

outFName = "output/output.nc"
ncStat  = nf90_create ( outFName, NF90_CLOBBER, ncid )
ncStat  = nf90_def_dim ( ncid, 'nP', H%nP, np_id )
ncStat  = nf90_def_dim ( ncid, 'x_nBins', x%nBins, x_nBins_id )
ncStat  = nf90_def_dim ( ncid, 't', NF90_UNLIMITED, t_id )
ncStat  = nf90_def_dim ( ncid, 'scalar', 1, scalar_id )

ncStat  = nf90_def_var ( ncid, 'x', NF90_REAL, (/np_id,t_id/), x_id )
ncStat  = nf90_def_var ( ncid, 'vx', NF90_REAL, (/np_id,t_id/), vx_id )
ncStat  = nf90_def_var ( ncid, 'weight', NF90_REAL, (/np_id,t_id/), weight_id )

ncStat  = nf90_def_var ( ncid, 'rho', NF90_REAL, (/x_nBins_id,t_id/), rho_id )
ncStat  = nf90_def_var ( ncid, 'phi', NF90_REAL, (/x_nBins_id,t_id/), phi_id )
ncStat  = nf90_def_var ( ncid, 'Ex', NF90_REAL, (/x_nBins_id,t_id/), Ex_id )
ncStat  = nf90_def_var ( ncid, 'dx', NF90_REAL, scalar_id, dx_id )
ncStat  = nf90_enddef ( ncid )

nc_start   = (/ 1, 1 /)
nc_count   = (/ H%nP, 1 /) 
nc_count_grid   = (/ x%nBins, 1 /)
ncStat  = nf90_put_var ( ncid, dx_id, x%step )

outCnt = 0

call clock1%start_timer()

time_loop: &
do t=1,1
    !write(*,'(i4)', advance = 'no') t 
    !flush ( output_unit )
 
    !   calculate charge density
    !   ------------------------
   
    call clockHist%start_timer()

    rhoNGP  = 0.0
    do i=1,H%nP 
    
        H%p(i)%ii   = & 
                nInt ( &
                    (H%p(i)%x-x%binCenters(1)) / x%rng * x%nBins &
                     +1) 
        if (H%p(i)%ii .eq. x%nBins+1) H%p(i)%ii = 1

        rhoNGP(H%p(i)%ii) = rhoNGP(H%p(i)%ii) + 1
    
    enddo

    rhoNGP = (-1.0) * e * weight / x%step * rhoNGP 

    
    write(*,'(a30,2x,f5.1)') 'hist time: ', clockHist%elapsed_time()

    !   enforce charge neutrality
    !   -------------------------

    total_charge    = sum ( rhoNGP )
    rhoNGP  = rhoNGP - total_charge / x%nBins 
    !write(*,*) 'CPU total:'
    !write(*,'(f15.7)') total_charge 
    !! compare the GPU and CPU versions

    !do i=1,x%nBins
    !    write(*,*) rhoNGP(i), rhoNGP_cuda(i)
    !enddo
    !write(*,*) '----------------------------' 
    !write(*,*) sum(rhoNGP), sum(rhoNGP_cuda)

   
    
    !   sovle laplace eqn
    !   -----------------

    lapack_a    = laplaceOp
    lapack_b(:,1)    = -rhoNGP(1:x%nBins) / e0 * x%step**2
    
    call sgesv ( lapack_n, &
                lapack_nRhs, &
                lapack_a, &
                lapack_lda, &
                lapack_ipiv, &
                lapack_b, &
                lapack_ldb, &
                lapack_info )
    
    check_lapack_result: &
    if ( lapack_info .ne. 0 ) then
        write(*,*) 'picRF.f90 [153]: lapack error', lapack_info
        stop
    endif check_lapack_result
    
    phi = lapack_b(:,1)
  
    
    !!   try jacobi iterative method
    !!   ---------------------------

    !phi_jacobi  = 0
    !jStat = jacobi_iter ( laplaceOp, -rhoNGP(1:x%nBins) / e0 * x%step**2, x = phi_jacobi ) 


    !   calculate Ex
    !   ------------
    
    Ex(1)   = ( phi(x%nBins) - phi(2) ) / ( 2.0 * x%step )
    Ex(x%nBins)   = ( phi(x%nBins-1) - phi(1) ) / ( 2.0 * x%step )
    
    calculate_Ex: &
    do i=2,x%nBins-1
        Ex(i)   = ( phi(i-1) - phi(i+1) ) / ( 2.0 * x%step )
    enddo calculate_Ex


!    !   print debug info
!    !   ----------------
!    
!    write(*,'(4a10)'), 'x', 'rho', 'phi', 'Ex'
!    do i=1,x%nBins
!        write(*,'(4f10.2)') x%binCenters(i), rhoNGP(i)*1e6, phi(i), Ex(i)
!    enddo
!    write(*,*) 'Total charge: ', sum ( rhoNGP )

    !   move particles
    !   --------------
    
    move_particles: &
    do i=1,H%nP
    
        H%p(i)%vx   = H%p(i)%vx + dt / ( me ) * H%p(i)%q * Ex(H%p(i)%ii)
        H%p(i)%x    = H%p(i)%x + H%p(i)%vx * dt

        !   periodic particle boundaries
        if ( H%p(i)%x < x%mn ) H%p(i)%x = H%p(i)%x + x%rng
        if ( H%p(i)%x .ge. x%mx ) H%p(i)%x = H%p(i)%x - x%rng

        dx_too_large: & 
        if ( H%p(i)%x < x%mn .or. H%p(i)%x .gt. x%mx ) then
              write(*,*) 'picRF.f90 [275]:' 
              write(*,*) 'ERROR: reduce dt', &
                (H%p(i)%x-x%mn)/x%rng, &
                dt / ( me ) * H%p(i)%q * Ex(H%p(i)%ii)
              stop
        endif dx_too_large

    enddo move_particles
    
    !stop
     !   write data
    !   ----------

    !write(*,*) H%p%vx

    if (mod(t,output_freq) .eq. 0) then

        write(*,*) t
        tmpFileTime = clock1%elapsed_time()

        outCnt  = outCnt + 1
        nc_start(2) = outCnt 
        ncStat  = nf90_put_var ( ncid, x_id, H%p%x, &
            start = nc_start, count = nc_count )
        ncStat  = nf90_put_var ( ncid, vx_id, H%p%vx, &
            start = nc_start, count = nc_count )
        ncStat  = nf90_put_var ( ncid, weight_id, H%p%w, &
            start = nc_start, count = nc_count )

        ncStat  = nf90_put_var ( ncid, rho_id, rhoNGP, &
            start = nc_start, count = nc_count_grid )
        ncStat  = nf90_put_var ( ncid, phi_id, phi, &
            start = nc_start, count = nc_count_grid )
        ncStat  = nf90_put_var ( ncid, Ex_id, Ex, &
            start = nc_start, count = nc_count_grid )

        tmpFileTime2 = clock1%elapsed_time()
        fileTime = fileTime + (tmpFileTime2 - tmpFileTime)
    endif

enddo time_loop

! call cuda version

call clockHistCuda%start_timer()
call cudahist(cuda_x,cuda_vx,H%nP,rhoNGP_cuda,x%nBins,x%binCenters(1),&
    x%mx,x%rng,weight,x%step,Ex_cuda,dt)
write(*,'(a30,2x,f5.1)') 'hist cuda time: ', clockHistCuda%elapsed_time()


do i=1,x%nBins
    write(*,*) phi(i), rhoNGP_cuda(i), -rhoNGP(i) / e0 * x%step**2
enddo
write(*,*) '----------------------------' 
write(*,*) sum ( phi ), sum ( phi_jacobi )

do i=1,x%nBins
    write(*,*) Ex(i), Ex_cuda(i)
enddo
write(*,*) '----------------------------' 
write(*,*) sum ( Ex ), sum ( Ex_cuda )


ncStat  = nf90_close ( ncid )

write(*,'(a30,2x,f5.1)') 'Elapsed time: ', clock1%elapsed_time()
write(*,'(a30,2x,f5.1)') 'netCDF time: ',fileTime 

end program picRF
