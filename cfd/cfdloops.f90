program cfd

  use mpi

  implicit none

  double precision, allocatable ::  psi(:,:), zet(:,:), vel(:,:,:)
  double precision, allocatable ::  psitmp(:,:), zettmp(:,:)

  integer, allocatable :: rgb(:,:,:)

  integer :: scalefactor,  numiter

  integer, parameter :: maxline = 32
  character(len=maxline) :: tmparg

!  Basic sizes of simulation

  integer, parameter :: bbase = 10
  integer, parameter :: hbase = 15
  integer, parameter :: wbase =  5
  integer, parameter :: mbase = 32
  integer, parameter :: nbase = 32

  double precision :: re
  logical :: irrotational

!  re = 3.7 seems to be stability limit with Jacobi

!  Some auxiliary parameters and variables

  integer :: m, n, ln, b, h, w
  integer :: iter, i, j

  double precision :: tstart, tstop, ttot, titer, modvsq, hue

!  Variables needed for parallelisation

  integer :: comm, rank, size, ierr

!  Parallel initialisation

  comm = MPI_COMM_WORLD

  call mpi_init(ierr)

  call mpi_comm_rank(comm, rank, ierr)
  call mpi_comm_size(comm, size, ierr)

!  Read in parameters

  if (command_argument_count() /= 2 .and. command_argument_count() /= 3) then

     if (rank == 0) write(*,*) 'Usage: cfd <scale> <numiter> [reynolds]'
     call mpi_finalize(ierr)
     stop

  end if

  if (rank == 0) then 

     call get_command_argument(1, tmparg)
     read(tmparg,*) scalefactor
     call get_command_argument(2, tmparg)
     read(tmparg,*) numiter

     if (command_argument_count() == 3) then

        irrotational = .false.
        call get_command_argument(3, tmparg)
        read(tmparg,*) re
        
     else

        re = -1.0
        irrotational = .true.

     end if

     write(*,fmt='('' Scale factor = '',i3,'', iterations = '', i6)') &
           scalefactor, numiter

     if (irrotational) then
        
        write(*,*) 'Irrotational flow'
        
     else

        write(*,fmt='('' Reynolds number = '', f6.3)') re
        
     end if

  end if

!  Broadcast runtime parameters to all the other processors
  
  call mpi_bcast(scalefactor,  1, MPI_INTEGER, 0, comm, ierr)
  call mpi_bcast(numiter,      1, MPI_INTEGER, 0, comm, ierr)
  call mpi_bcast(re, 1,  MPI_DOUBLE_PRECISION, 0, comm, ierr)
  call mpi_bcast(irrotational, 1, MPI_LOGICAL, 0, comm, ierr)

!  Calculate b, h & w and m & n
        
  b = bbase*scalefactor 
  h = hbase*scalefactor
  w = wbase*scalefactor 
  m = mbase*scalefactor
  n = nbase*scalefactor

  re = re / dble(scalefactor)

!  Calculate local size

  ln = n/size

!  Consistency check

  if (size*ln /= n) then

     if (rank == 0) write(*,*) 'ERROR: n = ', n, ' does not divide onto ', &
          size, ' processes'

     call mpi_finalize(ierr)
     stop

  end if

  if (rank == 0) then
     write(*,fmt='('' Running CFD on '', i4, '' x '', i4, '' grid using '', &
                  &i4, '' process(es)'')') m, n, size
  end if

!  Allocate arrays, including halos on psi and tmp

  allocate(psi(0:m+1, 0:ln+1))
  allocate(zet(0:m+1, 0:ln+1))

  allocate(psitmp(0:m+1, 0:ln+1))

  if (.not. irrotational) then

     allocate(zettmp(0:m+1, 0:ln+1))

  end if

  allocate(vel(2,m,ln))
  allocate(rgb(3,m,ln))

!  Zero the psi array

  psi(:,:) = 0.0
  zet(:,:) = 0.0

!  Set the psi boundary condtions which are constant

   call boundarypsi(psi, m, ln, b, h, w, comm)

!  Do a boundary swap of psi

   call haloswap(psi, m, ln, comm)

   if (.not. irrotational) then

!    Update the zeta boundary condtions which depend on psi

     call boundaryzet(zet, psi, m, ln, comm)

!  Do a boundary swap of zeta

     call haloswap(zet, m, ln, comm)

  end if

!  Begin iterative Jacobi loop

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Starting main loop ...'
     write(*,*)
  end if
     
!  Barriers purely to get consistent timing - not needed for correctness

  call mpi_barrier(comm, ierr)

  tstart = mpi_wtime()

  do iter = 1, numiter

!  Compute the new psi based on the old one

     if (irrotational) then

        psitmp(1:m, 1:ln) = 0.25d0*(psi(2:m+1, 1:ln) + psi(0:m-1, 1:ln) + &
                                    psi(1:m, 2:ln+1) + psi(1:m, 0:ln-1)     )
                                    

     else


! Not clear if we should fuse these by hand or not. Probably best to
! split them as I already have "do i, do j" which is naive so why not
! be 100% naive? Not fusing will also be more like the array syntax,
! and will give more scope for compiler optimisation.

        do i = 1, m
           do j = 1, ln

              psitmp(i,j) =  0.25d0*(psi(i+1,j) + psi(i-1,j) + &
                                     psi(i,j+1) + psi(i,j-1) - &
                                     zet(i,j))

              zettmp(i,j) =  0.25d0*(zet(i+1,j) + zet(i-1,j) +     &
                                     zet(i,j+1) + zet(i,j-1)   ) - &
                            re/16.0*((psi(i,j+1)-psi(i,j-1)) *     &
                                     (zet(i+1,j)-zet(i-1,j)) -     &
                                     (psi(i+1,j)-psi(i-1,j)) *     &
                                     (zet(i,j+1)-zet(i,j-1))   )

           end do
        end do

     end if

!  Copy back

     psi(1:m, 1:ln) = psitmp(1:m, 1:ln)

     if (.not. irrotational) then

        zet(1:m, 1:ln) = zettmp(1:m, 1:ln)

     end if

!  Do a boundary swap

     call haloswap(psi, m, ln, comm)

     if (.not. irrotational) then

        call haloswap(zet, m, ln, comm)

!    Update the zeta boundary condtions which depend on psi

        call boundaryzet(zet, psi, m, ln, comm)
        
     end if

!  End iterative Jacobi loop

     if (mod(iter,1000) == 0) then
        if (rank == 0) write(*,*) 'completed iteration ', iter
     end if

  end do

  call mpi_barrier(comm, ierr)

  tstop = mpi_wtime()

  ttot  = tstop-tstart
  titer = ttot/dble(numiter)

  if (rank == 0) then

     write(*,*) 
     write(*,*) '... finished'
     write(*,*)
     write(*,fmt='('' Time for '', i6, '' iterations was '',&
          &g11.4, '' seconds'')') numiter, ttot
     write(*,fmt='('' Each individual iteration took '', g11.4, '' seconds'')') &
          titer
     write(*,*)
     write(*,*) 'Writing output file ...'

  end if

  do i = 1, m
     do j = 1, ln

        vel(1,i,j) =   (psi(i,j+1)-psi(i,j-1)) / 2.0
        vel(2,i,j) = - (psi(i+1,j)-psi(i-1,j)) / 2.0

        modvsq = vel(1,i,j)**2 + vel(2,i,j)**2
        hue = modvsq**0.4

        call hue2rgb(hue, rgb(1,i,j), rgb(2,i,j), rgb(3,i,j))

     end do
  end do

!  Output result

  call writedatafiles(rgb, vel, m, ln, scalefactor, size, rank, comm)

  if (rank == 0) then

     call writeplotfile(m, n, scalefactor)

  end if

! Finish

  call mpi_finalize(ierr)

  if (rank == 0) then
     write(*,*) ' ... finished'
     write(*,*)
     write(*,*) 'CFD completed'
     write(*,*)
  end if

end program cfd


subroutine boundarypsi(psi, m, ln, b, h, w, comm)

  use mpi

  implicit none

  integer :: m, ln, b, h, w, comm
  
  double precision, dimension(0:m+1, 0:ln+1) :: psi

  integer, dimension(MPI_STATUS_SIZE) :: status

  integer :: size, rank, ierr, i, j
  integer :: jstart, jstop

  call mpi_comm_size(comm, size, ierr)
  call mpi_comm_rank(comm, rank, ierr)

  jstart = ln*rank + 1
  jstop  = jstart + ln - 1

!  Set the boundary conditions on the bottom edge

  if (rank == 0) then

     do i = b+1, b+w-1
        psi(i, 0) = float(i-b)
     end do

     do i = b+w, m
        psi(i, 0) = float(w)
     end do

  end if

!  Set the boundary conditions on the right hand side

  do j = 1, h

     if (j .ge. jstart .and. j .le. jstop) then
        psi(m+1,j-jstart+1) = float(w)
     end if

  end do

  do j = h+1, h+w-1

     if (j .ge. jstart .and. j .le. jstop) then
        psi(m+1,j-jstart+1) = float(w-j+h)
     end if

  end do

end subroutine boundarypsi

subroutine boundaryzet(zet, psi, m, ln, comm)

  use mpi

  implicit none

  integer :: m, ln, comm
  
  double precision, dimension(0:m+1, 0:ln+1) :: zet, psi

  integer, dimension(MPI_STATUS_SIZE) :: status

  integer :: size, rank, ierr, i, j

  call mpi_comm_size(comm, size, ierr)
  call mpi_comm_rank(comm, rank, ierr)

! Set the zeta boundary conditions which depend on psi

  do j = 1, ln

     zet(0,  j) = 2.0*(psi(1,j) - psi(0,  j))
     zet(m+1,j) = 2.0*(psi(m,j) - psi(m+1,j))

  end do

  if (rank == 0) then

     do i = 1, m
        zet(i,0) = 2.0*(psi(i,  1)-psi(i,0))
     end do

  end if

  if (rank == size-1) then

     do i = 1, m
        zet(i,ln+1) = 2.0*(psi(i,ln)-psi(i,ln+1))
     end do

  end if

end subroutine boundaryzet

subroutine haloswap(x, m, n, comm)

  use mpi

  implicit none

  integer :: m, n, comm, rank, uprank, dnrank, size, ierr

  double precision :: x(0:m+1, 0:n+1)

  integer, parameter :: tag = 1

  integer, dimension(MPI_STATUS_SIZE) :: status

  call mpi_comm_size(comm, size, ierr)
  call mpi_comm_rank(comm, rank, ierr)

!  No need to halo swap for serial code

  if (size > 1) then

!  Send upper boundaries and receive lower ones
  
     if (rank == 0) then

        call mpi_send(x(1,n), m, MPI_DOUBLE_PRECISION, rank+1, &
                     tag, comm, status, ierr)

     else if (rank == size-1) then

        call mpi_recv(x(1,0), m, MPI_DOUBLE_PRECISION, rank-1, &
                      tag, comm, status, ierr)
     else

        call mpi_sendrecv(x(1,n), m, MPI_DOUBLE_PRECISION, rank+1, tag, &
                          x(1,0), m, MPI_DOUBLE_PRECISION, rank-1, tag, &
                          comm, status, ierr)

     end if
  
!  Send lower boundaries and receive upper ones

     if (rank == 0) then

        call mpi_recv(x(1,n+1), m, MPI_DOUBLE_PRECISION, rank+1, &
                      tag, comm, status, ierr)

     else if (rank == size-1) then

        call mpi_send(x(1,1), m, MPI_DOUBLE_PRECISION, rank-1, &
                     tag, comm, status, ierr)

     else

        call mpi_sendrecv(x(1,1)  , m, MPI_DOUBLE_PRECISION, rank-1, tag, &
                         x(1,n+1), m, MPI_DOUBLE_PRECISION, rank+1, tag, &
                         comm, status, ierr)
     end if

end if

end subroutine haloswap


subroutine writedatafiles(rgb, vel, m, n, scale, size, rank, comm)

  use mpi

  implicit none

  integer :: scale, comm, size, rank, m, n

  integer, dimension(3,m,n) :: rgb(3,m,n)
  double precision, dimension(2,m,n) :: vel

  integer :: irank, i, j, k,  ierr, ix, iy

  integer, parameter :: tag = 1, iounitvel = 10, iounitcol = 11

  integer, dimension(MPI_STATUS_SIZE) :: status

!  Receive data

  if (rank == 0) then

     open(unit=iounitcol, file='colourmap.dat', form='formatted')
     open(unit=iounitvel, file='velocity.dat',  form='formatted')

     do irank = 0, size-1

        if (irank /= 0) then

           call mpi_recv(rgb, 3*m*n, MPI_INTEGER, irank, tag, &
                         comm, status, ierr)

           call mpi_recv(vel, 2*m*n, MPI_DOUBLE_PRECISION, irank, tag, &
                         comm, status, ierr)

        end if

!  Write out

        do j = 1, n

           iy = irank*n + j

           do i = 1, m

              ix = i

!  Write colour map of velocity magnitude at every point

              write(iounitcol,fmt='(i4,1x,i4,1x,i3,1x,i3,1x,i3)') &
                    ix, iy, rgb(1,i,j), rgb(2,i,j), rgb(3,i,j)

!  Only write velocity vectors every "scale" points

              if (mod(ix-1,scale) == (scale-1)/2 .and. &
                  mod(iy-1,scale) == (scale-1)/2         ) then

                 write(iounitvel,fmt='(i4,1x,i4,1x,g12.5,1x,g12.5)') &
                       ix, iy, vel(1,i,j), vel(2,i,j)
              end if
             
           end do
        end do
     end do

     close(unit=iounitcol)
     close(unit=iounitvel)

  else

     call mpi_ssend(rgb, 3*m*n, MPI_INTEGER, 0, tag, comm, ierr)
     call mpi_ssend(vel, 2*m*n, MPI_DOUBLE_PRECISION, 0, tag, comm, ierr)

  end if

end subroutine writedatafiles

subroutine writeplotfile(m, n, scale)

  implicit none

  integer :: m, n, scale
  integer, parameter :: iounit = 10

  open(unit=iounit, file='cfd.plt', form='formatted')

  write(iounit,*) 'set size square'
  write(iounit,*) 'set key off'
  write(iounit,*) 'unset xtics'
  write(iounit,*) 'unset ytics'

  write(iounit,fmt='('' set xrange ['',i4,'':'',i4, '']'')') 1-scale, m+scale
  write(iounit,fmt='('' set yrange ['',i4,'':'',i4, '']'')') 1-scale, n+scale

  write(iounit,fmt='('' plot "colourmap.dat" w rgbimage, "velocity.dat" u 1:2:&
       &('',i2,''*0.75*$3/sqrt($3**2+$4**2)):&
       &('',i2,''*0.75*$4/sqrt($3**2+$4**2)) &
       &with vectors  lc rgb "#7F7F7F"'')') scale, scale

  close(unit=iounit)

end subroutine writeplotfile


subroutine hue2rgb(hue, r, g, b)

  implicit none
  
  double precision :: hue, colfunc

  integer :: r, g, b
  integer, parameter :: rgbmax = 255

  r = rgbmax*colfunc(hue-1.0)
  g = rgbmax*colfunc(hue-0.5)
  b = rgbmax*colfunc(hue    )

end subroutine hue2rgb


double precision function colfunc(x)

  implicit none

  double precision :: x, absx, val

  double precision, parameter :: x1 = 0.2, x2 = 0.5

  absx = abs(x)

  if (absx .gt. x2) then
     val = 0.0
  else if (absx .lt. x1) then
     val = 1.0
  else
     val = 1.0 - ((absx-x1)/(x2-x1))**2
  end if

  colfunc = val
      
end function colfunc
