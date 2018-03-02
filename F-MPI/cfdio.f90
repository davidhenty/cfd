module cfdio

  use mpi

  implicit none

contains

subroutine writedatafiles(psi, m, n, scale, comm)

  integer :: m, n, scale, comm
  double precision ::  psi(0:m+1, 0:n+1)

  double precision, allocatable :: vel(:,:,:)
  integer, allocatable :: rgb(:,:,:)

  double precision :: modvsq, hue
  integer :: size, rank, irank, i, j, k,  ierr, ix, iy

  integer, parameter :: tag = 1, iounitvel = 10, iounitcol = 11

  integer, dimension(MPI_STATUS_SIZE) :: status

  call mpi_comm_size(comm, size, ierr)
  call mpi_comm_rank(comm, rank, ierr)

! Compute local velocities and colours

  allocate(rgb(3,m,n))
  allocate(vel(2,m,n))

  do i = 1, m
     do j = 1, n

        vel(1,i,j) =   (psi(i,j+1)-psi(i,j-1)) / 2.0
        vel(2,i,j) = - (psi(i+1,j)-psi(i-1,j)) / 2.0

        modvsq = vel(1,i,j)**2 + vel(2,i,j)**2
        hue = modvsq**0.4

        call hue2rgb(hue, rgb(1,i,j), rgb(2,i,j), rgb(3,i,j))

     end do
  end do

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

  double precision :: hue

  integer :: r, g, b
  integer, parameter :: rgbmax = 255

  r = rgbmax*colfunc(hue-1.0)
  g = rgbmax*colfunc(hue-0.5)
  b = rgbmax*colfunc(hue    )

end subroutine hue2rgb


double precision function colfunc(x)

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

end module cfdio
