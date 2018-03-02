module boundary

  use mpi

  implicit none

contains

subroutine boundarypsi(psi, m, n, b, h, w, comm)

  integer :: m, n, b, h, w, comm
  
  double precision, dimension(0:m+1, 0:n+1) :: psi

  integer :: size, rank, ierr, i, j
  integer :: jstart, jstop

  call mpi_comm_size(comm, size, ierr)
  call mpi_comm_rank(comm, rank, ierr)

  jstart = n*rank + 1
  jstop  = jstart + n - 1

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

subroutine boundaryzet(zet, psi, m, n, comm)

  integer :: m, n, comm
  
  double precision, dimension(0:m+1, 0:n+1) :: zet, psi

  integer :: size, rank, ierr, i, j

  call mpi_comm_size(comm, size, ierr)
  call mpi_comm_rank(comm, rank, ierr)

! Set the zeta boundary conditions which depend on psi

  do j = 1, n

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
        zet(i,n+1) = 2.0*(psi(i,n)-psi(i,n+1))
     end do

  end if

end subroutine boundaryzet

subroutine haloswap(x, m, n, comm)

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
                     tag, comm, ierr)

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
                     tag, comm, ierr)

     else

        call mpi_sendrecv(x(1,1)  , m, MPI_DOUBLE_PRECISION, rank-1, tag, &
                         x(1,n+1), m, MPI_DOUBLE_PRECISION, rank+1, tag, &
                         comm, status, ierr)
     end if

end if

end subroutine haloswap

end module boundary
