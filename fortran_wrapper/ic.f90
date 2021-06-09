subroutine write_file(nsteps, pname, omega1, filename)

  integer                           :: nsteps
  character(8)                      :: stepString
  character(200)                    :: pname

  integer                           :: nx,ny
  character(200)                     :: fname, filename
  real(8), dimension(62001)             :: omega1

  ! auxilary variables
  integer                     :: i,j,lev
  integer                     :: ngrid,timestep,numPoints
  real(8), allocatable        :: q(:,:),omega(:,:,:),f(:)
  real(8)                     :: dx, x0, y0, time

  write(stepString, '(I5.5)') nsteps

  fname = trim(adjustl(pname))//trim(adjustl(stepString))//'.bin'

  open(101, file=trim(adjustl(fname)), access='stream', form='unformatted')

  read(101) nx
  read(101) ny
  read(101) ngrid
  read(101) dx
  read(101) x0
  read(101) y0
  read(101) numPoints

  allocate(q(ngrid,2*nx*ny+nx+ny))
  allocate(omega(ngrid,nx-1,ny-1),f(2*numPoints))

  ! read the flow field
  do lev=1,ngrid
    do i=1,2*nx*ny+nx+ny
      read(101) q(lev,i)
    enddo
  enddo

  do lev=1,ngrid
    do i=1,nx-1
      do j=1,ny-1
        read(101) omega(lev,i,j)
      enddo
    enddo
  enddo

  do i=1,2*numPoints
    read(101) f(i)
  enddo

  read(101) timestep
  read(101) time

  close(101)

  print*, 'read file ', fname
  print*, 'time and timestep are ', time, timestep

  ! change the level-1 omega values
   do j=1,ny-1
     do i=1,nx-1
       omega(1,i,j) = omega1(i + (nx-1)*(j-1))
     enddo
   enddo

  ! change the time and timestep
  timestep = 0

  ! write the file back with new file name
  open(101,file=adjustl(trim(filename)),access='stream',form='unformatted')

  write(101) nx
  write(101) ny
  write(101) ngrid
  write(101) dx
  write(101) x0
  write(101) y0
  write(101) numPoints

  do lev=1,ngrid
  do i=1,2*nx*ny+nx+ny
    write(101) q(lev,i)
  enddo
  enddo

  do lev = 1,ngrid
    do i=1,nx-1
      do j=1,ny-1
        write(101) omega(lev,i,j)
      enddo
    enddo
  enddo

  do i=1,2*numPoints
    write(101) f(i)
  enddo

  write(101) timestep
  write(101) time

  close(101)

  print*, 'written file complete ', adjustl(trim(filename))

  deallocate(omega,f,q)

end subroutine write_file

subroutine read_file(filename,omega1)
  ! Inputs: nx,ny,filename
  ! Output: omega1

  integer                           :: nx,ny
  character(200)                     :: filename
  real(8), dimension(62001)        :: omega1

  ! auxilary variables
  integer                     :: i,j,lev
  integer                     :: ngrid,timestep,numPoints
  real(8), allocatable        :: q(:,:),omega(:,:,:),f(:)
  real(8)                     :: dx, x0, y0, time

  !f2py intent(in)        :: filename
  !f2py intent(out)       :: omega1

  ! open the file to be read
  open(101,file=adjustl(trim(filename)),access='stream',form='unformatted')

  ! read the essential parameters
  read(101) nx
  read(101) ny
  read(101) ngrid
  read(101) dx
  read(101) x0
  read(101) y0
  read(101) numPoints

  ! allocate variables based on said parameters
  allocate(q(ngrid,2*nx*ny+nx+ny))
  allocate(omega(ngrid,nx-1,ny-1),f(2*numPoints))

  ! read the flow field
  do lev=1,ngrid
    do i=1,2*nx*ny+nx+ny
      read(101) q(lev,i)
    enddo
  enddo

  do lev=1,ngrid
    do i=1,nx-1
      do j=1,nx-1
        read(101) omega(lev,i,j)
      enddo
    enddo
  enddo

  do i=1,2*numPoints
    read(101) f(i)
  enddo

  read(101) timestep
  read(101) time

  close(101)

  ! we only need the level-1 of omega and q and f
  ! put it in 1-D array
  do j=1,ny-1
    do i=1,nx-1
      omega1(i + (nx-1)*(j-1)) = omega(1,i,j)
    enddo
  enddo

  print*, 'read file ', adjustl(trim(filename))

  deallocate(omega,f,q)

end subroutine read_file
