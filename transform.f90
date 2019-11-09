! 2019-11-08 by Dianmo 
!
! This small program is used to transform the specific fragment of data
! in the unformatted file "md.out" into formatted file "md.dat".
! Data in "md.dat" will be in this order:
! x(i) y(i) z(i) vx(i) vy(i) vz(i), i stands for the i-th atom.

program transform
  implicit none
  integer, parameter :: natom = 43
  real, parameter :: dt = 0.002
  integer :: nstep, i, j
  real :: tstart, tend
  real(kind=8) :: rv(6,natom)

  open(01,file="md.out",form="unformatted",status="unknown")
  open(02,file="md.dat",form="formatted",status="unknown")
  rewind(01)
  rewind(02)

  write(*,*) "Starts at time(ps):"
  read(*,*) tstart
  write(*,*) "Ends at time(ps):"
  read(*,*) tend
  
  if (tstart > tend) then
    write(*,*) "Error: End time must be later than start time!"
    stop
  endif
  
  nstep = (tend - tstart) / dt
  
  skip: do i = 2, tstart/dt
    read(01)
  enddo skip

  do i = 0, nstep
    read(01) j, rv
    write(02,'(F10.3)') tstart
    do j = 1, natom
      write(02,"(6(SP,F12.6))") rv(1:6,j)
    enddo
    write(02,*) ''
    tstart = tstart + dt
  enddo

end program transform
        
