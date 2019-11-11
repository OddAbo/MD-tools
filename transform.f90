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
  logical :: iexist
  integer :: nstep, i, j
  real :: tstart, tend
  real(kind=8) :: rv(6,natom)


  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    write(*,*) "Error: ""md.out"" doesn't exist!"
    stop
  endif

  inquire(file="md.dat",exist=iexist)
  if(iexist) then
    write(*,*) "Error: ""md.dat"" already exists!"
    stop
  endif

  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="md.dat", form="formatted",status="new")
  rewind(10)
  rewind(20)

  write(*,*) "Starts at time(ps):"
  read(*,*) tstart
  write(*,*) "Ends at time(ps):"
  read(*,*) tend
  
  if (tstart > tend) then
    write(*,*) "Error: End time must be later than start time!"
    stop
  endif
  
  nstep = floor((tend - tstart) / dt)
  
  skip: do i = 2, floor(tstart/dt)
    read(10)
  enddo skip

  do i = 1, nstep
    read(10) j, rv
    do j = 1, natom
      write(20,"(6(SP,F12.6))") rv(1:6,j)
      ! write(20,"(3(SP,F12.6))",advance="no") rv(1:3,j)
    enddo
    write(20,"")
  enddo

end program transform
