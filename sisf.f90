! 2019-11-21 by Dianmo
!
! The space Fourier transform of the self part of the van Hove correlation
! function G_s(r,t) gives the self-intermediate-scattering function F_s(q,t)
! (also called auto-correlation function).

program sisf
  implicit none
  logical :: iexist
  real, parameter :: pi = 3.14159265, dt = 0.002
  integer, parameter :: natom = 43
  integer :: iend, istep, iatom
  real(kind=8) :: rv_0(6,natom), rv_new(6,natom), qmax
  real(kind=8) :: fs, dist
  
  ! 2.7 is the radius of the first peak in rdf
  qmax = 2 * pi / 2.7

  inquire(file="md.out",exist=iexist)
  if(.not.iexist) then
    stop "Error: ""md.out"" doesn't exist!"
  endif

  inquire(file="sisf.dat",exist=iexist)
  if(iexist) then
    stop "Error: ""sisf.dat"" already exists!&
    & Rename or remove it before running this program again."
    stop
  endif

  open(10,file="md.out",form="unformatted",status="old")
  open(20,file="sisf.dat", form="formatted",status="new")
  rewind(10)
  rewind(20)
  write(20,*) "t      Fs(q,t)"
  
  read(10) istep, rv_0
  rewind(10)

  loop: do while (.true.)
    read(10,iostat=iend) istep, rv_new
    if(is_iostat_end(iend)) exit loop
    fs = 0
    
    do iatom = 1, natom
      ! TODO: periodic boundary condition
      dist = sum((rv_new(1:3,iatom)-rv_0(1:3,iatom))**2)
      dist = sqrt(dist)
      fs = fs + cos(qmax*dist)
    enddo
    fs = fs / float(natom)
    write(20,"(2f10.6)") (istep-1)*dt, fs
    if (fs < 0) exit loop
  enddo loop

end program sisf
